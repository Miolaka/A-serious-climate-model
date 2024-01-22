#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ascii.h"
#include <unistd.h>
#include "repwvl_thermal.h"

#define k_B 1.3806503e-23
#define h 6.62607015e-34
#define c 299792458.0
#define sigma 5.670374419e-8
#define g 9.80665
#define albedo 0.30933137398971344
#define cp 1004.0


// ##########################################################
// ############ functions ###################################
// ##########################################################

void sort_vertically(double Theta[], int nlay, int unstable)
{
    unstable = 1;
    while (unstable)
    {
        unstable = 0;
        for (int i = nlay - 2; i >= 0; i--)
        { // 9,8,7,6,5,...,0
            if (Theta[i] < Theta[i + 1])
            {
                unstable = 1;
                double tmp = Theta[i + 1];
                Theta[i + 1] = Theta[i];
                Theta[i] = tmp;
            }
        }
    }
}

void get_Theta_from_T(double Theta[], double T[], int nlay, double p[])
{
    for (int i = nlay - 1; i >= 0; i--)
    {
        Theta[i] = T[i] * pow((100000.0 / ((p[i] + p[i + 1]) / 2)), 2.0 / 7.0);
    }
}

void get_T_from_Theta(double Theta[], double T[], int nlay, double p[])
{
    for (int i = nlay - 1; i >= 0; i--)
    {
        T[i] = Theta[i] * pow((((p[i] + p[i + 1]) / 2.0) / 100000.0), 2.0 / 7.0);
    }
}

double calc_epsilon(double tau, double mu)
{
    return 1.0 - exp(-tau / mu);
}

double calc_B(double T, double lambda, double wavelength)
{

    wavelength = wavelength * 1e-9;                                                                                            /* convert from nm to m */
    return (2.0 * h * c * c / (lambda * lambda * lambda * lambda * lambda)) / (exp(h * c / (lambda * k_B * T)) - 1.0) / 1.0e9; /* W/(m2 nm sterad) */
}

// ##########################################################
// ############ variables ###################################
// ##########################################################

void geometic_optic(double incident, double transmittance, double reflectance){

}
/*
use double addington to determine T and R
*/
void double_addington(){
    //get molecular scattering 
    //get molecular absorption
    //get omega from Rayleigh
    delta_E[]
}

int main(int argc, char **argv)
{
    double E_solar = 340.25 * (1 - albedo);
    double lambda;
    double dlambda;
    double R = 287.0;
    double wavelength;

    int unstable;
    int convection = 1;
    int convection_counter = 0;
    int nlay = 20;
    int nlev = nlay + 1;
    double T[nlay];
    double Theta[nlay];
    double p[nlev];
    double E_up[nlev];
    double E_down[nlev];
    double L_up[nlev];
    double L_down[nlev];
    double delta_E[nlay];
    double dT[nlay];
    double epsilon[nlay];
    double alpha[nlay];
    double init_srf_T = 288;
    double init_top_T = 100;
    double delta_T = (init_srf_T - init_top_T) / nlay;
    double B[nlay];
    double mu;

    double dt = 1000;
    double dTmax = 5;
    double du = 0.1;
    double delta_p = 100000 / nlay;

    // ##########################################################
    // ############ get the tau(lambda, höhe) ###################
    // ############ höhe entspricht layer (i) ###################
    // ##########################################################

    int nwvl = 0, status = 0;

    double *plev = NULL, *plevPa = NULL, *Tlev = NULL, *rholev = NULL;
    double *wvl = NULL;
    double *tmp = NULL;
    double Bg = 0;

    double *H2O_VMR = NULL;
    double *CO2_VMR = NULL;
    double *O3_VMR = NULL;
    double *N2O_VMR = NULL;
    double *CO_VMR = NULL;
    double *CH4_VMR = NULL;
    double *O2_VMR = NULL;
    double *HNO3_VMR = NULL;
    double *N2_VMR = NULL;

    /* first read atmospheric profile */
    status = read_9c_file("test.atm", &tmp, &plev, &Tlev, &rholev, &H2O_VMR, &O3_VMR, &CO2_VMR, &CH4_VMR, &N2O_VMR, &nlev);
    if (status != 0)
    {
        fprintf(stderr, "Error %d reading %s\n", status, "test.atm");
        return status;
    }

    O2_VMR = calloc(nlev, sizeof(double));
    N2_VMR = calloc(nlev, sizeof(double));
    HNO3_VMR = calloc(nlev, sizeof(double));
    CO_VMR = calloc(nlev, sizeof(double));
    plevPa = calloc(nlev, sizeof(double));

    for (int ilev = 0; ilev < nlev - 1; ilev++)
    {
        /* convert from ppm to absolute concentrations */
        H2O_VMR[ilev] = (H2O_VMR[ilev] + H2O_VMR[ilev + 1]) / 2 * 1E-6;
        O3_VMR[ilev] = (O3_VMR[ilev] + O3_VMR[ilev + 1]) / 2 * 1E-6;
        CO2_VMR[ilev] = (CO2_VMR[ilev] + CO2_VMR[ilev + 1]) / 2 * 1E-6;
        CH4_VMR[ilev] = (CH4_VMR[ilev] + CH4_VMR[ilev + 1]) / 2 * 1E-6;
        N2O_VMR[ilev] = (N2O_VMR[ilev] + N2O_VMR[ilev + 1]) / 2 * 1E-6;

        /* define missing species */
        O2_VMR[ilev] = 0.2095;
        N2_VMR[ilev] = 0.7808;
        HNO3_VMR[ilev] = 0.0;
        CO_VMR[ilev] = 0.0;
    }

    /* convert pressure from hPa to Pa */
    for (int ilev = 0; ilev < nlev; ilev++)
    {
        p[ilev] = plev[ilev] * 100.0;
    }

    /* obtain spectral absorption coefficients */

    double **tau = NULL;
    double *weight = NULL;

    // check if tautot array is filled with values
    // printf("%7.9f", tautot[0][0]);

    // ##########################################################
    // ############ Initialize the vertical profile #############
    // ##########################################################

    //  layers
    //  use Tlev as starting profile
    for (int i = 0; i < nlev - 1; i++)
    {
        T[i] = (Tlev[i] + Tlev[i + 1]) / 2;
    }

    get_Theta_from_T(Theta, T, nlay, p);

    printf("Unsorted vertical layers\n");
    for (int i = 0; i < nlay; i++)
    {
        printf("%d %7.1f %7.1f %7.1f\n", i, (p[i] + p[i + 1]) / 2, T[i], Theta[i]);
    }

    // #############################################################
    //  Now sort the Theta profile, so it is increasing with height#
    //  (Theta_0 > Theta_1 > ... > Theta_10 (ground level))#########
    // #############################################################

    sort_vertically(Theta, nlay, unstable);
    get_T_from_Theta(Theta, T, nlay, p); // Fragen - whrsl. egal

    // Print the sorted Theta values
    printf("Sorted vertical layers\n");
    for (int i = 0; i < nlay; i++)
    {
        printf("%d %7.1f %7.1f %7.1f\n", i, (p[i] + p[i + 1]) / 2, T[i], Theta[i]);
    }

    // #############################################################
    //  Start with the time loop####################################
    // #############################################################

    // determine optical thickness from quantities at levels
    // include in time loop because the temperature values change during each run
    read_tau("./ReducedLookupFile_thermal_50wvls_corrected.nc", nlev, p, T,
             H2O_VMR, CO2_VMR, O3_VMR, N2O_VMR, CO_VMR, CH4_VMR, O2_VMR, HNO3_VMR, N2_VMR, // pressure p is always level quantity
             &tau, &wvl, &weight, &nwvl, 0);                                               // nlev = 20 (layer length)

    /*
        alpha = 0.000;
        epsilon = 0.000;
    for (double h = 1; h<= 1000; h++ ) {

        alpha = h/1000;
        epsilon = alpha;


        printf("%7.5f %7.5f %7.5f\n", epsilon, T[9], Theta[9]);
    }
    */

    // time loop
    for (int j = 0; j < 8000; j++)
    {

        // Initialize E_up and E_down arrays to zero
        for (int i = 0; i < nlev; i++)
        {
            E_up[i] = 0.0;
            E_down[i] = 0.0;
        }

        // #############################################################
        //  Start with wavelength loop over lambda #####################
        // #############################################################

        for (int w = 0; w < nwvl; w++)
        {

            lambda = wvl[w] * 1e-9;

            // Calculate B plancksche Strahlung:

            for (int i = nlev - 2; i >= 0; i--)
            {
                B[i] = weight[w] * calc_B(T[i], lambda, wavelength);
            }

            // #############################################################
            // ############ Start loop over angles (mu) ####################
            // #############################################################

            for (double u = 0.05; u < 1; u += 0.1)
            {

                // printf("%7.2f, %7.7f""\n", u,  alpha);

                // #############################################################
                // ############ Calc E_up, E_down and Delta_E###################
                // #############################################################

                // calc epsilon-array:

                for (int i = 0; i < nlay; i++)
                {
                    epsilon[i] = calc_epsilon(tau[w][i], u);
                    alpha[i] = epsilon[i];
                }

                // ################ E_Up #######################################
                // upward radiation
                // ground lvl upward radiation fixed

                // Randbedingung 1: surface level
                L_up[nlev - 1] = B[nlev - 2]; // epsilon=1

                // remaining lvls
                for (int i = nlev - 2; i >= 0; i--)
                {

                    L_up[i] = L_up[i + 1] * (1 - alpha[i]) + epsilon[i] * B[i];
                    E_up[i] = E_up[i] + L_up[i] * u * du * 2.0 * M_PI;
                }

                // ################ E_DOWN #####################################
                // downward radiation
                // Randbedingung 2: uppermost level fixed
                L_down[0] = 0;

                // remaining lvls
                for (int i = 0; i <= (nlev - 2); i++)
                {

                    L_down[i + 1] = L_down[i] * (1.0 - alpha[i]) + epsilon[i] * B[i];
                    E_down[i + 1] = E_down[i + 1] + L_down[i + 1] * u * du * 2.0 * M_PI;
                }

            } // mu loop closes
        }     // lambda loop closes

        // ################ Delta_E #####################################
        // Radiation differences
        // lowermost layer radiation difference
        delta_E[nlay - 1] = E_down[nlev - 2] - E_up[nlev - 2] + E_solar;

        for (int i = nlay - 2; i >= 0; i--)
        {
            delta_E[i] = E_up[i + 1] + E_down[i] - E_down[i + 1] - E_up[i];
        }

        // #############################################################
        // ###### Temperature and Theta increase -> convection #########
        // #############################################################

        // get dt from Max-heating rate

        double dT_dt_max = 0;
        double dT_dt = 0;

        for (int i = 0; i <= nlay - 1; i++)
        {
            dT_dt = (delta_E[i] * g / ((p[i + 1] - p[i]) * cp));

            printf("%7.10f\n", dT_dt);

            if (fabs(dT_dt) > dT_dt_max)
            {
                dT_dt_max = fabs(dT_dt);
            }
        }

        dt = dTmax / dT_dt_max;

        if (dt > 30000)
        {
            dt = 30000;
        }

        printf("new dt value is:"
               "\n");
        printf("%7.1f %7.10f\n ", dt, dT_dt_max);

        // calculate temperature differences
        for (int i = 0; i <= nlay - 1; i++)
        {
            dT[i] = (delta_E[i] * g / ((p[i + 1] - p[i]) * cp)) * dt;
            T[i] = T[i] + dT[i];
            // Theta[i] = Theta[i] + dT[i];
        }

        // calculate new Theta values, convection, get new T
        get_Theta_from_T(Theta, T, nlay, p);
        sort_vertically(Theta, nlay, unstable);
        get_T_from_Theta(Theta, T, nlay, p);

        // print new values
        printf("new values after dt"
               "\n");
        for (int k = 0; k < nlay; k++)
        {
            printf("%d %d %7.1f %7.1f %7.1f %7.1f %7.1f\n", j, k, (p[k] + p[k + 1]) / 2, T[k], Theta[k], E_down[k], E_up[k]);
        }
    }

    double q;
    double r_droplet;
    double E_pc = q/4/pi/epsilon/r_droplet/r_droplet;
    double L_ray;
    L_ray [i] = delta_E[i]*delta_E[i];
    double tau_sca[];
    double tau_abs[];
    for(int i=0; i<nlay; i++){
        double beta_sca[i]
        tau_sca[i] = beta_sca[i]*delta_z[i];

        double omega_0 = tau_sca[i]/(tau_sca[i]+tau_abs[i]); 

    }

    double rho; //density
    double M_a = 28.97; //Molar mass
    double u = 1.66e-27
    double n_n = rho/M_a/n_n;



    for (int i = 0; i < nwvl; i++)
    {

        free(tau[i]);
    }

    free(tau);

    return 0;
}

