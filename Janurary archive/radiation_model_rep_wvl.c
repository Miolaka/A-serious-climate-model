#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ascii.h"
#include <unistd.h>
#include "repwvl_thermal.h"
#include "repwvl_solar.h"

#define k_B 1.3806503e-23
#define h 6.62607015e-34
#define c 299792458.0

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

// ############ functions ###################################
// ##########################################################

/* calculate layer properties t, r, rdir, sdir, and tdir from            */
/* layer optical thickness dtau, asymmetry parameter g,                  */
/* single scattering albedo omega0, and cosine of solar zenith angle mu0 */

void eddington_v2(double dtau, double g, double omega0, double mu0,
                  double *t, double *r, double *rdir, double *sdir, double *tdir)
{
    double alpha1 = 0, alpha2 = 0, alpha3 = 0, alpha4 = 0, alpha5 = 0, alpha6 = 0;
    double a11 = 0, a12 = 0, a13 = 0, a23 = 0, a33 = 0;
    double lambda = 0, b = 0, A = 0;
    double denom = 0;

    /* first, avoid omega0=1 because of instability */
    if (omega0 > 0.999999)
        omega0 = 0.999999;

    alpha1 = (1.0 - omega0) + 0.75 * (1.0 - omega0 * g);
    alpha2 = -(1.0 - omega0) + 0.75 * (1.0 - omega0 * g);

    lambda = sqrt(alpha1 * alpha1 - alpha2 * alpha2);

    A = 1.0 / (alpha2 / (alpha1 - lambda) * exp(lambda * dtau) - alpha2 / (alpha1 + lambda) * exp(-lambda * dtau));

    a11 = A * 2.0 * lambda / alpha2;
    a12 = A * (exp(lambda * dtau) - exp(-lambda * dtau));

    b = 0.5 - 0.75 * g * mu0;
    alpha3 = -omega0 * b;
    alpha4 = omega0 * (1 - b);

    denom = (1.0 / mu0 / mu0 - lambda * lambda);
    alpha5 = ((alpha1 - 1.0 / mu0) * alpha3 - alpha2 * alpha4) / denom;
    alpha6 = (alpha2 * alpha3 - (alpha1 + 1.0 / mu0) * alpha4) / denom;

    a33 = exp(-dtau / mu0);

    a13 = alpha5 * (1.0 - (a11) * (a33)) - alpha6 * (a12);
    a23 = -(a12)*alpha5 * (a33) + alpha6 * ((a33) - (a11));

    *t = a11;
    *r = a12;
    *tdir = a33;
    *rdir = a13 / mu0;
    *sdir = a23 / mu0;
}

// ##########################################################
// ############ variables ###################################
// ##########################################################

int main(int argc, char **argv)
{

    double sigma = 5.670374419e-8;
    double sigma_sca_mol = 4.84e-31; // unit is 1/m^2
    double g = 9.80665;
    double albedo = 0.30933137398971344;
    double cp = 1004.0; //[J/kgK]
    double E_solar = 340.25 * (1 - albedo);
    double lambda;
    double dlambda;
    double M_air = 0.02896; // unit is kg/mol
    double rho_air = 1.204; // unit is kg/m^3
    double N_A = 6.02214086e23;

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

    double g_asym_param_mol; //=0

    // ##########################################################

    // ##########################################################
    // ############ vertikales Profil und Spurengase einlesen ###
    // ############ Variablen für thermal einlesen ##############
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
    double *NO2_VMR = NULL;

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
    NO2_VMR = calloc(nlev, sizeof(double));

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
        NO2_VMR[ilev] = 10E-12;
    }

    /* convert pressure from hPa to Pa */
    for (int ilev = 0; ilev < nlev; ilev++)
    {
        p[ilev] = plev[ilev] * 100.0;
    }

    /* obtain spectral absorption coefficients */

    double *weight = NULL;

    // check if tautot array is filled with values
    // printf("%7.9f", tautot[0][0]);

    // #########################################################
    // ############ Variablen declaren für #####################
    // ############  representative wvl. solar #################
    // #########################################################

    int nwvl_s = 0;
    double *wvl_s = NULL;
    double *edir = NULL, *edir_lambda = NULL;

    edir = calloc(nlev, sizeof(double));
    edir_lambda = calloc(nlev, sizeof(double));

    /* obtain spectral absorption coefficients */

    double *weight_s;
    double *E0_s;

    char repwvlfilename[FILENAME_MAX] = "./ReducedLookupFile_solar_50wvls.nc";

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

    for (int j = 0; j < 8000; j++)
    {

        // SOLAR WAVELENGTH PART
        double **tau_SW_abs_mol = NULL; // read in

        // double **tau_SW_abs_cloud[nlay] = NULL;
        // double **tau_SW_sca_cloud[nlay] = NULL;
        // double **tau_SW_ext_cloud[nlay] = NULL;

        // include also function for the solar tau values and short wavelengths
        read_tau_solar(repwvlfilename, nlev, plevPa, T,
                       H2O_VMR, CO2_VMR, O3_VMR, N2O_VMR, CO_VMR, CH4_VMR, O2_VMR, HNO3_VMR, N2_VMR, NO2_VMR,
                       &tau_SW_abs_mol, &wvl_s, &weight_s, &E0_s, &nwvl_s);

        double tau_SW_sca_mol[nwvl_s];       // calculate by lambda^4 formula
        double tau_SW_ext_mol[nwvl_s][nlay]; // calculate via single scattering albedo formula

        // start with short wavelength loop
        g_asym_param_mol = 0;
        // mu0 = 0.7;
        double t[nlay];
        double r[nlay];
        double rdir[nlay];
        double sdir[nlay];
        double tdir[nlay];

        double R[nlay];
        double Trans[nlay];
        double T_dir[nlay];
        double S_dir[nlay];
        double sigma_sca_mol;
        double omega0[nwvl_s][nlay];

        double E_down_SW[nlev];
        double E_up_SW[nlev];
        double E_down_SW_lambda[nlev];
        double E_up_SW_lambda[nlev];

        // calculate the solar irradiance dependant on the wavelengths and absorption + weight
        // reset edir, E_down_SW, E_up_SW

        double mu0 = 1; // 60 degree solar zenith angle

        for (int i = 0; i < nlev; i++)
        {
            E_down_SW[i] = 0.0;
            E_up_SW[i] = 0.0;
            edir[i] = 0;
        }

        for (int w = 0; w < nwvl_s; w++)
        {

            lambda = wvl_s[w];
            sigma_sca_mol = 4.84e-31 * pow(550 / lambda, 4);
            tau_SW_sca_mol[w] = (sigma_sca_mol * delta_p * N_A) / (M_air * g);

            printf("%d %e \n", w, tau_SW_sca_mol[0]);

            for (int i = 0; i <= nlay; i++)
            {
                tau_SW_ext_mol[w][i] = tau_SW_sca_mol[w] + tau_SW_abs_mol[w][i];
                omega0[w][i] = tau_SW_sca_mol[w] / tau_SW_ext_mol[w][i];

                printf("%e \n", tau_SW_ext_mol[0][i]);

                eddington_v2(tau_SW_ext_mol[w][i], g_asym_param_mol, omega0[w][i], mu0,
                             &(t[i]), &(r[i]), &(rdir[i]), &(sdir[i]), &(tdir[i])); // is array filled by this? yes

                // calculate Ri+1, Ti+1, Sdiri+1

                R[0] = r[0];
                R[i + 1] = r[i + 1] + (R[i] * pow(t[i + 1], 2)) / (1.0 - R[i] * r[i + 1]);

                Trans[0] = t[0];
                Trans[i + 1] = (Trans[i] * t[i + 1]) / (1.0 - R[i] * r[i + 1]);

                T_dir[0] = 1;
                T_dir[i] = T_dir[i] * tdir[i];

                S_dir[0] = sdir[0];
                S_dir[i + 1] = (t[i + 1] * S_dir[0] + T_dir[i] * rdir[i + 1] * R[i] * t[i + 1] + T_dir[i] * sdir[i + 1]) / (1.0 - R[i] * r[i + 1]);
            }

            /*   for (int k = 0; k < nlay; k++) {
               printf("%d %e %e %e %e %e\n", w, t[k], r[k], rdir[k], sdir[k], tdir[k]);
           }

           */

            // calculate irradiances at levels starting from bottom (i=nlev) to top (i=0)
            // at the surface n ( = nlev): lowermost level
            // with Ag being the ground Albedo defined for total = direct + diffuse irradiance
            double Ag = 0.3;

            edir_lambda[0] = E0_s[w] * mu0 * weight_s[w];

            for (int i = 0; i <= nlev - 1; i++)
            {
                edir_lambda[i + 1] = edir_lambda[i] * exp(-tau_SW_abs_mol[w][i] / mu0);
            }

            E_down_SW_lambda[nlev] = (S_dir[nlev - 1] + T_dir[nlev - 1] * R[nlev - 1] * Ag) * edir_lambda[0] / (1.0 - R[nlev - 1] * Ag);
            E_up_SW_lambda[nlev] = Ag * (E_down_SW_lambda[nlev] + T_dir[nlev - 1] * edir_lambda[0]);

            for (int i = nlev; i > 1; i--)
            {

                E_down_SW_lambda[i - 1] = (R[i - 2] * E_down_SW_lambda[i] + edir_lambda[0] * S_dir[i - 2] + edir_lambda[i - 1] * rdir[i - 1] * R[i - 2]) / (1.0 - R[i - 2] * r[i - 1]);
                E_up_SW_lambda[i - 1] = (t[i - 1] * E_up_SW_lambda[i] + edir_lambda[0] * S_dir[i - 2] * r[i - 1] + edir_lambda[i - 1] * rdir[i - 1]) / (1.0 - R[i - 2] * r[i - 1]);
            }

            // boundary conditions: calculate irradiances at the surface (n = nlev)
            // top of the atmosphere: uppermost level

            E_down_SW_lambda[0] = 0;
            E_up_SW_lambda[0] = t[0] * E_up_SW_lambda[1] + rdir[0] * edir_lambda[0];

            for (int i = nlev; i >= 0; i--)
            {

                E_down_SW[i] = E_down_SW[i] + E_down_SW_lambda[i];
                E_up_SW[i] = E_up_SW[i] + E_up_SW_lambda[i];
                edir[i] += edir_lambda[i];
            }
        }

        // THERMAL WAVELENGTH PART

        // determine optical thickness from quantities at levels
        // include in time loop because the temperature values change during each run
        double **tau_LW_abs_mol = NULL; // read in
        double tau_LW_sca_mol[nlay];    // = 0
        double tau_LW_ext_mol[nlay];    // equals absorption

        // double **tau_LW_abs_cloud[nlay] = NULL;
        // double **tau_LW_sca_cloud[nlay] = NULL;
        // double **tau_LW_ext_cloud[nlay] = NULL;

        read_tau("./ReducedLookupFile_thermal_50wvls_corrected.nc", nlev, p, T,
                 H2O_VMR, CO2_VMR, O3_VMR, N2O_VMR, CO_VMR, CH4_VMR, O2_VMR, HNO3_VMR, N2_VMR, // pressure p is always level quantity
                 &tau_LW_abs_mol, &wvl, &weight, &nwvl, 0);                                    // nlev = 20 (layer length)

        // Initialize E_up and E_down arrays to zero
        for (int i = 0; i < nlev; i++)
        {
            E_up[i] = 0.0;
            E_down[i] = 0.0;
        }

        /*
            alpha = 0.000;
            epsilon = 0.000;
        for (double h = 1; h<= 1000; h++ ) {

            alpha = h/1000;
            epsilon = alpha;


            printf("%7.5f %7.5f %7.5f\n", epsilon, T[9], Theta[9]);
        }
        */

        // #############################################################
        //  Start with wavelength loop over lambda (long wave)##########
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
                    epsilon[i] = calc_epsilon(tau_LW_abs_mol[w][i], u);
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

        double E_down_total[nlev];
        double E_up_total[nlev];

        for (int i = 0; i <= (nlev - 2); i++)
        {
            E_down_total[i] = E_down[i] + E_down_SW[i] + edir[i];
            E_up_total[i] = E_up[i] + E_up_SW[i];
        }

        // ################ Delta_E #####################################
        // Radiation differences
        // lowermost layer radiation difference

        delta_E[nlay - 1] = E_down_total[nlev - 2] - E_up_total[nlev - 2];

        for (int i = nlay - 2; i >= 0; i--)
        {
            delta_E[i] = E_up_total[i + 1] + E_down_total[i] - E_down_total[i + 1] - E_up_total[i];
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

        /*

             //print new values
            printf("new values after dt" "\n");
            for (int k = 0; k < nlay; k++) {
                printf("%d %d %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f\n", j, k, (p[k] + p[k+1])/2, T[k], Theta[k], E_down_SW[k], E_up_SW[k], edir[k]);
            }

         */
    }

    /*

    for (int i=0; i<nwvl; i++){

        free(tau_LW_abs_mol[i]);
    }


    free(tau_LW_abs_mol);

    */

    return 0;
}
