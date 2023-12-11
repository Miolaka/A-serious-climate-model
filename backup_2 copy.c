#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ascii.h"




//##########################################################
//############ variables ###################################
//##########################################################




int main(int argc, char **argv) {
    int unstable;
    int convection = 1;
    int convection_counter = 0;
    int nlayer = 20;
    int nlevels = nlayer + 1;
    float T[nlayer];
    float Theta[nlayer];
    float p[nlevels];
    float E_up[nlevels];
    float E_down[nlevels];
    float L_up[nlevels];
    float L_down[nlevels];
    float delta_E[nlayer];
    float dT[nlayer];
    float g=9.80665;
    float E_solar=235;
    float sigma = 5.670374419e-8;
    float lambda;
    float epsilon[nlayer];
    float alpha[nlayer];
    float cp=1005; //[J/kgK]
    float dt;
    float delta_p = 100000 / nlayer;
    float init_srf_T = 288;
    float init_top_T = 100;
    float delta_T = (init_srf_T - init_top_T) / nlayer;
    float B[nlayer];
    float tau;
    float mu;
    float dTmax=5;
    float du=0.1;




    float h=6.62e-34;
    float c=3e8;
    float k_B=1.38e-23;
    float dlambda;
   


   
//##########################################################
//############ functions ###################################
//##########################################################




void sort_vertically(float Theta[], int nlayer) {
    unstable = 1;
    while (unstable) {
        unstable = 0;
        for (int i = nlayer - 2; i >= 0; i--) { // 9,8,7,6,5,...,0
            if (Theta[i] < Theta[i + 1]) {
                unstable = 1;
                float tmp = Theta[i + 1];
                Theta[i + 1] = Theta[i];
                Theta[i] = tmp;    //temperature swap/heating/cooling?
            }
        }
    }
}




void get_Theta_from_T(float Theta[], float T[]) {
    for (int i = nlayer - 1; i >= 0; i--) {
        Theta[i] = T[i] * pow((100000 / ((p[i]+p[i+1])/2 )), 2.0 / 7.0);
        }
    }




void get_T_from_Theta(float Theta[], float T[]) {
    for (int i = nlayer - 1; i >= 0; i--) {
        T[i] = Theta[i] * pow((((p[i]+p[i+1])/2 )/ 100000), 2.0 / 7.0);
        }
    }




inline float calc_epsilon(float tau, float mu) {
    return 1 - exp(- tau / mu);
}








float calc_B(float T, float lambda) {
    return (2*h*c*c/(lambda*lambda*lambda*lambda*lambda)) * 1/(exp(h*c/(lambda*k_B*T))-1);
}












//##########################################################
//############ get the tau(lambda, höhe) ###################
//############ höhe entspricht layer (i) ###################
//##########################################################




  int nwvl=0, nlyr=0;
  int status=0;
   
  char tauCO2filename[FILENAME_MAX]="./lbl.co2.asc";
  char tauH2Ofilename[FILENAME_MAX]="./lbl.h2o.asc";
  char tauCH4filename[FILENAME_MAX]="./lbl.ch4.asc";
  char tauN2Ofilename[FILENAME_MAX]="./lbl.n2o.asc";
  char tauO3filename[FILENAME_MAX]="./lbl.o3.asc";
   


  double *wvl=NULL;        // 1D array
  double **tauCO2=NULL;    // 2D array
  double **tauH2O=NULL;    // 2D array
  double **tauCH4=NULL;    // 2D array
  double **tauN2O=NULL;    // 2D array
  double **tauO3=NULL;    // 2D array
  double **tautot=NULL;    // 2D array


  /* read CO2 optical thickness profile */
  status = ASCII_file2xy2D (tauCO2filename,  
                &nwvl, &nlyr,
                &wvl, &tauCO2);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, tauCO2filename);
    return status;
  }


  fprintf (stderr, " ... read %d wavelengths and %d layers from %s\n",
       nwvl, nlyr, tauCO2filename);


 
  /* read H2O optical thickness profile */
  status = ASCII_file2xy2D (tauH2Ofilename,  
                &nwvl, &nlyr,
                &wvl, &tauH2O);


  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, tauH2Ofilename);
    return status;
  }


  fprintf (stderr, " ... read %d wavelengths and %d layers from %s\n",
       nwvl, nlyr, tauH2Ofilename);
   
    /* read CH4  optical thickness profile */
  status = ASCII_file2xy2D (tauCH4filename,  
                &nwvl, &nlyr,
                &wvl, &tauCH4);


  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, tauCH4filename);
    return status;
  }


  fprintf (stderr, " ... read %d wavelengths and %d layers from %s\n",
       nwvl, nlyr, tauCH4filename);
   
   
    /* read N2O optical thickness profile */
  status = ASCII_file2xy2D (tauN2Ofilename,  
                &nwvl, &nlyr,
                &wvl, &tauN2O);


  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, tauN2Ofilename);
    return status;
  }


  fprintf (stderr, " ... read %d wavelengths and %d layers from %s\n",
       nwvl, nlyr, tauN2Ofilename);
   
     /* read O3 optical thickness profile */
  status = ASCII_file2xy2D (tauO3filename,  
                &nwvl, &nlyr,
                &wvl, &tauO3);


  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, tauO3filename);
    return status;
  }


  fprintf (stderr, " ... read %d wavelengths and %d layers from %s\n",
       nwvl, nlyr, tauO3filename);






  /* add all gases to obtain total optical thickness */


  /* first allocate memory */
  tautot=calloc(nwvl, sizeof(double));
  for (int iv=0; iv<nwvl; iv++)
    tautot[iv] = calloc(nlyr, sizeof(double));
   
  for (int iv=0; iv<nwvl; iv++)
    for (int ilyr=0; ilyr<nlyr; ilyr++)
      tautot[iv][ilyr] = tauCO2[iv][ilyr] + tauH2O[iv][ilyr] + tauCH4[iv][ilyr] + tauN2O[iv][ilyr] + tauO3[iv][ilyr] ;


//check if tautot array is filled with values
//printf("%7.9f", tautot[0][0]);




//##########################################################
//############ Initialize the vertical profile #############
//##########################################################
   


//  levels
    for(int i=0; i<nlevels; i++){ //0,1,2,...,10
        p[i] = i*delta_p;
    }
       




//  layers
    for (int i = 0; i < nlayer; i++) {
        T[i] = 288.00 - (nlayer-i) * delta_T;
    }




    get_Theta_from_T(Theta,T);




printf("Unsorted vertical layers\n");
    for (int i = 0; i < nlayer; i++) {
        printf("%d %7.1f %7.1f %7.1f\n", i, (p[i] + p[i+1])/2, T[i], Theta[i]);
    }




//#############################################################
// Now sort the Theta profile, so it is increasing with height#
// (Theta_0 > Theta_1 > ... > Theta_10 (ground level))#########
//#############################################################




    sort_vertically(Theta,  nlayer);
    get_T_from_Theta(Theta, T); //Fragen - whrsl. egal
   
    // Print the sorted Theta values
    printf("Sorted vertical layers\n");
    for (int i = 0; i < nlayer; i++) {
        printf("%d %7.1f %7.1f %7.1f\n", i, (p[i] + p[i+1])/2, T[i], Theta[i]);
    }












//#############################################################
// Start with the time loop####################################
//#############################################################




/*
    alpha = 0.000;
    epsilon = 0.000;
for (float h = 1; h<= 1000; h++ ) {
   
    alpha = h/1000;
    epsilon = alpha;




    printf("%7.5f %7.5f %7.5f\n", epsilon, T[9], Theta[9]);
}
*/




//time loop
for(int j=0; j<800; j++){




// Initialize E_up and E_down arrays to zero
for (int i = 0; i < nlevels; i++) {
    E_up[i] = 0.0;
    E_down[i] = 0.0;
}








//#############################################################
// Start with wavelength loop over lambda #####################
//#############################################################




for (int w = 0; w<nwvl-1; w++) {


    lambda = wvl[w+1]*1e-9;
    dlambda = (wvl[w+1] - wvl[w])*1e-9;




//#############################################################
//############ Start loop over angles (mu) ####################
//#############################################################




//Calculate B plancksche Strahlung:


    for(int i=nlevels-2; i>=0; i--){
        B[i] = calc_B(T[i], lambda);
     }




for (float u=0.05; u <= 0.95; u += 0.1){


   
    //printf("%7.2f, %7.7f""\n", u,  alpha);


//#############################################################
//############ Calc E_up, E_down and Delta_E###################
//#############################################################
//################ E_Up #######################################


//calc epsilon-array:


    for(int i=nlevels-2; i>=0; i--){
        epsilon[i] = calc_epsilon(tautot[w][i], u);
        alpha[i] = epsilon[i];
     }




//upward radiation
//ground lvl upward radiation fixed , nlevels=11, nlayer=10








    //Randbedingung 1: surface level
    L_up[nlevels-1] = B[nlevels-2]*dlambda; //epsilon=1




    //remaining lvls
    for(int i=nlevels-2; i>=0; i--){


        L_up[i] = L_up[i+1] * (1 - alpha[i]) + epsilon[i] * B[i] *dlambda;
        E_up[i] = E_up[i] + L_up[i] * u * du * 2 * M_PI;
     }




//################ E_DOWN #####################################
    //downward radiation
    //Randbedingung 2: uppermost level fixed
    L_down[0] = 0;




    //remaining lvls  
    for (int i=0; i<=(nlevels-2); i++){


        L_down[i+1] = L_down[i] * (1 - alpha[i]) +  epsilon[i] * B[i] * dlambda;
        E_down[i+1] = E_down[i+1] + L_down[i+1] * u * du * 2 * M_PI;
    }




} // mu loop closes
} // lambda loop closes




//################ Delta_E #####################################
    //Radiation differences
    //lowermost layer radiation difference
    delta_E[nlayer-1] = E_down[nlevels-2] - E_up[nlevels-2] + E_solar;




    for (int i= nlayer-2; i >= 0 ; i--){
        delta_E[i] = E_up[i+1] + E_down[i] - E_down[i+1] - E_up[i];




    }




//#############################################################
//###### Temperature and Theta increase -> convection #########
//#############################################################


    // get dt from Max-heating rate
   
    float dT_dt_max = 0;
    float dT_dt = 0;


    for (int i=0; i<=nlayer-1; i++){
        dT_dt = (delta_E[i] * g / ((p[i+1]-p[i])*cp));


        printf("%7.10f\n", dT_dt);
       
 
        if (fabs(dT_dt) > dT_dt_max ){
            dT_dt_max = fabs(dT_dt);
        }
    }


    dt = dTmax / dT_dt_max;


    printf("new dt value is:" "\n");
    printf("%7.1f %7.10f\n", dt, dT_dt_max);


   






    //calculate temperature differences
    for (int i=0; i<=nlayer-1; i++){
        dT[i] = (delta_E[i] * g / ((p[i+1]-p[i])*cp)) * dt;
        T[i] = T[i] + dT[i];
        Theta[i] = Theta[i] + dT[i];




    }




    //calculate new Theta values
    sort_vertically(Theta, nlayer);
    get_T_from_Theta(Theta, T);












     //print new values
    printf("new values after dt" "\n");
    for (int k = 0; k < nlayer; k++) {
        printf("%d %d %7.1f %7.1f %7.1f\n", j, k, (p[k] + p[k+1])/2, T[k], Theta[k]);
    }
 
 
 }  




for (int i=0; i<nwvl; i++){
   
    free(tautot[i]);
}




free(tautot);


    return 0;
}