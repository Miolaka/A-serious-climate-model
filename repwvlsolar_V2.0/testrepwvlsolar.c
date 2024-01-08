#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "ascii.h"
#include "repwvl_solar.h"

int main ()
{
  int nlev=0, nwvl=0;

  char atmfilename[FILENAME_MAX]="./test.atm";

  int status = 0;
  double *plev=NULL, *plevPa=NULL, *Tlev=NULL, *Tlyr=NULL, *rholev=NULL;
  double *wvl=NULL;
  double *zlev=NULL;

  double *H2O_VMR=NULL,  *H2O_VMR_lyr=NULL;
  double *CO2_VMR=NULL,  *CO2_VMR_lyr=NULL;
  double *O3_VMR=NULL,   *O3_VMR_lyr=NULL;
  double *N2O_VMR=NULL,  *N2O_VMR_lyr=NULL;
  double *CO_VMR=NULL,   *CO_VMR_lyr=NULL;
  double *CH4_VMR=NULL,  *CH4_VMR_lyr=NULL;
  double *O2_VMR=NULL,   *O2_VMR_lyr=NULL;
  double *HNO3_VMR=NULL, *HNO3_VMR_lyr=NULL;
  double *N2_VMR=NULL,   *N2_VMR_lyr=NULL;
  double *NO2_VMR=NULL,  *NO2_VMR_lyr=NULL;
    
  double *edir=NULL, *edir_lambda=NULL;

  /* first read atmospheric profile */
  status = read_9c_file (atmfilename, &zlev, &plev, &Tlev, &rholev, &H2O_VMR, &O3_VMR, &CO2_VMR, &CH4_VMR, &N2O_VMR, &nlev);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, atmfilename);
    return status;
  }
    
  fprintf (stderr, " ... read %d levels from %s\n", nlev, atmfilename);

  O2_VMR   = calloc (nlev, sizeof(double));
  N2_VMR   = calloc (nlev, sizeof(double));
  NO2_VMR  = calloc (nlev, sizeof(double));
  HNO3_VMR = calloc (nlev, sizeof(double));
  CO_VMR   = calloc (nlev, sizeof(double));

  plevPa   = calloc (nlev, sizeof(double));

  /* scaling atmospheric trace gases and defining extra ones */
  for (int ilev=0; ilev<nlev; ilev++) {
    H2O_VMR [ilev] *= 1E-6;
    O3_VMR  [ilev] *= 1E-6;
    CO2_VMR [ilev] *= 1E-6;
    CH4_VMR [ilev] *= 1E-6;
    N2O_VMR [ilev] *= 1E-6;
    O2_VMR  [ilev] = 0.2095;
    N2_VMR  [ilev] = 0.7808;
    NO2_VMR [ilev] = 10E-12;
    HNO3_VMR[ilev] = 0.0;
    CO_VMR  [ilev] = 0.0;

    plevPa  [ilev] = plev[ilev] * 100.0;
  }

    
  /* layer quantities */
  int nlyr = nlev-1;

  Tlyr  = calloc (nlyr, sizeof(double));
  
  H2O_VMR_lyr  = calloc (nlyr, sizeof(double));
  CO2_VMR_lyr  = calloc (nlyr, sizeof(double));
  O3_VMR_lyr   = calloc (nlyr, sizeof(double));
  N2O_VMR_lyr  = calloc (nlyr, sizeof(double));
  CO_VMR_lyr   = calloc (nlyr, sizeof(double));
  CH4_VMR_lyr  = calloc (nlyr, sizeof(double));
  O2_VMR_lyr   = calloc (nlyr, sizeof(double));
  HNO3_VMR_lyr = calloc (nlyr, sizeof(double));
  N2_VMR_lyr   = calloc (nlyr, sizeof(double));
  NO2_VMR_lyr  = calloc (nlyr, sizeof(double));
  
  for (int ilyr=0; ilyr<nlyr; ilyr++)  {
    Tlyr[ilyr] = (Tlev[ilyr] + Tlev[ilyr+1])/2.0;
    
    H2O_VMR_lyr [ilyr] = (H2O_VMR[ilyr]  + H2O_VMR[ilyr+1])/2.0;
    CO2_VMR_lyr [ilyr] = (CO2_VMR[ilyr]  + CO2_VMR[ilyr+1])/2.0;
    O3_VMR_lyr  [ilyr] = (O3_VMR[ilyr]   + O3_VMR[ilyr+1])/2.0;
    N2O_VMR_lyr [ilyr] = (N2O_VMR[ilyr]  + N2O_VMR[ilyr+1])/2.0;
    CO_VMR_lyr  [ilyr] = (CO_VMR[ilyr]   + CO_VMR[ilyr+1])/2.0;
    CH4_VMR_lyr [ilyr] = (CH4_VMR[ilyr]  + CH4_VMR[ilyr+1])/2.0;
    O2_VMR_lyr  [ilyr] = (O2_VMR[ilyr]   + O2_VMR[ilyr+1])/2.0;
    HNO3_VMR_lyr[ilyr] = (HNO3_VMR[ilyr] + HNO3_VMR[ilyr+1])/2.0;
    N2_VMR_lyr  [ilyr] = (N2_VMR[ilyr]   + N2_VMR[ilyr+1])/2.0;
    NO2_VMR_lyr [ilyr] = (NO2_VMR[ilyr]  + NO2_VMR[ilyr+1])/2.0;
  }

  edir          = calloc (nlev, sizeof(double));
  edir_lambda   = calloc (nlev, sizeof(double));
  

  /* obtain spectral absorption coefficients */
  double **tau;
  double *weight;
  double *E0;

  char repwvlfilename[FILENAME_MAX] = "./ReducedLookupFile_solar_50wvls.nc";

  fprintf (stderr, " ... read_tau_solar (\"%s\")\n", repwvlfilename);
  fflush (stderr);
  
  fprintf (stderr, " ... calculating tau from %d layers\n", nlyr);
  read_tau_solar (repwvlfilename, nlev, plevPa, Tlyr,
		  H2O_VMR_lyr, CO2_VMR_lyr, O3_VMR_lyr, N2O_VMR_lyr, CO_VMR_lyr, CH4_VMR_lyr, O2_VMR_lyr, HNO3_VMR_lyr, N2_VMR_lyr, NO2_VMR_lyr,
		  &tau, &wvl, &weight, &E0, &nwvl);

  
  /* reset edir */
  for (int ilev=0; ilev<nlev; ilev++) 
    edir[ilev] = 0;

  double mu0 = 1;  // 60 degree solar zenith angle
  
  /* radiative transfer calculation; calculate direct solar irradiance */
  for (int iv=0; iv<nwvl; iv++) {   // wavelength loop */
    edir_lambda[0] = E0[iv] * mu0;
    
    for (int ilev=0; ilev<nlev-1; ilev++)
      edir_lambda[ilev+1] = edir_lambda[ilev] * exp(-tau[iv][ilev]/mu0);

    /* sum irradiance */
    for (int ilev=0; ilev<nlev; ilev++)
      edir[ilev] += edir_lambda[ilev] * weight[iv];
  }
  
  /* print result for initial atmosphere and stop */
  for (int ilev=0; ilev<nlev; ilev++)
    fprintf (stdout, "%2d %7.1f %7.1f %5.1f %7.1f\n", 
	     ilev, zlev[ilev], plev[ilev], Tlev[ilev], edir[ilev]);
  
  return 0;
}
