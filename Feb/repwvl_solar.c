#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include "repwvl_solar.h"


// function prototypes
void read_solar_irradiance(const char *reducedLkpPath, int nWvl, double **solar_irradiance);

void read_cloud_values(const char *reducedLkpPath, int nWvl, double **g, double **omega0, double **tau);

void get_no2_vmr(const char *reducedLkpPath, int nLev_ref, double *p_ref, double *rho_ref, double **NO2_VMR);

void get_weight(const char *reducedLkpPath, double **weight);



/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int sgn(double x) {
  if (x > 0.0) return 1;
  if (x < 0.0) return -1;
  return 0;
}

size_t LowerPos(double *tempsOnLayer, double currT, int P_REFSIZE) {
  double delta;
  size_t lowerPos=0;
  int sign;
  
  sign = sgn(tempsOnLayer[0] - currT);
  
  for (size_t p_ref_pos = 1; p_ref_pos != P_REFSIZE; ++ p_ref_pos ) {
    
    delta = tempsOnLayer[p_ref_pos] - currT;
    
    if (sign != sgn(delta)) {
      lowerPos = p_ref_pos - 1;
      return lowerPos;
    }
    else {
      if (p_ref_pos == P_REFSIZE - 1) {
	lowerPos = p_ref_pos - 1;
	return lowerPos;
      }
    }
    sign = sgn(delta);
  }

  return lowerPos;
}


void read_tau_solar (const char *reducedLkpPath, int nLev, double *pLev, double *T, double *H2O_VMR,
		     double *CO2_VMR, double *O3_VMR, double *N2O_VMR, double *CO_VMR, double *CH4_VMR,
		     double *O2_VMR, double *HNO3_VMR, double *N2_VMR, double *NO2_VMR,
		     double ***tau, double **wvl, double **weight, double **Esolar, int *nWvl) {

  const int numOfSpecies = 10;
  int nLyr=0;
  
  nLyr = nLev - 1;
  
  size_t xsec_nbooks, xsec_npages, xsec_nrows, xsec_ncols, vmrs_ref_nrows, vmrs_ref_ncols, t_pert_nelem, num_of_nodes;
  
  const double avog = 6.02214076e23;   // mol^⁻1
  const double molMassAir = 0.0289647; // (kg/mol)
  const double earthAccel = 9.80665;   // (m/s²)
  //  const double K_B = 1.38E-23; //Boltzmann const
  
  int retval;
  int ncid, xsec_nbooksid, xsec_npagesid, xsec_nrowsid, xsec_ncolsid, vmrs_ref_nrowsid, vmrs_ref_ncolsid, t_pert_nelemid, num_of_nodesid, xsecid, vmrs_refid, t_refid, p_gridid, t_pertid, ChosenWvlsid, ChosenWeightsid, c0_O3id, c1_O3id, c2_O3id, c0_NO2id, c1_NO2id, c2_NO2id;
  
  /* Open the file. NC_NOWRITE tells netCDF we want read-only access to the file */
  if ((retval = nc_open(reducedLkpPath, NC_NOWRITE, &ncid)))
    ERR(retval);
    
  if ((retval = nc_inq_dimid(ncid, "xsec_nbooks", &xsec_nbooksid)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, xsec_nbooksid, &xsec_nbooks)))
    ERR(retval);
  if ((retval = nc_inq_dimid(ncid, "xsec_npages", &xsec_npagesid)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, xsec_npagesid, &xsec_npages)))
    ERR(retval);
  if ((retval = nc_inq_dimid(ncid, "xsec_nrows", &xsec_nrowsid)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, xsec_nrowsid, &xsec_nrows)))
    ERR(retval);
  if ((retval = nc_inq_dimid(ncid, "xsec_ncols", &xsec_ncolsid)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, xsec_ncolsid, &xsec_ncols)))
    ERR(retval);
  if ((retval = nc_inq_dimid(ncid, "vmrs_ref_nrows", &vmrs_ref_nrowsid)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, vmrs_ref_nrowsid, &vmrs_ref_nrows)))
    ERR(retval);
  if ((retval = nc_inq_dimid(ncid, "vmrs_ref_ncols", &vmrs_ref_ncolsid)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, vmrs_ref_ncolsid, &vmrs_ref_ncols)))
    ERR(retval);
  if ((retval = nc_inq_dimid(ncid, "t_pert_nelem", &t_pert_nelemid)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, t_pert_nelemid, &t_pert_nelem)))
    ERR(retval);
  if ((retval = nc_inq_dimid(ncid, "f_grid_nelem", &num_of_nodesid)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, num_of_nodesid, &num_of_nodes)))
    ERR(retval);
    
    
  *nWvl = num_of_nodes;

  // !!! VMRS is upside down !!!
  double VMRS[numOfSpecies][nLyr];
  for (int i=0; i<nLyr; i++) {
    VMRS[0][nLyr-1-i] = H2O_VMR[i];
    VMRS[1][nLyr-1-i] = H2O_VMR[i];
    VMRS[2][nLyr-1-i] = CO2_VMR[i];
    VMRS[3][nLyr-1-i] = O3_VMR[i];
    VMRS[4][nLyr-1-i] = N2O_VMR[i];
    VMRS[5][nLyr-1-i] = CO_VMR[i];
    VMRS[6][nLyr-1-i] = CH4_VMR[i];
    VMRS[7][nLyr-1-i] = O2_VMR[i];
    VMRS[8][nLyr-1-i] = HNO3_VMR[i];
    VMRS[9][nLyr-1-i] = N2_VMR[i];
  } 
  
  size_t startp[4], countp[4];
  
  startp[0] = 0;
  startp[1] = 0;
  startp[2] = 0;
  startp[3] = 0;
  countp[0] = xsec_nbooks;
  countp[1] = xsec_npages;
  countp[2] = 1;
  countp[3] = xsec_ncols;
  
  double xsec[xsec_nbooks][xsec_npages][xsec_ncols];
  double vmrs_ref[vmrs_ref_nrows][vmrs_ref_ncols];
  double t_pert[t_pert_nelem];
  double press[vmrs_ref_ncols];
  double t_ref[vmrs_ref_ncols];
  double wvl_array[num_of_nodes];
  double weights_array[num_of_nodes];
  double c0_O3_array[xsec_nrows];
  double c1_O3_array[xsec_nrows];
  double c2_O3_array[xsec_nrows];
  double c0_NO2_array[xsec_nrows];
  double c1_NO2_array[xsec_nrows];
  double c2_NO2_array[xsec_nrows];
    
  if ((retval = nc_inq_varid(ncid, "xsec", &xsecid)))
    ERR(retval);
    
    
  if ((retval = nc_inq_varid(ncid, "vmrs_ref", &vmrs_refid)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, vmrs_refid, &vmrs_ref[0][0])))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "t_ref", &t_refid)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, t_refid, &t_ref[0])))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "p_grid", &p_gridid)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, p_gridid, &press[0])))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "t_pert", &t_pertid)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, t_pertid, &t_pert[0])))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "ChosenWvls", &ChosenWvlsid)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, ChosenWvlsid, &wvl_array[0])))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "ChosenWeights", &ChosenWeightsid)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, ChosenWeightsid, &weights_array[0])))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "c0_O3", &c0_O3id)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, c0_O3id, &c0_O3_array[0])))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "c1_O3", &c1_O3id)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, c1_O3id, &c1_O3_array[0])))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "c2_O3", &c2_O3id)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, c2_O3id, &c2_O3_array[0])))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "c0_NO2", &c0_NO2id)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, c0_NO2id, &c0_NO2_array[0])))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "c1_NO2", &c1_NO2id)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, c1_NO2id, &c1_NO2_array[0])))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "c2_NO2", &c2_NO2id)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, c2_NO2id, &c2_NO2_array[0])))
    ERR(retval);
  
  
  *wvl    = (double *)  calloc (*nWvl, sizeof(double));
  *weight = (double *)  calloc (*nWvl, sizeof(double));
  *tau    = (double **) calloc (*nWvl, sizeof(double *));
  for (int iwvl=0; iwvl<*nWvl; iwvl++)
    (*tau)[iwvl] = (double *) calloc (nLyr, sizeof(double));
  
  
  for (int iwvl=0; iwvl<num_of_nodes; iwvl++) {
    (*wvl)[iwvl] = wvl_array[iwvl];
    (*weight)[iwvl] = weights_array[iwvl];
  } 
  
  double tempsOnLayer[t_pert_nelem];
  double numDens[nLyr];
  
  double tempXsec;
  
  size_t lowPosT;
  size_t lowPosP;
  double c0, cP, cT, cPT;
  double delP, delT;
  
  double midT, midP, midVMR;
  
  double tauMat[nLyr][xsec_nrows];
  
  for (int i = 0; i<nLyr; ++i)
    for (int j = 0; j!= xsec_nrows; ++j)
      tauMat[i][j] = 0;

  // !!! numDens is upside down !!! 
  for (int lyr = 0; lyr<nLyr; ++lyr)
    numDens[lyr] = (pLev[nLev - 1 - lyr] - pLev[nLev - 2 - lyr]) * avog / molMassAir / earthAccel;
  
  //Calculation of optical thickness
  for (int iwvl = 0; iwvl != xsec_nrows; ++iwvl) {
    startp[2] = iwvl;
    
    if ((retval = nc_get_vara_double (ncid, xsecid, startp, 
				      countp, &xsec[0][0][0])))
      ERR(retval);
    
    for (int lyr = 0; lyr < nLyr; ++lyr) {
      midT = T[nLyr - 1 - lyr];

      double sigma_O3 = (c0_O3_array[iwvl] + c1_O3_array[iwvl]*(midT - 273.15) + c2_O3_array[iwvl]*(midT - 273.15)*(midT - 273.15)) * 1E-20;

      // hydrostatic eqilibrium
      double tau_O3 = O3_VMR[nLev - 2 - lyr] * sigma_O3 * (pLev[nLev - 1 - lyr] - pLev[nLev - 2 - lyr]) * avog / molMassAir / earthAccel * 1E-4; 

      // original formula by Sophie Meier
      // double numDensExact =  (pLev[nLev - 1 - lyr] + pLev[nLev - 2 - lyr])/(2*K_B * T[nLyr - 1 - lyr])*1E-6;
      // double tau_O3_org = numDensExact * O3_VMR[nLev - 2 - lyr] * sigma_O3 * (zLev[nLev - 2 - lyr] - zLev[nLev - 1 - lyr])*1E5; 
      // fprintf (stdout, "SIGMA %d %d %e %e %e    %e %e\n", iwvl, lyr, sigma_O3, tau_O3, tau_O3_org);
      
      tauMat[nLev - 2 - lyr][iwvl] += tau_O3;
        
      double sigma_NO2 = (c0_NO2_array[iwvl] + c1_NO2_array[iwvl]*(midT - 273.15) + c2_NO2_array[iwvl]*(midT - 273.15)*(midT - 273.15)) * 1E-20;
      double tau_NO2 = NO2_VMR[nLev - 2-lyr] * sigma_NO2 * (pLev[nLev - 1 - lyr]- pLev[nLev - 2 - lyr]) * avog / molMassAir / earthAccel * 1E-4;
            
      tauMat[nLev - 2 - lyr][iwvl] += tau_NO2;
      
      for (int spec = 0; spec != numOfSpecies; ++spec) {
	
        int spec2 = spec;
        midP = (pLev[nLev - 2 - lyr] + pLev[nLev - 1 - lyr]) / 2;

	midVMR = VMRS[spec][lyr];
        
        lowPosP = LowerPos(press, midP, vmrs_ref_ncols);
        
        for (int pertNum = 0; pertNum != t_pert_nelem; ++pertNum)
	  tempsOnLayer[pertNum] = t_ref[lowPosP] + t_pert[pertNum];

        lowPosT = LowerPos(tempsOnLayer,midT, t_pert_nelem);
        
        
        c0 = xsec[lowPosT][spec][lowPosP];
        cT = xsec[lowPosT + 1][spec][lowPosP] - c0;
        cP = xsec[lowPosT][spec][lowPosP + 1] - c0;
        cPT = xsec[lowPosT + 1][spec][lowPosP + 1] -cP - cT - c0;
        delT = (midT- tempsOnLayer[lowPosT]) / (tempsOnLayer[lowPosT + 1] - tempsOnLayer[lowPosT]);
        delP = (midP - press[lowPosP]) / (press[lowPosP + 1] - press[lowPosP]);
        
        tempXsec = (c0 + cT * delT + cP * delP + cPT * delT * delP)*midVMR*numDens[lyr];
        
        if(spec == 0) {
        
	  c0 = xsec[lowPosT][spec2][lowPosP]*vmrs_ref[0][lowPosP]*numDens[lyr]*(midVMR*midVMR)/(vmrs_ref[0][lowPosP]*vmrs_ref[0][lowPosP]);
	  cT = xsec[lowPosT + 1][spec2][lowPosP]*vmrs_ref[0][lowPosP]*numDens[lyr]*(midVMR*midVMR)/(vmrs_ref[0][lowPosP]*vmrs_ref[0][lowPosP]) - c0;
	  cP = xsec[lowPosT][spec2][lowPosP + 1]*vmrs_ref[0][lowPosP+1]*numDens[lyr]*(midVMR*midVMR)/(vmrs_ref[0][lowPosP+1]*vmrs_ref[0][lowPosP+1]) - c0;
	  cPT = xsec[lowPosT + 1][spec2][lowPosP + 1]*vmrs_ref[0][lowPosP+1]*numDens[lyr]*(midVMR*midVMR)/(vmrs_ref[0][lowPosP+1]*vmrs_ref[0][lowPosP+1])  -cP - cT - c0;
	  delT = (midT- tempsOnLayer[lowPosT]) / (tempsOnLayer[lowPosT + 1] - tempsOnLayer[lowPosT]);
	  delP = (midP - press[lowPosP]) / (press[lowPosP + 1] - press[lowPosP]);
            
	  tempXsec = (c0 + cT * delT + cP * delP + cPT * delT * delP);
        }

        tauMat[nLev - 2 - lyr][iwvl] += tempXsec;
      }
    }
  }
  
  for (int iwvl=0; iwvl<*nWvl; iwvl++)
    for (int ilyr=0; ilyr<nLyr; ilyr++)
      (*tau)[iwvl][ilyr] = tauMat[ilyr][iwvl];

  read_solar_irradiance(reducedLkpPath, *nWvl, Esolar);
}


void read_solar_irradiance(const char *reducedLkpPath, int nWvl, double **solar_irradiance) {
    
  int ncid, retval, wvlSolarid, solIrrid;
  size_t numWvlSolar;
    
  if ((retval = nc_open(reducedLkpPath, NC_NOWRITE, &ncid)))
    ERR(retval);
    
  if ((retval = nc_inq_dimid(ncid, "solar_irradiance_nelem", &wvlSolarid)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, wvlSolarid, &numWvlSolar)))
    ERR(retval);
    
  double irradianceArray[numWvlSolar];
    
  /* Get the varid of the data variable, based on its name. */
  if ((retval = nc_inq_varid(ncid, "solar_irradiance", &solIrrid)))
    ERR(retval);
  
  /* Read the data. */
  if ((retval = nc_get_var_double(ncid, solIrrid, &irradianceArray[0])))
    ERR(retval);
  
    
    
  if(numWvlSolar != nWvl) {
    printf("Error, number of wavelengths inconsistent %d : %ld", nWvl, numWvlSolar);
  }
    
  *solar_irradiance    = (double *)  calloc (numWvlSolar, sizeof(double));
    
  for (int iwvl=0; iwvl<nWvl; iwvl++)
    (*solar_irradiance)[iwvl] = irradianceArray[iwvl];  
}




void read_cloud_values(const char *reducedLkpPath, int nWvl, double **g, double **omega0, double **tau) {
    
  int ncid, retval;
    
  if ((retval = nc_open(reducedLkpPath, NC_NOWRITE, &ncid)))
    ERR(retval);
    
  int nuWvlsid;
  size_t numWvlTau;
    
  if ((retval = nc_inq_dimid(ncid, "number_tau_values", &nuWvlsid)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, nuWvlsid, &numWvlTau)))
    ERR(retval);
    
  if(numWvlTau != nWvl) {
    printf("Error, number of wavelengths inconsistent; %d : %ld", nWvl, numWvlTau);
  }
    
    
  int gid, omega0id, tauid;
  double tau_array[numWvlTau];
  double gArray[numWvlTau];
  double omegArray[numWvlTau];
    
  if ((retval = nc_inq_varid(ncid, "gValue_cloud", &gid)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, gid, &gArray[0])))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "omega0_cloud", &omega0id)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, omega0id, &omegArray[0])))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "tau_cloud", &tauid)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, tauid, &tau_array[0])))
    ERR(retval);

    
  *g      = (double *)  calloc (numWvlTau, sizeof(double));
  *omega0 = (double *)  calloc (numWvlTau, sizeof(double));
  *tau    = (double *)  calloc (numWvlTau, sizeof(double));
  
    
  for (int iwvl=0; iwvl<numWvlTau; iwvl++) {
    (*tau)[iwvl] = tau_array[iwvl]/156.3361;
    (*g)[iwvl] = gArray[iwvl];
    (*omega0)[iwvl] = omegArray[iwvl];
  }
    
    
    
}



void get_no2_vmr(const char *reducedLkpPath, int nLev_ref, double *p_ref, double *rho_ref, double **NO2_VMR) {
    
  int ncid, retval;
    
  if ((retval = nc_open(reducedLkpPath, NC_NOWRITE, &ncid)))
    ERR(retval);
    
  int levid, numDensid, pressid;
  size_t nlev_atm;
    
  if ((retval = nc_inq_dimid(ncid, "reference_NO2_numberDensity_nLev", &levid)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, levid, &nlev_atm)))
    ERR(retval);
    
  double press_no2[nlev_atm];
  double numDens_no2[nlev_atm];
    
  if ((retval = nc_inq_varid(ncid, "pressure_reference_NO2", &pressid)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, pressid, &press_no2[0])))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "reference_NO2_numberDensity", &numDensid)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, numDensid, &numDens_no2[0])))
    ERR(retval);
   
    
    
  double numDens_interpolated[nLev_ref];
  *NO2_VMR   = (double *)  calloc (nLev_ref -1, sizeof(double));
    
  for (int ilev_ref=0; ilev_ref < nLev_ref; ilev_ref++) {
    int ilev_atm = 0;
    while (ilev_atm < nlev_atm && press_no2[ilev_atm] < p_ref[ilev_ref]) {
      ilev_atm += 1;
    }
    numDens_interpolated[ilev_ref] = (numDens_no2[ilev_atm-1] + (numDens_no2[ilev_atm]-numDens_no2[ilev_atm-1])/(press_no2[ilev_atm]-press_no2[ilev_atm-1])*(p_ref[ilev_ref]-press_no2[ilev_atm-1]));
        
  }
    
  for (int ilev_ref=0; ilev_ref < nLev_ref-1; ilev_ref++) {
    (*NO2_VMR)[ilev_ref] = (numDens_interpolated[ilev_ref] + numDens_interpolated[ilev_ref+1])/(2*rho_ref[ilev_ref]);
  }
    
}

void get_weight(const char *reducedLkpPath, double **weight) {
    
  int ncid, retval;
    
  if ((retval = nc_open(reducedLkpPath, NC_NOWRITE, &ncid)))
    ERR(retval);
    
  int levid, weightid;
  size_t num_of_nodes;
    
  if ((retval = nc_inq_dimid(ncid, "numOfNodes", &levid)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, levid, &num_of_nodes)))
    ERR(retval);
    
  double weights_array[num_of_nodes];
    
  if ((retval = nc_inq_varid(ncid, "ChosenWeights", &weightid)))
    ERR(retval);
  if ((retval = nc_get_var_double(ncid, weightid, &weights_array[0])))
    ERR(retval);

  *weight = (double *)  calloc (num_of_nodes, sizeof(double));

  for (int iwvl=0; iwvl<num_of_nodes; iwvl++)
    (*weight)[iwvl] = weights_array[iwvl];
}
