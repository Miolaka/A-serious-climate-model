#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include "repwvl_thermal.h"

/* Handle errors by printing an error message and exiting with a
  * non-zero status. */
 #define ERRCODE 2
 #define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

static int sgn(double x) {
  if (x > 0.0) return 1;
  if (x < 0.0) return -1;
  return 0;
}

static size_t LowerPos(double *tempsOnLayer, double currT, int P_REFSIZE){
  double delta;
  size_t lowerPos=0;
  int sign;
  
  sign = sgn(tempsOnLayer[0] - currT);
  
  for(size_t p_ref_pos = 1; p_ref_pos != P_REFSIZE; ++ p_ref_pos ) {
    
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



void read_tau (const char *reducedLkpPath, int nLev, double *pLev, double *T, double *H2O_VMR,
	       double *CO2_VMR, double *O3_VMR, double *N2O_VMR, double *CO_VMR, double *CH4_VMR,
	       double *O2_VMR, double *HNO3_VMR, double *N2_VMR,
	       double ***tau, double **wvl, double **weight, int *nWvl,
	       int prop_at_Lev) {
  //Optionally we could use default values for some trace gases and reduce input.
  const int numOfSpecies = 10;
  int nProp=0;
  
  if (prop_at_Lev==0)  // temperature and mixing ratios for layers
    nProp = nLev - 1;
  else                 // temperature and mixing ratios for levels
    nProp = nLev;
  
  
  size_t xsec_nbooks, xsec_npages, xsec_nrows, xsec_ncols, vmrs_ref_nrows, vmrs_ref_ncols, t_pert_nelem, num_of_nodes, nls_pert_nelem;
  
  const double avog = 6.02214076e23; // mol^⁻1
  const double molMassAir = 0.0289647; // (kg/mol)
  const double earthAccel = 9.80665; // (m/s²)
  
  
  int retval;
  int ncid, xsec_nbooksid, xsec_npagesid, xsec_nrowsid, xsec_ncolsid, vmrs_ref_nrowsid, vmrs_ref_ncolsid, t_pert_nelemid, num_of_nodesid, nls_pert_nelemid, nls_pertid, xsecid, vmrs_refid, t_refid, p_gridid, t_pertid, ChosenWvlsid, ChosenWeightsid;
  
  /* Open the file. NC_NOWRITE tells netCDF we want read-only access
     * to the file.*/
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
    if ((retval = nc_inq_dimid(ncid, "numOfNodes", &num_of_nodesid)))
        ERR(retval);
    if ((retval = nc_inq_dimlen(ncid, num_of_nodesid, &num_of_nodes)))
        ERR(retval);
    if ((retval = nc_inq_dimid(ncid, "nls_pert_nelem", &nls_pert_nelemid)))
        ERR(retval);
    if ((retval = nc_inq_dimlen(ncid, nls_pert_nelemid, &nls_pert_nelem)))
        ERR(retval);
    
    *nWvl = num_of_nodes;
        
    
  double currentPress[nLev];
  double currentTemp[nProp];
  /*
  double H2O_VMR_VEC[nProp];
  double CO2_VMR_VEC[nProp];
  double O3_VMR_VEC[nProp];
  double N2O_VMR_VEC[nProp];
  double CO_VMR_VEC[nProp];
  double CH4_VMR_VEC[nProp];
  double O2_VMR_VEC[nProp];
  double HNO3_VMR_VEC[nProp];
  double N2_VMR_VEC[nProp];
  */
  
  for (int i=0; i<nLev; i++){
    currentPress[nLev-1-i] = pLev[i];  
  } 
  for (int i=0; i<nProp; i++){
    currentTemp[nProp-1-i] = T[i];
    /*
    H2O_VMR_VEC[nProp-1-i] = H2O_VMR[i];
    CO2_VMR_VEC[nProp-1-i] = CO2_VMR[i];
    O3_VMR_VEC[nProp-1-i] = O3_VMR[i];
    N2O_VMR_VEC[nProp-1-i] = N2O_VMR[i];
    CO_VMR_VEC[nProp-1-i] = CO_VMR[i];
    CH4_VMR_VEC[nProp-1-i] = CH4_VMR[i];
    O2_VMR_VEC[nProp-1-i] = O2_VMR[i];
    HNO3_VMR_VEC[nProp-1-i] = HNO3_VMR[i];
    N2_VMR_VEC[nProp-1-i] = N2_VMR[i];
    */
  }
    
  double VMRS[numOfSpecies][nProp];
  for(int i=0; i<nProp; i++){
    VMRS[0][nProp-1-i] = H2O_VMR[i];
    VMRS[1][nProp-1-i] = H2O_VMR[i];
    VMRS[2][nProp-1-i] = CO2_VMR[i];
    VMRS[3][nProp-1-i] = O3_VMR[i];
    VMRS[4][nProp-1-i] = N2O_VMR[i];
    VMRS[5][nProp-1-i] = CO_VMR[i];
    VMRS[6][nProp-1-i] = CH4_VMR[i];
    VMRS[7][nProp-1-i] = O2_VMR[i];
    VMRS[8][nProp-1-i] = HNO3_VMR[i];
    VMRS[9][nProp-1-i] = N2_VMR[i];
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
  
  double nls_pert[nls_pert_nelem];
  double xsec[xsec_nbooks][xsec_npages][xsec_ncols];
  double vmrs_ref[vmrs_ref_nrows][vmrs_ref_ncols];
  double t_pert[t_pert_nelem];
  double press[vmrs_ref_ncols];
  double t_ref[vmrs_ref_ncols];
  double wvl_array[num_of_nodes];
  double weights_array[num_of_nodes];
    
  
  /* Get the varid of the data variable, based on its name. */
    if ((retval = nc_inq_varid(ncid, "nls_pert", &nls_pertid)))
       ERR(retval);
  
    /* Read the data. */
    if ((retval = nc_get_var_double(ncid, nls_pertid, &nls_pert[0])))
       ERR(retval);
  
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
  
  
  *wvl    = (double *)  calloc (*nWvl, sizeof(double));
  *weight = (double *)  calloc (*nWvl, sizeof(double));
  *tau    = (double **) calloc (*nWvl, sizeof(double *));
  for (int iwvl=0; iwvl<*nWvl; iwvl++)
    (*tau)[iwvl] = (double *) calloc (nLev-1, sizeof(double));
  
  
  for(int iwvl=0; iwvl<num_of_nodes; iwvl++){
    (*wvl)[iwvl] = wvl_array[iwvl];
    (*weight)[iwvl] = weights_array[iwvl];
  } 
  
  double tempsOnLayer[t_pert_nelem];
  double numDens[nLev-1];
  
  double tempXsec;
  
  size_t lowPosT;
  size_t lowPosP;
  double c0, cP, cT, cPT;
  double delP, delT;
  
  double midT, midP, midVMR;
  
  double tauMat[nLev-1][xsec_nrows];
  
  for (int i = 0; i != nLev - 1; ++i)
    for (int j = 0; j!= xsec_nrows; ++j)
      tauMat[i][j] = 0;
  
  for(int lyr = 0; lyr != nLev - 1; ++lyr){
    numDens[lyr] = (currentPress[lyr] - currentPress[lyr + 1]) * avog / molMassAir / earthAccel;
  }
  
  
  
  //Calculation of optical thickness
  for (int wvl = 0; wvl != xsec_nrows; ++wvl) {
    startp[2] = wvl;
    
    if ((retval = nc_get_vara_double(ncid, xsecid, startp, 
                       countp, &xsec[0][0][0])))
      ERR(retval);
    
    for(int lyr = 0; lyr < nLev - 1; ++lyr){
      
      for(int spec = 0; spec != numOfSpecies; ++spec){
	
        int spec2 = spec;
        midP = (currentPress[lyr + 1] + currentPress[lyr]) / 2;

        if(prop_at_Lev==0) {
            midT = currentTemp[lyr];
            midVMR = VMRS[spec][lyr];
        }
        else {
            midT = (currentTemp[lyr + 1] + currentTemp[lyr]) / 2;
            midVMR = (VMRS[spec][lyr + 1] + VMRS[spec][lyr]) / 2;
        }
        
        lowPosP = LowerPos(press, midP, vmrs_ref_ncols);
        
        
        for(int pertNum = 0; pertNum != t_pert_nelem; ++pertNum){
            tempsOnLayer[pertNum] = t_ref[lowPosP] + t_pert[pertNum];
        }
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
        

        tauMat[nLev - 2 - lyr][wvl] += tempXsec;
    }
    }
    
      
    
  }
  
  
  
  
  
  for (int iwvl=0; iwvl<*nWvl; iwvl++)
    for (int ilyr=0; ilyr<nLev-1; ilyr++)
      (*tau)[iwvl][ilyr] = tauMat[ilyr][iwvl];
  
  
    
  
}
  
  
  
  
  
