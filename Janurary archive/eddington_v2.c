
/* calculate layer properties t, r, rdir, sdir, and tdir from            */
/* layer optical thickness dtau, asymmetry parameter g,                  */
/* single scattering albedo omega0, and cosine of solar zenith angle mu0 */

void eddington_v2 (double dtau, double g, double omega0, double mu0,
		   double *t, double *r, double *rdir, double *sdir, double *tdir)
{
  double alpha1=0, alpha2=0, alpha3=0, alpha4=0, alpha5=0, alpha6=0;
  double a11=0, a12=0, a13=0, a23=0, a33=0;
  double lambda=0, b=0, A=0;
  double denom=0;

  /* first, avoid omega0=1 because of instability */
  if (omega0 > 0.999999)
    omega0=0.999999;

  alpha1= (1.0-omega0)+0.75*(1.0-omega0*g);
  alpha2=-(1.0-omega0)+0.75*(1.0-omega0*g);
  
  lambda=sqrt(alpha1*alpha1-alpha2*alpha2);
  
  A=1.0/(alpha2/(alpha1-lambda)*exp(lambda*dtau)-alpha2/(alpha1+lambda)*exp(-lambda*dtau));
  
  a11=A*2.0*lambda/alpha2;
  a12=A*(exp(lambda*dtau)-exp(-lambda*dtau));
  
  b=0.5-0.75*g*mu0;
  alpha3=-omega0*b; 
  alpha4=omega0*(1-b);
  
  denom = (1.0/mu0/mu0-lambda*lambda);
  alpha5=((alpha1-1.0/mu0)*alpha3-alpha2*alpha4)/denom;
  alpha6=(alpha2*alpha3-(alpha1+1.0/mu0)*alpha4)/denom;
  
  a33=exp(-dtau/mu0);
  
  a13=alpha5*(1.0-(a11)*(a33))-alpha6*(a12);
  a23=-(a12)*alpha5*(a33)+alpha6*((a33)-(a11));

  *t    = a11;
  *r    = a12;
  *tdir = a33;
  *rdir = a13 / mu0;
  *sdir = a23 / mu0;
}
