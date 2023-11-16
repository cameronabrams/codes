
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_sf_bessel.h>
      
double * Rlam;  /* array of R*lambda_i, where J0(R*lambda_i) = 0 */
double * cfacs; /* array of factors in summation which do not 
		   depend on r or t, and can therefore be precomputed */
int n_evalues;  /* number of eigenvalues */

/* Initializes both the eigenvalue array and the precomputed factor array */
void init ( double R, double ra ) {
  int i=0;

  Rlam = (double*)malloc(n_evalues*sizeof(double));
  cfacs = (double*)malloc(n_evalues*sizeof(double));
  for (i=0;i<n_evalues;i++){
    Rlam[i] = gsl_sf_bessel_zero_J0(i+1);
    cfacs[i] = gsl_sf_bessel_J1(ra/R*Rlam[i]) / (Rlam[i] / R) 
      / pow(gsl_sf_bessel_J1(Rlam[i]),2); 
  }
}

/* Computes and returns the analytic c(r,t) */
double c_analytic ( double r, double t, double ra, 
		    double R, double D, double c_0) {
  int j;
  double c = 0.0;
  for (j=0;j<n_evalues;j++){
    c += exp(-Rlam[j]*Rlam[j]*D*t/R/R) * gsl_sf_bessel_J0(r/R*Rlam[j]) 
      * cfacs[j];
  }
  c *= 2*c_0*ra/R/R;
  return c;
}


int main (int argc, char * argv[] ) 
{
  int i,j,k;
  double r;          // radial coordinate (m)
  double t;          // time (s)
  double t_0,t_f,dt; // initial and final time values, and time step
  double r_i,r_f,dr; // inner and outer radii, and radial increment
  double ra,c_0;     // initial condition: c(r<ra,0) = c_0 [=] g/m^3
  double R;          // radial position of zero-concentration boundary (m)
  double D;          // diffusivity (m^2/s) 
  int log=0;         // flag for time update; if log, t = t*dt 
                     //   instead of t = t+dt


  /* Default values */
  D = 1.e-10;  /* 1.e-6 cm^2/s */
  R = 1.e-3;   /* 1 mm */
  ra = 1.e-5;  /* 0.01 mm  = 10 micron */
  c_0 = 1.e6;  /* g/m^3 = 1.e3 mg/cm^3 */
  r_i = r_f = 0.005;   /* 0.5 mm only */
  dr = 0.001;
  t_i = 0.0;
  t_f = 3.6e3;
  dt = 1.0;

  n_evalues=1000000;
  
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-t"))
      sscanf(argv[++i],"%lf,%lf,%lf",&t_0,&t_f,&dt);
    else if (!strcmp(argv[i],"-r"))
      sscanf(argv[++i],"%lf,%lf,%lf",&r_i,&r_f,&dr);
    else if (!strcmp(argv[i],"-log")) log = 1;   
    else if (!strcmp(argv[i],"-ne")) n_evalues = atoi(argv[++i]);   
    else if (!strcmp(argv[i],"-ra")) ra = atof(argv[++i]);
    else if (!strcmp(argv[i],"-D")) D = atof(argv[++i]);
    else if (!strcmp(argv[i],"-R")) R = atof(argv[++i]);
    else if (!strcmp(argv[i],"-c_0")) c_0 = atof(argv[++i]);
  }
  
  fprintf(stdout,"#---------------SI UNITS----------------\n");  
  fprintf(stdout,"# Number of eigenvalues is %i\n",n_evalues);
  fprintf(stdout,"# D %.5le m^2/s, R %.5le m, "
	  "ra %.5le m, c_0 %.5le g/m^3\n#\n",
	  D,R,ra,c_0);
  
  init(R,ra);

  fprintf(stdout,"# radius(m) time(s) concentration(g/m^3)\n");

  for (t=t_0;t<=t_f;t=log?(t*dt):(t+dt)) {
    for (r=r_i;r<=r_f;r+=dr) {
      fprintf(stdout,"%.5le %.5le %.5le\n",r,t,c_analytic(r,t,ra,R,D,c_0));
    }
  }
}

