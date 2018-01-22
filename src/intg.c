#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
/*
 * intg.c -- uses the trapezoidal rule to integrate a function
 * supplied in the discrete form x_i, f(x_i).  Computes the maximum
 * value of f, the integral of f, `thickness' of f (integral/max),
 * and the average value of x, <x>, over the interval.
 * 
 * The data is assumed to be cartesian, such that the differential
 * element is simply d(x).  The user can specify that the data is
 * spherically-symmetric polar, such that the differential element
 * is 4*pi*r^2*d(r), by using the command-line switch `-p'.
 *
 * Alternatively, the user can specify that the x-data is the 
 * polar angle theta, in which case the volume element is 
 * sin(theta)d(theta), by using the command-line switch `-pad' (if
 * the xdata is angle in degrees) or `-par' (if in RADIANS).
 *
 * Two columns of input data can be read from stdin (no arguments)
 * or from an input file (single argument).  There is no limit to the 
 * number of data points that can be considered.  Output is to stdout, and
 * can be either `short' (default) or `long' format.  Short format is
 * a single number, the integral of f(x).  Long output includes all computed
 * quantities.  Long output is specified by the command-line switch `-l'.
 *
 * Integration limits can be chosen with the switch `-limits x0,x1' where
 * x0 is the lower limit and x1 is the upper limit.  If no limits are given,
 * the integration is performed over the entire domain.
 * 
 * 
 * (c) 2000-2004 cameron abrams
 */

void enq ( double A[], double x, int n )
{
  int i;
  for (i=0;i<n-1;i++) A[i]=A[i+1];
  A[i]=x;
}

int nData;
char ln[255];
int main (int argc, char * argv[])
{
    FILE *fp=NULL;
    char *fn=NULL;
    char sfmt[255];
    char * nfmt="%.5lf";
    int i,d=2;
    double s0, s1, s2, s3, fx, max=1.e-9, x, y;
    double X[3]={0,0,0},Y[3]={0,0,0};
    double x0,x1;  // limits
    short polar=0, pangle=0, longout=0, limits=0;
    double fac=4*M_PI;
    double xspan=0.0;
    
    for (i=1;i<argc;i++) {
      if (argv[i][0]!='-') fn=argv[i];
      else if (!strcmp(argv[i],"-p")) polar=1;
      else if (!strcmp(argv[i],"-pad")) pangle=1;
      else if (!strcmp(argv[i],"-par")) pangle=2;
      else if (!strcmp(argv[i],"-l")) longout=1;
      else if (!strcmp(argv[i],"-nf")) nfmt=argv[++i];
      else if (!strcmp(argv[i],"-limits")) {
	limits=1;
	sscanf(argv[++i],"%lf,%lf",&x0,&x1);
      }
    }

    if (polar&&pangle) {
      fprintf(stderr,"Error -- specify *either* -p or -pa, but not both.\n");
      exit(-1);
    }
    
    if (!fn) fp=stdin;
    else fp=fopen(fn, "r");
    
    i=0;
    s0=s1=s2=s3=0.0;
    while (fgets(ln, 255, fp)) {
      if (ln[0]!='#') {
	sscanf(ln, "%le %le", &x, &y);
	if (!limits||(limits&&(x<=x1&&x>=x0))) {
	  //printf("data %.5le %.5le\n", x,y);
	  if (!i) xspan = x; 
	  enq(X,x,d);
	  enq(Y,y,d);
	  if (i>=(d-1)) {
	    /* trapezoidal rule kernel */
	    fx=(X[1]-X[0])*(Y[1]+Y[0]);
	    if (polar) fx*=fac*0.25*(X[0]*X[0]+X[1]*X[1]+2*X[0]*X[1]);
	    else if (pangle==1) fx*=180.0/M_PI*sin(M_PI/180.*0.5*(X[0]+X[1]));
	    else if (pangle==2) fx*=sin(0.5*(X[0]+X[1]));
	    s0+=fx;
	    s1+=0.5*(X[0]+X[1])*fx;
	    s2+=0.5*(X[0]*X[0]+X[1]*X[1])*fx;
	    s3+=(X[1]-X[0])*(Y[1]*Y[1]+Y[0]*Y[0]); // y -> y^2
	  }
	  max=y>max?y:max;
	  i++;
	}
      }
    }
    xspan = x-xspan;
    nData=i;
    if (!nData) {
      fprintf(stderr,"ERROR -- no data\n");
      exit(-1);
    }
    s0*=0.5;
    s3*=0.5;
    if (s0) s1*=0.5/s0;
    if (s0) s2*=0.5/s0;

    if (longout) sprintf(sfmt,"# nData %%i I %s I/max %s <x> %s <x^2> %s <y> %s <(Delta-y)^2>^(1/2) %s\n",
			 nfmt,nfmt,nfmt,nfmt,nfmt,nfmt);
    else sprintf(sfmt,"%s\n",nfmt);

    if (longout) printf(sfmt,nData, s0, s0/max, s1, s2, 
			s0/xspan, sqrt(s3/xspan-(s0*s0/xspan/xspan)));
    else printf(sfmt,s0);

    if (fp!=stdin) fclose(fp);

}
