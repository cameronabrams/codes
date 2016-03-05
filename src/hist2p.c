#include <stdio.h>
#include <math.h>
/*
 * hist2p.c:  
 * 
 * operates on an input histogram supplied in the discrete form 
 * x_i, H(x_i).
 * 
 * Computes the probability density function, P(x), given information
 * about the volume element d(x), such that
 * 
 * 1 = sum_i[H(x_i)] = integral[P(x)d(x)]
 *
 * Optionally, (a) computes the `boltzmann inversion' of this probability
 * density as
 *
 * U(x) ~ -kT * ln [P(x)],
 *
 * and, (b) computes the relevant force factor by 2nd-order finite difference:
 *
 * f(x_i) = -d/dx[U(x)]|x_i
 *        = [U(x(i+1)) - U(x(i-1))]/(x_(i+1) - x(i-1))
 *
 * Note on units:  the natural units for angles are radians.  
 * The input histogram may be a function of angle in degrees,
 * as therefore can the discretization.  The differential element
 * d(x) when x is an angle is sin(x)dx, where dx MUST be in radians.
 * Hence, when the finite difference (x_(i+1) - x(i-1)) is used for
 * dx, it is converted to radians.
 * Then, the units of f(x_i) are always energy per RADIAN.
 *
 * Two columns of input data can be read from stdin (no arguments)
 * or from an input file (single argument).  Output is to stdout.
 * MAXPTS can be redefined (below) if more than 5000 points are 
 * being considered.
 * 
 * (c) 2001, 2002 cameron abrams
 * max-planck-institute for polymer research
 * mainz, germany
 */

#define MAXPTS 5000
double x[MAXPTS];
double H[MAXPTS];
double Hsm[MAXPTS];
double P[MAXPTS];
double U[MAXPTS];
double Usm[MAXPTS];
double F[MAXPTS];

/* setf -- allows user to set force to an explicit value over
   an explicit interval.  */
#define MAXSETF 10

#define sgn(a) ((a)<0.0?-1.0:1.0)

double v_rad ( double x ) {
  return x*x;
}
double v_pol ( double x ) {
  return sin(x);
}
double v_azi ( double x ) {
  return 1.0;
}

enum aunits {RADIANS,DEGREES};
enum vtype {RAD,POL,AZI};
char * Velem[3] = {"rad","pol","azi"};
double (*v_elem[3])( double );

double Umax=20.0;
double Fmax=-1;

int nData;
char ln[255];
int main (int argc, char * argv[])
{
    FILE *fp=NULL;
    char *fn=NULL;
    int i, vtyp=AZI,au=RADIANS, nslav, lo, hi, j, nsetf=0;
    double diff, tot, slav, *u, *h;
    double deg2rad=M_PI/180.0;
    double setfmin[MAXSETF],setfmax[MAXSETF],setf[MAXSETF];
    char *fmt=NULL;
    char prfmt[50];
    double (*ve)(double);
    int Usmooth=0, Hsmooth=0, window=10;

    v_elem[RAD]=&v_rad;
    v_elem[POL]=&v_pol;
    v_elem[AZI]=&v_azi;

    for (i=1;i<argc;i++) {
      if (argv[i][0]!='-') fn=argv[i];
      else if (!strcmp(argv[i],"-fmt")) fmt=argv[++i];
      else if (!strcmp(argv[i],"-umax")) Umax=atof(argv[++i]);
      else if (!strcmp(argv[i],"-fmax")) Fmax=atof(argv[++i]);
      else if (!strcmp(argv[i],"-usm")) Usmooth=1;
      else if (!strcmp(argv[i],"-hsm")) Hsmooth=1;
      else if (!strcmp(argv[i],"-win")) window=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-v")) {
	int j;
	i++;
	for (j=0;j<3&&strcmp(argv[i],Velem[j]);j++);
	if (j==3) {
	  fprintf(stderr,"Error: vol elem typ `%s' not recognized.\n",
		  argv[i]);
	  fprintf(stderr,"Please specify one of ");
	  for (j=0;j<3;j++) 
	    fprintf(stderr," %s%s",Velem[j],j<2?",":".\n");
	  exit(-1);
	}
	vtyp=j;
	ve=v_elem[j];
      }
      else if (!strcmp(argv[i],"-au")) {
	i++;
	if (!strcmp(argv[i],"deg")) au=DEGREES;
	else if (!strcmp(argv[i],"rad")) au=RADIANS;
	else {
	  fprintf(stderr,"Error: angle units type `%s' not recognized.\n",
		  argv[i]);
	  fprintf(stderr,"Please specify `deg' or `rad'.\n");
	  exit(-1);
	}
      }
      else if (!strcmp(argv[i],"-setf")) {
	sscanf(argv[++i],"%lf,%lf,%lf",
	       &setfmin[nsetf],&setfmax[nsetf],&setf[nsetf]);
	nsetf++;
      }
    }

    if (!fmt) sprintf(prfmt,"%.5lf %.5lf\n");
    else sprintf(prfmt,"%s %s\n", fmt, fmt);

    if (!fn) fp=stdin;
    else fp=fopen(fn, "r");
    
    i=0;
    tot=0.0;
    while (fgets(ln, 255, fp)) {
      if (ln[0]!='#') {
	sscanf(ln, "%lf %lf", &x[i], &H[i]);
	i++;
	if (i==MAXPTS) {
	  fprintf(stderr,"Error: MAXPTS(%i) too small.\n",MAXPTS);
	  exit(-1);
	}
      }
    }
    nData=i;
    //    fprintf(stderr,"# ndata %i\n",nData);fflush(stderr);
    if (Hsmooth) { // smooth the histogram if requested
      for (i=0;i<nData;i++) {
	if (i>window/2&&i<(nData-window/2)) {
	  lo=i-window/2;
	  hi=i+window/2;
	}
	else if (i<=window/2) {
	  lo=0;
	  hi=i+window/2;
	}
	else if (i>=(nData-window/2)) {
	  lo=i-window/2;
	  hi=nData-1;
	}
	nslav=0;
	slav=0.0;
	for (j=lo;j<=hi;j++) {
	  slav+=H[j];
	  nslav++;
	}
	Hsm[i]=slav/nslav;
      }
    }

    h=(Hsmooth?Hsm:H);

    /* renormalize if the input is not normalized already or
       if it was smoothed */
    for (i=0;i<nData;i++) tot+=h[i];
    if (tot!=1.0) {
      fprintf(stdout,"# renormalizing by %.15f\n",tot);
      for (i=0;i<nData;i++) h[i]/=tot;
      tot=0.0;
      for (i=0;i<nData;i++) tot+=h[i];
      fprintf(stdout,"# now sums to %.15f\n",tot);
    }

    /* compute the probability density distribution P(x), accounting
     * for the proper volume element for this type of coordinate */
    for (i=0;i<nData;i++) {
      tot=(*ve)((au==DEGREES?deg2rad:1.0)*x[i]);
      P[i]=tot?h[i]/tot:0.0;
      U[i]=P[i]>0.0?-log(P[i]):1.0e99;
      U[i]=U[i]>Umax?Umax:U[i];
      Usm[i]=0.0; // just an initialization
      //      fprintf(stderr,"%.10f %.10f %.10f %.10f\n",
      //           x[i],H[i],P[i],U[i]>Umax?Umax:U[i]);fflush(stderr);
    }
    
    /* smooth the potential */
    if (Usmooth) {
      for (i=0;i<nData;i++) {
	if (i>window/2&&i<nData-window/2) {
	  lo=i-window/2;
	  hi=i+window/2;
	}
	else if (i<=window/2) {
	  lo=0;
	  hi=i+window/2;
	}
	else if (i>=nData-window/2) {
	  lo=i-window/2;
	  hi=nData-1;
	}
	nslav=0;
	slav=0.0;
	for (j=lo;j<=hi;j++) {
	  slav+=U[j];
	  nslav++;
	}
	Usm[i]=slav/nslav;
      }
      //      for (i=0;i<nData;i++) {
      //	U[i]=Usm[i];
      //      }
    }

    /* Differentiate and negate potential, and output; note that
     * angles are converted to natural units of radians! */
    printf("# x   force   Norm.histogram   Boltz.potential\n");
    diff=0.0;
    u=(Usmooth?Usm:U);
    for (i=0;i<nData;i++)
    {
      if (i==0) 
	diff=(u[i+1]-u[i])/((x[i+1]-x[i])*(au==DEGREES?deg2rad:1.0));
      else if (i==(nData-1)) 
	diff=(u[i]-u[i-1])/((x[i]-x[i-1])*(au==DEGREES?deg2rad:1.0));
      else 
	diff=(u[i+1]-u[i-1])/((x[i+1]-x[i-1])*(au==DEGREES?deg2rad:1.0));
      F[i]=-diff;
      if (Fmax!=-1&&fabs(F[i])>Fmax) F[i]=sgn(F[i])*Fmax;
      if (nsetf) {
	for (j=0;j<nsetf;j++) {
	  if (x[i]>=setfmin[j]&&x[i]<=setfmax[j])
	    F[i]=setf[j];
	}
      }
      // printf("%.10f %.10f %.10f %.10f %.10f\n", x[i], F[i], P[i], U[i], Usm[i]);
      printf("%.10lf %.13lf %.13lf %.13lf\n", x[i], F[i], h[i], u[i]);
      fflush(stdout);
    }
    
/*     printf("# nData = %i\n", nData); */
    
    if (fp!=stdin) fclose(fp);

}
