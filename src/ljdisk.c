/* 
   Metropolis Monte Carlo simulation of Lennard-Jones disks in a square periodic domain

   Cameron F. Abrams

   Written for the fun

   compile using "gcc -o ljdisk ljdisk.c -lm -lgsl"
   (assumes the GNU Scientific Library is installed)

   runs as "./ljdisk -N <number_of_particles> -rho <density> \
                     -nc <numcycles(1e6)>  \
		     -s <seed(?)> -maxdr maximum displacement length"

   You must have the GNU Scientific Library installed; see
   the coursenotes to learn how to do this.

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2019
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>

/* minimum image convention */
double MIC ( double d, double h ) {
  if (d<-h) d+=2*h;
  else if (d>h) d-=2*h;
  return d;
}

double LJ_E ( double * rx, double * ry, int n, double bhs, int * h, int nb, double hrmax ) {
  int i, j, b;
  double e = 0.0, dx, dy, dr2, dr_6, dr_12;
  for (i=0;i<n;i++) {
     for (j=i+1;j<n;j++) {
	 dx=MIC(rx[i]-rx[j],bhs);
	 dy=MIC(ry[i]-ry[j],bhs);
         dr2=dx*dx+dy*dy;
	 if (h) {
            b=(int)(nb*sqrt(dr2)/hrmax);
	    if (b<nb) h[b]+=2;
	 }
	 dr_6=pow(dr2,-3);
	 dr_12=dr_6*dr_6;
	 e+=dr_12-dr_6;
     }
  }
  return 4*e;
}

double LJ_deltaE ( double * rx, double * ry, int N, int i, double idx, double idy, double bhs ) {
   int j;
   double e = 0.0, dx, dy, dr2, dr_6, dr_12;
   double newx=MIC(rx[i]+idx,bhs);
   double newy=MIC(ry[i]+idy,bhs);
   for (j=0;j<N;j++) {
     if (j!=i) {
        dx=MIC(newx-rx[j],bhs);
	dy=MIC(newy-ry[j],bhs);
	dr2=dx*dx+dy*dy;
	dr_6=pow(dr2,-3);
	dr_12=dr_6*dr_6;
	e+=dr_12-dr_6;
        dx=MIC(rx[i]-rx[j],bhs);
	dy=MIC(ry[i]-ry[j],bhs);
	dr2=dx*dx+dy*dy;
	dr_6=pow(dr2,-3);
	dr_12=dr_6*dr_6;
	e-=dr_12-dr_6;
     }
   }
   return 4*e;
}

void out ( FILE * fp, double * rx, double * ry, int n ) {
  int i;
  for (i=0;i<n;i++) {
    fprintf(fp,"%.5lf %.5lf\n",rx[i],ry[i]);
  }
}

/* Initialize particle positions by assigning them
   randomly while avoiding overlaps. */
void init ( double * rx, double * ry, 
	    int n, double bhs, gsl_rng * r ) {
  int i,j,reject;
  double L = 2*bhs, r2, sx, sy;
  /* nMax is maximum number of insertion trials */
  int nMax = 10000, nTrials=0;
  
  FILE * tmp_fp;

  for (i=0;i<n;i++) {
    reject = 1;
    nTrials = 0;
    while (reject&&nTrials<nMax) {
      reject = 0;
      /* randomly locate particle i in an LxL domain with (0,0) dead-center */
      rx[i] = bhs*(1.0-2*gsl_rng_uniform(r));
      ry[i] = bhs*(1.0-2*gsl_rng_uniform(r));
      /* check for overlaps; reject if an overlap is found */
      for (j=0;j<n;j++) {
	if (j!=i) {
	  sx = rx[i]-rx[j];
	  sy = ry[i]-ry[j];
	  sx = MIC(sx,bhs);
	  sy = MIC(sy,bhs);
	  r2 = sx*sx+sy*sy;
	  if (r2 < 1.0) reject=1;
	}
      }
      nTrials++;
    }
    if (nTrials==nMax) {
      fprintf(stderr,"# Error -- could not initialize position "
	      "of particle %i in %i trials\n",i,nTrials);
      fprintf(stderr,"# program ends.\n");
      exit(-1);
    }
  }
  tmp_fp=fopen("ldtmpcfg_init.xy","w");
  out(tmp_fp,rx,ry,n);
  fclose(tmp_fp);
}

int main ( int argc, char * argv[] ) {

  double * rx, * ry;
  int N=-1;
  double L=-1.0, bhs;
  double rho=-1.0;
  double T = 1.0;
  double dr=0.5,dx,dy;
  double arg,test_x,bf,EE,E0,dE,phi,this_dr;
  int i,j,a,c;
  int nCycles = 10;
  int nAcc, accept;
  int short_out = 0;
  FILE * tmp_fp;

  int *h;
  int nbins=100;
  double hrmax=bhs;

  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned long int Seed = 23410981;

  /* Here we parse the command line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-L")) L=atof(argv[++i]);
    else if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i],"-rho")) rho=atof(argv[++i]);
    else if (!strcmp(argv[i],"-maxdr")) dr=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nc")) nCycles = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-so")) short_out = 1;
    else if (!strcmp(argv[i],"-seed")) Seed = (unsigned long)atoi(argv[++i]);
  }

  /* If N was not set by the user, calculate it based on 
     the given values of density and box length */
  if (N==-1) {
    N = (int)(rho*L*L);
  } else if (rho==-1) {
  }
  /* Otherwise, recompute the box length. */
  else {
    L = sqrt(((double)N)/rho);
  }

  bhs = L/2;
  hrmax = bhs;

  if (!short_out) 
    fprintf(stdout,"# L = %.5lf; rho = %.5lf; N = %i\n",
	  L,rho,N);

  /* Seed the random number generator */
  gsl_rng_set(r,Seed);

  /* allocate the position arrays */
  rx = (double*)malloc(N*sizeof(double));
  ry = (double*)malloc(N*sizeof(double));
  h = (int*)malloc(nbins*sizeof(int));
  
  /* generate initial positions that fit inside 
     periodic box with side-length and guarantee no 
     particles overlap. */
  init(rx,ry,N,bhs,r);
  for (i=0;i<nbins;i++) h[i]=0;

  EE=E0=LJ_E(rx,ry,N,bhs,NULL,0,0);
  fprintf(stdout,"# initial energy %.5lf\n",EE);

  nAcc = 0;
  for (c=0;c<nCycles;c++) {
    /* Make N displacement attempts */
    for (a=0;a<N;a++) {
      /* randomly select a particle */
      i=(int)gsl_rng_uniform_int(r,N);
      /* calculate displacement */
      phi=2*M_PI*gsl_rng_uniform(r);
      this_dr=dr*gsl_rng_uniform(r);
      dx = this_dr*cos(phi);
      dy = this_dr*sin(phi);
      /* displace particle */
      dE = LJ_deltaE(rx,ry,N,i,dx,dy,bhs);
      accept=0;
      arg=-dE/T;
      if (arg>0) accept=1;
      else {
	 bf=exp(arg);
	 test_x=gsl_rng_uniform(r);
	 if (test_x<bf) accept=1;
      }
      if (accept) {
         EE+=dE;
	 rx[i]+=dx;
	 ry[i]+=dy;
         nAcc++;
      }
    }
    fprintf(stdout,"# cycle %d energy check: %.5lf tally %.5lf\n",c,LJ_E(rx,ry,N,bhs,h,nbins,hrmax),EE);
  }
  tmp_fp=fopen("hdtmpcfg","w");
  out(tmp_fp,rx,ry,N);
  fclose(tmp_fp);
  tmp_fp=fopen("rdf.dat","w");
  for (i=0;i<nbins;i++) 
      fprintf(tmp_fp,"%.5f %.5f\n",i*hrmax/nbins,h[i]/(N*nCycles*M_PI*((i+1)*(i+1)-i*i)*(hrmax/nbins)*(hrmax/nbins)*rho));
  fclose(tmp_fp);
  fprintf(stdout,"Results:\n"
	    "Number of Trial Moves:          %i\n"
	    "Maximum Displacement Length:    %.5lf\n"
	    "Acceptance Ratio:               %.5lf\n",
	    N*nCycles,dr,((double)nAcc)/(N*nCycles));
}
