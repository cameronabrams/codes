/* 
   Metropolis Monte Carlo simulation of hard disks confined in a square

   Cameron F. Abrams

   Written for the course CHE 800-002, Molecular Simulation
   Spring 0304

   compile using "gcc -o hdisk hdisk.c -lm -lgsl"
   (assumes the GNU Scientific Library is installed)

   runs as "./hdisk -N <number_of_particles> -rho <density> \
                    -R <radius_of_circle> \
                    -nc <numcycles(1e6)>  \
		    -s <seed(?)>"

   You must have the GNU Scientific Library installed; see
   the coursenotes to learn how to do this.

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2004-2016
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

void out_fig ( FILE * fp, double * rx, double * ry, int n, double L2 ) {
  int i;
  double L = sqrt(L2);
  int iL = (int)(L*300);
  fprintf(fp,"#FIG 3.2\n" "Landscape\n" "Center\n" "Inches\n" "Letter\n"
	  "100.00\n" "Single\n" "-2\n" "1200 2\n");
  fprintf(fp,"2 2 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 5\n                    %i %i %i %i %i %i %i %i %i %i \n",0,0,0,iL,iL,iL,iL,0,0,0);
  for (i=0;i<n;i++) {
    int x = (int)((rx[i])*300);
    int y = (int)((ry[i])*300);
    fprintf(fp,"1 3 0 1 0 7 50 -1 -1 0.000 1 0.0000 %i %i "
	    "150 150 %i %i %i %i\n",x,y,x,y,x,y);
  }
}

void out ( FILE * fp, double * rx, double * ry, int n ) {
  int i;
  for (i=0;i<n;i++) {
    fprintf(fp,"%.5lf %.5lf\n",rx[i],ry[i]);
  }
}

/* Initialize particle positions by assigning them
   randomly while avoiding overlaps. */
void init ( double * rx, double * ry, int * state,
	    int n, double L2, double s2, gsl_rng * r ) {
  int i,j,reject;
  double L = sqrt(L2), r2, sx, sy;
  /* nMax is maximum number of insertion trials */
  int nMax = 10000, nTrials=0;
  
  FILE * tmp_fp;

  for (i=0;i<n;i++) {
    reject = 1;
    nTrials = 0;
    while (reject&&nTrials<nMax) {
      reject = 0;
      /* randomly locate particle in square (0,0)->(L,L) */
      rx[i] = L*gsl_rng_uniform(r);
      ry[i] = L*gsl_rng_uniform(r);
      state[i]=0; // inside the square
      /* check for overlaps; reject if an overlap is found */
      for (j=0;j<n;j++) {
	if (j!=i) {
	  sx  = rx[i]-rx[j];
	  sy  = ry[i]-ry[j];
	  r2  = sx*sx + sy*sy;
	  if (r2 < s2) reject=1;
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
  tmp_fp=fopen("0.fig","w");
  out_fig(tmp_fp,rx,ry,n,L2);
  fclose(tmp_fp);
}

int main ( int argc, char * argv[] ) {

  double * rx, * ry;
  int N=-1,c,a;
  double L=5, L2, s2=1.0;
  double rho=0.5;
  double dr=0.1,dx,dy,theta,sx,sy,r2;
  int i,j;
  int nCycles = 10, fSamp = 1;
  int nAcc, reject;
  int noob=0,novl=0;
  int short_out = 0;
  FILE * tmp_fp;
  char fn[20];
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned long int Seed = 23410981;
  int * state;

  /* Here we parse the command line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-L")) L=atof(argv[++i]);
    else if (!strcmp(argv[i],"-rho")) rho=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dr")) dr=atof(argv[++i]);
    else if (!strcmp(argv[i],"-s")) s2=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nc")) nCycles = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-fs")) fSamp = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-so")) short_out = 1;
    else if (!strcmp(argv[i],"-seed")) Seed = (unsigned long)atoi(argv[++i]);
  }

  /* precompute square of particle diameter */
  s2 = s2*s2;

  /* compute area of domain */
  L2 = L*L;

  /* If N was not set by the user, calculate it based on 
     the given values of density and radius */
  if (N==-1) {
    N = (int)(rho*L2);
  }
  /* Otherwise, recompute the radius. */
  else {
   L2 = ((double)N)/rho;
   L=sqrt(L2);
  }

  if (!short_out) 
    fprintf(stdout,"# L = %.5lf; rho = %.5lf; N = %i\n",
	  sqrt(L2),rho,N);

  /* Seed the random number generator */
  gsl_rng_set(r,Seed);

  /* allocate the position arrays */
  rx = (double*)malloc(N*sizeof(double));
  ry = (double*)malloc(N*sizeof(double));

  /* allocate state array */
  state = (int*)malloc(N*sizeof(int));

  /* generate initial positions that fit inside 
     circle with radius R and guarantee no 
     particles overlap. */
  init(rx,ry,state,N,L2,s2,r);

  nAcc = 0;
  for (c=0;c<nCycles;c++) {
    /* Make N displacement attempts */
    for (a=0;a<N;a++) {
      /* randomly select a particle */
      i=(int)gsl_rng_uniform_int(r,N);
      /* calculate displacement */
      dx = dr*(0.5-gsl_rng_uniform(r));
      dy = dr*(0.5-gsl_rng_uniform(r));
      /* displace particle */
      rx[i]+=dx;
      ry[i]+=dy;

      /* reject if outside domain */
      
      reject=rx[i]<0||ry[i]<0||rx[i]>L||ry[i]>L;
      if (reject) noob++;

      /* check for overlaps with other particles */
      if (!reject) {
	for (j=0;j<N;j++) {
	  if (j!=i) {
	    sx  = rx[i]-rx[j];
	    sy  = ry[i]-ry[j];
	    r2  = sx*sx + sy*sy;
	    if (r2 < s2) reject=1;
	  }
	}
	if (reject) novl++;
      }

      /* if move is rejected, undo displacement */
      if (reject) {
	rx[i]-=dx;
	ry[i]-=dy;
      }
      else nAcc++;
    }
    if (!(c%fSamp)) {
      sprintf(fn,"%i.fig",c);
      tmp_fp=fopen(fn,"w");
      out_fig(tmp_fp,rx,ry,N,L2);
      fclose(tmp_fp);
    }
  }
  if (short_out) {
    fprintf(stdout,"%.6lf %.6lf %.6lf %.6lf\n",
	    dr,((double)nAcc)/(N*nCycles),
	    ((double)noob)/(N*nCycles-nAcc),
	    ((double)novl)/(N*nCycles-nAcc)
	    );
  }
  else {
    fprintf(stdout,"Results:\n"
	    "Number of Trial Moves:          %i\n"
	    "Maximum Displacement Length:    %.5lf\n"
	    "Acceptance Ratio:               %.5lf\n"
	    "Reject Fraction Out-of-bounds:  %.5lf\n"
	    "Reject Fraction Overlap:        %.5lf\n",
	    N*nCycles,dr,((double)nAcc)/(N*nCycles),
	    ((double)noob)/(N*nCycles-nAcc),
	    ((double)novl)/(N*nCycles-nAcc)
	    );
  }
}
