/* 
   Metropolis Monte Carlo simulation of a Lennard-Jones fluid
   in a periodic boundary

   Cameron F. Abrams

   Written for the course CHE 800-002, Molecular Simulation
   Spring 0304

   compile using "gcc -o mclj mclj.c -lm -lgsl"
   (assumes the GNU Scientific Library is installed)

   runs as "./mclj -N <number_of_particles> -rho <density> \
                    -nc <numcycles(1e6)>  -dr <delta-r> \
		    -s <seed(?)> -ne <#equil.cycles(100)>"

   You must have the GNU Scientific Library installed; see
   the coursenotes to learn how to do this.

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2004
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>


/* An N^2 algorithm for computing the total energy.  The virial
   is also computed and returned in *vir. */
double total_e ( double * rx, double * ry, double * rz, int N, double L,
		 double rc2, int tailcorr, double ecor, 
		 int shift, double ecut, double * vir ) {
   int i,j;
   double dx, dy, dz, r2, r6i;
   double e = 0.0, hL=L/2.0;

   *vir=0.0;
   for (i=0;i<(N-1);i++) {
     for (j=i+1;j<N;j++) {
	dx  = (rx[i]-rx[j]);
	dy  = (ry[i]-ry[j]);
	dz  = (rz[i]-rz[j]);
	/* Periodic boundary conditions: Apply the minimum image
	   convention; note that this is *not* used to truncate the
	   potential as long as there an explicit cutoff. */
	if (dx>hL)       dx-=L;
	else if (dx<-hL) dx+=L;
	if (dy>hL)       dy-=L;
	else if (dy<-hL) dy+=L;
	if (dz>hL)       dz-=L;
	else if (dz<-hL) dz+=L;
	r2 = dx*dx + dy*dy + dz*dz;
	if (r2<rc2) {
	  r6i   = 1.0/(r2*r2*r2);
	  e    += 4*(r6i*r6i - r6i) - (shift?ecut:0.0);
	  *vir += 48*(r6i*r6i-0.5*r6i);
	}
     }
   }
   return e+(tailcorr?(N*ecor):0.0);
}

/* Initialize particle positions by assigning them
   on a cubic grid, then scaling positions 
   to achieve a given box size and thereby, volume,
   and density */
void init ( double * rx, double * ry, double * rz,
	    int n, double L, gsl_rng * r ) {
  int i,ix,iy,iz;
  
  int n3=2;
  /* Find the lowest perfect cube, n3, greater than or equal to the
     number of particles */
  while ((n3*n3*n3)<n) n3++;

  ix=iy=iz=0;
  /* Assign particle positions */
  for (i=0;i<n;i++) {
    rx[i] = ((double)ix+0.5)*L/n3;
    ry[i] = ((double)iy+0.5)*L/n3;
    rz[i] = ((double)iz+0.5)*L/n3;
    ix++;
    if (ix==n3) {
      ix=0;
      iy++;
      if (iy==n3) {
	iy=0;
	iz++;
      }
    }
  }
}

int main ( int argc, char * argv[] ) {

  double * rx, * ry, * rz;
  int N=216,c,a;
  double L=0.0;
  double rho=0.5, T=1.0, rc2 = 1.e20, vir, vir_old, vir_sum, pcor, V;
  double E_new, E_old, esum, rr3, ecor, ecut;
  double dr=0.1,dx,dy,dz;
  double rxold,ryold,rzold;
  int i,j;
  int nCycles = 10, nSamp, nEq=1000;
  int nAcc;
  int short_out=0;
  int shift=0;
  int tailcorr=1;

  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned long int Seed = 23410981;

  /* Here we parse the command line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-rho")) rho=atof(argv[++i]);
    else if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dr")) dr=atof(argv[++i]);
    else if (!strcmp(argv[i],"-rc")) rc2=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nc")) nCycles = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-ne")) nEq = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-so")) short_out=1;
    else if (!strcmp(argv[i],"+tc")) tailcorr=0;
    else if (!strcmp(argv[i],"-sh")) shift=1;
    else if (!strcmp(argv[i],"-seed")) 
      Seed = (unsigned long)atoi(argv[++i]);
    else {
      fprintf(stderr,"Error.  Argument '%s' is not recognized.\n",
		 argv[i]);
      exit(-1);
    }
  }

  /* Compute the side-length */
  L = pow((V=N/rho),0.3333333);

  /* Compute the tail-corrections; assumes sigma and epsilon are both 1 */
  rr3 = 1.0/(rc2*rc2*rc2);
  ecor = 8*M_PI*rho*(rr3*rr3*rr3/9.0-rr3/3.0);
  pcor = 16.0/3.0*M_PI*rho*rho*(2./3.*rr3*rr3*rr3-rr3);
  ecut = 4*(rr3*rr3*rr3*rr3-rr3*rr3);

  /* Compute the *squared* cutoff, reusing the variable rc2 */
  rc2*=rc2;

  /* For computational efficiency, use reciprocal T */
  T = 1.0/T;

  /* compute box volume */
  V = L*L*L;
  
  /* Output some initial information */
  fprintf(stdout,"# NVT MC Simulation of a Lennard-Jones fluid\n");
  fprintf(stdout,"# L = %.5lf; rho = %.5lf; N = %i; rc = %.5lf\n",
	  L,rho,N,sqrt(rc2));
  fprintf(stdout,"# nCycles %i, nEq %i, seed %d, dR %.5lf\n",
	  nCycles,nEq,Seed,dr);
  
  /* Total number of cycles is number of "equilibration" cycles plus
     number of "production" cycles */
  nCycles+=nEq;

  /* Seed the random number generator */
  gsl_rng_set(r,Seed);

  /* Allocate the position arrays */
  rx = (double*)malloc(N*sizeof(double));
  ry = (double*)malloc(N*sizeof(double));
  rz = (double*)malloc(N*sizeof(double));

  /* Generate initial positions on a cubic grid, 
     and measure initial energy */
  init(rx,ry,rz,N,L,r);
  E_old = total_e(rx,ry,rz,N,L,rc2,tailcorr,ecor,shift,ecut,&vir_old);

  nAcc = 0;
  esum = 0.0;
  nSamp = 0;
  vir_sum = 0.0;
  for (c=0;c<nCycles;c++) {

    /* Randomly select a particle */
    i=(int)gsl_rng_uniform_int(r,N);
    /* calculate displacement */
    dx = dr*(0.5-gsl_rng_uniform(r));
    dy = dr*(0.5-gsl_rng_uniform(r));
    dz = dr*(0.5-gsl_rng_uniform(r));

    /* Save the current position of particle i */
    rxold=rx[i];
    ryold=ry[i];
    rzold=rz[i];

    /* Displace particle i */
    rx[i]+=dx;
    ry[i]+=dy;
    rz[i]+=dz;

    /* Apply periodic boundary conditions */
    if (rx[i]<0.0) rx[i]+=L;
    if (rx[i]>L)   rx[i]-=L;
    if (ry[i]<0.0) ry[i]+=L;
    if (ry[i]>L)   ry[i]-=L;
    if (rz[i]<0.0) rz[i]+=L;
    if (rz[i]>L)   rz[i]-=L;

    /* Get the new energy */
    E_new = total_e(rx,ry,rz,N,L,rc2,tailcorr,ecor,shift,ecut,&vir);
      
    /* Conditionally accept... */
    if (gsl_rng_uniform(r) < exp(-T*(E_new-E_old))) {
      E_old=E_new;
      vir_old=vir;
      nAcc++;
    }
    /* ... or reject the move; reassign the old positions */
    else {
      rx[i]=rxold;
      ry[i]=ryold;
      rz[i]=rzold;
    }

    /* Sample: default frequency is once per trial move; We must
       include results of a move regardless of whether the move is
       accepted or rejected. */
    if (c>nEq) {
      esum+=E_old;
      vir_sum+=vir_old;
      nSamp++;
    }
  }

  /* Output delta-r, the acceptance ratio, 
     and the average energy/particle */
  if (short_out)
    fprintf(stdout,"%.6lf %.5lf %.5lf %.5lf\n",
	    dr,((double)nAcc)/(N*nCycles),
	    esum/nSamp/N,vir_sum/3.0/nSamp/V+rho*T+pcor);
  else
    fprintf(stdout,"NVT Metropolis Monte Carlo Simulation"
	    " of the Lennard-Jones fluid.\n"
	    "---------------------------------------------\n"
	    "Number of particles:              %i\n"
	    "Number of cycles:                 %i\n"
	    "Cutoff radius:                    %.5lf\n"
	    "Maximum displacement:             %.5lf\n"
	    "Density:                          %.5lf\n"
	    "Temperature:                      %.5lf\n"
	    "Tail corrections applied?         %s\n"
	    "Shifted potential?                %s\n"
	    "Results:\n"
	    "Potential energy tail correction: %.5lf\n"
	    "Pressure tail correction:         %.5lf\n"
	    "Potential energy shift at cutoff  %.5lf\n"
	    "Acceptance ratio:                 %.5lf\n"
	    "Energy/particle:                  %.5lf\n"
	    "Ideal gas pressure:               %.5lf\n"
	    "Virial:                           %.5lf\n"
	    "Total Pressure:                   %.5lf\n"
	    "Program ends.\n",
	    N,nCycles,sqrt(rc2),dr,rho,1.0/T,
	    tailcorr?"Yes":"No",shift?"Yes":"No",
	    ecor,pcor,ecut,
	    ((double)nAcc)/(N*nCycles),
	    esum/nSamp/N,
	    rho/T,vir_sum/3.0/nSamp/V,
	    vir_sum/3.0/nSamp/V+rho/T+(tailcorr?pcor:0.0));
}
