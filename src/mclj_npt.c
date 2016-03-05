/* 
   Metropolis Monte Carlo simulation of a Lennard-Jones fluid
   in a periodic boundary:  NPT implementation

   Cameron F. Abrams

   Written for the course CHE 800-002, Molecular Simulation
   Spring 0304

   compile using "gcc -o mclj_npt mclj_npt.c -lm -lgsl"
   (assumes the GNU Scientific Library is installed)

   runs as "./mclj_npt -N <number_of_particles> -rho <density> \
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
#include <string.h>
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

int mcmove ( double * rx, double * ry, double * rz, int i, int N,
	     int tailcorr, double ecor, int shift, double ecut, 
	     double L, double rc2, double dr, 
	     double T, double * E_old, double * vir_old,
	     gsl_rng * r ) {
  double dx, dy, dz;
  double rxold, ryold, rzold, E_new, vir;

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

  E_new = total_e(rx,ry,rz,N,L,rc2,tailcorr,ecor,shift,ecut,&vir);

  /* Conditionally accept... */
  if (gsl_rng_uniform(r) < exp(-T*(E_new-(*E_old)))) {
    *E_old=E_new;
    *vir_old=vir;
    return 1;
  }
  /* ... or reject the move; reassign the old positions */
  else {
    rx[i]=rxold;
    ry[i]=ryold;
    rz[i]=rzold;
    return 0;
  }
}

int mcvol ( double * rx, double * ry, double * rz, int N,
	    double * L, int tailcorr, double * ecor, double * pcor,
	    int shift, double * ecut, double * rc2,
	    double dlnV, double T, double P, 
	    double * E_old, double * vir_old, gsl_rng * r ) {

  int i;
  double V0 = (*L)*(*L)*(*L);
  double f, lnV, lnV0, V, L1, E_new, vir, rho, rr3, arg, rc;


  lnV0 = log(V0);
  lnV = lnV0 + dlnV*(gsl_rng_uniform(r)-0.5);
  V=exp(lnV);
  L1 = pow(V,0.333333333);
  f = L1/(*L);
  for (i=0;i<N;i++) {
    rx[i]*=f; // rx[i] = rx[i] * f;
    ry[i]*=f;
    rz[i]*=f;
  }
  (*L) = L1;
  /* compute the new cutoff and tail corrections, if necessary */
  (*rc2) *= f*f;
  rc=sqrt(*rc2);
  rho = N/V;
  rr3 = 1.0/(rc*rc*rc);
  *ecor = 8*M_PI*rho*(rr3*rr3*rr3/9.0-rr3/3.0);
  *pcor = 16.0/3.0*M_PI*rho*rho*(2./3.*rr3*rr3*rr3-rr3);
  *ecut = 4*(rr3*rr3*rr3*rr3-rr3*rr3);
  
  E_new = total_e(rx,ry,rz,N,*L,*rc2,tailcorr,*ecor,shift,*ecut,&vir);

  /* Conditionally accept... */
  arg = -T*(P*(V-V0)-(N+1)*(lnV-lnV0)/T+E_new-(*E_old));
  if (gsl_rng_uniform(r) < exp(arg)) {
    *E_old=E_new;
    *vir_old=vir;
    return 1;
  }
  /* ... or reject the move; reassign the old positions */
  else {
    f = 1.0/f;
    for (i=0;i<N;i++) {
      rx[i]*=f;
      ry[i]*=f;
      rz[i]*=f;
    }
    (*L)*=f;
    /* recompute the old cutoff and tail corrections, if necessary */
    (*rc2) *= f*f;
    rho = N/V0;
    rc = sqrt(*rc2);
    rr3 = 1.0/(rc*rc*rc);
    *ecor = 8*M_PI*rho*(rr3*rr3*rr3/9.0-rr3/3.0);
    *pcor = 16.0/3.0*M_PI*rho*rho*(2./3.*rr3*rr3*rr3-rr3);
    *ecut = 4*(rr3*rr3*rr3*rr3-rr3*rr3);
    return 0;
  }
  
}

int main ( int argc, char * argv[] ) {

  double * rx, * ry, * rz;
  int N = 216,c;
  double L = 0.0;
  double rho = 0.5, T = 1.0, P = 1.0, beta = 1.0;
  double rc = 2.5, rc2 = rc*rc, vir, p_sum, rho_sum, pcor, V;
  double E, esum, rr3, ecor, ecut;
  double dr = 0.1, dlnV = 0.1;
  int i;
  int nCycles = 10, nSamp, nEq = 1000;
  int acc, nAcc, nTransAcc, nVolAcc, nTransAtt, nVolAtt;
  int short_out = 0;
  int tailcorr = 1, shift = 0;

  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned long int Seed = 23410981;

  /* Here we parse the command line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-rho")) rho=atof(argv[++i]);
    else if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i],"-P")) P=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dr")) dr=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dlnV")) dlnV=atof(argv[++i]);
    else if (!strcmp(argv[i],"-rc")) rc=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nc")) nCycles = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-ne")) nEq = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-so")) short_out=1;
    else if (!strcmp(argv[i],"+tc")) tailcorr=0;
    else if (!strcmp(argv[i],"-sh")) shift=1;
    else if (!strcmp(argv[i],"-seed")) 
      Seed = (unsigned long)atoi(argv[++i]);
    else {
      fprintf(stderr,"Error.  '%s' not recognized.\n",argv[i]);
      exit(-1);
    }
  }

  /* Compute the *squared* cutoff */
  rc2 = rc*rc;

  /* For computational efficiency, use reciprocal T */
  beta = 1.0/T;

  /* Compute the side-length */
  L = pow(N/rho,0.3333333);

  /* Compute the tail-corrections; assumes sigma and epsilon are both 1 */
  rr3 = 1.0/(rc*rc*rc);
  ecor = 8*M_PI*rho*(rr3*rr3*rr3/9.0-rr3/3.0);
  pcor = 16.0/3.0*M_PI*rho*rho*(2./3.*rr3*rr3*rr3-rr3);
  ecut = 4*(rr3*rr3*rr3*rr3-rr3*rr3);

  /* Output some initial information */
  fprintf(stdout,"# NPT MC Simulation of a Lennard-Jones fluid\n");
  fprintf(stdout,"# P = %.5lf, T = %.5lf, L0 = %.5lf;"
	  " rho0 = %.5lf; N = %i; rc = %.5lf\n",
	  P,T,L,rho,N,rc);
  fprintf(stdout,"# nCycles %i, nEq %i, seed %lu, dR %.5lf, dlnV %.5lf\n",
	  nCycles,nEq,Seed,dr,dlnV);
  
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
  E = total_e(rx,ry,rz,N,L,rc2,tailcorr,ecor,shift,ecut,&vir);

  nAcc = 0;
  nTransAcc = 0;
  nTransAtt = 0;
  nVolAcc = 0;
  nVolAtt = 0;
  esum = 0.0;
  nSamp = 0;
  rho_sum = 0.0;
  p_sum = 0.0;
  for (c=0;c<nCycles;c++) {

    /* Get a random number between 0 and N */
    i=(int)gsl_rng_uniform_int(r,N+1);
    if (i<N) { /* perform a displacement move */
      nTransAtt++;
      acc = mcmove(rx,ry,rz,i,N,tailcorr,ecor,shift,ecut,
		   L,rc2,dr,beta,&E,&vir,r);
      nTransAcc+=acc;
    }
    else { /* perform a volume scaling move */
      nVolAtt++;
      acc = mcvol(rx,ry,rz,N,           /* Coordinates */
		  &L,                   /* box size may change! */
		  tailcorr,&ecor,&pcor, /* tail corrections may change */
		  shift,&ecut,&rc2,     /* shift and cutoff may change */
		  dlnV,beta,P,          /* volume delta, temp, pressure */
		  &E,&vir,              /* Energy and virial are returned */
		  r);                   /* Random number generator */
      nVolAcc+=acc;
    }
    nAcc+=acc;

    /* Sample: default frequency is once per trial move; We must
       include results of a move regardless of whether the move is
       accepted or rejected. */
    if (c>nEq) {
      esum+=E;
      p_sum+=vir/3.0/(L*L*L)+pcor;
      rho_sum+=N/(L*L*L);
      nSamp++;
    }
  }

  if (short_out)
    fprintf(stdout,"%.3lf %.3lf | %.6lf %.6lf : %.5lf %.5lf [%.5lf] %.5lf\n",
	    T,P,dr,dlnV,((double)nAcc)/(N*nCycles),
	    esum/nSamp/N,p_sum/nSamp+rho_sum/nSamp*T,rho_sum/nSamp);
  else
    fprintf(stdout,"NPT Metropolis Monte Carlo Simulation"
	    " of the Lennard-Jones fluid.\n"
	    "---------------------------------------------\n"
	    "Number of particles:              %i\n"
	    "Number of cycles:                 %i\n"
	    "Maximum particle displacement:    %.5lf\n"
	    "Maximum log-volume displacement:  %.5lf\n"
	    "Temperature:                      %.5lf\n"
	    "Pressure:                         %.5lf\n"
	    "Tail corrections used?            %s\n"
	    "Shifted potentials used?          %s\n"
	    "Results:\n"
	    "Displacement attempts:            %i\n"
	    "Volume change attempts:           %i\n"
	    "Acceptance ratio, ptcl displ.     %.5lf\n"
	    "Acceptance ratio, vol moves       %.5lf\n"
	    "Overall acceptance ratio:         %.5lf\n"
	    "Energy/particle:                  %.5lf\n"
	    "Density:                          %.5lf\n"
	    "Computed Pressure:                %.5lf\n"
	    "Program ends.\n",
	    N,nCycles,dr,dlnV,T,P,
	    tailcorr?"Yes":"No",shift?"Yes":"No",
	    nTransAtt,nVolAtt,
	    ((double)nTransAcc)/(nTransAtt?nTransAtt:1),
	    ((double)nVolAcc)/(nVolAtt?nVolAtt:1),
	    ((double)nAcc)/nCycles,
	    esum/nSamp/N,rho_sum/nSamp,
	    p_sum/nSamp+rho_sum/nSamp*T);
  return 0;
}
