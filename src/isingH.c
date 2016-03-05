/* 
   Metropolis Monte Carlo simulation of a 2D Ising system

   Cameron F. Abrams

   Written for the course CHE 614 Chem. Engr. Thermodynamics II
   Winter 1516

   compile using "gcc -o ising ising.c -lm"

   runs as "./ising -L <sidelength(20)> -T <temperature(1.0)> \
                    -nc <numcycles(1e6)> -fs <samplefreq(100)> \
		    -s <seed(?)>"

   For example, to run a 20x20 system at T = 0.5 for 1e7 cycles
   sampling every 100 cycles, the command looks like
           
          ./ising -L 20 -T 0.5 -nc 1e7 -fs 100
   
   The default values are shown in parentheses above.

   The Hamiltonian is 

   H = -J sum_<ij> s_i * s_j + h * s_i,

   where "sum_<ij>" is the sum over all unique
   nearest neighbor pairs, s_i = +/- 1, and J 
   is an arbitrary "coupling" parameter having
   units of energy.  We work in a unit system
   where energy is measured in units of J and
   temperature is nondimensionalized as (k_B*T)/J.

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2011-2016
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* This function computes and returns the change in
   system energy when spin (i,j) is flipped.  The
   modulo arithmetic (the % operator) ensures 
   periodic boundaries. The syntax "i?(i-1):(L-1)"
   performs the following check:  If i is non-zero,
   return i-1, otherwise return L-1.  This also 
   ensures periodic boundaries.  */
double E ( int ** F, double ** h, int L, int i, int j ) {
  return -2.0*(F[i][j])*((double)(F[i?(i-1):(L-1)][j]+F[(i+1)%L][j]+
				F[i][j?(j-1):(L-1)]+F[i][(j+1)%L])+h[i][j]);
}

/* Sample the system; compute the average magnetization
   and the average energy per spin */
double samp ( int ** F, double ** h, int L, double * s, double * e, double * sisj, int cc ) {
  int i,j,d;

  *s=0.0;
  *e=0.0;
  for (i=0;i<L/2;i++) sisj[i]=0.0;
  /* Visit each position (i,j) in the lattice */
  for (i=0;i<L;i++) {
    for (j=0;j<L;j++) {
      *e-=(double)(F[i][j])*(F[i][(j+1)%L]+F[(i+1)%L][j]+h[i][j]);
      if (j==cc) {
	for (d=0;d<L/2;d++) 
	  sisj[d]+=F[i][j]*(F[(i+d)%L][j]+F[i][(j+d)%L]);
	*s+=(double)F[i][j];
      }
    }
  }
  *s/=(L);
  *e/=(L*L);
  for (d=0;d<L/2;d++) 
    sisj[d]/=(2*L);
}

/* Randomly assigns all spins */
void init ( int ** F, int L ) {
  int i,j;

  /* Visit each position (i,j) in the lattice */
  for (i=0;i<L;i++) {
    for (j=0;j<L;j++) {
      /* 2*x-1, where x is randomly 0,1.
         lrand48() returns a non-negative integer, and
         %2 takes its modulo-2; if even, this zero, if odd, 1. */
      F[i][j]=2*(int)(lrand48()%2)-1;
    }
  }
}

void init_field ( double ** h, int L, double H ) {
  int i,j;
  for (i=0;i<L;i++) {
    for (j=0;j<L;j++) {
      if (!i && j<L/2) h[i][j]=H;
      else if (!i && j>=L/2) h[i][j]=-H;
      else h[i][j]=0.0;
    }
  }
}

int main (int argc, char * argv[]) {

  /* System parameters */
  int ** F;       /* The 2D array of spins; i.e., the "magnet" */
  double ** h, H=0;    /* The field; H is the magnitude and h[][] is
			computed at each spin position */
  int L = 20;     /* The sidelength of the magnet */
  int N;          /* The total number of spins = L*L */
  double T = 1.0; /* Dimensionless temperature = (T*k)/J */

  /* Run parameters */
  int nCycles = 1000000; /* number of MC cycles to run; one cycle is N 
			    consecutive attempted spin flips */
  int fSamp = 1000;      /* Frequency with which samples are taken */

  /* Computational variables */
  int nSamp;      /* Number of samples taken */
  int de;         /* energy change due to flipping a spin */
  double b;       /* Boltzman factor */
  double x;       /* random number */
  int i,j,a,c;    /* loop counters */

  /* Observables */
  double s=0.0, ssum=0.0;    /* average magnetization */
  double e=0.0, esum=0.0;    /* average energy per spin */
  double * sisj, * sisjsum;  /* sisj array */
  int cc = 10;               /* column in which correlation function is computed */

  unsigned long int Seed = 23410981;

  /* Here we parse the command line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-L")) L=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i],"-H")) H=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nc")) nCycles = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-fs")) fSamp = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-s")) Seed = (unsigned long)atoi(argv[++i]);
    else if (!strcmp(argv[i],"-cc")) cc=atoi(argv[++i]);
  }
  
  /* Seed the pseudorandom number generator */
  srand48(Seed);

  /* Output some initial information */
  fprintf(stdout,"# Ising simulation, NVT Metropolis Monte Carlo -- cfa 2016\n");
  fprintf(stdout,"# Field: inhomogeneous, +H for spins on left half of top row; -H for spins on right half\n");
  fprintf(stdout,"# L = %i, T = %.3lf, H = %.5lf, nCycles %i, fSamp %i, Seed %u\n",
	  L,T,H,nCycles,fSamp,Seed);

  /* Compute the number of spins */
  N=L*L;

  /* Allocate memory for the magnet */
  F=(int**)malloc(L*sizeof(int*));
  for (i=0;i<L;i++) F[i]=(int*)malloc(L*sizeof(int));

  /* Allocate memory for the field */
  h=(double**)malloc(L*sizeof(double*));
  for (i=0;i<L;i++) h[i]=(double*)malloc(L*sizeof(double));

  /* Allocate memory for sisj array */
  sisj=(double*)malloc(L/2*sizeof(double));
  sisjsum=(double*)malloc(L/2*sizeof(double));

  /* Set up the initial field */
  init_field(h,L,H);

  /* Generate an initial state */
  init(F,L);

  /* For computational efficiency, convert T to reciprocal T */
  T=1.0/T;

  s = 0.0;
  e = 0.0;
  for (i=0;i<L/2;i++) sisj[i]=sisjsum[i]=0.0;
  nSamp = 0;
  for (c=0;c<nCycles;c++) {
    /* Make N flip attempts */
    for (a=0;a<N;a++) {
      /* randomly select a spin */
      i=(int)(drand48()*L);
      j=(int)(drand48()*L);
      /* get the "new" energy as the incremental change due
         to flipping spin (i,j) */
      de = E(F,h,L,i,j);
      /* compute the Boltzmann factor; recall T is now
         reciprocal temperature */
      // fprintf(stdout,"de*T %.10le\n",de*T);fflush(stdout);
      b = exp(de*T);
      /* pick a random number between 0 and 1 */
      x = drand48();
      /* accept or reject this flip */
      if (x<b) { /* accept */
	/* flip it */
	F[i][j]*=-1;
      }
    }
    /* Sample and accumulate averages */
    if (!(c%fSamp)) {
      samp(F,h,L,&s,&e,sisj,cc);
      ssum+=s;
      esum+=e;
      for (i=0;i<L/2;i++) sisjsum[i]+=sisj[i];
      nSamp++;
    }
  }
  fprintf(stdout,"# The average magnetization is in column %i %.5lf\n",cc,ssum/nSamp);
  fprintf(stdout,"# The average energy per spin in is %.5lf\n",esum/nSamp);
  fprintf(stdout,"# Spin-spin correlation function, column %i\n# dist  <sisj> cij\n",cc);
  for (i=0;i<L/2;i++)
    fprintf(stdout,"%i  %.5lf %.5lf\n",i,sisjsum[i]/nSamp,sisjsum[i]/nSamp-ssum*ssum/(nSamp*nSamp));
  fprintf(stdout,"# Program ends.\n");
}
