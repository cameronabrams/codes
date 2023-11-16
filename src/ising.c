/* 
   Metropolis Monte Carlo simulation of a 2D Ising system

   Cameron F. Abrams

   Written for the course CHE T580, Modern Molecular Simulations
   Spring 2021

   compile using "gcc -o ising ising.c -lm -lgsl"
   (assumes the GNU Scientific Library is installed)

   runs as "./ising -L <sidelength(20)> -T <temperature(1.0)> \
                    -nc <numcycles(1e6)> -fs <samplefreq(100)> \
		    -s <seed(?)>"

   For example, to run a 20x20 system at T = 0.5 for 1e7 cycles
   sampling every 100 cycles, the command looks like
           
          ./ising -L 20 -T 0.5 -nc 1e7 -fs 100
   
   The default values are shown in parentheses above.

   You must have the GNU Scientific Library installed.

   The Hamiltonian is 

   H = -J sum_<ij> s_i * s_j,

   where "sum_<ij>" is the sum over all unique
   nearest neighbor pairs, s_i = +/- 1, and J 
   is an arbitrary "coupling" parameter having
   units of energy.  We work in a unit system
   where energy is measured in units of J and
   temperature is nondimensionalized as (k_B*T)/J.

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2004-2021
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

/* dE() returns the change in energy when spin (i,j) 
   is flipped. Modulo arithmetic (the % operator) ensures 
   periodic boundaries. The syntax "i?(i-1):(L-1)"
   performs the following check:  If i is non-zero,
   return i-1, otherwise return L-1.  This also 
   ensures periodic boundaries.  
   
   Parameters:
    M: 2D array of spins (the magnet)
    L: side-length of magnet
    i,j: coordinates of proposal spin
   
   */
int dE ( int ** M, int L, int i, int j ) {
  return -2*(M[i][j])*(M[i?(i-1):(L-1)][j]+M[(i+1)%L][j]+
		      M[i][j?(j-1):(L-1)]+M[i][(j+1)%L]);
}

/* Sample the system; compute the average magnetization
   and the average energy per spin */
double samp ( int ** M, int L, double * s, double * e ) {
  int i,j;

  *s=0.0;
  *e=0.0;
  /* Visit each position (i,j) in the lattice */
  for (i=0;i<L;i++) {
    for (j=0;j<L;j++) {
      *s+=(double)M[i][j];
      *e-=(double)(M[i][j])*(M[i][(j+1)%L]+M[(i+1)%L][j]);
    }
  }
  *s/=(L*L);
  *e/=(L*L);
}

/* Randomly assigns all spins */
void init ( int ** M, int L, gsl_rng * r ) {
  int i,j;

  /* Visit each position (i,j) in the lattice */
  for (i=0;i<L;i++) {
    for (j=0;j<L;j++) {
      /* 2*x-1, where x is randomly 0,1 */
      M[i][j]=2*(int)gsl_rng_uniform_int(r,2)-1;
    }
  }
}

int main (int argc, char * argv[]) {

  /* System parameters */
  int ** M;       /* The 2D array of spins */
  int L = 20;     /* The sidelength of the array */
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

  /* This line creates a random number generator
     of the "Mersenne Twister" type, which is much
     better than the default random number generator. */
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned long int Seed = 23410981;

  /* Here we parse the command line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-L")) L=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nc")) nCycles = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-fs")) fSamp = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-s")) Seed = (unsigned long)atoi(argv[++i]);
  }
  
  /* Output some initial information */
  fprintf(stdout,"# command: ");
  for (i=0;i<argc;i++) fprintf(stdout,"%s ",argv[i]);
  fprintf(stdout,"\n");
  fprintf(stdout,"# ISING simulation, NVT Metropolis Monte Carlo\n");
  fprintf(stdout,"# L = %i, T = %.3lf, nCycles %i, fSamp %i, Seed %u\n",
	  L,T,nCycles,fSamp,Seed);

  /* Seed the random number generator */
  gsl_rng_set(r,Seed);

  /* Compute the number of spins */
  N=L*L;

  /* Allocate memory for the system */
  M=(int**)malloc(L*sizeof(int*));
  for (i=0;i<L;i++) F[i]=(int*)malloc(L*sizeof(int));

  /* Generate an initial state */
  init(M,L,r);

  /* For computational efficiency, convert T to reciprocal T */
  T=1.0/T;

  s = 0.0;
  e = 0.0;
  nSamp = 0;
  for (c=0;c<nCycles;c++) {
    /* Make N flip attempts */
    for (a=0;a<N;a++) {
      /* randomly select a spin */
      i=(int)gsl_rng_uniform_int(r,L);
      j=(int)gsl_rng_uniform_int(r,L);
      /* get the "new" energy as the incremental change due
         to flipping spin (i,j) */
      de = dE(F,L,i,j);
      /* compute the Boltzmann factor; recall T is now
         reciprocal temperature */
      b = exp(de*T);
      /* pick a random number between 0 and 1 */
      x = gsl_rng_uniform(r);
      /* accept or reject this flip */
      if (x<b) { /* accept */
	      /* flip it */
	      M[i][j]*=-1;
      }
    }
    /* Sample and accumulate averages */
    if (!(c%fSamp)) {
      samp(M,L,&s,&e);
      fprintf(stdout,"%i %.5le %.5le\n",c,s,e);
      fflush(stdout);
      ssum+=s;
      esum+=e;
      nSamp++;
    }
  }
  fprintf(stdout,"# The average magnetization is %.5lf\n",ssum/nSamp);
  fprintf(stdout,"# The average energy per spin is %.5lf\n",esum/nSamp);
  fprintf(stdout,"# Program ends.\n");
}
