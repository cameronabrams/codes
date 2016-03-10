/* 
   Metropolis Monte Carlo simulation of a 2D lattice gas

   Cameron F. Abrams

   Written for the course CHE 800-002, Adv. Chem. Engr. Thermodynamics
   Winter 1112

   compile using "gcc -o lgT lgT.c -lm"

   runs as "./lgT -L <sidelength(20)> -T <temperature(1.0)> \
                  -nc <numcycles(1e6)> -fs <samplefreq(100)> \
	          -s <seed(?)> -mu <mu(0)>"

   For example, to run a 20x20 system at T = 0.5, mu = -2, for 1e7 cycles
   sampling every 100 cycles, the command looks like
           
          ./lgT -L 20 -T 0.5 -mu -2 -nc 1e7 -fs 100

   Append the flag "-novis" if you do NOT want to watch the lattice gas!
   
   The default values are shown in parentheses above.

   The Hamiltonian is 

   H = E_\nu - \mu N_\nu = -\epsilon\sum_{<ij>} n_i n_j - \mu\sum_i n_i,

   where "sum_<ij>" is the sum over all unique
   nearest neighbor pairs, n_i = 0,1, and \epsilon
   it the strength of the attractive interaction between
   occupied neighboring sites and \mu is the chemical potential.

   The argument of the exponential in the Boltzmann factor:

   -\beta H = \beta\epsilon\sum_{ij} n_i n_j +\beta\mu\sum_i n_i

   Take the unit of energy as \epsilon and let temperature be
   non-dimensionalized as T^* = T (k_B / \epsilon), for convenience.
   Then

   \beta^* = \beta (\epsilon), and

   -\beta H = -(\beta^*) H/\epsilon 
            = (\beta^*)[\sum_{ij} n_i n_j + (\mu/\epsilon)\sum_i n_i]  

   Therefore, we specify a dimensionless temperature T* (whose
   reciprocal is \beta^*) and dimensionless chemical potential,
   \mu/\epsilon, for a system where \epsilon is the unit of energy.

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2012
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* This function computes and returns the change in (E - mu N) when
   site (i,j) occupancy is flipped.  The modulo arithmetic (the %
   operator) ensures periodic boundaries. The syntax "i?(i-1):(L-1)"
   performs the following check: If i is non-zero, return i-1,
   otherwise return L-1.  This also ensures periodic boundaries.  */
double D_E_muN ( int ** F, int L, int i, int j, double mu ) {
  int dN = -(2*F[i][j]-1); // -1 if proposed move is 1->0, 1 if move is 0->1
  double dE = -dN*(F[i?(i-1):(L-1)][j]+F[(i+1)%L][j]+
		   F[i][j?(j-1):(L-1)]+F[i][(j+1)%L]);
  double mu_dN = dN*mu;
  return (dE-mu_dN);
}

/* Sample N, N^2, E, and E^2 */
double samp ( int ** F, int L, double * n, double * n2, double * e, double * e2 ) {
  int i,j;

  *n=*n2=0.0;
  *e=*e2=0.0;

  for (i=0;i<L;i++) {
    for (j=0;j<L;j++) {
      *n+=(double)F[i][j];
      *e-=(double)(F[i][j])*(F[i][(j+1)%L]+F[(i+1)%L][j]);
    }
  }
  *n2=(*n)*(*n);
  *e2=(*e)*(*e);
}

/* Randomly assigns all occupancies */
void init ( int ** F, int L ) {
  int i,j;

  /* Visit each position (i,j) in the lattice */
  for (i=0;i<L;i++) {
    for (j=0;j<L;j++) {
      F[i][j]=(int)(lrand48()%2);
    }
  }
}

void draw ( int ** F, int L, int c, double n, double e,
	    char o, char u, int visual, struct timespec * frm ) {
  int i,j;
  double L2=L*L;

  fprintf(stdout,"%c[3;0H",27);
  fprintf(stdout,"cycle % 8d  <n> % 8.5lf <E>/(L^2) % 8.5lf\n",
	  c,n/L2,e/L2);
  if (!visual) return;
  for (i=0;i<L;i++) {
    for (j=0;j<L;j++) {
      fprintf(stdout,"% 2c",(F[i][j]==0)?u:o);
    }
    fprintf(stdout,"\n");
  }
  if (frm) nanosleep(frm,NULL);
}

int main (int argc, char * argv[]) {

  /* System parameters */
  int ** F;       /* The 2D array of sites; i.e., the "gas" */
  int L = 20;     /* The sidelength of the lattice */
  int N;          /* The total number of sites = L*L */
  double T = 1.0; /* Dimensionless temperature = (T*k)/epsilon */

  double mu = 0.0; /* Dimensionless chemical potential = mu/epsilon */

  /* Run parameters */
  int nCycles = 1000000; /* number of MC cycles to run; one cycle is N 
			    consecutive attempted occupancy flips */
  int fSamp = 1000;      /* Frequency with which samples are taken */

  /* Computational variables */
  int nSamp;      /* Number of samples taken */
  double de;      /* energy change due to flipping an occupancy */
  double b;       /* Boltzman factor */
  double x;       /* random number; used in Metropolis criterion */
  int i,j,a,c;    /* loop counters */

  /* Observables */
  double n=0.0, nsum=0.0;    /* average number of particles */
  double e=0.0, esum=0.0;    /* average energy per site */
  double n2=0.0, n2sum=0.0;    /* average squared number of particles */
  double e2=0.0, e2sum=0.0;    /* average squared energy per site */

  unsigned long int Seed = 23410981;

  char occu_char = '#';
  char empt_char = ' ';

  int visual = 1;

  struct timespec frm={0,50000000};

  /* Seed the pseudorandom number generator */
  srand48(Seed);

  /* Parse command-line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-L")) L=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i],"-mu")) mu=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nc")) nCycles = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-fs")) fSamp = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-s")) Seed = (unsigned long)atoi(argv[++i]);
    else if (!strcmp(argv[i],"-oc")) occu_char = argv[++i][0];
    else if (!strcmp(argv[i],"-ec")) empt_char = argv[++i][0];
    else if (!strcmp(argv[i],"-novis")) visual = 0;
  }
  
  /* Output some initial information */
  fprintf(stdout,"%c[2J",27);
  fprintf(stdout,"%c[0;0H",27);
  fprintf(stdout,"# Lattice gas simulation, mu-VT Metropolis Monte Carlo -- cfa 2012\n");
  fprintf(stdout,"# L = %i, T = %.3lf, mu/ep = %.3lf, nCycles %i, fSamp %i, Seed %u\n",
	  L,T,mu,nCycles,fSamp,Seed);

  /* Compute the number of spins */
  N=L*L;

  /* Allocate memory for the magnet */
  F=(int**)malloc(L*sizeof(int*));
  for (i=0;i<L;i++) F[i]=(int*)malloc(L*sizeof(int));

  /* Generate an initial state */
  init(F,L);

  /* For computational efficiency, convert T to reciprocal T */
  T=1.0/T;

  n = n2 = 0.0;
  e = e2 = 0.0;
  nSamp = 0;
  for (c=0;c<nCycles;c++) {
    /* Make N flip attempts */
    for (a=0;a<N;a++) {
      /* randomly select a site */
      i=(int)(drand48()*L);
      j=(int)(drand48()*L);
      /* get the change in (E - mu N) due to flipping occupancy (i,j) */
      de = D_E_muN(F,L,i,j,mu);
      /* compute the Boltzmann factor; recall T is now reciprocal
         temperature */
      if (de<0.0) {/* accept */
	/* flip it */
	F[i][j]=(F[i][j]+1)%2;
      } else {
	b = exp(-de*T);
	/* pick a random number between 0 and 1 */
	x = drand48();
	/* accept or reject this flip */
	if (x<b) { /* accept */
	  /* flip it */
	  F[i][j]=(F[i][j]+1)%2;
	}
      }
    }
    /* Sample and accumulate averages */
    if (!(c%fSamp)) {
      samp(F,L,&n,&n2,&e,&e2);
      nsum+=n;
      esum+=e;
      n2sum+=n2;
      e2sum+=e2;
      nSamp++;
      draw(F,L,c,nsum/nSamp,esum/nSamp,occu_char,empt_char,visual,&frm);
      fflush(stdout);
    }
  }
  fprintf(stdout,"# The average occupancy is % .5lf +/- % .5lf\n",nsum/nSamp/(L*L), sqrt((n2sum-nsum*nsum/nSamp)/nSamp)/(L*L));
  fprintf(stdout,"# The average energy per site is % .5lf +/- % .5lf\n",esum/nSamp/(L*L),sqrt((e2sum-esum*esum/nSamp)/nSamp)/(L*L));
  fprintf(stdout,"# Program ends.\n");
}
