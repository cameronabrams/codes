/* 
   Metropolis Monte Carlo simulation of a 2D Ising system

   Cameron F. Abrams

   Written for the CHE 614, Chem. Engr. Thermo. II
   Winter 1516

   compile using "gcc -o isingU isingU.c -lm"

   runs as "./isingU -L <sidelength(20)> -T <temperature(1.0)> \
                    -nc <numcycles(1e6)> -fs <samplefreq(100)> \
		    -s <seed(?)> -window <M_min,M_max>"

   For example, to run a 20x20 system at T = 0.5 for 1e7 cycles
   sampling every 100 cycles, confined to a spin window 
   between -98 and -100, the command looks like
           
          ./isingU -L 20 -T 0.5 -nc 1e7 -fs 100 -window -100,-98
   
   The default values are shown in parentheses above.

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
double E ( int ** F, int L, double h, int i, int j ) {
  return -2.0*(F[i][j])*((double)(F[i?(i-1):(L-1)][j]+F[(i+1)%L][j]+
				F[i][j?(j-1):(L-1)]+F[i][(j+1)%L])+h);
}

/* Sample the system; compute the magnetization M
   and the average energy per spin e */
void samp ( int ** F, int L, double h, int * M, double * e ) {
  int i,j;
  *M=0;
  *e=0.0;
  for (i=0;i<L;i++) {
    for (j=0;j<L;j++) {
      *M+=F[i][j];
      *e-=(F[i][j])*(F[i][(j+1)%L]+F[(i+1)%L][j]+h);
    }
  }
}

/* Randomly assigns all spins such that initial M is between M_min and
   M_max */
void init ( int ** F, int L, int M_min, int M_max ) {
  int i,j;

  /* compute required number of up spins */
  int Nup=0.5*(0.5*(M_max+M_min)+L*L);
  int nup=0;
  /* flip every spin down */
  for (i=0;i<L;i++) {
    for (j=0;j<L;j++) {
      F[i][j]=-1;
    }
  }

  /* randomly flip the required number of spins up */
  while (nup<Nup) {
    i=(int)(drand48()*L);
    j=(int)(drand48()*L);
    if (F[i][j]==-1) {
      F[i][j]=1;
      nup++;
    }
  }  
}

int runMC_and_update_hist ( int ** F, int L, double h, int m, int ws, double *mag_hist, int nbin, int lbin, double T, int nCycles, int fSamp ) {

  int M_min = m;
  int M_max = m+ws;

  int ntb = ws/2 + 1;
  int * th;

  double sf = 1.0;
  int N = L*L;


  /* Computational variables */
  int nSamp;      /* Number of samples taken */
  double de;         /* energy change due to flipping a spin */
  double b;       /* Boltzman factor */
  double x;       /* random number */
  int i,j,a,c;    /* loop counters */

  /* Observables */
  int M;  // number of up spins minus number of down spins
  int msum=0;
  double e=0.0, esum=0.0;    /* average energy per spin */

  th = (int*) malloc(ntb*sizeof(int));
  
  for (i=0;i<ntb;i++) th[i]=0;

  /* Generate an initial state */
  init(F,L,M_min,M_max);

  /* For computational efficiency, convert T to reciprocal T */
  T=1.0/T;

  M = 0;
  e = 0.0;
  nSamp = 0;
  samp(F,L,h,&M,&e);
  //if (!quiet) fprintf(stderr,"# initial mag %i\n",M);
  //  exit(-1);
  for (c=0;c<nCycles;c++) {
    /* Make N flip attempts */
    for (a=0;a<N;a++) {
      /* randomly select a spin */
      i=(int)(drand48()*L);
      j=(int)(drand48()*L);
      /* compute the change in M due to flipping this spin */
      M+=-F[i][j]*2;
      /* if the trial value for new M is within the window, proceed */
      if (M>=M_min && M<=M_max) {
	/* get the "new" energy as the incremental change due
	   to flipping spin (i,j) */
	de = E(F,L,h,i,j);
	/* compute the Boltzmann factor; recall T is now
	   reciprocal temperature */
	b = exp(de*T);
	/* pick a random number between 0 and 1 */
	x = drand48();
	/* accept or reject this flip */
	if (x<b) { /* accept */
	  /* flip it */
	  F[i][j]*=-1;
	}
	else { /* undo the trial change in spin */
	  M+=F[i][j]*2;
	}
      }
      else { /* undo the trial change in spin */
	M+=F[i][j]*2;
      }
    }
    /* update the histogram */
    th[(M-M_min)/2]++;
    /* Sample and accumulate averages */
    if (!(c%fSamp)) {
      int tmpM;
      samp(F,L,h,&tmpM,&e);
      if (tmpM!=M) {
	fprintf(stderr,"ERROR: mismatch %i %i\n",M,tmpM);
	exit(-1);
      }
      msum+=M;
      esum+=e;
      nSamp++;
      //fprintf(stdout,"% 8.5lf % 8.5lf\n",(double)msum/(L*L),e);
      //fflush(stdout);
    }
  }

  // stitch th onto mag_hist
  if (!lbin) { // this is first
    for (i=0;i<ntb;i++) {
      mag_hist[i]=th[i];
    }
  } else {
    // match th[0] to mag_hist[lbin]
    sf = mag_hist[lbin]/th[0];
    for (i=1;i<ntb;i++) {
      mag_hist[lbin+i]=th[i]*sf;
    }
  }

  free(th);
  return 0;
}


int main ( int argc, char * argv [] ) {

    /* System parameters */
  int ** F;       /* The 2D array of spins; i.e., the "magnet" */
  int L = 20;     /* The sidelength of the magnet */
  int N;          /* The total number of spins = L*L */
  double T = 1.0; /* Dimensionless temperature = (T*k)/J */
  double h = 0.0; /* field */

  /* Run parameters */
  int nCycles = 1000000; /* number of MC cycles to run; one cycle is N 
			    consecutive attempted spin flips */
  int fSamp = 1000;      /* Frequency with which samples are taken */

  /* Umbrella sampling parameters*/
  int M_min, M_max; /* magnetization window boundaries */
  double * mag_hist;
  int nbin;
  char * my_filename="hist.dat";
  FILE *fp;
  int i;
  int ws;
  int m;
  double sum;

  unsigned long int Seed = 23410981;

  int quiet=0;

  /* Here we parse the command line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-L")) L=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i],"-h")) h=atof(argv[++i]);
    else if (!strcmp(argv[i],"-hf")) my_filename=argv[++i];
    else if (!strcmp(argv[i],"-nc")) nCycles = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-q")) quiet=1;
    else if (!strcmp(argv[i],"-fs")) fSamp = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-s")) Seed = (unsigned long)atoi(argv[++i]);
    else if (!strcmp(argv[i],"-windowsize")) ws = atoi(argv[++i]);
  }

  /* half-domain histogram */
  M_min=0;
  M_max=L*L;

  /* allocate the histogram */
  nbin=(M_max-M_min)/2+1;
  mag_hist=(double*)malloc(nbin*sizeof(double));
  for (i=0;i<nbin;i++) mag_hist[i]=0;

  /* Seed the pseudorandom number generator */
  srand48(Seed);
  
  /* Output some initial information */
  if (!quiet) {
    fprintf(stdout,"# Ising simulation, NVT Metropolis Monte Carlo Umbrella Sampling -- cfa 2016\n");
    fprintf(stdout,"# L = %i, T = %.3lf, h = %.3lf, nCycles %i, fSamp %i, Seed %u, Windowsize %i\n",
	    L,T,h,nCycles,fSamp,Seed,ws);
  }
  /* Compute the number of spins */
  N=L*L;

  /* Allocate memory for the magnet */
  F=(int**)malloc(L*sizeof(int*));
  for (i=0;i<L;i++) F[i]=(int*)malloc(L*sizeof(int));
  
  for (m = 0; m < M_max; m+=ws) {
    fprintf(stdout,"# Running window [%i,%i]...\n",m,m+ws);
    runMC_and_update_hist(F,L,h,m,ws,mag_hist,nbin,(m+ws)/2-1,T,nCycles,fSamp);
  }

  sum=0.0;
  for (i=0;i<nbin;i++) sum+=mag_hist[i];
  fp=fopen(my_filename,"w");
  for (i=0;i<nbin;i++) {
    fprintf(fp,"%i %.5lf %.5le\n",i,mag_hist[i]/sum/2,-T*log(mag_hist[i]/sum/2));
  }
  fclose(fp);
  fprintf(stdout,"Program ends. %s created.\n",my_filename);

}
