// 
// tourney.c 
// Cameron F Abrams cfa22@drexel.edu
//
// Tournament simulator
//
// computes the number of underdogs (type-0) remaining
// in each round of a single-elimination tournament
// initialized with equal numbers of favorites (type-1) and 
// underdogs.  The variable "usp" is the "upset probability"
// which is the probability that an underdog defeats a favorite
// in any one game.  The number of rounds can be specified by "-nr"
// and the number of individual simulations by "-ns".
//
// (c) 201i8 cfa
//
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// does one round of the tournament, storing winning state 
// in S2 for the result of the game between adjacent teams
// in S1
void do_round ( int * S1, int * S2, int N, double usp ) {
    double x;
    int i;
    for (i=0;i<N;i+=2) {
	if (S1[i]==S1[i+1]) S2[i/2]=S1[i];
	else {
           x=drand48();
	   if (x<usp) S2[i/2]=0;
	   else S2[i/2]=1;
	}
    }
}

// increments the nz[] array at position g by the number
// of zeros in the the S1[] array
void samp ( int * S1, int N, int * nz, int g ) {
   int i;
   for (i=0;i<N;i++) if (S1[i]==0) nz[g]++;
}

// initializes the first round randomly with statistically
// equally likely underdogs and favorites in each slot
void do_init_a ( int * S, int N ) {
   int i;
   for (i=0;i<N;i++) {
      if (drand48()<0.5) S[i]=0;
      else S[i]=1;
   }
}

int main ( int argc, char * argv[] ) {
   int * S1;
   int * S2;
   int * Stmp;
   int N0, N;
   int nr=10;
   unsigned long int Seed=30479;
   int * nz;
   int ns=1000;
   int i,j;
   double usp=0.0;

   for (i=1;i<argc;i++) {
     if (!strcmp(argv[i],"-nr")) nr=atoi(argv[++i]);
     else if (!strcmp(argv[i],"-ns")) ns=atoi(argv[++i]);
     else if (!strcmp(argv[i],"-usp")) usp=atof(argv[++i]);
     else if (!strcmp(argv[i],"-seed")) Seed=(unsigned long int)atoi(argv[++i]);
   }

   srand48(Seed);

   N0=(int)pow(2,nr-1);
   N=N0;
   S1=(int*)malloc(N*sizeof(int));
   S2=(int*)malloc(N*sizeof(int));
   nz=(int*)malloc(nr*sizeof(int));
   for (i=0;i<nr;i++) nz[i]=0;

   fprintf(stdout,"# nsim %i usp %.5lf N0 %i\n",ns,usp,N0);fflush(stdout);

   for (i=0;i<ns;i++) {
     N=N0;
     do_init_a(S1,N);
     samp(S1,N,nz,0);
     for (j=1;j<nr;j++) {
	do_round(S1,S2,N,usp);
	Stmp=S1;
	S1=S2;
	S2=Stmp;
	N/=2;
        samp(S1,N,nz,j);
     }
   }

   N=N0;
   for (i=0;i<nr;i++) {
     fprintf(stdout,"%i %i %.8lf\n",i,N,((double)nz[i])/ns/(double)N);
     N/=2;
   }
}
