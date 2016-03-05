/* 
   Blocked averages computation a la Flyvbjerg and Petersen.
   [J. Chem. Phys., 91:461-466,1989]

   Reads in a Y vs. X data, and recursively block adjacent 
   elements, computing variance.  Details:  See Appendix D
   in F&S.

   Outputs column-oriented data as
   <generation> <mean> <variance>

   Cameron F. Abrams

   Written for the course CHE 800-002, Molecular Simulation
   Spring 0304

   compile using "gcc -o flyvbjerg.c flyvbjerg -lm"

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2004
*/
#include <stdio.h>
#include <stdlib.h>


/* Perform one blocking operation on the array A[] */
void block (double * A, int * N) {
  int i;
  for (i=0;i<(*N)/2;i++) {
    A[i] = 0.5*(A[2*i]+A[2*i+1]);
  }
  for (;i<(*N);i++) A[i]=0.0;
  (*N) = (*N) / 2;
}

/* Compute the mean (m), 2nd moment (s2) and 4th moment (s4) of the
   array A[] */
void moments (double * A, int N, double * m, double * s2, double * s4) {
  double ss=0.0,ssss=0.0;
  int i;

  *m=0.0;
  for (i=0;i<N;i++) {
    *m+=A[i];
    ss+=A[i]*A[i];
    ssss+=A[i]*A[i]*A[i]*A[i];
  }
  *m/=N; ss/=N; ssss/=N;
  *s2 = ss - *m**m;
  *s4 = ssss - *m**m**m**m;
}

int main (int argc, char * argv[] ) {
  int i,j,dum;
  int nData;
  int ngen;
  int nrml=1;
  int data_limit=0;
  double * A;
  double mn,s2,s4,av,sav;
  char * fn=NULL, ln[255];
  FILE * fp;

  for (i=1;i<argc;i++) {  
    if (!strcmp(argv[i],"-ng")) ngen=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-na")) nrml=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-lim")) data_limit=atoi(argv[++i]);
    else fn=argv[i];
  }
  fp=fopen(fn,"r");
  fprintf(stderr,"# determining datasize\n");
  nData=0;
  while (fgets(ln,255,fp)&&(!data_limit||nData<data_limit))
    if (ln[0]!='#') nData++;
  fprintf(stderr,"# %d samples = %d bytes\n",nData,nData*sizeof(double));
  fflush(stderr);
  rewind(fp);
  A = (double*)calloc(nData,sizeof(double));
  fprintf(stderr,"# scanning...");fflush(stderr);
  i=0;
  for (i=0;i<nData;i++) {
    fgets(ln,255,fp);
    if (ln[0]!='#') {
      sscanf(ln,"%lf\n",&A[i]);
      A[i]/=nrml;
    }
  }
  fprintf(stderr,"done\n");
  j=0;
  while (nData>4) {
    moments(A,nData,&mn,&s2,&s4);
    av=sqrt(s2/(nData-1));
    sav=av/sqrt(2*(nData-1));
    fprintf(stdout,"%d %.5lf %.5lf %.5lf\n",j++,mn,av,sav);
    block(A,&nData);
  }
  fprintf(stderr,"# program ends\n");
}
