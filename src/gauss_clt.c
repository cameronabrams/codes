#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

double GaussRand_CLT (int nTerms, double var, int even) {
  int i;
  double sum;

  sum = 0.0;
  for (i = 0; i < nTerms; i++)
    sum += ((double)rand())/((double)RAND_MAX);
  sum=var*sqrt(12*nTerms)*(sum/nTerms - 0.5);

  return even?sum<0.0?-sum:sum:sum;
}

int main (int argc, char * argv[]) {

  int nSamples=0, nTerms=0;
  double var=0,dx=0.0,x;
  int i=0, N=0;
  int even=0;

  int seed=186091;
  
  for (i=0;i<argc;i++) {
    if (!strcmp(argv[i],"-s")) nSamples=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-t")) nTerms=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-v")) var=atof(argv[++i]);
    else if (!strcmp(argv[i],"-seed")) seed=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-e")) even=1;
  }

  srand(seed);

  for (i=0;i<nSamples;i++) {
    x=GaussRand_CLT(nTerms,var,even);
    fprintf(stdout,"%.5lf\n",x);
  }

  return 0; 

}
