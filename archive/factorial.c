#include <stdio.h>
#include <gsl/gsl_sf.h>
int main (int argc, char * argv[]) {
  int N=atoi(argv[1]);
  int M=2;
  int M1, M2;
  double O=gsl_sf_pow_int(M,N);
  int i;

  for (i=0;i<=N;i++) {
    M1=i;
    M2=N-i;
    fprintf(stdout,"%.3f %.6le\n",((double)i-N/2)/N,gsl_sf_choose(N,M1)/O*N);
  }
}
