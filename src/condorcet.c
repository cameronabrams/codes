#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>

/* vp is the probability that a voter selects correctly. */
int vote ( double vp, gsl_rng * r ) {
  return (int)(gsl_rng_uniform(r)<vp);
}

int voteround ( double vp, int nv, gsl_rng * r ) {
  int iv;
  double result=0.0;
  for (iv=0; iv<nv; iv++) result+=vote(vp,r);
  return (int)(result/nv>0.5);
}

int main ( int argc, char * argv[] ) {
  double vp=0.51, p;
  int i, rnd;
  int nv=100, iv;
  int nr=100, ir;
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  double pav=0.0;

  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-vp")) vp =atof(argv[++i]);
    else if (!strcmp(argv[i],"-nv")) nv =atoi(argv[++i]);
    else if (!strcmp(argv[i],"-nr")) nr =atoi(argv[++i]);
  }
  
  pav=0.0;
  for (ir=0;ir<nr;ir++) {
    rnd=voteround(vp,nv,r);
    //fprintf(stderr,"round %i vote %i\n",ir,rnd);
    pav+=rnd;
  }
  pav/=nr;

  fprintf(stdout,"# Result for %.10lf, %i, %i is %.5lf\n",vp,nv,nr,pav);

}
