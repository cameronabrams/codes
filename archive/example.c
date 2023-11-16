#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rotmov.h"


/* example main program for demonstrating crankshaft and generalized pivot moves

   written by Cameron Abrams
   2004
*/


/* populates the array of points r with a random walk */
void make_random_walk ( double ** r, int N ) {

  int i;
  double phi, costheta, sintheta, cosphi, sinphi, dx, dy, dz;
  r[0][0]=r[0][1]=r[0][2]=0.0;
  

  for (i=1;i<N;i++) {
    phi=2*M_PI*(double)rand()/RAND_MAX;
    costheta = 1.0 - 2*(double)rand()/RAND_MAX;
    sintheta = sqrt(1-costheta*costheta);
    cosphi=cos(phi);
    sinphi=sin(phi);
    dx=sintheta*cosphi;
    dy=sintheta*sinphi;
    dz=costheta;
    r[i][0]=r[i-1][0]+dx;
    r[i][1]=r[i-1][1]+dy;
    r[i][2]=r[i-1][2]+dz;
  }
}

void echo_bond_lengths ( double ** r, int N ) {

  int i;
  double rabs;

  for (i=0;i<N-1;i++) {
    rabs=sqrt((r[i][0]-r[i+1][0])*(r[i][0]-r[i+1][0])+
	      (r[i][1]-r[i+1][1])*(r[i][1]-r[i+1][1])+
	      (r[i][2]-r[i+1][2])*(r[i][2]-r[i+1][2]));
    fprintf(stderr,"# bond %i-%i: %.3lf\n",i,i+1,rabs);      
  }

}

/* outputs the random walk; monomers between a and b are "highlighted" */
void xyz_out_highlight_subchain (FILE * fp, double **r, int N, int a, int b ) {
  
  int i,typ;

  if (a>b) {
    int tmp=a;
    a=b;
    b=tmp;
  }

  fprintf(fp,"%i 0 0\n\n",N);
  for (i=0;i<N;i++) {
    typ=1;
    if (i==a || i==b) typ=5;
    if (i>a && i<b) typ=6;
    fprintf(fp,"%i %.4lf %.4lf %.4lf\n", typ,r[i][0], r[i][1], r[i][2]);
  }
}

/* outputs the random walk; monomers between a and b are "highlighted" */
void xyz_out_highlight_subchain_plus1 (FILE * fp, double **r, int N, int a, int b, int c ) {
  
  int i,typ;

  if (a>b) {
    int tmp=a;
    a=b;
    b=tmp;
  }

  fprintf(fp,"%i 0 0\n\n",N);
  for (i=0;i<N;i++) {
    typ=1;
    if (i==a || i==b) typ=5;
    if (i>a && i<b) typ=6;
    if (i==c) typ=7;
    fprintf(fp,"%i %.4lf %.4lf %.4lf\n", typ,r[i][0], r[i][1], r[i][2]);
  }
}


int main (int argc, char * argv[]) {

  int i,seed =101;
  int N = 100;
  int a = 80;
  int b = 90;
  int c = -1;
  double g1=M_PI, g2=M_PI;
  double ** r;
  FILE * fp;

  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-seed")) seed=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
    /* a->b crank by g1 */
    else if (!strcmp(argv[i],"-ab")) sscanf(argv[++i],"%i,%i,%lf",&a,&b,&g1);
    /* a->b->c generalized pivots by g1 and g2 */
    else if (!strcmp(argv[i],"-abc")) sscanf(argv[++i],"%i,%i,%i,%lf,%lf",&a,&b,&c,&g1,&g2);
  }

  /* seed the random number generator */
  srand(seed);

  /* setup the rotation matrices */
  setup_rotation_matrices();

  /* allocate the array of points */
  r = (double **)malloc(N*sizeof(double*));
  for (i=0;i<N;i++) r[i]=(double*)calloc(3,sizeof(double));

  /* create the random walk */
  make_random_walk(r,N);
  //echo_bond_lengths(r,N);

  /* perform the desired move; if c is set by the user, assume it
     is to be a generalized pivot; otherwise it is a crankshaft */
  if (c==-1) {
    fprintf(stderr,"# crank\n");

    /* output original random walk */
    fp=fopen("tmp0.xyz","w");
    xyz_out_highlight_subchain(fp,r,N,a,b);
    fclose(fp);

    /* crank it */
    crankshaft(r,N,a,b,g1);
    echo_bond_lengths(r,N);

    /* output the resulting random walk */
    fp=fopen("tmp1.xyz","w");
    xyz_out_highlight_subchain(fp,r,N,a,b);
    fclose(fp);
  }
  else {
    fprintf(stderr,"# pivot\n");

    fp=fopen("tmp0.xyz","w");
    xyz_out_highlight_subchain_plus1(fp,r,N,b,c,a);
    fclose(fp);

    /* pivot it */
    generalized_pivot(r,N,a,b,c,g1,g2);
    //echo_bond_lengths(r,N);

    fp=fopen("tmp1.xyz","w");
    xyz_out_highlight_subchain_plus1(fp,r,N,b,c,a);
    fclose(fp);
  }

}
