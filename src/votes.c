#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* 
 Simulation of random fractional voting among seven voters for three candidates.  
 Objective is to determine the likelihood that all three candidates have scores greater than 2.

(c) 2018, cameron f abrams
cfa22@drexel.edu
*/

int main ( int argc, char * argv[] ) {
  double *x, *y; // array of fractional votes for candidate X and Y, respectively.
  // vote for candidate Z by voter i is 1-x[i]-y[i].
  int np=7;  // number of voters
  int i,j;
  int ng=0;  // number of simulations in which all three candidates have scores > 2
  int n;
  int nsim=10000000; // number of independent simulations to run
  double vx, vy, vz;
  unsigned long int Seed = 2948572;

  /* command-line parsing */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-seed")) Seed=(unsigned long int)atoi(argv[++i]);
    else if (!strcmp(argv[i],"-nsim")) nsim=atoi(argv[++i]);
  }

  x=(double*)malloc(np*sizeof(double));
  y=(double*)malloc(np*sizeof(double));
  
  srand48(Seed);

  ng=0;
  for (n=0;n<nsim;n++) {
    vx=0.0;
    vy=0.0;
    vz=0.0;
    for (i=0;i<np;i++) {
      x[i]=drand48();
      y[i]=drand48()*(1-x[i]);
      vx+=x[i];
      vy+=y[i];
      vz+=1-x[i]-y[i];
    }
    // if all of X, Y, and Z have scores greater than two, increment ng
    if (vx>2.0 && vy>2.0 && vz>2.0) {
//     fprintf(stdout,"%.3lf %.3lf %.3lf = %.3lf\n",
//	vx,vy,vz,vx+vy+vz);
     ng++;
    }
  }
  fprintf(stdout,"%.5lf\n",((double)ng)/nsim);
}
