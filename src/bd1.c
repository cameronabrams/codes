/* 1-D Brownian dynamics with Umbrella Sampling
   
   D R E X E L   U N I V E R S I T Y
   Department of Chemical and Biological Engineering
   CHE 614 -- Thermodynamics II
   C. F. Abrams -- Winter 1516

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void f ( double x, double a, double b, double c, double * e, double * f ) {
  double x2=x*x;
  double x3=x2*x;
  double x4=x2*x2;
  double bb=b-2*a;
  *e = a*x4+bb*x2 + c*x + a;
  *f = -(4*a*x3+2*bb*x + c);
}

// compute the integral of the Boltzmann factor of f(x)
// T is dimensionless temperature
double Q ( double T, double xmin, double xmax, double a, double b, double c ) {
  int N=1000; // default number of integration points
  double dx = (xmax-xmin)/N;
  int i;
  double x=xmin;
  double e, dum;
  double s=0;
  for (i=0;i<N;i++) {
    f(x,a,b,c,&e,&dum);
    s+=((i==0||i==(N-1))?1:2)*exp(-e/T);
    x+=dx;
  }
  return dx/2.0 * s;
}

int main  ( int argc, char * argv [] ) {
  double x;
  double x0 = -10.;
  double xt, dx;
  double h1, h2;

  double xlo= -3;
  double xhi = -2;

  double dt = 0.0001;

  int I; // interval
  double dd=1.0;//0.25;
  int nI;// number of intervals of size dd

  int i;
  int nsteps = 80000000;
  double gam = 5.0;
  double T = 1.0;
  double force, e;
  //  double a = 0.01;
  //  double b = -0.25;
  //  double c = 0.3;
  double a = 0.02;
  double b = -1.0;
  double c = 0.0;

  int seed = 1243;

  // histogram
  int * h;
  int N=25; // size of histogram
  double xmin = -8;
  double xmax = 8;
  double bs;// = (xmax-xmin)/N;
  int bin;
  int minbin,maxbin;
  int n=0;
  double q;
  double buf;
  FILE * fp;
  char fn[20];

  srand(seed);
  
  h=(int*)malloc(N*sizeof(int));

  nI = (int)((xmax-xmin)/dd);

  buf = dd/(N-1);
  h1=dt/gam;
  h2=sqrt(6*T*h1);

  for (I=0;I<nI;I++) {
    minbin=9999;
    maxbin=-9999;
    xlo = xmin + I*dd;
    xhi = xlo + dd + buf;
    bs = (xhi-xlo)/N;
    for (i=0;i<N;i++) h[i]=0;
    n=0;
    x0=0.5*(xhi+xlo);
    for (i=1;i<=nsteps;i++) {
      f(x,a,b,c,&e,&force);
      dx = h1*force+h2*(1-2*((double)rand())/RAND_MAX);
      xt = dx + x0;
      if (xt > xhi || xt < xlo) x = x0;
      else x = xt;
      x0 = x;
      bin = (int)((x-xlo)/bs);
      //fprintf(stdout,"%i\n",bin);
      if (bin<minbin) minbin=bin;
      else if (bin>maxbin) maxbin=bin;
      h[bin]++;
      n++;
    }
    sprintf(fn,"d%i",I);
    fp=fopen(fn,"w");
    fprintf(fp,"#x  rho(x)  boltz(x)\n");
    for (i=0;i<N;i++) {
      x=xlo+i*bs;
      f(x,a,b,c,&e,&force);
      fprintf(fp,"%.5lf %.5le %.5le %.5le\n",xlo+i*bs,((double)h[i])/(n*bs),e);
    }
    fclose(fp);
    fprintf(stderr,"# interval [%.3lf,%.3lf] [%i,%i] complete.  Histogram in %s\n",xlo,xhi,minbin,maxbin,fn);
  }


}
