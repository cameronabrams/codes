#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct VEC {
  double x, y, z;
} vec;


double diffnorm ( vec * a, vec * b ) {
   double dx=a->x-b->x;
   double dy=a->y-b->y;
   double dz=a->z-b->z;
   return sqrt(dx*dx+dy*dy+dz*dz);
} 

void random_step ( vec * L, int i, double R ) {
   double phi, theta, cphi, sphi, ctheta, stheta;
   phi=M_PI*(1-2*((double)rand())/RAND_MAX);
   cphi=cos(phi);
   sphi=sin(phi);
   ctheta=1-2*((double)rand())/RAND_MAX;
   theta=acos(ctheta);
   stheta=sin(stheta);
   L[i].x=L[i-1].x+R*stheta*cphi;
   L[i].y=L[i-1].y+R*stheta*sphi;
   L[i].z=L[i-1].z+R*cphi;
}

int main (int argc, char * argv[]) {
  FILE * fp=NULL;
  int nsteps = 10000;
  vec * L=NULL;
  double * dist;
  int i, w;
  int nwalkers = 100;

  L=(vec*)malloc(nsteps*sizeof(vec));
  dist=(double*)malloc(nsteps*sizeof(vec));

  for (i=1;i<nsteps;i++) {
    dist[i]=0.0;
  }

  L[0].x=L[0].y=L[0].z=0.0;
  for (w=0;w<nwalkers;w++) {
    fprintf(stdout,"Walker %d...\n",w);
    for (i=1;i<nsteps;i++) {
      random_step(L,i,1.0);
      dist[i]+=diffnorm(&L[i],&L[0]);
    }
  }

  for (i=1;i<nsteps;i++) {
    dist[i]/=nwalkers;
  }

  fp=fopen("out.dat","w");
  for (i=0;i<nsteps;i++) {
    fprintf(fp,"%i %.5lf\n",i,dist[i]);
  }
  fclose(fp);

}
