#include <stdio.h>
#include <stdlib.h>

int main (int argc, char * argv[])
{
  int ny, nx;
  float dy, dx;
  int i,j;
  float x,y;

  ny=atoi(argv[1]);
  nx=atoi(argv[2]);
  dy=atof(argv[3]);
  dx=atof(argv[4]);

  //printf("%i %i %f %f\n",ny,nx,dy,dx);

  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      printf("%.2f %.2f\n", (j+1)*dx, (i+1)*dy);
    }
  }

}
