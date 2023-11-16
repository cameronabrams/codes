#include <stdio.h>
#include <string.h>
#include <math.h>

#define sqr(a) ((a)*(a))
#define cub(a) ((a)*(a)*(a))

double vsph ( double R )
{
  return 4./3.*M_PI*cub(R);
}

double vseg ( double R, double g ) 
{
  if (g>R) return -1;
  return 2./3.*M_PI*cub(R)*(1.-g/R-0.5*g*(sqr(R)-sqr(g))/cub(R));
}

double voverlap ( double R1, double R2, double l )
{
  double g1,g2;

  if (l>(R1+R2)) return 0.0;

  if (l<R1||l<R2) return -1;

  g1=(sqr(R1)-sqr(R2)+sqr(l))/(2*l);
  g2=l-g1;

  return vseg(R1,g1)+vseg(R2,g2);

}

void usage ( char * huh )
{
  fprintf(stderr,"`%s'? read the source\n",huh);
}

int main ( int argc, char * argv[])
{
  int i=0;
  double R1=0, R2=0, l=0;

  double f=1.0;

  for (i=1;i<argc;i++)
  {
    if (!strcmp(argv[i],"-r1")) R1=atof(argv[++i]);
    else if (!strcmp(argv[i],"-r2")) R2=atof(argv[++i]);
    else if (!strcmp(argv[i],"-l")) l=atof(argv[++i]);
    else usage(argv[i]);
  }

/*   fprintf(stdout,"Sphere 1: %.5lf  Sphere 2: %.5lf\n", */
/* 	  vsph(R1),vsph(R2)); */

/*   fprintf(stdout,"Sum: %.5lf  Overlap: %.5lf  Total: %.5lf\n", */
/* 	  vsph(R1)+vsph(R2),voverlap(R1,R2,l), */
/* 	  vsph(R1)+vsph(R2)-voverlap(R1,R2,l)); */


  for (f=1.0;f<2.0;f+=0.01) {
    fprintf(stdout,"%.3f %.3f\n",
	    f,77.08*cub(f)+2*voverlap(f*R1,f*R2,l));
  }

}
