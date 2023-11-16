/* xyz_sep.c
 *
 * For a file containing lines of xyz data, output all 
 * separations.
 */

#include <stdio.h>
#include <math.h>

typedef struct XYZPT
{
   double x,y,z;
} xyz;

xyz data[100000];

char scr_ln[255];
#define ANINT(a) (double) (long)  ( ((a) < 0.0) ? ((a) - 0.5) : ((a) + 0.5) )

void main (int argc, char **argv)
{
   FILE *fp = fopen(argv[1],"r");
   int i,j,N;
   int pbc=argc>2;
   double box=atof(argv[2]);
   xyz d;
   double D;

   if (!fp) exit(1);

   i=0;
   while (fgets(scr_ln,255,fp))
   {
     sscanf(scr_ln,"%lf %lf %lf", &(data[i].x), &(data[i].y), &(data[i].z));
     i++;
   }
   fprintf(stderr,"%i points\n",N=i);
   fprintf(stderr,"%i argc %i pbc\n",argc,pbc);
   if (pbc) fprintf(stderr,"%.5lf boxsize\n",box);
   for (i=0;i<N;i++){
     for (j=i+1;j<N;j++){
       d.x=data[i].x-data[j].x;
       d.y=data[i].y-data[j].y;
       d.z=data[i].z-data[j].z;
       if (pbc){
	 d.x-=box*ANINT(d.x/box);
	 d.y-=box*ANINT(d.y/box);
	 d.z-=box*ANINT(d.z/box);
       }
       D=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
       fprintf(stdout,"%i %i %.6lf\n",i+1,j+1,D);
     }
   }
}
