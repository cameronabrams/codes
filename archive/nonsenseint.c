#include <stdio.h>
#include <stdlib.h>
#include <float.h>

unsigned getbits ( unsigned long x , int p, int n ) 
{
  return (x >> (p+1-n)) & ~(~0 << n);
}
typedef unsigned int       uint32;
int main (int argc, char * argv[])
{

  uint32 a[2];
  double * pd;
  int i;

  double D=0.0;//(double)atof(argv[1]);
  unsigned long *d=(long*)&D;
  unsigned long L=0,*l=&L;

  sscanf(argv[1],"%u",&a[0]);
  sscanf(argv[2],"%u",&a[1]);
  memcpy(d,&a[0],8);
  printf("%s\n",memcmp(d,&a[0],8)?"bad":"ok");

  //  *l=1<<1;

  printf("%e(%u:%u)\n",D,sizeof(D),sizeof(a));
  printf("is d nan? %s\n",d!=d?"yes":"no");

  for (i=8*sizeof(double)-1;i>-1;i--) {
    printf("%s%s",getbits(*d,i,1)?"1":"0",!((i)%4)?" ":"");
  }
  printf("\n");

/*   a[0]=a[1]=0; */

/*   pd=(double*)&a[0]; */

/*   (*pd)=d; */

/*   printf("%e %u %u\n",(*pd),a[0],a[1]); */

/*   a[0]=atoi(argv[2]); a[1]=atoi(argv[3]); */

/*   printf("%8x %u %u\n",(*pd),a[0],a[1]); */

}
