/* ptest.c -- this program will test two methods of accessing
 * array elements.  The test is the following:  given an array
 * of N double-precision numbers, compute the average value.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ADIM 1000000

double A[ADIM];

typedef struct ARRAY_ELEM * aep;

typedef struct ARRAY_ELEM 
{
  double x;
  aep next;
} _ae_ ;

_ae_ L[ADIM];

void main (int argc, char **argv)
{
  int i=0;
  int N=atoi(argv[1]);
  aep p=NULL;
  double sum=0.0;
  
  if (N>ADIM) {printf("Whoops! N %i is too big.\n", N);exit(0);}

  /* initialize the overlay array/fill the array */
  for (i=0;i<N;i++)
  {
    L[i].next=i<N-1?&(L[i+1]):NULL;
    A[i]=((double)rand()/RAND_MAX);
    L[i].x=A[i];
  }

  /* compute the average */
  sum=0.0;
#ifdef LL
  for (p=L;p;p=p->next) sum+=p->x;
#else
  for (i=0;i<N;i++) sum+=A[i];
#endif
  sum/=N;

  printf("%.5lf\n",sum);

}
