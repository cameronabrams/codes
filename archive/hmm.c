#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "bits.h"


int main ( int argc, char * argv[] )
{
   uint p;
   int anIndex,aSegment;

   int i=0,s;

   p=atoi(argv[1]);

   printf("0x%8x at 0x%8x %s\n",p,&p,sprint_as_bits(&p,sizeof(p)));
   //flipbit(&p,atoi(argv[2]),8*sizeof(p));
   //   printf("0x%8x %s\n",p,sprint_as_bits(&p,8*sizeof(p)));

   exit(0);


}
