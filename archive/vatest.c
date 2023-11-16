#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

void theFunct ( char * a, ... )
{
  va_list ap;
  int i,j=-1;
  va_start(ap,a);

  i = va_arg(ap,int);
  j = va_arg(ap,int);
  fprintf(stdout,"%s %i %i\n",a,i,j);

  va_end(ap);
}

void main ( int argc, char * argv[] )
{
  int a=0,b=0,c=0;
  void (*func)( char * a, ... );
  
  sscanf(argv[1],"%i,%i,%i",&a,&b,&c);
  printf("%i %i %i\n",a,b,c);
  func=theFunct;

  (*func)("Here is the first:",a,b);

}
