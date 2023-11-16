#include <stdio.h>
void main (int argc, char * argv[])
{
  int n = atoi(argv[1]);

  if ( n & 1 ) printf("odd\n");

}
