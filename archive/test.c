#include <stdio.h>
#include <stdlib.h>

int main ( int argc, char * argv[] ) {
  fprintf(stdout,"float %d double %d int %d\n",
	  sizeof(float), sizeof(double), sizeof(int));
}
