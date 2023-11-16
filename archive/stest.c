/* stest -- just playing around */

/* lesson learned:  if an instance of dt A is declared in
   b.c, while the function to initialize it is defined in
   a.c, b.c calls this init function.  if a.c passes in the
   parameter list the address of a static global in a.c... */

#include <stdio.h>
#include "ss.h"

static int K_ = -3;

static int testFunction1 ( double x, double y )
{
  int i=-K_;
  if (x>y) i=K_;
  return i;
}

int main (int argc, char * argv[])
{
   myType A1, * a=&A1;
   
   InitMyType(a,atof(argv[1]),atof(argv[2]),testFunction1);

   printf("%i\n",a->map(a->e1,a->e2));

}
