/* stest -- just playing around */

#include <stdio.h>
#include "ss.h"

void InitMyType ( myType * mt, double x, double y, 
		  int (*func)(double,double)) {
  mt->e1=x;
  mt->e2=y;
  mt->map=func;
}

