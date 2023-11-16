/* stest -- just playing around */

#include <stdio.h>

typedef struct {
  double e1;
  double e2;
  int (*map)(double,double);
} myType;

void InitMyType ( myType * mt, double x, double y, 
		  int (*func)(double,double));

