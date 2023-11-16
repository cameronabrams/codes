#ifndef POINT_H
#define POINT_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct POINT * ptPtr;
typedef struct POINT
{
    double x;
    double y;
    double z;
    int g, s;
    ptPtr chi, sib;
} pt;


ptPtr ptPtr_Init (double x, double y, double z);
double ptPtr_abs (ptPtr p);
ptPtr ptPtr_cross (ptPtr res, ptPtr a, ptPtr b);
void ptPtr_out(FILE * fp, ptPtr p);
void freeptPtr (ptPtr p);

#endif

