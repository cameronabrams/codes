#include "point.h"

ptPtr ptPtr_Init (double x, double y, double z)
{
    ptPtr p = malloc(sizeof(pt));
    if (p)
    {
	p->x = x;
	p->y = y;
	p->z = z;
	p->chi = p->sib = NULL;
	p->g = 0;
	p->s = 0;
    }
    return p;
}

double ptPtr_abs (ptPtr p)
{
    if (!p) return 0.0;
    return sqrt(p->x*p->x + p->y*p->y + p->z*p->z);
}

ptPtr ptPtr_cross (ptPtr res, ptPtr a, ptPtr b)
{
    if (!a || !b) return res;
    if (!res) res = ptPtr_Init(0.0, 0.0, 0.0);
    res->x = a->y*b->z - b->y*a->z;
    res->y = a->z*b->x - b->z*a->x;
    res->z = a->x*b->y - b->x*a->y;
    return res;
}

void ptPtr_out(FILE * fp, ptPtr p)
{
    if (!fp || !p) return;
    fprintf(fp, "#%.5le, %.5le, %.5le\n", p->x, p->y, p->z);
}

void freeptPtr (ptPtr p)
{
    if (!p) return;
    freeptPtr(p->chi);
    freeptPtr(p->sib);
    free((ptPtr)p);
}
