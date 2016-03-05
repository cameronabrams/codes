#include <stdio.h>
#include <stdlib.h>

/** Code01  (c) 2004
 *
 *  Cameron Abrams
 *  Drexel University
 *  Department of Chemical Engineering
 * 
 *  Philadelphia
 *
 *  This code will compute the (x,y) positions of the nodes
 *  of a Harter-Heightway Dragon.  The user specifies
 *  a generation for the Dragon, and a single data file
 *  of Y vs X is generated, suitable for plotting.
 *
 *  This code serves as the first test of "download/compile/run"
 *  for a graduate class on molecular simulation.
 *  
 */

#define RIGHT 1
#define LEFT -1

typedef struct pointNode * ptp;
typedef struct pointNode
{
  int g;
  double x;
  double y;
  ptp next;
} pt;

ptp pointNode_Init(int g, double x, double y) {
  ptp p = malloc(sizeof(pt));
  if (p) {
    p->g=g;
    p->x = x;
    p->y = y;
    p->next=NULL;
    return p;
  }
  else
    return NULL;
}

void pointNode_Free (ptp P) {
  if (!P) return;
  free((ptp)(P->next));
  free((ptp)P);
}

void pts_out (FILE * fp, ptp P, int g) {
  ptp p=P;
  if (!fp||!P) return;
  for (p=P;p;p=p->next) {
      fprintf(fp, "%.10lf %.10lf %i\n", p->x, p->y, p->g);
  }
}

int main (int argc, char * argv[]) {
  ptp P = NULL;
  char fn[20];
  int n = atoi(argv[1]);	    /* number of generations */
  ptp p = NULL, q = NULL;
  int i = 0, d = RIGHT;
  double xm, ym, xn, yn;
  FILE * fp = NULL;
  
  P = pointNode_Init(0, -0.5, 0);
  P->next = pointNode_Init(0, 0.5, 0);  /* set up generation 0 */
  
  /* For each generation */
  for (i = 0; i < n; i++) {
    d = RIGHT;
    for (p = P; p->next; p = p->next->next) {
      q = p->next; /* hold the endpoint at q */
      /* find p-q segment midpoint */
      xm = (p->x + q->x) / 2.0;
      ym = (p->y + q->y) / 2.0;
      if (i == 0 || d == LEFT) {
	xn = p->y - ym + xm;	/* "left" displacement */
	yn = xm - p->x + ym;
      }
      else {
	xn = ym - p->y + xm;	/* "right" displacement */
	yn = p->x - xm + ym;
      }
      /* create the new point */
      p->next = pointNode_Init(i+1, xn, yn);
      /* insert this new point between p and q */
      p->next->next = q;
      d*=-1; /* flip the displacement direction */
    }

    /* generate output filename, open it, and output points */
    sprintf(fn,"%i.hhd",i);
    fp=fopen(fn,"w");
    pts_out(fp, P, i);
    fclose(fp);
  }
  pointNode_Free(P);
}
