#include <stdio.h>
#include <math.h>
/*
 * linint.c -- simple linear interpolator.
 *
 * given exactly 2 points (x_0,y_0) and (x_1,y_1), and
 * either a value for x *or* y, interpolates to
 * find the other.
 *
 * the points are read in from two lines of stdin.
 * output is to stdout. 
 * 
 * (c) 2001 cam abrams
 */
typedef struct POINT
{
    float x, y;
} pt;

char ln[255];
int main (int argc, char * argv[])
{
  FILE * fp;
  int i,nData;
  char *fmt=NULL;
  char PFMT[50],*pfmt=PFMT;
  float x=0,y=0,m,b;
  int xs=0,ys=0;
  pt d[2];
    
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-fmt")) fmt=argv[++i];
    else if (argv[i][0]=='x') {x=atof(&argv[i][2]);xs=1;}
    else if (argv[i][0]=='y') {y=atof(&argv[i][2]);ys=1;}
  }

  if (!xs&&!ys) {
    fprintf(stderr,"error: must specify either an x or y\n");
    exit(-1);
  }

  if (!fmt) strcpy(pfmt,"%.5lf %.5lf\n");
  else sprintf(pfmt,"%s %s\n", fmt, fmt);

  fp=stdin;
    
  i=0;
  while (fgets(ln, 255, fp)) {
    if (ln[0]!='#') {
      sscanf(ln, "%f %f", &(d[i].x), &(d[i].y));
      i++;
    }
  }
  nData=i;
  if (nData != 2) {
    fprintf(stderr,"error: need exactly 2 points\n");
    exit(-1);
  }
    

  m=(d[1].y-d[0].y)/(d[1].x-d[0].x);
  b=d[1].y-m*d[1].x;

  if (xs) fprintf(stdout,pfmt,x,m*x+b);
  else if (ys) fprintf(stdout,pfmt,(y-b)/m,y);
}
