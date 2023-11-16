#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*
 * derv.c -- uses 2nd order finite difference to differentiate a function
 * supplied in the discrete form x_i, f(x_i).
 * 
 * Two columns of input data can be read from stdin (no arguments)
 * or from an input file (single argument).  Output is 2 columns, 
 * [x_i, f'(x_i)], to stdout.
 * MAXPTS can be redefined (below) if more than 5000 points are 
 * being considered.
 * 
 * (c) 1999-2010 cameron abrams
 */
typedef struct POINT
{
    double x, y;
} pt;


void boxcar_data (pt * data, int n, int b) {
  int i,j,nn;
  double x;
  pt * y=(pt*)malloc(n*sizeof(pt));
  for (i=0;i<n;i++) {
    x=0.0;
    nn=0;
    for (j=i-b/2;j<=i+b/2;j++) {
      if (j>-1&&j<n) {
	x+=data[j].y;
	nn++;
      }
    }
    x/=nn;
    printf("%.5le %.5le %i\n",data[i].y,x,nn);
    y[i].y=x;
  }
  for (i=0;i<n;i++) {
    data[i].y=y[i].y;
  }
  free(y);
} 


#ifndef MAXPTS
#define MAXPTS 5000
#endif
pt data[MAXPTS];
int nData;
char ln[255];
int main (int argc, char * argv[])
{
  FILE *fp=NULL;
  char *fn=NULL;
  int i;
  double diff, max=1.e-9;
  char *fmt=NULL;
  char prfmt[50];
  int boxsize=0;
  double lastdiff;
  int diffonly=0;

  for (i=1;i<argc;i++) {
    if (argv[i][0]!='-') fn=argv[i];
    else if (!strcmp(argv[i],"-fmt")) fmt=argv[++i];
    else if (!strcmp(argv[i],"-box")) boxsize=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-d")) diffonly=1;
  }

  if (!fmt) strcpy(prfmt,"% 10.5lf% 10.5lf\n");
  else sprintf(prfmt,"%s %s\n", fmt, fmt);
  
  if (!fn) fp=stdin;
  else fp=fopen(fn, "r");
    
  i=0;
  while (fgets(ln, 255, fp)) {
    if (ln[0]!='#') { 
      sscanf(ln, "%lf %lf", &(data[i].x), &(data[i].y));
      max=(data[i].y>max?data[i].y:max);
      i++;
    }
  }
  if (fp!=stdin) fclose(fp);
  nData=i;

  if (boxsize) boxcar_data(data,nData,boxsize);
    
  diff=0.0;
  for (i=0;i<nData;i++) {
    if (diffonly) {
      if (i) {
	printf("% 10.5lf\n",data[i].y-data[i-1].y);
      }
    }
    else {
      if (i==0) diff=(data[i].y-data[i+1].y)/(data[i].x-data[i+1].x);
      else if (i==(nData-1))
	diff=(data[i-1].y-data[i].y)/(data[i-1].x-data[i].x);
      else diff=(data[i-1].y-data[i+1].y)/(data[i-1].x-data[i+1].x);
      if (i) {
	if (diff>0.0&&lastdiff<0.0) {
	  fprintf(stdout,"# local min at %.4le %.4le\n",data[i].x,data[i].y);
	}
	else if (diff<0.0&&lastdiff>0.0) {
	  fprintf(stdout,"# local max at %.4le %.4le\n",data[i].x,data[i].y);
	}
      }
      lastdiff=diff;
      printf(prfmt, data[i].x, diff);
    }
  }
    
  fprintf(stderr,"# nData = %i\n", nData);
    

}
