#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXDATA 100000

typedef struct THEPOINT {
  double x;
  double y;
} pt;



int main ( int argc, char * argv [] ) {

  pt data[MAXDATA];
  int nData=0;
  double plattime[100];
  int nPlat=0;
  int plat= 0;
  char ln[255];
  FILE *fp=NULL;
  char *fn=NULL;
  int i;
  double diff, max=1.e-9;
  char *fmt=NULL;
  char prfmt[50];
  int trans;
  
  for (i=1;i<argc;i++) {
    if (argv[i][0]!='-') fn=argv[i];
    else if (!strcmp(argv[i],"-fmt")) fmt=argv[++i];
  }

  if (!fmt) strcpy(prfmt,"%.5lf %.5lf\n");
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
  nData=i;
  
  plat=0;
  plattime[plat]=0.0;
  trans=0;
  for (i=0;i<nData;i++) {
    if (data[i].y>2.0) trans=1;
    else {
      if (trans) { // end of a plateau 
	plat++;
	plattime[plat]=data[i].x-plattime[plat-1];
      }
      trans=0;
    }
  }
  nPlat=plat;


  for (i=0;i<nPlat;i++) {
    fprintf(stdout,"%.5lf\n",plattime[i]);
  }

  

}
