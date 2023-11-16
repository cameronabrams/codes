#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXDATA 100000

typedef struct THEPOINT {
  double x;
  double y;
  double z;
} pt;



int main ( int argc, char * argv [] ) {

  pt data[MAXDATA];
  int nData=0;
  double plattime[100],platvalue[100];
  int nPlat=0,nThisPlat=0;
  int plat= 0;
  double cutoff=2.0;
  char ln[255];
  FILE *fp=NULL;
  char *fn=NULL;
  int i;
  double diff, max=1.e-9;
  char *fmt=NULL;
  char prfmt[50];
  int trans;
  int intervals=0;
  double diff2,maxdiff2;
  int quiet=0;

  for (i=1;i<argc;i++) {
    if (argv[i][0]!='-') fn=argv[i];
    else if (!strcmp(argv[i],"-fmt")) fmt=argv[++i];
    else if (!strcmp(argv[i],"-c")) cutoff=atof(argv[++i]);
    else if (!strcmp(argv[i],"-i")) intervals=1;
    else if (!strcmp(argv[i],"-q")) quiet=1;
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

  diff=0.0;
  for (i=0;i<nData;i++) {
    if (i==0) diff=(data[i].y-data[i+1].y)/(data[i].x-data[i+1].x);
    else if (i==(nData-1))
      diff=(data[i-1].y-data[i].y)/(data[i-1].x-data[i].x);
    else diff=(data[i-1].y-data[i+1].y)/(data[i-1].x-data[i+1].x);
    data[i].z=diff;
    diff2=diff*diff;
    maxdiff2+=diff2;
  }
  maxdiff2/=nData;
  diff=cutoff*sqrt(maxdiff2);
  if (!quiet) fprintf(stderr,"# using diff=%.5lf\n",diff);

  plat=0;
  plattime[plat]=0.0;
  trans=0;
  nThisPlat=0;
  plat++;
  for (i=0;i<nData;i++) {
    if (data[i].z<-diff) trans=1;
    else {
      platvalue[plat]+=data[i].y;
      nThisPlat++;
      if (trans) { // end of a plateau 
	plattime[plat]=data[i].x;
	platvalue[plat]/=nThisPlat;
	trans=0;
	plat++;
	nThisPlat=0;
      }
    }
  }
  nPlat=plat;

  // fprintf(stderr,"%i\n",intervals);
  for (i=1;i<nPlat;i++) {
    fprintf(stdout,"%.5lf %.5lf\n",plattime[i]-plattime[i-1],platvalue[i]);
  }

  

}
