/* sliding average */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXDATA 10000
char ln[255];
int main (int argc, char * argv[]) {
  double x[MAXDATA];
  int i;
  int k;
  int n;
  int w=5;
  int nData;
  FILE *fp;
  char * fn=NULL;
  double avg;

  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-w")) w=atoi(argv[++i]);
    else fn=argv[i];
  }
  //fprintf(stderr,"%i %x\n",w,fn);
  if (fn) fp=fopen(fn,"r");
  else fp=stdin;
  i=0;
  while (fgets(ln,255,fp)) {
    sscanf(ln,"%lf",&x[i]);
    i++;
  }
  nData = i;
  if (fp!=stdin) fclose(fp);
 
  for (i=0;i<nData;i++) {
    avg=0;
    n=0;
    for (k=-w;k<=w;k++) {
      if ((i+k)>=0&&(i+k)<nData) {
	avg+=x[i+k];
	n++;
      }
    }
    printf("%.5lf\n",avg/=n);
  }
}
