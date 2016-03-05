#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXLINE 255
#define MAXDATA 100000

double p2 ( double x) {
  return 0.5*(3*x*x-1);
}

double dot ( double x1, double y1, double z1, double x2, double y2, double z2 ) {
  return x1*x2+y1*y2+z1*z2;
}

int main ( int argc, char * argv [] ) {
  char * fn=NULL;
  FILE * fp = stdin;
  char scr_line_[MAXLINE];
  int nData =0;
  double t[MAXDATA], x[MAXDATA], y[MAXDATA], z[MAXDATA];
  double runp2[MAXDATA], runp2sq[MAXDATA];
  double this_p2;
  int nins[MAXDATA];
  int i,j;

  if (argc>1) fn=argv[1];
  if (fn) {
    fp=fopen(fn, "r");
    if (!fp) {printf("Error. Could not open %s.\n", fn);exit(-1);}
  }

  nData=0;
  while (fgets(scr_line_,MAXLINE,fp)) {
    sscanf(scr_line_,"%lf %lf %lf %lf\n",&t[nData],&x[nData],&y[nData],&z[nData]);
    //    fprintf(stderr,"%lf %lf %lf %lf\n",t[nData],x[nData],y[nData],z[nData]);
    nData++;
  }
  printf("# nData: %i\n",nData);
  fclose(fp);

  for (i=0;i<nData;i++) {
    runp2[i]=0.0;
    runp2sq[i]=0.0;
    nins[i]=0;
  }

  for (i=0;i<nData;i++) {
    for (j=i;j<nData;j++) {
      this_p2=p2(dot(x[i],y[i],z[i],x[j],y[j],z[j])/(sqrt(dot(x[i],y[i],z[i],x[i],y[i],z[i]))*sqrt(dot(x[j],y[j],z[j],x[j],y[j],z[j]))));
      runp2[j-i]+=this_p2;
      runp2sq[j-i]+=(this_p2*this_p2);
      nins[j-i]++;
    }
  }
  fp=fopen("p2corr.dat", "w");
  for (i=0;i<nData;i++) {
    fprintf(fp,"%.5lf %.5lf %.5lf %i\n",t[i],0.4*runp2[i]/nins[i],0.4*sqrt((runp2sq[i]-runp2[i]*runp2[i]/nins[i])/nins[i]),nins[i]);
  }
  fclose(fp);


}
