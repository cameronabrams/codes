#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main ( int argc, char * argv[] ) {
  char ln[255];
  double s=0.0,ss=0.0,x=0.0;
  int n = 0;
  while (fgets(ln,255,stdin)) {
    if (ln[0]!='#') {
      n++;
      x=atof(ln);
      s+=x;
      ss+=x*x;
    }
  }
  printf("%i %.4lf %.4lf\n",n,s/n,sqrt((ss-s*s/n)/n));
}
