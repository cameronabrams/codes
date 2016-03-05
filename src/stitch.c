#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int readfile ( FILE * fp, double * X, double * FX ) {
    int i;
    char ln[255];
    while (fgets(ln,255,fp)) {
      if (ln[0]!='#') {
        sscanf(ln,"%lf %lf",&X[i],&FX[i]);
        i++;
      }
    }
    return i;
}

int main ( int argc, char * argv [] ) {

    double * x1, * fx1;
    double * x2, * fx2;
    double Y;
    
    int N1=10000;
    int N2=1000;
    int i;
    
    FILE * fp;
    
    x1=(double*)malloc(N1*sizeof(double));
    x2=(double*)malloc(N2*sizeof(double));
    fx1=(double*)malloc(N1*sizeof(double));
    fx2=(double*)malloc(N2*sizeof(double));
    
    fp=fopen(argv[1],"r");
    N1=readfile(fp,x1,fx1);
    fclose(fp);

    fp=fopen(argv[2],"r");
    N2=readfile(fp,x2,fx2);
    fclose(fp);
    
    // assuming same size
    Y=fx1[N1-1]/fx2[0];
    
    for (i=0;i<N2;i++) fx2[i]*=Y;
    
    for (i=0;i<(N1-1);i++) fprintf(stdout,"%.5lf %.5le\n",x1[i],fx1[i]);
    for (i=0;i<N2;i++) fprintf(stdout,"%.5lf %.5le\n",x2[i],fx2[i]);
}
