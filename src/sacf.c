/* 
   Scalar autocorrelation function calculation, straightforward
   algorithm.  Reads in an arbitrary number of time-series data,
   one item per line.

   Cameron F. Abrams, cfa22@drexel.edu

   compile using "gcc -o sacf sacf.c -lm"

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2004-2018
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Prints usage information */
void usage ( void ) {
  fprintf(stdout,"sacf usage:\n");
  fprintf(stdout,"sacf [options]\n\n");
  fprintf(stdout,"Options:\n");
  fprintf(stdout,"\t -dt [real(0.001)]\t\ttime increment between data\n");
  fprintf(stdout,"\t -h           \t\tPrint this info.\n");
}

/* Reads in scalar time series */
int x_in (FILE * fp, double * x, int * N) {
  int i;
  char ln[255];
  *N=0;
  while (fgets(ln,255,fp)) {  
    sscanf(ln,"%lf",&x[(*N)]);
    (*N)++;
  }
  return 0;
}

int main ( int argc, char * argv[] ) {

  double * x, * ax, ssum=0.0;
  int ndata=0, NDATA=100000;
  int t, dt, * cnt, tmax=0;
  double time_step = 0.001;
  int start=0,stop=0,step=1,ngr=0;
  int i;

  char * fnf;
  FILE * fp = stdin;

  char * fn = NULL;

  /* Here we parse the command line arguments;  If
   you add an option, document it in the usage() function! */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-dt")) dt=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-tmax")) tmax=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-i")) fn=argv[++i];
    else if (!strcmp(argv[i],"-h")) {
      usage(); exit(0);
    }
    else {
      fprintf(stderr,"Error: Command-line argument '%s' not recognized.\n",
	      argv[i]);
      exit(-1);
    }
  }
  
  x=(double*)calloc(NDATA,sizeof(double));
  ax=(double*)calloc(NDATA,sizeof(double));
  cnt=(int*)calloc(NDATA,sizeof(int));

  if (fn) fp=fopen(fn,"r");
  x_in(fp,x,&ndata);
  if (fn) fclose(fp);

  if (tmax>ndata) {
    fprintf(stderr,"ERROR: not enough data for tmax %i\n",tmax);
    exit(-1);
  }

  /* Compute autocorrelation using
     the straightforward algorithm */
  fprintf(stderr,"# computing...\n");fflush(stderr);
  for (t=0;t<ndata;t++) {
    ssum+=x[t]*x[t];
    for (dt=0;(t+dt)<tmax;dt++) {
      cnt[dt]++;  /* number of origins for interval length dt  */
      ax[dt] += x[t+dt]*x[t];
    }
  }
  for (t=0;t<tmax;t++) {
    ax[t] /= cnt[t]?(cnt[t]):1;
    fprintf(stdout,"%.5lf %.8lf\n",
	    t*time_step,ax[t]/ssum*ndata);
  }

}
