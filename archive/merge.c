/*
 * ydist.c
 *
 * (c) 1999-2009 cfa
 * 
 * ydist.c computes a histogram and accompanying discrete distribution
 * from a column of numbers.  The input may contain two columns: the
 * first is the data to be binned, denoted by x_i, and the second (if
 * present) is taken as a single value of some function f(x_i). If
 * present, the functions <f(x)>, <del-f(x)> (std.dev.), <min(f(x))>,
 * and <max(f(x))> are computed from this data.
 * 
 * There is no restriction to the number of input data lines.
 * The maximum number of histogram bins is set as 5000, but may
 * be changed at compile time with the flag -DMAXBINS=100000, for 
 * example.
 *
 * usage:
 * 
 * ydist [-bs #] [-ymin #] [-ymax #] [-f #] [-n #] [autoscale #] [-| |file]
 * 
 * -ymin    allows user to specify minimum bound on histogram;
 *	    if not supplied by user, it is **assumed to be 0.0**.
 * -ymax    allows user to specify maximum bound on histogram;
 *	    if not supplied by user, it is set to the maximum value of
 *	    of the data read in.
 * -bs	    allows user to specify bin size of histogram; default is 1.0.
 *          future: have the code automatically compute an optimal binsize
 *          based on the min,max,stddev,mean of the input data?
 * -n	    allows user to specify the *histogram* normalization constant, 
 *          if it is different from the count of data points.  The normalized
 *          histogram is converted into a discrete distribution function by
 *          dividing by the bin size.
 * -f	    bin growth factor, for binning logarithmic data; default is 1.0
 *	    (no growth).  Logarithmic data is *not* converted into a 
 *          distribution function.
 * -autoscale turns on autoscaling; requires an argument which is the
 *          number of bins desired
 * - (or nothing)
 *          tells ydist.c that the input is a column of numbers from stdin, 
 *          terminated by ^D or EOF.
 * (filename) 
 *          tells ydist.c to open and read a single column from the file
 *          named 'filename'.
 * 
 * compile this as 'cc -o ydist ydist.c -lm'
 * 
 * 15Apr1999 Berkeley, California
 * 07Jul2000 Mainz, Germany
 * 13Mar2009 Philadelphia, PA
 *
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define MAXLINE 255
char scr_line_[MAXLINE];
#ifndef MAXBINS
#define MAXBINS 5000
#endif
int bins[MAXBINS];

#ifndef DOUBLE_PRECISION
typedef float _fl_;
#else
typedef double _fl_;
#endif

_fl_ hist[MAXBINS], secsum[MAXBINS], secmin[MAXBINS], secmax[MAXBINS], secssum[MAXBINS];

typedef struct DATA {
  _fl_ y;
  _fl_ y2;
} data_t;

#ifndef MAXDATA
#define MAXDATA 100000000
#endif

void usage (void)
{
  printf("usage:\n");
  printf("ydist {- | <fileName>} (options) \n");
  printf("Options:\n");
  printf("\t-ymin # -ymax #\n");
  printf("\t-n <norm const(calculated)>\n");
  printf("\t-bs <bin size (1.0)>\n");
  printf("\t-f <bin size growth factor (1.0=const. bin size)>\n");
  exit(-1);
}

int bin (_fl_ x, _fl_ bs, _fl_ Ymin, _fl_ f)
{
  int i=0;
  _fl_ sum=Ymin;
    
  while (i<MAXBINS&&x>sum) sum+=bs*pow(f,i++);
  if (i==MAXBINS) {fprintf(stderr,"error bigbin y=%.5lf\n",x);exit(0);}
  
  return i;
}

_fl_ binv (int b, _fl_ bs, _fl_ Ymin, _fl_ f)
{
  int i=0;
  _fl_ sum=0.0;
  
  for (i=0;i<b;i++) sum+=bs*pow(f,i);
  sum-=0.5*bs*pow(f,--i);
  return sum+Ymin;
}

char y1str[50];
void main (int argc, char * argv[])
{
  int i=0, b=0;
  char * fn=NULL;
  char * p=NULL;
  _fl_ sum_a=0.0, sumxx_a=0.0;
  _fl_ y1=0.0, y2=0.0;
  _fl_ totaly=0.0;
  _fl_ ymin=1.e9, ymax=-1.e9, bval=0.0;
  _fl_ Ymin=0.0, Ymax=-1.e9;
  int nData=0;
  _fl_ binsize=1.0;
  _fl_ f=1.0;
  int nbins=10;
  FILE * fp=stdin;
  data_t * data;
  int autoscale=0;
  int q=1;
  int autoscale_pad=0;
  data = (data_t*)calloc(MAXDATA,sizeof(data_t));

  for (i=0;i<MAXBINS;i++) {
    bins[i]=0;
    hist[i]=secsum[i]=secssum[i]=0.0;
    secmax[i]=-1.e10;
    secmin[i]=1.e10;
  }
  
  for (i=1;i<argc;i++) {
    if (argv[i][0] != '-') fn=argv[i];
    else if (!strcmp(argv[i], "-q")) q=1;
    else if (!strcmp(argv[i], "-")) fp=stdin; /* not necessary */
    else if (!strcmp(argv[i], "-bs")) binsize=atof(argv[++i]);
    else if (!strcmp(argv[i], "-f")) f=atof(argv[++i]);
    else if (!strcmp(argv[i], "-ymin")) Ymin=atof(argv[++i]);
    else if (!strcmp(argv[i], "-ymax")) Ymax=atof(argv[++i]);
    else if (!strcmp(argv[i], "-n")) totaly=atof(argv[++i]);
    else if (!strcmp(argv[i], "-autoscale")) {autoscale=1; nbins=atoi(argv[++i]);}
    else if (!strcmp(argv[i], "-autoscale_pad")) {autoscale_pad=atoi(argv[++i]);}
    else usage();
  }

  if (!q) printf("# YDIST cfa 2009\n"); 

  if (fn) {
    fp=fopen(fn, "r");
    if (!fp) {printf("Error. Could not open %s.\n", fn);exit(-1);}
  }
    
  nData=0;
  while (fgets(scr_line_,MAXLINE,fp)) {
    p=scr_line_;
    if (isspace(*p)) while(isspace(*(++p)));
    /* read the first number on the line */
    sscanf(p,"%s",y1str);
    y1=atof(y1str);
    sum_a+=y1;
    sumxx_a+=y1*y1;
    ymax=(y1>ymax?y1:ymax);
    ymin=(y1<ymin?y1:ymin);
    /* advance the pointer until we hit a second number or EOL */
    p+=strlen(y1str);
    while(isspace(*(++p)));
    y2=0.0;
    if (p) y2=atof(p);
    data[nData].y=y1;
    data[nData].y2=y2;
    nData++;
    //fprintf(stderr,"%.5le\n",y2);
  }

  if (autoscale) {
    binsize = (ymax-ymin)/nbins;
    Ymin=ymin-autoscale_pad*2*binsize;
    Ymax=ymax+autoscale_pad*2*binsize;
    if (!q) printf("# autoscale count %i %.5lf %.5lf %.5lf %i\n",nData,ymin,ymax,binsize,nData);
  }

  for (i=0;i<nData;i++) {
    /* determine which bin to increment */
    //    fprintf(stderr,"y[%i] %.5lf\n",i,data[i].y);
    b=bin(data[i].y,binsize,Ymin,f);
    bins[b]++;
    secsum[b]+=data[i].y2;
    secssum[b]+=data[i].y2*data[i].y2;
    if (data[i].y2>secmax[b]) secmax[b]=data[i].y2;
    if (data[i].y2<secmin[b]) secmin[b]=data[i].y2;
  }
  if (fp!=stdin) fclose(fp);
  if (!totaly) totaly=nData;
  if (Ymax==-1.e9) Ymax=(_fl_)((int)ymax + 1);
  printf("# mean: %.5lf",sum_a/nData);
  printf(" sd: %.5lf", sqrt((sumxx_a-sum_a*sum_a/nData)/nData));
  printf(" min: %.5lf max :%