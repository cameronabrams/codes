/*
 * 2ddist.c
 *
 * (c) 1999-2016 cfa
 * 
 * 2ddist.c computes a 2-dimensional histogram h(x1,x2), from a 
 * two-column file of sample data, called x1 and x2.  
 * 
 * usage:
 * 
 * 2ddist [-r1 x1min,x1max,x1binsize] [-r2 ...] [-|file]
 * 
 * 
 * 15Apr1999 Berkeley, California
 * 07Jul2000 Mainz, Germany
 * 27Feb2004-2016 Philadelphia, Pennsylvania
 *
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define MAXLINE 255
char scr_line_[MAXLINE];
#ifndef MAXBINS
#define MAXBINS 5001
#endif

#ifndef MAXDATA
#define MAXDATA 80000000
#endif

#define DOUBLE_PRECISION

#ifndef DOUBLE_PRECISION
typedef float _fl_;
#else
typedef double _fl_;
#endif

int bins[MAXBINS][MAXBINS];
_fl_ vals[MAXBINS][MAXBINS];

void usage (void)
{
  printf("usage:\n");
  printf("2ddist {- | <fileName>} (options) \n");
  printf("Options:\n");
  printf("\t-r1 min,max,incr -r2 min,max,incr\n");
  printf("\t-wt use weights that appear in third column of input\n");
  printf("\t-kT kBT (if included, free energy is also output)\n");
  exit(-1);
}

int bin ( _fl_ x, _fl_ xmin, _fl_ xbs ) 
{
  return (int)((x-xmin)/xbs);
}

void two_bin (int *i, int *j, _fl_ x1, _fl_ x2, _fl_ x1min, _fl_ x2min, _fl_ x1bs, _fl_ x2bs)
{
  *i=bin(x1,x1min,x1bs);
  *j=bin(x2,x2min,x2bs);
}

_fl_ binv (int b, _fl_ xbs, _fl_ xmin)
{
  return b*xbs+xmin;
}

typedef struct DATA {
  _fl_ x1;
  _fl_ x2;
  _fl_ x3;
  _fl_ x4;
} data_t;

void main (int argc, char * argv[])
{
  int i=0, j=0, b=0, k=0;
  char * fn=NULL;
  char * p=NULL;
  _fl_ x1, x1min=1.e9, x1max=-1.e9, x1bs;
  _fl_ x2, x2min=1.e9, x2max=-1.e9, x2bs;
  FILE * fp=stdin;
  int autoscale=0; 
  int nb1, nb2;
  int nData=0;
  double nNormData=0.0;
  data_t * data;
  double kT = 0.0, fe, BIGFE=1.e6;
  int weighted=0, ns;

  fprintf(stderr,"# 2DDIST (c) 2016 CFA\n");
  fflush(stderr);
  
  for (i=1;i<argc;i++) {
    if (argv[i][0] != '-') fn=argv[i];
    else if (!strcmp(argv[i], "-")) fp=stdin; /* not necessary */
    else if (!strcmp(argv[i], "-r1")) 
#ifndef DOUBLE_PRECISION
      sscanf(argv[++i],"%f,%f,%f",&x1min,&x1max,&x1bs);
#else
      sscanf(argv[++i],"%lf,%lf,%lf",&x1min,&x1max,&x1bs);
#endif
    else if (!strcmp(argv[i], "-r2"))
#ifndef DOUBLE_PRECISION 
      sscanf(argv[++i],"%f,%f,%f",&x2min,&x2max,&x2bs);
#else
      sscanf(argv[++i],"%lf,%lf,%lf",&x2min,&x2max,&x2bs);
#endif
    else if (!strcmp(argv[i], "-autoscale")) {
      autoscale=1;
      sscanf(argv[++i],"%i,%i",&nb1,&nb2);
    }
    else if (!strcmp(argv[i], "-kT")) kT = atof(argv[++i]);
    else if (!strcmp(argv[i], "-wt")) weighted=1;
    else usage();
  }

  nb1=(int)((x1max-x1min)/x1bs);
  nb2=(int)((x2max-x2min)/x2bs);

  data=(data_t*)calloc(MAXDATA,sizeof(data_t));
  
  memset(bins,0,MAXBINS*MAXBINS);
  memset(vals,0,MAXBINS*MAXBINS);

  if (nb1>MAXBINS||nb2>MAXBINS) {
    fprintf(stderr,"ERROR: Apparent number of bins %i,%i is too large (max %i).\n",nb1,nb2,MAXBINS);
    fprintf(stderr,"Recompile with larger MAXBINS\n");
    exit(-1);
  }
  if (fn) {
    fp=fopen(fn, "r");
    if (!fp) {printf("Error. Could not open %s.\n", fn);exit(-1);}
  }
  nData=0;
  while (fgets(scr_line_,MAXLINE,fp)) {
#ifndef DOUBLE_PRECISION
    ns=sscanf(scr_line_,"%f %f %f %f",&data[nData].x1,&data[nData].x2,&data[nData].x3,&data[nData].x4);
#else
    ns=sscanf(scr_line_,"%lf %lf %lf %lf",&data[nData].x1,&data[nData].x2,&data[nData].x3,&data[nData].x4);
#endif
    if (ns==2 && weighted) {
      data[nData].x3=1.0;
    } 
    if (autoscale) {
      if (data[nData].x1>x1max) x1max=data[nData].x1;
      if (data[nData].x1<x1min) x1min=data[nData].x1;
      if (data[nData].x2>x2max) x2max=data[nData].x2;
      if (data[nData].x2<x2min) x2min=data[nData].x2;
    }
    nData++;
  }

  if (autoscale) {
    x1bs = (x1max-x1min)/nb1;
    x2bs = (x2max-x2min)/nb2;
    x1min = x1min-2*x1bs;
    nb1+=2;
    x2min = x2min-2*x1bs;
    nb2+=2;
  }

  fflush(stderr);
  nNormData=0.0;
  for (k=0;k<nData;k++) {
    two_bin(&i,&j,data[k].x1,data[k].x2,x1min,x2min,x1bs,x2bs);
    //fprintf(stderr,"%i %i %i (%i %i)\n",k,i,j,MAXBINS,MAXBINS);
    if (i>=0&&j>=0) {
      if (weighted) {
         bins[i][j]+=data[k].x3;
         nNormData+=data[k].x3;
         //fprintf(stdout,"%i %.5le\n",k,nNormData);fflush(stdout);
      }
      else {
         bins[i][j]++;
         nNormData++;
      }
      vals[i][j]+=sqrt(data[k].x3*data[k].x3+data[k].x4*data[k].x4);
    }
  }
  fprintf(stderr,"#nData %i nNormData %.4lf: x1 [%.3lf,%.3lf,%.5lf]; x2 [%.3lf,%.3lf,%.5lf]\n",
	  nData,nNormData,x1min,x1max,x1bs,x2min,x2max,x2bs); 

  for (i=0;i<MAXBINS;i++) {
    if (binv(i,x1bs,x1min) >= x1min && binv(i,x1bs,x1min) <= x1max) {
      for (j=0;j<MAXBINS;j++) {
	if (binv(j,x2bs,x2min) >= x2min && binv(j,x2bs,x2min) <= x2max) {
	  printf("%.3lf %.3lf %.6le %.6le %i %.6le ",binv(i,x1bs,x1min),
		 binv(j,x2bs,x2min),((_fl_)bins[i][j])/nNormData,
		 ((_fl_)bins[i][j])/nNormData/(x1bs*x2bs),bins[i][j],bins[i][j]?vals[i][j]/bins[i][j]:0.0);
	  if (kT>0.0) {
	    fe=bins[i][j]?-kT*log(((_fl_)bins[i][j])/nNormData/(x1bs*x2bs)):BIGFE;
	    if (fe!=fe || fe>BIGFE) fe=BIGFE;
	    printf("%.6le",fe);
	  }
	  printf("\n");
	}
      }
      printf("\n");
    }
  }
}

