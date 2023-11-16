#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

/*
 * psclu2.c -- processes simple cluster info from my polysoap runs
 * 
 */


#define sqr(a) ((a)*(a))

void usage (void) {fprintf(stderr,"no help for you!\n");}

typedef struct FRAME {
  int id;
  double time;
  int n;
  int *agg;
  int **mem;
  double **dist;
  double *com[3];
} frame_t;

void frame_out ( frame_t * F )
{
  int i,j;
  fprintf(stdout,"frame time %lf, clusters %i\n",F->time,F->n);
  for (i=0;i<F->n;i++) {
    fprintf(stdout,"    cluster %i pop %i com %lf %lf %lf\n",
	    i,F->agg[i],F->com[0][i],F->com[1][i],F->com[2][i]);
    fprintf(stdout,"    members: ");
    for (j=0;j<F->agg[i];j++) {
      fprintf(stdout," %i",F->mem[i][j]);
    }
    fprintf(stdout,"\n");
  }
}

void intarr_sort ( int * A, int n ) {

  int i, sorted=0,swp;
  while (!sorted) {
    sorted=1;
    for (i=0;i<n-1;i++)
      if (A[i]>A[i+1]) {swp=A[i];A[i]=A[i+1];A[i+1]=swp;sorted=0;}
  }
}

enum { tDist, realSpcCorr };

int main (int argc, char * argv[])
{
  unsigned int toDo = tDist;
  frame_t * trj;
  int i=0, j=0, k=0, dtmax=100, nFrames=0, sideChainLength=5, dummy=0;
  FILE * fp=stdin;
  char * fn=NULL, *p=NULL;
  char scrln[1000];
  int bin[20];
  double sums[6];
  int count[6];

  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-h")) usage();
    else if (!strcmp(argv[i],"-dtmax")) dtmax=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-scl")) sideChainLength=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-tdist")) toDo = tDist;
    else if (!strcmp(argv[i],"-rsc")) toDo = realSpcCorr;
    else fn=argv[i];
  }

  memset(sums,0,6*sizeof(double));
  memset(count,0,6*sizeof(int));

  if (fn) fp=fopen(fn,"r");
  if (!fp) {fprintf(stderr,"Cannot open %s.\n",fn);exit(-1);}
  fscanf(fp,"%i",&nFrames);
  trj=(frame_t*)malloc(nFrames*sizeof(frame_t));
  if (!trj) {
    fprintf(stderr,"Error: not enough memory for %d frames.\n",nFrames);
    exit(-1);
  }


  /* read in the trajectory */
  //fprintf(stderr,"# input\n");fflush(stderr);
  for (i=0;i<nFrames;i++) {
    memset(bin,0,20*sizeof(int));
    fscanf(fp,"%d %lf\n",&trj[i].n,&trj[i].time);
    //fprintf(stderr,"Frame %i has %i clusters at time %lf\n",i,trj[i].n,trj[i].time);fflush(stderr);
    for (j=0;j<3;j++) trj[i].com[j]=(double*)malloc(trj[i].n*sizeof(double));
    trj[i].agg=(int*)malloc(trj[i].n*sizeof(int));
    trj[i].mem=(int**)malloc(trj[i].n*sizeof(int*));
    trj[i].dist=(double**)malloc(trj[i].n*sizeof(double*));
    //    fscanf(fp,"%s",scrln); // skip the # line
    fgets(scrln,1000,fp); //fprintf(stderr,"scrln [%s]",scrln); fflush(stderr);
    for (j=0;j<trj[i].n;j++) {
      //      fgets(scrln,1000,fp);
      fscanf(fp,"%i %i %lf %lf %lf",&dummy,&trj[i].agg[j],&trj[i].com[0][j],
	     &trj[i].com[1][j],&trj[i].com[2][j]);
      bin[trj[i].agg[j]/5-1]++;
      trj[i].mem[j]=(int*)malloc(trj[i].agg[j]*sizeof(int));
      trj[i].dist[j]=(double*)malloc(trj[i].agg[j]*sizeof(double));
      if (!trj[i].mem[j]) {
	fprintf(stderr,"Error: out of memory\n");
	exit(-1);
      }
      for (k=0;k<trj[i].agg[j];k++) fscanf(fp,"%i",&trj[i].mem[j][k]); 
      for (k=0;k<trj[i].agg[j];k++) fscanf(fp,"%lf",&trj[i].dist[j][k]); 
      fscanf(fp,"\n");
      if (trj[i].agg[j]==6) {
	for (k=0;k<6;k++){ 
	  //	  fprintf(stdout,"%i %lf\n",k,trj[i].dist[j][k]); 
	  sums[k]+=trj[i].dist[j][k];
	  count[k]++;
	}
      }
    }
    if (toDo == tDist) 
      for (j=5;j<55;j+=5) printf("%.5lf %i %i\n",trj[i].time,j,bin[j/5-1]);
  }

  
  if (toDo == realSpcCorr)
    for (k=0;k<6;k++) { 
      fprintf(stdout,"%i %lf %i\n",k,sums[k],count[k]); 
    }

  return 0;
}
