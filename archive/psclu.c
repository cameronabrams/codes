#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

/*
 * psclu.c -- processes simple cluster info from my polysoap runs
 * 
 */
#ifndef MAXCLUSTERS
#define MAXCLUSTERS 255
#endif
#ifndef MAXCLUSTPOP
#define MAXCLUSTPOP 255
#endif
#ifndef BUFSIZE
#define BUFSIZE 701
#endif

void usage (void) {fprintf(stderr,"no help for you!\n");}

int main (int argc, char * argv[])
{
  int i=0, j=0, k=0, m=0, n=0, dummy=0, nFrames=0, nAtomPerFrame=0, thisAtomID=0;
  int thisClustId, nClust=0, oldNClust=0, found=0, dt=0;
  int sideChainLength=5, dtmax=100;
  //  int clustPop[MAXCLUSTERS];
  //int clustMembers[MAXCLUSTERS][MAXCLUSTPOP];
  int clustPopBuf[BUFSIZE][MAXCLUSTERS];
  int clustMembersBuf[BUFSIZE][MAXCLUSTERS][MAXCLUSTPOP];
  int bufSz=0, clust_evol=0;
  int clustLifetime[BUFSIZE][MAXCLUSTERS];
  FILE * fp=stdin;
  char * fn=NULL;
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-h")) usage();
    else if (!strcmp(argv[i],"-dtmax")) dtmax=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-scl")) sideChainLength=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-evol")) clust_evol=atoi(argv[++i]);
    else fn=argv[i];
  }

  if (dtmax>=BUFSIZE) {
    fprintf(stderr,"Please choose a max time interval less than %i\n",BUFSIZE);
    exit(-1);
  }

  if (fn) fp=fopen(fn,"r");
  if (!fp) {fprintf(stderr,"Cannot open %s.\n",fn);exit(-1);}
  fscanf(fp,"%i",&nFrames);
  if (nFrames>BUFSIZE) {
    fprintf(stderr,"# Too many frames; recompile with BUFSIZE > %i\n",BUFSIZE);
    exit(-1);
  }
  memset(clustPopBuf,0,sizeof(clustPopBuf));
  memset(clustMembersBuf,0,sizeof(clustMembersBuf));
  /* read in the trajectory */
  //fprintf(stderr,"# input\n");fflush(stderr);
  for (i=0;i<nFrames;i++) {
    fscanf(fp,"%i",&nAtomPerFrame);
    for (j=0;j<nAtomPerFrame;j++) {
      fscanf(fp,"%i %i",&thisAtomID,&thisClustId);
      if (thisClustId>=0&&thisClustId<MAXCLUSTERS) {
	clustMembersBuf[i][thisClustId][clustPopBuf[i][thisClustId]++]=thisAtomID;
      }
      else {fprintf(stderr,"Error: bad cluster id: %i\n",thisClustId);exit(-1);}
    }
    if (clust_evol) {
      fprintf(stdout,"%d ",i);
      j=0;
      while (clustPopBuf[bufSz-1][j]) fprintf(stdout,"%d ",clustPopBuf[bufSz-1][j++]);
      fprintf(stdout,"\n");
    }
  }

  // fprintf(stderr,"# calc\n");fflush(stderr);
  memset(clustLifetime,0,sizeof(clustLifetime));
  for (i=0;i<nFrames-1;i++) {
    for (m=0;clustPopBuf[i][m];m++) { /* all clusters `m' at time `i' */
      if (clustPopBuf[i][m]!=-1) { /* this cluster alive at an earlier time */
	for (j=i+1;j<nFrames;j++) {
	  dt=j-i;
	  found=0;
	  for (n=0;clustPopBuf[j][n];n++) { /* all clusters `n' at time j>i */
	    if ((clustPopBuf[i][m]==clustPopBuf[j][n])&&
		(!memcmp(clustMembersBuf[i][m],clustMembersBuf[j][n],
			 clustPopBuf[i][m]*sizeof(int)))) {
		clustPopBuf[j][n]=-1;
		found=1;
	    }
	  }
	  if (!found) { /* cluster m at time i not found at time j */
	    if (!clustLifetime[i][m]) clustLifetime[i][m]=dt;
	  }
	}
      }
    }
  }
  for (i=0;i<nFrames;i++) {
    for (j=0;j<MAXCLUSTERS;j++) {
      if (clustLifetime[i][j]) fprintf(stdout,"%i %i\n",clustPopBuf[i][j],clustLifetime[i][j]);
    }
  }
    

  return 0;
}
