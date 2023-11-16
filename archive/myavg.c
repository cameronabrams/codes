/* (c) 2009 Cameron F Abrams */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct POINT * ppt;

typedef struct POINT {
  int i;
  int nd;
  double * d;
  ppt next;
} pt;

ppt newpt ( int i, int nd, double * d ) {
  ppt p=malloc(sizeof(pt));
  int j;
  //fprintf(stderr,"newpt : %i %i\n",i,nd);fflush(stderr);
  p->i=i;
  p->nd=nd;
  p->d=calloc(nd,sizeof(double));
  if (d) memcpy(p->d,d,nd*sizeof(double));
  p->next=NULL;
/*   fprintf(stderr,"newpt rtn: %i %i : ",p->i,p->nd);fflush(stderr); */
/*   for (j=0;j<p->nd;j++) { */
/*     fprintf(stderr,"%lf ",p->d[j]);fflush(stderr); */
/*   } */
/*   fprintf(stderr,"\n");fflush(stderr); */
  //fprintf(stderr,"newpt rtn\n");
  return p;
}

double mydist2 ( ppt A, ppt B ) {
  double d=0.0;
  int i;
  // fprintf(stderr,"mydist2: A %i %i B %i %i\n",
  //	  A->i,A->nd,B->i,B->nd);fflush(stderr);
  /* first element is time */
  for (i=1;i<A->nd;i++) {
    // fprintf(stderr,"mydist2 %i %i ... %lf %lf\n",
    //	    A->i,B->i,A->d[i],B->d[i]);fflush(stderr);
    d+=(A->d[i]-B->d[i])*(A->d[i]-B->d[i]);
  }
  // fprintf(stderr,"mydist2 %i %i %lf\n",A->i,B->i,d);fflush(stderr);
  return d;
}

int main ( int argc, char * argv[] ) {
  char ln[255];
  char dstr[255];
  pt * P = NULL, * pp = NULL;
  double * s = NULL;
  double * td = NULL;
  int npi=0;
  int n = 0;
  int i =0; 
  int f = 0;
  int nf=0;
  char * p;
  char cum=0;
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-cum")) cum=1;
    if (!strcmp(argv[i],"-npi")) npi=1;
    if (!strcmp(argv[1],"-nf")) nf=atoi(argv[++i]);
  }
  while (fgets(ln,255,stdin)) {
    p=ln;
    if (ln[0]!='#') {
      if (!n) {
	f=0;
	while (*p && isprint(*p) && *p!='\n') {
	  if (isspace(*p)) while(isspace(*(++p)));
	  //printf("[%s]\n",p);
	  sscanf(p,"%s",dstr);
	  //printf("%s\n",dstr);
	  p+=strlen(dstr);
	  f++;
	}
	f--;
	//printf("%i fields\n",f);
	s=(double*)malloc(f*sizeof(double));
	td=(double*)malloc(f*sizeof(double));
	for (i=0;i<f;i++) s[i]=0.0;
      }
      p=ln;
      for (i=0;i<f;i++) {
	if (isspace(*p)) while(isspace(*(++p)));
	sscanf(p,"%s",dstr);
	p+=strlen(dstr);
	td[i]=atof(dstr);
	s[i]+=td[i];
	if (cum) {
	  if (!i) printf("%s ",dstr);
	  else printf("%.4lf ",s[i]/(n+1));
	}
      }
      if (cum) printf("\n");
      if (npi) {
	if (!P) {
	  P=newpt(n,f,td);
	  pp=P;
	}
	else {
	  pp->next=newpt(n,f,td);
	  pp=pp->next;
	}
      }
      n++;
    }
  }
  if (nf==0) nf=f;
  if (!cum && !npi) {
    for (i=0;i<nf;i++) {
      printf("%.4lf ",s[i]/n);
    }
    printf("\n");
  }

  if (npi) {
    int j;
    double mind,dd;
    ppt minp=P;
    ppt A;
    for (i=0;i<f;i++) {
      s[i]/=n;
    }
    A=newpt(-1,f,s);
 /*    fprintf(stderr,"newpt A: %i %i : ",A->i,A->nd);fflush(stderr); */
/*     for (j=0;j<A->nd;j++) { */
/*       fprintf(stderr,"%lf ",A->d[j]);fflush(stderr); */
/*     } */
/*     fprintf(stderr,"\n");fflush(stderr); */
    mind=1.e9;
    for (pp=P;pp;pp=pp->next) {
      //fprintf(stderr,"calling mydist2 on p %i %i A %i %i\n",pp->i,pp->nd,A->i,A->nd);fflush(stderr);
      dd=mydist2(pp,A);
      if (dd<mind) {
	mind=dd;
	minp=pp;
      }
    }
    printf("%i ",minp->i);
    for (j=1;j<A->nd;j++) {
      fprintf(stdout,"%lf ",A->d[j]);
    }
    fprintf(stdout,"\n");
/*     printf("INFO: avg %lf %lf pt %i %lf %lf\n", */
/* 	   A->d[1],A->d[2],minp->i,minp->d[1],minp->d[2]); */
  }

  return 0;
}
