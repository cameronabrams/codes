#include <stdio.h>
#include <math.h>
#include <string.h>
/*
 * collapse.c -- collapses blocks of contiguous lines in a column-
 * oriented data file into single lines.  A block of lines is
 * defined such that each member of the block has the same 
 * element values in at least one of the fields (by default the
 * first).
 * 
 * 
 * (c) 2001 cam abrams
 * max-planck-institute for polymer research
 * mainz, germany
 */

uint numFields ( char * ln, uint len )
{
  uint i=0;
  char *p=ln;
  while (*p) {
    while (*p&&isspace(*p)) p++;
    if (*p&&(isprint(*p)&&!isspace(*p))) {
      i++;
      while (*p&&(isprint(*p)&&!isspace(*p))) p++;
    }
  }
  return i;
}

void ln2buf ( float * X, char * ln, uint n )
{
  uint i=0;
  char str[80];
  char *p=ln;
  for (i=0;i<n;i++) {
    while (*p&&isspace(*p)) p++;
    sscanf(p,"%s",str);
    //fprintf(stdout,"{%s}[%i]\n",str,strlen(str));
    p+=strlen(str);
    X[i]=atof(str);
  }
}

void scalbuf ( float * X, float x, uint n )
{
  uint i=0;
  for (i=0;i<n;i++) X[i]*=x;
}

void copybuf ( float * D, float * S, uint n )
{
  uint i=0;
  for (i=0;i<n;i++) D[i]=S[i];
}

void addbuf ( float * D, float * S, uint n )
{
  uint i=0;
  for (i=0;i<n;i++) D[i]+=S[i];
}

char buf_[255];
char * prbuf ( float * X, char * nfmt, uint n )
{
  uint i=0;
  char fmtStr[80];
  buf_[0]=(char)0;
  sprintf(fmtStr,"%s%%s",nfmt);
  for (i=0;i<n;i++) {
    sprintf(&buf_[strlen(buf_)],fmtStr,X[i],(i==(n-1)?"":" "));
  }
  return buf_;
}

uint nData, NF;
char ln[255];
int main (int argc, char * argv[])
{
    FILE *fp=NULL;
    char *fn=NULL;
    char sfmt[255];
    char * nfmt="%.5f";
    uint i,cnt;
    int com=0;
    float *OX=NULL,*NX=NULL,*S=NULL;
    
    for (i=1;i<argc;i++)
    {
	if (argv[i][0]!='-') fn=argv[i];
        else if (!strcmp(argv[i],"-nf")) nfmt=argv[++i];
        else if (!strcmp(argv[i],"-com")) com=atoi(argv[++i])-1;
    }

    if (!fn) fp=stdin;
    else fp=fopen(fn, "r");
    
    i=0;cnt=0;
    while (fgets(ln,255,fp))
    {
      if (ln[0]!='#' && !isspace(ln[0]))
      {
	if (!OX&&!NX) {
	  OX=(float*)calloc(NF=numFields(ln,255),sizeof(float));
	  NX=(float*)calloc(NF,sizeof(float));
	  S=(float*)calloc(NF,sizeof(float));
	  ln2buf(OX,ln,NF);
	  copybuf(S,OX,NF);
	  //fprintf(stdout,"[0] [%s]\n",prbuf(S,nfmt,NF));
	  cnt=1;
	}
	else {
	  ln2buf(NX,ln,NF);
	  //fprintf(stdout,"[%i] [%s]\n",i,prbuf(NX,nfmt,NF));
	  if (NX[com]!=OX[com]) {
	    scalbuf(S,1.0/cnt,NF);
	    fprintf(stdout,"%s # %i\n",prbuf(S,nfmt,NF),cnt);
	    copybuf(OX,NX,NF);
	    copybuf(S,OX,NF);
	    cnt=1;
	  }
	  else {
	    copybuf(OX,NX,NF);
	    addbuf(S,OX,NF);
	    cnt++;
	  }
	}
	i++;
      }
    }
    scalbuf(S,1.0/cnt,NF);
    fprintf(stdout,"%s # %i\n",prbuf(S,nfmt,NF),cnt);
    nData=i;
    if (fp!=stdin) fclose(fp);

}
