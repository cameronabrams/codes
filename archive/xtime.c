#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

enum {DAY,MONTH,YEAR,HOUR,MINUTE};
char * fstr[] = {"day","month","year","hour","minute"};

void blank_tm (struct tm * t) {
  if (!t) return;
  t->tm_sec=t->tm_min=t->tm_hour=t->tm_mday=t->tm_mon=
  t->tm_year=t->tm_wday=t->tm_yday=t->tm_isdst=0;
}

void copy_tm (struct tm * t1, struct tm * t2) {
  if (!t1||!t2) return;
  t1->tm_sec=t2->tm_sec;
  t1->tm_min=t2->tm_min;
  t1->tm_hour=t2->tm_hour;
  t1->tm_mday=t2->tm_mday;
  t1->tm_mon=t2->tm_mon;
  t1->tm_year=t2->tm_year;
  t1->tm_wday=t2->tm_wday;
  t1->tm_yday=t2->tm_yday;
  t1->tm_isdst=t2->tm_isdst;
}

int main ( int argc, char * argv[] ) {
  struct tm t0, * pt0; short init=0;
  struct tm t;
  time_t it,it0;
  FILE * fp=NULL;
  char * fn=NULL;
  char ln[255],dstr[100];
  char * dfmt=NULL;
  int ord[] = {DAY,MONTH,YEAR,HOUR,MINUTE};
  int *tdat;
  char fmt[255]="%d/%d/%d.%d:%d"; /*default*/
  float y,h;
  int i,p,s;
  int nf=sizeof(ord)/sizeof(int);

  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-f")) fn=argv[++i];
    else if (!strcmp(argv[i],"-dfmt")) dfmt=argv[++i];
  }

  if (dfmt) {
    for (p=0;fmt[p]!='\0';p++) fmt[p]='\0';
    i=0;
    for (p=0;dfmt[p]!='\0';p++) {
      if (dfmt[p]=='%') {
	fmt[p]=dfmt[p++];
	switch (dfmt[p]) {
	  case 'd':  ord[i]=DAY; break;
	  case 'm':  ord[i]=MONTH; break;
	  case 'Y':  ord[i]=YEAR;  break;
	  case 'M':  ord[i]=MINUTE; break;
	  case 'H':  ord[i]=HOUR; break;
	  default:     fprintf(stderr,"dfmt error\n");exit(-1);
	}
	fmt[p]='d';
	i++;
      }
      else fmt[p]=dfmt[p];
    }
    nf=i;
  }

  tdat=(unsigned int*)malloc((unsigned)nf*sizeof(int));
  for (i=0;i<nf;i++) tdat[i]=0;

  //fprintf(stderr,"dfmt %s fmt %s ",dfmt,fmt);
  //for (i=0;i<nf;i++) fprintf(stderr,"%s ",fstr[ord[i]]);
  //fprintf(stderr,"\n");
  
  fp=fn?fopen(fn,"r"):stdin;

  time(&it);
  pt0=localtime(&it);
  p=pt0->tm_isdst;
  blank_tm(&t0);
  t0.tm_isdst=p;
  blank_tm(&t);

  while (fgets(ln,255,fp)) {
    blank_tm(&t);
    sscanf(ln,"%s %f",dstr,&y);
    //fprintf(stderr,"d `%s' y `%f'\n",dstr,y);fflush(stderr);
    sscanf(dstr,fmt,tdat,tdat+1,tdat+2,tdat+3,tdat+4,tdat+5,tdat+6,tdat+7);
    for (i=0;i<nf;i++) {
      switch (ord[i]) {
        case DAY:    t.tm_mday=tdat[i]; break;
        case MONTH:  t.tm_mon=tdat[i]-1; break;
        case YEAR:   t.tm_year=tdat[i]-1900; break;
        case HOUR:   t.tm_hour=tdat[i]; break;
        case MINUTE: t.tm_min=tdat[i]; break;
      }
    }
    if (!init) {
      copy_tm(&t0,&t);
      it0=mktime(&t0);
      init=1;
    }
    it=mktime(&t);
    s=difftime(it,it0);
    h=(float)s/60.0/60.0;
    fprintf(stdout,"%.3f %.3f\n",h,y);
  }  
}
