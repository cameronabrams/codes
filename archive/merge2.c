#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

/*
 * merge2.c -- merges any number of {Y}-vs-X data files
 * to produce a single {Y}-vs-X output file, where Y_i is
 * an average of Y_i's from each file.
 *
 * {Y} is a set of columns, and each file is assumed to have the same
 * format of one {Y}-vs-X point per line as:
 *
 * x-value y1-value y2-value ...
 *
 * Input lines beginning with "%" or "#" are ignored.
 * All input fields (except the first column of x-data) are assumed to be 
 * floating-point values.
 *
 * Y-data output is optionally three columns for each column of y-data input:
 * (mean) "+/-" (std.dev.)
 * 
 * followed by the count of instances which produced the average in that row.
 *
 * Each file need not have the same number of lines; the average is 
 * performed only those values available.
 *
 * The first file in the argument list is by default the "pattern" file.
 * 
 * 
 * 
 * (c) 1999-2011 cameron abrams
 * 
 */

#define MAXFILES 100 
#define MAXLINES 5000
#define MAXCOLUMNS 100

int ln2row (char * p, double row[])
{
  char substr[255];
  int i=0;
  if (!p) return 0;
//    printf("row [%s]",p);
  while (*p&&sscanf(p,"%s",substr)) {
    p+=strlen(substr);
    while(isspace(*(++p)));
    if (substr[strlen(substr)-1]==',') substr[strlen(substr)-1]='\0';
//        printf("  -> convert [%s] --> [%.5le]\n", substr,atof(substr));
    if (strcmp(substr, "+/-")) row[i++]=atof(substr);
  }
  return i;
}

void ln2hln (char * hln, char * p)
{
  char substr[255], * q = hln;
  int i=0;
  if (!p||!q) return;
  while (*p&&sscanf(p,"%s",substr)) {
    p+=strlen(substr);
    while(isspace(*(++p)));
    if (q==hln) q+=sprintf(q,"%s\t",substr);
    else q+=sprintf(q,"+-- %s --+\t",substr);
  }
  *(--q)='\n'; *(++q)='\0';
}


typedef struct FILEDATA {
  int N; // number of lines
  int M; // number of columns
  double ** d; // pointer to data
} fdata;

typedef struct MRGDATA {
  int N; // number of lines
  int M; // number of columns
  int * c; // counts per column
  double ** d; // pointer to data
  double ** s; // pointer to stddev
} mdata;

fdata * new_filedata ( void ) {
  int i,j; 
  fdata * f = malloc(sizeof(fdata));
  f->N=f->M=0;
  f->d=(double**)malloc(MAXLINES*sizeof(double*));
  for (i=0;i<MAXLINES;i++) {
    f->d[i]=(double*)malloc(MAXCOLUMNS*sizeof(double));
  }
  for (i=0;i<MAXLINES;i++) {
    for (j=0;j<MAXCOLUMNS;j++) {
      f->d[i][j]=0.0;
    }
  }
  return f;
}

mdata * new_mrgdata ( void ) {
  int i,j; 
  mdata * f = malloc(sizeof(fdata));
  f->N=f->M=0;
  f->c=(int*)malloc(MAXCOLUMNS*sizeof(int));
  f->d=(double**)malloc(MAXLINES*sizeof(double*));
  f->s=(double**)malloc(MAXLINES*sizeof(double*));
  for (i=0;i<MAXLINES;i++) {
    f->d[i]=(double*)malloc(MAXCOLUMNS*sizeof(double));
    f->s[i]=(double*)malloc(MAXCOLUMNS*sizeof(double));
  }
  for (i=0;i<MAXLINES;i++) {
    for (j=0;j<MAXCOLUMNS;j++) {
      f->d[i][j]=0.0;
      f->s[i][j]=0.0;
    }
  }
  for (j=0;j<MAXCOLUMNS;j++) {
    f->c[j]=0;
  }
  return f;
}

fdata * get_file_data ( char * fn ) {
  int i,nf;
  char ln[1000];
  FILE * fp = fopen(fn,"r");
  if (fp) {
    fdata * f = new_filedata();
    
    //    fprintf(stderr,"Reading file %s...\n",fn);
    //    fflush(stderr);

    i=0;
    while (fgets(ln,1000,fp)) {
      if (ln[0]!='#'&&ln[0]!='%') {
	//	fprintf(stderr,"  -> reading data from file into line %i...\n",i);fflush(stderr);
	nf=ln2row(ln,f->d[i]);
	if (!i) f->M=nf;
	i++;
      }
    }
    f->N=i;
    close(fp);
    //   fprintf(stderr,"Finished read-in of %s: data is %i x %i\n",fn,f->N,f->M);
    return f;
  }
  else return NULL;
}

mdata * mrg_file_data (fdata ** file_data, int N) {
  int i,j,k,c;
  mdata * m = new_mrgdata();
  double s,ss;

  int maxN=0;
  for (i=0;i<N;i++) {
    if (file_data[i]->N>maxN) maxN=file_data[i]->N;
  }

  m->N=maxN;
  m->M=0;
  for (i=0;i<maxN;i++) {
    for (j=0;j<MAXCOLUMNS;j++) {
      c=0;
      s=ss=0.0;
      for (k=0;k<N;k++) {
	if (i<file_data[k]->N&&j<file_data[k]->M) {
	  s+=file_data[k]->d[i][j];
	  ss+=file_data[k]->d[i][j]*file_data[k]->d[i][j];
	  c++;
	}
	else if (!m->M) {
	  m->M=j;
	}
      }
      if (c) {
	m->d[i][j]=s/c;
	m->s[i][j]=sqrt((ss-s*s/c)/c);
      }
    }
  }

  return m;
}

void mrg_data_out ( FILE * fp, mdata * m, int isd ) {

  int i,j;
  
  for (i=0;i<m->N;i++) {
    for (j=0;j<m->M;j++) {
      fprintf(fp,"%.6lf ",m->d[i][j]);
      if (isd) fprintf(fp,"sd-> %.3lf <-sd  ",m->s[i][j]);
    }
    fprintf(fp,"\n");
  }

}

int main (int argc, char * argv[])
{
  fdata ** file_data=NULL;
  mdata * mrg;
  char * fn[MAXFILES];
  int includeStddev=0;
  char * fmt = "% .5lf";

  int nFiles=0,i;

  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-sd")) includeStddev=1;
    else if (!strcmp(argv[i],"-fmt")) fmt=argv[++i];
    else fn[nFiles++]=argv[i];
  }

  fprintf(stderr, "# merging %i file%s\n", nFiles, nFiles>1?"s":"");
  fflush(stderr);

  file_data = malloc(nFiles*sizeof(fdata*));
  for (i=0;i<nFiles;i++) {
    file_data[i]=get_file_data(fn[i]);
  }

  mrg=mrg_file_data(file_data,nFiles);

  mrg_data_out(stdout,mrg,includeStddev);

  fprintf(stderr,"Program ends.\n");
}
