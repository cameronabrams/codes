/* transpose.c
 * (c) 2001 cameron abrams
 * max-planck-institute for polymer research
 *
 * This program transposes any column oriented data file.
 * 
 * i.e., a file that looks like
 
        a b
        c d
 *
 * becomes
 *
        a c
        b d
 *
 * One limitation is that the maximum number of rows/columns is 1000.
 *
 */

#include <stdio.h>
#include <string.h>

char ln[25000];

char matrix[1000][1000][25];

int main (int argc, char * argv[])
{
  FILE * fp=stdin;
  char * p;
  int n,m,i,j,ncol=0;

  if (argc==2) fp=fopen(argv[1],"r");
  if (!fp) {fprintf(stderr,"Error. Can't open %s\n",argv[1]);exit(-1);}
  n=0;
  while (fgets(ln,25000,fp)) {
    p=ln;
    m=0;
    while (*p&&sscanf(p,"%s",matrix[n][m])) {
      //fprintf(stderr,"%s ",matrix[n][m]);
      m++;
      while (*p&&isgraph(*p))p++; 
      while (*p&&isspace(*p))p++; 
    }
    if (!ncol) ncol=m;
    //    fprintf(stderr,"\n");
    //    fprintf(stderr,"# ln %i elements %i\n",n,m);
    //    for (i=0;i<m;i++)
    //      fprintf(stderr,"%i (%s) ",i,matrix[n][i]);
    //    fprintf(stderr,"\n");
    n++;
  }
  for (i=0;i<ncol;i++) {
    for (j=0;j<n;j++) {
      fprintf(stdout,"%s ",matrix[j][i]);
    }
    fprintf(stdout,"\n");
  }
}
