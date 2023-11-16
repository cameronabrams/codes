#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

/*
 * merge.c -- merges an arbitrary number (MAX 25) of parallel data files 
 * to produce an average.  All files must (a) have the same
 * number of rows and (b) have the same number of fields per row.
 * All fields (except the first column of x-data) are assumed to be 
 * floating-point values.
 * 
 * Y-data output is three columns for each column of y-data input:
 * (mean) "+/-" (std.dev.)
 * 
 * (c) 1999 cameron abrams
 * 
 */
 
char ln[5000], hln[5000];
#define MAXFILES 1001
#define MAXFIELDS 50
FILE * fp[MAXFILES];
char * fn[MAXFILES];

int ln2row (char * p, double row[])
{
    char substr[255];
    int i=0;
    if (!p) return 0;
//    printf("row [%s]",p);
    while (*p&&sscanf(p,"%s",substr))
    {
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
    while (*p&&sscanf(p,"%s",substr))
    {
	p+=strlen(substr);
	while(isspace(*(++p)));
	if (q==hln) q+=sprintf(q,"%s\t",substr);
	else q+=sprintf(q,"+-- %s --+\t",substr);
    }
    *(--q)='\n'; *(++q)='\0';
}


int main (int argc, char * argv[])
{
    int i=0, j=0, nf=MAXFIELDS;
    short header;
    double buf[MAXFIELDS];
    double sums[MAXFIELDS], ssums[MAXFIELDS];
    int nFiles=argc-1;
    char * reading,  * p;
    char x[25];
    for (i=1;i<argc;i++) fn[i-1]=argv[i];
    fprintf(stderr, "# merging %i files\n", nFiles);
    
    for (i=0;i<nFiles;i++) fp[i]=fopen(fn[i], "r");
    reading=fn[0];
    header=0;
    while (reading)
    {
	for (j=0;j<nf;j++) sums[j]=ssums[j]=0.0;
	for (i=0;i<nFiles;i++)
	{
	    reading=fgets(ln,255,fp[i]);
	    if (reading&&ln[0]!='#'&&ln[0]!='%') 
	    {
//		if (i==0&&!header) {printf("h %s",hln);header=1;}
		p=ln;
		/* might be some white space in front of the first # */
		while(isspace(*(p++))); p--;
		/* extract the first substring as the x-value */
		sscanf(p,"%s",x);
//		printf("f[%i] p->[%c][%s]\n",i,p[0],p+1);
//		printf("  -> x-value(str) [%s] l[%i]\n", x,strlen(x));
		p+=strlen(x);
//		printf(" -> p->[%c][%s]\n",p[0],p+1);
		/* the next statement assumes there is at least ONE
		   character separating two columns. */
		while(isspace(*(++p)));

//		printf(" -> p->[%c][%s]\n",p[0],p+1);
		nf=ln2row(p, buf);
		for (j=0;j<nf;j++) {sums[j]+=buf[j];ssums[j]+=buf[j]*buf[j];}
	    }
//	    if (i==0&&(ln[0]=='#'||ln[0]=='%')) {ln2hln(hln,ln);}
	}
	if (reading&&ln[0]!='#'&&ln[0]!='%')
	{
      	    printf("%s", x);
	    for (j=0;j<nf;j++) 
	    {
		/* compute each average and standard deviation */
		sums[j]/=nFiles;
		ssums[j]=sqrt(ssums[j]/nFiles-sums[j]*sums[j]);
/*		printf("\t%.3lf +/- %.3lf", sums[j], ssums[j]);*/
		printf("\t%.3le", sums[j]);
	    }
	    printf("\n");
	}
    }
    for (i=0;i<nFiles;i++) fclose(fp[i]);
    
}
