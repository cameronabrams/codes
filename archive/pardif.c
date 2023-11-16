#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

/*
 * pardif.c -- subtracts two parallel data files to produce a third
 * parallel data file.  Two files are defined as `parallel' if they:
 * (a) have the same number of rows, and 
 * (b) have the same number of fields per row, and
 * (c) the first field of the same row in all files is the same.
 *
 * All fields (except the first column of x-data) are assumed to be 
 * floating-point values.
 * 
 * The output is to stdout, and consists of a single file where each
 * field `f_i' of each row `i' is the difference of these two fields
 * f_i(file1) - f_i(file2)
 * 
 * (c) 2001 cameron abrams
 * 
 */
 
char ln[5000], hln[5000];
#define MAXFIELDS 50
FILE * fp[2];
char * fn[2];

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

int main (int argc, char * argv[])
{
    int i=0, j=0, nf=MAXFIELDS;
    double buf[MAXFIELDS];
    double difs[MAXFIELDS];
    int nFiles=argc-1;
    char * reading,  * p;
    char x[25];
    for (i=1;i<argc;i++) fn[i-1]=argv[i];
    if (nFiles!=2) {
      fprintf(stderr,"Need two file names to pardif\n");
      exit(-1);
    }
    fprintf(stderr, "# pardiffing %i files\n", nFiles);
    
    for (i=0;i<nFiles;i++) fp[i]=fopen(fn[i], "r");
    reading=fn[0];
    while (reading)
    {
	for (j=0;j<nf;j++) difs[j]=0.0;
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
		for (j=0;j<nf;j++) {difs[j]=difs[j]+(i?-buf[j]:buf[j]);}
	    }
//	    if (i==0&&(ln[0]=='#'||ln[0]=='%')) {ln2hln(hln,ln);}
	}
	if (reading&&ln[0]!='#'&&ln[0]!='%')
	{
      	    printf("%s", x);
	    for (j=0;j<nf;j++) 
	    {
		printf("\t%.5le", difs[j]);
	    }
	    printf("\n");
	}
    }
    for (i=0;i<nFiles;i++) fclose(fp[i]);
    
}
