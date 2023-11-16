/*
 * fft.c    by cam abrams
 * (c) 1999,2000 cfa
 * 
 * computes a discrete fourier transform of a list of x_i, f(x_i) values.
 * 07Jul2000 Mainz, Germany
 *
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define MAXLINE 255
char scr_line_[MAXLINE];
#ifndef MAXDATA
#define MAXDATA 2048
#endif
float data[MAXDATA+1],x[MAXDATA+1];

void realft(float data[], unsigned long n, int isign); // recipes_c

void usage (void)
{
    printf("usage:\n");
    printf("fft filename \n");
    exit(-1);
}

int main (int argc, char * argv[])
{
    int i=0,sgn=1;
    char * fn=NULL;
    int nData=0;
    FILE * fp=stdin;
    char y1str[50], x1str[50];
    
    for (i=0;i<MAXDATA;i++) data[i]=x[i]=0.0;

    for (i=1;i<argc;i++)
    {
	if (argv[i][0] != '-') fn=argv[i];
        else usage();
    }
    if (fn)
    {
        fp=fopen(fn, "r");
        if (!fp) {printf("Error. Could not open %s.\n", fn);exit(-1);}
    }
    
    nData=0;
    while (fgets(scr_line_,MAXLINE,fp))
    {
	/* read the first number on the line */
	sscanf(scr_line_,"%s %s",x1str,y1str);
	data[++nData]=atof(y1str);
	x[nData]=atof(x1str);
    }
    if (fp!=stdin) fclose(fp);

    realft(data,nData-1,sgn);

    for (i=2;i<nData+1;i++)
      fprintf(stdout,"%.10le %.10le\n",1.0/x[i],data[i]);
    
}
