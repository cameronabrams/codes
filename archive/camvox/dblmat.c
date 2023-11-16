/*
 * MD series -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1996-1998 Cam Abrams, Dept. of Chemical Engineering
 * Regents of the University of California, Berkeley
 * 
 * Module Name:		dblmat
 * Module Class:	Auxiliary
 * Module Description:  floating point square matrices
 * 
 */
#include "dblmat.h"

void dblmat_zero(double mat[][DBLMAT_MAXDIM], int dim)
{
    int i = 0, j = 0;
    if (dim > DBLMAT_MAXDIM)
    {
	printf("Error -- matrix too big.\n");
	return;
    }
    for (i = 0; i < dim; i++)
    {
	for (j = 0; j < dim; j++)
	{
	    mat[i][j] = 0.0;
	}
    }
}

void dblmat_identity(double mat[][DBLMAT_MAXDIM], int dim)
{
    int i = 0, j = 0;
    if (dim > DBLMAT_MAXDIM)
    {
	printf("Error -- matrix too big.\n");
	return;
    }
    for (i = 0; i < dim; i++)
    {
	for (j = 0; j < dim; j++)
	{
	    mat[i][j] = 0.0;
	}
	mat[i][i] = 1.0;
    }
}

void dblmat_copy(double dest[][DBLMAT_MAXDIM],
	         double orig[][DBLMAT_MAXDIM], int dim)
{
    int i = 0, j = 0;
    if (dim > DBLMAT_MAXDIM)
    {
	printf("Error -- matrix too big.\n");
	return;
    }
    for (i = 0; i < dim; i++)
    {
	for (j = 0; j < dim; j++)
	{
	    dest[i][j] = orig[i][j];
	}
    }
}

void dblmat_add(double dest[][DBLMAT_MAXDIM], 
		double a[][DBLMAT_MAXDIM],
		double b[][DBLMAT_MAXDIM], int dim)
{
    int i = 0, j = 0;
    if (dim > DBLMAT_MAXDIM)
    {
	printf("Error -- matrix too big.\n");
	return;
    }
    for (i = 0; i < dim; i++)
    {
	for (j = 0; j < dim; j++)
	{
	    dest[i][j] = a[i][j] + b[i][j];
	}
    }
}

void dblmat_matmult(double dest[][DBLMAT_MAXDIM], 
		double a[][DBLMAT_MAXDIM],
		double b[][DBLMAT_MAXDIM], int dim)
{
    int i = 0, j = 0, k = 0;
    if (dim > DBLMAT_MAXDIM)
    {
	printf("Error -- matrix too big.\n");
	return;
    }
    
#if DIAG
    printf("A=\n");
    for (i = 0; i < dim; i++)
    {
	for (j = 0; j < dim; j++)
	{
	    printf("%.2lf\t", a[i][j]);
	}
	printf("\n");
    }
    printf("B=\n");
    for (i = 0; i < dim; i++)
    {
	for (j = 0; j < dim; j++)
	{
	    printf("%.2lf\t", b[i][j]);
	}
	printf("\n");
    }
#endif
    for (i = 0; i < dim; i++)
    {
	for (j = 0; j < dim; j++)
	{
	    dest[i][j] = 0.0;
	    for (k = 0; k < dim; k++)
	    {
		dest[i][j] += a[i][k]*b[k][j];
	    }
	}
    }
}

void dblmat_out(FILE * fp, double a[][DBLMAT_MAXDIM], int dim)
{
    int i = 0, j = 0;
    for (i = 0; i < dim; i++)
    {
	fprintf(fp, "#");
	for (j = 0; j < dim; j++)
	{
	   fprintf(fp, "%.2lf\t", a[i][j]);
	}
	fprintf(fp, "\n");
    }    
}

