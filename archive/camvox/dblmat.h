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
#ifndef DBLMAT_H
#define DBLMAT_H

#include <math.h>
#include <stdio.h>

#define DBLMAT_MAXDIM	4

void dblmat_zero(double mat[][DBLMAT_MAXDIM], int dim);
void dblmat_identity(double mat[][DBLMAT_MAXDIM], int dim);
void dblmat_copy(double dest [][DBLMAT_MAXDIM],
		 double orig[][DBLMAT_MAXDIM], int dim);
void dblmat_add(double dest[][DBLMAT_MAXDIM], 
		double a[][DBLMAT_MAXDIM],
		double b[][DBLMAT_MAXDIM], int dim);
void dblmat_matmult(double dest[][DBLMAT_MAXDIM], 
		double a[][DBLMAT_MAXDIM],
		double b[][DBLMAT_MAXDIM], int dim);
void dblmat_out(FILE * fp, double a[][DBLMAT_MAXDIM], int dim);

#endif
