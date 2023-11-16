#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/* 
   Routines for crankshaft and generalized pivot moves 
   in 3-space.

   Written by Cameron Abrams, 21-Oct-2004
   Drexel University
   Philadelphia
*/

/* Global rotation matrices declared.  setup_rotation_matrices() should be called from the 
   main program to allocate them. */
void setup_rotation_matrices ( void );

/* 
   crankshaft:  rotates a subset of points around
   an axis defined by the direction of a to b by an angle gam.
   The subset of points are those between indices a and b in the array of points.

   double ** r is an array of 3-component vectors, each representing a point in 3-space; 
               e.g., the x,y,z-components of point n are r[n][0], r[n][1], r[n][2].

   N is the number of points in r;
   a and b are indices in the array of points defining the range to be cranked;
   gam is the angle in radians.

   The order of a and b is irrelevant; a may be greater than or less than b.
   If a and b are equal, then an error results.

   IMPORTANT:  It is assumed that both a and b are between 0 and N-1!

   IMPORTANT:  This function uses three global matrices (R,R1,R2) 
   which must be allocated by a call to setup_rotation_matrices.
   
   written by Cameron Abrams
   Drexel University
   2004
*/
void crankshaft ( double ** r, int N, int a, int b, double gam );

/*
  generalized_pivot:  pivots a subset of points.  The "hinge" for the pivot
  is defined at point b, and the reference axis is the vector pointing from b to a.
  The hinge angle is defined by the vector from b to a and the vector from b to c.
  Two angles define the move:  g1 measures how much the hinge bends, and g2 a rotation
  of all points between b and c rotated azimuthally around the axis.

   double ** r is an array of 3-component vectors; e.g., the x,y,z-components
               of point n are r[n][0], r[n][1], r[n][2].

   N is the number of points in r;

   b is the index of the "hinge point" in the array r[].
   a is the index of the "reference point" in r[].
   c is the index of the "endpoint" in r[].
   
  

   IMPORTANT:  uses three global matrices which must have been previously
   allocated.
   IMPORTANT:  no error checking is done; it is assumed that a and b
   are chosen with 0,N-1.
   
   written by Cameron Abrams
   Drexel University
   2004
*/
void generalized_pivot ( double ** r, int N, int a, int b, int c, double g1, double g2 );
