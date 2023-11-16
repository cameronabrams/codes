#include "rotmov.h"

/* Routines for vector rotation in 3-space.
   Written by Cameron Abrams, 2004
*/


void matrix_identity (double ** R, int dim)
{
  int i,j;
  for (i=0;i<dim;i++) {
    for (j=0;j<dim;j++) {
      R[i][j]=(i==j)?1.0:0.0;
    }
  }
}

/* multiplies square matrices A and B to produce R */
void matrix_multiply (double ** R, double ** A, double ** B, int dim)
{

  int i,j,k;

  for (i=0;i<dim;i++) {
    for (j=0;j<dim;j++) {
      R[i][j]=0.0;
      for (k=0;k<dim;k++)
	R[i][j]+=A[i][k]*B[k][j];
    }
  }

}

/* performs matrix multiplication of vector r_old by A and puts results in r_new */
void matrix_vector_multiply ( double * r_new, double ** A, double * r_old, int dim )
{

  int i,j;

  for (i=0;i<dim;i++) {
    r_new[i]=0.0;
    for (j=0;j<dim;j++) {
      r_new[i]+=A[i][j]*r_old[j];
    }
  }

}

void matrix_make_rotation_z (double ** R, double gam)
{
  double s=sin(gam);
  double c=cos(gam);
  matrix_identity(R,3);
  R[0][0]=c;
  R[0][1]=-s;
  R[1][0]=s;
  R[1][1]=c;
}

void matrix_make_rotation_y (double ** R, double gam)
{
  double s=sin(gam);
  double c=cos(gam);
  matrix_identity(R,3);
  R[0][0]=c;
  R[0][2]=s;
  R[2][0]=-s;
  R[2][2]=c;
}


/* Global rotation matrices declared.  setup_rotation_matrices() should be called from the 
   main program to allocate them. */
static double ** R, ** R1, ** R2, ** R3;
void setup_rotation_matrices (void) {
  int i;
  int dim=3;
  R=(double**)malloc(dim*sizeof(double*));
  R1=(double**)malloc(dim*sizeof(double*));
  R2=(double**)malloc(dim*sizeof(double*));
  R3=(double**)malloc(dim*sizeof(double*));
  for (i=0;i<dim;i++) {
    R[i]=(double*)calloc(dim,sizeof(double));
    R1[i]=(double*)calloc(dim,sizeof(double));
    R2[i]=(double*)calloc(dim,sizeof(double));
    R3[i]=(double*)calloc(dim,sizeof(double));
  }
}

/* 
   crankshaft:  rotates a subset of points (from index a+1 to b-1) around
   an axis defined by the direction of a to b by an angle gam.

   double ** r is an array of 3-component vectors; e.g., the x,y,z-components
               of point n are r[n][0], r[n][1], r[n][2].

   N is the number of points in r;
   a and b are indices in the array of points defining the range to be cranked;
   gam is the angle in radians

   IMPORTANT:  uses three global matrices which must have been previously
   allocated.
   IMPORTANT:  no error checking is done; it is assumed that a and b
   are chosen with 0,N-1.
   
   written by Cameron Abrams
   Drexel University
   2004
*/
enum {NORMAL, REVERSE};
void crankshaft (double ** r, int N, int a, int b, double gam ) {

  int i, sense=NORMAL;
  double phi;
  double theta,cosT,sinT,cosP,sinP,rabs;
  double tv[3], sh[3];

  if (a<b) sense=NORMAL;
  else if (a>b) sense=REVERSE;
  else {
    fprintf(stderr,"Error in crankshaft:  sequence %i=%i is invalid.\n",a,b);
    exit(-1);
  }

  matrix_identity(R,3);
  matrix_identity(R1,3);
  matrix_identity(R2,3);

  /* save position of monomer a */
  sh[0]=r[a][0];
  sh[1]=r[a][1];
  sh[2]=r[a][2];

  /* shift origin of monomers a to b to origin at position a */
  if (sense==NORMAL) {
    for (i=a;i<=b;i++) {
      r[i][0]-=sh[0];
      r[i][1]-=sh[1];
      r[i][2]-=sh[2];
    }
  }
  else {
    for (i=a;i>=b;i--) {
      r[i][0]-=sh[0];
      r[i][1]-=sh[1];
      r[i][2]-=sh[2];
    }
  }

  /* compute polar coordinates of director */
  rabs=sqrt(r[b][0]*r[b][0]+r[b][1]*r[b][1]+r[b][2]*r[b][2]);
  cosT=r[b][2]/rabs;
  theta=acos(cosT);
  sinT=sqrt(1-cosT*cosT);
  cosP=r[b][0]/(rabs*sinT);
  sinP=r[b][1]/(rabs*sinT);
  phi=acos(cosP);
  if (sinP<0) phi=2*M_PI-phi;

  /* construct the rotation matrix */
  matrix_make_rotation_z(R1,-phi);
  matrix_make_rotation_y(R2,-theta);
  matrix_multiply(R,R2,R1,3);
  matrix_make_rotation_z(R1,gam);
  matrix_multiply(R2,R1,R,3);
  matrix_make_rotation_y(R1,theta);
  matrix_multiply(R,R1,R2,3);
  matrix_make_rotation_z(R2,phi);
  matrix_multiply(R1,R2,R,3);

  /* rotate all vectors between a+1 and b-1 inclusive */
  if (sense==NORMAL) {
    for (i=a+1;i<b;i++) {
      matrix_vector_multiply(tv,R1,r[i],3);
      r[i][0]=tv[0];  r[i][1]=tv[1];  r[i][2]=tv[2];
    }
  }
  else {
    for (i=a-1;i>b;i--) {
      matrix_vector_multiply(tv,R1,r[i],3);
      r[i][0]=tv[0];  r[i][1]=tv[1];  r[i][2]=tv[2];
    }
  }

  /* shift all points back to global origin */
  if (sense==NORMAL) {
    for (i=a;i<=b;i++) {
      r[i][0]+=sh[0];
      r[i][1]+=sh[1];
      r[i][2]+=sh[2];
    }
  }
  else {
    for (i=a;i>=b;i--) {
      r[i][0]+=sh[0];
      r[i][1]+=sh[1];
      r[i][2]+=sh[2];
    }
  }      
}

/*
  generalized_pivot:  pivots a subset of points (b+1 to c) and translates points
  (c+1 to N-1) by the amount by which c is translated.  The "hinge" for the pivot
  is defined at point b, and the axis is the vector pointing from b to a.  Two angles
  define the move:  g1 measures how much the hinge bends, and g2 measures a rotation
  of all points between b and c rotated azimuthally around the axis.

   double ** r is an array of 3-component vectors; e.g., the x,y,z-components
               of point n are r[n][0], r[n][1], r[n][2].

   N is the number of points in r;
   a and b are indices in the array of points defining the range to be cranked;
   gam is the angle in radians

   IMPORTANT:  uses three global matrices which must have been previously
   allocated.
   IMPORTANT:  no error checking is done; it is assumed that a and b
   are chosen with 0,N-1.
   
   written by Cameron Abrams
   Drexel University
   2004
*/
void generalized_pivot ( double ** r, int N, int a, int b, int c, double g1, double g2 ) 
{

  int i,sense=NORMAL;
  double phi0,phi1;
  double theta0,theta1,cosT,sinT,cosP,sinP,rabs;
  double tv[3], sh[3], shc[3];

  if (a < b && b < c) {
    sense=NORMAL;
  }
  else if (a > b && b > c) {
    sense=REVERSE;
  }
  else {
    fprintf(stderr,"ERROR in generalized_pivot:  sequence %i=%i=%i is not valid.\n",
	    a,b,c);
    exit(-1);
  }

  matrix_identity(R,3);
  matrix_identity(R1,3);
  matrix_identity(R2,3);

  /* save position of monomers b and c */
  sh[0]=r[b][0];
  sh[1]=r[b][1];
  sh[2]=r[b][2];

  shc[0]=r[c][0];
  shc[1]=r[c][1];
  shc[2]=r[c][2];

  /* shift origin of monomers a to c to origin at position b */
  if (sense==NORMAL) {
    for (i=a;i<=c;i++) {
      r[i][0]-=sh[0];
      r[i][1]-=sh[1];
      r[i][2]-=sh[2];
    }
  }
  else {
    for (i=a;i>=c;i--) {
      r[i][0]-=sh[0];
      r[i][1]-=sh[1];
      r[i][2]-=sh[2];
    }
  }

  /* the reference for this move is r(b->a) */
  /* compute polar coordinates of reference: theta0, phi0 */
  rabs=sqrt(r[a][0]*r[a][0]+r[a][1]*r[a][1]+r[a][2]*r[a][2]);
  cosT=r[a][2]/rabs;
  theta0=acos(cosT);
  sinT=sqrt(1-cosT*cosT);
  cosP=r[a][0]/(rabs*sinT);
  sinP=r[a][1]/(rabs*sinT);
  phi0=acos(cosP);
  if (sinP<0) phi0=2*M_PI-phi0;

  /* construct the rotation matrix */
  matrix_make_rotation_z(R1,M_PI-phi0);     // rotate reference into xz along -x
  matrix_make_rotation_y(R2,theta0-M_PI);   // align local -z along reference
  matrix_multiply(R,R2,R1,3);

  /* rotate all vectors between a and c, including c */
  if (sense==NORMAL) {
    for (i=b+1;i<=c;i++) {
      matrix_vector_multiply(tv,R,r[i],3);
      r[i][0]=tv[0];  r[i][1]=tv[1];  r[i][2]=tv[2];
    }
  }
  else {
    for (i=b-1;i>=c;i--) {
      matrix_vector_multiply(tv,R,r[i],3);
      r[i][0]=tv[0];  r[i][1]=tv[1];  r[i][2]=tv[2];
    }
  }
  
  /* the actor for this move is r(b->c) */
  /* compute polar coordinates of aligned actor: theta1,phi1 */
  rabs=sqrt(r[c][0]*r[c][0]+r[c][1]*r[c][1]+r[c][2]*r[c][2]);
  cosT=r[c][2]/rabs;
  theta1=acos(cosT);
  sinT=sqrt(1-cosT*cosT);
  cosP=r[c][0]/(rabs*sinT);
  sinP=r[c][1]/(rabs*sinT);
  phi1=acos(cosP);
  if (sinP<0) phi1=2*M_PI-phi1;

  matrix_identity(R,3);
  matrix_make_rotation_z(R1,-phi1);  // rotate actor into xz
  matrix_multiply(R2,R1,R,3);
  matrix_make_rotation_y(R1,g1);            // bend the hinge
  matrix_multiply(R,R1,R2,3);
  //matrix_identity(R2,3);
  matrix_make_rotation_z(R2,g2);            // rotate around local z by g2
  matrix_multiply(R1,R2,R,3);
  matrix_make_rotation_z(R2,phi1);   // rotate actor out of xz
  matrix_multiply(R,R2,R1,3);

  /* rotate */
  if (sense==NORMAL) {
    for (i=b+1;i<=c;i++) {
      matrix_vector_multiply(tv,R,r[i],3);
      r[i][0]=tv[0];  r[i][1]=tv[1];  r[i][2]=tv[2];
    }
  }
  else {
    for (i=b-1;i>=c;i--) {
      matrix_vector_multiply(tv,R,r[i],3);
      r[i][0]=tv[0];  r[i][1]=tv[1];  r[i][2]=tv[2];
    }
  }

  matrix_identity(R,3);
  matrix_make_rotation_y(R2,M_PI-theta0);    // unalign the reference
  matrix_multiply(R1,R2,R,3);
  matrix_make_rotation_z(R2,phi0-M_PI);           // rotate reference out of xz
  matrix_multiply(R,R2,R1,3);

  /* rotate */
  if (sense==NORMAL) {
    for (i=b+1;i<=c;i++) {
      matrix_vector_multiply(tv,R,r[i],3);
      r[i][0]=tv[0];  r[i][1]=tv[1];  r[i][2]=tv[2];
    }
  }
  else {
    for (i=b-1;i>=c;i--) {
      matrix_vector_multiply(tv,R,r[i],3);
      r[i][0]=tv[0];  r[i][1]=tv[1];  r[i][2]=tv[2];
    }
  }

  /* shift all back to global origin */
  if (sense==NORMAL) {
    for (i=a;i<=c;i++) {
      r[i][0]+=sh[0];
      r[i][1]+=sh[1];
      r[i][2]+=sh[2];
    }
  }
  else {
    for (i=a;i>=c;i--) {
      r[i][0]+=sh[0];
      r[i][1]+=sh[1];
      r[i][2]+=sh[2];
    }
  }

  /* compute shift of points above c */
  sh[0]=r[c][0]-shc[0];
  sh[1]=r[c][1]-shc[1];
  sh[2]=r[c][2]-shc[2];

  /* shift monomers of points above c */
  if (sense==NORMAL) {
    for (i=c+1;i<N;i++) {
      r[i][0]+=sh[0];
      r[i][1]+=sh[1];
      r[i][2]+=sh[2];
    }
  }
  else {
    for (i=c-1;i>=0;i--) {
      r[i][0]+=sh[0];
      r[i][1]+=sh[1];
      r[i][2]+=sh[2];
    }
  }
}
