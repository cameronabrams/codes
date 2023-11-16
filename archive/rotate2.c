#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/* Routines for vector rotation in 3-space.
   Written by Cameron Abrams, 2004
*/


void matrix_identity (double ** R, int dim)
{
  int i,j;
  for (i=0;i<dim;i++) {
    for (j=0;j<dim;j++) {
      R[i][j]=i==j?1.0:0.0;
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
void crankshaft (double ** r, int N, int a, int b, double gam ) {

  int i;
  double dTheta,phi;
  double theta,cosT,sinT,cosP,sinP,rabs;
  double tv[3], sh[3];

  matrix_identity(R,3);
  matrix_identity(R1,3);
  matrix_identity(R2,3);

  /* save position of monomer a */
  sh[0]=r[a][0];
  sh[1]=r[a][1];
  sh[2]=r[a][2];

  /* shift origin of monomers a to b to origin at position a*/
  for (i=a;i<b;i++) {
    r[i][0]-=sh[0];
    r[i][1]-=sh[1];
    r[i][2]-=sh[2];
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
  for (i=a+1;i<b;i++) {
    matrix_vector_multiply(tv,R1,r[i],3);
    r[i][0]=tv[0];  r[i][1]=tv[1];  r[i][2]=tv[2];
  }  

  /* shift all points back to global origin */
  for (i=a;i<b;i++) {
    r[i][0]+=sh[0];
    r[i][1]+=sh[1];
    r[i][2]+=sh[2];
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

  int i;
  double dTheta,phi;
  double theta,cosT,sinT,cosP,sinP,rabs;
  double tv[3], sh[3], shc[3];

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
  for (i=a;i<=c;i++) {
    r[i][0]-=sh[0];
    r[i][1]-=sh[1];
    r[i][2]-=sh[2];
  }

  /* the director for this move is r(b->a) */
  /* compute polar coordinates of director */
  rabs=sqrt(r[a][0]*r[a][0]+r[a][1]*r[a][1]+r[a][2]*r[a][2]);
  cosT=r[a][2]/rabs;
  theta=acos(cosT);
  sinT=sqrt(1-cosT*cosT);
  cosP=r[a][0]/(rabs*sinT);
  sinP=r[a][1]/(rabs*sinT);
  phi=acos(cosP);
  if (sinP<0) phi=2*M_PI-phi;
  
  /* construct the rotation matrix */
  matrix_make_rotation_z(R1,-phi);
  matrix_make_rotation_y(R2,-theta);
  matrix_multiply(R,R2,R1,3);
  matrix_make_rotation_z(R1,g2);
  matrix_multiply(R2,R1,R,3);
  matrix_make_rotation_y(R1,g1);
  matrix_multiply(R,R1,R2,3);
  matrix_make_rotation_z(R2,phi);
  matrix_multiply(R1,R2,R,3);

  /* rotate all vectors between a and b */
  for (i=b+1;i<=c;i++) {
    matrix_vector_multiply(tv,R1,r[i],3);
    r[i][0]=tv[0];  r[i][1]=tv[1];  r[i][2]=tv[2];
  }  

  /* shift back to global origin */
  for (i=b;i<=c;i++) {
    r[i][0]+=sh[0];
    r[i][1]+=sh[1];
    r[i][2]+=sh[2];
  }

  /* compute shift of points above c */
  sh[0]=r[c][0]-shc[0];
  sh[1]=r[c][1]-shc[1];
  sh[2]=r[c][2]-shc[2];

  /* shift monomers of points above c */
  for (i=c+1;i<N;i++) {
    r[i][0]+=sh[0];
    r[i][1]+=sh[1];
    r[i][2]+=sh[2];
  }


}

void make_random_walk ( double ** r, int N ) {

  int i;
  double phi, costheta, sintheta, cosphi, sinphi, dx, dy, dz;
  r[0][0]=r[0][1]=r[0][2]=0.0;
  

  for (i=1;i<N;i++) {
    phi=2*M_PI*(double)rand()/RAND_MAX;
    costheta = 1.0 - 2*(double)rand()/RAND_MAX;
    sintheta = sqrt(1-costheta*costheta);
    cosphi=cos(phi);
    sinphi=sin(phi);
    dx=sintheta*cosphi;
    dy=sintheta*sinphi;
    dz=costheta;
    r[i][0]=r[i-1][0]+dx;
    r[i][1]=r[i-1][1]+dy;
    r[i][2]=r[i-1][2]+dz;
  }
}

void xyz_out (FILE * fp, double **r, int N ) {
  
  int i;

  fprintf(fp,"%i 0 0\n\n",N);
  for (i=0;i<N;i++) fprintf(fp,"1 %.4lf %.4lf %.4lf\n", r[i][0], r[i][1], r[i][2]);

}

void xyz_out_highlight_subchain (FILE * fp, double **r, int N, int a, int b ) {
  
  int i,typ;

  fprintf(fp,"%i 0 0\n\n",N);
  for (i=0;i<N;i++) {
    typ=1;
    if (i==a || i==b) typ=5;
    if (i>a && i<b) typ=6;
    fprintf(fp,"%i %.4lf %.4lf %.4lf\n", typ,r[i][0], r[i][1], r[i][2]);
  }
}

int main (int argc, char * argv[]) {

  int i,seed =101;
  int N = 100;
  int a = 79;
  int b = 80;
  int c = 90;
  double ** r;
  FILE * fp;

  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-seed")) seed=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
  }

  srand(seed);

  setup(3);

  r = (double **)malloc(N*sizeof(double*));
  for (i=0;i<N;i++) r[i]=(double*)calloc(3,sizeof(double));

  make_random_walk(r,N);

  fp=fopen("tmp0.xyz","w");
  xyz_out_highlight_subchain(fp,r,N,b,c);
  fclose(fp);

  //crankshaft(r,N,a,b,M_PI);
  generalized_pivot(r,N,a,b,c,M_PI/2,M_PI);

  fp=fopen("tmp1.xyz","w");
  xyz_out_highlight_subchain(fp,r,N,b,c);
  fclose(fp);

}
