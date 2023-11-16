#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


/* Routines for vector rotation in 3-space.
   Uses GSL vectors and matrices.
   Written by Cameron Abrams, 2004
*/

void gsl_matrix_multiply (gsl_matrix * R, gsl_matrix * A, gsl_matrix * B)
{

  int i,j,k;

  /* check dimensions */
  if (R->size1 != A->size1 || R->size2 != B->size2) {
    fprintf(stderr,"gsl_matrix_multiply:  results is wrongly dimensioned:\n"
	    "expecting %ix%i, got %ix%i\n",A->size1,B->size2,
	    R->size1,R->size2);
    exit(-1);
  }
  if (A->size2 != B->size1) {
    fprintf(stderr,"gsl_matrix_multiply:\n\tNumber of columns in A (%i) "
	    "does not match number of rows in B (%i)\n",
	    A->size2,B->size1);
    exit(-1);
  }

  for (i=0;i<R->size1;i++) {
    for (j=0;j<R->size2;j++) {
      gsl_matrix_set(R,i,j,0.0);
      for (k=0;k<A->size1;k++)
	gsl_matrix_set(R,i,j,(gsl_matrix_get(R,i,j)+
			     (gsl_matrix_get(A,i,k)*gsl_matrix_get(B,k,j))));
    }
  }

}

void gsl_matrix_vector_multiply ( gsl_vector * r_new, gsl_matrix * T, gsl_vector * r_old )
{

  int i,j;

  /* check dimensions */
  if (r_new->size != r_old->size && T->size1 != T->size2) {
    fprintf(stderr,"gsl_matrix_vector_multiply: source and result vectors have different sizes.\n");
    exit(-1);
  }
  if (r_new->size != T->size2) {
    fprintf(stderr,"gsl_matrix_vector_multiply: Number of columns of matrix (%i) does "
	    "not match number of rows of vector (%i).\n",T->size2,r_new->size);
    exit(-1);
  }

  for (i=0;i<r_new->size;i++) {
    r_new->data[i]=0.0;
    for (j=0;j<T->size2;j++) {
      r_new->data[i]+=gsl_matrix_get(T,i,j)*r_old->data[j];
    }
  }

}

/* Matrices used in function; must be allocated in main pgm */
gsl_matrix * T1, * T2, * T3, * T, * Ttmp;
void gsl_vector_ThetaRotate_constPhi ( gsl_vector * r_new, gsl_vector * r_old, double dTheta )
{

  double r;
  double cosT, sinT, cosP, sinP;
  double cosT1,sinT1;
  double theta, phi;

  /* extract the cartesian coordinates from the gsl_vector structure */
  double x = gsl_vector_get(r_old,0);
  double y = gsl_vector_get(r_old,1);
  double z = gsl_vector_get(r_old,2);

  /* compute the polar coordinates */
  /* polar coordinates for given has 0 < theta < pi */
  r = sqrt(x*x+y*y+z*z);
  cosT=z/r;
  sinT=sqrt(1-cosT*cosT);
  cosP=x/(r*sinT);
  sinP=y/(r*sinT);
  theta = acos(cosT);
  phi = acos(cosP);
  if (sinP<0.0) phi = 2*M_PI-phi;

  /* theta rotation can be through 2pi */
  cosT1=cos(dTheta);
  sinT1=sin(dTheta);


  /* Build the three rotation matrices */
  /* T1 rotates around z-axis by -phi */
  gsl_matrix_set(T1,0,0,cosP);
  gsl_matrix_set(T1,0,1,sinP);
  gsl_matrix_set(T1,0,2,0.0);
  gsl_matrix_set(T1,1,0,-sinP);
  gsl_matrix_set(T1,1,1,cosP);
  gsl_matrix_set(T1,1,2,0.0);
  gsl_matrix_set(T1,2,0,0.0);
  gsl_matrix_set(T1,2,1,0.0);
  gsl_matrix_set(T1,2,2,1.0);

  /* T2 rotates around y-axis by dTheta */
  gsl_matrix_set(T2,0,0,cosT1);
  gsl_matrix_set(T2,0,1,0.0);
  gsl_matrix_set(T2,0,2,sinT1);
  gsl_matrix_set(T2,1,0,0.0);
  gsl_matrix_set(T2,1,1,1.0);
  gsl_matrix_set(T2,1,2,0.0);
  gsl_matrix_set(T2,2,0,-sinT1);
  gsl_matrix_set(T2,2,1,0.0);
  gsl_matrix_set(T2,2,2,cosT1);

  /* T3 rotates around z-axis by phi */
  gsl_matrix_set(T3,0,0,cosP);
  gsl_matrix_set(T3,0,1,-sinP);
  gsl_matrix_set(T3,0,2,0.0);
  gsl_matrix_set(T3,1,0,sinP);
  gsl_matrix_set(T3,1,1,cosP);
  gsl_matrix_set(T3,1,2,0.0);
  gsl_matrix_set(T3,2,0,0.0);
  gsl_matrix_set(T3,2,1,0.0);
  gsl_matrix_set(T3,2,2,1.0);

  /* Construct overall rotation matrix */
  gsl_matrix_multiply(Ttmp,T2,T1);
  gsl_matrix_multiply(T,T3,Ttmp);

  /* perform the rotation */
  gsl_matrix_vector_multiply(r_new,T,r_old);

}


void setup ( int dim ) {

  T1=gsl_matrix_calloc(dim,dim);
  T2=gsl_matrix_calloc(dim,dim);
  T3=gsl_matrix_calloc(dim,dim);
  Ttmp=gsl_matrix_calloc(dim,dim);
  T=gsl_matrix_calloc(dim,dim);

}

int main (int argc, char * argv[]) {

  int i;
  gsl_vector * r_new = gsl_vector_calloc(3);
  gsl_vector * r_old = gsl_vector_calloc(3);

  double th=0.0;

  setup(3);

  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-v")) {
      sscanf(argv[++i],"%lf,%lf,%lf",&r_old->data[0],&r_old->data[1],&r_old->data[2]);
    }
    else if (!strcmp(argv[i],"-t")) th=atof(argv[++i]);
    else if (!strcmp(argv[i],"-td")) th=M_PI/180.0*atof(argv[++i]);
  }

  gsl_vector_ThetaRotate_constPhi(r_new,r_old,th);

  fprintf(stdout,"# original: %.4lf %.4lf %.4lf (%.4lf); rotate theta by %.2lf rad:\n",
	  r_old->data[0],r_old->data[1],r_old->data[2],
	  sqrt(r_old->data[0]*r_old->data[0]+r_old->data[1]*r_old->data[1]+r_old->data[2]*r_old->data[2]),th);
  fprintf(stdout,"%.4lf %.4lf %.4lf (%.4lf)\n",
	  r_new->data[0],r_new->data[1],r_new->data[2],
	  sqrt(r_new->data[0]*r_new->data[0]+r_new->data[1]*r_new->data[1]+r_new->data[2]*r_new->data[2]));

}
