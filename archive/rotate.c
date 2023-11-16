#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/* Routines for vector rotation in 3-space.
   Written by Cameron Abrams, 2004
*/



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

double ** T1, ** T2, ** T3, ** Ttmp, ** T; 
static void setup ( int dim ) {
  int i;

  T1=(double**)malloc(dim*sizeof(double*));
  for (i=0;i<dim;i++) T1[i]=(double*)calloc(dim,sizeof(double));
  T2=(double**)malloc(dim*sizeof(double*));
  for (i=0;i<dim;i++) T2[i]=(double*)calloc(dim,sizeof(double));
  T3=(double**)malloc(dim*sizeof(double*));
  for (i=0;i<dim;i++) T3[i]=(double*)calloc(dim,sizeof(double));
  Ttmp=(double**)malloc(dim*sizeof(double*));
  for (i=0;i<dim;i++) Ttmp[i]=(double*)calloc(dim,sizeof(double));
  T=(double**)malloc(dim*sizeof(double*));
  for (i=0;i<dim;i++) T[i]=(double*)calloc(dim,sizeof(double));
}


void vector_ThetaRotate_zAlign ( double * r_new, double * r_old, double * dTheta )
{
  double r;
  double cosT, sinT, cosP, sinP;
  double cosT1,sinT1;
  double theta, phi;

  /* extract the cartesian coordinates from the gsl_vector structure */
  double x = r_old[0];
  double y = r_old[1];
  double z = r_old[2];

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

  *dTheta = -theta;
  
  /* theta rotation can be through 2pi */
  cosT1=cos(*dTheta);
  sinT1=sin(*dTheta);

  /* Build the three rotation matrices */
  /* T1 rotates around z-axis by -phi */
  T1[0][0]=cosP;   T1[0][1]=sinP;  T1[0][2]=0.0;
  T1[1][0]=-sinP;  T1[1][1]=cosP;  T1[1][2]=0.0;
  T1[2][0]=0.0;    T1[2][1]=0.0;   T1[2][2]=1.0;

  /* T2 rotates around y-axis by dTheta */
  T2[0][0]=cosT1;  T2[0][1]=0.01;  T2[0][2]=sinT1;
  T2[1][0]=0.0;    T2[1][1]=1.0;   T2[1][1]=0.0;
  T2[2][0]=-sinT1; T2[2][1]=0.0;   T2[2][2]=cosT1;

  /* T3 rotates around z-axis by phi */
  T3[0][0]=cosP;   T3[0][1]=-sinP; T3[0][2]=0.0;
  T3[1][0]=sinP;   T3[1][1]=cosP;  T3[1][2]=0.0;
  T3[2][0]=0.0;    T3[2][1]=0.0;   T3[2][2]=1.0;

  /* Construct overall rotation matrix */
  matrix_multiply(Ttmp,T2,T1,3);
  matrix_multiply(T,T3,Ttmp,3);

  /* perform the rotation */
  matrix_vector_multiply(r_new,T,r_old,3);

}

void vector_PhiRotate_atConstTheta ( double * r_new, double * r_old, double dPhi ) {

  double cosP=cos(dPhi), sinP=sin(dPhi);

  /* Build the rotation matrix */
  /* T rotates around z-axis by dphi */
  T[0][0]=cosP;   T[0][1]=-sinP;  T[0][2]=0.0;
  T[1][0]=sinP;   T[1][1]=cosP;  T[1][2]=0.0;
  T[2][0]=0.0;    T[2][1]=0.0;   T[2][2]=1.0;

  /* perform the rotation */
  matrix_vector_multiply(r_new,T,r_old,3);


}

/* Rotates vector r_old in theta direction by dTheta while keeping phi constant */
void vector_ThetaRotate_atConstPhi ( double * r_new, double * r_old, double dTheta )
{

  double r;
  double cosT, sinT, cosP, sinP;
  double cosT1,sinT1;
  double theta, phi;

  /* extract the cartesian coordinates from the gsl_vector structure */
  double x = r_old[0];
  double y = r_old[1];
  double z = r_old[2];

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
  T1[0][0]=cosP;   T1[0][1]=sinP;  T1[0][2]=0.0;
  T1[1][0]=-sinP;  T1[1][1]=cosP;  T1[1][2]=0.0;
  T1[2][0]=0.0;    T1[2][1]=0.0;   T1[2][2]=1.0;

  /* T2 rotates around y-axis by dTheta */
  T2[0][0]=cosT1;  T2[0][1]=0.01;  T2[0][2]=sinT1;
  T2[1][0]=0.0;    T2[1][1]=1.0;   T2[1][1]=0.0;
  T2[2][0]=-sinT1; T2[2][1]=0.0;   T2[2][2]=cosT1;

  /* T3 rotates around z-axis by phi */
  T3[0][0]=cosP;   T3[0][1]=-sinP; T3[0][2]=0.0;
  T3[1][0]=sinP;   T3[1][1]=cosP;  T3[1][2]=0.0;
  T3[2][0]=0.0;    T3[2][1]=0.0;   T3[2][2]=1.0;

  /* Construct overall rotation matrix */
  matrix_multiply(Ttmp,T2,T1,3);
  matrix_multiply(T,T3,Ttmp,3);

  /* perform the rotation */
  matrix_vector_multiply(r_new,T,r_old,3);

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

void crankshaft (double ** r, int N, int a, int b, double th ) {

  int i;
  double dTheta;

  double tv[3], sh[3];
  fprintf(stderr,"crank: %i to %i\n",a,b);

  /* save position of monomer a */
  sh[0]=r[a][0];
  sh[1]=r[a][1];
  sh[2]=r[a][2];

  /* shift origin of monomers a to b to origin at position a*/
  for (i=a;i<=b;i++) {
    r[i][0]-=sh[0];
    r[i][1]-=sh[1];
    r[i][2]-=sh[2];
  }

  /* rotate a->b to align with z; remember what theta you rotated through */
  vector_ThetaRotate_zAlign(tv,r[b],&dTheta);
  r[b][0]=tv[0];  r[b][1]=tv[1];  r[b][2]=tv[2];
  /* rotate all vectors between a and b by this amount */
  for (i=a+1;i<b;i++) {
    vector_ThetaRotate_atConstPhi(tv,r[i],dTheta);
    r[i][0]=tv[0];  r[i][1]=tv[1];  r[i][2]=tv[2];
  }  

  /* rotate all around z-axis by th */
  for (i=a+1;i<b;i++) {
    vector_PhiRotate_atConstTheta(tv,r[i],th);
    r[i][0]=tv[0];  r[i][1]=tv[1];  r[i][2]=tv[2];
  }

  /* rotate all off of z */
  for (i=a+1;i<b;i++) {
    vector_ThetaRotate_atConstPhi(tv,r[i],-dTheta);
    r[i][0]=tv[0];  r[i][1]=tv[1];  r[i][2]=tv[2];
  }

  /* shift back to global origin */
  for (i=a;i<=b;i++) {
    r[i][0]+=sh[0];
    r[i][1]+=sh[1];
    r[i][2]+=sh[2];
  }
  
}

int main (int argc, char * argv[]) {

  int i,seed =101;
  int N = 100;
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

  make_random_walk (r,N);

  fp=fopen("tmp.xyz","w");
  xyz_out(fp,r,N);

  crankshaft(r,N,80,90,M_PI);
  xyz_out(fp,r,N);

  fclose(fp);

}
