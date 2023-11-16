#include <stdio.h>
#include <math.h>
#define SDIM 3
typedef float _fl_;
typedef _fl_ TVector[SDIM];
typedef _fl_ TMatrix [SDIM][SDIM];
void mycross ( TVector r, TVector a, TVector b, int dim )
{
  /*uint i;
  for (i=0;i<dim;i++)
    r[i]=a[(i+1)%dim]*b[(i+2)%dim]-b[(i+1)%dim]*a[(i+2)%dim];
  */
  r[0]=a[1]*b[2]-b[1]*a[2];
  r[1]=a[2]*b[0]-b[2]*a[0];
  r[2]=a[0]*b[1]-b[0]*a[1];
}

void SMult (TVector a, float x, uint dim)
{
  int i;
  for (i=0;i<dim;i++) a[i]*=x;
}

float ScalarProd(TVector a, TVector b, uint sdim)
{
  int i;
  float s=0.0;
  for (i=0;i<sdim;i++) s+=a[i]*b[i];
  return s;
}

float NormSqr(TVector a, uint sdim)
{
  int i;
  float s=0.0;
  for (i=0;i<sdim;i++) s+=a[i]*a[i];
  return s;  
}

float Norm (TVector a, uint sdim)
{
  return sqrt(NormSqr(a,sdim));
}

float tor ( TVector b[3], int dim )
{
  TVector c1,c2;
  uint ref;
  mycross(c1,b[0],b[1],dim);
  mycross(c2,b[1],b[2],dim);
  SMult(c1,1.0/Norm(c1,dim),dim);
  SMult(c2,1.0/Norm(c2,dim),dim);
  ref=(ScalarProd(b[2],c1,dim)<0)?0:1;
  return ref*360.0 + (ref?-1:1)*180/M_PI*acos(ScalarProd(c1,c2,dim));
}
void makeRotMat ( TMatrix A, _fl_ ctheta, _fl_ phi )
{
  A[2][2] = -ctheta;                   /* rotate z-axis through (pi-theta) */
  A[2][0] = -sqrt(1 - ctheta*ctheta);  /* -sinTheta */
  A[0][1] = sin(phi);
  A[1][1] =  cos(phi);
  A[2][1] =  0;
  A[0][0] =  A[1][1]*A[2][2];
  A[0][2] = -A[1][1]*A[2][0];
  A[1][0] = -A[0][1]*A[2][2];
  A[1][2] =  A[0][1]*A[2][0];
}

void matMult ( TMatrix D, TMatrix A, TMatrix B, int dim)
{
  _fl_ s;
  int i,j,k;
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      s = 0.0;
      for (k = 0; k < dim; k++) {
	s += A[i][k]*B[k][j];
      }
      D[i][j]=s;
    }
  }
}

void matCopy ( TMatrix D, TMatrix A, int dim)
{
  int i,j;
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      D[i][j]=A[i][j];
    }
  }
}

void matPrint ( TMatrix A, int dim )
{
  int i,j;
  for (i=0;i<dim;i++) {
    for (j=0;j<dim;j++) printf("%.0f ",A[i][j]);
    printf("\n");
  }
}

int main (int argc, char * argv[])
{
  _fl_ phi,theta;
  TVector b[3],r[4];
  TMatrix R,Th,T={{1,0,0},{0,1,0},{0,0,1}};
  int i;

  r[0][0]=0; r[0][1]=0; r[0][2]=1;
  phi=0;
  theta=0;
  
  makeRotMat(R,cos(theta),phi);
  matMult(Th,T,R,SDIM);
  matCopy(T,Th,SDIM);
  
  printf("R2)\n");
  matPrint(R,SDIM);
  printf("T2)\n");
  matPrint(T,SDIM);

  r[1][0] = r[0][0]+T[0][2];
  r[1][1] = r[0][1]+T[1][2];
  r[1][2] = r[0][2]+T[2][2];

  phi=0;
  theta=M_PI_2;

  makeRotMat(R,cos(theta),phi);
  matMult(Th,T,R,SDIM);
  matCopy(T,Th,SDIM);
  
  printf("R3)\n");
  matPrint(R,SDIM);
  printf("T3)\n");
  matPrint(T,SDIM);

  r[2][0] = r[1][0]+T[0][2];
  r[2][1] = r[1][1]+T[1][2];
  r[2][2] = r[1][2]+T[2][2];

  phi=M_PI_2;
  theta=M_PI_2;

  makeRotMat(R,cos(theta),phi);
  matMult(Th,T,R,SDIM);
  matCopy(T,Th,SDIM);
  
  printf("R4)\n");
  matPrint(R,SDIM);
  printf("T4)\n");
  matPrint(T,SDIM);

  r[3][0] = r[2][0]+T[0][2];
  r[3][1] = r[2][1]+T[1][2];
  r[3][2] = r[2][2]+T[2][2];

  for (i=0;i<3;i++) {
    b[i][0]=r[i+1][0]-r[i][0];
    b[i][1]=r[i+1][1]-r[i][1];
    b[i][2]=r[i+1][2]-r[i][2];
  }

  for (i=0;i<3;i++) fprintf(stdout,"%i) %.2f %.2f %.2f\n",i,b[i][0],b[i][1],b[i][2]);
  fprintf(stdout,"%.5f\n",tor(b,SDIM));

}
