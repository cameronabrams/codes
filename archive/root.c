#include <stdio.h>
#include <string.h>
#include <math.h>

#define JMAX 100

float rtnewt(void (*funcd)(float, float *, float *), float x1, float x2,
        float xacc)
{
        void nrerror(char error_text[]);
        int j;
        float df,dx,f,rtn;

        rtn=0.5*(x1+x2);
        for (j=1;j<=JMAX;j++) {
                (*funcd)(rtn,&f,&df);
                dx=f/df;
                rtn -= dx;
                if ((x1-rtn)*(rtn-x2) < 0.0)
                        nrerror("Jumped out of brackets in rtnewt");
                if (fabs(dx) < xacc) return rtn;
        }
        nrerror("Maximum number of iterations exceeded in rtnewt");
        return 0.0;
}

float a=1.0,b=1.0,v=1.0;
void yukr12 ( float x, float * f, float * dfdx )
{
  float x12=x*x*x*x*x*x;
  x12*=x12;
  x12=1.0/x12;
  *f=4.0*x12-a/x*exp(-b*x)-v;
  *dfdx=-48.0*x12/x+(1+b*x)*exp(-b*x)*a/(x*x);
}

void r126 ( float x, float * f, float * dfdx )
{
  float x6=x*x*x*x*x*x;
  float x12=x6*x6;
  x6=1.0/x6;
  x12=1.0/x12;
  *f=x<1.12?4.0*x12-4.0*x6+1.0-v:-v;
  *dfdx=-48.0*x12/x+24.0*x6/x;
}

void polynom ( float x, float * f, float * dfdx )
{
   float x5=x*x*x*x*x;
   float x6=x5*x;
   float x9=x6*x*x*x;
   float x10=x5*x5;
   float a=(1./(2*M_PI)-0.6);

   *f=a*x10+x6-0.4;
   *dfdx=10*a*x9+6*x5;
}
   

float e; /* user-selected well depth */
float rm=1.122462; /* 2^(1/6) */
void b12_e ( float x, float * f, float * dfdx )
{
  float r12=rm*rm*rm*rm,r11;
  float ermc,A,B;
  r12=r12*r12*r12;
  r11=r12/rm;
  ermc=e*r12-1;
  A=r11*(ermc*(1+x*rm)+12);
  B=(ermc+1);
  *f=A*exp(-x*rm)-B*12*exp(-x);
  *dfdx=(-rm*A+r12*ermc)*exp(-x*rm)+B*12*exp(-x);
}

void usage ( char * huh )
{
  fprintf(stderr,"`%s'? read the source\n",huh);
}

int main ( int argc, char * argv[])
{
  float xacc=1.e-6,x,y,d,a,c;
  float x1=0.0, x2=1.0;
  float res=0.0;
  char * fmt="%.5f\n";
  short dat=0;
  int i=0;

  e=-1.0;
  for (i=1;i<argc;i++)
  {
    if (!strcmp(argv[i],"-a")) a=atof(argv[++i]);
    else if (!strcmp(argv[i],"-b")) b=atof(argv[++i]);
    else if (!strcmp(argv[i],"-v")) v=atof(argv[++i]);
    else if (!strcmp(argv[i],"-rm")) rm=atof(argv[++i]);
    else if (!strcmp(argv[i],"-e")) e=atof(argv[++i]);
    else if (!strcmp(argv[i],"-x2")) x2=atof(argv[++i]);
    else if (!strcmp(argv[i],"-x1")) x1=atof(argv[++i]);
    else if (!strcmp(argv[i],"-xacc")) xacc=atof(argv[++i]);
    else if (!strcmp(argv[i],"-fmt")) fmt=argv[++i];
    else if (!strcmp(argv[i],"-dat")) dat=1;
    else usage(argv[i]);
  }

  //fprintf(stderr,"Potential root finder\n");
  //fprintf(stderr,"a=%.3f b=%.3f v=%.3f\n",a,b,v);
  //fprintf(stderr,"x1=%.3f x2=%.3f xacc=%.3e\n",x1,x2,xacc);
//  if (a!=1.0) res=rtnewt(&yukr12,x1,x2,xacc);
//  else res=rtnewt(&r126,x1,x2,xacc); 
  if (!dat)
  {
   res=rtnewt(&polynom,x1,x2,xacc);
   printf("root = %.5f\n",res);
//   fprintf(stdout,fmt,res);
   //   a=12.0/(pow(rm,11)*exp(-res*rm)*(1+res*rm) - 12*exp(-res));
   //   c=1+a*exp(-res);
   //printf("# e=%.5f rm=%.5f --> c=%.5f a=%.5f b=%.5f\n",
   //  e,rm,c,a,res);
   //for (x=0.8;x<2.51;x+=0.01) 
   //printf("%.5f %.5f\n",x,c/pow(x,12)-a/x*exp(-res*x));
  }

  if (dat) 
   for(x=0.0;x<100;x+=0.1) 
   {
    polynom(x,&y,&d);
    fprintf(stdout,"%.5e %.5e %.5e\n",x,y,d);
   }

}
