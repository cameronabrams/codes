/* 
   Mean-squared displacement calculation, straightforward algorithm.
   Reads in an arbitrary number of configuration snapshots 
   in XYZ format, and computes msd(t).

   Cameron F. Abrams

   Written for the course CHE 800-002, Molecular Simulation
   Spring 0304

   compile using "gcc -o rdf rdf.c -lm"

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2004
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Prints usage information */
void usage ( void ) {
  fprintf(stdout,"msd usage:\n");
  fprintf(stdout,"msd [options]\n\n");
  fprintf(stdout,"Options:\n");
  fprintf(stdout,"\t -fnf [string]\t\tFile name format (%i.xyz)\n");
  fprintf(stdout,"\t -for start,stop,step\t\tloop control\n");
  fprintf(stdout,"\t -mddt [real(0.001)]\t\tMD time step\n");
  fprintf(stdout,"\t -h           \t\tPrint this info.\n");
}

/* Reads in a configuration snapshot in XYZ format (such
   as that created by mdlj.c) */
int xyz_in (FILE * fp, double * rx, double * ry, double * rz, 
	    double * vx, double * vy, double * vz, 
	    int * N) {
  int i;
  int has_vel, dum;
  fscanf(fp,"%i %i\n\n",N,&has_vel);
  
  for (i=0;i<(*N);i++) {
    fscanf(fp,"%i %lf %lf %lf ",&dum,&rx[i],&ry[i],&rz[i]);
    if (has_vel) fscanf(fp, "%lf %lf %lf\n",&vx[i],&vy[i],&vz[i]);
  }
  return has_vel;
}

void com ( double * rx, double * ry, double * rz, int N, 
	   double * cx, double * cy, double * cz) {
  int i;
  *cx=0.0;*cy=0.0;*cz=0.0;
  for (i=0;i<N;i++) {
    (*cx)+=rx[i];
    (*cy)+=ry[i];
    (*cz)+=rz[i];
  }
  (*cx)/=N;  (*cy)/=N;  (*cz)/=N;
}

void traj_in (char * fnf, int start, int stop, int step,
	      double ** rx, double ** ry, double ** rz,
	      double ** vx, double ** vy, double ** vz,
	      int * N, 
	      double * cx, double * cy, double * cz, int * M ) {
  int i,t;
  FILE * fp;
  char fn[50];

  fprintf(stderr,"# Reading in trajectory...");fflush(stderr);
  t=0;
  for (i=start;i<stop+1;i+=step) {
    /* Generate the file name */
    sprintf(fn,fnf,i);
    /* Open the file */
    fp=fopen(fn,"r");
    /* If the file is there... */
    if (fp) {
      /* Read in the data, ignoring velocities */
      xyz_in(fp,rx[t],ry[t],rz[t],vx[t],vy[t],vz[t],N);
      /* Compute the system center of mass from all coordinates */
      com(rx[t],ry[t],rz[t],*N,&cx[t],&cy[t],&cz[t]);
      /* Close the file */
      fclose(fp);
      /* Increment the frame index */
      t++;
    }
  }
  (*M)=t;
  fprintf(stderr,"done: %i frames.\n",t);fflush(stderr);
}

int main ( int argc, char * argv[] ) {

  double ** rx, ** ry, ** rz;
  double ** vx, ** vy, ** vz;
  int N=216, M=0;
  double * sdx, * sdy, * sdz;
  double * cx, * cy, * cz;
  double * sdcx, * sdcy, * sdcz;
  int t, dt, * cnt, tmax=0;
  double md_time_step = 0.001;
  double L=0.0, V=0.0;
  double dr=0.1, r, vb, nid;
  double rc2 = 10.0, rho=0.85;
  int start=0,stop=0,step=1,ngr=0;
  int i;

  char * fnf;
  FILE * fp;

  char fn[35];

  /* Here we parse the command line arguments;  If
   you add an option, document it in the usage() function! */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-for"))
      sscanf(argv[++i],"%i,%i,%i",&start,&stop,&step);
    else if (!strcmp(argv[i],"-fnf")) fnf=argv[++i];
    else if (!strcmp(argv[i],"-mddt")) md_time_step=atof(argv[++i]);
    else if (!strcmp(argv[i],"-tmax")) tmax=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-h")) {
      usage(); exit(0);
    }
    else {
      fprintf(stderr,"Error: Command-line argument '%s' not recognized.\n",
	      argv[i]);
      exit(-1);
    }
  }

  /* Compute the number of frames in the trajectory */
  M = (stop - start) / step + 1;
  if (!tmax) tmax=M;
  fprintf(stderr,"# Expecting %i frames...\n",M);

  /* Allocate the trajectory arrays */
  rx = (double**)malloc(M*sizeof(double*));
  for (i=0;i<M;i++) rx[i] = (double*)malloc(N*sizeof(double));
  ry = (double**)malloc(M*sizeof(double*));
  for (i=0;i<M;i++) ry[i] = (double*)malloc(N*sizeof(double));
  rz = (double**)malloc(M*sizeof(double*));
  for (i=0;i<M;i++) rz[i] = (double*)malloc(N*sizeof(double));
  vx = (double**)malloc(M*sizeof(double*));
  for (i=0;i<M;i++) vx[i] = (double*)malloc(N*sizeof(double));
  vy = (double**)malloc(M*sizeof(double*));
  for (i=0;i<M;i++) vy[i] = (double*)malloc(N*sizeof(double));
  vz = (double**)malloc(M*sizeof(double*));
  for (i=0;i<M;i++) vz[i] = (double*)malloc(N*sizeof(double));

  /* Allocate the system center of mass trajectory */
  cx = (double*)malloc(M*sizeof(double));
  cy = (double*)malloc(M*sizeof(double));
  cz = (double*)malloc(M*sizeof(double));

  /* Read in the trajectory */
  traj_in(fnf,start,stop,step,rx,ry,rz,vx,vy,vz,&N,cx,cy,cz,&M);
  fprintf(stderr,"# Frames read-in: %i\n",M);

  /* Allocate and initialize the squared-displacement teams */
  sdx=(double*)calloc(M,sizeof(double));
  sdy=(double*)calloc(M,sizeof(double));
  sdz=(double*)calloc(M,sizeof(double));
  sdcx=(double*)calloc(M,sizeof(double));
  sdcy=(double*)calloc(M,sizeof(double));
  sdcz=(double*)calloc(M,sizeof(double));
  cnt=(int*)calloc(M,sizeof(int));


  /* Compute the mean-squared displacement using
     the straightforward algorithm */
  fprintf(stderr,"# computing...\n");fflush(stderr);
  for (t=0;t<M;t++) {
    for (dt=1;(t+dt)<tmax;dt++) {
      cnt[dt]++;  /* number of origins for interval length dt  */
      for (i=0;i<N;i++) {
	sdx[dt] += (rx[t+dt][i] - rx[t][i])*(rx[t+dt][i] - rx[t][i]);
	sdy[dt] += (ry[t+dt][i] - ry[t][i])*(ry[t+dt][i] - ry[t][i]);
	sdz[dt] += (rz[t+dt][i] - rz[t][i])*(rz[t+dt][i] - rz[t][i]);
      }
      sdcx[dt] += (cx[t+dt] - cx[t])*(cx[t+dt] - cx[t]);
      sdcy[dt] += (cy[t+dt] - cy[t])*(cy[t+dt] - cy[t]);
      sdcz[dt] += (cz[t+dt] - cz[t])*(cz[t+dt] - cz[t]);
    }
  }
  for (t=0;t<M;t++) {
    sdx[t] /= cnt[t]?(N*cnt[t]):1;
    sdy[t] /= cnt[t]?(N*cnt[t]):1;
    sdz[t] /= cnt[t]?(N*cnt[t]):1;
    sdcx[t] /= cnt[t]?(cnt[t]):1;
    sdcy[t] /= cnt[t]?(cnt[t]):1;
    sdcz[t] /= cnt[t]?(cnt[t]):1;
    fprintf(stdout,"%.5lf %.8lf %.8lf %.8lf %.8lf"
	    " %.8lf %.8lf %.8lf %.8lf %i %.5lf %.5lf %.5lf\n",
	    t*step*md_time_step,sdx[t],sdy[t],sdz[t],
	    sdx[t]+sdy[t]+sdz[t],sdcz[t],sdcy[t],sdcz[t],
	    sdcz[t]+sdcy[t]+sdcz[t],N*cnt[t],cx[t],cy[t],cz[t]);
  }
}
