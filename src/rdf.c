/* 
   Radial distribution function calculation.
   Reads in an arbitrary number of configuration snapshots 
   in XYZ format, and computes g(r).

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
  fprintf(stdout,"rdf usage:\n");
  fprintf(stdout,"rdf [options]\n\n");
  fprintf(stdout,"Options:\n");
  fprintf(stdout,"\t -fnf [string]\t\tFile name format (%i.xyz)\n");
  fprintf(stdout,"\t -for start,stop,step\t\tloop control\n");
  fprintf(stdout,"\t -dr [real]\t\tRadial binsize (0.1)\n");
  fprintf(stdout,"\t -rc [real]\t\tCutoff radius\n");
  fprintf(stdout,"\t -rho [real]\t\tNumber density (0.85)\n");
  fprintf(stdout,"\t -h           \t\tPrint this info.\n");
}

/* Reads in a configuration snapshot in XYZ format (such
   as that created by mdlj.c) */
int xyz_in (FILE * fp, double * rx, double * ry, double * rz, 
	     int * N) {
  int i;
  int has_vel, dum;
  double vx,vy,vz;
  fscanf(fp,"%i %i\n\n",N,&has_vel);
  
  for (i=0;i<(*N);i++) {
    fscanf(fp,"%i %lf %lf %lf ",&dum,&rx[i],&ry[i],&rz[i]);
    if (has_vel) fscanf(fp, "%lf %lf %lf\n",&vx,&vy,&vz);
  }
  return has_vel;
}

/* An N^2 algorithm for computing interparticle separations
   and updating the radial distribution function histogram. 
   Implements parts of Algorithm 7 in F&S. */
void update_hist ( double * rx, double * ry, double * rz, 
		   int N, double L, double rc2, double dr, int * H ) {
   int i,j;
   double dx, dy, dz, r2;
   int bin;
   double hL = L/2;

   for (i=0;i<(N-1);i++) {
     for (j=i+1;j<N;j++) {
       dx  = (rx[i]-rx[j]);
       dy  = (ry[i]-ry[j]);
       dz  = (rz[i]-rz[j]);
       if (dx>hL)       dx-=L;
       else if (dx<-hL) dx+=L;
       if (dy>hL)       dy-=L;
       else if (dy<-hL) dy+=L;
       if (dz>hL)       dz-=L;
       else if (dz<-hL) dz+=L;
       r2 = dx*dx + dy*dy + dz*dz;
       if (r2<rc2) {
	 bin=(int)(sqrt(r2)/dr);
	 H[bin]+=2;
       }
     }
   }
}

int main ( int argc, char * argv[] ) {

  double * rx, * ry, * rz;
  int N=216, nB=0;
  int * H;
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
    if (!strcmp(argv[i],"-dr")) dr=atof(argv[++i]);
    else if (!strcmp(argv[i],"-rc")) rc2=atof(argv[++i]);
    else if (!strcmp(argv[i],"-for"))
      sscanf(argv[++i],"%i,%i,%i",&start,&stop,&step);
    else if (!strcmp(argv[i],"-rho")) rho=atof(argv[++i]);
    else if (!strcmp(argv[i],"-fnf")) fnf=argv[++i];
    else if (!strcmp(argv[i],"-h")) {
      usage(); exit(0);
    }
    else {
      fprintf(stderr,"Error: Command-line argument '%s' not recognized.\n",
	      argv[i]);
      exit(-1);
    }
  }

  /* compute the number of bins, and allocate the histogram */
  nB = (int)(rc2/dr) + 1;
  H = (int*)calloc(nB,sizeof(int));

  /* Compute the *squared* cutoff, reusing the variable rc2 */
  rc2*=rc2;

  /* Allocate the position arrays */
  rx = (double*)malloc(N*sizeof(double));
  ry = (double*)malloc(N*sizeof(double));
  rz = (double*)malloc(N*sizeof(double));

  ngr=0;
  for (i=start;i<stop+1;i+=step) {
    /* Generate the file name */
    sprintf(fn,fnf,i);
    /* Open the file */
    fp=fopen(fn,"r");
    /* If the file is there... */
    if (fp) {
      /* Update the number of samples */
      ngr++;
      /* Read in the data, ignoring velocities */
      xyz_in(fp,rx,ry,rz,&N);
      /* Close the file */
      fclose(fp);
      /* If this is the first snapshot, take this opportunity
	 to compute the side-length from the given density
	 and detected number of particles from the first
	 read-in. */
      if (i==start) {
	/* Compute the side-length */
	L = pow((V=N/rho),0.3333333);
	/* if the cutoff is longer than half the box-length, spit
	   out a warning message, then reassign the cutoff appropriately. */
	if (sqrt(rc2)>L/2) {
	  fprintf(stderr,"# Warning: the maximum radius you have requested,"
		  " %.5lf, is greater than one-half the box length, %.5lf\n",
		  sqrt(rc2),L/2);
	  fprintf(stderr,"# Reducing rc:\n");
	  rc2 = L/2;
	  rc2 *= rc2;
	}
      }
      /* Update the histogram by sampling the configuration */
      update_hist(rx,ry,rz,N,L,rc2,dr,H);
      fprintf(stderr,"# %s...\n",fn);
    }
  }
  /* Normalize and output g(r) to the terminal */
  for (i=0;i<nB;i++) {
    r=dr*(i+0.5);
    vb=((i+1)*(i+1)*(i+1)-i*i*i)*dr*dr*dr;
    nid=(4./3.)*M_PI*vb*rho;
    fprintf(stdout,"%.4lf %.4lf\n",i*dr,(double)(H[i])/(ngr*N*nid));
  }
}
