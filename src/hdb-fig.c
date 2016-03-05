/* 
   Metropolis Monte Carlo simulation of 2D hard dumbbells confined in a circle

   Cameron F. Abrams

   Written for the course CHE 800-002, Molecular Simulation
   Spring 0304

   compile using "gcc -o hdb hdb.c -lm -lgsl"
   (assumes the GNU Scientific Library is installed)

   runs as "./hdisk -N <number_of_particles> -rho <density> \
                    -R <radius_of_circle> \
                    -nc <numcycles(1e6)>  \
		    -seed <seed(?)> \
		    -dw <weight_of_displacement_trial_move(0.5)> \
		    -dr <max_displacement(1.0)> \
		    -da <max_angle_rotation(pi/2)> \
		    -r_0 <dumbbell bond-length(1.0)> \
		    -s <particle_diameter(1.0)>"

   You must have the GNU Scientific Library installed; see
   the coursenotes to learn how to do this.

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2004
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

void out_fig ( FILE * fp, double * rx, double * ry, int n, double R2 ) {
  int i;
  double R = sqrt(R2);
  int iR = (int)(R*300);
  fprintf(fp,"#FIG 3.2\n" "Landscape\n" "Center\n" "Inches\n" "Letter\n"
	  "100.00\n" "Single\n" "-2\n" "1200 2\n");
  fprintf(fp,"1 3 0 1 0 7 50 -1 -1 0.000 1 0.0000 %i %i %i %i %i %i %i %i\n",iR,iR,iR,iR,iR,iR,2*iR,iR);
  for (i=0;i<n;i++) {
    int x = (int)((rx[i] + R)*300);
    int y = (int)((ry[i] + R)*300);
    fprintf(fp,"1 3 0 1 0 7 50 -1 -1 0.000 1 0.0000 %i %i "
	    "150 150 %i %i %i %i\n",x,y,x,y,x+150,y);
    if (!(i%2)) {
      int x2 = (int)((rx[i+1] + R)*300);
      int y2 = (int)((ry[i+1] + R)*300);
      fprintf(fp,"2 1 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 2\n"
			"%i %i %i %i\n",x,y,x2,y2);
    }
  }
}

/* Initialize dumbbell positions by assigning them
   randomly while avoiding overlaps. */
void init ( double * rx, double * ry, 
	    int n, double R2, double s2, double r_02, gsl_rng * r,
	    int fig_out) {
  int i,j,reject;
  double R = sqrt(R2), r2, sx, sy, r_0 = sqrt(r_02), a;
  /* nMax is a maximum number of insertion trials */
  int nMax = 50000, nTrials=0, nPartTrials=0;
  FILE * tmp_fp;

  for (i=0;i<n;i++) {
    reject = 1;
    nTrials = 0;
    while (reject&&nTrials<nMax) {
      reject = 0;
      rx[i] = R*(1.0-2*gsl_rng_uniform(r));
      ry[i] = sqrt(R2-rx[i]*rx[i])*(1.0-2*gsl_rng_uniform(r));
      for (j=0;!reject&&j<i;j++) {
	sx  = rx[i]-rx[j];
	sy  = ry[i]-ry[j];
	r2  = sx*sx + sy*sy;
	if (r2 < s2) reject=1;
      }
      if (!reject) {
      /* position of i is accepted.  try to place its partner */
	i++;
	reject=1;
	nPartTrials=0;
	while (reject&&nPartTrials<nMax) {
	  reject=0;
	  a=gsl_rng_uniform(r)*M_PI*2;
	  rx[i] = r_0*cos(a);
	  ry[i] = r_0*sin(a);
	  rx[i] += rx[i-1];
	  ry[i] += ry[i-1];
	  if ((rx[i]*rx[i]+ry[i]*ry[i])>R2) reject=1;
	  for (j=0;(!reject)&&(j<i);j++) {
	    if (j!=i) {
	      sx  = rx[i]-rx[j];
	      sy  = ry[i]-ry[j];
	      r2  = sx*sx + sy*sy;
	      if (r2 < s2) reject=1;
	    }
	  }
	  nPartTrials++;
	}
	if (nPartTrials==nMax) {
	  reject=1;
	}
      }
      nTrials++;
    }
    if (nTrials==nMax) {
      fprintf(stderr,"# Error -- could not initialize position "
	      "of particle %i in %i trials\n",i,nTrials);
      fprintf(stderr,"# program ends.\n");
      exit(-1);
    }
  }
  if (fig_out) {
    tmp_fp = fopen("0.fig","w");
    out_fig(tmp_fp,rx,ry,n,R2);
    fclose(tmp_fp);
    fprintf(stdout,"# init complete.\n");
  }
}
enum { DISP, ROTATE };

int main ( int argc, char * argv[] ) {

  double * rx, * ry;
  int N=-1,c,a;
  double R2=16.0, s2=1.0, r_02=1.0, dA, cs, sn, xc, yc;
  double rho=0.5;
  double xs[2], ys[2], disp_wt = 0.5;
  double dr=0.1,dx,dy,theta,sx,sy,r2, da=M_PI/180.0, nx,ny;
  int i,j,k;
  int nCycles = 10;
  int reject;
  int nDispAttempt, nRotAttempt, nDispAcc, nRotAcc;
  int noob=0,novl=0, move;
  int fSamp=1000;
  int fig_out = 0;
  char fn[20];

  FILE * tmp_fp;

  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned long int Seed = 23410981;

  /* Here we parse the command line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-N"))        N  = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-R"))   R2 = atof(argv[++i]);
    else if (!strcmp(argv[i],"-rho")) rho= atof(argv[++i]);
    else if (!strcmp(argv[i],"-dr"))  dr = atof(argv[++i]);
    else if (!strcmp(argv[i],"-da"))  da = atof(argv[++i]);
    else if (!strcmp(argv[i],"-s"))   s2 = atof(argv[++i]);
    else if (!strcmp(argv[i],"-r_0")) r_02=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dw"))  disp_wt=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nc"))  nCycles = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-fs"))  fSamp = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-fo"))  fig_out = 1;
    else if (!strcmp(argv[i],"-seed"))Seed = (unsigned long)atoi(argv[++i]);
  }

  /* Some checks */
  if (disp_wt < 0.0 || disp_wt > 1.0) {
    fprintf(stderr,"# Error: please select a weight between 0 and 1\n");
    exit(-1);
  }

  s2=s2*s2;
  r_02=r_02*r_02;

  /* If N was not set by the user, calculate it based on 
     the given values of density and radius */
  if (N==-1) {
    N = (int)(rho*M_PI*R2);
    if (N%2) N--;
  }
  /* Otherwise, recompute the radius. */
  else {
    if (N%2) N--;
    R2 = ((double)N)/rho/M_PI;
  }

  fprintf(stdout,"# R = %.5lf; rho = %.5lf; N = %i; r_0 = %.5lf; s = %.5lf\n",
	  sqrt(R2),rho,N,sqrt(r_02),sqrt(s2));

  /* Seed the random number generator */
  gsl_rng_set(r,Seed);

  /* allocate the position arrays */
  rx = (double*)malloc(N*sizeof(double));
  ry = (double*)malloc(N*sizeof(double));

  /* generate initial positions that fit inside 
     circle with radius R and guarantee no 
     particles overlap. */
  init(rx,ry,N,R2,s2,r_02,r,fig_out);

  nDispAcc = nRotAcc = 0;
  nDispAttempt = nRotAttempt = 0;
  for (c=0;c<nCycles;c++) {
    /* Make N attempts */
    for (a=0;a<N;a++) {
      /* randomly select a molecule */
      i=(int)gsl_rng_uniform_int(r,N/2);
      i*=2;  /* N is assumed even */
      /* save its position */
      for (k=0;k<2;k++) {
	xs[k] = rx[i+k];
	ys[k] = ry[i+k];
      }
      /* decide whether to move this molecule, or to rotate it. */
      if (gsl_rng_uniform(r) < disp_wt) {
	move = DISP;
	/* calc displacement */
	dx = dr*(0.5-gsl_rng_uniform(r));
	dy = dr*(0.5-gsl_rng_uniform(r));
	/* displace molecule */
	for (k=0;k<2;k++) {
	  rx[i+k]+=dx;
	  ry[i+k]+=dy;
	}
	nDispAttempt++;
      }
      else {
	move = ROTATE;
	/* compute random angle */
	dA = da*(0.5-gsl_rng_uniform(r));
	cs=cos(dA);
	sn=sin(dA);
	/* translate COM to origin (convert to internal coords) */
	xc=0.5*(rx[i]+rx[i+1]);
	yc=0.5*(ry[i]+ry[i+1]);
	for (k=0;k<2;k++) {
	  rx[i+k]-=xc;
	  ry[i+k]-=yc;
	}
	/* rotate and transform back to global coords*/
	for (k=0;k<2;k++) {
	  nx  = rx[i+k]*cs - ry[i+k]*sn;
	  ny  = rx[i+k]*sn + ry[i+k]*cs;
	  rx[i+k] = nx+xc;
	  ry[i+k] = ny+yc;
	}
	nRotAttempt++;
      }

      /* Determine whether we should reject this move */
      /* First, check to make sure molecule stays inside circle */
      reject=0;
      for (k=0;!reject&&k<2;k++) {
	/* compute new distance to origin */
	r2=rx[i+k]*rx[i+k]+ry[i+k]*ry[i+k];
	/* reject move if outside circle */
	reject=(r2>R2);
	if (reject) noob++;
      }

      /* Second, check for overlaps with other particles */
      for (j=0;(!reject)&&(j<N);j++) {
	if (j!=i&&j!=(i+1)) {
	  sx  = rx[i]-rx[j];
	  sy  = ry[i]-ry[j];
	  r2  = sx*sx + sy*sy;
	  reject = r2<s2;
	  /* check the partner */
	  if (!reject) {
	    sx  = rx[i+1]-rx[j];
	    sy  = ry[i+1]-ry[j];
	    r2  = sx*sx + sy*sy;
	    reject = r2<s2;
	  }
	}
	if (reject) novl++;
      }
      

      /* if move is rejected, undo trial move */
      if (reject) {
	for (k=0;k<2;k++) {
	  rx[i+k]=xs[k];
	  ry[i+k]=ys[k];
	}
      }
      else {
	if (move == DISP) nDispAcc++;
	else if (move == ROTATE) nRotAcc++;
      }
    }
    if (!(c%fSamp)) {
      if (fig_out) {
	sprintf(fn,"%i.fig",c);
	tmp_fp = fopen(fn,"w");
	out_fig(tmp_fp,rx,ry,N,R2);
	fclose(tmp_fp);
      }
    }
  }
  fprintf(stdout,"Results:\n"
	  "Number of trial moves:             %i\n"
	  "Maximum displacement length:       %.3lf\n"
	  "Number of displacement attempts:   %i\n"
	  "Maximum rotation angle (radians):  %.3lf\n"
	  "Number of rotation attempts:       %i\n"
	  "Displacement acceptance ratio:     %.3lf\n"
	  "Rotation acceptance ratio:         %.3lf\n"
	  "Reject Fraction Out-of-bounds:     %.5lf\n"
	  "Reject Fraction Overlap:           %.5lf\n",
	  N*nCycles,dr,nDispAttempt,da,nRotAttempt,
	  nDispAttempt?((double)nDispAcc)/nDispAttempt:0,
	  nRotAttempt?((double)nRotAcc)/nRotAttempt:0,
	  ((double)noob)/(N*nCycles-nRotAcc-nDispAcc),
	  ((double)novl)/(N*nCycles-nRotAcc-nDispAcc)
	  );
}
