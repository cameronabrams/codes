/* 
   Molecular Dynamics simulation of a Lennard-Jones fluid
   in a periodic boundary using the Nose-Hoover Thermostat

   Cameron F. Abrams

   Written for the course CHE 800-002, Molecular Simulation
   Spring 0304

   compile using "gcc -o mdlj_nhc mdlj_nhc.c -lm -lgsl"
   (assumes the GNU Scientific Library is installed)

   You must have the GNU Scientific Library installed; see
   the coursenotes to learn how to do this.

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2004
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Prints usage information */
void usage ( void ) {
  fprintf(stdout,"mdlj_and usage:\n");
  fprintf(stdout,"mdlj_and [options]\n\n");
  fprintf(stdout,"Options:\n");
  fprintf(stdout,"\t -N [integer]\t\tNumber of particles\n");
  fprintf(stdout,"\t -rho [real]\t\tNumber density\n");
  fprintf(stdout,"\t -dt [real]\t\tTime step\n");
  fprintf(stdout,"\t -rc [real]\t\tCutoff radius\n");
  fprintf(stdout,"\t -ns [real]\t\tNumber of integration steps\n");
  fprintf(stdout,"\t -so       \t\tShort-form output (unused)\n");
  fprintf(stdout,"\t -T0 [real]\t\tSetpoint temperature\n");
  fprintf(stdout,"\t -Q [real,real]\t\tThermostat masses (2)\n");  
  fprintf(stdout,"\t -Tjump [integer,real]\t\tTemperature jump at a given time-step to a given temperature\n");
  fprintf(stdout,"\t -fs [integer]\t\tSample frequency\n");
  fprintf(stdout,"\t -sf [a|w]\t\tAppend or write config output file\n");
  fprintf(stdout,"\t -icf [string]\t\tInitial configuration file\n");
  fprintf(stdout,"\t -seed [integer]\tRandom number generator seed\n");
  fprintf(stdout,"\t -uf          \t\tPrint unfolded coordinates in output files\n");
  fprintf(stdout,"\t -h           \t\tPrint this info\n");
}

/* Writes the coordinates in XYZ format to the output stream fp.
   The integer "z" is the atomic number of the particles, required
   for the XYZ format. The array ix contains the number of x-dir 
   periodic boundary crossings a particle has performed; thus,
   the "unfolded" coordinate is rx[i]+ix[i]*L. */
void xyz_out (FILE * fp, 
	      double * rx, double * ry, double * rz, 
	      double * vx, double * vy, double * vz, 
	      int * ix, int * iy, int * iz, double L,
	      int N, 
	      double * xi, double * vxi, double * Q, int nl,
	      int z, int put_vel, int unfold) {
  int i;

  fprintf(fp,"%i %i\n\n",N,put_vel);
  for (i=0;i<N;i++) {
    fprintf(fp,"%i %.8lf %.8lf %.8lf ",z,
	    rx[i]+(unfold?(ix[i]*L):0.0),
	    ry[i]+(unfold?(iy[i]*L):0.0),
	    rz[i]+(unfold?(iz[i]*L):0.0));
    if (put_vel)
      fprintf(fp,"%.8lf %.8lf %.8lf",vx[i],vy[i],vz[i]);
    fprintf(fp,"\n");
  }
}

int xyz_in (FILE * fp, double * rx, double * ry, double * rz, 
	     double * vx, double * vy, double * vz, 
	     int * N) {
  int i;
  int has_vel, dum;
  fscanf(fp,"%i %i\n\n",N,&has_vel);
  for (i=0;i<(*N);i++) {
    fscanf(fp,"%i %lf %lf %lf ",&dum,&rx[i],&ry[i],&rz[i]);
    if (has_vel) { // read velocities
      fscanf(fp,"%lf %lf %lf",&vx[i],&vy[i],&vz[i]);
    }
  }
  return has_vel;
}

/* An N^2 algorithm for computing forces and potential energy.  The virial
   is also computed and returned in *vir. */
double total_e ( double * rx, double * ry, double * rz, 
		 double * fx, double * fy, double * fz, 
		 int N, double L,
		 double rc2, double ecor, double ecut, double * vir ) {
   int i,j;
   double dx, dy, dz, r2, r6i;
   double e = 0.0, hL=L/2.0,f;

   /* Zero the forces */
   for (i=0;i<N;i++) {
     fx[i]=fy[i]=fz[i]=0.0;
   }
   
   *vir=0.0;
   for (i=0;i<(N-1);i++) {
     for (j=i+1;j<N;j++) {
	dx  = (rx[i]-rx[j]);
	dy  = (ry[i]-ry[j]);
	dz  = (rz[i]-rz[j]);
	/* Periodic boundary conditions: Apply the minimum image
	   convention; note that this is *not* used to truncate the
	   potential as long as there an explicit cutoff. */
	if (dx>hL)       dx-=L;
	else if (dx<-hL) dx+=L;
	if (dy>hL)       dy-=L;
	else if (dy<-hL) dy+=L;
	if (dz>hL)       dz-=L;
	else if (dz<-hL) dz+=L;
	r2 = dx*dx + dy*dy + dz*dz;
	if (r2<rc2) {
	  r6i   = 1.0/(r2*r2*r2);
	  e    += 4*(r6i*r6i - r6i) - ecut;
	  f     = 48*(r6i*r6i-0.5*r6i);
	  fx[i] += dx*f/r2;
	  fx[j] -= dx*f/r2;
	  fy[i] += dy*f/r2;
	  fy[j] -= dy*f/r2;
	  fz[i] += dz*f/r2;
	  fz[j] -= dz*f/r2;
	  *vir += f;
	}
     }
   }
   return e+N*ecor;
}

/* Initialize particle positions by assigning them
   on a cubic grid, then scaling positions 
   to achieve a given box size and thereby, volume,
   and density */
void init ( double * rx, double * ry, double * rz,
	    double * vx, double * vy, double * vz,
	    int * ix, int * iy, int * iz,
	    int n, 
	    double * xi, double * vxi,
	    int nl, double L, gsl_rng * r, double T0,
	    double * KE, char * icf) {
  int i,iix,iiy,iiz;
  double cmvx=0.0,cmvy=0.0,cmvz=0.0;
  double T, fac;
  int n3=2;
  int vel_ok=0;
  
  /* If icf has a value, assume it is the name of a file containing
     the input configuration in XYZ format */
  if (icf) {
    FILE * fp = fopen(icf,"r");
    if (fp) vel_ok = xyz_in(fp,rx,ry,rz,vx,vy,vz,&n);
    else {
      fprintf(stderr,"# error: could not read %s\n",icf);
      exit(-1);
    }
  }
  /* Assign particles on a cubic lattice */
  else {

    /* Find the lowest perfect cube, n3, greater than or equal to the
       number of particles */
    while ((n3*n3*n3)<n) n3++;
  
    iix=iiy=iiz=0;
    /* Assign particle positions */
    for (i=0;i<n;i++) {
      rx[i] = ((double)iix+0.5)*L/n3;
      ry[i] = ((double)iiy+0.5)*L/n3;
      rz[i] = ((double)iiz+0.5)*L/n3;
      iix++;
      if (iix==n3) {
	iix=0;
	iiy++;
	if (iiy==n3) {
	  iiy=0;
	  iiz++;
	}
      }
    }
  }
  /* If no velocities yet assigned, randomly pick some */
  if (!vel_ok) {
    for (i=0;i<n;i++) {
      vx[i]=gsl_ran_gaussian(r,1.0);
      vy[i]=gsl_ran_gaussian(r,1.0);
      vz[i]=gsl_ran_gaussian(r,1.0);
    }
  }
  /* Take away any center-of-mass drift; compute initial KE */
  for (i=0;i<n;i++) {
    cmvx+=vx[i];
    cmvy+=vy[i];
    cmvz+=vz[i];
  }
  (*KE)=0;
  for (i=0;i<n;i++) {
    vx[i]-=cmvx/n;
    vy[i]-=cmvy/n;
    vz[i]-=cmvz/n;
    (*KE)+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
  }
  (*KE)*=0.5;
  /* if T0 is specified, scale velocities */
  if (T0>0.0) {
    T=(*KE)/n*2./3.;
    fac=sqrt(T0/T);
    (*KE)=0;
    for (i=0;i<n;i++) {
      vx[i]*=fac;
      vy[i]*=fac;
      vz[i]*=fac;
      (*KE)+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
    }
    (*KE)*=0.5;
  }
  /* Initialize periodic boundary crossing counter arrays */
  memset(ix,0,n*sizeof(int));
  memset(iy,0,n*sizeof(int));
  memset(iz,0,n*sizeof(int));

  memset(xi,0,nl*sizeof(double));
  memset(vxi,0,nl*sizeof(double));
}

/* Chain half integration step (Algorithm 30) */
void chain ( double * KE, double dt, double dt_2,
	     double dt_4, double dt_8,
	     double * Q, double * xi, double * vxi,
	     double * vx, double * vy, double * vz,
	     int nl, int N, double T) {
  double G1, G2, s;
  int i;

  G2= (Q[0]*vxi[0]*vxi[0]-T);
  vxi[1]+=G2*dt_4;
  vxi[0]*=exp(-vxi[1]*dt_8);
  G1=(2*(*KE)-3*N*T)/Q[0];
  vxi[0]+=G1*dt_4;
  vxi[0]*=exp(-vxi[1]*dt_8);
  xi[0]+=vxi[0]*dt_2;
  xi[1]+=vxi[1]*dt_2;
  s=exp(-vxi[0]*dt_2);
  for (i=0;i<N;i++) {
    vx[i]*=s; vy[i]*=s; vz[i]*=s;
  }
  (*KE)*=(s*s);
  vxi[0]*=exp(-vxi[1]*dt_8);
  G1=(2*(*KE)-3*N*T)/Q[0];
  vxi[0]+=G1*dt_4;
  vxi[0]*=exp(-vxi[1]*dt_8);
  G2=(Q[0]*vxi[0]*vxi[0]-T)/Q[1];
  vxi[1]+=G2*dt_4;
}

int main ( int argc, char * argv[] ) {

  double * rx, * ry, * rz;
  double * vx, * vy, * vz;
  double * fx, * fy, * fz;
  int * ix, * iy, * iz;
  int nl;
  double * xi, * vxi, * Q;
  int N=216,c,a;
  double L=0.0;
  double rho=0.5, Tb = 1.0, nu=1.0, rc2 = 1.e20;
  double vir, vir_sum, pcor, V;
  double PE, KE, TE, ecor, ecut, T0=0.0, TE0;
  double rr3,dt=0.001, dt2, dt_2, dt_4, dt_8, sigma;
  int tj_ts=-1;
  double tj_Tb=1.0;
  int i,j,s;
  int nSteps = 10, fSamp=100;
  int short_out=0;
  int use_e_corr=0;
  int unfold = 0;

  char fn[20];
  FILE * out;
  char * wrt_code_str = "w";
  char * init_cfg_file = NULL;

  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned long int Seed = 23410981;

  /* Allocate arrays for the NH Chain */
  nl = 2;  /* for now... */
  xi = (double*)malloc(nl*sizeof(double));
  vxi = (double*)malloc(nl*sizeof(double));
  Q = (double*)malloc(nl*sizeof(double));
  /* Default masses */
  Q[0] = Q[1] = 0.1;

  /* Here we parse the command line arguments;  If
   you add an option, document it in the usage() function! */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-rho")) rho=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nu")) nu=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dt")) dt=atof(argv[++i]);
    else if (!strcmp(argv[i],"-rc")) rc2=atof(argv[++i]);
    else if (!strcmp(argv[i],"-ns")) nSteps = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-so")) short_out=1;
    else if (!strcmp(argv[i],"-T0")) T0=atof(argv[++i]);
    else if (!strcmp(argv[i],"-Tb")) Tb=atof(argv[++i]);
    else if (!strcmp(argv[i],"-Q")) sscanf(argv[++i],"%lf,%lf",
					   &Q[0],&Q[1]);
    else if (!strcmp(argv[i],"-Tjump")) sscanf(argv[++i],"%i,%lf",
					       &tj_ts,&tj_Tb);
    else if (!strcmp(argv[i],"-fs")) fSamp=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-sf")) wrt_code_str = argv[++i];
    else if (!strcmp(argv[i],"-icf")) init_cfg_file = argv[++i];
    else if (!strcmp(argv[i],"-ecorr")) use_e_corr = 1;
    else if (!strcmp(argv[i],"-seed")) Seed = (unsigned long)atoi(argv[++i]);
    else if (!strcmp(argv[i],"-uf")) unfold = 1;
    else if (!strcmp(argv[i],"-h")) {
      usage(); exit(0);
    }
    else {
      fprintf(stderr,"Error: Command-line argument '%s' not recognized.\n",
	      argv[i]);
      exit(-1);
    }
  }

  /* Compute the side-length */
  L = pow((V=N/rho),0.3333333);

  /* Compute the tail-corrections; assumes sigma and epsilon are both 1 */
  rr3 = 1.0/(rc2*rc2*rc2);
  ecor = use_e_corr?8*M_PI*rho*(rr3*rr3*rr3/9.0-rr3/3.0):0.0;
  pcor = use_e_corr?16.0/3.0*M_PI*rho*rho*(2./3.*rr3*rr3*rr3-rr3):0.0;
  ecut = 4*(rr3*rr3*rr3*rr3-rr3*rr3);

  /* Compute the *squared* cutoff, reusing the variable rc2 */
  rc2*=rc2;

  /* compute the squared time step */
  dt2  = dt*dt;
  dt_2 = 0.5*dt;
  dt_4 = 0.5*dt_2;
  dt_8 = 0.5*dt_4;  // thanks, jgrime@uchicago.edu

  /* Compute sigma */
  sigma = sqrt(Tb);

  /* Output some initial information */
  fprintf(stdout,"# Nose-Hoover-Chain-Thermostat MD Simulation"
	  " of a Lennard-Jones fluid\n");
  fprintf(stdout,"# L = %.5lf; rho = %.5lf; N = %i; rc = %.5lf\n",
	  L,rho,N,sqrt(rc2));
  fprintf(stdout,"# nSteps %i, seed %d, dt %.5lf, "
	  "Tb %.5lf, Q0 %.5lf Q1 %.5lf\n",
	  nSteps,Seed,dt,Tb,Q[0],Q[1]);
  
  /* Seed the random number generator */
  gsl_rng_set(r,Seed);
  
  /* Allocate the position arrays */
  rx = (double*)malloc(N*sizeof(double));
  ry = (double*)malloc(N*sizeof(double));
  rz = (double*)malloc(N*sizeof(double));

  /* Allocate the boundary crossing counter arrays */
  ix = (int*)malloc(N*sizeof(int));
  iy = (int*)malloc(N*sizeof(int));
  iz = (int*)malloc(N*sizeof(int));

  /* Allocate the velocity arrays */
  vx = (double*)malloc(N*sizeof(double));
  vy = (double*)malloc(N*sizeof(double));
  vz = (double*)malloc(N*sizeof(double));

  /* Allocate the force arrays */
  fx = (double*)malloc(N*sizeof(double));
  fy = (double*)malloc(N*sizeof(double));
  fz = (double*)malloc(N*sizeof(double));

  /* Generate initial positions on a cubic grid, 
     and measure initial energy */
  init(rx,ry,rz,vx,vy,vz,ix,iy,iz,N,xi,vxi,nl,L,r,T0,&KE,init_cfg_file);
  sprintf(fn,"%i.xyz",0);
  out=fopen(fn,"w");
  xyz_out(out,rx,ry,rz,
	      vx,vy,vz,ix,iy,iz,
	      L,N,
	      xi,vxi,Q,nl,
	      16,1,unfold);
  fclose(out);

  PE = total_e(rx,ry,rz,fx,fy,fz,N,L,rc2,ecor,ecut,&vir);
  TE0=PE+KE;
  
  fprintf(stdout,"# step PE KE TE drift T P\n");

  /* Nose-Hoover-Chain (Algorithms 30, 31, 32) */
  for (s=0;s<nSteps;s++) {

    /* do a temperature jump at the prescribed time */
    if (s==tj_ts) Tb = tj_Tb;

    chain(&KE,dt,dt_2,dt_4,dt_8,Q,xi,vxi,vx,vy,vz,nl,N,Tb);
    
    /* First integration half-step */
    KE = 0.0;
    for (i=0;i<N;i++) {
      rx[i]+=vx[i]*dt_2;
      ry[i]+=vy[i]*dt_2;
      rz[i]+=vz[i]*dt_2;
      /* Apply periodic boundary conditions */
      if (rx[i]<0.0) { rx[i]+=L; ix[i]--; }
      if (rx[i]>L)   { rx[i]-=L; ix[i]++; }
      if (ry[i]<0.0) { ry[i]+=L; iy[i]--; }
      if (ry[i]>L)   { ry[i]-=L; iy[i]++; }
      if (rz[i]<0.0) { rz[i]+=L; iz[i]--; }
      if (rz[i]>L)   { rz[i]-=L; iz[i]++; }
    }
    /* Calculate forces */
    PE = total_e(rx,ry,rz,fx,fy,fz,N,L,rc2,ecor,ecut,&vir);
      
    /* Second integration half-step */
    for (i=0;i<N;i++) {
      vx[i]+=dt*fx[i];
      vy[i]+=dt*fy[i];
      vz[i]+=dt*fz[i];
      rx[i]+=vx[i]*dt_2;
      ry[i]+=vy[i]*dt_2;
      rz[i]+=vz[i]*dt_2;
      KE+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
    }
    KE*=0.5;

    chain(&KE,dt,dt_2,dt_4,dt_8,Q,xi,vxi,vx,vy,vz,nl,N,Tb);

    TE=PE+KE;
    fprintf(stdout,"%i %.5lf %.5lf %.5lf %.5lf %.5le %.5lf %.5lf\n",
	    s,s*dt,PE,KE,TE,(TE-TE0)/TE0,KE*2/3./N,rho*KE*2./3./N+vir/3.0/V);
    if (!(s%fSamp)) {
      sprintf(fn,"%i.xyz",!strcmp(wrt_code_str,"a")?0:s);
      out=fopen(fn,wrt_code_str);
      xyz_out(out,rx,ry,rz,vx,vy,vz,ix,iy,iz,L,N,xi,vxi,Q,nl,16,1,unfold);
      fclose(out);
    }
  }
}
