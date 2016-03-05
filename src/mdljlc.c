/* 
   Microcanonical Molecular Dynamics simulation of a Lennard-Jones fluid
   in a periodic boundary

   Cameron F. Abrams

   Written for the course CHE 800-002, Molecular Simulation
   Spring 0304
   (with link-cells)

   compile using "gcc -o mdlj mdlj.c -lm -lgsl"
   (assumes the GNU Scientific Library is installed)

   You must have the GNU Scientific Library installed; see
   the coursenotes to learn how to do this.

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2006
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Prints usage information */
void usage ( void ) {
  fprintf(stdout,"mdlj usage:\n");
  fprintf(stdout,"mdlj [options]\n\n");
  fprintf(stdout,"Options:\n");
  fprintf(stdout,"\t -N [integer]\t\tNumber of particles\n");
  fprintf(stdout,"\t -rho [real]\t\tNumber density\n");
  fprintf(stdout,"\t -dt [real]\t\tTime step\n");
  fprintf(stdout,"\t -rc [real]\t\tCutoff radius\n");
  fprintf(stdout,"\t -ns [real]\t\tNumber of integration steps\n");
  fprintf(stdout,"\t -so       \t\tShort-form output (unused)\n");
  fprintf(stdout,"\t -T0 [real]\t\tInitial temperature\n");
  fprintf(stdout,"\t -fs [integer]\t\tSample frequency\n");
  fprintf(stdout,"\t -sf [a|w]\t\tAppend or write config output file\n");
  fprintf(stdout,"\t -icf [string]\t\tInitial configuration file\n");
  fprintf(stdout,"\t -seed [integer]\tRandom number generator seed\n");
  fprintf(stdout,"\t -uf          \t\tPrint unfolded coordinates in output files\n");
  fprintf(stdout,"\t -h           \t\tPrint this info\n");
}


/* Returns the cell index of the 'index'th neighbor of cell ic */
uint nborcell ( uint ic, uint index, uint ncx, uint ncy, uint ncz ) {
  uint icz=ic/(ncx*ncy);
  uint icy=(ic-icz*ncx*ncy)/ncx;
  uint icx=ic-icz*ncx*ncy-icy*ncx;
  int jcx,jcy,jcz;

  if (index==0) {
    jcx=icx-1;  if (jcx<0) jcx=ncx-1;
    jcy=icy-1;  if (jcy<0) jcy=ncy-1;
    jcz=icz+1;  if (jcz==ncz) jcz=0;
  }
  else  if (index==1) {
    jcx=icx-1;  if (jcx<0) jcx=ncx-1;
    jcy=icy;
    jcz=icz+1;  if (jcz==ncz) jcz=0;
  }
  else  if (index==2) {
    jcx=icx-1;  if (jcx<0) jcx=ncx-1;
    jcy=icy+1;  if (jcy==ncy) jcy=0;
    jcz=icz+1;  if (jcz==ncz) jcz=0;
  }
  else  if (index==3) {
    jcx=icx;  
    jcy=icy-1;  if (jcy<0) jcy=ncy-1;
    jcz=icz+1;  if (jcz==ncz) jcz=0;
  }
  else  if (index==4) {
    jcx=icx;  
    jcy=icy;
    jcz=icz+1;  if (jcz==ncz) jcz=0;
  }
  else  if (index==5) {
    jcx=icx;
    jcy=icy+1;  if (jcy==ncy) jcy=0;
    jcz=icz+1;  if (jcz==ncz) jcz=0;
  }
  else  if (index==6) {
    jcx=icx+1;  if (jcx==ncx) jcx=0;
    jcy=icy-1;  if (jcy<0) jcy=ncy-1;
    jcz=icz+1;  if (jcz==ncz) jcz=0;
  }
  else  if (index==7) {
    jcx=icx+1;  if (jcx==ncx) jcx=0;
    jcy=icy;
    jcz=icz+1;  if (jcz==ncz) jcz=0;
  }
  else  if (index==8) {
    jcx=icx+1;  if (jcx==ncx) jcx=0;
    jcy=icy+1;  if (jcy==ncy) jcy=0;
    jcz=icz+1;  if (jcz==ncz) jcz=0;
  }
  else  if (index==9) {
    jcx=icx-1; if (jcx<0) jcx=ncx-1;
    jcy=icy-1; if (jcy<0) jcy=ncy-1;
    jcz=icz;  if (jcz==ncz) jcz=0;
  }
  else  if (index==10) {
    jcx=icx-1; if (jcx<0) jcx=ncx-1;
    jcy=icy;
    jcz=icz;
  }
  else  if (index==11) {
    jcx=icx-1; if (jcx<0) jcx=ncx-1;
    jcy=icy+1; if (jcy==ncy) jcy=0;
    jcz=icz;
  }
  else  if (index==12) {
    jcx=icx;
    jcy=icy-1; if (jcy<0) jcy=ncy-1;
    jcz=icz;
  }
  else  if (index==13) {
    jcx=icx;
    jcy=icy;
    jcz=icz;
  }
  else  if (index==14) {
    jcx=icx;
    jcy=icy+1; if (jcy==ncy) jcy=0;
    jcz=icz;
  }
  else  if (index==15) {
    jcx=icx+1; if (jcx==ncx) jcx=0;
    jcy=icy-1; if (jcy<0) jcy=ncy-1;
    jcz=icz;
  }
  else  if (index==16) {
    jcx=icx+1; if (jcx==ncx) jcx=0;
    jcy=icy;
    jcz=icz;
  }
  else  if (index==17) {
    jcx=icx+1; if (jcx==ncx) jcx=0;
    jcy=icy+1; if (jcy==ncy) jcy=0;
    jcz=icz;
  }
  else  if (index==18) {
    jcx=icx-1; if (jcx<0) jcx=ncx-1;
    jcy=icy-1; if (jcy<0) jcy=ncy-1;
    jcz=icz-1; if (jcz<0) jcz=ncz-1;
  }
  else  if (index==19) {
    jcx=icx-1; if (jcx<0) jcx=ncx-1;
    jcy=icy;
    jcz=icz-1; if (jcz<0) jcz=ncz-1;
  }
  else  if (index==20) {
    jcx=icx-1; if (jcx<0) jcx=ncx-1;
    jcy=icy+1; if (jcy==ncy) jcy=0;
    jcz=icz-1; if (jcz<0) jcz=ncz-1;
  }
  else  if (index==21) {
    jcx=icx;
    jcy=icy-1; if (jcy<0) jcy=ncy-1;
    jcz=icz-1; if (jcz<0) jcz=ncz-1;
  }
  else  if (index==22) {
    jcx=icx;
    jcy=icy;
    jcz=icz-1; if (jcz<0) jcz=ncz-1;
  }
  else  if (index==23) {
    jcx=icx;
    jcy=icy+1; if (jcy==ncy) jcy=0;
    jcz=icz-1; if (jcz<0) jcz=ncz-1;
  }
  else  if (index==24) {
    jcx=icx+1; if (jcx==ncx) jcx=0;
    jcy=icy-1; if (jcy<0) jcy=ncy-1;
    jcz=icz-1; if (jcz<0) jcz=ncz-1;
  }
  else  if (index==25) {
    jcx=icx+1; if (jcx==ncx) jcx=0;
    jcy=icy;
    jcz=icz-1; if (jcz<0) jcz=ncz-1;
  }
  else  if (index==26) {
    jcx=icx+1; if (jcx==ncx) jcx=0;
    jcy=icy+1; if (jcy==ncy) jcy=ncy-1;
    jcz=icz-1; if (jcz<0) jcz=ncz-1;
  } 

  //fprintf(stderr,"# nborcell ic %d : %d %d %d | index %d nbor %d : %d %d %d\n",
  //	  ic,icx,icy,icz,index,jcz*ncx*ncy+jcy*ncx+jcx,jcx,jcy,jcz);


  return (uint)(jcz*ncx*ncy+jcy*ncx+jcx);
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
	      int N, int z, int put_vel, int unfold) {
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

void setup_linklist ( int N, double L, double rc2, 
		      int * ncx, int * ncy, int * ncz, int * nc,
		      double * rcx, double * rcy, double * rcz,
		      int ** ll, int ** head ) {

  (*ll) = (int*)malloc(N*sizeof(int));
  *ncx = (int)(L/sqrt(rc2));
  *ncy = (int)(L/sqrt(rc2));
  *ncz = (int)(L/sqrt(rc2)); 
  *nc=(*ncx)*(*ncy)*(*ncz);
  *rcx=L/(*ncx);
  *rcy=L/(*ncy);
  *rcz=L/(*ncz);
  (*head)=(int*)malloc((*nc)*sizeof(int));
  fprintf(stderr,"# setup_linklist: rc %.5lf rcx %.5lf rcy %.5lf rcz %.5lf\n",
	  sqrt(rc2),*rcx,*rcy,*rcz);
}

double total_e ( double * rx, double * ry, double * rz, 
		 double * fx, double * fy, double * fz, 
		 int N, double L,
		 double rc2, double ecor, double ecut, double * vir,
		 int * ll, int * head, 
		 int nc, int ncx, int ncy, int ncz,
		 double rcx, double rcy, double rcz, int use_n2 ) {
   int i;
   double dx, dy, dz, r2, r6i;
   double e = 0.0, hL=L/2.0,f;
   uint ix,iy,iz,ic,jc,jjc;
   int ii,jj;

   /* Zero the forces */
   for (i=0;i<N;i++) {
     fx[i]=fy[i]=fz[i]=0.0;
   }

   /* enlist the atoms */
   //fprintf(stderr,"# enlist %i atoms: nc[xzy]: %i %i %i nc %i rc[xyz]: %.3lf %.3lf %.3lf\n",N,
   //   	   ncx,ncy,ncz,nc,rcx,rcy,rcz);
   for (i=0;i<N;i++) ll[i] = -1;
   for (i=0;i<nc;i++) head[i]=-1;
   for (i=0;i<N;i++) {
      ix=(uint)(rx[i]/rcx);
      iy=(uint)(ry[i]/rcy);
      iz=(uint)(rz[i]/rcz);
      ic = iz*ncx*ncy+iy*ncx+ix;
      ll[i]=head[ic];
      head[ic]=i;
      //fprintf(stderr,"# enlist %i %.5lf %.5lf %.5lf : %i %d %d %d\n",
      //     i,rx[i],ry[i],rz[i],ic,ix,iy,iz);
    }

   if (use_n2) {
     for (ii=0;ii<N;ii++) {
       for (jj=ii+1;jj<N;jj++) {
	 dx  = (rx[ii]-rx[jj]);
	 dy  = (ry[ii]-ry[jj]);
	 dz  = (rz[ii]-rz[jj]);
	 if (dx>hL)       dx-=L;
	 else if (dx<-hL) dx+=L;
	 if (dy>hL)       dy-=L;
	 else if (dy<-hL) dy+=L;
	 if (dz>hL)       dz-=L;
	 else if (dz<-hL) dz+=L;
	 r2 = dx*dx + dy*dy + dz*dz;
	 if (r2<rc2) {
	   //fprintf(stderr,"#ccpp %d %d %d %d\n",ic,jjc,ii>jj?jj:ii,ii>jj?ii:jj);
           //fprintf(stderr,"# ii %i jj %i r %.5lf\n",ii,jj,sqrt(r2));
	   r6i   = 1.0/(r2*r2*r2);
	   e    += 4*(r6i*r6i - r6i) - ecut;
	   f     = 48*(r6i*r6i-0.5*r6i);
	   fx[ii] += dx*f/r2;
	   fx[jj] -= dx*f/r2;
	   fy[ii] += dy*f/r2;
	   fy[jj] -= dy*f/r2;
	   fz[ii] += dz*f/r2;
	   fz[jj] -= dz*f/r2;
	   *vir += f;
	 }
       }
     }
     return e+N*ecor;
   }

   *vir=0.0;
    for (ic=0;ic<nc;ic++) {
      for (ii=head[ic];ii>-1;ii=ll[ii]) {
	for (jc=0;jc<14;jc++) {
	  jjc=nborcell(ic,jc,ncx,ncy,ncz);
	  for (jj=head[jjc];jj>-1;jj=ll[jj]) {
	    if ((ic!=jjc)||(ic==jjc&&ii<jj)) {
	      dx  = (rx[ii]-rx[jj]);
	      dy  = (ry[ii]-ry[jj]);
	      dz  = (rz[ii]-rz[jj]);
	      if (dx>hL)       dx-=L;
	      else if (dx<-hL) dx+=L;
	      if (dy>hL)       dy-=L;
	      else if (dy<-hL) dy+=L;
	      if (dz>hL)       dz-=L;
	      else if (dz<-hL) dz+=L;
	      r2 = dx*dx + dy*dy + dz*dz;
	      if (r2<rc2) {
                double thise;
		//fprintf(stderr,"#ccpp %d %d %d %d\n",ic,jjc,ii>jj?jj:ii,ii>jj?ii:jj);
                r6i   = 1.0/(r2*r2*r2);
		e    += thise=4*(r6i*r6i - r6i) - ecut;
		f     = 48*(r6i*r6i-0.5*r6i);
		//fprintf(stderr,"# ii %i jj %i r %.8lf e %.8lf f %.8lf\n",ii,jj,sqrt(r2),thise,f/r2);
                fx[ii] += dx*f/r2;
		fx[jj] -= dx*f/r2;
		fy[ii] += dy*f/r2;
		fy[jj] -= dy*f/r2;
		fz[ii] += dz*f/r2;
		fz[jj] -= dz*f/r2;
		*vir += f;
	      }
	    }
	  }
	}
      }
    }
//fprintf(stderr,"# done\n");
    return e;//+N*ecor;
}

/* Initialize particle positions by assigning them
   on a cubic grid, then scaling positions 
   to achieve a given box size and thereby, volume,
   and density */
void init ( double * rx, double * ry, double * rz,
	    double * vx, double * vy, double * vz,
	    int * ix, int * iy, int * iz,
	    int n, double L, gsl_rng * r, double T0,
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
      vx[i]=gsl_ran_exponential(r,1.0);
      vy[i]=gsl_ran_exponential(r,1.0);
      vz[i]=gsl_ran_exponential(r,1.0);
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
}

int main ( int argc, char * argv[] ) {

  double * rx, * ry, * rz;
  double * vx, * vy, * vz;
  double * fx, * fy, * fz;
  int * ix, * iy, * iz;
  int N=216,c,a;
  double L=0.0;
  double rho=0.5, T=0.0, rc2 = 1.e20, vir, vir_old, vir_sum, pcor, V;
  double PE, KE, TE, ecor, ecut, T0=1.0, TE0;
  double rr3,dt=0.001, dt2;
  int i,j,s;
  int nSteps = 10, fSamp=100, fOut=1;
  int short_out=0;
  int use_e_corr=0;
  int unfold = 0;

  /* timing variables */
  clock_t t0, t1;
  double sec;

  /* link-list variables */
  int ncx, ncy, ncz, nc;
  int * ll, * head;
  double rcx, rcy, rcz;
  int use_n2=0;

  char fn[20];
  FILE * out;
  char * wrt_code_str = "w";
  char * init_cfg_file = NULL;
  
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned long int Seed = 23410981;

  /* Here we parse the command line arguments;  If
   you add an option, document it in the usage() function! */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-rho")) rho=atof(argv[++i]);
    else if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dt")) dt=atof(argv[++i]);
    else if (!strcmp(argv[i],"-rc")) rc2=atof(argv[++i]);
    else if (!strcmp(argv[i],"-ns")) nSteps = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-so")) short_out=1;
    else if (!strcmp(argv[i],"-T0")) T0=atof(argv[++i]);
    else if (!strcmp(argv[i],"-fs")) fSamp=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-fo")) fOut=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-sf")) wrt_code_str = argv[++i];
    else if (!strcmp(argv[i],"-icf")) init_cfg_file = argv[++i];
    else if (!strcmp(argv[i],"-ecorr")) use_e_corr = 1;
    else if (!strcmp(argv[i],"-seed")) Seed = (unsigned long)atoi(argv[++i]);
    else if (!strcmp(argv[i],"-uf")) unfold = 1;
    else if (!strcmp(argv[i],"-n2")) use_n2=1;
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

  /* Setup the link-list */
  setup_linklist(N,L,rc2,&ncx,&ncy,&ncz,&nc,&rcx,&rcy,&rcz,&ll,&head);

  /* compute the squared time step */
  dt2=dt*dt;

  /* Output some initial information */
  fprintf(stdout,"# NVE MD Simulation of a Lennard-Jones fluid\n");
  if (use_n2) fprintf(stdout,"# Using N2 loop\n");
  else fprintf(stdout,"# Using the link-cell algorithm; nc = %d\n", nc);
  fprintf(stdout,"# L = %.5lf; rho = %.5lf; N = %i; rc = %.5lf\n",
	  L,rho,N,sqrt(rc2));
  fprintf(stdout,"# nSteps %i, seed %d, dt %.5lf\n",
	  nSteps,Seed,dt);
  fflush(stdout);
  
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
  init(rx,ry,rz,vx,vy,vz,ix,iy,iz,N,L,r,T0,&KE,init_cfg_file);
  sprintf(fn,"%i.xyz",0);
  out=fopen(fn,"w");
  xyz_out(out,rx,ry,rz,vx,vy,vz,ix,iy,iz,L,N,16,1,unfold);
  fclose(out);

  t0=clock();

  PE = total_e(rx,ry,rz,fx,fy,fz,N,L,rc2,ecor,ecut,&vir_old,
	       ll,head,nc,ncx,ncy,ncz,rcx,rcy,rcz,use_n2);
  TE0=PE+KE;
  
  fprintf(stdout,"# step PE KE TE drift T P\n");

  for (s=0;s<nSteps;s++) {

    /* First integration half-step */
    for (i=0;i<N;i++) {
      rx[i]+=vx[i]*dt+0.5*dt2*fx[i];
      ry[i]+=vy[i]*dt+0.5*dt2*fy[i];
      rz[i]+=vz[i]*dt+0.5*dt2*fz[i];
      vx[i]+=0.5*dt*fx[i];
      vy[i]+=0.5*dt*fy[i];
      vz[i]+=0.5*dt*fz[i];
      /* Apply periodic boundary conditions */
      if (rx[i]<0.0) { rx[i]+=L; ix[i]--; }
      if (rx[i]>L)   { rx[i]-=L; ix[i]++; }
      if (ry[i]<0.0) { ry[i]+=L; iy[i]--; }
      if (ry[i]>L)   { ry[i]-=L; iy[i]++; }
      if (rz[i]<0.0) { rz[i]+=L; iz[i]--; }
      if (rz[i]>L)   { rz[i]-=L; iz[i]++; }
    }
    /* Calculate forces */
    PE = total_e(rx,ry,rz,fx,fy,fz,N,L,rc2,ecor,ecut,&vir,
		 ll,head,nc,ncx,ncy,ncz,rcx,rcy,rcz,use_n2);
      
    /* Second integration half-step */
    KE = 0.0;
    for (i=0;i<N;i++) {
      vx[i]+=0.5*dt*fx[i];
      vy[i]+=0.5*dt*fy[i];
      vz[i]+=0.5*dt*fz[i];
      KE+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
    }
    KE*=0.5;
    TE=PE+KE;
    if (fOut && !(s%fOut))
      fprintf(stdout,"%i %.5lf %.5lf %.5lf %.5lf %.5le %.5lf %.5lf\n",
	      s,s*dt,PE,KE,TE,(TE-TE0)/TE0,KE*2/3./N,rho*KE*2./3./N+vir/3.0/V);
    if (fSamp && !(s%fSamp)) {
      sprintf(fn,"%i.xyz",!strcmp(wrt_code_str,"a")?0:s);
      out=fopen(fn,wrt_code_str);
      xyz_out(out,rx,ry,rz,vx,vy,vz,ix,iy,iz,L,N,16,1,unfold);
      fclose(out);
    }
  }
  t1=clock();
  sec= difftime(t1,t0)/CLOCKS_PER_SEC;
  fprintf(stdout,"# timing: %.1lf sec, %.3lf p-u/s\n",
	  sec,N*s/sec);
}
