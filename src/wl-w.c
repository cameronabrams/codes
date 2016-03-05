#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <gsl/gsl_rng.h>

/***  wanglandau.c 
 *
 *
 *    (c) 2004 cameron abrams
 *    drexel university
 *    department of chemical engineering
 *
 *    philadelphia
 * 
 *    Implements the Wang-Landau flat histogram method for either
 *    an Ising or Potts model in two dimensions, a la 
 *    Wang&Landau, PRL 86:2050 (2001) and PRE 64:056101 (2001).
 *
 *
 ***/

enum { POTTS, ISING, NULL_SCHEME };
char * scheme_name [NULL_SCHEME] = {"POTTS", "ISING"};

int get_ebin ( int e, int n, int Scheme ) {
  /* Both:  e_0 = 2*n
     Ising has n-1 occupiable levels:  
              -e0,-e0+8,-e0+12,...,-4,0,4,...e0-12,e0-8,e0
     Potts has n occupiable levels: 
              as in -e0,-e0+4,-e0+6,...,-6,-4,-2,0
  */
  int e_0 = 2*n;
  int ee; 
  ee = e+e_0;
  if (ee) {ee/=2; if (Scheme==ISING) ee/=2; ee-=1;}
  if (Scheme==ISING && e==e_0) ee--;
  return ee;
}

int get_energy ( int i, int n, int Scheme ) {
  int e_0 = 2*n, e;
  if (!i) return -e_0;
  else {
    e=(i+1)*(Scheme==ISING?4:2)-e_0;
    if (Scheme==ISING && i==(n-2)) e=e_0;
  }
  return e;
}

int E ( int ** F, int L, int Scheme ) {
  int i,j;
  int energy=0;
  for (i=0;i<L;i++) {
    for (j=0;j<L;j++) {
      if (Scheme==POTTS) energy-=((int)(F[i][j]==F[i][(j+1)%L])+
				  (int)(F[i][j]==F[(i+1)%L][j]));
      else energy+=F[i][j]*F[i][(j+1)%L]+F[i][j]*F[(i+1)%L][j];
    }
  }
  return energy;
}

int dE ( int ** F, int L, int Scheme, int i, int j, int Q ) {
  if (Scheme==POTTS) {
    double energy = (int)(F[i][j]==F[i][(j+1)%L])
      + (int)(F[i][j]==F[(i+1)%L][j])
      + (int)(F[i][j]==F[i][j?(j-1):(L-1)])
      + (int)(F[i][j]==F[i?(i-1):(L-1)][j]);
    energy -= (int)(Q==F[i][(j+1)%L])
      + (int)(Q==F[(i+1)%L][j])
      + (int)(Q==F[i][j?(j-1):(L-1)])
      + (int)(Q==F[i?(i-1):(L-1)][j]);
    return energy;
  }
  else  return -2*(F[i][j])*(F[i?(i-1):(L-1)][j]+F[(i+1)%L][j]+
		      F[i][j?(j-1):(L-1)]+F[i][(j+1)%L]);
}

int HistFlatnessTest ( int * H, int n, double tol, int h_min ) {
  int i;

  double mean;
  double min,diff,rm;

  /* if no tolerance, use the minimal condition that every bin
     has at least one visit */
  if (tol==0.0) {
    rm=0;
    for (i=0;!rm&&i<n;i++) rm = rm || (H[i]==0);
    if (rm) fprintf(stderr," e(%i) is empty.\n",i);
    else fprintf(stderr," all bins hit.\n");
    return rm;
  }
  if (h_min) {
    rm=0;
    for (i=0;!rm&&i<n;i++) rm = rm || (H[i]<h_min);
    //if (rm) fprintf(stderr," e(%i) has less than %i hits.\n",i,h_min);
    return rm;
  }

  min = 1.e9;
  mean=0.0;
  for (i=0;i<n;i++) mean+=H[i];
  mean/=n;
  for (i=0;i<n;i++) {
    diff=((double)(H[i]))/mean;
    if (diff<min) min=diff;
  }

  fprintf(stderr,"av %.5le mn %.5le\n",mean,min);fflush(stderr);
  if (min < (1.0-tol)) return 1;
  else return 0;
  
}

void read_dosfile( char * idosfn, int * L, int * Scheme, 
		   int * nelevels, int * segment,
		   double * lnf, double ** lng) {
  FILE * ifp = fopen(idosfn,"r");
  double en, lng_i, h, f;
  int i, e, N;
  char ln[255];
  char * scheme_str;
  int argc=0;
  char * argv[255], *p;
  char tpln[255];

  if (ifp) {
    fgets(tpln,255,ifp);
    p=tpln;
    while (*p) {
      argv[argc++]=p;
      while (!isspace(*p)&&*p!='\n') p++;
      if (*p) {
	*p='\0';
	p++;
	while (isspace(*p)&&*p!='\n') p++;
      }
    }      
    for (i=0;i<argc;i++) {
      if (!strcmp(argv[i],"f")) sscanf(argv[++i],"%le",&f);
      else if (!strcmp(argv[i],"segment")) *segment=atoi(argv[++i]);
      else if (!strcmp(argv[i],"Scheme")) scheme_str=argv[++i];
      else if (!strcmp(argv[i],"L")) *L=atoi(argv[++i]);
    }
    N=(*L)*(*L);
    *Scheme=NULL_SCHEME;
    for (i=0;i<NULL_SCHEME;i++) {
      if (!strcmp(scheme_str,scheme_name[i])) *Scheme=i;
    }
    if (*Scheme==NULL_SCHEME) {
      fprintf(stderr,"# error: Scheme name %s is not recognized\n",scheme_str);
      exit(-1);
    }
    *nelevels=*Scheme==POTTS?(N):(N-1);
    
    *lng=(double*)calloc(*nelevels,sizeof(double));

    fprintf(stderr,"# restart from %s:  segment %i f %.10le"
	    " Scheme %s L %i\n",
	    idosfn,*segment,f,scheme_name[*Scheme],*L);

    *lnf=log(f);
    while (fgets(ln,255,ifp)) {
      sscanf(ln,"%i %i %lf %le %le\n",&i,&e,&en,&h,&lng_i);
      if (i>=*nelevels) {
	fprintf(stderr,"# Error;  file %s contains too many entries\n" 
		"Expecting %i\n",
		idosfn,*nelevels);
	exit(-1);
      }
      (*lng)[i]=lng_i;
    }
    fclose(ifp);
  }
  else {
    fprintf(stderr,"# error: cannot read %s\n",idosfn);
    exit(-1);
  }
}

void write_data (char * fn, int segment, double lnf, int N,
		 int Scheme, int * ehist, double * lng, int nelevels) {
  FILE * fp=fopen(fn,"w");
  int i;

  if (fp) {
    fprintf(fp,"# wl segment %i f %.10le Scheme %s L %i\n",
	    segment,exp(lnf),scheme_name[Scheme],sqrt(N));
   for (i=0;i<nelevels;i++) {
      fprintf(fp,"%i %i %.5lf %.5le %.5le\n",
	      i,get_energy(i,N,Scheme),
	      ((double)get_energy(i,N,Scheme))/N,
	      (double)(ehist[i]),
	      lng[i]);
    }
    fclose(fp);
  }
  else {
    fprintf(stderr,"# error: could not open output %s\n",
	    fn);
    exit(-1);
  }
}

int main ( int argc, char * argv[] ) {

  int L=20; // size of one side of square domain
  int Q=10; // number of possible values per spin
  int ** F = NULL, F_save;
  int N;
  int i,j;

  int energy, new_energy, old_energy, d_energy;
  int e_min, e_max, new_spin;

  int nelevels, * ehist;
  double * lng;
  double lnf = 1.0, prob, x;
  double f_tol = 1.0e-8;
  double h_tol = 0.20;
  double fac, update_factor = 0.5, T=1.0;

  char * idosfn = NULL;

  int segment, old_ebin, new_ebin, nFlips, otf_write=0, ground=0;
  
  char fn[255], * tail_fn = NULL;
  FILE * fp, * tail_fp = NULL;

  int nSweeps, nAcc, inSegment, keep_going, nTotalSweeps=0;
  int sampleInterval = 10000, h_min = 0;

  int Scheme=POTTS;
  int rfac,rshf;

  int quiet = 0;
  int doWangLandau = 1;
  int doNVT = 0, maxSweeps = 1000;

  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);

  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-L")) L=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-Q")) Q=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-si")) sampleInterval=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-I")) Scheme=ISING;
    else if (!strcmp(argv[i],"-htol")) h_tol = atof(argv[++i]);
    else if (!strcmp(argv[i],"-hmin")) h_min = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-q")) quiet=1;
    else if (!strcmp(argv[i],"-restart")) idosfn=argv[++i];
    else if (!strcmp(argv[i],"-otf_write")) otf_write=1;
    else if (!strcmp(argv[i],"-tail_fn")) tail_fn=argv[++i];
    else if (!strcmp(argv[i],"-f")) update_factor = atof(argv[++i]);
    else if (!strcmp(argv[i],"-g")) ground = 1;
    else if (!strcmp(argv[i],"-nvt")) {
      doWangLandau=0;
      doNVT=1;
      T = atof(argv[++i]);
      maxSweeps = atoi(argv[++i]);
    }
  }

  if (Scheme==ISING&&Q!=2) {
    fprintf(stderr,"# Warning:  Ising scheme, setting Q = 2\n");
    Q=2;
  }

  if (tail_fn) {
    tail_fp=fopen(tail_fn,"w");
    fprintf(tail_fp,"# sweep energy\n");
  }

  // establish the initial density of states 
  if (idosfn) read_dosfile(idosfn,&L,&Scheme,&nelevels,&segment,&lnf,&lng);
  else {
    N=L*L;
    nelevels=Scheme==POTTS?(N):(N-1);
    segment=0;
    lng=(double*)calloc(nelevels,sizeof(double));
    memset(lng,0,nelevels*sizeof(double));
  }
  N=L*L;

  /* Allocate the lattice, and initialize */
  F=(int**)malloc(L*sizeof(int*));
  for (i=0;i<L;i++) F[i]=(int*)calloc(L,sizeof(int));
  
  /* Allocate the histogram */
  ehist=(int*)calloc(nelevels,sizeof(int));

  rfac=(Scheme==ISING)?2:1;
  rshf=(Scheme==ISING)?-1:0;

  fprintf(stderr,"# %s Sampling of the %s System\n"
	  "# N = %i, Q = %i, sample_interval = %i\n",
	  doWangLandau?"Wang-Landau":"NVT",
	  scheme_name[Scheme],N,Q,sampleInterval);

  if (doWangLandau) {
    if (!h_min)
      fprintf(stderr,"# flatness at %.3lf%%\n",100.0*(1.0-h_tol));
    else
      fprintf(stderr,"# flatness after a minimum of %i hits\n",h_min);
  }
  if (doNVT)
    fprintf(stderr,"# Temperature %.5lf, maxSweeps %i\n",T,maxSweeps);

  keep_going=1;
  while (keep_going) {

    if (doWangLandau) {
      fprintf(stderr,"# segment %i, f %.10le...\n",segment,exp(lnf));
      fflush(stderr);
    }

    // zero the histogram
    memset(ehist,0,nelevels*sizeof(int));

    // assign all spins randomly to initialize this run segment
    for (i=0;i<L;i++) {
      for (j=0;j<L;j++) {
	F[i][j]=ground?1:(rfac*(int)gsl_rng_uniform_int(r,Q)+rshf);
      }
    }
    // update the energy bin and DOS bin for this new state
    old_energy=E(F,L,Scheme);
    old_ebin=get_ebin(old_energy,N,Scheme);
    ehist[old_ebin]++;
    lng[old_ebin]+=lnf;

    // perform sweeps to execute random walk in energy
    // and stop when a predefined histogram flatness
    // criterion is met
    nSweeps=0;
    inSegment=1;
    while (inSegment) {

      nFlips=0;
      nAcc=0;
      while (nFlips < N) { // one sweep is N flip attempts
	// randomly select a spin
	i=(int)gsl_rng_uniform_int(r,L);
	j=(int)gsl_rng_uniform_int(r,L);
	if (Scheme==ISING) 
	  new_spin = -F[i][j];
	else
	  new_spin = (int)gsl_rng_uniform_int(r,Q);
	d_energy = dE(F,L,Scheme,i,j,new_spin);
	new_ebin=get_ebin(old_energy+d_energy,N,Scheme);
	// calc acceptance criterion
	if (doWangLandau)
	  prob=exp(lng[old_ebin]-lng[new_ebin]);
	else if (doNVT)
	  prob=exp(-d_energy/T);
	if (prob>1.0) prob=1.0;

	// select a uniform random variate between 0 and 1
	x=gsl_rng_uniform(r);
	if (x<prob) { // accept the move
	  ehist[new_ebin]++;
	  lng[new_ebin]+=lnf;
	  nAcc++;
	  old_energy+=d_energy;
	  old_ebin=new_ebin;
	  F[i][j]=new_spin;
	}
	else {
	  ehist[old_ebin]++;
	  lng[old_ebin]+=lnf;
	}

	nFlips++;
      }
      // perform test
      if (doNVT) inSegment=nSweeps<maxSweeps;
      if (nSweeps && !(nSweeps%sampleInterval)) {
	//fprintf(stderr,"# sw %i ",nSweeps);fflush(stderr);
	if (doWangLandau)
	  inSegment=HistFlatnessTest(ehist,nelevels,h_tol,h_min);
	else if (doNVT) {
	  sprintf(fn,"nvt%i.dat",nSweeps);
	  write_data(fn,segment,lnf,N,Scheme,ehist,lng,nelevels);
	}
	if (otf_write) 
	  write_data("tmp_hist",segment,lnf,N,Scheme,ehist,lng,nelevels);
      }
      else {
	if (doWangLandau) inSegment=1;
      }
      if (tail_fp) fprintf(tail_fp,"%i %.10le\n",nSweeps,old_energy);
      nSweeps++;
    }
    nTotalSweeps+=nSweeps-1;
    fprintf(stderr,"# end of segment %i : %i total sweeps\n",
	    segment,nTotalSweeps);
    // end of segment
    // renormalize density of states; same for Ising or Potts
    fac=log((double)Q)-lng[0];
   
    for (i=0;i<nelevels;i++) lng[i]+=fac;

    // output histogram and density of states
    sprintf(fn,"wl%i.dat",segment);
    write_data(fn,segment,lnf,N,Scheme,ehist,lng,nelevels);

    // update factor 
    lnf *= update_factor;
    if (doWangLandau && lnf > f_tol) keep_going=1;
    else keep_going=0;
    segment++;
  }
  if (tail_fp) {
    fprintf(tail_fp,"# end of data\n");
    fclose(tail_fp);
  }
  fprintf(stderr,"# Program ends.\n");

}
