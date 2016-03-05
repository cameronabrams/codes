/* 
   Microcanonical Molecular Dynamics simulation of 
   a diamond-cubic silicon crystal using the
   Stillinger-Weber potential in a periodic boundary

   Reference:  Stillinger & Weber, "Computer simulation of local order in 
   condensed phases of silicon", Phys. Rev. B, V31, No. 8, 5262, 1985.

   Cameron F. Abrams

   Written for the course CHE 800-002, Molecular Simulation
   Spring 0304

   compile using "gcc -o mdswsi mdswsi.c -lm -lgsl"
   (assumes the GNU Scientific Library is installed)

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
#include <gsl/gsl_randist.h>

/* Prints usage information */
void usage ( void ) {
  fprintf(stdout,"mdswsi usage:\n");
  fprintf(stdout,"mdswdi [options]\n\n");
  fprintf(stdout,"Options:\n");
  fprintf(stdout,"\t -N [integer]\t\tNumber of atoms\n");
  fprintf(stdout,"\t -rho [real]\t\tNumber density\n");
  fprintf(stdout,"\t -dt [real]\t\tTime step\n");
  fprintf(stdout,"\t -rc [real]\t\tCutoff radius\n");
  fprintf(stdout,"\t -ns [real]\t\tNumber of integration steps\n");
  fprintf(stdout,"\t -so       \t\tShort-form output (unused)\n");
  fprintf(stdout,"\t -T0 [real]\t\tInitial temperature\n");
  fprintf(stdout,"\t -Tb [real]\t\tBath temperature, Andersen thermostat\n");
  fprintf(stdout,"\t -nu [real]\t\tCollision frequency, Andersen thermostat\n");
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
	      int * ix, int * iy, int * iz, 
	      double Lx, double Ly, double Lz,
	      int N, int z, int put_vel, int unfold) {
  int i;
  double cx=0.0,cy=0.0,cz=0.0;

  fprintf(fp,"%i %i\n\n",N,put_vel);
  for (i=0;i<N;i++) {
    fprintf(fp,"%i %.8lf %.8lf %.8lf ",z,
	    rx[i]+(unfold?(ix[i]*Lx-cx):0.0),
	    ry[i]+(unfold?(iy[i]*Ly-cy):0.0),
	    rz[i]+(unfold?(iz[i]*Lz-cz):0.0));
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


#define ANGSTROMS_PER_SIGMA	2.0951
#define EV_PER_EPSILON		2.1678
#define PICOSECS_PER_TAU	0.076634   
#define G_PER_M			4.66362659e-23 
#define AMU_PER_M		28.085 
#define GperCC_PER_MperCUSIG	5.071185526
#define KELVIN_PER_SWT		25156.73798
/* 
Subset 1.a:  Two-body model potential parameters 
*/ 
#define A_		7.0495562770
#define B_		0.6022245584
#define P_		4
#define Q_		0
#define a_		1.8

/* Subset 1.b:  Three-body model potential parameters */
#define LAM_		21
#define GAM_ 		1.2

/* Set 2:  convenient expressions involving parameters */
#define NP_		-4
#define NP_1_		-5
#define NQ_		0
#define NQ_1_		-1
#define BP_		2.408898232
#define a2_		3.24
#define ONE_THIRD       0.3333333333333333333333

double total_e ( double * rx, double * ry, double * rz, 
		 double * fx, double * fy, double * fz, 
		 int N, double Lx, double Ly, double Lz) {

   int i,j,k;
   double dxij, dyij, dzij, r2ij, rij;
   double dxik, dyik, dzik, r2ik, rik;
   double dxjk, dyjk, dzjk, r2jk, rjk;
   double vcuij, evcuij, vlj, res, v2, C2;
   double rhoij, rij_a, rik_a, rjk_a;
   double vcuik, rhoik, vcujk, rhojk;
   double vcui, vcuj, vcuk;
   double hjik, hijk, hikj;
   double hij, hik, hjk;
   double ri, rj, rk;
   double si, sj, sk;
   double ei, ej, ek;
   double cosjik, cosijk, cosikj;
   double cosjik_3, cosijk_3, cosikj_3;
   double eiNjik_x, eiNjik_y, eiNjik_z;
   double eiNkij_x, eiNkij_y, eiNkij_z;
   double ejNijk_x, ejNijk_y, ejNijk_z;
   double ejNkji_x, ejNkji_y, ejNkji_z;
   double ekNikj_x, ekNikj_y, ekNikj_z;
   double ekNjki_x, ekNjki_y, ekNjki_z;
   double hRij_x, hRij_y, hRij_z;
   double hRik_x, hRik_y, hRik_z;
   double hRjk_x, hRjk_y, hRjk_z;
   double dF3i_x, dF3i_y, dF3i_z;
   double dF3j_x, dF3j_y, dF3j_z;
   double dF3k_x, dF3k_y, dF3k_z;
   double e2 = 0.0, e3 = 0.0;
   double hLx=Lx/2.0, hLy=Ly/2.0, hLz=Lz/2.0; 
   double f;

   /* Zero the forces */
   for (i=0;i<N;i++) {
     fx[i]=fy[i]=fz[i]=0.0;
   }
   
   for (i=0;i<(N-1);i++) {
     for (j=i+1;j<N;j++) {
       rij=vcuij=0.0;
       dxij  = (rx[i]-rx[j]);
       dyij  = (ry[i]-ry[j]);
       dzij  = (rz[i]-rz[j]);
       if (dxij>hLx)       dxij-=Lx;
       else if (dxij<-hLx) dxij+=Lx;
       if (dyij>hLy)       dyij-=Ly;
       else if (dyij<-hLy) dyij+=Ly;
       if (dzij>hLz)       dzij-=Lz;
       else if (dzij<-hLz) dzij+=Lz;
       r2ij = dxij*dxij + dyij*dyij + dzij*dzij;
       rij = sqrt(r2ij);
       if (r2ij<a2_) {
	 vcuij = 1.0/(rij - a_);
	 evcuij = 0.0;
	 if (vcuij > -30.0) evcuij = exp(vcuij);
	 if (evcuij != evcuij) evcuij = 0.0;
	 res = 1.0;
	 if (NQ_) res = pow(rij,NQ_);
	 vlj = B_ * pow(rij,NP_) - res;
	 v2 = A_ * vlj * evcuij;
	 C2  = BP_*pow(rij,NP_1_);
	 if (Q_) C2 -= Q_*pow(rij,NQ_1_);
	 C2 *= A_*evcuij;
	 C2 += v2*vcuij*vcuij;
	 C2 /= rij;

	 fx[i] += dxij*C2;
	 fx[j] -= dxij*C2;
	 fy[i] += dyij*C2;
	 fy[j] -= dyij*C2;
	 fz[i] += dzij*C2;
	 fz[j] -= dzij*C2;

	 e2+=v2;
       }
#ifdef THREE_BODY
       rij_a=rhoij=vcuij=0.0;
       if (rij<a_) {
	 rij_a = rij - a_;
	 rhoij = rij*rij_a*rij_a;
	 vcuij = GAM_/rij_a;
	 if (vcuij > -30) vcuij = exp(vcuij);
	 else vcuij = 0.0;
	 if (vcuij != vcuij) vcuij = 0.0;
       }
       for (k=j+1;k<N;k++) { 

	 hjik = hijk = hikj = 0.0;
	 hij = hik = hjk = 0.0;
	 ri = rj = rk = 0.0;
	 si = sj = sk = 0.0;
	 ei = ej = ek = 0.0;
	 cosjik = cosijk = cosikj = 0.0;
	 cosjik_3 = cosijk_3 = cosikj_3 = 0.0;

	 eiNjik_x = 0.0; eiNjik_y = 0.0; eiNjik_z = 0.0;
	 eiNkij_x = 0.0; eiNkij_y = 0.0; eiNkij_z = 0.0;
	 ejNijk_x = 0.0; ejNijk_y = 0.0; ejNijk_z = 0.0;
	 ejNkji_x = 0.0; ejNkji_y = 0.0; ejNkji_z = 0.0;
	 ekNikj_x = 0.0; ekNikj_y = 0.0; ekNikj_z = 0.0;
	 ekNjki_x = 0.0; ekNjki_y = 0.0; ekNjki_z = 0.0;
	 
	 rik_a = vcuik = rhoik = 0.0;
	 rjk_a = vcujk = rhojk = 0.0;

	 dxik  = (rx[i]-rx[k]);
	 dyik  = (ry[i]-ry[k]);
	 dzik  = (rz[i]-rz[k]);
	 if (dxik>hLx)       dxik-=Lx;
	 else if (dxik<-hLx) dxik+=Lx;
	 if (dyik>hLy)       dyik-=Ly;
	 else if (dyik<-hLy) dyik+=Ly;
	 if (dzik>hLz)       dzik-=Lz;
	 else if (dzik<-hLz) dzik+=Lz;
	 r2ik = dxik*dxik + dyik*dyik + dzik*dzik;

	 dxjk  = (rx[j]-rx[k]);
	 dyjk  = (ry[j]-ry[k]);
	 dzjk  = (rz[j]-rz[k]);
	 if (dxjk>hLx)       dxjk-=Lx;
	 else if (dxjk<-hLx) dxjk+=Lx;
	 if (dyjk>hLy)       dyjk-=Ly;
	 else if (dyjk<-hLy) dyjk+=Ly;
	 if (dzjk>hLz)       dzjk-=Lz;
	 else if (dzjk<-hLz) dzjk+=Lz;
	 r2jk = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk;
	 
	 if (r2ik<a2_) {
	   rik = sqrt(r2ik);
	   rik_a = rik - a_;
	   rhoik = rik*rik_a*rik_a;
	   vcuik = GAM_/rik_a;
	   if (vcuik > -30.0) vcuik = exp(vcuik);
	   else vcuik = 0.0;
	   if (vcuik != vcuik) vcuik = 0.0;
	 }
	 if (r2jk<a2_) {
	   rjk = sqrt(r2jk);
	   rjk_a = rjk - a_;
	   rhojk = rjk*rjk_a*rjk_a;
	   vcujk = GAM_/rjk_a;
	   if (vcujk > -30.0) vcujk = exp(vcujk);
	   else vcujk = 0.0;
	   if (vcujk != vcujk) vcujk = 0.0;
	 }
	 vcui = vcuij*vcuik;
	 if (vcui) {
	   ri = rij*rik;
	   cosjik = (dxij*dxik+dyij*dyik+dzij*dzik)/ri;
	   cosjik_3 = cosjik + ONE_THIRD;
	   si = LAM_*vcui*cosjik_3;
	   ei = 2*si/ri;
	   hjik = si*cosjik_3;
	 }
	 vcuj = vcuij*vcujk;
	 if (vcuj) {
	   rj = rij*rjk;
	   cosijk = -(dxij*dxjk+dyij*dyjk+dzij*dzjk)/rj;
	   cosijk_3 = cosijk + ONE_THIRD;
	   sj = LAM_*vcuj*cosijk_3;
	   ej = 2*sj/rj;
	   hijk = sj*cosijk_3;
	 }
	 vcuk = vcuik*vcujk;
	 if (vcuk) {
	   rk = rik*rjk;
	   cosikj = (dxjk*dxik+dyjk*dyik+dzjk*dzik)/rk;
	   cosikj_3 = cosikj + ONE_THIRD;
	   sk = LAM_*vcuk*cosikj_3;
	   ek = 2*sk/rk;
	   hikj = sk*cosikj_3;
	 }

#if DIAG
	 printf("  i %i j %i k %i\n",i,j,k);
	 printf("  rij %.5lf rik %.5lf rjk %.5lf rc %.5lf\n",
		rij,rik,rjk,a_);
	 printf("  vcuij %.5lf vcuik %.5lf vcujk %.5lf\n",
		vcuij,vcuik,vcujk);
	 printf("  jik {%i|%i/%i}\n", j,i,k);
	 printf("          vcu%i: %.10lf\n",i,vcui);
	 if (!vcui) printf("          NONE were calculated:\n");
	 printf("            r%i: %.10lf\n",i,ri);
	 printf("  cos(%i,%i,%i): %.10lf\n",j,i,k,cosjik);
	 printf("cos(%i,%i,%i)_3: %.10lf\n",j,i,k,cosjik_3);
	 printf("           si%i: %.10lf\n",i,si);
	 printf("           ei%i: %.10lf\n",i,ei);
	 printf("  ->h(%i,%i,%i): %.10lf\n",j,i,k,hjik);
	 printf("  ijk {%i|%i/%i}\n", i,j,k);
	 printf("          vcu%i: %.10lf\n",j,vcuj);
	 if (!vcuj) printf("          NONE were calculated:\n");
	 printf("            r%i: %.10lf\n",j,rj);
	 printf("  cos(%i,%i,%i): %.10lf\n",i,j,k,cosijk);
	 printf("cos(%i,%i,%i)_3: %.10lf\n",i,j,k,cosijk_3);
	 printf("           sj%i: %.10lf\n",j,sj);
	 printf("           ej%i: %.10lf\n",j,ej);
	 printf("  ->h(%i,%i,%i): %.10lf\n",i,j,k,hijk);
	 printf("  ikj {%i|%i/%i}\n", i,k,j);
	 printf("          vcu%i: %.10lf\n",k,vcuk);
	 if (!vcuk) printf("          NONE were calculated:\n");
	 printf("            r%i: %.10lf\n",k,rk);
	 printf("  cos(%i,%i,%i): %.10lf\n",i,k,j,cosikj);
	 printf("cos(%i,%i,%i)_3: %.10lf\n",i,k,j,cosikj_3);
	 printf("           sk%i: %.10lf\n",k,sk);
	 printf("           ek%i: %.10lf\n",k,ek);
	 printf("  ->h(%i,%i,%i): %.10lf\n",i,k,j,hikj);                   
#endif

	 eiNjik_x = ei*(dxik - rik/rij*cosjik*dxij);
	 eiNjik_y = ei*(dyik - rik/rij*cosjik*dyij);
	 eiNjik_z = ei*(dzik - rik/rij*cosjik*dzij);
	 
	 eiNkij_x = ei*(dxij - rij/rik*cosjik*dxik);
	 eiNkij_y = ei*(dyij - rij/rik*cosjik*dyik);
	 eiNkij_z = ei*(dzij - rij/rik*cosjik*dzik);
	 	 
	 ejNijk_x = ej*(dxjk + rjk/rij*cosijk*dxij);
	 ejNijk_y = ej*(dyjk + rjk/rij*cosijk*dyij);
	 ejNijk_z = ej*(dzjk + rjk/rij*cosijk*dzij);
	 
	 ejNkji_x = ej*(-dxij - rij/rjk*cosijk*dxjk);
	 ejNkji_y = ej*(-dyij - rij/rjk*cosijk*dyjk);
	 ejNkji_z = ej*(-dzij - rij/rjk*cosijk*dzjk);
	 
	 ekNikj_x = ek*(-dxjk + rjk/rik*cosikj*dxik);
	 ekNikj_y = ek*(-dyjk + rjk/rik*cosikj*dyik);
	 ekNikj_z = ek*(-dzjk + rjk/rik*cosikj*dzik);

	 ekNjki_x = ek*(-dxik + rik/rjk*cosikj*dxjk);
	 ekNjki_y = ek*(-dyik + rik/rjk*cosikj*dyjk);
	 ekNjki_z = ek*(-dzik + rik/rjk*cosikj*dzjk);
	 
	 if (rhoij) hij = GAM_*(hijk+hjik)/rhoij;
	 if (rhoik) hik = GAM_*(hjik+hikj)/rhoik;
	 if (rhojk) hjk = GAM_*(hikj+hijk)/rhojk;

	 hRij_x = hij * dxij;
	 hRij_y = hij * dyij;
	 hRij_z = hij * dzij;
	 hRik_x = hik * dxik;
	 hRik_y = hik * dyik;
	 hRik_z = hik * dzik;
	 hRjk_x = hjk * dxjk;
	 hRjk_y = hjk * dyjk;
	 hRjk_z = hjk * dzjk;

	 dF3i_x = hRij_x + hRik_x - eiNjik_x - eiNkij_x + ejNijk_x + ekNikj_x;
	 dF3i_y = hRij_y + hRik_y - eiNjik_y - eiNkij_y + ejNijk_y + ekNikj_y;
	 dF3i_z = hRij_z + hRik_z - eiNjik_z - eiNkij_z + ejNijk_z + ekNikj_z;
	 
	 dF3j_x = -hRij_x + hRjk_x - ejNijk_x - ejNkji_x + eiNjik_x + ekNjki_x;
	 dF3j_y = -hRij_y + hRjk_y - ejNijk_y - ejNkji_y + eiNjik_y + ekNjki_y;
	 dF3j_z = -hRij_z + hRjk_z - ejNijk_z - ejNkji_z + eiNjik_z + ekNjki_z;

	 dF3k_x = -hRjk_x - hRik_x - ekNikj_x - ekNjki_x + eiNkij_x + ejNkji_x;
	 dF3k_y = -hRjk_y - hRik_y - ekNikj_y - ekNjki_y + eiNkij_y + ejNkji_y;
	 dF3k_z = -hRjk_z - hRik_z - ekNikj_z - ekNjki_z + eiNkij_z + ejNkji_z;

	 fx[i] += dF3i_x;
	 fy[i] += dF3i_y;
	 fz[i] += dF3i_z;
	 fx[j] += dF3j_x;
	 fy[j] += dF3j_y;
	 fz[j] += dF3j_z;
	 fx[k] += dF3k_x;
	 fy[k] += dF3k_y;
	 fz[k] += dF3k_z;
	 
	 e3 += (hijk + hjik + hikj);
       }
#endif
     }
   }
   //   fprintf(stderr,"# e2 %.5lf e3 %.5lf\n",e2,e3);
   return e2 + e3;
}

/* Initialize particle positions by assigning them
   on a diamond cubic lattice, ncx unit cells in x-dir,
   ncy unit cells in y-dir, ncz unit cells in z-dir,
   then scale positions to achieve a given box size 
   and thereby, volume, and density */
void init ( double * rx, double * ry, double * rz,
	    double * vx, double * vy, double * vz,
	    int * ix, int * iy, int * iz,
	    int ncx, int ncy, int ncz, double rho,
	    int * N, double * Lx, double * Ly, double * Lz,
	    gsl_rng * r, double T0,
	    double * KE, char * icf) {
  int i,iix,iiy,iiz;
  double a0, b0, c0;
  int ic, il;
  double cl_x, cl_y, cl_z;
  double hcl_x, hcl_y, hcl_z;
  double qcl_x, qcl_y, qcl_z;
  double tqcl_x, tqcl_y, tqcl_z;

  double cmvx=0.0,cmvy=0.0,cmvz=0.0,f3;
  double T, fac;
  int n3=2;
  int vel_ok=0;
  
  /* If icf has a value, assume it is the name of a file containing
     the input configuration in XYZ format */
  if (icf) {
    FILE * fp = fopen(icf,"r");
    if (fp) vel_ok = xyz_in(fp,rx,ry,rz,vx,vy,vz,N);
    else {
      fprintf(stderr,"# error: could not read %s\n",icf);
      exit(-1);
    }
  }
  /* Assign particles on a diamond cubic lattice */
  else {

    a0 = 1.0;
    b0 = 1.0;
    c0 = 1.0;
    f3 = 8/(a0*b0*c0)/rho;
    f3 = pow(f3,ONE_THIRD);
    a0*=f3;
    b0*=f3;
    c0*=f3;
    fprintf(stderr,"# a0 %.5lf\n",a0);

    cl_x = a0;
    cl_y = b0;
    cl_z = c0;
    hcl_x = cl_x*0.5;
    hcl_y = cl_y*0.5;
    hcl_z = cl_z*0.5;
    qcl_x = hcl_x*0.5;
    qcl_y = hcl_y*0.5;
    qcl_z = hcl_z*0.5;
    tqcl_x = qcl_x*3.0;
    tqcl_y = qcl_y*3.0;
    tqcl_z = qcl_z*3.0;

    /* assign the unit cell */

    rx[0] = 0.0;    ry[0] = 0.0;    rz[0] = 0.0;
    rx[1] = hcl_x;  ry[1] = hcl_y;  rz[1] = 0.0;
    rx[2] = 0.0;    ry[2] = hcl_y;  rz[2] = hcl_z;
    rx[3] = hcl_x;  ry[3] = 0.0;    rz[3] = hcl_z;
    rx[4] = qcl_x;  ry[4] = qcl_y;  rz[4] = qcl_z;
    rx[5] = tqcl_x; ry[5] = tqcl_y; rz[5] = qcl_z;
    rx[6] = qcl_x;  ry[6] = tqcl_y; rz[6] = tqcl_z;
    rx[7] = tqcl_x; ry[7] = qcl_y;  rz[7] = tqcl_z;
    
    (*Lx) = ncx*a0;
    (*Ly) = ncy*b0;
    (*Lz) = ncz*c0;

    (*N) = 8 * ncx*ncy*ncz;
    il=0;
    for (iiz=0;iiz<ncz;iiz++) {
      for (iiy=0;iiy<ncy;iiy++) {
	for (iix=0;iix<ncx;iix++) {
	  for (ic = 0; ic < 8; ic++) {
	    rx[ic+il] = rx[ic] + a0*iix;
	    ry[ic+il] = ry[ic] + b0*iiy;
	    rz[ic+il] = rz[ic] + c0*iiz;
	  }
	  il+=8;
	}
      }
    }
  }
  /* If no velocities yet assigned, randomly pick some */
  if (!vel_ok) {
    for (i=0;i<(*N);i++) {
      vx[i]=gsl_ran_gaussian(r,1.0);
      vy[i]=gsl_ran_gaussian(r,1.0);
      vz[i]=gsl_ran_gaussian(r,1.0);
    }
  }
  /* Take away any center-of-mass drift; compute initial KE */
  for (i=0;i<(*N);i++) {
    cmvx+=vx[i];
    cmvy+=vy[i];
    cmvz+=vz[i];
  }
  (*KE)=0;
  for (i=0;i<(*N);i++) {
    vx[i]-=cmvx/(*N);
    vy[i]-=cmvy/(*N);
    vz[i]-=cmvz/(*N);
    (*KE)+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
  }
  (*KE)*=0.5;
  /* if T0 is specified, scale velocities */
  if (T0>0.0) {
    T=(*KE)/(*N)*2./3.;
    fac=sqrt(T0/T);
    (*KE)=0;
    for (i=0;i<(*N);i++) {
      vx[i]*=fac;
      vy[i]*=fac;
      vz[i]*=fac;
      (*KE)+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
    }
    (*KE)*=0.5;
  }
  /* Initialize periodic boundary crossing counter arrays */
  memset(ix,0,(*N)*sizeof(int));
  memset(iy,0,(*N)*sizeof(int));
  memset(iz,0,(*N)*sizeof(int));
}

int main ( int argc, char * argv[] ) {

  double * rx, * ry, * rz;
  double * vx, * vy, * vz;
  double * fx, * fy, * fz;
  int * ix, * iy, * iz;
  int N=216,c,a;
  double Lx=0.0,Ly=0.0,Lz=0.0;
  double rho=0.138, T0=1.0, Tb=1.0, nu = 1.0, vir, vir_old, vir_sum, pcor, V;
  double PE, KE, TE, ecor, ecut, TE0, sigma;
  double rr3,dt=0.001, dt2;
  int i,j,s;
  int nSteps = 10, fSamp=100;
  int short_out=0;
  int use_e_corr=0;
  int unfold = 0;
  int ncx = 3, ncy = 3, ncz = 3;
  char fn[20];
  FILE * out;
  char * wrt_code_str = "w";
  char * init_cfg_file = NULL;
  int veloc = 1;
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned long int Seed = 23410981;

  /* Here we parse the command line arguments;  If
   you add an option, document it in the usage() function! */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-rho")) rho=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dt")) dt=atof(argv[++i]);
    else if (!strcmp(argv[i],"-ns")) nSteps = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-so")) short_out=1;
    else if (!strcmp(argv[i],"-T0")) T0=atof(argv[++i]);
    else if (!strcmp(argv[i],"-Tb")) Tb=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nu")) nu=atof(argv[++i]);
    else if (!strcmp(argv[i],"-fs")) fSamp=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-sf")) wrt_code_str = argv[++i];
    else if (!strcmp(argv[i],"-icf")) init_cfg_file = argv[++i];
    else if (!strcmp(argv[i],"-ecorr")) use_e_corr = 1;
    else if (!strcmp(argv[i],"-seed")) Seed = (unsigned long)atoi(argv[++i]);
    else if (!strcmp(argv[i],"-uf")) unfold = 1;
    else if (!strcmp(argv[i],"+v")) veloc=0;
    else if (!strcmp(argv[i],"-nc")) {
      sscanf(argv[++i],"%i,%i,%i",&ncx,&ncy,&ncz);
    }
    else if (!strcmp(argv[i],"-h")) {
      usage(); exit(0);
    }
    else {
      fprintf(stderr,"Error: Command-line argument '%s' not recognized.\n",
	      argv[i]);
      exit(-1);
    }
  }

  /* compute the squared time step */
  dt2=dt*dt;

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
  init(rx,ry,rz,vx,vy,vz,
       ix,iy,iz,
       ncx,ncy,ncz,
       rho,&N,&Lx,&Ly,&Lz,r,T0,&KE,init_cfg_file);
  fprintf(stderr,"# init done\n");
  /* Output some initial information */
  fprintf(stdout,"# NVE/T MD Simulation of Stillinger-Weber Silicon\n");
  fprintf(stdout,"# L = (%.5lf,%.5lf,%.5lf); rho = %.5lf; N = %i\n",
	  Lx,Ly,Lz,rho,N);
  fprintf(stdout,"# nSteps %i, seed %d, dt %.5lf\n",
	  nSteps,Seed,dt);
  fprintf(stdout,"# Andersen thermostat, nu %.5lf\n",nu); 

  /* Compute sigma for Andersen thermostat */
  sigma = sqrt(Tb);

  sprintf(fn,"init.xyz");
  out=fopen(fn,"w");
  xyz_out(out,rx,ry,rz,vx,vy,vz,ix,iy,iz,Lx,Ly,Lz,N,16,veloc,unfold);
  fclose(out);

  PE = total_e(rx,ry,rz,fx,fy,fz,N,Lx,Ly,Lz);
  TE0=PE+KE;
  
  fprintf(stdout,"# step PE KE TE drift T\n");

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
      if (rx[i]<0.0)  { rx[i]+=Lx; ix[i]--; }
      if (rx[i]>Lx)   { rx[i]-=Lx; ix[i]++; }
      if (ry[i]<0.0)  { ry[i]+=Ly; iy[i]--; }
      if (ry[i]>Ly)   { ry[i]-=Ly; iy[i]++; }
      if (rz[i]<0.0)  { rz[i]+=Lz; iz[i]--; }
      if (rz[i]>Lz)   { rz[i]-=Lz; iz[i]++; }
    }
    /* Calculate forces */
    PE = total_e(rx,ry,rz,fx,fy,fz,N,Lx,Ly,Lz);
      
    /* Second integration half-step */
    KE = 0.0;
    for (i=0;i<N;i++) {
      vx[i]+=0.5*dt*fx[i];
      vy[i]+=0.5*dt*fy[i];
      vz[i]+=0.5*dt*fz[i];
      KE+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
    }
    KE*=0.5;
    /* Andersen thermostat */
    KE = 0.0;
    for (i=0;i<N;i++) {
      if (gsl_rng_uniform(r) < nu*dt) {
	vx[i]=gsl_ran_gaussian(r,sigma);
	vy[i]=gsl_ran_gaussian(r,sigma);
	vz[i]=gsl_ran_gaussian(r,sigma);
      }
      KE+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
    }
    KE*=0.5;
    TE=PE+KE;
    fprintf(stdout,"%i %.5lf %.5lf %.5lf %.5lf %.5le %.5lf\n",
	    s,s*dt,PE/N,KE/N,TE/N,(TE-TE0)/TE0,KE*2/3./N);
    if (!(s%fSamp)) {
      sprintf(fn,"%i.xyz",!strcmp(wrt_code_str,"a")?0:s);
      out=fopen(fn,wrt_code_str);
      xyz_out(out,rx,ry,rz,vx,vy,vz,ix,iy,iz,Lx,Ly,Lz,N,16,veloc,unfold);
      fclose(out);
    }
  }
}
