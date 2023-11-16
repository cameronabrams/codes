#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

/* Computes and outputs fourier transform of y vs x data
   reads from stdin
   cameron abrams
 */

int main (int argc, char * argv[]) {
  int i, n = 5001;
  double data[n];  // y
  double r[n];     // x

  gsl_fft_real_wavetable * real;
  gsl_fft_halfcomplex_wavetable * hc;
  gsl_fft_real_workspace * work;

  /* read in data */
  char ln[255];
  i=0;
  while (i<5000&&fgets(ln,255,stdin)) {
    if (ln[0]!='#') {
      sscanf(ln,"%lf %lf",&r[i],&data[i]);
    }
    i++;
  }
  n = i;
  //  fprintf(stderr,"# got %i datapoints\n",i);fflush(stderr);

  work = gsl_fft_real_workspace_alloc (n);
  real = gsl_fft_real_wavetable_alloc (n);

  gsl_fft_real_transform (data, 1, n,
			  real, work);

  gsl_fft_real_wavetable_free (real);
  

  /* data is now transformed; output */
  for (i=0;i<n;i++) {
    fprintf(stdout,"%i %.5le\n",i,data[i]);
  }
  return 0;
}


