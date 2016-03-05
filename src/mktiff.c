/* mktiff.c
 *
 * generates a single tiff image of a field z(x,y)
 * 
 * the tiff library must be installed
 *
 * compile as "gcc -o mktiff mktiff.c -lm -ltiff
 *
 * Cameron Abrams
 * 2005-2006
 */

#include <stdio.h>
#include <stdlib.h>
#include <tiffio.h>
#include <string.h>
#include <math.h>

#ifndef MAX_ROWS
#define MAX_ROWS 1000
#endif
#ifndef MAX_COLUMNS
#define MAX_COLUMNS 1000
#endif

typedef struct VMDCOLORMAP {
  double x;
  char r,g,b;
} vmdcm;

const vmdcm VMDCOLORS[6]={{-10,0,0,255},{-11,255,0,0},{-12,128,128,128},
			  {-37,205,0,205},{-20,0,205,205},{-15,0xbb,0x77,0x77}};

int getVMDCOLOR(double x) {
  int i;
  for (i=0;i<6;i++) {
    if (x==VMDCOLORS[i].x) {
      fprintf(stderr,"found vmd color %lf in database\n",x);
      return i;
    }
  }
  return -1;
}

typedef struct FIELD {
  double * x;
  int N;
  double * y;
  int M;
  double ** f;
} field;

field * new_field (int N, int M) {
  int i;
  field * f = (field*)malloc(sizeof(field));
  f->N=N;
  f->M=M;
  f->y=(double*)malloc(N*sizeof(double));
  f->x=(double*)malloc(M*sizeof(double));
  f->f=(double**)malloc(N*sizeof(double*));
  for (i=0;i<N;i++) {
    f->f[i]=(double*)malloc(M*sizeof(double));
  }
  return f;
}

void free_field ( field * f ) {
  int i;
  free((double*)f->y);
  free((double*)f->x);
  for (i=0;i<f->N;i++) {
    free((double*)f->f[i]);
  }
  free((double**)f->f);
  free((field*)f);
}

field * realloc_field ( field * f, int N, int M ) {
  int i;
  f->N=N;
  f->M=M;
  f->y=(double*)realloc((double*)f->y,M*sizeof(double));
  f->x=(double*)realloc((double*)f->x,N*sizeof(double));
  f->f=(double**)realloc((double**)f->f,N*sizeof(double*));
  for (i=0;i<N;i++) {
     f->f[i]=(double*)realloc((double*)f->f[i],M*sizeof(double));
  }
  return f;
}


field * a_leftof_b_composite ( field * a, field * b ) {
  int i,oldM;
  if (a->N!=b->N) {
    fprintf(stderr,"ERROR: cannot left-of composite fields with different numbers of rows.\n");
    return NULL;
  }
  oldM=a->M;
  a=realloc_field(a,a->N,a->M+b->M);
  for (i=0;i<a->N;i++) {
    memcpy(&a->f[i][oldM],&b->f[i][0],b->M*sizeof(double));
  }
  memcpy(&a->x[oldM],&b->x[0],b->M*sizeof(double));
  return a;
}

field * horizontal_grade ( field * f, int dir ) {
  if (f) {
    int i,j;
    int h=f->N; // number of rows, i
    int w=f->M; // number of columns, j

    for (i=0;i<h;i++) {
      for (j=0;j<w;j++) {
	f->f[i][j]=((double)(j))/w;
      }
    }
    return f;
  }     
  return NULL;
}

field * vertical_grade ( field * f, int dir ) {
  if (f) {
    int i,j;
    int h=f->N; // number of rows, i
    int w=f->M; // number of columns, j

    for (i=0;i<h;i++) {
      for (j=0;j<w;j++) {
	if (dir==1) f->f[i][j]=((double)(i))/h;
	if (dir==-1) f->f[h-i-1][j]=((double)(i))/h;
      }
    }
    return f;
  }     
  return NULL;
}

field * uniform_field ( field * f, double v ) {
  if (f) {
    int i,j;
    for (i=0;i<f->N;i++) {
      for (j=0;j<f->M;j++) {
	f->f[i][j]=v;
      }
    }
    return f;
  }
  return NULL;
}

field * read_field_gnuplot ( char * fn ) {
  FILE * fp = fopen(fn,"r");
  if (fp) {
    field * f = new_field(MAX_ROWS,MAX_COLUMNS);
    char ln[500];
    int i,j,maxj;
    /* gnuplot format: data is blocked with blank lines between "rows" of points; i.e.,
       each block is all points for a row */
    i=0;
    j=0;
    maxj=0;
    while (fgets(ln,500,fp)) {
      if (ln[0]!='%'&&ln[0]!='#') {
	if (ln[0]=='\n') {
	  i++;
	  if (j>maxj) maxj=j;
	  j=0;
	} else {
	  sscanf(ln,"%lf %lf %lf",&f->x[i],&f->y[j],&f->f[i][j]);
	  j++;
	}
      }
    }
    fclose(fp);
    f=realloc_field(f,i,maxj);
    fprintf(stderr,"INFO: read_field_gnuplot read %i rows and %i columns\n",f->N,f->M);
    return f;
  }
  return NULL;

}

field * unit_scale_field ( field * f, double * desired_fmin, double * desired_fmax ) {
  int i,j;
  double minf=1.e10,maxf=-1.e10,span;

  if (desired_fmin && desired_fmax) {
    minf=*desired_fmin;
    maxf=*desired_fmax;
  } else {

    for (i=0;i<f->N;i++) {
      for (j=0;j<f->M;j++) {
	if (f->f[i][j]>maxf) maxf=f->f[i][j];
	else if (f->f[i][j]<minf) minf=f->f[i][j];
      }
    }
  }

  span=maxf-minf;
  for (i=0;i<f->N;i++) {
    for (j=0;j<f->M;j++) {
      f->f[i][j]=(f->f[i][j]-minf)/span;
    }
  }
  return f;
}

int make_tiff_image ( char * image, int w, int h, int spp, char * fn ) {

  TIFF *out = TIFFOpen(fn,"w");
  tsize_t linebytes = spp*w;
  int row;
  unsigned char * buf = NULL;

  TIFFSetField(out,TIFFTAG_IMAGEWIDTH,w);
  TIFFSetField(out,TIFFTAG_IMAGELENGTH,h);
  TIFFSetField(out,TIFFTAG_SAMPLESPERPIXEL,spp);
  TIFFSetField(out,TIFFTAG_BITSPERSAMPLE,8);
  TIFFSetField(out,TIFFTAG_ORIENTATION,ORIENTATION_TOPLEFT);
  TIFFSetField(out,TIFFTAG_PLANARCONFIG,PLANARCONFIG_CONTIG);
  TIFFSetField(out,TIFFTAG_PHOTOMETRIC,PHOTOMETRIC_RGB);

  if (TIFFScanlineSize(out)<linebytes) 
    buf=(unsigned char*)_TIFFmalloc(linebytes);
  else
    buf=(unsigned char*)_TIFFmalloc(TIFFScanlineSize(out));

  for (row=0;row<h;row++) {
    memcpy(buf,&image[row*linebytes],linebytes);
    if (TIFFWriteScanline(out,buf,row,0)<0) break;
  }

  TIFFClose(out);
  if (buf) _TIFFfree(buf);
  fprintf(stderr,"INFO: created %s.\n",fn);
  return 0;
}

int get_hsv ( double x, double hsv[], double hmin, double hmax ) {
  double y=1-x;
  if (x>(1-hmin)) {
    hsv[1]=1; hsv[2]=1.0;
    hsv[0]=y<hmax?hmax:y;
  }
  else {
    hsv[0]=hmin;//0.8;
    hsv[1]=1;
    hsv[2]=x/(1-hmin);
  }
  return 0;
}

int hsv2rgb ( double hsv[], char rgb[] ) {
  int var_i;
  double var_1,var_2,var_3;
  double var_r,var_g,var_b,var_h;

  if ( hsv[1]==0 ) {
    rgb[0] = hsv[2] * 255;
    rgb[1] = hsv[2] * 255;
    rgb[2] = hsv[2] * 255;
  }
  else {
    var_h = hsv[0] * 6;
    var_i = (int)( var_h );             //Or ... var_i = floor( var_h )
    var_1 = hsv[2] * ( 1 - hsv[1] );
    var_2 = hsv[2] * ( 1 - hsv[1] * ( var_h - var_i ) );
    var_3 = hsv[2] * ( 1 - hsv[1] * ( 1 - ( var_h - var_i ) ) );

    if      ( var_i == 0 ) { var_r = hsv[2]     ; var_g = var_3 ; var_b = var_1; }
    else if ( var_i == 1 ) { var_r = var_2 ; var_g = hsv[2]     ; var_b = var_1; }
    else if ( var_i == 2 ) { var_r = var_1 ; var_g = hsv[2]     ; var_b = var_3; }
    else if ( var_i == 3 ) { var_r = var_1 ; var_g = var_2 ; var_b = hsv[2];     }
    else if ( var_i == 4 ) { var_r = var_3 ; var_g = var_1 ; var_b = hsv[2];     }
    else                   { var_r = hsv[2]     ; var_g = var_1 ; var_b = var_2; }
 
    rgb[0] = var_r * 255;
    rgb[1] = var_g * 255;
    rgb[2] = var_b * 255;
  }
  return 0;
}


int get_rgb ( char * c, double x, double hmin, double hmax ) {

  double hsv[3];
  c[0]=c[1]=c[2]=c[3]=0;

  if (x==-1) { // white
    c[0]=c[1]=c[2]=c[3]=255;
    return 0;
  }
  else if (x<0.0) { // VMD COLOR
    int i;
    fprintf(stderr,"detected vmd color index %lf\n",x);
    i=getVMDCOLOR(x);
    c[0]=VMDCOLORS[i].r;
    c[1]=VMDCOLORS[i].g;
    c[2]=VMDCOLORS[i].b;
    c[3]=255;
    return 0;
  }

  if (x<0.0) x=0.0;
  if (x>1.0) x=1.0;

  get_hsv(x,hsv,hmin,hmax);
  hsv2rgb(hsv,c);

  return 0;
}

int make_raw_image ( char * image, int lb, int spp, 
		     field * f, double hmin, double hmax ) {
  int i,j;
  char clr[4];
  int h=f->N;
  int w=f->M;
  
  memset(&image[0],0,lb*h);

  for (i=0;i<h;i++) {
    for (j=0;j<w;j++) {
      //      fprintf(stderr,"# %i %i %.5lf\n",i,j,f->f[i][j]);fflush(stderr);
      if (!get_rgb(clr,f->f[i][j],hmin,hmax))
	memcpy(&image[i*lb+j*spp],clr,spp);
      else {
	fprintf(stderr,"#error mapping color from %.4f at %i %i\n",f->f[i][j],i,j);
	exit(-1);
      }
    }
  }
  return 0;
}


struct scat {
  char fn[50];
  char clr[3];
  double * x, * y;
  int n;
};
typedef struct scat TScat;

int read_scat (TScat * s) {
  FILE * fp = fopen(s->fn,"r");
  s->n=0;
  s->x=s->y=NULL;
  if (fp) {
    int nl=0;
    char ln[255];
    int i;

    while (fgets(ln,255,fp)) if (ln[0]!='#') nl++;
    s->x=(double*)malloc(nl*sizeof(double));
    s->y=(double*)malloc(nl*sizeof(double));
    rewind(fp);
    //fprintf(stderr,"# %s %i lines\n",s->fn,nl);
    i=0;
    while (fgets(ln,255,fp)) {
      if (ln[0]!='#') {
	sscanf(ln,"%lf %lf",&s->x[i],&s->y[i]);
	//fprintf(stderr,"#in %i %lf %lf\n",i,s->x[i],s->y[i]);
	//fflush(stderr);
	i++;
      }
    }
    s->n=i;
    fclose(fp);
    //fprintf(stderr,"# read %s; %i pts\n",s->fn,i);
    fflush(stderr);
  }
  return 0;
}

int raw_image_scat ( char * image, int lb, int spp, int w, int h, TScat * s, int ns) {
  int i,j;
  int ypx,xpx;

  for (i=0;i<ns;i++) {
    //fprintf(stderr,"# drawing scat %i...\n",i);fflush(stderr);
    for (j=0;j<s[i].n;j++) {
      ypx=(int)rint((float)s[i].y[j]);
      xpx=(int)rint(s[i].x[j]);
      //fprintf(stderr,"# %i %lf %lf %i %i (%i)\n",i,s[i].x[j],s[i].y[j],xpx,ypx,(h-ypx-1));
      memcpy(&image[(h-ypx-1)*lb+xpx*spp],s[i].clr,spp);
    }
  }
  return 0;
}

int main ( int argc, char * argv[] ) {

  int sampleperpixel = 3;
  int height=401,width=401;
  int i,p;
  tsize_t linebytes = 0;
  int logscale=0;
  char * ifn=NULL;
  char * tfn="my.tif";
  char *image=NULL;
  field * f=NULL;
  field * g=NULL, *gg=NULL;
  double ** field, maxv, minv;
  double a_max, a_min;
  double hmin = 0.8, hmax = 0.0;
  int a_set=0;
  int make_colorbar=0;
  int nscat=0;
  TScat * S=(TScat*)malloc(10*sizeof(TScat));

  for (i=1;i<argc;i++) {
    if (!strcmp("-i",argv[i])) ifn=argv[++i];
    else if (!strcmp("-o",argv[i])) tfn=argv[++i];
    else if (!strcmp("-r",argv[i])) {
      sscanf(argv[++i],"%lf,%lf",&a_min,&a_max);
      a_set=1;
    }
    else if (!strcmp("-hue",argv[i])) {
      sscanf(argv[++i],"%lf,%lf",&hmin,&hmax);
    }
    else if (!strcmp("-log",argv[i])) logscale=1;
    else if (!strcmp("-mcb",argv[i])) {
      make_colorbar=1;
    }
    else if (!strcmp("-geom",argv[i])) sscanf(argv[++i],"%i,%i",&width,&height);
    else if (!strcmp("-scat",argv[i])) {
      if (nscat==10) {
	fprintf(stderr,"#error; cannot handle more than %i"
		" scatter plots. Skipping.\n",10);
      }
      else {
	int r,g,b;
	sscanf(argv[++i],"%i,%i,%i,%s",&r,&g,&b,S[nscat].fn);
	S[nscat].clr[0]=(char)r;
	S[nscat].clr[1]=(char)g;
	S[nscat].clr[2]=(char)b;
	read_scat(&S[nscat]);
	nscat++;
      }
    }
  }


  if (ifn) {
    f=read_field_gnuplot(ifn);
    
    f=unit_scale_field(f,a_set?&a_min:NULL,a_set?&a_max:NULL);

    g=new_field(f->N,2);
    g=uniform_field(g,-1.0);
    g=a_leftof_b_composite(g,f);
    free_field(f);
    f=g;

    g=new_field(f->N,2);
    g=read_field_gnuplot("sequencebar.dat");
    //    g=vertical_grade(g,-1);
    g=a_leftof_b_composite(g,f);
    free_field(f);
    f=g;

    //    fprintf(stderr,"#INFO: unit_scale_field returns.\n");fflush(stderr);
    image=(char*)malloc(f->M*f->N*sampleperpixel*sizeof(char)); 
    linebytes = sampleperpixel*f->M;
    make_raw_image(image,linebytes,sampleperpixel,f,hmin,hmax);
    //    fprintf(stderr,"#INFO: make_raw_image returns.\n");fflush(stderr);
    make_tiff_image(image,f->M,f->N,sampleperpixel,tfn);
    free_field(f);

/*     image=(char*)malloc(width*height*sampleperpixel*sizeof(char)); */
/*     field=(double**)malloc(height*sizeof(double*)); */
/*     for (i=0;i<height;i++) field[i]=(double*)malloc(width*sizeof(double)); */

/*     if (read_field(ifn,width,height,field,&minv,&maxv)) { */
/*       fprintf(stderr,"# error reading %s\n",ifn); */
/*       exit(-1); */
/*     } */
/*     if (a_set) { */
/*       minv=a_min; */
/*       maxv=a_max; */
/*     } */

/*     make_raw_image(image,linebytes,sampleperpixel,field,width,height,minv,maxv,logscale,hmin,hmax); */
    
/*     if (nscat) raw_image_scat(image,linebytes,sampleperpixel,width,height,S,nscat); */

/*     make_tiff_image(image,width,height,sampleperpixel,tfn); */
    
/*     free(image); */
/*     free(field); */
/*   } */

    if (make_colorbar) {
      int width=20;
      int height=f->N;
      linebytes = sampleperpixel*width;
      image=(char*)malloc(width*height*sampleperpixel*sizeof(char));
      //    cbfield=make_colorbar_field(width,height);
      
      //make_raw_image(image,linebytes,sampleperpixel,cbfield,hmin,hmax);
      
      make_tiff_image(image,width,height,sampleperpixel,"colorbar.tif");
      
      free(image);
      
    }
    fprintf(stderr,"INFO: program ends.\n");
  }
  else {
    int h=200;
    int w=400;
    linebytes = sampleperpixel*w;
    image=(char*)malloc(h*w*sampleperpixel*sizeof(char));
    f=new_field(h,w);
    fprintf(stderr,"New field %i %i\n",h,w);fflush(stderr);
    f=horizontal_grade(f,1);
    fprintf(stderr,"horizontal grade %i %i\n",h,w);fflush(stderr);
    make_raw_image(image,linebytes,sampleperpixel,f,hmin,hmax);
    fprintf(stderr,"raw image %i %i\n",h,w);fflush(stderr);
    make_tiff_image(image,w,h,sampleperpixel,"horizontal.tif");
    fprintf(stderr,"tiff image %i %i\n",h,w);fflush(stderr);
    f=vertical_grade(f,1);
    make_raw_image(image,linebytes,sampleperpixel,f,hmin,hmax);
    make_tiff_image(image,w,h,sampleperpixel,"vertical.tif");
    fprintf(stderr,"INFO: test ends.\n");    
  }
}
