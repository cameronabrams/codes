#include <tiffio.h>
#include <stdio.h>
#include <stdlib.h>

/* Image data is stored 4 bytes ("samples") per pixel in a character array.
   An image of w pixels by h pixels has 4wh bytes, and to address
   the j'th pixel of row i (both starting with 0), you 
   go to element [i*w*4+j*4].  The data you write at that locations
   is 4 bytes: red (0-254), green (0-254), blue (0-254), and a zero.
   Greys have all the same values for r, g, b; black is (0,0,0); white
   is (1,1,1).

   Let's say you want to generate a 100 x 100 image.

   int h = 100, w = 100;
   char * image;
   char pxl[4];
   ...

   image = (char*)malloc(4*h*w*sizeof(char));

   Make image white to start
   pxl[0]=pxl[1]=pxl[2]=(char)1;  pxl[4]=0;
   for (i=0;i<h;i++) for (j=0;j<w;j++) memcpy(image+i*w*4+j*4,clr,4*sizeof(char));
   
   ...
   
   (code to generate image)

   
   Your program will have to generate the byte array for an image; then
   a simple call to make_tiff_image will generate the image file.

   Cameron 2006
*/

/* make_tiff_image 

   image : character array of image data
   w, h  : image width and height in pixels
   spp   : samples (bytes) per pixel (usually 4; 1 each for r, g, b, and channel)
   fn    : name of output file
*/
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
  return 0;
}
