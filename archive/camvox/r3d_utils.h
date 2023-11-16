#ifndef R3D_UTILS_H
#define R3D_UTILS_H

#include <stdio.h>
#include "point.h"
#include "colormap.h"
#include "dblmat.h"
#define MAX_ARC 50
#define MAX_ARC_ARROW_CONE 20
#ifndef SQRT_3
#define SQRT_3	    1.732050807568877293527446341505
#endif
#ifndef SQRT_3_2
#define SQRT_3_2    0.866025403784438646763723170752
#endif
enum {U, X, Y, Z};
#define MAXROTS 15

void r3d_header (FILE * fp, double rots[], int rot_order[], 
		     double zoom, ptPtr shift, double span);
void r3d_tri (FILE * fp, ptPtr p0, ptPtr p1, ptPtr p2, color_t color);
void r3d_cone (FILE * fp, ptPtr b, ptPtr t, double r, int n, color_t color);
void r3d_cyl (FILE * fp, ptPtr p0, ptPtr p1, double r, color_t color);
void r3d_3darrow (FILE * fp, ptPtr p0, ptPtr p1, double r, 
		double ahl, double ahw, int n, color_t color);
void r3d_sphere (FILE * fp, ptPtr p0, double r, color_t color);
void r3d_disc (FILE * fp, ptPtr c, ptPtr axis, double r, 
	    int n, color_t color, ptPtr arc[]);
void r3d_ring (FILE * fp, ptPtr c, ptPtr axis, double r, double t, int n, 
		color_t color, ptPtr arc[]);
void r3d_halfring (FILE * fp, ptPtr c, ptPtr axis, double r, double t, int n, 
		color_t color, ptPtr arc[]);
void r3d_quad (FILE * fp, ptPtr p0, ptPtr p1, ptPtr p2, ptPtr p3, color_t color);
void r3d_arc (FILE * fp, ptPtr p0, ptPtr p1, ptPtr p2,  
	     int n, double liner, color_t color, ptPtr arc[]);
void r3d_calcArcPoints (ptPtr arc[], int n, ptPtr c, ptPtr axis, double r);

void r3d_volumeElement_cubic (FILE * fp, ptPtr c, 
			      double d, color_t color, double clarity);
void r3d_volumeElement_polar (FILE * fp, ptPtr c, 
			      double dr, double dt, double dp,
			      color_t color, double clarity);

void r3d_transparency_on (FILE * fp, double clarity);
void r3d_transparency_off (FILE * fp);




#endif
