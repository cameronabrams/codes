/* colormap.h (c) 1997 Cam Abrams
 * University of California, Berkeley
 * 
 * colormap module for Geomview 1.6.1
 * ==================================
 *  Uniquely maps a real number in some given range to
 *  an rgb color in some given rgb range.
 */
 
#ifndef COLORMAP_H
#define COLORMAP_H

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef ONE_SIXTH
#define ONE_SIXTH 0.16666666666666666666666666667
#endif
#ifndef ONE_THIRD
#define ONE_THIRD 0.33333333333333333333333333333
#endif
#ifndef TWO_THIRDS
#define TWO_THIRDS 0.66666666666666666666666666666
#endif
#ifndef NULL
#define NULL 0
#endif

typedef struct rgb_color
{
    float r;
    float g;
    float b;
} color_t;

typedef enum { C_R, C_G, C_B } hue_t;

hue_t color_maxHue ( color_t *c );
/* returns which hue is maximum in c */

hue_t color_minHue ( color_t *c );
/* returns which hue is minimum in c */

float color_maxVal ( color_t *c );
/* returns maximum value in c */

float color_minVal ( color_t *c );
/* returns minimum value in c */

hue_t color_nextHue ( hue_t Hue, int flow);
/* returns next hue from H in direction flow */

float color_hueVal ( color_t *c,  hue_t Hue );
/* returns value of hue H in color c */

color_t *color_init ( color_t *res, float r, float g, float b );
/* initializes color pointed to by res to rgb values given */

float color_color2x ( color_t *c );
/* converts color pointed to by c into a unique fraction [0,1);
 * returns the fraction
 */

color_t *colorPtr (char * colorStr);
char * colorStr (color_t * c);

color_t *color_x2color (color_t *c,  float x);
/* converts fraction [0,1) into a unique color c;
 * returns pointer to color.
 */

void color_printf (color_t *c);
/* outputs rgb values of color c */

color_t *color_scalMultColor (color_t *c, double x);
/* multiplies each value in c by scalar value x;
 * returns c.
 */

color_t *color_mapFloat ( float f0, color_t *c0, float f1, color_t *c1, 
			 float f, color_t *c, int flow);
/* maps float f (in range f0,f1) to color c (in range c0,c1) with
 * flow direction 'flow' (+1 == r->g->b->r);
 * returns pointer to color.
 */

#endif
