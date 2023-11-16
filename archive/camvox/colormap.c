/* colormap.c	(c) 1997 Cam Abrams
 * University of California,  Berkeley
 * Dept. of Chemical Engineering
 * Graves Group
 */

#include "colormap.h"

float hue_max_ = 1.0;
float hue_min_ = 0.0;


#include "cam_colors.cf"

hue_t color_maxHue ( color_t *c )
{
    if (c->r >= c->b && c->r >= c->g) return C_R;
    if (c->g >= c->b && c->g >= c->r) return C_G;
    if (c->b >= c->r && c->b >= c->g) return C_B;
    return C_R;
}

hue_t color_minHue ( color_t *c )
{
    if (c->r <= c->b && c->r <= c->g) return C_R;
    if (c->g <= c->b && c->g <= c->r) return C_G;
    if (c->b <= c->r && c->b <= c->g) return C_B;
    return C_R;
}

float color_maxVal ( color_t *c )
{
    return color_hueVal(c, color_maxHue(c));
}

float color_minVal ( color_t *c )
{
    return color_hueVal(c, color_minHue(c));
}

color_t *color_init ( color_t *res, float r, float g, float b )
{
    if (!res) return NULL;
    res->r = r;
    res->g = g;
    res->b = b;
    return res;
}

color_t *colorPtr (char * colorStr)
{
    if (!colorStr) return NULL;
    #include "cam_str2clr.cf"
    return NULL;
}

char tmp_colorStr[255];
char * colorStr (color_t * c)
{
    char * rv = tmp_colorStr;
    rv[0] = '\0';
    #include "cam_clr2str.cf"
    return rv;
}

hue_t color_nextHue ( hue_t Hue, int flow)
{
    if (flow != 1 && flow != -1) return Hue;
    if (Hue != C_R && Hue != C_B && Hue != C_G) return Hue;
    if (Hue == C_R && flow == 1) return C_G;
    if (Hue == C_R && flow == -1) return C_B;
    if (Hue == C_G && flow == 1) return C_B;
    if (Hue == C_G && flow == -1) return C_R;
    if (Hue == C_B && flow == 1) return C_R;
    if (Hue == C_B && flow == -1) return C_G;
    return C_G;
}

float color_hueVal ( color_t *c,  hue_t Hue )
{
    if (Hue != C_R && Hue != C_B && Hue != C_G) return 0.0;
    if (Hue == C_R) return c->r;
    if (Hue == C_G) return c->g;
    if (Hue == C_B) return c->b;
    return C_G;
}

float color_color2x ( color_t *c )
{
    if (color_maxHue(c) == C_R)
    {
	float f = c->g - c->b;
	if (f < 0.0) return (1.0 + ONE_SIXTH * f) * (hue_max_ - hue_min_);
	return ONE_SIXTH * f * (hue_max_ - hue_min_);  
    }
    if (color_maxHue(c) == C_G)
	return (ONE_THIRD + ONE_SIXTH * (c->b - c->r)) * (hue_max_ - hue_min_);
  
    if (color_maxHue(c) == C_B) 
	return (TWO_THIRDS + ONE_SIXTH * (c->r - c->g)) * (hue_max_ - hue_min_);

    return 0.0;  
}

color_t *color_x2color (color_t *c,  float x)
{
    double i,  *ip = &i;
    color_init(c, 0.0, 0.0, 0.0);
    x = modf(x, ip);
    if (x < 0.0) x += 1.0;
    if (x > (1.0 - ONE_SIXTH)) x -= 1.0;
    if (x <= ONE_SIXTH || x >= (1.0 - ONE_SIXTH)) {
	c->r = hue_max_;
	if (x < 0.0) {
	    c->g = hue_min_;
	    c->b = -6*x*(hue_max_ - hue_min_);
	}
	else {
	    c->b = hue_min_;
	    c->g = 6*x*(hue_max_ - hue_min_);
	}
    }
    else if (x <= 0.5) {
	c->g = hue_max_;
	if (x < ONE_THIRD) {
	    c->b = hue_min_;
	    c->r = 6*(ONE_THIRD - x)*(hue_max_ - hue_min_);
	}
	else {
	    c->r = hue_min_;
	    c->b = 6*(x - ONE_THIRD)*(hue_max_ - hue_min_);
	}
    }
    else {
	c->b = hue_max_;
	if (x < TWO_THIRDS) {
	    c->r = hue_min_;
	    c->g = 6*(TWO_THIRDS - x)*(hue_max_ - hue_min_);
	}
	else {
	    c->g = hue_min_;
	    c->r = 6*(x - TWO_THIRDS)*(hue_max_ - hue_min_);
	}
    }
    return c;
}

void color_printf (color_t *c) {
    if (!c) printf("0.0 0.0 0.0");
    printf("%.3f %.3f %.3f", c->r, c->g, c->b);
}

color_t *color_scalMultColor (color_t *c, double x)
{
    return color_init(c, c->r*x, c->g*x, c->b*x);
}

color_t *color_mapFloat ( float f0, color_t *c0, float f1, color_t *c1, 
			 float f, color_t *c, int flow)
{
    float M1, m1, M0, m0;
    float m = (f1 - f)/(f1 - f0);
    float x0 = color_color2x(c0);
    float x1 = color_color2x(c1);
    float dx = x1 - x0;
    float xc = 0.0;
    
    M1 = color_maxVal(c1);
    M0 = color_maxVal(c0);
    m1 = color_minVal(c1);
    m0 = color_minVal(c0);
    
    /* check constraints on c0 and c1 */
    if ((M1 != M0) || (m1 != m0)) return color_init(c, 0.0, 0.0, 0.0);
    
    /* if the float is not within the range (f0,f1) then assign the slope
     * so that if f < f0, xc = x0, and if f > f1, xc = x1.
     */
    
    if (f > f0 && f < f1)
    {
	if (flow == 1) xc = m*dx + x0;
	else xc = x1 - m*dx;
    }
    else if (f <= f0)
	xc = x0;
    else if (f >= f1)
	xc = x1;
    
    
    return color_x2color (c, xc);

}
