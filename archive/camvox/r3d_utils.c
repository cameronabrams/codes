#include "r3d_utils.h"
#include "cam_colors.ext"

char * rot_label[4] = 
{
    "U", "X", "Y", "Z"
};

/* render Header parameters and their default values */
double ImageBufferSize_ = 5.0;
int xTiles_ = 32, yTiles_ = 32;
int xPixPerTile_ = 18, yPixPerTile_ = 18;
int Scheme_ = 3;
color_t * bgClr = &black;
int Shadows_ = 1;
int Phong_ = 25;
double SecLtContr_ = 0.15;
double AmbLtContr_ = 0.05;
double SpecRefComp_ = 0.25;
double Eyepos_ = 4.0;
pt MLP = {1, 1, 1, 0, 0, NULL, NULL};
const ptPtr MainLtPos = &MLP;
short transparency_ = 0;
short RotateLightWithObject_ = 0;
double tmpMat1[4][4];
double tmpMat2[4][4];

void r3d_header (FILE * fp, double rots[], int rot_order[],
		 double zoom, ptPtr shift, double span)
{
    int i = 0, j = 0, r = 0;
    double trans[4][4] = 
    {
	1.0, 0.0, 0.0, 0.0, 
	0.0, 1.0, 0.0, 0.0, 
	0.0, 0.0, 1.0, 0.0, 
	0.0, 0.0, 0.0, 1.0
    };
    ptPtr tpt = ptPtr_Init(0, 0, 0);
    double cosT = 0.0, sinT = 0.0;
    if (!fp) return;
    if (!span) span = 1.0;
    for (r = 0; r < MAXROTS && rot_order[r]; r++)
    {
	dblmat_identity(tmpMat1, 4);
	dblmat_identity(tmpMat2, 4);
	sinT = sin(M_PI/180.0*rots[r]);
	cosT = cos(M_PI/180.0*rots[r]);
	switch (rot_order[r])
	{
	    case X: /* Rotation about X-axis:  yy, yz, zy, zz elements assigned */
		    tmpMat1[1][1] = cosT;	/* yy */
		    tmpMat1[1][2] = -sinT;	/* yz */
		    tmpMat1[2][1] = sinT;	/* zy */
		    tmpMat1[2][2] = cosT;	/* zz */
		    break;
	    case Y: /* Rotation about Y-axis:  xx, xz, zx, zz elements assigned */
		    tmpMat1[0][0] = cosT;	/* xx */
		    tmpMat1[0][2] = -sinT;	/* xz */
		    tmpMat1[2][0] = sinT;	/* zx */
		    tmpMat1[2][2] = cosT;	/* zz */
		    break;
	    case Z: /* Rotation about Z-axis:  xx, xy, yx, yy elements assigned */
		    tmpMat1[0][0] = cosT;	/* xx */
		    tmpMat1[0][1] = -sinT;	/* xy */
		    tmpMat1[1][0] = sinT;	/* yx */
		    tmpMat1[1][1] = cosT;	/* yy */
		    break;
	    default:
		    printf("ERROR!  Invalid rotation.\n");
		    break;
	}
	dblmat_copy(tmpMat2, trans, 4);
	dblmat_matmult(trans, tmpMat2, tmpMat1, 4);
    }
    if (shift)
    {
	trans[3][0] = shift->x/span;
	trans[3][1] = shift->y/span;
	trans[3][2] = shift->z/span;
    }
    if (zoom) trans[3][3] = 1.0/zoom;
    fprintf(fp, "%s\n", "tmp1");
    fprintf(fp, "%i %i     tiles in x,y\n", xTiles_, yTiles_);
    fprintf(fp, "%i %i     pixels (x,y) per tile\n", xPixPerTile_, yPixPerTile_);
    fprintf(fp, "%i        Scheme\n", Scheme_);
    fprintf(fp, "%.3lf %.3lf %.3lf  Background color\n", 
			    bgClr->r, bgClr->g, bgClr->b);
    fprintf(fp, "%s        shadows? T[yes] F[no]\n", (Shadows_ ? "T" : "F"));
    fprintf(fp, "%i        Phong power\n", Phong_);
    fprintf(fp, "%.2lf     secondary light contribution\n", SecLtContr_);
    fprintf(fp, "%.2lf     ambient light contribution\n", AmbLtContr_);
    fprintf(fp, "%.2lf     specular reflection component\n", SpecRefComp_);
    fprintf(fp, "%.2lf     eye position\n", Eyepos_);
    if (MainLtPos)
    {
	if (RotateLightWithObject_)
	{
	    fprintf(stderr, "Orbiting...");
	    tpt->x = MainLtPos->x*trans[0][0] + 
		     MainLtPos->y*trans[1][0] +
		     MainLtPos->z*trans[2][0];
	    tpt->y = MainLtPos->x*trans[0][1] + 
		     MainLtPos->y*trans[1][1] +
		     MainLtPos->z*trans[2][1];
	    tpt->z = MainLtPos->x*trans[0][2] + 
		     MainLtPos->y*trans[1][2] +
		     MainLtPos->z*trans[2][2];
	    MainLtPos->x = tpt->x;
	    MainLtPos->y = tpt->y;
	    MainLtPos->z = tpt->z;
	}
	fprintf(fp, "%.2lf %.2lf %.2lf  main light source position\n", 
		MainLtPos->x, MainLtPos->y, MainLtPos->z);	
    }
    else fprintf(fp, "1 1 1\tmain light source position\n");
    for (i = 0; i < 4; i++)
    {
	for (j = 0; j < 4; j++)
	{
	    fprintf(fp, "%.4lf ", trans[i][j]);
	}
	fprintf(fp, "\n");
    }
    fprintf(fp, "3         mixed objects\n");
    fprintf(fp, "*\n*\n*\n");
    
    
}

void r3d_tri (FILE * fp, ptPtr p0, ptPtr p1, ptPtr p2, color_t color)
{
    if (!fp) return;
    if (!p0 || !p1 || !p2) return;
    
    fprintf(fp, "1\n");
    fprintf(fp, "%.4lf %.4lf %.4lf  %.4lf %.4lf %.4lf  %.4lf %.4lf %.4lf  %.4lf %.4lf %.4lf \n", 
	p0->x, p0->y, p0->z, 
	p1->x, p1->y, p1->z, 
	p2->x, p2->y, p2->z, 
	color.r, color.g, color.b);
}

void r3d_quad (FILE * fp, ptPtr p0, ptPtr p1, ptPtr p2, ptPtr p3, color_t color)
{
    if (!fp) return;
    if (!p0 || !p1 || !p2 || !p3) return;
    r3d_tri(fp, p1, p0, p3, color);
    r3d_tri(fp, p1, p2, p3, color);
}

void r3d_sphere (FILE * fp, ptPtr p0, double r, color_t color)
{
    if (!fp) return;
    if (!p0 || !r) return;
    
    fprintf(fp, "2\n");
    fprintf(fp, "%.4lf %.4lf %.4lf  %.4lf %.4lf %.4lf %.4lf \n", 
	p0->x, p0->y, p0->z, 
	r, 
	color.r, color.g, color.b);
}

void r3d_cyl (FILE * fp, ptPtr p0, ptPtr p1, double r, color_t color)
{
    if (!fp) return;
    if (!p0 || !p1 || !r) return;
    
    fprintf(fp, "5\n");
    fprintf(fp, "%.4lf %.4lf %.4lf  %.4lf  %.4lf %.4lf %.4lf  %.4lf  %.4lf %.4lf %.4lf\n", 
	p0->x, p0->y, p0->z, r, 
	p1->x, p1->y, p1->z, 9.99, 
	color.r, color.g, color.b);
}


void r3d_disc (FILE * fp, ptPtr c, ptPtr axis, double r, int n, 
		color_t color, ptPtr arc[])
{
    int i = 0, j = 0;

    r3d_calcArcPoints(arc, n, c, axis, r);

    for (i = 0; i < n; i++)
    {
	r3d_tri(fp, c, arc[i], arc[(i+1)%n], color);
    }
}

void r3d_ring (FILE * fp, ptPtr c, ptPtr axis, double r, double t, int n, 
		color_t color, ptPtr arc[])
{
    int i = 0, j = 0;
    ptPtr arc_inner[MAX_ARC];
    for (i = 0; i < n; i++) arc_inner[i] = ptPtr_Init(0.0, 0.0, 0.0);
    
    r3d_calcArcPoints(arc, n, c, axis, r);
    r3d_calcArcPoints(arc_inner, n, c, axis, r-t);
    for (i = 0; i < n; i++)
    {
	r3d_quad(fp, arc[i], arc[(i+1)%n],
	    arc_inner[(i+1)%n], arc_inner[i], color);
    }
    for (i = 0; i < n; i++) freeptPtr(arc_inner[i]);
}

void r3d_halfring (FILE * fp, ptPtr c, ptPtr axis, double r, double t, int n, 
		color_t color, ptPtr arc[])
{
    int i = 0, j = 0;
    ptPtr arc_inner[MAX_ARC];
    for (i = 0; i < n; i++) arc_inner[i] = ptPtr_Init(0.0, 0.0, 0.0);
    
    r3d_calcArcPoints(arc, n, c, axis, r);
    r3d_calcArcPoints(arc_inner, n, c, axis, r-t);
    for (i = 0; i < (n%2 ? (n-1)/2 : n/2); i++)
    {
	r3d_quad(fp, arc[i], arc[(i+1)%n],
	    arc_inner[(i+1)%n], arc_inner[i], color);
    }
    for (i = 0; i < n; i++) freeptPtr(arc_inner[i]);
}

void r3d_calcArcPoints (ptPtr arc[], int n, ptPtr c, ptPtr axis, double r)
{
    int i = 0, j = 0;
    double wedge_rad = 0.0, axmag = 0.0, naz, nax;
    double x = 0.0, y = 0.0, z = 0.0, sina = 0.0, cosa = 0.0, theta = 0.0;
    double ang_x = 0.0, ang_y = 0.0, ang_z = 0.0;
    if (!c || !axis) return;
    if (!n || !r) return;
    if (n > MAX_ARC) return;
    
    axmag = sqrt(axis->x*axis->x+axis->y*axis->y+axis->z*axis->z);
     
    
    /* draw a circle of radius r centered at origin with z-hat as the
     * axis vector, then rotate the disc until the axis vector aligns
     * with the axis argument. */
    theta = 0.0;
    wedge_rad = M_PI*2/n;
    for (i = 0; i < n; i++)
    {
	arc[i]->x = r*cos(theta);
	arc[i]->y = r*sin(theta);
	arc[i]->z = 0.0;
	/* draws radius-r arc centered at origin */
	theta += wedge_rad;
    }
    
    /* To align the disc along the requested axis, at most two 
     * rotations are required.  The first can either be around the
     * x or y axis (in the current plane of the disc).
     * The second is around the z-axis.  */
    ang_y = 0.0; 
    ang_x = 0.0;
    if (axis->x || axis->y) 
    {
	if (!(axis->y)) /* need rotate along y axis only */
	{
	    ang_y = acos(axis->z/axmag);
	    cosa = cos(ang_y);
	    sina = sin(ang_y);
	    for (i = 0; i < n; i++)
	    {
		x = arc[i]->x*cosa;
	        z = -arc[i]->x*sina;
		arc[i]->x = x;
		arc[i]->z = z;
	    }
	}
	else if (!(axis->x)) /* need rotate along x axis only */
	{
	    ang_x = acos(axis->z/axmag);
	    cosa = cos(ang_x);
	    sina = sin(ang_x);
	    for (i = 0; i < n; i++)
	    {
		y = arc[i]->y*cosa;
	        z = arc[i]->y*sina;
		arc[i]->y = y;
		arc[i]->z = z;
	    }
	}
	else /* two-step rotation required */
	{
	    ang_y = acos(axis->z/axmag);
	    cosa = cos(ang_y);
	    sina = sin(ang_y);
	    for (i = 0; i < n; i++)
	    {
		x = arc[i]->x*cosa;
	        z = -arc[i]->x*sina;
		arc[i]->x = x;
		arc[i]->z = z;
	    }
	    ang_z = acos(axis->x/(sqrt(axis->x*axis->x+axis->y*axis->y)));
	    if (axis->y < 0.0) ang_z *= -1;
	    cosa = cos(ang_z);
	    sina = sin(ang_z);
	    for (i = 0; i < n; i++)
	    {
		x = arc[i]->x*cosa - arc[i]->y*sina;
	        y = arc[i]->x*sina + arc[i]->y*cosa;
		arc[i]->x = x;
		arc[i]->y = y;
	    }
	}
    }
    
    /* Now, translate the center of the disc */
    for (i = 0; i < n; i++)
    {
 	arc[i]->x += c->x;
	arc[i]->y += c->y;
	arc[i]->z += c->z;
    }
    
}

void r3d_arc (FILE * fp, ptPtr p0, ptPtr p1, ptPtr p2,
	     int n, double liner, color_t color, ptPtr arc[])
{
    int i = 0, j = 0;
    double wedge_rad = 0.0, axmag = 0.0, naz, nax;
    double x = 0.0, y = 0.0, z = 0.0, sina = 0.0, cosa = 0.0, theta = 0.0;
    double ang_x = 0.0, ang_y = 0.0, ang_z = 0.0, r = 0.0;
    ptPtr cross = ptPtr_Init(0.0, 0.0, 0.0);
    ptPtr mp01 = ptPtr_Init(0.0, 0.0, 0.0);
    ptPtr mp12 = ptPtr_Init(0.0, 0.0, 0.0);
    ptPtr disp01 = ptPtr_Init(0.0, 0.0, 0.0);
    ptPtr disp12 = ptPtr_Init(0.0, 0.0, 0.0);
    ptPtr disp02 = ptPtr_Init(0.0, 0.0, 0.0);
    if (!fp) return;
    if (!p0 || !p1 || !p2) return;
    if (!n) return;
    
    mp01->x = (p0->x + p1->x)/2.0;
    mp01->y = (p0->y + p1->y)/2.0;
    mp01->z = (p0->z + p1->z)/2.0;
    mp12->x = (p1->x + p2->x)/2.0;
    mp12->y = (p1->y + p2->y)/2.0;
    mp12->z = (p1->z + p2->z)/2.0;

    disp01->x = p1->x - p0->x;
    disp01->y = p1->y - p0->y;
    disp01->z = p1->z - p0->z;
    disp12->x = p2->x - p1->x;
    disp12->y = p2->y - p1->y;
    disp12->z = p2->z - p1->z;
    disp02->x = p2->x - p0->x;
    disp02->y = p2->y - p0->y;
    disp02->z = p2->z - p0->z;

    cross->x = (disp01->y*disp12->z - disp12->y*disp12->z);
    cross->y = (disp01->z*disp12->x - disp12->z*disp12->x);
    cross->z = (disp01->x*disp12->y - disp12->x*disp12->y);

    
    /* Now, each point on the arc of this disc belongs to exactly
     * one chord.  So, draw each chord as a cylinder. */
    for (i = 0; i < n-1; i++)
    {
	r3d_cyl(fp, arc[i], arc[(i+1)], liner, color);
    }
}


void r3d_cone (FILE * fp, ptPtr b, ptPtr t, double r, int n, color_t color)
{
    ptPtr arc[MAX_ARC];
    int i;
    ptPtr axis = ptPtr_Init(0.0, 0.0, 0.0);
    if (!fp) return;
    if (!t || !b || !r) return;
    
    axis->x = t->x - b->x;
    axis->y = t->y - b->y;
    axis->z = t->z - b->z;
    
    for (i = 0; i < MAX_ARC; i++) arc[i] = ptPtr_Init(0.0, 0.0, 0.0);
    
    fprintf(fp, "#cone -- n = %i, MAX_ARC = %i\n", n, MAX_ARC);
    fprintf(fp, "#cone -- draw base:\n");
    fflush(fp);
    r3d_disc(fp, b, axis, r, n, color, arc);
    fprintf(fp, "#cone -- draw sides:\n");
    fflush(fp);
    for (i = 0; i < n; i++)
    {
	r3d_tri(fp, arc[i], t, arc[(i+1)%n], color);
    }
    fprintf(fp, "#cone -- done\n");
    fflush(fp);
    for (i = 0; i < MAX_ARC; i++) freeptPtr(arc[i]);
}

void r3d_3darrow (FILE * fp, ptPtr p0, ptPtr p1, double r, 
		  double ahl, double ahw, int n, color_t color)
{
    ptPtr ahb = ptPtr_Init(0.0, 0.0, 0.0);
    ptPtr disp = NULL;
    double tl = 0.0, fac = 0.0;
    
    if (!fp || !p0 || !p1) return;
    
    disp = ptPtr_Init(p1->x-p0->x, p1->y-p0->y, p1->z-p0->z);
    tl = sqrt(disp->x*disp->x + disp->y*disp->y + disp->z*disp->z);
    fac = (tl - ahl)/tl;
    if (fac < 0.0) 
    {
	fprintf(fp, "#error -- requested arrowhead length(%.5le)\n", ahl);
	fprintf(fp, "#greater than length of arrow(%.5le).\n", tl);
	return;
    }
    ahb->x = p0->x + fac*(p1->x-p0->x);
    ahb->y = p0->y + fac*(p1->y-p0->y);
    ahb->z = p0->z + fac*(p1->z-p0->z);
    
    fprintf(fp, "#arrow -- draw shaft\n");
    fflush(stdout);
    r3d_cyl(fp, p0, ahb, r, color);
    fprintf(fp, "#arrow -- draw head\n");
    fflush(stdout);
    r3d_cone(fp, ahb, p1, ahw/2.0, 
	(n > MAX_ARC_ARROW_CONE ? MAX_ARC_ARROW_CONE : n), color);
    fprintf(fp, "#arrow -- done\n");
    fflush(stdout);
    
}

void r3d_volumeElement_cubic (FILE * fp, ptPtr c, double d, color_t color, double clarity)
{
    ptPtr v[8];
    int i = 0;
    
    if (!c || !fp) return;
    
    if (transparency_ && !clarity) r3d_transparency_off(fp);
    if (clarity > 0.0) r3d_transparency_on(fp, clarity);
    for (i = 0; i < 8; i++) v[i] = ptPtr_Init(0.0, 0.0, 0.0);

/****** UNDER CONSTRUCTION 9/2/98 ***************************/
/****** UNDER CONSTRUCTION 9/2/98 ***************************/
/****** UNDER CONSTRUCTION 9/2/98 ***************************/
    v[0]->x = c->x + d/2.0;
    v[0]->y = c->y + d/2.0;
    v[0]->z = c->z + d/2.0;
    v[1]->x = c->x - d/2.0;
    v[1]->y = v[0]->y;
    v[1]->z = v[0]->z;
    v[2]->x = v[1]->x;
    v[2]->y = c->y - d/2.0;
    v[2]->z = v[0]->z;
    v[3]->x = v[0]->x;
    v[3]->y = v[2]->y;
    v[3]->z = v[0]->z;
    v[4]->x = v[0]->x;
    v[4]->y = v[0]->y;
    v[4]->z = c->z - d/2.0;
    v[5]->x = v[1]->x;
    v[5]->y = v[1]->y;
    v[5]->z = v[4]->z;
    v[6]->x = v[2]->x;
    v[6]->y = v[2]->y;
    v[6]->z = v[4]->z;
    v[7]->x = v[3]->x;
    v[7]->y = v[3]->y;
    v[7]->z = v[4]->z;
    
    r3d_quad(fp, v[0], v[1], v[2], v[3], color);  /* +z (top) */
    r3d_quad(fp, v[3], v[7], v[4], v[0], color);  /* +x (right) */
//    r3d_quad(fp, v[0], v[4], v[5], v[1], color);  /* +y (back) */
//    r3d_quad(fp, v[1], v[5], v[6], v[2], color);  /* -x (left) */
    r3d_quad(fp, v[2], v[6], v[7], v[3], color);  /* -y (front) */
//    r3d_quad(fp, v[4], v[5], v[6], v[7], color);  /* -z (bottom) */

    for (i = 0; i < 8; i++) freeptPtr(v[i]);
}

void r3d_volumeElement_polar (FILE * fp, ptPtr c, 
			      double dr, double dt, double dp,
			      color_t color, double clarity)
{
    ptPtr v[8];
    int i = 0;
    
    double r_p, r_m, t_p, t_m, p_p, p_m;
    if (!c || !fp) return;
     
    if (transparency_ && !clarity) r3d_transparency_off(fp);
    if (clarity > 0.0) r3d_transparency_on(fp, clarity);
    for (i = 0; i < 8; i++) v[i] = ptPtr_Init(0.0, 0.0, 0.0);

    r_p = c->x + 0.5*dr;
    r_m = c->x - 0.5*dr;
    t_p = c->y + 0.5*dt;
    t_m = c->y - 0.5*dt;
    p_p = c->z + 0.5*dp;
    p_m = c->z - 0.5*dp;

    /* Angles are given in degrees, so convert to radians here */
    t_p *= M_PI/180.0;
    t_m *= M_PI/180.0;
    p_p *= M_PI/180.0;
    p_m *= M_PI/180.0;
    r_p *= 0.5;
    r_m *= 0.5;  /* Scaling for the viewing window */
    
    /* 0 --> -+- */
    v[0]->x = r_m*sin(t_p)*cos(p_m);
    v[0]->y = r_m*sin(t_p)*sin(p_m);
    v[0]->z = r_m*cos(t_p);
    
    /* 1 --> --- */
    v[1]->x = r_m*sin(t_m)*cos(p_m);
    v[1]->y = r_m*sin(t_m)*sin(p_m);
    v[1]->z = r_m*cos(t_m);
    
    /* 2 --> --+ */
    v[2]->x = r_m*sin(t_m)*cos(p_p);
    v[2]->y = r_m*sin(t_m)*sin(p_p);
    v[2]->z = r_m*cos(t_m);
    
    /* 3 --> -++ */
    v[3]->x = r_m*sin(t_p)*cos(p_p);
    v[3]->y = r_m*sin(t_p)*sin(p_p);
    v[3]->z = r_m*cos(t_p);
    
    /* 4 --> ++- */
    v[4]->x = r_p*sin(t_p)*cos(p_m);
    v[4]->y = r_p*sin(t_p)*sin(p_m);
    v[4]->z = r_p*cos(t_p);

    /* 5 --> +-- */
    v[5]->x = r_p*sin(t_m)*cos(p_m);
    v[5]->y = r_p*sin(t_m)*sin(p_m);
    v[5]->z = r_p*cos(t_m);
    
    /* 6 --> +-+ */
    v[6]->x = r_p*sin(t_m)*cos(p_p);
    v[6]->y = r_p*sin(t_m)*sin(p_p);
    v[6]->z = r_p*cos(t_m);
    
    /* 7 --> +++ */
    v[7]->x = r_p*sin(t_p)*cos(p_p);
    v[7]->y = r_p*sin(t_p)*sin(p_p);
    v[7]->z = r_p*cos(t_p);
    
    
/* ****** UNDER CONSTRUCTION 9/30/98 ********** */
/* ****** UNDER CONSTRUCTION 9/30/98 ********** */
/* ****** UNDER CONSTRUCTION 9/30/98 ********** */
/* ****** UNDER CONSTRUCTION 9/30/98 ********** */
 
//    r3d_quad(fp, v[0], v[1], v[2], v[3], color);  /* -r (back) */
//    r3d_quad(fp, v[3], v[7], v[4], v[0], color);  /* +t (bottom) */
    r3d_quad(fp, v[0], v[4], v[5], v[1], color);  /* +p (right) */
    r3d_quad(fp, v[1], v[5], v[6], v[2], color);  /* -t (top) */
    r3d_quad(fp, v[2], v[6], v[7], v[3], color);  /* -p (left) */
    r3d_quad(fp, v[4], v[5], v[6], v[7], color);  /* +r (front) */
   
    for (i = 0; i < 8; i++) freeptPtr(v[i]);
}

void r3d_transparency_on (FILE * fp, double clarity)
{
    if (clarity > 1.0 || clarity < 0.0)
    {
	fprintf(fp, "# Warning!  Invalid clarity requested.  Setting to 0.8.\n");
	clarity = 0.8;
    }
    if (!transparency_)
    {
	fprintf(fp, "# Transparency (Clarity = %.2lf) ON.\n", clarity);
	fprintf(fp, "8\n");
	fprintf(fp, " 17.  0.6       -1.0 -1.0 -1.0     %0.2lf   0 0 0 0\n", 
	    clarity);
	transparency_ = 1;
    }
}

void r3d_transparency_off (FILE * fp)
{
    if (transparency_)
    {
	fprintf(fp, "# Transparency off.\n");
	fprintf(fp, "9\n");
	transparency_ = 0;
    }
}
