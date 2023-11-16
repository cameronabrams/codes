/* camvox.c 
 * 
 * (c) 1997-2000 cameron abrams
 * 
 * start_doc
 *
 * Name
 *    camvox
 *
 * Usage
 *    camvox -f datafile [data options] \
 *				  [Raster3d options] > output.r3d
 *                                            **or** | render [render options]
 *
 * Description
 *    camvox is used to render a volume histogram (VH).  A volume histogram
 *    is a set of bins, called `voxels', each of which is indexed by three 
 *    independent coordinates, such as (x,y,z) or (r,theta,phi).  Each voxel
 *    has associated with it a single `value'.  So a volume histogram can
 *    be thought of as nothing more than a discretized function of 
 *    three independent variables.  camvox reads in volume histogram
 *    data from an ASCII text input file of the following form:
 *
 *    #HEADERS
 *    x_1 y_1 z_1 f_1
 *    ...
 *    x_i y_i z_i f_i
 *
 *    Headers include keyword/value pairs.  The following keyword
 *    value pairs are recognized in an input data file:
 *    #POLAR  (no value); specifies that the VH is in polar coordinates
 *                        (default is cartesian)
 *    #[XYZ]_BS value; specifies [XYZ] binsize
 *    #[RTP]_BS value; specifies [r,theta,phi] binsize
 *
 *    Lines that do not begin with `#' must contain data.
 *    Each line of date contains four numbers:  the
 *    first are the three independent coordinates indexing a single
 *    voxel, and f_i is the voxel's value.  Voxel values are assumed
 *    to be in the range [0:1].
 *    Voxel values are used to set the voxel's color in the final image.
 * 
 *    camvox's primary purpose is to translate this data file into
 *    a stream of commands that is interpreted by the Raster3D program
 *    `render', which actually constructs the image of the histogram.
 *    This stream can be redirected into a file, or piped directly
 *    to render itself (the preferred mode of operation for me anyway).
 *
 *    camvox can operate in two coordinate systems:
 *      1. Cartesian: each voxel is indexed by (x,y,z) and voxels
 *         are cubic.  Each of (x,y,z) must lie in the domain [-0.5:0.5].
 *      2. Spherical: each voxel is indexed by (r,theta,phi), and
 *         voxels are discretized as dr, r*sin(theta)*d(theta), r*d(phi).
 *         The domain of r is [0:1], theta [0:180], phi [0:360].
 *
 *     
 *
 * Command-line Options (dreadfully incomplete)
 *    -f <name of data file>
 *    -voxd <#>; dimension of cubic voxel; discretization of cartesian VH
 *               May appear in the data file as a line of the format:
 *               #X_BS <binsize>
 *    -voxdr <#>; r-dimension of spherical coordinate voxel
 *               May appear in the data file as a line of the format:
 *               #R_BS <binsize>
 *    -voxdt <#>; theta-dimension of spherical coordinate voxel
 *               May appear in the data file as a line of the format:
 *               #T_BS <binsize>
 *    -voxdp <#>; phi-dimension of spherical coordinate voxel
 *               May appear in the data file as a line of the format:
 *               #P_BS <binsize>
 *    -voxmin <#>; minimum voxel value; voxels with values below this
 *                  are not rendered.
 *
 * Imported
 *    point
 *    dblmat
 *    colormap
 *    r3d_utils
 *    cam_colors.ext
 *
 * Package
 *    md srs 1
 *
 * Module
 *    camvox
 *
 * Author
 *    Cameron Abrams
 *
 * Copyright
 *    Copyright 1997-2000 Cameron Abrams and the Regents of the University
 *    of California, Berkeley
 *    
 * end_doc
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "dblmat.h"
#include "point.h"
#include "colormap.h"
#include "r3d_utils.h"
#include "cam_colors.ext"

extern char * rot_label[];
extern double ImageBufferSize_;
extern int xTiles_, yTiles_;
extern int xPixPerTile_, yPixPerTile_;
extern int Scheme_;
extern color_t * bgClr;
extern int Shadows_;
extern int Phong_;
extern double SecLtContr_;
extern double AmbLtContr_;
extern double SpecRefComp_;
extern double Eyepos_;
extern ptPtr MainLtPos;

void usage (FILE * fp)
{
    if (!fp) return;
    fprintf(fp, "Usage: camvox [control params] [> <r3dName>]\n");
    fprintf(fp, "Raster3D render control parameters:\n");
    fprintf(fp, "\t [-xt <xTiles(32)>] [-yt <yTiles(32)>]\n");
    fprintf(fp, "\t [-xp <xPixelsPerTile(18)>] [-yp <yPixelsPerTile(18)>]\n");
    fprintf(fp, "\t [-scheme <scheme(3)>]\n");
    fprintf(fp, "\t [-bg <r> <g> <b> (0 0 0)]\n");
    fprintf(fp, "\t [-noshadows] [-phong <phongPower(25)>]\n");
    fprintf(fp, "\t [-slc <secondaryLightContributionFraction(0.15)>]\n");
    fprintf(fp, "\t [-alc <ambientLightContributionFraction(0.05>]\n");
    fprintf(fp, "\t [-src <specularReflectionComponentFraction(0.25)>]\n");
    fprintf(fp, "\t [-eyepos <EYEPOS(4.0)>]\n");
    fprintf(fp, "\t [-mlpos <x> <y> <z> (1 1 1)] (Main Light Position)\n");
}

enum { MALCOLM_X, ARCHIE_BUNKER };
char scr_ln[255];
void main (int argc, char * argv[])
{
    FILE * fp = NULL;
    char * d_h = "camvox::main", cmt = '#';
    char * data3file = "data.3dist";
    char keyword[255];
    double keyval = 0.0;
    int i = 0,  j = 0;
    double scale = 1.0, span = 1.0;
    ptPtr shift = ptPtr_Init(0.0, 0.0, 0.0);
    ptPtr pt = ptPtr_Init(0.0, 0.0, 0.0);
    ptPtr tpt[4] = {NULL, NULL, NULL, NULL};
    ptPtr zpax = ptPtr_Init(0.0, 0.0, 0.5);
    ptPtr zmax = ptPtr_Init(0.0, 0.0, -0.5);
    ptPtr xpax = ptPtr_Init(0.5, 0.0, 0.0);
    ptPtr xmax = ptPtr_Init(-0.5, 0.0, 0.0);
    ptPtr ypax = ptPtr_Init(0.0, 0.5, 0.0);
    ptPtr ymax = ptPtr_Init(0.0, -0.5, 0.0);
    color_t tc = {0.0, 0.0, 0.0};
    int rot_order[MAXROTS];
    double rots[MAXROTS];
    ptPtr origin = ptPtr_Init(0.0, 0.0, 0.0);
    int arcsegs = 20, nmin, nmax, nsegs;
    double liner = 0.0025, axisr = 0.005;
    double x = 0.0, y = 0.0, z = 0.0, c = 0, rmin, rmax;
    double voxd = 0.01, vox_minv = 0.0, clarity = 0.0;
    double voxdr = 0.04, voxdt = 2.0,  voxdp = 5.0;
    double theta_i = -1.0, dTheta = 0.0, tickTheta = 0.0;
    short greyScaleDirection = MALCOLM_X;
    short nThetaTicks = 0;
    ptPtr arc[MAX_ARC];
    short twin = 0;
    short use_clarity = 0, inColor = 1, spec_sphere = 0, just_draw_legend = 0;
    short nodata = 0, polar_ = 0;
    double tickLen = 0.03, tickWidDeg = 2.0;
    double dimFactor = 0.75;
    color_t * axClr = &grey2;
    
    /* Initializations */
    Shadows_=0;
    for (i=0; i<MAXROTS; i++) rots[i]=0.0;
    for (i=0; i<MAXROTS; i++) rot_order[i]=0;
    for (i=0; i<MAX_ARC; i++) arc[i]=ptPtr_Init(0.0, 0.0, 0.0);
    for (i=0; i<4; i++) tpt[i]=ptPtr_Init(0.0, 0.0, 0.0);

    /* Handle command-line arguments (old way) */
    for (i=1; i<argc; i++)
    {
	if (!strcmp(argv[i], "-xr"))
	{
	    if (j==MAXROTS)
	    {
		printf("Error: too many rotations requested.\n");
	    }
	    else
	    {
		rots[j]=atof(argv[++i]);
		rot_order[j]=X;
		j++;
	    }
	}
	else if (!strcmp(argv[i], "-yr"))
	{
	    if (j==MAXROTS)
	    {
		printf("Error: too many rotations requested.\n");
	    }
	    else
	    {
		rots[j]=atof(argv[++i]);
		rot_order[j]=Y;
		j++;
	    }
	}
	else if (!strcmp(argv[i], "-zr"))
	{
	    if (j==MAXROTS)
	    {
		printf("Error: too many rotations requested.\n");
	    }
	    else
	    {
		rots[j]=atof(argv[++i]);
		rot_order[j]=Z;
		j++;
	    }
	}
	else if (!strcmp(argv[i], "-f")) data3file = argv[++i];
	else if (!strcmp(argv[i], "-voxd")) voxd = atof(argv[++i]);
	else if (!strcmp(argv[i], "-voxdr")) voxdr = atof(argv[++i]);
	else if (!strcmp(argv[i], "-voxdt")) voxdt = atof(argv[++i]);
	else if (!strcmp(argv[i], "-voxdp")) voxdp = atof(argv[++i]);
	else if (!strcmp(argv[i], "-voxc")) use_clarity = 1;
	else if (!strcmp(argv[i], "-voxmin")) vox_minv = atof(argv[++i]);
	else if (!strcmp(argv[i], "-spec_sphere")) spec_sphere = 1;
	else if (!strcmp(argv[i], "-legend_only")) 
	    just_draw_legend = (short)atoi(argv[++i]);
	else if (!strcmp(argv[i], "-white_hi")) 
	    greyScaleDirection = ARCHIE_BUNKER;
 	else if (!strcmp(argv[i], "-arcsegs")) arcsegs = atoi(argv[++i]);
	else if (!strcmp(argv[i], "-liner")) liner = atof(argv[++i]);
	else if (!strcmp(argv[i], "-axisr")) axisr = atof(argv[++i]);
	else if (!strcmp(argv[i], "-nocolor")) inColor = 0;
	else if (!strcmp(argv[i], "-dimfac")) dimFactor = atof(argv[++i]);
	else if (!strcmp(argv[i], "-twin")) twin = 1;
	else if (!strcmp(argv[i], "-nodata")) nodata = 1;
	else if (!strcmp(argv[i], "-theta_i")) theta_i = atof(argv[++i]);
	else if (!strcmp(argv[i], "-nthetaticks")) 
	    nThetaTicks=(short)atoi(argv[++i]);
	else if (!strcmp(argv[i], "-ticklen")) tickLen = atof(argv[++i]);
	else if (!strcmp(argv[i], "-tickwiddeg")) tickWidDeg = atof(argv[++i]);
	else if (!strcmp(argv[i], "-xt")) xTiles_ = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-yt")) yTiles_ = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-xp")) xPixPerTile_ = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-yp")) yPixPerTile_ = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-scheme")) Scheme_ = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-phong")) Phong_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-slc")) SecLtContr_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-alc")) AmbLtContr_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-src")) SpecRefComp_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-eyepos")) Eyepos_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-ibs")) ImageBufferSize_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-zoom")) scale = atof(argv[++i]);
	else if (!strcmp(argv[i], "-bg"))
	{
	    if (isdigit(argv[++i][0]))
	    {
		bgClr = color_init(bgClr, 
			atof(argv[i]), atof(argv[++i]), atof(argv[++i]));
	    }
	    else bgClr = colorPtr(argv[i]);
	}
 	else if (!strcmp(argv[i], "-mlpos"))
	{
	    x = atof(argv[++i]);
	    y = atof(argv[++i]);
	    z = atof(argv[++i]);
	    MainLtPos = ptPtr_Init(x, y, z);
	}
 	else if (!strcmp(argv[i], "-trans"))
	{
	    shift->x = atof(argv[++i]);
	    shift->y = atof(argv[++i]);
	    shift->z = atof(argv[++i]);
	}
 	else if (!strcmp(argv[i], "-shadows")) 
	{
	    Shadows_ = 1;
	}
	else
	{
	    if (strcmp(argv[i], "-h"))
		printf("Error:  Keyword %s not recognized.\n", argv[i]);
	    usage(stdout);
	    exit(-1);
	}
    }
    
    /* First, output the Raster3D header, which contains 
     * rendering instructions, scaling and rotation information, etc. */
    r3d_header(stdout, rots, rot_order, scale, shift, span);
    
    /* Now begin the description of the scene. */
    pt->x = 0.0;
    pt->y = 0.0;
    pt->z = 0.0;
    if (!inColor) axClr = &grey4;
    if (!just_draw_legend)
    {
      /* first draw the coordinate axes -- drawn as spherical */
	rmin = 0.1;
	rmax = 0.5;
	nmin = arcsegs;
	nmax = MAX_ARC;
	for (x = rmin; x <= rmax; x += 0.1)
	{
	    nsegs = arcsegs + (int)((x-rmin)/(rmax-rmin)*(nmax-nmin));
	    if (nsegs%2) nsegs++;
	    if (!twin)
		r3d_halfring(stdout, origin, zmax, x, 2*liner, nsegs, *axClr, arc);
	    else
		r3d_ring(stdout, origin, zpax, x, 2*liner, nsegs, *axClr, arc);
	    r3d_halfring(stdout, origin, ypax, x, 2*liner, nsegs, *axClr, arc);
	}
	/* theta-ticks on the outer arc:
	 *  - divide the arc into nThetaTicks+1 segments
	 */
        if (nThetaTicks)
	{
	    dTheta = M_PI/2.0/(nThetaTicks+1);
	    for (i = 0; i < nThetaTicks+1; i++)
	    {
		tpt[0]->x = (0.5-tickLen)*sin(tickTheta)*cos(0.1/180*M_PI);
		tpt[0]->y = (0.5-tickLen)*sin(tickTheta)*sin(-0.1/180*M_PI);
		tpt[0]->z = (0.5-tickLen)*cos(tickTheta);
		tpt[1]->x = (0.5)*sin(tickTheta-tickWidDeg/180.0*M_PI/2.0)*cos(0.1/180*M_PI);
		tpt[1]->y = (0.5)*sin(tickTheta-tickWidDeg/180.0*M_PI/2.0)*sin(-0.1/180*M_PI);
		tpt[1]->z = (0.5)*cos(tickTheta-tickWidDeg/180.0*M_PI/2.0);
		tpt[2]->x = (0.5)*sin(tickTheta+tickWidDeg/180.0*M_PI/2.0)*cos(0.1/180*M_PI);
		tpt[2]->y = (0.5)*sin(tickTheta+tickWidDeg/180.0*M_PI/2.0)*sin(-0.1/180*M_PI);
		tpt[2]->z = (0.5)*cos(tickTheta+tickWidDeg/180.0*M_PI/2.0);
		tickTheta += dTheta;
		r3d_tri(stdout, tpt[0], tpt[1], tpt[2], *axClr);
	    }
	    
	}
	tpt[0]->z = tpt[1]->z = tpt[2]->z = 0.0;
	r3d_cyl(stdout, origin, zpax, axisr, *axClr);
	tpt[0]->x = xpax->x;
	tpt[0]->y = xpax->y-axisr;
	tpt[1]->x = xpax->x;
	tpt[1]->y = xpax->y+axisr;
	tpt[2]->x = xmax->x;
	tpt[2]->y = xmax->y-axisr;
	tpt[3]->x = xmax->x;
	tpt[3]->y = xmax->y+axisr;
	r3d_quad(stdout, tpt[0], tpt[1], tpt[3], tpt[2], *axClr);
	tpt[0]->x = ypax->x-axisr;
	tpt[0]->y = ypax->y;
	tpt[1]->x = ypax->x+axisr;
	tpt[1]->y = ypax->y;
	tpt[2]->x = (twin ? ymax->x-axisr : -axisr);
	tpt[2]->y = (twin ? ymax->y : 0);
	tpt[3]->x = (twin ? ymax->x+axisr : axisr);
	tpt[3]->y = (twin ? ymax->y : 0);
	r3d_quad(stdout, tpt[0], tpt[1], tpt[3], tpt[2], *axClr);
	r3d_3darrow(stdout, zpax, ptPtr_Init(0.0, 0.0, zpax->z+0.1), axisr, 
	    0.07, 0.03, arcsegs, *axClr);
	if (theta_i != -1.0)
	{
	    r3d_3darrow(stdout, 
		ptPtr_Init(-0.5*sin(M_PI/180.0*theta_i), 
		0.0, 0.5*cos(M_PI/180.0*theta_i)), origin, axisr, 
		0.07, 0.03, arcsegs, (inColor ? yellow : grey3));
	    if (spec_sphere)
		r3d_sphere(stdout, 
		    ptPtr_Init(0.5*sin(M_PI/180.0*theta_i), 
			0.0, 0.5*cos(M_PI/180.0*theta_i)), 
			1.5*axisr, yellow);
	}

	/* Now read the data file, and for each line, draw a voxel */
	fp = fopen(data3file, "r");
	if (fp && !nodata)
	{
	    while (fgets(scr_ln, 255, fp))
	    {
	      /* Line is data if it does NOT begin with `#' */
		if (scr_ln[0] != '#') 
		{
		    sscanf(scr_ln, "%lf %lf %lf %lf", &x, &y, &z, &c);	
		    pt->x = x;
		    pt->y = y;
		    pt->z = z;
		    tc.r = tc.g = tc.b = 
			dimFactor*(greyScaleDirection == ARCHIE_BUNKER ? 
			     c : (1.0 - c));
		    clarity = 0.6;
		    if (inColor)
		    {
			color_mapFloat(0.0, &blue, 1.0, &red, c, &tc, -1);
		    }
		    use_clarity = 0;
		    if (c < 0.25) use_clarity = 1;
		    if (c >= vox_minv)
		    {
			if (!polar_)
			{
			    r3d_volumeElement_cubic(stdout, pt, voxd, tc, 
				(use_clarity ? clarity : 0.0));
			    /* the twin cube has same x and z but -y */
			    if (twin)
			    {
				pt->y = -pt->y;
				r3d_volumeElement_cubic(stdout, pt, voxd, tc, 
				    (use_clarity ? clarity : 0.0));
			    }
			}
			else
			{
			    /* volume element is polar, so 
			     *	    pt->x is "r", pt->y is "theta", 
			     *	    pt->z is "phi"; angles in degrees.
			     */
			    r3d_volumeElement_polar(stdout, pt, 
				voxdr, voxdt, voxdp, tc, 
				(use_clarity ? clarity : 0.0));
			    if (twin)
			    {
				pt->z = -pt->z;
				r3d_volumeElement_polar(stdout, pt, 
				    voxdr, voxdt, voxdp, tc, 
				    (use_clarity ? clarity : 0.0));
			    }
			}
		    }
		}
		/* Line contains keyword/value pairs if it DOES begin
		   with `#' */
		else 
		{
		    sscanf(&(scr_ln[1]), "%s %lf", keyword, &keyval);
		    if (!(strcmp(keyword, "POLAR")))
		    {
			polar_ = 1;
		    }
		    else if (!strcmp(keyword, "X_BS"))
			voxd = keyval/2;
		    else if (!strcmp(keyword, "Y_BS"))
			voxd = keyval/2;
		    else if (!strcmp(keyword, "Z_BS"))
			voxd = keyval/2;
		    else if (!strcmp(keyword, "R_BS"))
			voxdr = keyval;
		    else if (!strcmp(keyword, "T_BS"))
			voxdt = keyval;
		    else if (!strcmp(keyword, "P_BS"))
			voxdp = keyval;
		}
	    }
	    fclose(fp);
	}
	else
	{
	    printf("# ERROR -- cannot read from %s.\n", data3file);
	    printf("# Program exits.\n");
	    exit(-1);
	}
    }
    else if (just_draw_legend)
    {
	printf("# Just drawing legend, size = %i\n", just_draw_legend);
	pt->x = 0.0;
	pt->z = 0.0;
	c = 1.0;
	voxd = 0.90/just_draw_legend;
	clarity = 0.6;
	for (pt->y = 0.45; pt->y > -0.45; pt->y -= voxd)
	{
	    color_mapFloat(0.0, &blue, 1.0, &red, c, &tc, -1);
	    use_clarity = 0;
	    if (c < 0.25) use_clarity = 1;
	    r3d_volumeElement_cubic(stdout, pt, voxd, tc, (use_clarity ? clarity : 0.0));
	    c = (0.45+pt->y-voxd)/0.90;
	}
    }
    printf("#Program Ends.\n");  
}
