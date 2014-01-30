#ifndef fitsim_h
#define fitsim_h

#include "cpgplot.h"
#include "structdef.h"
#include "coords.h"
#include "fitswcs.h"
#include "fitsfuncs.h"

/*.......................................................................
 *
 * Definitions
 *
 */

#define PIXDSS   0.02     /* arcmin/pix on DSS plates from skyview */
#define PIXNIRC  0.15     /* NIRC pixel scale (arcsec/pix) */
#define PIXCCD13 0.37166  /* CCD13 pixel scale (arcsec/pix) */
#define PIXLRIS  0.215    /* LRIS pixel scale (arcsec/pix) */
#define PIXWF    0.0996   /* WFPC2 WF pixel scale (arcsec/pix) */
#define PIXPC    0.04554  /* WFPC2 PC pixel scale (arcsec/pix) */
#define PIXNIC1  0.043    /* NICMOS NIC1 pixel scale (arcsec/pix) */
#define PIXDEF   0.05     /* Default pixel scale (arcsec/pix) */
#define MAXINST 7         /* Number of instruments in enumeration */

/*.......................................................................
 *
 * Enumeration for PGPLOT color
 */

enum {
  WHITE,
  BLACK,
  RED,
  GREEN,
  BLUE,
  CYAN,
  MAGENTA,
  YELLOW,
  ORANGE
};

/*.......................................................................
 *
 * Enumeration for instrument value
 *
 */

enum {
  QUIT = -1,
  OTHER, 
  DSS, 
  NIRC, 
  CCD13, 
  LRIS,
  RADIO,
  WF,
  PC,
  NIC1
};

/*.......................................................................
 *
 * Enumeration for setup parameter values
 *
 */

enum {
  UNSET = -1,
  FALSE,
  TRUE,
  AUTO
};

/*.......................................................................
 *
 * Enumeration for setup file parameters
 *
 */

enum {
  INST,
  SOURCE,
  TITLE,
  PIXSCALE,
  IMSIZE,
  IMCENT,
  AXISCENT,
  ZOOM,
  FRANGE,
  CIRCLE,
  CMUL,
  CONT,
  LABH,
  LABW,
  LAB,
  COMPASS,
  MARK,
  NOCONTOUR,
  NOLABEL,
  SLITMASK,
  BORDER,
  TRANS,
  DRAW,
  LWEIGHT,
  PLOTWFPC,
  DEFAULT,
  SETUPERR
};

/*.......................................................................
 *
 * Enumeration for axis label
 *
 */

enum {
  ARCSEC,
  ARCMIN,
  PIXEL
};

/*.......................................................................
 *
 * Structure definitions
 *
 */

typedef struct {
  int instrument;     /* Instrument (see enumeration above) */
  char source[80];    /* Source name */
  char title[80];     /* Title for plot */
  float cheight;      /* Character scale height for external labels */
  int pixset;         /* Flag set to 1 when pixel scale has been set */
  float pixscale;     /* Pixel scale for plot */
  int centset;        /* Flag set to 1 when central pixel has been set */
  float centx, centy; /* Central pixel */
  int sizeset;        /* Flag set to 1 when image size has been set */
  int xsize, ysize;   /* Size of image in pixels */
  float rxstart,rxend;  /* Start and end pixel for displayed window */
  float rystart,ryend;  /* Start and end pixel for displayed window */
  int axcentset;      /* Flag set to 1 when axis center has been set */
  float axcentx;      /* Axis center x position (in pixels) */
  float axcenty;      /* Axis center y position (in pixels) */
  float zoom[2];      /* Zoom factors for the two axes */
  float hi,lo;        /* High and low values to display */
  int docontour;      /* Flag set to FALSE for no contours */
  int ncont;          /* Number of contour levels */
  float *cont;        /* Contour levels */
  float cmul;         /* Multiplier of contour levels */
  int contcolor;      /* Color of contours */
  int dolabel;        /* Flag set to FALSE for no internal labels */
  int nintlab;        /* Number of internal labels */
  float labheight;    /* Text height for internal labels */
  int labwidth;       /* Text line width for internal labels */
  Labinfo *intlabs;   /* Internal labels */
  int ncircle;        /* Number of circles to plot */
  Labinfo *circle;    /* Circle positions, sizes, and labels */
  int drawcompass;    /* Flag set to 1 if a compass is requested */
  int drawmarker;     /* Flag set to 1 if a object marker is requested */
  float markx,marky;  /* Location of object marker */
  char slitfile[80];  /* Name of slitmask input file */
  Secat *slitmask;    /* Actual slitmask info */
  int nmask;          /* Number of members of slitmask structure */
  float rcirc;        /* Radius of circles, used in slitmask plots */
  int border;         /* Flag set if plot boundary and labels are requested */
  int axislab;        /* Flag for axis label scale (arsec or arcmin) */
  char xlab[80];      /* Label for x-axis */
  char ylab[80];      /* Label for y-axis */
  int trans;          /* Flag for linear, log or sqrt scaling */
  int ndraw;          /* Number of lines to draw */
  int lweight;        /* Width of line segments (must be between 1 and 201) */
  int ltype;          /* Line-drawing styles */
  Pos *drawstart;     /* Array of line starting positions */
  Pos *drawend;       /* Array of line ending positions */
  int plotwfpc;       /* Flag set if WFPC2 overlay is desired */
  Pos wfpccent;       /* Central position for WFPC overlay */
  float wfpcpa;       /* PA corresponding to PA_V3 for WFPC */
} Setup;

/*.......................................................................
 *
 * Function declarations
 *
 */

void fitsplt_copyright_notice();
int plot_fitsfile_init(char *fitsfile, char *setupfile, Image *image,
		       Setup *setup, int usefile, int autopars,
		       int screenflag);
int display_image(Image *image, Setup *setup);
int adjust_display(Image *image, Setup *setup);
int radplt(Image *image, Pos centpos, float rmax);
int calc_ibound(float rval, int imax);
void doborder(Setup *setup);
int draw_compass(Image *image, Setup *setup);
int compass_angle_keck(Image *image, int instrument, double *alpha);
int compass_angle_wfpc2(Image *image, double *alpha);
int compass_angle_nicmos(Image *image, double *alpha);
void draw_locator(Setup *setup);
int plot_circles(Setup *setup);
int plot_secat(Setup *setup, Secat *secat, int ncat, int format, int color,
	       float scale);
int plot_slitmask(Image *image, Setup *setup);
int draw_lines_interactive(Setup *setup);
int set_line_attrib(Setup *setup);
int draw_lines_setup(Setup *setup);
void draw_line(Pos begin, Pos end);
void plot_wfpc(Setup *setup);
int label_ra_dec(Image *image);
int print_image_data(Image *image);
int get_plotname(Image *image, Setup *setup, char *plottype, char *plotname);

/*
 * Functions in setup_fitsim.c
 */

Setup *new_setup(int size);
Setup *del_setup(Setup *setup);
void setup_help();
Setup *fill_setup(Image *image, char *setupfile, int usefile, int autopars);
int init_setup(Image *image, Setup *setup, int autopars);
int set_setup_defaults(Image *image, Setup *setup);
int setup_file(Image *image, Setup *setup, char *inname);
int read_setup_line(char *line, char *keyword);
int setup_menu(Image *image, Setup *setup);
int setup_interact(Image *image, Setup *setup);
int read_slitmask(Image *image, Setup *setup);

/*
 * Functions in get_params.c
 */

int get_instrument(Image *image, Setup *setup);
void read_instrument(Image *image, Setup *setup);
int get_object(Image *image, Setup *setup, char *filename);
void read_object(Image *image, Setup *setup);
int get_pixscale(Image *image, Setup *setup);
void read_dss_pixscale(Image *image, Setup *setup);
int read_file_pixscale(Image *image, float scale);
int get_imsize(Image *image, Setup *setup);
int get_zooms(Setup *setup);
int check_imsize(Image *image, Setup *setup, char *line, int interactive);
void get_title(Image *image, Setup *setup);
void toggle_compass(Setup *setup);
void get_marker_info(Setup *setup);
int get_frange(Image *image, Setup *setup);
int get_epoch(Image *image, float *epoch);
int get_transfunc(Setup *setup);
int get_intlab(Image *image, Setup *setup);
int get_labscale(Setup *setup);
int get_border_info(Image *image, Setup *setup);
int get_wfpc_overlay(Setup *setup);

/*
 * Functions in contours_fitsim.c
 */

int get_contours(Setup *setup);
int pcontour(float *cont, int ncont, float max);
int mcontour(float *cont, int ncont, float max);
int lcontour(float *cont, int ncont, float max);

/*
 * Slicing functions
 */

int imslice(Image *image);
Pos *make_vslice(Image *image, int x0, int ymin, int ymax, float *smin,
		 float *smax, double *ysmax);
int plot_vslice(Pos *slice, double ymin, double ymax, float smin,
		float smax, int dopoints);
void find_vslice_max(Image *image, int x, int ylo, int yhi, 
		     float *slmax, int *ymax);
int fit_continuum(Pos *slice, int npoints, double bklo_1, double bkhi_1, 
		  double bklo_2, double bkhi_2, float *m, float *b);
int pos_lsf(Pos *pos, int npoints, float *m, float *b);

#endif
