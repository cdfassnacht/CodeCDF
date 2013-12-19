#ifndef lc_setup_h
#define lc_setup_h

/*.......................................................................
 *
 * Enumeration for setup file parameters
 *
 */

enum {
  DOCHI,
  DOXCORR,
  DOACORR,
  DODISP,
  DODCF,
  DOCURVEFIT,
  DISPCHOICE,
  D2DELTA,
  OUTFILE,
  ACHIFILE,
  CCHIFILE,
  DCHIFILE,
  CHILOG,
  XCLOG,
  DOOVERLAP,
  MU0,
  NMU,
  TAU0,
  DTAU,
  NTAU,
  FLAGBAD,
  MEANCHOICE,
  DOSMOOTH,
  BOXCAR,
  MEDIAN,
  TRIANGLE,
  VARBOX,
  VARTRI,
  GAUSS,
  NINTERP,
  INTSTEP,
  INTSTART,
  ASKSTART,
  ROOT,
  DEFAULT,
  SETUPERR
};

/*.......................................................................
 *
 * Enumeration for type of smoothing or interpolation (setup->dosmooth)
 *
 */

enum {
  SMUNSET = -2, /* Value is not set */
  NOSMOOTH,     /* No smoothing or interpolation performed */
  INTONLY,      /* Piecewise linear interpolation */
  SMONLY,       /* Smoothing and interpolating by using a running mean */
  SMINPLACE,    /* Replace each point by its smoothed value - NO interp. */
};

/*.......................................................................
 *
 * Enumeration for smoothing functions (setup->smtype)
 *
 */

enum {
  SMBOXCAR,
  SMMEDIAN,
  SMTRIANGLE,
  SMGAUSS,
  SMVARBOX,
  SMVARTRI
};

/*.......................................................................
 *
 * Enumeration for dispersion analysis method (setup->dispchoice)
 *
 */

enum {
  D21,      /* Pelt et al. D^2_1 (non-parametric) method */
  D21M,     /* Pelt et al. D^2_1 method with all curves simultaneously */
  D22,      /* Pelt et al. D^2_2 method */
  DLOVELL   /* Lovell et al. modification of the D^2_2 method */
};

/*.......................................................................
 *
 * Enumeration for speed flag.
 *
 */

enum {
  SLOW,
  FAST
};

/*.......................................................................
 *
 * Verbosity enumeration
 */

enum {
  NONE,
  BESTFIT,
  GRID
};

/*.......................................................................
 *
 * Yes/no enumeration
 *
 */

enum {
  UNSET = -1,
  NO,
  YES
};

/*.......................................................................
 *
 * Structure definitions
 *
 */

typedef struct {
  int ncurves;          /* Number of curves being compared */
  int nfiles;           /* Number of input files */
  char *infile[4];      /* Names of the input light curve files */
  int npoints[4];       /* Number of points in each input light curve */
  char *setupfile;      /* Name of the optional setup file */
  char *outfile;        /* Name of output file (used for simulations) */
  int tauset;           /* Flag set to YES if tau0 was set in input file */
  int ntau;             /* Number of delay steps on either side of tau0 */
  double tau0[4];       /* Initial guess for delays */
  double dtau;          /* Stepsize to use in grid search for delays */
  int dooverlap;        /* Set to YES for using ratios from overlap region */
  int nmu;              /* Number of mag steps on either side of mu0 */
  double mu0[4];        /* Initial guesses for magnifications */
  double dmu;           /* Stepsize for magnification */
  int dochi;            /* Flag set to YES for chisq analysis */
  int doxcorr;          /* Flag set to YES for cross-corr analysis */
  int doacorr;          /* Flag set to YES for auto-corr analysis */
  int dodisp;           /* Flag set to YES for dispersion analysis */
  int dodcf;            /* Flag set to YES for DCF analysis */
  int docurvefit;       /* Flag set to YES for curve fitting analysis */
  int dispchoice;       /* Choice of dispersion analysis method */
  float d2delta;        /* delta parameter for D^2_2 dispersion method */
  char achifile[MAXC];  /* File for B-A chisq minimization output */
  char cchifile[MAXC];  /* File for B-C chisq minimization output */
  char dchifile[MAXC];  /* File for B-D chisq minimization output */
  char chilog[MAXC];    /* Name of chisq log file */
  char xclog[MAXC];     /* Name of cross-correlation log file */
  int dosmooth;         /* Type of smoothing or interp. (see enum above) */
  int smtype;           /* Type of smoothing function (see enum above) */
  float smwidth;        /* Width of smoothing window in days */
  int ninterp;          /* Number of points in interpolated curve */
  float intstep;        /* Step size for added interpolation */
  float intstart;       /* Starting day for interpolation */
  int askstart;         /* Ask user about starting day for interpolation? */
  int nvar;             /* Number of points in variable boxcar box */
  int flagbad;          /* Flag set to 1 to flag bad days */
  int meanchoice;       /* Method of normalizing secondary flux cals */
  char root[MAXC];      /* Root name for files */
} Setup;

/*.......................................................................
 *
 * Function declarations
 *
 */

Setup *new_setup(int size);
Setup *del_setup(Setup *setup);
Setup *setup_from_command_line(char *argv[], int narg);
int get_setup_params(Setup *setup, Fluxrec **lc);
int setup_file(Setup *setup, char *inname);
int read_setup_line(char *line, char *keyword);
int setup_interp(Setup *setup);
int setup_delays(Setup *setup);
void setup_interp_summary(Setup *setup);
void setup_delays_summary(Setup *setup);
void get_meanchoice(Setup *setup);
void get_interp_choice(Setup *setup);
void get_interp_step(Setup *setup);
void get_smooth_fn(Setup *setup);
void get_smooth_width(Setup *setup);
void get_nvar(Setup *setup);
void set_grid_params(Setup *setup, Fluxrec *raw, int nraw);

#endif
