#ifndef lc_funcs_h
#define lc_funcs_h

#include "structdef.h"
#include "lc_setup.h"

#define NFLUXSTEP 50   /* Num. of steps on each side of central flux ratio */
#define FLUXSTEP 0.001 /* Size of steps in flux ratio */
#define NFIT 2         /* Number of varying parameters in the chisq fitting */
#define N08 4          /* Number of 1608 light curves */
#define ALAG 31.5      /* t_BA for fake light curves */
#define CLAG 36.5      /* t_BC for fake light curves */
#define DLAG 80.5      /* t_BD for fake light curves */
#define F34_0 823.859  /* Total flux in clean components for 1634 model */
#define F35_0 221.225  /* Total flux in clean components for 1634 model */
#define RAB 2.0244     /* Best-fit season 1 flux density ratio btwn A and B */
#define RCB 1.0376     /* Best-fit season 1 flux density ratio btwn C and B */
#define RDB 0.3482     /* Best-fit season 1 flux density ratio btwn D and B */
#define AMEAN 34.15    /* Mean flux density of component A light curve */
#define BMEAN 16.65    /* Mean flux density of component B light curve */
#define CMEAN 17.24    /* Mean flux density of component C light curve */
#define DMEAN 5.92     /* Mean flux density of component D light curve */
#define ASTART 	374    /* Starting date for A configuration (MJD - 50000) */
#define AEND 461       /* Ending date for A configuration (MJD - 50000) */
#define BSTART 494     /* Starting date for B configuration (MJD - 50000) */
#define FLATCHOICE 3   /* Determines how to flatten the 1608 light curves */
#define WINFRAC 1.0    /* Fraction used to determine number of ind. points */

/*.......................................................................
 *
 * Structure definitions
 *
 */

typedef struct {
  float val0;       /* Central value for parameter search */
  float dval;       /* Step size for parameter search */
  int nval;         /* Number of steps for parameter search */
  int minstep;      /* Minimum step number for parameter search */
  int maxstep;      /* Maximum step number for parameter search */
} Prange;

/*.......................................................................
 *
 * Function declarations
 *
 */

Prange *new_prange(int size);
Prange *del_prange(Prange *prange);
Fluxrec *new_fluxrec(int size);
Fluxrec *del_fluxrec(Fluxrec *fluxrec);
Fluxrec **load_light_curves(Setup *setup, int *fresult);
Fluxrec *read_fluxrec_1curve(char *inname, char comment, int *nlines);
int write_fluxrec(Fluxrec *fluxrec, int npoints, char *outfile,
		  int no_dup_days, float dtol);

void print_log(FILE *chifp, FILE *xcfp, char *logstring);
int read_data(Fluxrec *fl34, Fluxrec *fl35, Fluxrec *fl08[], 
	      int nlines, FILE *lfp, FILE *cfp);
int read_1608(Fluxrec *fl08[], int nlines, FILE *lfp, char *filename);
int write_1608(Fluxrec *flux[], int npoints, char *filename, int verbose);
int read_bad_ext(char *inname, char comment, Fluxrec *bada, Fluxrec *badb,
		 Fluxrec *badc, Fluxrec *badd);
int flag_bad(Fluxrec *fl34, Fluxrec *fl35, Fluxrec *fl08[], int nlines,
	     int *flagbad, char *badfilename, int *nbad);
float *make_flat(Fluxrec *fl34, Fluxrec *fl35, int nlines, int meanchoice);
Fluxrec *flat_field(Fluxrec *raw, float *flat, int nlines, float fracrms,
		    char *outname, int doprint);
Fluxrec *norm_constant(Fluxrec *raw, int nlines, float constant);
Fluxrec *norm_curve(Fluxrec *raw, int nlines, char *source, int doprint);
Fluxrec *norm_zero_mean(Fluxrec *raw, int nlines);
Fluxrec *norm_config(Fluxrec *raw, int nlines, char *source, int doprint);
Fluxrec *calc_flrat(Fluxrec *flux1, Fluxrec *flux2, int size);
int set_mu0(Fluxrec *fl08[], int nlines, Setup *setup);
int set_mu_grid(Fluxrec *lc[], int *npoints, Setup *setup);
int set_tau_grid(Fluxrec *lc[], int *npoints, int *index, Setup *setup);
int ratio_err(Fluxrec *flux1, Fluxrec *flux2, int nlines, float *fracrms);
int calc_mean(Fluxrec *flux, int nlines, float *mean, float *rms);
int calc_mean_dt(Fluxrec *flux, int nlines, float *mean, float *rms,
		 float startday, float endday, float dayshift);
int calc_swmean(Fluxrec *flux, int nlines, float *mean);
int fit_1608(Fluxrec *flux[], int npoints, int nbad[], float dt);
int fit_delays(Fluxrec *flux[], int npoints, float mean[], float d0[], 
	       float dtmax, float dt, int ndeg, Fluxrec *compos, int ncompos, 
	       int *first, float *bestchi, Fluxrec bestpars[]);
float *fit_poly(Fluxrec *flux, int npoints, int ndeg, int verbose);
void poly(float x, float p[], int ndeg);
void sinpoly(float x, float p[], int ndeg);
float chisq_fit(Fluxrec *flux, int npoints, float *a, int ncoeff, int nhidden,
		int verbose);
int fit_parab(float x1, float y1, float x2, float y2, float x3, float y3,
	      float *a, float *b, float *c, int doprint);
int daycmp(const void *v1, const void *v2);
int fluxcmp(const void *v1, const void *v2);
Fluxrec *make_compos(Fluxrec *flux[], int ncurves, int *npoints, 
		     int *index, float *lag, float *mu, int *ncompos,
		     int verbose);

#endif
