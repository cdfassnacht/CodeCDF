#ifndef correlate_h
#define correlate_h

#include "lc_funcs.h"
#include "lc_interp.h"

#define CORRFUNC 3     /* Choice of method for doing cross-correlation */
#define ZEROFAC 1      /* Ratio of sizes of zero-padded and input arrays */

/*.......................................................................
 *
 * Function declarations
 *
 */

int call_xcorr_fft(Fluxrec *flux[], int size, char *filename, int nbad, 
		   FILE *logfp, int doprint);
Fluxrec *do_corr(Fluxrec *flux1, Fluxrec *flux2, int size, int *corsize,
		 int nbad);
Fluxrec *cross_corr_fft(float *flux1, float *flux2, int size, float dx);
Fluxrec *cross_corr_nr(float *flux1, float *flux2, int npoints, float dx);
int call_xcorr_time(Fluxrec *raw[], Fluxrec *interp[], int nraw, int ninterp, 
		    Setup *setup, FILE *logfp, char *filename, int doprint);
Fluxrec *cross_corr_time(Fluxrec *flux1, Fluxrec *flux2, int n1, int n2, 
			 float dx, Setup *setup);
int corr_bevington(Fluxrec *flux1, Fluxrec *flux2, int npoints, Fluxrec *r,
		   int nind);
double bevington_prob(Fluxrec maxcorr, int nind);
int corr_white_peterson(Fluxrec *flux1, Fluxrec *flux2, int npoints, float *r);
int call_acorr(Fluxrec *interp[], int size, char *filename, int nbad);
float *zero_pad(Fluxrec *flux, int npoints);
float *wt_hanning(Fluxrec *flux, int npoints);

#endif
