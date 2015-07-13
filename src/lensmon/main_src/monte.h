#ifndef monte_h
#define monte_h

#include "structdef.h"
#include "lc_setup.h"
#include "monte_setup.h"

/*.......................................................................
 *
 * Function declarations
 *
 */

int rand_curves(Fluxrec *fl08[], int nlines, int nbad, int ncurves);
int gauss_noise(Fluxrec *noise, int npoints, float sigma, long *randseed);
float *find_spacings(Fluxrec *obsflux, int npoints);
float *get_ideal_multfac(Fluxrec *mflux[], MC_Setup *setup, int ncomp);
int split_ideal(Fluxrec *ideal, Fluxrec *iflux[], MC_Setup *setup);
int make_monte(Fluxrec *ideal, Fluxrec *mflux[], MC_Setup *setup, 
	       int ncurves);
int rand_sampling(Fluxrec *obsflux, int npoints, float *samplings,
		     float *randsamp, int *nsamp, long *randseed);

#endif
