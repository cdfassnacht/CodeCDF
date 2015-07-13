/*
 * mk_monte_lcurve.c
 *
 * Usage: mk_monte_lcurve [1608 filename] [ideal filename]
 *
 * Creates fake light curves based on the idealized 1608 light curve.
 *
 * 08Jul98 CDF,  Modification of lcurve.c
 * v21Jul98 CDF, Added flat-fielding of the raw curves, which gives them
 *               the proper error bars for the Monte Carlo curve generation.
 * v05Sep98 CDF, Inserted improved bad-day flagging.
 *               Made number of curves interactive.
 * v23Feb99 CDF, Made sigma for adding random noise to curves interactive.
 * v22Jul99 CDF, Changed input to be previously flat-fielded curve (produced
 *                by fluxplot.sm) and took out the flat-fielding section
 *                of this program.
 * v02Apr02 CDF, Put some of the functionality into new functions in the
 *                monte_setup.c library.
 *               Put many variables into the newly-defined MC_Setup structure.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "nr.h"
#include "lc_funcs.h"
#include "monte_setup.h"
#include "monte.h"

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;                      /* Looping variable */
  int no_error=1;             /* Flag set to 0 on error */
  char lensfile[MAXC];        /* Input file name for 1608 data */
  char idealfile[MAXC];       /* Name of ideal curve file */
  char setupfile[MAXC];       /* Name of setup file */
  MC_Setup *setup=NULL;       /* Container for setup information */
  Fluxrec *fl08[N08]={NULL};  /* 1608 light curves */
  Fluxrec *ideal=NULL;        /* Ideal light curve */
  FILE *lfp=NULL;             /* 1608 data file pointer */

  /*
   * Check input line
   */

  if(argc < 3) {
    fprintf(stderr,"\nUsage: mk_monte_lcurve [1608 filename] ");
    fprintf(stderr,"[ideal curve filename] ([setup filename])\n\n");
    fprintf(stderr," The setup file is optional.\n\n");
    return 1;
  }

  /*
   * Initialize the Setup container
   */

  if(!(setup = new_MC_Setup(1))) {
    fprintf(stderr,"ERROR. Exiting program.\n");
    return 1;
  }

  /*
   * Put setup parameters into setup structure from setup file
   */

  if(argc == 4) {
    strcpy(setupfile,argv[3]);
    if(MC_Setup_file(setup,setupfile))
      no_error = 0;
  }

  /*
   * Fill in parts of the setup structure that weren't filled in
   *  from setup file.
   */

  if(no_error)
    if(setup_monte(setup))
      no_error = 0;

  /*
   * Summarize setup parameters
   */

  if(no_error)
    setup_monte_summary(setup);

  /*
   * Open 1608 and  files.
   */

  strcpy(lensfile,argv[1]);

  printf("\n");
  if(!(lfp = open_readfile(lensfile)))
    no_error = 0;

  /*
   * First count the number of data input lines in order to set
   *  sizes of arrays.
   */

  setup->nobs = n_lines(lfp,'#');

  if(setup->nobs == 0) {
    fprintf(stderr,"ERROR. No valid data lines.\n");
    no_error = 0;
  }
  else {
    printf("\n%d data lines in %s.\n\n",setup->nobs,lensfile);
    rewind(lfp);
  }

  /*
   * Allocate arrays
   */

  if(no_error)
    for(i=0; i<N08; i++) {
      if(!(fl08[i] = new_fluxrec(setup->nobs)))
	no_error = 0;
    }

  /*
   * Read in data
   */

  if(no_error)
    if(read_1608(fl08,setup->nobs,lfp,lensfile))
      no_error = 0;

  /*
   * Fill the idealized light curve container
   */

  strcpy(idealfile,argv[2]);
  if(!(ideal = read_fluxrec(idealfile,'#',&setup->nideal)))
    no_error = 0;

  /*
   * Get number of realizations to create
   */

  if(no_error && setup->nout == 0) {
    setup->nout = 1000;   /* Default value */
    if(get_nout(setup))
      no_error = 0;
  }

  /*
   * Generate fake curves with make_monte
   */

  if(no_error) {
    setup->ftol = 0.05;
    if(make_monte(ideal,fl08,setup,setup->nout))
      no_error = 0;
  }

  /*
   * Clean up
   */

  printf("\nCleaning up\n\n");
  ideal = del_fluxrec(ideal);
  for(i=0; i<N08; i++) {
    fl08[i] = del_fluxrec(fl08[i]);
  }

  if(lfp)
    fclose(lfp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: Exiting program\n");
    return 1;
  }
}

