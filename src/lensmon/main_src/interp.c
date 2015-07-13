/*
 * interp.c
 *
 * Usage: interp [1608 filename]  [output filename] ([setup filename]),
 *   where the output file contains the smoothed and/or interpolated
 *   light curve.  The setup file is optional.
 *
 * Interpolates 1608 light curves that have already been "flat-fielded".
 *  The flat-fielding is done in the fluxplot.sm macros.
 *
 * 20Jul99 CDF,  A modification of lcurve_mc.c.
 * v16Sep99 CDF, Slight modifications to deal with new version of lc_setup.c
 * v06Sep00 CDF, Changed output format slightly.
 * v13Sep00 CDF, Changed call to interp_switch to deal with changes in
 *                lc_interp.c.
 * v04Oct00 CDF, Moved printing of output file into new write_1608 function
 *                in lc_funcs.c.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "nr.h"
#include "lc_setup.h"
#include "lc_funcs.h"
#include "lc_interp.h"

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;                      /* Looping variable */
  int no_error=1;             /* Flag set to 0 on error */
  int nlines;                 /* Number of data lines in lens input file */
  int nbad[N08]={0,0,0,0};    /* Number of bad days in bad day file */
  char lensfile[MAXC];        /* Input file name for 1608 data */
  char outname[MAXC];         /* Name for smoothed curve output file */
  char setupfile[MAXC];       /* Name of setup file */
  Fluxrec *fl08[N08]={NULL};  /* 1608 light curves */
  Fluxrec *gr08[N08]={NULL};  /* Gridded 1608 curves (smoothed) */
  Setup *setup=NULL;          /* Container for setup information */
  FILE *lfp=NULL;             /* 1608 data file pointer */
 
  /*
   * Check input line
   */

  if(argc < 3) {
    fprintf(stderr,"\nUsage: interp [1608 filename] [output filename]");
    fprintf(stderr," ([setup filename]),\n\n");
    fprintf(stderr," The input 1608 file has already been");
    fprintf(stderr," flat-fielded (in fluxplot.sm) and\n");
    fprintf(stderr," where the output curve contains the smoothed and/or");
    fprintf(stderr," interpolated light curves.\n");
    fprintf(stderr," The setup file is optional.\n\n");
    return 1;
  }

  /*
   * Initialize the Setup container
   */

  if(!(setup = new_setup(1))) {
    fprintf(stderr,"ERROR. Exiting program.\n");
    return 1;
  }

  /*
   * Put setup parameters into setup structure from setup file
   */

  if(argc == 4) {
    strcpy(setupfile,argv[3]);
    if(setup_file(setup,setupfile))
      no_error = 0;
  }

  /*
   * Fill in parts of the setup structure that weren't filled in
   *  from setup file.
   */

  if(no_error)
    if(setup_interp(setup))
      no_error = 0;

  /*
   * Summarize setup parameters
   */

  if(no_error)
    setup_interp_summary(setup);

  /*
   * Open lens light curve file 
   */

  strcpy(lensfile,argv[1]);
  if(!(lfp = open_readfile(lensfile)))
    no_error = 0;

  /*
   * First count the number of data input lines in order to set
   *  sizes of arrays.
   */

  nlines = n_lines(lfp,'#');

  if(nlines == 0) {
    fprintf(stderr,"ERROR. No valid data lines.\n");
    no_error = 0;
  }
  else {
    printf("\n%d data lines in %s.\n\n",nlines,lensfile);
    rewind(lfp);
  }

  /*
   * Allocate arrays
   */

  if(no_error)
    for(i=0; i<N08; i++) {
      if(!(fl08[i] = new_fluxrec(nlines)))
	no_error = 0;
    }

  /*
   * Read in data
   */

  if(no_error)
    if(read_1608(fl08,nlines,lfp,lensfile))
      no_error = 0;

  /*
   * Set parameters for interpolated curves.
   */

  if(no_error)
    set_grid_params(setup,fl08[0],nlines);

  /*
   * Smooth or linearly interpolate the flattened curves, depending on
   *  flag in setup.
   */

  if(setup->dosmooth != NOSMOOTH && no_error) {
    printf("\nSmoothing and/or interpolating light curves...\n");
    for(i=0; i<N08; i++) {
      if(!(gr08[i] = interp_switch(fl08[i],nlines,setup,nbad[i])))
	no_error = 0;
    }
  }

  /*
   * Print out smoothed and/or interpolated curves
   */

  if(no_error && setup->dosmooth != NOSMOOTH) {
    printf("\nPrinting out smoothed light curves.\n");
    strcpy(outname,argv[2]);
    printf("  ");
    if(write_1608(gr08,setup->ninterp,outname,1))
      no_error = 0;
  }

  /*
   * Clean up
   */

  if(no_error)
    printf("\nCleaning up\n");

  for(i=0; i<4; i++) {
    fl08[i] = del_fluxrec(fl08[i]);
    gr08[i] = del_fluxrec(gr08[i]);
  }
  setup = del_setup(setup);

  if(lfp)
    fclose(lfp);

  if(no_error) {
    printf("\nFinished with program interp.c.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"ERROR: Exiting program\n");
    return 1;
  }
}
