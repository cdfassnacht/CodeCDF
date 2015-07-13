/*
 * best_cmod.c
 *
 * This program takes the output from the difmap modelfitting scripts for
 *  the lens monitoring program and finds the fluxes associated with the
 *  best fit models.  The best fit models are written to model files.
 *
 * Usage: best_cmod [input date] [source_root]
 *
 * 02Apr97 CDF,  A modification of the old getflux.c
 * v14Apr97 CDF, Accounting for post 21Mar97 1635 models.
 * v10May01 CDF, A big re-write to make more general.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dataio.h"
#include "modfuncs.h"

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int no_error=1;     /* Flag set to 0 on error */
  float scal;         /* Best-fit scale factor */
  char obsdate[MMAX]; /* Epoch of observation */
  char root[MMAX];    /* Source root name */
  char inmod[MMAX];   /* Input "scalmod" filename */
  char outmod[MMAX];  /* Output model filename */
  char basemod[MMAX]; /* Base model filename */
  FILE *mfp=NULL;     /* Model output file */

  /*
   * Check command line format
   */

  if(argc != 3) {
    fprintf(stderr,"\nUsage: best_cmod [obsdate] [source_root]\n\n");
    fprintf(stderr," e.g. For observations of 1633+741 on Feb 23, 1998\n");
    fprintf(stderr,"      type \"best_cmod 19980223 1633\"\n\n");
    return 1;
  }

  /*
   * Get relevant info from command line
   */

  strcpy(obsdate,argv[1]);
  strcpy(root,argv[2]);

  /*
   * Set filenames for input "scalmod" file and output model file.
   */

  sprintf(inmod,"%s_%s.scalmod",obsdate,root);
  sprintf(outmod,"%s_%s.mod",obsdate,root);

  /*
   * Get best-fit scale
   */

  if((cal_flux(inmod,&scal)) == 1) {
    fprintf(stderr,"ERROR: best_cmod\n");
    return 1;
  }

  /*
   * Open base model file
   */


  /*
   * Create best-fit model
   */

  sprintf(basemod,"%s%s_mod_0",MODDIR,root);
  if(scale_data(scal,scal,1.0,basemod))
      no_error = 0;

  /*
   * Clean up and exit
   */

  if(mfp)
    fclose(mfp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: best_cmod\n");
    return 1;
  }
}
