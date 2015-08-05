/*
 * find_delays.c
 *
 * Usage: find_delays [image_pos_file] [model_params_file]
 *
 * Description:  Given the lens parameters, the program calculates
 *                the predicted time delays between all the components.
 *
 * 28Feb96 CDF
 * v23Apr96 CDF, Include GL file io functions
 * v08Jul97 CDF, Consolidate I/O function into a modinfo I/O.  Better
 *                error checking and documentation.
 * v07Oct97 CDF, Calculate time delays at MODEL image positions rather
 *                than at observed image positions.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structdef.h"
#include "file_io.h"
#include "list_tools.h"
#include "lensprogs.h"

#define PI 3.141592653589793
#define D2R (PI/180)
#define MAXLINE 100

int main(int argc, char *argv[])
{
  int i;                    /* Looping variable */
  int no_error=1;           /* Flag set to 0 on error */
  float *nzfact=NULL;       /* Array of redshift factors */
  float *del_t=NULL;        /* Array of time delays */
  float zl,zs,omega0;       /* Lens and source redshifts, Omega_0 */
  char line[MAXLINE];       /* General string for reading info */
  char name[MAXLINE];       /* System name */
  Posinfo *piptr;           /* Pointer to navigate modinfo->posobs */
  Posinfo *pmptr;           /* Pointer to navigate modinfo->posimmod */
  Modinfo *modinfo=NULL;    /* Model, image and source info */
  List imlist;              /* List of images produced by lens model */

  /*
   * Check input on command line.
   */

  if (argc < 3) {
    fprintf(stderr,"\nUsage: find_delays [image_pos_file] ");
    fprintf(stderr,"[model_params_file] [source_info_file]\n\n");
    return 1;
  }

  /*
   * Fill modinfo structure with data from input files
   */

  if(!(modinfo = modinfo_io(argv[1],argv[2],argv[3],name))) {
    fprintf(stderr,"Exiting program\n");
    return 1;
  }

  /*
   * Check to make sure number of redshift factors in image_pos file matches
   *  the number of sources in the source file
   */

  if(!(nzfact = check_nsrc(modinfo)))
    no_error = 0;

  /*
   * Get information about pixel scale, errors and chisq calculation
   *  methods.
   */

  if(no_error)
    if(get_fitinfo(modinfo))
      no_error = 0;

  /*
   * If there is no source input file, calculate rough source position
   *  as initial starting point.
   */

  if(!(argv[3]) && no_error)
    if(source_guess(modinfo,nzfact))
      no_error = 0;

  /*
   * Get cosmological information
   */

  if(no_error) {
    printf("\nEnter redshift of lens (just hit return if unknown):  ");
    gets(line);
    if(strcmp(line,"") == 0)
      zl = 0.0;
    else
      sscanf(line,"%f",&zl);
    printf("Enter redshift of source (just hit return if unknown):  ");
    gets(line);
    if(strcmp(line,"") == 0)
      zs = 0.0;
    else
      sscanf(line,"%f",&zs);
    printf("Enter Omega_0 (Hit return for Omega_0 = 1.0):  ");
    gets(line);
    if(strcmp(line,"") == 0)
      omega0 = 1.0;
    else
      sscanf(line,"%f",&omega0);
  }

  /*
   * The model image positions need to be calculated.  Do this by starting 
   *  at observed image positions and letting Newton's method find optimal 
   *  image postions for the model.
   *
   * Copy observed image positions into model image positions, to be
   *  used as a starting guess for Newton's method.
   */

  if(no_error) {
    printf("\nUsing Newton's method to find model image positions\n");
    modinfo->nimmod = modinfo->nobs;
    if(!(modinfo->posimmod = new_posinfo(modinfo->nobs,1)))
      no_error = 0;
  }

  if(no_error) {
    for(i=0,piptr=modinfo->posobs,pmptr=modinfo->posimmod; i<modinfo->nobs;
	i++,piptr++,pmptr++)
      *pmptr = *piptr;
  }
  
  /*
   * Use Newton's method for finding image positions, given the lens 
   *  and starting with the observed image positions
   */

  if(no_error) {
    printf("\nFinding image positions predicted by lens model\n");
    if(!(imlist = find_images_newton(modinfo,1)))
      no_error = 0;
  }

  /*
   * Transfer info from list into modinfo->posimmod
   */

  if(no_error)
    if(condense_imlist(imlist,modinfo,1))
      no_error = 0;

  /*
   * Calculate delays
   */

  if(no_error) {
    printf("\n");
    for(i=0,pmptr=modinfo->posimmod; i<modinfo->nimmod; i++,pmptr++)
      printf("Model image %d located at %7.4f %7.4f\n",i+1,pmptr->x,
	     pmptr->y);
    printf("\n");
    if(!(del_t = predict_delay(modinfo,zl,zs,omega0)))
      no_error = 0;
  }

  /*
   * Clean up
   */

  free_list(&imlist);
  del_t = del_array(del_t);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR.  Exiting program\n\n");
    return 1;
  }
}

int MAIN_(void)
{
  return 0;
}
