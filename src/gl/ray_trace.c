/*
 * ray_trace.c
 *
 * Usage: ray_trace [image_pos_file] [model_params_file] ([source_file])
 *
 * Description:  Sets up a potential and traces rays through it to
 *                determine image position, using a known model.
 *                Given lens parameters and observed image positions,
 *                the program calculates the potential on a grid, the 
 *                time delay contours and finds the image positions 
 *                determined by the potential.
 *
 * 08Feb96, CDF First working version.
 * v15Feb96, CDF Calculate image magnifications
 * v18Feb96, CDF Put functions into libraries for greater versatility.
 *               Finds critical curves and caustics
 * v03Mar96, CDF Put model type into a function table for versatility.
 * v12Mar96, CDF Changed library structure for type definitions
 * v13Mar96, CDF Changed potentials to a more modular form, ie. instead
 *                of using an SIS+external or SIS+internal, set up
 *                SIS, internal and external so they can be combined in
 *                a general way.
 * v14Mar96, CDF Added the possibility of having multiple background
 *                sources at different redshifts being lensed.  Included
 *                a "redshift factor" which scales the potential depending
 *                on the redshift of the background source.
 * v15Mar96, CDF Changed list of model image positions to a list instead
 *                of a pointer array
 * v21Mar96, CDF Slight debugging
 * v09Apr96, CDF Add interactive plotting choices
 * v19Apr96, CDF Added possibility of computing flux-weighted chisq
 * v20Apr96, CDF Change to GL file_io functions
 * v24Apr96, CDF More sophisticated time-delay contouring
 * v14Jul97, CDF Massive re-write putting almost all calculations and plotting
 *                into functions in the gl library.
 *               More efficient I/O.
 * v14Nov97, CDF Moved plot_boxes call into lensprogs.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "file_io.h"
#include "list_tools.h"
#include "potentials.h"
#include "lensprogs.h"
#include "plotfuncs.h"
#include "cpgplot.h"

#define PI 3.141592653589793
#define D2R (PI/180)
#define DEFGRIDSIZE 128    /* Default grid size (pixels) */

int main(int argc, char *argv[])
{
  int i;                  /* Looping variable */
  int no_error=1;         /* Flag set to 0 on error */
  int size;               /* Size of grid */
  float *nzfact;          /* Redshift scaling factor */
  float *hess;            /* Hessian array */
  float *fptr;            /* General pointers to float arrays */
  char line[MAXC];        /* General character array for reading info */
  char name[MAXC];        /* Source name */
  Lensparams *parptr;     /* General pointer to parameters */
  Poten *poten=NULL;      /* Gravitational potential grid */
  Poten *impot=NULL;      /* Potential at image positions */
  Posinfo *grid=NULL;     /* The coordinate grid for the program */
  Posinfo *piptr;
  Modinfo *modinfo=NULL;  /* Container for model, image and source info */
  List imlist=NULL;       /* List containing model image positions */

  /*
   * Check input on command line.
   */

  if (argc < 3 || argc > 4) {
    fprintf(stderr,"Usage:  ray_trace [image_pos_file] [model_params_file]");
    fprintf(stderr," ([source_file])\n");
    return 1;
  }

  /*
   * Fill get model, image and source data from input files
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
   * Print out chisq for this model
   */

  if(no_error) {
    printf("\n** Chisq calculation for this model **\n");
    if(sp_chisq(modinfo->posobs->spos,modinfo,1) < 0.0)
      no_error = 0;
  }

  /* 
   * Set up image size
   */

  if(no_error) {
    size = DEFGRIDSIZE;
    printf("\nEnter size of field in pixels.\n");
    printf(" One value only, the field will be square [%d]:  ",size);
    gets(line);
    if(strcmp(line,"") != 0) {
      while(sscanf(line,"%d",&size) != 1 || size < 0) {
	fprintf(stderr,"ERROR: Invalid input\n");
	fprintf(stderr,"Enter size of field in pixels:  ");
	gets(line);
      }
    }
  }

  /*
   * Set up grid for potential calculations
   */

  if(no_error) {
    printf("\nSetting up coordinate grid...\n");
    if(!(grid = make_grid(size,modinfo->pixscale)))
      no_error = 0;
  }

  if(no_error) {
    printf("\nAllocating arrays....\n");
    if(!(poten = new_poten(size,size)))
      no_error = 0;
  }

  /* 
   * Calculate potential on grid for quick checks.
   */

  if(no_error) {
    printf("\nCalculating Potential Grid...\n");
    for(i=0,parptr=modinfo->lp; i<modinfo->nmod; i++,parptr++) 
      poten_table[parptr->pottyp].fn(poten,grid,parptr,size,size,0);
  }

  /*
   * Find rough image positions using grid potential.  Note that
   *  if modinfo->posimmod has been calculated already in finding
   *  the (unknown) source position, it needs to be freed.
   */

  if(no_error) {
    printf("\nFinding rough image positions...\n");
    modinfo->posimmod = del_posinfo(modinfo->posimmod);
    if(find_images_grid(modinfo,nzfact,size,grid,poten,
			0.1*modinfo->poserr,1,0))
      no_error = 0;
    else {
      printf("Rays through %d grid points hit within %f of source\n",
	     modinfo->nimmod,modinfo->poserr);
      printf(" position in source plane\n");
    }
  }

  /*
   * Use Newton's method for finding better image positions
   */

  if(no_error) {
    printf("\nFinding better image positions...\n");
    if(!(imlist = find_images_newton(modinfo,0)))
      no_error = 0;
  }

  /*
   * Eliminate duplicates from the list of better image positions and
   *  put final list into modinfo->posimmod.
   */

  if(no_error)
    if(condense_imlist(imlist,modinfo,1))
      no_error = 0;

  /*
   * Calculate potential at model image positions
   */

  if(no_error) {
    printf("\nCalculating potential at model image positions\n");
    if(!(impot = calc_poten(modinfo->posimmod,modinfo->nimmod,modinfo)))
      no_error = 0;
  }

  /*
   * Find image magnifications from the inverse Hessian
   */

  if(no_error) {
    printf("\nCalculating image magnifications....\n");
    if(!(hess = calc_hess(impot,modinfo->posimmod,modinfo->nimmod,1)))
      no_error = 0;
  }

  if(no_error) {
    printf(" Image      x        y     magnif.   flux\n");
    printf("-------  -------  -------  -------  -------\n");
    for(i=0,piptr=modinfo->posimmod,fptr=hess; i<modinfo->nimmod; 
	i++,piptr++,fptr++)
      printf("%5d    %7.4f  %7.4f  %7.3f  %7.2f\n",i+1,piptr->x,piptr->y,
	     1/(*fptr),piptr->sflux/fabs(*fptr));
  }

  /*
   * Plot things
   */

  printf("\n");
  plot_open(1.5);
  plot_labs("",size,modinfo->pixscale);

  /*
   * Check with the user about various plotting options.
   */

  if(no_error)
    if(plot_poten(poten,size,modinfo))
      no_error = 0;

  if(no_error) {
    printf("Plot time delay contours (y/n)? [y]  ");
    gets(line);
    if (strcmp(line,"") == 0 || line[0] == 'y' || line[0] == 'Y')
      if(plot_time_delay(grid,poten,size,modinfo,hess,impot))
	 no_error = 0;
  }

  /*
   * Find caustics using the grid
   */

  if(no_error) {
    printf("Plot caustics and critical curves (y/n)? [y]  ");
    gets(line);
    if (strcmp(line,"") == 0 || line[0] == 'y' || line[0] == 'Y') {
      printf("\nFinding caustic positions, please be patient...\n");
      if(find_caustics(grid,poten,size,modinfo))
	no_error = 0;
    }
  }

  /*
   * Plot model and observed image positions
   */

  if(no_error) {
    if(plot_modim(modinfo,grid,poten,size,nzfact))
      no_error = 0;
  }

  /*
   * Clean up and exit
   */

  plot_close();

  grid = del_posinfo(grid);
  poten = del_poten(poten);
  impot = del_poten(impot);
  hess = del_array(hess);
  nzfact = del_array(nzfact);
  modinfo = del_modinfo(modinfo);
  free_list(&imlist);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"\nERROR.  Exiting program.\n\n");
    return 1;
  }
}

int MAIN_(void)
{
  return 0;
}
