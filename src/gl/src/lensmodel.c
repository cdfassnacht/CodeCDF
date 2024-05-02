/*
 * model_fit.c
 *
 * Usage: model_fit [model_params_file] [image_pos_file] [source_flux_file]
 *  ** NB: The source flux file is optional **
 *
 * Description:  
 *
 * 06Apr96 CDF, First working version
 * v08Apr96 CDF, Added a rough grid search to improve initial guess
 * v19Apr96 CDF, Added possibility of computing flux-weighted chisq
 * v20Apr96 CDF, Change to GL file_io functions
 * v06Mar96 CDF, Added fitting to image fluxes
 * v28Jan97 CDF, Make pixel scale dependent on positional error
 * v27May97 CDF, Split off many steps into independent functions, which
 *                are also called by ray_trace and the other programs.
 *               Better documentation.
 * v06Aug97 CDF, Changed output file format in writemod
 * v08Aug97 CDF, Changed default gamma step size to 0.01
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structdef.h"
#include "file_io.h"
#include "list_tools.h"
#include "lensprogs.h"
#include "fitting.h"

#define PI 3.141592653589793
#define D2R (PI/180)
#define MAX 200
#define MAXNSTEP 9
#define DEFPOSVAR 0.05        /* Default dvar size for positions (arcsec) */
#define DEFGAMVAR 0.01        /* Default dvar size for gamma */
#define DEFTHETAVAR (5*D2R)  /* Default dvar size for theta (rad) */
#define DEFFLUXVAR 0.2        /* Default dvar size for source flux (mJy) */
#define FTOL 1.0e-7        /* Tolerance for downhill simplex fitting */

/*.......................................................................
 *
 * Function declarations
 *
 */

float *fill_stepsize(char *lensname, char *srcname);
int get_stepinfo(char *filename, float *sptr, int nsteps);
int fill_var(float *modpars, float *errarray, float *stepsize,
	     Modinfo *modinfo);
int struct_to_modpars(float valptr, int *varptr, int *maxvar, float *mpptr,
		      float *erraptr, float *stepsize, int j);
int writemod(char *name, Modinfo *modinfo, float *chisqlist);

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;
  int no_error=1;          /* Flag set to 0 on error */
  int niter;               /* Number of downhill simplex iterations */
  int nvar;                /* Number of varying parameters */
  int maxvar;
  int *intptr=NULL;        /* General pointer to integer arrays */
  float *nzfact=NULL;      /* Array of redshift factors */
  float *modpars=NULL;     /* Array of varying model parameters */
  float *stepsize=NULL;    /* Array of simplex spacings */
  float *errarray=NULL;    /* Array of errors for varying parameters */
  float *chisqlist=NULL;   /* Array of chisq values for the simplex points */
  char line[MAX];          /* String variable used to read input */
  char name[MAX];          /* Source name */
  Modinfo *modinfo=NULL;   /* Model, image and source info */

  /*
   * Check input on command line.
   */

  if (argc < 3 || argc > 4) {
    fprintf(stderr,"Usage:  lensmodel [image_pos_file] [model_params_file]");
    fprintf(stderr," ([source_file])\n");
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
   * Determine number of varying parameters
   */

  if(no_error) {
    intptr = modinfo->var;
    nvar = 0;
    maxvar = 1;
    for(i=0; i<(5*modinfo->nmod+4*modinfo->nfsrc); i++,intptr++)
      if (*intptr==1)
	nvar++;
      else if (*intptr>maxvar) {
	nvar++;
	maxvar = *intptr;
      }
    modinfo->nlinks = maxvar;
    modinfo->nvar = nvar;
    if(nvar == 0) {
      fprintf(stderr,"\nERROR: 0 varying components.  No model fitting\n");
      no_error = 0;
    }
  }

  /*
   * Allocate memory for arrays of varying parameters and errors on those
   *  parameter values.
   */

  if(no_error)
    if(!(modpars = new_array(nvar,nvar+1)))
      no_error = 0;

  if(no_error)
    if(!(errarray = new_array(nvar,1)))
      no_error = 0;

  /*
   * Fill the step size array
   */

  if(!(stepsize = fill_stepsize(argv[2],argv[3]))) {
    fprintf(stderr,"Exiting program\n");
    return 1;
  }

  if(fill_var(modpars,errarray,stepsize,modinfo)) {
    fprintf(stderr,"Exiting program\n");
    return 1;
  }

  /*
   * Do an optional grid of initial guesses, each followed by a
   *  downhill simplex fitting.  The values of modpars returned
   *  by mod_grid_search are the best-fit initial-guess values from 
   *  the entire grid of intial guesses.  If the grid search is done
   *  it still must be followed by the downhill simplex fitting to
   *  get the final best-fit model
   */

  if(no_error) {
    printf("Do you want to do a preliminary grid search (y/n)? [y] ");
    gets(line);
    if(strcmp(line,"")==0 || line[0] == 'Y' || line[0] == 'y') {
      if(mod_grid_search(modpars,errarray,modinfo,1))
	no_error = 0;
    }
  }

  /*
   * Set up simplex for final model-fitting
   */

  if(no_error)
    if(!(chisqlist = simplex_setup(modpars,nvar,errarray,mod_chisq,
				   (void *) modinfo,1)))
      no_error = 0;

  /*
   * Do the simplex search through model space.
   */

  if(no_error) {
    printf("\nStarting downhill simplex....\n\n");
    amoeba(modpars,chisqlist,modinfo->nvar,FTOL,mod_chisq,&niter,
	   (void *) modinfo,1);
  }

  /* 
   * Get information from fitting routine into Modinfo 
   */

  if(no_error)
    if(modpars_to_modinfo(modpars,modinfo))
      no_error = 0;

  /*
   * Write out the results
   */

  if(no_error)
    if(writemod(name,modinfo,chisqlist))
      no_error = 0;

  /*
   * Clean up
   */

  nzfact = del_array(nzfact);
  stepsize = del_array(stepsize);
  errarray = del_array(errarray);
  modpars = del_array(modpars);
  chisqlist = del_array(chisqlist);
  modinfo = del_modinfo(modinfo);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: Exiting program\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function fill_var
 *
 * Fills in the containers for the varying parameters of the model as well
 *   as the "error" arrays that set the size of the simplex and the
 *   stepsize in the (optional) grid search.
 *
 * Inputs:  float *modpars     list of varying parameter values (set by
 *                               this function)
 *          float *errarray    array that sets simplex basis vector sizes
 *                               (set by this function)
 *          float *stepsize    array containing step size information
 *          Modinfo *modinfo   structure holding all other model information
 *
 * Output:  int (0 or 1)       0 ==> successful, 1 ==> error
 *
 */

int fill_var(float *modpars, float *errarray, float *stepsize,
	     Modinfo *modinfo)
{
  int i,j;
  int maxvar;
  int *intptr;
  float *fptr,*fptr2;
  Lensparams *parptr;
  Srcinfo *srcptr;

  intptr = modinfo->var;
  fptr = modpars;
  fptr2 = errarray;
  maxvar = 1;
  parptr = modinfo->lp;
  for(i=0; i<modinfo->nmod; i++,parptr++) {
    j = 0;
    switch(struct_to_modpars(parptr->xl.par,intptr++,&maxvar,fptr,
			      fptr2,stepsize,j++)) {
    case 0:
      break;
    case 1:
      fptr++;
      fptr2++;
      break;
    default:
      fprintf(stderr,"Error: fill_var\n");
      return 1;
    }
    switch(struct_to_modpars(parptr->yl.par,intptr++,&maxvar,fptr,
			     fptr2,stepsize,j++)) {
    case 0:
      break;
    case 1:
      fptr++;
      fptr2++;
      break;
    default:
      fprintf(stderr,"Error: fill_var\n");
      return 1;
    }
    switch(struct_to_modpars(parptr->b.par,intptr++,&maxvar,fptr,
			     fptr2,stepsize,j++)) {
    case 0:
      break;
    case 1:
      fptr++;
      fptr2++;
      break;
    default:
      fprintf(stderr,"Error: fill_var\n");
      return 1;
    }
    switch(struct_to_modpars(parptr->gamma.par,intptr++,&maxvar,fptr,
			     fptr2,stepsize,j++)) {
    case 0:
      break;
    case 1:
      fptr++;
      fptr2++;
      break;
    default:
      fprintf(stderr,"Error: fill_var\n");
      return 1;
    }
    switch(struct_to_modpars(parptr->theta_g.par,intptr++,&maxvar,fptr,
			     fptr2,stepsize,j++)) {
    case 0:
      break;
    case 1:
      fptr++;
      fptr2++;
      break;
    default:
      fprintf(stderr,"Error: fill_var\n");
      return 1;
    }
  }

  srcptr = modinfo->fsource;
  for(i=0; i<modinfo->nfsrc; i++,srcptr++) {
    j = 5;
    switch(struct_to_modpars(srcptr->x.par,intptr++,&maxvar,fptr,
			     fptr2,stepsize,j++)) {
    case 0:
      break;
    case 1:
      fptr++;
      fptr2++;
      break;
    default:
      fprintf(stderr,"Error: fill_var\n");
      return 1;
    }
    switch(struct_to_modpars(srcptr->y.par,intptr++,&maxvar,fptr,
			     fptr2,stepsize,j++)) {
    case 0:
      break;
    case 1:
      fptr++;
      fptr2++;
      break;
    default:
      fprintf(stderr,"Error: fill_var\n");
      return 1;
    }
    switch(struct_to_modpars(srcptr->flux.par,intptr++,&maxvar,fptr,
			     fptr2,stepsize,j++)) {
    case 0:
      break;
    case 1:
      fptr++;
      fptr2++;
      break;
    default:
      fprintf(stderr,"Error: fill_var\n");
      return 1;
    }
    switch(struct_to_modpars(srcptr->zf.par,intptr++,&maxvar,fptr,
			     fptr2,stepsize,j++)) {
    case 0:
      break;
    case 1:
      fptr++;
      fptr2++;
      break;
    default:
      fprintf(stderr,"Error: fill_var\n");
      return 1;
    }
  }

  printf("\nThere are %d variable parameters.\n",modinfo->nvar);
  fptr = modpars;
  fptr2 = errarray;
  for (i=0; i<modinfo->nvar; i++,fptr++,fptr2++) {
    printf("%d %8.4f  err=%f\n",i+1,*fptr,*fptr2);
  }
  printf("\n");

  return 0;
}

/*.......................................................................
 *
 * Function struct_to_modpars
 *
 * Takes the varying parameters of a structure and puts them into the
 *   list of values to be passed to the fitting routine
 *
 * Inputs:  float *valptr      pointer to the structure component
 *          int *varptr        tells whether the component was modified
 *                               by the model fitting
 *          int *maxvar        keeps track of varying parameters
 *          float *mpptr       pointer to the model parameters
 *          float *erraptr     pointer to the "errors" associated with
 *                               the varying parameters
 *          float *stepsize    values of the "errors" used to determine
 *                               the size of the simplex bases
 *          int j              index for stepsize
 *
 * Output:  int (0, 1 or -1)   0 ==> don't increment the pointers
 *                             1 ==> increment the pointers
 *                             -1 ==> error
 *
 */

int struct_to_modpars(float valptr, int *varptr, int *maxvar, float *mpptr,
		      float *erraptr, float *stepsize,int j)
{
  if (*varptr==1) {
    *mpptr = valptr;
    *erraptr = stepsize[j];
    return 1;
  }
  else if (*varptr > *maxvar) {
    *maxvar = *varptr;
    *mpptr = valptr;
    *erraptr = stepsize[j];
    return 1;
  }
  else
    return 0;

  return 0;
}

/*.......................................................................
 *
 * Function fill_stepsize
 *
 * Reads the model parameter input files (one for the lens, one for the
 *  source -- if the source file exists), and extracts from them the
 *  step sizes to be used in constructing any simplex and or grid
 *  in parameter space.  The step sizes will be on a separate line in
 *  the input files, with a '#' at the beginning of the line, followed
 *  by the word 'steps'.
 *
 * Inputs: char *lensname      name of lens input file
 *         char *srcname       name of source input file (NULL if no file)
 *
 * Output: float *stepsize     array containing step sizes, NULL on error
 *
 */

float *fill_stepsize(char *lensname, char *srcname)
{
  int lensread=0;          /* Flag set to 1 when lens info has been read */
  int srcread=0;           /* Flag set to 1 when source info has been read */
  float *stepsize=NULL;    /* Array of step sizes */
  float *sptr;             /* Pointer to navigate stepsize */

  /*
   * Allocate memory for the step sizes
   */

  if(!(stepsize = new_array(MAXNSTEP,1))) {
    fprintf(stderr,"ERROR: fill_stepsize\n");
    return NULL;
  }

  /*
   * Read in the lens info.
   * If the "steps" line was not found, fill in the lens section of
   *  the stepsize array with the default parameters
   */

  sptr = stepsize;
  switch(lensread = get_stepinfo(lensname,sptr,5)) {
  case 1:
    *(stepsize+4) *= D2R;
    break;
  case 0:
    printf("\nUsing default step size values for lens model\n");
    sptr = stepsize;
    *sptr++ = DEFPOSVAR;
    *sptr++ = DEFPOSVAR;
    *sptr++ = DEFPOSVAR;
    *sptr++ = DEFGAMVAR;
    *sptr++ = DEFTHETAVAR;
    break;
  default:
    fprintf(stderr,"ERROR: fill_stepsize\n");
    return del_array(stepsize);
  }

  /*
   * Now read in the source information (if any)
   */

  if(srcname) {
    switch(srcread = get_stepinfo(srcname,sptr,4)) {
    case 1:
      break;
    case 0:
      printf("\nUsing default step size values for source model\n");
      sptr = stepsize+5;
      *sptr++ = DEFPOSVAR;
      *sptr++ = DEFPOSVAR;
      *sptr++ = DEFFLUXVAR;
      *sptr++ = DEFFLUXVAR;
      break;
    default:
      fprintf(stderr,"ERROR: fill_stepsize\n");
      return del_array(stepsize);
    }
  }
  else {
    printf("\nUsing default step size values for source model\n");
    sptr = stepsize+5;
    *sptr++ = 0.3;
    *sptr++ = 0.3;
    *sptr++ = 0.0;
    *sptr++ = 0.0;
  }

  return stepsize;
}

/*.......................................................................
 *
 * Function get_stepinfo
 *
 * Given an input filename (either for the lens info or the source info)
 *  this function opens the file and reads in the step size information
 *  from it.  The step sizes will be on a separate line in the input 
 *  file, with a '#' at the beginning of the line, followed by the 
 *  word 'steps'.
 *
 * Inputs: char *filename      name of input file
 *         float *sptr         pointer to current position in stepsize
 *                              array
 *         int nsteps          number of steps expected on the input line.
 *
 * Output: int inforead        flag set to 1 if info found and read
 *                             on error, -1 is returned
 *
 */

int get_stepinfo(char *filename, float *sptr, int nsteps)
{
  int inforead=0;      /* Flag set to 1 if info found and read */
  char line[MAX];      /* General string for reading input */
  char keyword[MAX];   /* First word on line of input file */
  FILE *ifp=NULL;      /* Input file pointer */

  /*
   * Open the lens file
   */

  strcpy(line,filename);
  while(!(ifp = fopen(line,"r"))) {
    printf("%s is not a valid file name.  Enter the file name again:  ",
	   line);
    gets(line);
  }

  /*
   * Read the input lines
   */

  while(fgets(line,MAX,ifp) != NULL && !inforead) {
    if(line[0] == '#') {

      /*
       * Check to see if there is more than just a '#' on the input line
       */

      if(sscanf(line,"#%s",keyword) == 1) {

	/*
	 * See if keyword matches 'steps'
	 */

	if(strcmp(keyword,"steps") == 0)
	  if(nsteps == 5) {
	    if(sscanf(line,"#%s %f %f %f %f %f",keyword,sptr++,sptr++,
		      sptr++,sptr++,sptr++) != 6) {
	      fprintf(stderr,"ERROR:  get_stepinfo.  Bad input line\n");
	      return -1;
	    }
	    inforead = 1;
	  }
	  else if(nsteps == 4) {
	    if(sscanf(line,"#%s %f %f %f %f",keyword,sptr++,sptr++,
		      sptr++,sptr++) != 5) {
	      fprintf(stderr,"ERROR:  get_stepinfo.  Bad input line\n");
	      return -1;
	    }
	    inforead = 1;
	  }
      }
    }
  }

  return inforead;
}

/*.......................................................................
 *
 * Function writemod
 *
 * Writes the best-fit model information to two output files, one for
 *  the lens parameters and one for the source parameters.
 *
 * Inputs: char *name          source name
 *         Modinfo *modinfo    model information
 *         float *chisqlist    array containing chisq of best-fit model
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int writemod(char *name, Modinfo *modinfo, float *chisqlist)
{
  int i;                /* Looping variable */
  char outname[MAX];    /* String to get output file names */
  Lensparams *parptr;   /* Pointer to navigate modinfo->lp */
  Srcinfo *srcptr;      /* Pointer to navigate modinfo->fsource */
  FILE *mfp=NULL;       /* Model file pointer */
  FILE *sfp=NULL;       /* Source file pointer */

  /*
   * Get output file names.
   */

  printf("\nOutput file name for model parameters:  ");
  gets(outname);
  while(!(mfp = fopen(outname,"w"))) {
    fprintf(stderr,"Error in opening output file.  Enter filename again:  ");
    gets(outname);
  }

  printf("Output file name for source flux(es):  ");
  gets(outname);
  while(!(sfp = fopen(outname,"w"))) {
    fprintf(stderr,"Error in opening output file.  Enter filename again:  ");
    gets(outname);
  }

  /*
   * Write output
   */

  printf("Writing output files.\n\n");
  fprintf(mfp,"%s\n",name);
  printf("\nFinal model for %s.  chisq = %f\n\n",name,chisqlist[0]);
  printf("Type     x_l      y_l      b       gamma    theta_g\n");
  printf("-----  -------  ------- -------- --------- ---------\n");
  fprintf(mfp,"#Type      x_l        y_l       b         gamma    ");
  fprintf(mfp," theta_g\n");
  fprintf(mfp,"#-----  ---------  --------- ---------- -----------");
  fprintf(mfp," ---------\n");
  parptr = modinfo->lp;
  for(i=0;i<modinfo->nmod;i++,parptr++) {
    fprintf(mfp,"%4d %10.4f %10.4f %10.5f %10.5f %10.1f\n",
	    parptr->pottyp,parptr->xl.par,parptr->yl.par,parptr->b.par,
	    parptr->gamma.par,(parptr->theta_g.par)/D2R);
    printf("%3d %10.4f %8.4f %8.5f %8.5f %8.1f\n",
	    parptr->pottyp,parptr->xl.par,parptr->yl.par,parptr->b.par,
	    parptr->gamma.par,(parptr->theta_g.par)/D2R);
  }

  printf("\n   x        y     Flux     beta\n");
  printf("-------  ------- -------  --------\n");
  fprintf(sfp,"# %s\n",name);
  fprintf(sfp,"#    x          y       Flux       beta\n");
  fprintf(sfp,"#---------  --------- ---------  ---------\n");
  srcptr = modinfo->fsource;
  for(i=0;i<modinfo->nfsrc;i++,srcptr++) {
    fprintf(sfp,"%7.4f %10.4f %9.2f %9.4f\n",srcptr->x.par,srcptr->y.par,
	    srcptr->flux.par,srcptr->zf.par);
    printf("%7.4f %8.4f %7.2f %8.4f\n",srcptr->x.par,srcptr->y.par,
	   srcptr->flux.par,srcptr->zf.par);
  }

  /*
   * Clean up
   */

  if(mfp)
    fclose(mfp);
  if(sfp)
    fclose(sfp);

  return 0;
}

/*.......................................................................
 *
 * Dummy function for pgplot
 *
 */

 int MAIN_(void)
{
  return 0;
}
