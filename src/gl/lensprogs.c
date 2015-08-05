/*
 * lensprogs.c
 *
 * A collection of functions used in modelling and doing ray_tracing
 *  through gravitational lens potentials.
 *
 * 18Feb96 CDF, First working version
 * v14Jul97 CDF, Major rewrite, putting many lens info parameters into
 *                modinfo structures to make function calls easier.
 *               Found bugs in rough image and fine image location functions.
 * v05Aug97 CDF, Deleted an unnecessary temporary array in calc_sp
 * v07Aug97 CDF, Modified grid_iter to print out model parameters every
 *                time the minimum chisq improves.
 * v07Oct97 CDF, Changed predict_delay to use modinfo structure rather
 *                than individual lensparams, posinfo, etc. structures.
 *                Also put in better error checking and documentation.
 * v08Oct97 CDF, Have source_guess put the guessed source position(s)
 *                into the modinfo->posobs structure.
 *               Added function calc_poten to calculate potentials at
 *                positions contained in a Posinfo structure array.
 *               Change sp_chisq to always start with the observed
 *                image positions and then to use Newton's method
 *                to calculate model image positions.
 *               Split calc_chisq out from sp_chisq.
 *               Gave sp_chisq a verbose option
 * v14Nov97 CDF, Moved plot_boxes call from ray_trace into plot_modim
 *                function.
 *               Fixed calculation of step sizes computed in newton_step.
 * v17Nov97 CDF, Corrected RMS calculation in calc_chisq
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "list_tools.h"
#include "potentials.h"
#include "plotfuncs.h"
#include "fitting.h"
#include "cosmo.h"
#include "lensprogs.h"

#define IMPOSERR 1.0e-4
#define CERR 0.01
#define MAXITER 500
#define RAD2ASEC 206265.0  /* To convert radians to arcsec */
#define TINY 1.0e-18       /* To avoid division by zero */
#define MINCHFLG 1
#define MAXCHFLG 3
#define SAMPLES 5.0
#define FTOL 1.0e-7        /* Tolerance for downhill simplex fitting */
#define NDX 1              /* Grid points on either side of x_0 */
#define DEFBEAM 0.25       /* Default beam size (arcsec) */
#define DEFFERR 0.4        /* Default flux error (mJy) */
#define SCALEUP 2          /* Factor by which to make grids finer */

/*.......................................................................
 *
 * Function get_fitinfo
 *
 * Queries the user to get rms error values on the fluxes and image
 *  positions, as well as getting the method for computing chisq.
 *
 * Input:  Modinfo *modinfo    model and image information
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int get_fitinfo(Modinfo *modinfo)
{
  char line[MAXC];   /* String variable used for getting input */

  /* 
   * Set pixel scale and errors
   */

  modinfo->poserr = DEFBEAM;
  printf("Enter the beam size in arcsec [%4.2f]:  ",
	 modinfo->poserr);
  gets(line);
  if(strcmp(line,"") != 0) {
    while(sscanf(line,"%f",&modinfo->poserr) != 1 || modinfo->poserr < 0.0) {
      fprintf(stderr,"ERROR: Input is not valid\n");
      fprintf(stderr,"Enter the rms error in image positions in arcsec:  ");
      gets(line);
    }
  }
  modinfo->pixscale = modinfo->poserr / SAMPLES;

  modinfo->ferr = DEFFERR;
  printf("Enter the rms error in image fluxes in mJy [%4.2f]:  ",
	 modinfo->ferr);
  gets(line);
  if(strcmp(line,"") != 0) {
    while(sscanf(line,"%f",&modinfo->ferr) != 1 || modinfo->ferr < 0.0) {
      fprintf(stderr,"ERROR: Input is not valid\n");
      fprintf(stderr,"Enter the rms error in image flux in mJy:  ");
      gets(line);
    }
  }

  /*
   * Determine type of chisq calculation to perform
   */

  modinfo->chisqflag = get_chitype();

  return 0;
}

/*......................................................................
 *
 * Function mod_grid_search
 *
 * Does a rough brute-force search for the best fit model parameters
 *  by computing the best chisq for given image positions.
 *
 * Inputs: float *modpars      varying model parameters
 *         float *step         step sizes for each grid dimension
 *         Modinfo *modinfo    non-varying model information
 *         int verbose         if verbose=1 print out progress information
 *
 * Output: int (0 or 1)        0 ==> successful run, 1 ==> error
 */

int mod_grid_search(float *modpars, float *step, Modinfo *modinfo, 
		    int verbose)
{
  int i;               /* Looping variable */
  int no_error=1;      /* Flag set to 0 on error */
  int niter=0;         /* Number of iterations taken in loop */
  float chimin;        /* Min. value of chisq found by the grid search */
  float *minpars=NULL; /* Grid point giving min. chisq */
  float *tmppars=NULL; /* Dummy variable used by grid_iter */
  float *fptr,*fptr2;  /* Pointers to navigate minpars and modpars */
  void *info;          /* Void pointer to modinfo */

  printf("\nDoing a grid search ... Please be patient.\n\n");

  /*
   * Allocate arrays
   */

  if(!(minpars = new_array(modinfo->nvar,1))) {
    fprintf(stderr,"Error: mod_grid_search\n");
    return 1;
  }

  if(!(tmppars = new_array(modinfo->nvar,(modinfo->nvar)+1)))
    no_error = 0;

  /*
   * Calculate chisq at central grid point and put it in chimin, also
   *  put central values of the model parameters into minpars
   */

  if(no_error) {
    info = (void *) modinfo;
    if((chimin = mod_chisq(modpars,info)) < 0)
      no_error = 0;
    else {
      printf("\n** Chisq at grid center = %e\n",chimin);
      for(i=0,fptr=minpars,fptr2=modpars; i<modinfo->nvar; 
	  i++,fptr++,fptr2++) {
	*fptr = *fptr2;
      }
    }
  }

  /*
   * Run the grid
   */

  if(no_error) {
    i = 1;
    if(grid_iter(&i,modinfo->nvar,modpars,tmppars,step,run_amoeba,
		 &chimin,minpars,info,&niter))
      no_error = 0;
    else
      printf("\nGrid search produces best chisq = %f at\n",chimin);
  }


  /*
   * Print out best-fit values and put them into modpars array
   */

  if(no_error) {
    for(i=0,fptr=minpars,fptr2=modpars; i<modinfo->nvar; 
	i++,fptr++,fptr2++) {
      printf("  %f\n",*fptr);
      *fptr2 = *fptr;
    }
  }

  /*
   * Clean up
   */

  minpars = del_array(minpars);
  tmppars = del_array(tmppars);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"Error: mod_grid_search\n");
    return 1;
  }
}

/*......................................................................
 *
 * Function grid_iter
 *
 * Calculates a value of chisq (using the passed function 'func') at
 *  each point of a grid.  The minimum value of chisq is stored in
 *  the variable chisq, and the parameter values associated with the
 *  minimum-chisq model are stored in holdpars.
 * NB: Some value of chisq must be passed to the function when it is
 *  first called.  The best thing to do is to calculate chisq at the
 *  central point of the grid (x0) before calling the function.
 *
 * Inputs: int *n              number of axes completed -- modified
 *                              by the function
 *         int naxes           number of axes
 *         float *x0           point at the center of the grid
 *         float *x            location of the current grid point
 *         float *dx           array of step sizes for grid
 *         float (*func)       function used to evaluate chisq
 *         float *chimin       minimum value of chisq computed
 *         float *xmin         grid point associated with min chisq
 *         void *info          other model info not contained in x
 *         int *niter          number of iterations
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int grid_iter(int *n, int naxes, float *x0, float *x, float *dx, 
	      float (*func)(float *, float *, void *), float *chimin, 
	      float *xmin, void *info, int *niter)
{
  int i,j;        /* Looping variables */
  int step;       /* Number used to determine array members used */
  float tmpchisq; /* Calculated value of chisq */
  float *fptr;    /* Pointer to navigate xmin */
  float *fptr2;   /* Pointer to navigate x */

  for(i=-NDX; i<=NDX; i++) {
    step = *n - 1;
    if (*n == naxes) {
      (*niter)++;
      if(*niter % 100 == 0)
	printf("At iteration %d\n",*niter);
      *(x+step) = *(x0+step) + *(dx+step) * i;
      if((tmpchisq =(*func)(x,dx,info)) < 0) {
	fprintf(stderr,"ERROR: grid_iter\n");
	return 1;
      }
      fptr = x;

      if(tmpchisq < *chimin) {
	*chimin = tmpchisq;
	printf("** Chisq now = %e\n",tmpchisq);
	for(j=0,fptr=xmin,fptr2=x; j<naxes; j++,fptr++,fptr2++) {
	  *fptr = *fptr2;
	  printf("  %f\n",*fptr);
	}
      }
    }
    else if (*n < naxes) {
      *(x+step) = *(x0+step) + *(dx+step) * i;
      (*n)++;
      if(grid_iter(n,naxes,x0,x,dx,func,chimin,xmin,info,niter)) {
	fprintf(stderr,"ERROR: grid_iter\n");
	return 1;
      }
    }
  }

  (*n)--;

  return 0;
}

/*.......................................................................
 *
 * Function run_amoeba
 *
 * This function, called by grid_iter, sets up a simplex and then runs
 *  the downhill simplex process.  The simplex is set up by calling
 *  the function simplex_setup, and the downhill simplex is set up
 *  by calling the function amoeba.
 *
 * Inputs: float *x            starting vertex for the simplex
 *         float *dx           step sizes for the simplex
 *         void *info          pointer to the Modinfo structure containing
 *                              the rest of the model info.
 *
 * Output: float chisq         final value of chisq found by the
 *                              the downhill simplex.  -999.9 on error.
 *
 */

float run_amoeba(float *x, float *dx, void *info)
{
  int niter;         /* Number of simplex iterations */
  int verbose=0;     /* If verbose = 1 then print output */
  float *chisqlist;  /* List of chisq values associated with the simplex */
  float chisq;       /* Min. chisq from fitting */
  Modinfo *modinfo;  /* Un-voided version of info */

  /*
   * Cast info back to a modinfo structure 
   */

  modinfo = (Modinfo *) info;

  /*
   * Set up simplex
   */

  if(!(chisqlist = simplex_setup(x,modinfo->nvar,dx,mod_chisq,info,
				 verbose))) {
    fprintf(stderr,"ERROR: run_amoeba\n");
    return -999.9;
  }

  /*
   * Do the simplex search through model space.
   */

  amoeba(x,chisqlist,modinfo->nvar,FTOL,mod_chisq,&niter,info,verbose);

  chisq = chisqlist[0];

  chisqlist = del_array(chisqlist);

  return chisq;
}

/*......................................................................
 *
 * Function mod_chisq
 *
 * Takes input from model-fitting downhill simplex search and calculates
 *   the value of chisq associated with that particular model.  This
 *   function uses the find_source function to calculate chisq.
 *
 * Inputs:  float *modpars     model parameters
 *          void *minfo        container for the non-variable
 *                             model information
 *			       
 * Output:  float chisq        value of chisq for this model
 *
 */

float mod_chisq(float *modpars, void *minfo)
{
  float chisq;             /* Value of chisq calculated for this model */
  Pos sp={0,0};            /* Source position */
  Modinfo *modinfo=NULL;   /* Model info contained in minfo */

  /*
   * Recast void pointer to Modinfo
   */

  modinfo = (Modinfo *) minfo;

  /* 
   * Get information from fitting routine into Modinfo 
   */

  if(modpars_to_modinfo(modpars,modinfo)) {
    fprintf(stderr,"ERROR:  mod_chisq\n");
    return -999.0;
  }

  sp.x = modinfo->fsource->x.par;
  sp.y = modinfo->fsource->y.par;

  if((chisq = sp_chisq(sp,modinfo,0)) < 0) {
    fprintf(stderr,"ERROR: mod_chisq\n");
    return -999.0;
  }

  return chisq;
}

/*.......................................................................
 *
 * Function source_guess
 *
 * In the case of no source input file, calculate a source position and
 *  transfer it to the modinfo->fsource container
 *
 * Inputs: Modinfo *modinfo    model and image information
 *         float *nzfact       redshift factors
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int source_guess(Modinfo *modinfo, float *nzfact)
{
  int i,j;          /* Looping variables */
  Pos *sp=NULL;     /* Source position returned from calc_sp */
  Pos *xptr;        /* Pointer used to navigate sp */
  Posinfo *pptr;    /* Pointer used to navigate modinfo->posobs */
  Srcinfo *srcptr;  /* Pointer used to navigate modinfo->fsource */

  /*
   * Calculate sp by calling calc_sp
   */

  if(!(sp = calc_sp(modinfo,nzfact)))  {
    fprintf(stderr,"ERROR: source_guess\n");
    return 1;
  }

  /*
   * Transfer sp into modinfo->fsource and into modinfo->posobs
   */

  for(i=0,xptr=sp,srcptr=modinfo->fsource; i<modinfo->nfsrc; 
      i++,xptr++,srcptr++) {
    srcptr->x.par = xptr->x;
    srcptr->y.par = xptr->y;
    for(j=0,pptr=modinfo->posobs; i<modinfo->nobs; i++,pptr++)
      if(pptr->zf == srcptr->zf.par)
	pptr->spos = *xptr;
  }

  /*
   * Clean up
   */

  sp = del_pos(sp);

  return 0;
}


/*.......................................................................
 *
 * Function calc_sp
 *
 * Takes the observed image positions and fluxes and passes rays back
 *  through the model potential at those positions to find the source
 *  position (sp)
 *
 * Inputs: Modinfo *modinfo    model and image information
 *         float *nzfact       redshift factors
 *
 * Output: Pos *sp             the list of source positions
 */

Pos *calc_sp(Modinfo *modinfo, float *nzfact)
{
  int i,j;               /* Looping variables */
  int count;             /* Keeps track of number of images per source */
  Pos *sp;               /* Source position(s) (returned by function) */
  Pos *xptr;             /* Pointer for navigating sp */
  Posinfo *pptr1,*pptr2; /* Pointers for navigating Posinfo arrays */

  /*
   * Allocate memory for sp and pitmp
   */

  if(!(sp=new_pos(modinfo->nfsrc,1))) {
    fprintf(stderr,"ERROR: calc_sp\n");
    return NULL;
  }

  /*
   * Calculate the source position for each redshift factor
   */

  for(i=0,xptr=sp; i<modinfo->nfsrc; i++,xptr++) {
    count=0;

    /*
     * Count the number of images associated with this source
     */

    for(j=0,pptr1=modinfo->posobs; j<modinfo->nobs; j++,pptr1++)
      if (pptr1->zf == *(nzfact+i)) {
	count++;
      }

    /*
     * Now allocate memory for modinfo->posimmod (for this source) and 
     *  fill it with the appropriate position information
     */

    if(!(modinfo->posimmod = new_posinfo(count,1))) {
      fprintf(stderr, "ERROR: calc_sp\n");
      return NULL;
    }

    pptr2 = modinfo->posimmod;
    for(j=0,pptr1=modinfo->posobs; j<modinfo->nobs; j++,pptr1++)
      if (pptr1->zf == *(nzfact+i)) {
	*pptr2++ = *pptr1;
      }

    modinfo->nimmod = count;

    /*
     * Find the source position
     */

    if(find_source(xptr,modinfo,1)) {
      fprintf(stderr,"ERROR: calc_sp\n");
      return NULL;
    }

    /*
     * Free up pimod
     */

    modinfo->posimmod = del_posinfo(modinfo->posimmod);
  }

  return sp;
}

/*.......................................................................
 *
 * Function find_source
 *
 * Locates the source position for a given model potential by shooting
 *   rays at the observed image positions back through the potential and
 *   finding the average "hit" position on the source plane
 *
 * Inputs:  Pos source         source position (modified by the function)
 *          Modinfo *info      container for the image and
 *                               model information
 *          int verbose        if 0 supress output, otherwise print
 *
 * Output:  int (0 or 1)       0 ==> successful, 1 ==> error
 *
 */

int find_source(Pos *sp, Modinfo *info, int verbose)
{
  Pos source;

  /* 
   * Find rough source position by averaging positions found by
   *   passing rays back through the potential at the image
   *   positions.
   */

  if(rough_source(&source,info->posimmod,info->nimmod,info->lp,info->nmod,1)) {
    fprintf(stderr,"ERROR: find_source\n");
    return 1;
  }

  /* 
   * Get a better source position by minimizing chisq found by
   *   passing rays from source back through the potential
   *   and comparing the positions to the observed image positions.
   */

  if(better_source(&source,info,verbose)) {
    fprintf(stderr,"ERROR: find_source\n");
    return 1;
  }

  if(verbose)
    printf("\nBest source position is %7.4f %7.4f\n",source.x,source.y);

  *sp = source;

  return 0;
}

/*......................................................................
 *
 * Function rough_source
 *
 * Locates the source position for a given model potential by shooting
 *   rays at the observed image positions back through the potential and
 *   finding the average "hit" position on the source plane
 *
 * Inputs:  Pos source         source position (modified by this function)
 *          Posinfo *posinfo   list of image positions
 *          int nobs           number of image positions
 *          Lensparams *lp     lens potential  parameters
 *          int nmod           number of model types included in
 *                              the lens potential
 *          int verbose        if 0 supress output, otherwise print
 *
 * Output:  int (0 or 1)       0 ==> successful, 1 ==> error
 *
 */

int rough_source(Pos *source, Posinfo *posinfo, int nobs, Lensparams *lp, 
		 int nmod, int verbose)
{
  int i;              /* Looping variable */
  float sumx,sumy;    /* Used in computing average source position */
  Posinfo *piptr;     /* Pointer to navigate posinfo */
  Poten *spoten=NULL; /* Array of potential values at the image positons */
  Poten *potptr;      /* Pointer to navigate spoten */
  Lensparams *parptr; /* Pointer to navigate lp */

  printf("\nDoing grid search for source position...");

  /*
   * Find the potential at each of the image positions and put the results
   *  in the spoten array.
   */

  if(!(spoten = new_poten(nobs,1))) {
    fprintf(stderr,"\nERROR: rough_source\n");
    return 1;
  }

  for(i=0,parptr=lp; i<nmod; i++,parptr++)
    poten_table[parptr->pottyp].fn(spoten,posinfo,parptr,nobs,1,1);

  /* 
   * Find rough source position by averaging positions found by
   *   passing rays back through the potential at the image
   *   positions.
   */

  sumx = sumy = 0.0;
  for(i=0,potptr=spoten,piptr=posinfo; i<nobs; i++,piptr++,potptr++) {
    if (piptr->zf != posinfo->zf) {
      fprintf(stderr,"\nERROR:  sending multiple sources to rough_source\n");
      spoten = del_poten(spoten);
      return 1;
    }
    sumx += (piptr->x - (piptr->zf * potptr->dphdx));
    sumy += (piptr->y - (piptr->zf * potptr->dphdy));
  }
  source->x = sumx/nobs;
  source->y = sumy/nobs;

  if(verbose)
    printf("\nSource position from rough_source is %7.4f %7.4f\n",
	   source->x,source->y);

  spoten = del_poten(spoten);

  return 0;
}

/*......................................................................
 *
 * Function better_source
 *
 * Given a rough source position, this function finds a better position
 *   by shifting the source position to minimize the chisq value
 *   calculated by comparing the image positions produced by the source
 *   to the observed image positions.
 *
 */

int better_source(Pos *init, Modinfo *info, int verbose)
{
  int i;                /* Looping variable */
  int no_error=1;       /* Flag set to 0 on error */
  int nsteps;           /* Number of steps taken in fitting */
  float initchisq=0.0;  /* Chisq for initial guess */
  float *sparray=NULL;  /* Array of simplex points used in the fitting */
  float *simpstep=NULL; /* Simplex step sizes */
  float *chisq=NULL;    /* Chisq values associated with simplex points */
  Posinfo *tmpinfo;     /* Pointer to navigate info->posimmod */
  List calclist;        /* List of image positions calculated using model */
  List tmplist;         /* Pointer to navigate calclist */

  /*
   * Put initial guess of source position into the source position
   *  container of the posimmod array
   */

  for(i=0,tmpinfo=info->posimmod; i<info->nimmod; i++,tmpinfo++)
    tmpinfo->spos = *init;

  /*
   * Calculate the value of chisq produced by using the initial guess of
   *  source position.  The model image positions will be calculated
   *  by starting with the observed image positions and using Newton's 
   *  method.  The value of chisq is then computed by comparing the model 
   *  image positions (and fluxes) to the observed values.
   */

  if((initchisq = sp_chisq(*init,info,verbose)) < 0)
    no_error = 0;

  /*
   * Print out the rough guess, if running in verbose mode
   */

  if (no_error && verbose)
    printf("Initial chisq at %7.4f %7.4f is chisq = %f \n",
	   init->x,init->y,initchisq);

  /*
   * First do a rough grid search of positions around the initial
   *  guess source position to get a better initial guess
   */

  if(no_error)
    if(grid_search(init,&initchisq,info,verbose))
      no_error = 0;

  /*
   * Now do a downhill simplex search
   */

  if(no_error) 
    if(!(sparray = new_array(3,2)))
      no_error = 0;

  if(no_error)
    if(!(simpstep = new_array(2,1)))
      no_error = 0;

  /*
   * Set up simplex, first assigning the proper initial values to
   *  sparray and simpstep
   */

  if(no_error) {
    sparray[0] = init->x;
    sparray[1] = init->y;
    simpstep[0] = simpstep[1] = 5*info->pixscale;

    if(!(chisq = simplex_setup(sparray,2,simpstep,pos_cast,(void *) info,
			       verbose)))
      no_error = 0;
  }

  /*
   * Do the downhill simplex fitting
   */

  if(no_error) {
    if(verbose)
      printf("\nStarting Downhill Simplex\n\n");
    amoeba(sparray,chisq,2,1.0e-7,pos_cast,&nsteps,(void *) info,verbose);
  }

  /*
   * If running in verbose mode, figure out the new image positions
   *  using the best-fit source position
   */

  if (verbose && no_error) {

    /*
     * Put new source info into posimmod
     */

    for(i=0,tmpinfo=info->posimmod; i<info->nimmod; i++,tmpinfo++) {
      tmpinfo->spos.x = sparray[0];
      tmpinfo->spos.y = sparray[1];
    }

    /*
     * Calculate image positions
     */

    if(!(calclist = find_images_newton(info,0)))
      no_error = 0;

    /*
     * Print out results
     */

    if(no_error) {
      printf("Downhill simplex finds at %7.4f %7.4f, chisq = %f \n\n",
	     sparray[0],sparray[1],chisq[0]);
      printf("   x        y     calcx   calcy\n");
      for(i=0,tmpinfo=info->posimmod,tmplist=calclist; i<info->nimmod; 
	  i++,tmpinfo++,tmplist=tmplist->next)
	printf("%7.4f %7.4f %7.4f %7.4f\n",tmpinfo->x,tmpinfo->y,
	       tmplist->impos.x,tmplist->impos.y);
      printf("\n");
    }

    free_list(&calclist);
  }

  /*
   * Assign best-fit source position
   */

  if(no_error) {
    init->x = sparray[0];
    init->y = sparray[1];
  }

  /*
   * Clean up
   */

  chisq = del_array(chisq);
  sparray = del_array(sparray);
  simpstep = del_array(simpstep);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: better_source\n");
    return 1;
  }
}

/*......................................................................
 *
 * Function grid_search
 *
 * Does a rough brute-force search for the source position giving
 *  the best chisq for a given potential and image positions.
 *
 */

int grid_search(Pos *init, float *chisq, Modinfo *info, int verbose)
{
  int i,j;
  int nsteps=10;
  float stepsize;
  float tmpchisq;
  Pos tmppos,holdpos;

  stepsize = info->pixscale/5;

  holdpos = *init;
  for(i=-nsteps;i<nsteps;i++)
    for(j=-nsteps;j<nsteps;j++) {
      tmppos.x = init->x + i*stepsize;
      tmppos.y = init->y + j*stepsize;
      if((tmpchisq = sp_chisq(tmppos,info,0)) < 0) {
	fprintf(stderr,"ERROR: grid_search\n");
	return 1;
      }
      if (tmpchisq < *chisq) {
	*chisq = tmpchisq;
	holdpos = tmppos;
      }
    }

  if (verbose)
    printf("\nGrid search finds at %7.4f %7.4f chisq = %f \n",
	   holdpos.x,holdpos.y,*chisq);
  *init = holdpos;

  return 0;
}


/*......................................................................
 *
 * Function pos_cast
 *
 * Takes the passed values from amoeba, in the case of trying to
 *   find the best source position, casts them to the appropriate
 *   variables and then calls sp_chisq to get a chisq value.
 *
 * Inputs: float *p,        the source position, as a 2x1 float array
 *         void  *info,     all the other model information
 *
 * Output: float chisq,     the value returned from sp_chisq
 *
 */

float pos_cast(float *p, void *info)
{
  float chisq;
  Pos sp;
  Modinfo *trueinfo;

  /*
   * Cast info to the proper form
   */

  trueinfo = (Modinfo *) info;

  /*
   * Get the source position into a position variable
   */

  sp.x = *p;
  sp.y = *(p+1);

  /* Calculate chisq */

  chisq = sp_chisq(sp,trueinfo,0);

  return chisq;
}

/*.......................................................................
 *
 * Function get_chitype
 *
 * Determines the type of chisq calculation to be done:
 *   1. Fitting to image positions only
 *   2. Fitting to flux-weighted image positions
 *   3. Fitting to fluxes and flux-weighted image positions
 *
 * Inputs: none
 *
 * Output: int chisqflag       flag that sets the chisq calculation
 *
 */

int get_chitype()
{
  int chisqflag;    /* Chisq method flag */
  char line[MAXC];  /* General string for reading input */

  printf("\nChoose the type of chisq calculation to perform:\n");
  printf("   1. Fit to image positions only\n");
  printf("   2. Fit to flux-weighted image positions\n");
  printf("   3. Fit to fluxes and flux-weighted image positions\n");
  printf(" Default is case 1:  ");
  gets(line);
  if(strcmp(line,"") == 0)
    chisqflag = 1;
  else {
    sscanf(line,"%d",&chisqflag);
    if(chisqflag < MINCHFLG || chisqflag > MAXCHFLG) {
      fprintf(stderr,"Outside accepted values.  Taking chisqflag = 1\n");
      chisqflag = 1;
    }
  }

  return chisqflag;
}

/*.......................................................................
 *
 * Function sp_chisq
 *
 * Calculates a chisq value for a given potential and a given source
 *  position.  The flag determines the type of chisq determined:
 *    1. Image position only
 *        chisq = sum( (x_mod - x_obs)^2 / poserr^2 )
 *    2. Image positions weighted by fluxes
 *        chisq = sum( f^2 * (x_mod - x_obs)^2 / (poserr^2 * ferr^2)
 *    3. Fluxes + flux-weighted image positions
 *                   ( f^2 * (x_mod - x_obs)^2      (f_mod - f_obs)^2  )
 *        chisq = sum(------------------------- +  ------------------- )
 *                   (   poserr^2 * ferr^2               ferr^2        )
 *
 * Inputs: Pos sp              position of the source
 *         Modinfo *info       model, grid, source and calculation info
 *         int verbose         flag set to 1 for verbose output
 *
 * Output: float chisq         calculated value of chisq
 *
 */

float sp_chisq(Pos sp, Modinfo *info, int verbose)
{
  int i,j;                 /* Looping variables */
  int no_error=1;          /* Flag set to 0 on error */
  int noposimmod=0;        /* Flag set to 1 if no info->posimmod */
  float chisq;             /* Chisq for a given set of model & obs. values */
  float *hess=NULL;        /* Inverse magnification matrix */
  List calclist=NULL;      /* Calculated ist of model image positions */
  List tmplist;            /* Pointer to navigate calclist */
  Posinfo *pimod=NULL;     /* Model image positions */
  Posinfo *tmpinfo;        /* Pointer to navigate Posinfo arrays */
  Posinfo *tmpinfo2;       /* Pointer to navigate Posinfo arrays */
  Srcinfo *srcptr;         /* Pointer to navigate Srcinfo arrays */
  Poten *tmppoten=NULL;    /* Potential at the model image positions */

  /*
   *  Copy over info->posobs into info->posimmod.
   */

  if(!info->posimmod) {
    noposimmod = 1;
    if(!(info->posimmod = new_posinfo(info->nobs,1))) {
      fprintf(stderr, "ERROR: sp_chisq\n");
      return -999.0;
    }

    for(j=0,tmpinfo=info->posobs,tmpinfo2=info->posimmod; j<info->nobs; 
	j++,tmpinfo++,tmpinfo2++)
      *tmpinfo2 = *tmpinfo;

    info->nimmod = info->nobs;
  }

  /*
   * Get the proper source position and flux into the list
   */

  for(i=0,tmpinfo=info->posimmod; i<info->nimmod; i++,tmpinfo++) {
    tmpinfo->spos = sp;
    for(j=0,srcptr=info->fsource; j<info->nfsrc; j++,srcptr++)
      if(tmpinfo->zf == srcptr->zf.par) {
	tmpinfo->sflux = srcptr->flux.par;
	break;
      }
  }

  if(!(calclist = find_images_newton(info,0)))
    no_error = 0;

  /*
   * Check against division by zero
   */

  if(info->poserr == 0.0)
    info->poserr = TINY;
  if(info->ferr == 0.0)
    info->ferr = TINY;
  
  /*
   * If fitting to fluxes, compute the magnifications at the model image
   *    positions
   */

  if(no_error && info->chisqflag == 3) {

    /*
     * Transfer image positions from calclist to pimod, because the
     *  function that calculates potential at image positions expects the
     *  positions to be in a Posinfo array.
     */

    if(!(pimod = new_posinfo(info->nobs,1)))
      no_error = 0;

    if(no_error) {
      tmpinfo = pimod;
      for(tmplist=calclist; tmplist; tmplist=tmplist->next,tmpinfo++)
	*tmpinfo = tmplist->impos;
    }

    /*
     * Calculate the potentials (one at each image position)
     */

    if(no_error)
      if(!(tmppoten = calc_poten(pimod,info->nobs,info)))
	no_error = 0;

    /*
     * Find value of magnification tensor at model image positions
     */

    if(no_error)
      if(!(hess = calc_hess(tmppoten,pimod,info->nobs,1)))
	no_error = 0;
  }

  /*
   * Actually calculate chisq
   */

  chisq = 0.0;
  if(no_error)
    if(calc_chisq(info,calclist,hess,&chisq,verbose))
      no_error = 0;

  /*
   * Clean up and exit
   */

  free_list(&calclist);
  hess = del_array(hess);
  tmppoten = del_poten(tmppoten);
  pimod = del_posinfo(pimod);
  if(noposimmod)
    info->posimmod = del_posinfo(info->posimmod);

  if(no_error)
    return chisq;
  else {
    fprintf(stderr,"ERROR: sp_chisq\n");
    return -999.0;
  }
}

/*.......................................................................
 *
 * Function calc_chisq
 *
 * The function that actually does the chisq calculation.  The data that
 *  go into the calculation are printed out if the verbose flag is set.
 *
 * Inputs: Modinfo *modinfo    model and image info
 *         List calclist       list of model image positions
 *         float *hess         hessian values at image positions
 *         float *chisq        chisq (set by this function)
 *         int verbose         flag set to 1 for verbose output
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int calc_chisq(Modinfo *modinfo, List calclist, float *hess, float *chisq, 
	       int verbose)
{
  int i;             /* Looping variable */
  float err;         /* Error to be used in chisq calculation */
  float mag;         /* Component inverse magnification */
  float tmpchisq;    /* Chisq associated with individual position */
  double delx,dely;  /* Position difference in x and y */
  double del2;       /* Position difference squared (dx^2+dy^2) */
  double rms=0.0;    /* RMS error in position differences */
  double delf;       /* Flux difference */
  Posinfo *tmpinfo;  /* Pointer for navigating modinfo-posobs */
  List tmplist;      /* Pointer for navigating calclist */

  /*
   * Print out table header if running in verbose mode
   */

  if(verbose) {
    printf("\n   x        y     modx   mody    |del_x|   del_f  err_x ");
    printf("  err_f  chisq\n");
    printf("------- -------- ------ ------- --------- ------- -------");
    printf(" ------ --------\n");
  }

  /*
   * Loop through the images computing chisq for each one
   */

  for(i=0,tmpinfo=modinfo->posobs,tmplist=calclist; i<modinfo->nobs; 
      i++,tmpinfo++,tmplist=tmplist->next) {

    /*
     * Calculate position differences
     */

    delx = tmpinfo->x - tmplist->impos.x;
    dely = tmpinfo->y - tmplist->impos.y;
    del2 = delx*delx + dely*dely;
    switch(modinfo->chisqflag) {
    case 1:
      err = modinfo->poserr;
      delf = 0.0;
      break;
    case 2:
      err = modinfo->poserr * modinfo->ferr / tmpinfo->flux;
      delf = 0.0;
      break;
    case 3:
      err = modinfo->poserr * modinfo->ferr / tmpinfo->flux;
      mag = *(hess+i);
      delf = (tmpinfo->flux - tmpinfo->sflux/mag);
      break;
    default:
      fprintf(stderr,"ERROR: calc_chisq.  Not a valid option for calculating ");
      fprintf(stderr,"chisq.\n");
      return 1;
    }
    tmpchisq = del2/(err*err) + delf*delf/(modinfo->ferr*modinfo->ferr);
    *chisq += tmpchisq;
    rms += del2;
    if(verbose)
      printf("%7.4f %7.4f %7.4f %7.4f %9.5f %7.2f %7.5f %6.4f %f\n",tmpinfo->x,
	     tmpinfo->y,tmplist->impos.x,tmplist->impos.y,sqrt(del2),
	     delf,err,modinfo->ferr,tmpchisq);
  }

  if(verbose) {
    printf("\nRMS position error = %9.5f\n",sqrt(rms/(modinfo->nobs - 1)));
    printf("Total chisq = %f\n",*chisq);
  }

  return 0;
}

/*.......................................................................
 *
 * Function make_grid
 *
 * Creates a grid of positions centered on (0,0), having dimensions
 *   (size x size) and having pixel scale pixscale.
 *
 * Inputs: int size            one dimension of the grid (the grid
 *                               will be square)
 *         float pixscale      pixel scale of the grid
 *
 * Output: Posinfo *posinfo    array of positions
 *
 */

Posinfo *make_grid(int size,float pixscale)
{
  int ix,iy;
  Posinfo *posinfo,*tmp;

  if(!(posinfo = new_posinfo(size,size))) {
    fprintf(stderr,"Error: make_grid\n");
    return NULL;
  }
  tmp = posinfo;
  for(iy= 0; iy < size; iy++)
    for(ix = 0; ix < size; ix++) {
      tmp->y = (iy-(size/2.0)+0.5)*pixscale;
      tmp->x = (ix-(size/2.0)+0.5)*pixscale;
      tmp->zf = 1.0;
      tmp->flux = 0.0;
      tmp++;
    }

  return posinfo;
}

/*.......................................................................
 *
 * Function calc_poten
 *
 * Given an array of positions in the Posinfo format, this function
 *  calculates the lensing potential and its derivatives at each position.
 *
 * Inputs: Posinfo *posinfo    array of positions at which potential is
 *                               calculated
 *         int npos            number of image positions
 *         Modinfo *modinfo    image, lens and source info
 *
 * Output: Poten *poten        array of potentials and derivatives
 *
 */

Poten *calc_poten(Posinfo *posinfo, int npos, Modinfo *modinfo)
{
  int i;               /* Looping variable */
  Poten *poten=NULL;   /* Potential at the input positions */
  Lensparams *parptr;  /* Pointer to navigate modinfo->lp */

  /*
   * Allocate array for potentials
   */

  if(!(poten=new_poten(npos,1))) {
    fprintf(stderr,"ERROR: calc_poten\n");
    return NULL;
  }

  /*
   * Calculate potentials at image positions
   */

  for(i=0,parptr=modinfo->lp; i<modinfo->nmod; i++,parptr++) 
    poten_table[parptr->pottyp].fn(poten,posinfo,parptr,npos,1,0);

  return poten;
}

/*.......................................................................
 *
 * Function find_images_grid
 *
 * Finds rough image positions using the pre-calculated grid of potential
 *   values.  Images are formed at positions where a ray shot back
 *   through the potential hits the source position in the source plane.
 *   All images within a given error (poserr - passed to the function)
 *   are considered, but if rough=0 those within a given distance of a 
 *   previously found image are rejected if the ray shot back through 
 *   the potential hits farther from the true source position than the 
 *   ray associated with the previously found image.
 *
 * Inputs:  Modinfo *modinfo   model, source and image info (modified by
 *                              this function (modinfo->posimmod and
 *                              modinfo->nimmod get set))
 *          float *nzfact      list of redshift factors
 *          int size           size of the grid to check (it will
 *                                 have dimensions size x size)
 *          Posinfo *grid      grid of image positions
 *          Poten *poten       previously calculated grid of 
 *                                  potential values
 *          float poserr       effective size of "source hit" region
 *          int rough          if set to 1 accept all rays within err of
 *                                 the source position
 *          int verbose        set to 1 for verbose output
 *
 *
 * Output:  int (0 or 1)       0 ==> success, 1 ==> error
 *
 */

int find_images_grid(Modinfo *modinfo, float *nzfact, int size, Posinfo *grid, 
		     Poten *poten, float poserr, int rough, int verbose)
{
  int i,j;               /* Looping variables */
  int n=0;               /* Number of images found */
  int flag;              /* Set to 1 if neighbor point is on image list */
  float matchlength;     /* Minimum distance for a unique "hit" */
  double betax,betay;    /* Calculated position of source */
  double delx,dely,del;  /* Dist. btwn calculated and actual posn of source */
  float *zptr;           /* Pointer to navigate nzfact */
  List list;             /* List of found images (--> posinfo) */
  List tmplist;          /* Pointer to navigate list */
  Poten *potptr;         /* Pointer to navigate poten */
  Posinfo *gridptr;      /* Pointer to navigate grid */
  Posinfo *tmpinfo;      /* Pointer to navigate posinfo */
  Srcinfo *srcptr;       /* Pointer to navigate modinfo->fsource */

  /*
   * Initialize list
   */

  init_list(&list);

  /*
   * Start finding images
   */

  matchlength = 4.0*modinfo->pixscale;
  for(j=0,zptr=nzfact,srcptr=modinfo->fsource; j<modinfo->nfsrc;
      j++,zptr++,srcptr++) {
    for(i=0,potptr=poten,gridptr=grid; i<size*size; i++,gridptr++,potptr++) {

      /*
       * For each grid point, calculate the associated source position 
       *  assuming that an image is at that point.  Then calculate the
       *  distance (del) between the calculated source position and the
       *  actual source position.  If it is less than poserr,
       *  then store the position in a list (if it meets further criteria).
       */

      betax = gridptr->x - (*zptr * potptr->dphdx);
      betay = gridptr->y - (*zptr * potptr->dphdy);
      delx = betax - srcptr->x.par;
      dely = betay - srcptr->y.par;
      del = sqrt(delx*delx+dely*dely);

      if (del < 0.25)  {

	/*
	 * If rough=1, accept the grid point with no further criteria and
	 *  add it to the list.
	 */

	if(rough) {
	  n++;
	  gridptr->zf = *zptr;
	  gridptr->spos.x = srcptr->x.par;
	  gridptr->spos.y = srcptr->y.par;
	  gridptr->sflux = srcptr->flux.par;
	  append_list(&list,*gridptr,del);
	}

	/*
	 * If rough=0, compare the del_b for the grid point, with the
	 *  del_b calculated for all nearby points (within 
	 *  matchlength) that are already on the list.  If the
	 *  current grid point has a smaller del_b, then replace the
	 *  point on the list with the current point.
	 */

	else {
	  flag=0;
	  for(tmplist=list; tmplist; tmplist=tmplist->next)
	    if(fabs(gridptr->x-tmplist->impos.x)<matchlength && 
	       fabs(gridptr->y-tmplist->impos.y)<matchlength) {
	      flag++;
	      if(del < tmplist->delb) {
		gridptr->zf = *zptr;
		gridptr->spos.x = srcptr->x.par;
		gridptr->spos.y = srcptr->y.par;
		gridptr->sflux = srcptr->flux.par;
		tmplist->impos = *gridptr;
		tmplist->delb = del;
	      }
	    }

	  /* 
	   * If no nearby points are found, then add the grid point 
	   *  to the list.
	   */

	  if(!flag) {
	    n++;
	    gridptr->zf = *zptr;
	    gridptr->spos.x = srcptr->x.par;
	    gridptr->spos.y = srcptr->y.par;
	    gridptr->sflux = srcptr->flux.par;
	    append_list(&list,*gridptr,del);
	  }
	}
      }
    }
  }

  /*
   * Allocate memory for final image array (modinfo->posimmod) and transfer 
   *  info from list to this array
   */

  modinfo->nimmod = n;
  if(!(modinfo->posimmod = new_posinfo(n,1))) {
    fprintf(stderr,"\nERROR: find_images_grid\n");
    free_list(&list);
    return 1;
  }

  for(i=0,tmplist=list,tmpinfo=modinfo->posimmod; i<n;
      i++,tmpinfo++,tmplist=tmplist->next) {
    if(!rough && verbose) {
      printf("Image %d at %7.4f, %7.4f. Delta_b = %7.4f.  ",i+1,
	     tmplist->impos.x,tmplist->impos.y,tmplist->delb);
      printf("Redshift factor = %7.4f.\n",tmplist->impos.zf);
    }
    *tmpinfo = tmplist->impos;
  }

  /*
   * Clean up and exit
   */

  free_list(&list);
  return 0;
}

/*......................................................................
 *
 * Function find_images_newton
 *
 * Uses Newton's method to refine the list of rough image positions to
 *   a list of more accurate image positions.  Newton's method is a
 *   way to find roots of an equation.  In one dimension, f(x) is the
 *   function for which you want to find the root, i.e. you want to
 *   find the value of x that makes f(x) = 0.  You supply an initial
 *   guess x_0 and then each successive approximation for x is derived
 *   from the fact that 
 *
 *       f(x + del_x) = f(x) + f'(x) del_x + ...
 *
 *   implying that
 *
 *       x_(k+1) = x_k - (f(x)/f'(x)).
 *
 *   The extension to two dimensions is described in the comments for
 *   newton_step.
 *
 * Inputs:  Modinfo *info      container holding the model parameters
 *          int verbose        if 0 supress output
 *
 * Output:  List list          list of more accurate image 
 *                               positions and redshift factors
 *
 */

List find_images_newton(Modinfo *info, int verbose)
{
  int i,j;                      /* Looping variables */
  int count=0;                  /* Counts number of iterations */
  float delbx,delby;            /* Diff in computed and actual source pos */
  double dx,dy;                 /* Shifts computed by Netwon's method */
  Posinfo xk;                   /* Position found with Newton's method */
  Posinfo *tmpinfo;             /* Pointer to navigate posinfo */
  Lensparams *parptr;           /* Pointer to navigate info->lp */
  List list;                    /* List of image positions returned by fn */
  Poten tmppot;                 /* Potential at final position */
  Poten potzero={0.0,0.0,0.0,0.0,0.0,0.0};  /* Used to zero tmppot */

  /*
   * Initialize the list
   */

  init_list(&list);

  for(j=0,tmpinfo=info->posimmod; j<info->nimmod; j++,tmpinfo++) {

    /*
     * Set initial guess to image position passed to this function
     *  and set count = 0
     */

    xk = *tmpinfo;
    count = 0;

    /*
     * Compute the step-size to the next guess at the image position
     */

    if(newton_step(xk,info,&dx,&dy,tmpinfo)) {
      fprintf(stderr,"ERROR: find_images_newton\n");
      free_list(&list);
      return NULL;
    }

    /*
     * Loop using Newton's method until the computed stepsize is smaller
     *  than IMPOSERR or MAXITER iterations have been done.
     */

    while (sqrt(dx*dx+dy*dy) > IMPOSERR && count < MAXITER) {

      /*
       * Increment the iteration counter
       */

      count++;

      /*
       * Find the next image position using the computed dx and dy values
       */

      xk.x += dx;
      xk.y += dy;

      /*
       * Compute the step to the next image position
       */

      if(newton_step(xk,info,&dx,&dy,tmpinfo)) {
	fprintf(stderr,"ERROR: find_images_newton\n");
	free_list(&list);
	return NULL;
      }
    }

    /*
     * Compute the source position associated with the given image position
     *  and lens model, and find its distance to the "true" source position
     *  given by the model.
     */

    tmppot = potzero;
    for(i=0,parptr=info->lp; i<info->nmod; i++,parptr++)
      poten_table[parptr->pottyp].fn(&tmppot,&xk,parptr,1,1,0);
    delbx = (xk.x - (tmpinfo->zf * tmppot.dphdx)) - tmpinfo->spos.x;
    delby = (xk.y - (tmpinfo->zf * tmppot.dphdy)) - tmpinfo->spos.y;

    /*
     * Append the position to the list
     */

    if(verbose)
      printf("Image at  %6.4f, %6.4f.  %d iterations.\n",xk.x,xk.y,count);
    append_list(&list,xk,sqrt(delbx*delbx+delby*delby));
  }

  return list;
}

/*.......................................................................
 *
 * Function newton_step
 *
 * Computes the step to the next image using Newton's Method.  In the
 *  two-dimensional case, think of things in terms of vectors and
 *  matrices.  Once again we have
 *
 *    x_(k+1) = x_k + del_x
 *
 *  where
 *
 *    del_x = -J^{-1} . F
 *
 *  In the above equation, F is the vector of functions for which we are
 *  are trying to find roots, and J is the Jacobian matrix,
 *
 *    J_ij = dF_i/dx_j
 *
 *  Since we are dealing with the 2D case, we can easily find the inverse
 *  of the Jacobian matrix by remembering that the inverse of a 2D matrix
 *
 *         ( a  b )                    1    (  d  -b )
 *     A = (      )    is   A^{-1} = -----  (        )
 *         ( c  d )                  det A  ( -c   a )
 *
 *  For the Jacobian matrix,
 *    a = 1 - dphixx
 *    b = c = -dphixy
 *    d = 1 - dphiyy.
 *  So just plug and chug.
 *
 * Inputs: Posinfo xk          current image position
 *         Modinfo *modinfo    model information
 *         double *dx          step size in x direction (set by function)
 *         double *dy          step size in y direction (set by function)
 *         Posinfo *tmpinfo    info associated with current image position
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int newton_step(Posinfo xk, Modinfo *modinfo, double *dx, double *dy,
		Posinfo *tmpinfo)
{
  int i;                        /* Looping variable */
  float dphixx,dphixy,dphiyy;   /* Second derivatives of the potential */
  float delbx,delby;            /* Diff in computed and actual source pos */
  float detj;                   /* Determinant of the Jacobian matrix */
  Lensparams *parptr;           /* Pointer to navigate modinfo->lp */
  Poten tmppot;                 /* Potential calculated each iteration */
  Poten potzero={0.0,0.0,0.0,0.0,0.0,0.0};  /* Used to zero tmppot */

  /*
   * Zero out the potential and compute it at the current position.
   */

  tmppot = potzero;
  for(i=0,parptr=modinfo->lp; i<modinfo->nmod; i++,parptr++)
    poten_table[parptr->pottyp].fn(&tmppot,&xk,parptr,1,1,0);

  /*
   * Compute the redshift-corrected potential derivatives
   */

  dphixx = tmpinfo->zf * tmppot.dphdxx;
  dphiyy = tmpinfo->zf * tmppot.dphdyy;
  dphixy = tmpinfo->zf * tmppot.dphdxy;

  /*
   * Compute the determinant of the Jacobian matrix
   */

  detj = (1-dphixx) * (1 - dphiyy) - (dphixy * dphixy);

  /*
   * Compute the difference between the computed and actual source
   *  positions.
   */

  delbx = xk.x - (tmpinfo->zf * tmppot.dphdx) - tmpinfo->spos.x;
  delby = xk.y - (tmpinfo->zf * tmppot.dphdy) - tmpinfo->spos.y;

  /*
   * Using Newton's method, compute the steps to the next image position
   */

  *dx = -1.0 * ((1 - dphiyy)*delbx + dphixy*delby) / detj;
  *dy = -1.0 * (dphixy*delbx + (1 - dphixx)*delby) / detj;

  return 0;
}

/*.......................................................................
 *
 * Function condense_imlist
 *
 * Takes a list of image positions and eliminates all positions within
 *  one pixel (modinfo->pixscale) of each other.  The positions in the
 *  resulting list are transferred into the modinfo->posimmod array.
 *
 * Inputs: List imlist         list to be pared down
 *         Modinfo *modinfo    source, image and model info
 *         int verbose         flag set to 1 for verbose output
 *
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *          NB: modinfo gets modified by this function
 *
 */

int condense_imlist(List imlist, Modinfo *modinfo, int verbose)
{
  int nim=1;         /* Number of images in pared-down list */
  int flag;          /* Flag set to 1 if duplicate position is found */
  float matchlength; /* Minimum distance for a unique hit */
  List finallist;    /* Pared-down list */
  List tmplist;      /* Pointer to navigate imlist and finallist */
  List tmplist2;     /* Pointer to navigate imlist and finallist */
  Posinfo *piptr;    /* Pointer to navigate modinfo->posimmod */

  printf("\nEliminating duplicates from final image list...\n\n");

  /*
   * Initialize finallist and put first node of imlist into first node
   *  of finallist
   */

  init_list(&finallist);
  append_list(&finallist,imlist->impos,imlist->delb);

  /*
   * Loop through the rest of imlist, comparing image positions to those
   *  already in finallist.  If there are no positions in finallist within
   *  matchlength of the imlist position (flag=0), add the imlist 
   *  position to finallist.  Otherwise (flag=1), compare the two matching 
   *  positions and choose the one that produces a ray closest to the source 
   *  position inthe source plane (by comparing the delb parameter).
   */

  matchlength = 2.0 * modinfo->pixscale;
  for(tmplist=imlist->next; tmplist; tmplist=tmplist->next) {
    flag = 0;
    for(tmplist2=finallist; tmplist2; tmplist2=tmplist2->next)
      if (tmplist->impos.zf == tmplist2->impos.zf &&
	  fabs(tmplist->impos.x - tmplist2->impos.x) < matchlength &&
	  fabs(tmplist->impos.y - tmplist2->impos.y) < matchlength ) {
	flag++;
	if(tmplist->delb < tmplist2->delb) {
	  tmplist2->impos = tmplist->impos;
	  tmplist2->delb = tmplist->delb;
	}
      }
    if(!flag) {
      append_list(&finallist,tmplist->impos,tmplist->delb);
      nim++;
    }
  }

  /*
   * Check number of images found.
   */

  if(nim > modinfo->nimmod) {
    fprintf(stderr,"\nERROR: condense_imlist.  More images found than in ");
    fprintf(stderr,"original list (%d vs. %d)\n",nim,modinfo->nimmod);
    return 1;
  }

  /*
   * Transfer finallist information into modinfo->nimmod and 
   *  modinfo->posimmod
   */

  modinfo->nimmod = nim;
  for(piptr=modinfo->posimmod,tmplist=finallist; tmplist; 
      piptr++,tmplist=tmplist->next) {
    *piptr = tmplist->impos;
    if(verbose)
      printf("  Image at %7.4f, %7.4f\n",piptr->x,piptr->y);
  }

  return 0;
}

/*.......................................................................
 *
 * Function predict_delay
 *
 * Calculates the predicted time delay at each position given to the
 *   function as input.  Calculates the gradient of the potential at
 *   each position and returns an array of time delay values, one  
 *   for each position.
 *
 * Inputs: Modinfo *modinfo    lens model information, including image
 *                              positions
 *         float zl            redshift of the lens
 *         float zs            redshift of the source
 *         float omega0        Omega_0
 *
 * Output: float *delays       list of time delay values
 *
 */

float *predict_delay(Modinfo *modinfo, float zl, float zs, float omega0)
{
  int i,j;               /* Looping variables */
  int no_error=1;        /* Flag set to 0 on error */
  int zflag=0;           /* Flag for number of redshifts known */
  int nim;               /* Number of images */
  float delx,dely,del2;	 
  float k;               /* Scaling factor */
  float *delays=NULL;    /* Time delays between images */
  float *abstime=NULL;   /* Absolute delays associated with images */
  float *tptr;           /* Pointer for navigating delays and abstime */
  Poten *tmppoten=NULL;  /* Potential at image positions */
  Poten *potptr;         /* Pointer for navigating tmppoten */
  Posinfo *piptr;        /* Pointer for navigating modinfo->posimmod */

  /*
   * Allocate memory for arrays
   */

  nim = modinfo->nimmod;

  if(!(delays = new_array(nim*(nim-1)/2,1))) {
    fprintf(stderr,"ERROR: predict_delay\n");
    return NULL;
  }

  if(!(abstime = new_array(nim,1)))
    no_error = 0;

  /*
   * Calculate potential at image positions
   */

  if(no_error)
    if(!(tmppoten = calc_poten(modinfo->posimmod,nim,modinfo)))
      no_error = 0;

  /*
   * Calculate the absolute time delay at each of the image positions
   *  using angular distances.  This uses the "scaled" version of the
   *  time delay equation:
   *
   *       Delta_t = T_0 * [ 1/2 * (theta - beta)^2 - phi]
   */

  for(i=0,potptr=tmppoten,piptr=modinfo->posimmod,tptr=abstime; i<nim; 
      i++,piptr++,tptr++,potptr++) {

    /*
     * Do the calculations
     */

    delx = piptr->x - modinfo->fsource->x.par;
    dely = piptr->y - modinfo->fsource->y.par;
    del2 = delx*delx + dely*dely;
    *tptr = 0.5*del2 - ((piptr->zf) * potptr->ph);

    /*
     * Print out the results of the calculations for error-checking
     *  purposes
     */

    printf("Pot = %f, dpdx = %f, dpdy = %f\n",potptr->ph,potptr->dphdx,
	   potptr->dphdy);
    printf("Image %d -- Geom part = %g, Grav part = %g\n",
	   i+1,0.5*del2, -(piptr->zf) * potptr->ph);
    printf("x - beta = %f, del(phi) = %f\n",sqrt(del2),
	   sqrt(potptr->dphdx*potptr->dphdx + potptr->dphdy*potptr->dphdy));
  }

  /* 
   * Calculating the absolute delays in physical units.
   *
   *  It should be noted that del_t = T_0 * (0.5(x - beta)^2 - phi)
   *  where
   *         T_0 = (1+z_l) * D_l * D_s / (c * D_ls) / (206265)^2
   *
   *             = (1+z_l) * D_l * D_s * (sec/3e10cm) * (3e27cm/Gpc) 
   *               ---------------------------------------------------
   *                       (206265)^2 * D_ls
   *
   *             = 2.4e6 * (1+z_l) * (D / 1Gpc) sec
   *                     where D = D_l * D_s / D_ls
   *
   *             = 28 * (1+z_l) * (D / 1Gpc) days
   *
   *  or
   *         T_0 = 7.8e-22 * (1+z_l) * (D / 1cm) sec
   *
   */

  if(zl > 0.0) {
    if(zs > 0.0)
      zflag = 3;
    else
      zflag = 1;
  }
  else if(zs > 0.0)
    zflag = 2;

  printf("\n");
  switch(zflag) {
  case 0:
    k = 28;
    printf("T_0 = 28 * (1+z_l) * (D/1Gpc) days\n\n");
    printf("Multiply all numerical factors below by K = (1+z_l)*(D/1Gpc)\n");
    printf("  where D = D_l * D_s / D_ls.\n\n");
    printf("NB: the factor of h^{-1} is contained in D.\n\n");
    break;
  case 1:
    k = (1/C) * (1+zl) * ang_dist(0.0,zl,omega0,0.0) / (RAD2ASEC*RAD2ASEC);
    printf("T_0 = %6.4f * D_s / D_ls h^{-1} days\n\n",k/3600.0/24.0);
    printf("D_l = %6.1f h^{-1} Mpc\n\n",ang_dist(0.0,zl,omega0,0.0)/MPC2CM);
    printf("Multiply all numerical factors below by K = h^{-1} D_s / D_ls.\n\n");
    break;
  case 2:
    k = (1/C) * ang_dist(0.0,zs,omega0,0.0) / (RAD2ASEC*RAD2ASEC);
    printf("T_0 = %6.4f * (1+z_l) * D_l / D_ls h^{-1} days\n\n",k/3600.0/24.0);
    printf("D_s = %6.1f h^{-1} Mpc\n\n",ang_dist(0.0,zs,omega0,0.0)/MPC2CM);
    printf("Multiply all numerical factors below by ");
    printf("K = h^{-1} (1+z_l) D_l / D_ls.\n\n");
    break;
  case 3:
    k = (1/C) * (1+zl) * ang_dist(0.0,zl,omega0,0.0) * ang_dist(0.0,zs,omega0,0.0)
      / (ang_dist(zl,zs,omega0,0.0) * RAD2ASEC * RAD2ASEC);
    printf("T_0 = %6.4f h^{-1} days\n\n",k/3600.0/24.0);
    printf("D_l  = %6.1f h^{-1} Mpc\n",ang_dist(0.0,zl,omega0,0.0)/MPC2CM);
    printf("D_s  = %6.1f h^{-1} Mpc\n",ang_dist(0.0,zs,omega0,0.0)/MPC2CM);
    printf("D_ls = %6.1f h^{-1} Mpc\n\n",ang_dist(zl,zs,omega0,0.0)/MPC2CM);
    printf("Multiply all numerical factors below by K = h^{-1}.\n\n");
    break;
  default:
    k = 1.0e17;
    printf("Multiply all numerical factors below by K = (1+z_l)*(D/1Gpc)\n");
    printf("  where D = D_l * D_s / D_ls.\n\n");
    printf("NB: the factor of h^{-1} is contained in D.\n\n");
  }

  /*
   * Calculate relative time delays between components
   */

  if(no_error) {
    for(i=0,tptr=delays; i<nim; i++)
      for(j=i+1; j<nim; j++) {
	tptr++;
	*tptr = *(abstime+j) - *(abstime+i);
	printf("Time delay: delta_t_%d - delta_t_%d = ",j+1,i+1);
	printf("%6.2f*T0  = ",*tptr);
	printf("%6.2f*K days\n",*tptr * k / 3600.0 / 24.0);
      }
  }

  /*
   * Clean up and exit
   */

  abstime = del_array(abstime);
  tmppoten = del_poten(tmppoten);

  if(no_error)
    return delays;
  else {
    fprintf(stderr,"ERROR: predict_delay\n");
    return NULL;
  }
}

/*......................................................................
 *
 * Function time_delay
 *
 * Calculates the predicted time delay at each position given to the
 *   function as input.  This particular function is a specialized
 *   case of the function predict_delay in that it takes as input the
 *   grid of input positions and already calculated potentials.  Returns 
 *   an array of time delays, one for each position, which will be used
 *   to plot the time delay contours on the final output plot.
 *
 * Inputs:  Posinfo *grid      coordinate grid
 *          Poten *poten       previously calculated grid of 
 *                               potential values
 *          int size           one dimension of the grid (the grid
 *                              is square)
 *          Modinfo *modinfo   contains number of sources and source
 *                              positions
 *
 * Output:  float *time        grid of time delay values
 *
 */

float *time_delay(Posinfo *grid, Poten *poten, int size, Modinfo *modinfo)
{
  int i,j;           /* Looping variables */
  float *time=NULL;  /* Array of arrival times */
  float *tptr;       /* Pointer to navigate time array */
  float dx,dy;       /* Differences between (x,y) and (beta_x,beta_y) */
  Posinfo *gridptr;  /* Pointer to navigate grid */
  Poten *potptr;     /* Pointer to navigate poten */
  Srcinfo *srcptr;   /* Pointer to navigate modinfo->fsource */

  if(!(time = new_array(size*size,modinfo->nfsrc))) {
    fprintf(stderr,"Error: time_delay\n");
    return NULL;
  }

  tptr = time;
  for(j=0,srcptr=modinfo->fsource; j<modinfo->nfsrc; j++,srcptr++) {
    for(i=0,gridptr=grid,potptr=poten; i<size*size; 
	i++,gridptr++,potptr++,tptr++) {
      dx = gridptr->x - srcptr->x.par;
      dy = gridptr->y - srcptr->y.par;
      *tptr = 0.5 * ((dx*dx) + (dy*dy)) - potptr->ph;
    }
  }

  return time;
}

/*......................................................................
 *
 * Function calc_hess
 *
 * Calculates the Hessian for a potential.  The Hessian is defined
 *   as the determinant of the Jacobian matrix | d_beta / d_theta |, ie.
 *
 *              |                                  |
 *              |  1 - d2phi/dx2      -d2phi/dxdy  |
 *              |                                  |
 *              |   -d2phi/dxdy      1 - d2phi/dy2 |
 *              |                                  |
 *
 *   The Hessian is important because its inverse gives the magnification
 *   matrix.
 *
 * Inputs:  Poten *impot       potential at each of the image positions
 *          Posinfo *pinfo     list of image positions and redhshift factors
 *          int nim            number of image positions
 *          int flag           if flag=0 apply a single redshift factor
 *                               to every value of the potential
 *                             if flag=1 increment the redshift factor
 *                               array along with the potential array
 *
 * Output:  float *magptr      list of Hessian values
 *
 */

float *calc_hess(Poten *impot, Posinfo *pinfo, int nim, int flag)
{
  int i;
  float *magptr,*tmpptr;
  float hxx,hxy,hyy;             /* Components of the Hessian matrix */
  Posinfo *piptr;
  Poten *impotptr;

  if(!(magptr=new_array(nim,1))) {
    fprintf(stderr,"Error: calc_hess\n");
    return NULL;
  }

  impotptr = impot;
  tmpptr = magptr;
  piptr = pinfo;
  for(i=0;i<nim;i++,tmpptr++,impotptr++) {
    hxx = 1 - (piptr->zf * impotptr->dphdxx);
    hyy = 1 - (piptr->zf * impotptr->dphdyy);
    hxy = -(piptr->zf) * impotptr->dphdxy;
    *tmpptr = hxx*hyy - hxy*hxy;
    if(flag)
      piptr++;
  }

  return magptr;
}
      

/*......................................................................
 *
 * Function find_caustics
 *
 * Finds critical curves and caustics by looking for points where the
 *   Hessian matrix is zero ==> infinite magnifications.  The function
 *   uses the pre-calculated potential.  Plots these curves.
 *
 * Inputs:  Posinfo *grid      coordinate grid
 *          Poten *poten       pre-calculated grid of potential values
 *          int size           dimensions of the grid
 *          Modinfo *modinfo   model, source and image info
 *
 * Output:  int = 0,           for successful completion
 *
 */

int find_caustics(Posinfo *grid, Poten *poten, int size, Modinfo *modinfo)
{
  int i,j;             /* Looping variables */
  int no_error=1;      /* Flag set to 0 on error */
  int count;           /* Number of critical curve points found */
  float *hess=NULL;    /* Hessian matrix */
  float *hessptr;      /* Pointer to navigate hess */
  Posinfo *gptr;       /* Pointer to navigate grid */
  Posinfo *crit=NULL;  /* Array of critical curve points */
  Posinfo *critptr;    /* Pointer to navigate crit */
  Posinfo *caus=NULL;  /* Array of caustic curve points */
  Posinfo *causptr;    /* Pointer to navigate caus */
  List list;           /* List of critical curve points */
  List tmplist;        /* Pointer to navigate list */
  Poten *tmppot=NULL;  /* Potential at each critical curve point */
  Poten *tmpptr;       /* Pointer to navigate tmppot */
  Srcinfo *srcptr;     /* Pointer to modinfo->fsource */

  /*
   * Calculate the hessian associated with this potential
   */

  hess = calc_hess(poten,grid,size*size,0);

  /*
   * Loop through the list of sources, calculating critical curves and
   *  caustics for each one
   */

  for(j=0,srcptr=modinfo->fsource; j<modinfo->nfsrc; j++,srcptr++) {

    /*
     * Initialize the list that will contain the critical curve points
     */

    init_list(&list);

    /*
     * Find all points on the critical curves, ie. points for which the
     *  Hessian value is within CERR of 0.  Put these values into a
     *  linked list and count the number of values.
     */

    count=0;
    for(i=0,hessptr=hess,gptr=grid; i<size*size; i++,hessptr++,gptr++) {
      if (fabs(*hessptr) < CERR) {
	count++;
	append_list(&list,*gptr,0.0);
      }
    }
    printf("%d critical curve points found.\n",count);
    if(count == 0) {
      fprintf(stderr,"ERROR: find_caustics -- No critical points found\n\n");
      no_error = 0;
    }

    /*
     * Transfer the critical curve points from the linked list to a
     *  posinfo array and plot the points
     */

    if(no_error)
      if(!(crit=new_posinfo(count,1)))
	no_error = 0;

    if(no_error) {
      for(critptr=crit,tmplist=list; tmplist->next; 
	  tmplist=tmplist->next,critptr++)
	*critptr = tmplist->impos;

      /*
       * Plot the critical curves
       */

      plot_cont(hess,size,0.0,0.0,1,3);
    }
    
    /*
     * Calculate the value of the potential at the critical curve points
     */

    if(no_error)
      if(!(tmppot = calc_poten(crit,count,modinfo)))
	no_error = 0;

    /*
     * Calculate the caustic points associated with the critical curve
     *  points  and plot them.
     */

    if(no_error)
      if(!(caus = new_posinfo(count,1)))
	no_error = 0;

    if(no_error) {
      for(i=0,causptr=caus,critptr=crit,tmpptr=tmppot; i<count; 
	  i++,critptr++,tmpptr++,causptr++) {
	causptr->x = (critptr->x - tmpptr->dphdx);
	causptr->y = (critptr->y - tmpptr->dphdy);
      }

      /*
       * Plot the caustic points
       */

      plot_pts(caus,size,modinfo->pixscale,count,2,1);
    }

    /*
     * Clean up from this round
     */

    crit = del_posinfo(crit);
    caus = del_posinfo(caus);
    tmppot = del_poten(tmppot);
    free_list(&list);

    if(no_error)
      break;
  }

  /*
   * Clean up
   */

  hess = del_array(hess);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: find_caustics\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function plot_poten
 *
 * Plots a greyscale or contour plot of the gravitational potential
 *  associated with a lens model.
 *
 * Inputs: Poten *poten        grid of grav. potential
 *         int size            size of (square) grid
 *         Modinfo *modinfo    model info
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_poten(Poten *poten, int size, Modinfo *modinfo)
{
  int i;             /* Looping variable */
  int plotoption=3;  /* Flag that sets what to plot */
  float *phi=NULL;   /* Float version of potential */
  float *fptr;       /* Pointer to navigate phi */
  char line[MAXC];   /* General string for reading input */
  Poten *potptr;     /* Pointer to navigate poten */

  /*
   * Find out what the user wants plotted and set plotoption flag
   *  accordingly (0 ==> no plotting, 1 ==> grey, 2 ==> contour,
   *  3 ==> grey + contour).
   */

  printf("\nPlot greyscale of potential (y/n)? [n]  ");
  gets(line);
  if (strcmp(line,"") == 0 || line[0] == 'n' || line[0] == 'N')
    plotoption -=1;

  printf("Plot contours of potential (y/n)? [n]  ");
  gets(line);
  if (strcmp(line,"") == 0 || line[0] == 'n' || line[0] == 'N')
    plotoption -=2;

  /*
   * If we will be plotting, assign the potential values to phi
   */

  if(plotoption > 0) {
    if(!(phi = new_array(size,size))) {
      fprintf(stderr,"ERROR: plot_poten\n");
      return 1;
    }

    for(i=0,potptr=poten,fptr=phi; i<size*size; i++,potptr++,fptr++)
      *fptr = potptr->ph;
  }

  /*
   * Plot what has been requested
   */

  switch(plotoption) {
  case 0:
    break;
  case 1:
    plot_gray(phi,size);
    break;
  case 2:
    plot_cont(phi,size,0.16*modinfo->lp->b.par,0.0,30,1);
    break;
  case 3:
    plot_gray(phi,size);
    plot_cont(phi,size,0.16*modinfo->lp->b.par,0.0,30,1);
    break;
  default:
    fprintf(stderr,"\nERROR: plot_poten.  No potential being plotted\n");
  }

  /*
   * Clean up and exit
   */

  phi = del_array(phi);

  return 0;
}

/*.......................................................................
 *
 * Function plot_time_delay
 *
 * Plots contours of the time delay surface, if requested.  The function
 *  finds the location of the negative parity image positions (saddle
 *  points).  If there are more than one, the contour spacing is chosen
 *  to pass through at least two of the saddle points.  This makes the
 *  plot look snazzier.
 *
 * Inputs: Posinfo *grid       coordinate grid
 *         Poten *poten        gravitational potential grid
 *         int size            size of grids
 *         Modinfo *modinfo    model, image and source info
 *         float *hess         hessian values at image positions
 *         Poten *impot        grav. potential at image positions
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_time_delay(Posinfo *grid, Poten *poten, int size, Modinfo *modinfo, 
		    float *hess, Poten *impot)
{
  int i;               /* Looping variable */
  int no_error=1;      /* Flag set to 0 on error */
  int nneg=0;          /* Number of negative parity images */
  float *del_t=NULL;   /* Time delay grid */
  float *imdelt=NULL;  /* Time delays at neg. parity image positions */
  float *hessptr;      /* Pointer to navigate hess */
  float *dtptr;        /* Pointer to navigate imdelt */
  float dx,dy;         /* Differences between (x,y) and (beta_x,beta_y) */
  float contzero;      /* Zero point for contouring */
  float contstep;      /* Step size for contouring */
  Poten *potptr;       /* Pointer to navigate impot */
  Posinfo *piptr;      /* Pointer to navigate modinfo->posimmod */

  /*
   * Find time delay contours using the potential grid.
   */

  printf("\nCalculating the time delay surface....\n");
  if(!(del_t = time_delay(grid,poten,size,modinfo))) {
    fprintf(stderr,"\nERROR: plot_time_delay\n");
    return 1;
  }

  /*
   * Count the number of saddle points
   */

  for(i=0,hessptr=hess; i<modinfo->nimmod; i++,hessptr++)
    if (*hessptr < 0)
      nneg++;



  /*
   * If there is at least one saddle point, allocate memory for array
   *  of time delays at the saddle points and calculate the delays.
   */

  if(nneg > 0) {

    /*
     * Allocate memory for array of saddle point delays
     */

    if(!(imdelt = new_array(nneg,1)))
      no_error = 0;

    /*
     * Go through the Hessian and other arrays again, assigning the
     *  calculated value of the delay at each saddle point to imdelt
     */

    if(no_error) {
      dtptr = imdelt;
      for(i=0,piptr=modinfo->posimmod,hessptr=hess,potptr=impot; 
	  i<modinfo->nimmod; i++,piptr++,hessptr++,potptr++)
	if (*hessptr < 0) {

	  /*
	   * Calculate d_t at saddle point and assign value to imdelt
	   */

	  dx = piptr->x - piptr->spos.x;
	  dy = piptr->y - piptr->spos.y;
	  *dtptr++ = 0.5 * (dx*dx + dy*dy) - potptr->ph;
	}
    }

    /*
     * Now choose contour spacing based on number of saddle points
     */

    if(nneg == 1 && no_error) {
      contstep = 0.1*(modinfo->lp->b.par)*(modinfo->lp->b.par);
      contzero = *imdelt - 10*contstep;
    }
    else if (nneg >= 2 && no_error) {
      contstep = fabs(*imdelt - *(imdelt+1))/5;
      contzero = *imdelt - 10*contstep;
    }
  }

  /*
   * If there are no saddle points, base contour spacing on critical radius
   *  of first model in modinfo->lp
   */

  else if (no_error) {
    contstep = 0.1*(modinfo->lp->b.par)*(modinfo->lp->b.par);
    contzero = -999.9;
  }

  plot_cont(del_t,size,contstep,contzero,30,1);

  /*
   * Clean up
   */

  del_t = del_array(del_t);
  imdelt = del_array(imdelt);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"\nERROR: plot_time_delay\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function plot_modim
 *
 * Plots image positions found by using the lens potential and source
 *  position contained in the model.  Also plot observed image positions.
 *
 * Inputs: Modinfo *modinfo    model and image info
 *         Posinfo *grid       coordinate grid
 *         Poten *poten        gravitational potential grid
 *         int size            size of grids
 *         float *nzfact       array of redshift factors
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int plot_modim(Modinfo *modinfo, Posinfo *grid, Poten *poten, int size,
	       float *nzfact)
{
  int i;               /* Looping variable */
  int imtype=1;        /* Choice of images to plot */
  float sourcesize;    /* Size of background source */
  char line[MAXC];     /* General string for reading input */
  Lensparams *parptr;  /* Pointer to navigate modinfo->lp */

  /*
   * Choose whether to plot images corresponding to an unresolved 
   *  background source (choice 1) or a background source with finite
   *  size (choice 2)
   */

  printf("Image plotting choices:\n");
  printf("  1. Best image locations\n");
  printf("  2. Ray-shoot with finite source size\n");
  printf("Enter choice:  [1] ");
  gets(line);
  if(strcmp(line,"") != 0) {
    while(sscanf(line,"%d",&imtype) != 1) {
      fprintf(stderr,"ERROR: input not valid.  Enter choice again\n");
      gets(line);
    }
  }

  /*
   * Now do the plotting based on the choice above
   */
#if 1
  switch(imtype) {
  case 1:

    /*
     * If choosing an unresolved background source, just use image
     *  positions already in modinfo->posimmod
     */

    plot_pts(modinfo->posimmod,size,modinfo->pixscale,modinfo->nimmod,
	      6,18);
    break;

  case 2:

    /*
     * If ray-shooting back to a source of finite size, make the
     *  coordinate grid finer and then call find_images_grid with
     *  rough = 1, which returns a list of all grid points associated
     *  with rays that hit inside the source.
     * First free the memory associated with the coordinate and grav.
     *  potential grids.
     */

    grid = del_posinfo(grid);
    poten = del_poten(poten);

    /*
     * Now set the grid parameters to a finer scale
     */

    modinfo->pixscale /= SCALEUP;
    size *= SCALEUP;

    /*
     * Set up the grids with the new parameters
     */

    printf("Setting up new coordinate grid\n");
    if(!(grid = make_grid(size,modinfo->pixscale))) {
      fprintf(stderr,"Exiting program\n");
      return 1;
    }
    printf("Setting up new potential grid\n");
    if(!(poten = new_poten(size,size))) {
      fprintf(stderr,"Exiting program\n");
      return 1;
    }
    printf("Calculating potential\n");
    for(i=0,parptr=modinfo->lp; i<modinfo->nmod; i++,parptr++)
      poten_table[parptr->pottyp].fn(poten,grid,parptr,size,size,1);

    /*
     * Get information about the (circular) source
     */

    printf("Enter source size in arcsec:  ");
    gets(line);
    sscanf(line,"%f",&sourcesize);

    /*
     * Find the images and plot them
     */

    printf("Finding images\n");
    if(find_images_grid(modinfo,nzfact,size,grid,poten,sourcesize,1,0)) {
      fprintf(stderr,"Exiting program\n");
      return 1;
    }
    plot_pts(modinfo->posimmod,size,modinfo->pixscale,modinfo->nimmod,6,17);
    break;
   
  default:

    /*
     * Error message for bad choice
     */

    fprintf(stderr,"Not a valid choice -- continuing with program\n");
  }
#else
  plot_obsim(modinfo->posimmod,size,modinfo);

#endif
  /*
   * Plot observed image positions
   */


  plot_obsim(modinfo->posobs,size,modinfo);

  return 0;
}
