/*
 * 18Feb96 CDF, First working version
 * v14Jul97 CDF, Major rewrite
 * v07Oct97 CDF, Fixed bug in modinfo_io in which the source position
 *                was not getting transferred to the observed posinfo
 *                structure.
 * v08Oct98 CDF, Link to dataio.h library.
 *               Improvements to posinfo_io
 */

#include <stdio.h>
#include <stdlib.h>
#include "structdef.h"
#include "dataio.h"
#include "file_io.h"

#define PI 3.141592653589793
#define D2R (PI/180.0)
#define NPAR 5

/*.......................................................................
 *
 * Function: lp_io
 *
 * Gets lens parameter information from an input file.
 *
 * Inputs:  char *inpname      name of input file
 *          char *name         name of the source - set by function
 *          int *nmod          number of model types - set by function
 *
 * Output:  Lensparams *lp     lens parameters - set by function
 *
 */

Lensparams *lp_io(char *inpname, char *name, int *nmod)
{
  char line[MAXC];      /* General string for reading input */
  char modstr[MAXC];    /* String passed to get_modpars */
  char *cptr;           /* Pointer used to navigate strings */
  Lensparams *lp=NULL;  /* Array of lens parameters filled by this function */
  Lensparams *parptr;   /* Pointer for navigating lp */
  Modparam *modp=NULL;  /* Temporary array with parameter values */
  FILE *inp;            /* Input file pointer */

  /*
   * Check the filename
   */

  strcpy(line,inpname);
  while(!(inp = fopen(line,"r"))) {
    printf("%s is not a valid file name.  Enter the file name again:  ",
	   line);
    gets(line);
  }

  /*
   * Read in the first line, which contains the source name
   */

  fgets(line,MAXC,inp);
  sscanf(line,"%s",name);

  /*
   * Count the number of model elements
   */

  while(fgets(line,MAXC,inp) != NULL)
    if (line[0] != '#') 
      (*nmod)++;
  printf("There are %d model types.\n",*nmod);
  rewind(inp);

  /*
   * Allocate memory for the the lp array
   */

  if(!(lp = new_lp(*nmod))) {
    fprintf(stderr,"Error:  lp_io\n");
    return NULL;
  }

  /*
   * Put the info into the array
   */

  parptr = lp;
  fgets(line,MAXC,inp);
  while(fgets(line,MAXC,inp) != NULL) {
    if (line[0] != '#') {

      /*
       * Read in potential type
       */

      if(sscanf(line,"%d",&parptr->pottyp) != 1) {
	  fprintf(stderr,"Error: lp_io.  Could not read input from %s",
		  inpname);
	  return NULL;
      }

      /*
       * Then set modstr to be the remainder of the line
       */

      cptr = line;
      while(*cptr == ' ')
	cptr++;
      while(*cptr != ' ' && *cptr != '\n' && *cptr != '\0')
	cptr++;
      strcpy(modstr,cptr);
      if(!(modp = get_modpars(modstr,NPAR))) {
	fprintf(stderr,"ERROR: lp_io\n");
	return NULL;
      }
      parptr->xl = *modp;
      parptr->yl = *(modp+1);
      parptr->b  = *(modp+2);
      parptr->gamma = *(modp+3);
      parptr->theta_g = *(modp+4);
      parptr->theta_g.par *= D2R;
      parptr++;
      modp = del_modparam(modp);
    }
  }

  if(inp)
    fclose(inp);

  return lp;
}

/*.......................................................................
 *
 * Function get_modpars
 *
 * Takes an input string and reads the npar parameter values from it.
 *  If any of the parameter values is followed immediately by a 'v' or
 *  'V', then that parameter will vary in the modelfitting and the
 *  associated varflag is set to 1.
 *
 * Inputs: char *modstr        string containing parameter values
 *         int npar            number of parameter values
 *
 * Output: Modparam *modpar    array of parameters and varflags
 *
 */

Modparam *get_modpars(char *modstr, int npar)
{
  int go=1;        /* Flag set to 0 when loop should stop */
  int count=0;     /* Number of parameters read in */
  char junk[MAXC]; /* String used for reading parameters */
  char *cptr;      /* Pointer used to step through the string */
  Modparam *modp;  /* Array filled by the function */
  Modparam *mptr;  /* Pointer to navigate modp */

  /*
   * Allocate memory for the Modparam array
   */

  if(!(modp = new_modparam(npar))) {
    fprintf(stderr,"ERROR: get_modpars\n");
    return NULL;
  }

  /*
   * Initialize
   */

  mptr = modp;
  cptr = modstr;
  while(go) {

    /*
     * Move cptr over blank spaces
     */

    while(*cptr == ' ')
      cptr++;

    /*
     * Check for end of line
     */

    if(*cptr == '\n' || *cptr == '\0') {
      go = 0;
      break;
    }

    /*
     * If OK, read in the parameter value
     */

    else {
      count++;
      if(count <= npar) {
	strcpy(junk,cptr);
	if(sscanf(junk,"%f",&mptr->par) != 1) {
	  fprintf(stderr,"ERROR: get_modpars. Problem reading input\n");
	  return (del_modparam(modp));
	}

	/*
	 * Move cptr to the end of the float
	 */

	while(*cptr != 'v' && *cptr != 'V' && *cptr != ' ' && *cptr != '\n'
	      && *cptr != '\0')
	  cptr++;

	/*
	 * Check to see if this is a free parameter.  If it is, read its
	 *  varflag (which follows the 'v').
	 * NB: If nothing follows the 'v', then the varflag is set to 1
	 */

	if(*cptr == 'v' || *cptr == 'V') {
	  cptr++;
	  if(*cptr == ' ' || *cptr == '\n' || *cptr == '\0')
	    mptr->varflag = 1;
	  else {
	    strcpy(junk,cptr);
	    if(sscanf(junk,"%d",&mptr->varflag) != 1) {
	      fprintf(stderr,"ERROR: get_modpars. Problem reading input\n");
	      return (del_modparam(modp));
	    }
	    cptr++;
	  }
	}
	else
	  mptr->varflag = 0;
	mptr++;
      }
      else
	go = 0;
    }
  }

  if(count < npar) {
    fprintf(stderr,"ERROR: get_modpars. Unexpected number of paramters\n");
    return del_modparam(modp);
  }
  else if(count > npar) {
    fprintf(stderr,"ERROR: get_modpars. \n");
    fprintf(stderr," More parameters in input file than expected.\n");
    fprintf(stderr," Only using first %d parameters on each model line\n",
	    npar);
  }

  return modp;
}

/*.......................................................................
 *
 * Function:  posinfo_io
 *
 * Gets image positions, fluxes and redshift factors from an input file.
 *
 * Inputs:  char *inpname,     the name of the input file
 *          int  *nobs,        the number of observed images - set
 *                               by this function
 *
 * Output:  Posinfo *pi,       the position information
 *
 * v08Oct98 CDF, Added confirmation of input values. 
 *               Better error checking and documentation.
 */

Posinfo *posinfo_io(char *inpname, int *nobs)
{
  int no_error=1;       /* Flag set to 0 on error */   
  char line[MAXC];      /* General string for reading input */
  Posinfo *pi=NULL;     /* Container for position info */
  Posinfo *piptr;       /* Pointer to navigate pi */
  FILE *ifp=NULL;       /* Input file pointer */

  if (!(ifp = fopen(inpname,"r"))) {
    printf("%s is not a valid file name.  Enter the file name again.\n",
	   inpname);
    while(!(ifp = fopen(gets(line),"r"))) 
      printf(" Still not correct.  Enter the file name again.\n");
  }

  /*
   * Get number of sources in input file
   */

  if((*nobs = n_lines(ifp,'#')) == 0) {
    fprintf(stderr,"ERROR: posinfo_io.  No valid lines in input file\n");
    return NULL;
  }
  else {
    printf("posinfo_io: %d lines in input file.\n",*nobs);
    rewind(ifp);
  }

  /*
   * Allocate memory for image information
   */

  if(!(pi = new_posinfo(*nobs,1))) {
    fprintf(stderr,"ERROR:  posinfo_io\n");
    return NULL;
  }

  /*
   * Loop through input file, reading image inforation
   */

  piptr = pi;
  printf("\n--------------------------------------------------\n");
  printf(" Image data\n\n");
  printf("      x       y      flux    zfactor\n");
  printf("  -------- -------- -------- -------\n");
  while(fgets(line,MAXC,ifp) != NULL && no_error) {
    if (line[0] != '#') {
      if(sscanf(line,"%f %f %f %f",
		&piptr->x,&piptr->y,&piptr->flux,&piptr->zf) != 4) {
	fprintf(stderr,"ERROR: posinfo_io.  Bad input format.\n");
	no_error = 0;
      }
      else {
	printf("  %8.4f %8.4f %7f  %5.2f\n",piptr->x,piptr->y,piptr->flux,
	       piptr->zf);
	piptr++;
      }
    }
  }
  if(no_error)
    printf("\n--------------------------------------------------\n\n");
  
  /*
   * Clean up and return
   */

  if(ifp)
    fclose(ifp);

  if(no_error)
    return pi;
  else {
    fprintf(stderr,"ERROR: posinfo_io.\n");
    return del_posinfo(pi);
  }
}


/*.......................................................................
 *
 * Function:  fsource_io
 *
 * Gets source fluxes from an input file.
 *
 * Inputs:  char *inpname      name of the input file
 *          int  *nfsrc        number of observed images - set
 *                               by this function
 *          int *nofile        if *inpname == NULL set to 1, otherwise 0
 *
 * Output:  Srcinfo *fsource,  the source fluxes and redshift factors
 *
 */

Srcinfo *fsource_io(char *inpname, int *nfsrc, int *nofile)
{
  char line[MAXC];        /* General string variable for reading input */
  Srcinfo *fsource=NULL;  /* Source information array -- filled by function */
  Srcinfo *srcptr;        /* Pointer to navigate fsource */
  Modparam *modp=NULL;    /* Source parameters returned by get_modpars */
  FILE *inp;              /* Input file pointer */

  *nofile = 0;

  /*
   * If inpname == NULL, then there is no source input file.  In this case,
   *  give the source the default parameters
   */

  if(!inpname) {
    *nfsrc = 1;
    *nofile = 1;
    if(!(fsource = new_srcinfo(1,1))) {
      fprintf(stderr,"Error: fsource_io\n");
      return NULL;
    }
    fsource->x.par = 0.0;
    fsource->y.par = 0.0;
    fsource->flux.par = 1.0;
    fsource->zf.par = 1.0;
    fsource->x.varflag = fsource->y.varflag = 1;
    fsource->flux.varflag = fsource->zf.varflag = 0;
    printf("\nNo source flux file.  Setting source flux = 1.0\n\n");
    return fsource;
  }

  /*
   * Otherwise read source parameters from input file
   *
   * First check existence of input file
   */

  else {
    if (!(inp = fopen(inpname,"r"))) {
      printf("%s is not a valid file name.  Enter the file name again.\n",
	     inpname);
      while(!(inp = fopen(gets(line),"r"))) 
	printf(" Still not correct.  Enter the file name again.\n");
    }
  }

  /*
   * Count the number of sources in the input file
   */

  while(fgets(line,MAXC,inp) != NULL)
    if (line[0] != '#')
      (*nfsrc)++;
  if(!(fsource = new_srcinfo(*nfsrc,1))) {
    fprintf(stderr,"Error:  fsource_io\n");
    return NULL;
  }
  rewind(inp);

  /*
   * Read the input information
   */

  srcptr = fsource;
  while(fgets(line,MAXC,inp) != NULL) {
    if (line[0] != '#') {
      if(!(modp = get_modpars(line,4))) {
	fprintf(stderr,"ERROR: fsource_io.\n");
	return del_srcinfo(fsource);
      }
      srcptr->x = *modp;
      srcptr->y = *(modp+1);
      srcptr->flux = *(modp+2);
      srcptr->zf = *(modp+3);
      modp = del_modparam(modp);
      srcptr++;
    }
  }

  /*
   * Clean up
   */

  if(inp)
    fclose(inp);

  return fsource;
}


/*.......................................................................
 *
 * Function: modinfo_io
 *
 * Reads information into a modinfo structure from a pair of input files
 *
 * Inputs:  char *imposname    file containing observed image info
 *          char *datname      file containing the model information
 *          char *fsourcename  file containing the source info
 *          char *name         source name - set by this function
 *
 * Ouput:   Modinfo *modinfo   filled structure
 *
 */

Modinfo *modinfo_io(char *imposname, char *datname, char *fsourcename, 
		    char *name)
{
  int i,j;                 /* Looping variables */
  int nobs=0;              /* Number of observed images */
  int nmod=0;              /* Number of potentials */
  int nfsrc=0;             /* Number of sources */
  int nofile=0;            /* Flag set to 1 if there is no source file */
  int *var=NULL;           /* Array of flags for model fitting */
  int *intptr;             /* Pointer to navigate var */
  Posinfo *posinfo=NULL;   /* Image positions and fluxes */
  Posinfo *piptr;          /* Pointer to navigate posinfo */
  Srcinfo *fsource=NULL;   /* Source position and fluxes */
  Srcinfo *srcptr;         /* Pointer to navigate fsource */
  Lensparams *lpinit=NULL; /* Initial guess for lens parameters */
  Lensparams *parptr;      /* Pointer used to navigate lp */
  Modinfo *modinfo;

  /*
   * Fill position information from imposname
   */

  if(!(posinfo = posinfo_io(imposname,&nobs))) {
    fprintf(stderr,"Error:  modinfo_io\n");
    return NULL;
  }

  /*
   * Get number of source fluxes from fsourcename (if it exists)
   */

  if(!(fsource = fsource_io(fsourcename,&nfsrc,&nofile))) {
    fprintf(stderr,"Error: modinfo_io/n");
    return NULL;
  }

  /*
   * Read source name and initial guess of lens parameters from 
   *   data input file.
   */

  if(!(lpinit = lp_io(datname,name,&nmod))) {
    fprintf(stderr,"ERROR: modinfo_io\n");
    return NULL;
  }

  /*
   * Put appropriate source flux info into posinfo
   */

  piptr = posinfo;
  for(i=0; i<nobs; i++,piptr++) {
    srcptr = fsource;
    for(j=0; j<nfsrc; j++,srcptr++)
      if(piptr->zf == srcptr->zf.par) {
	piptr->sflux = srcptr->flux.par;
	piptr->spos.x = srcptr->x.par;
	piptr->spos.y = srcptr->y.par;
	break;
      }
  }

  /*
   * Create a new Modinfo structure and fill it with the information read
   *  in so far.
   */

  if(!(modinfo = new_modinfo())) {
    fprintf(stderr,"Error: modinfo_io\n");
    return NULL;
  }

  modinfo->nmod = nmod;
  modinfo->nobs = nobs;
  modinfo->nfsrc = nfsrc;
  modinfo->lp = lpinit;
  modinfo->posobs = posinfo;
  modinfo->fsource = fsource;
  modinfo->posimmod = NULL;

  /*
   * Print out the model parameters and, at the same time, put the
   *  varflags from the parameters into the var array in modinfo.
   */

  if(!(var = new_intarray(5*nmod+4*nfsrc,1))) {
    fprintf(stderr,"ERROR: modinfo_io\n");
    return NULL;
  }
  intptr = var;
  parptr = modinfo->lp;
  printf("\nModel parameters\n");

  /*
   * First the lens parameters 
   */

  for(i=0; i<nmod; i++,parptr++) {
    printf("%8.4f %d %8.4f %d %8.4f %d %8.4f %d %8.4f %d\n",
	   parptr->xl.par,parptr->xl.varflag,parptr->yl.par,
	   parptr->yl.varflag,parptr->b.par,parptr->b.varflag,
	   parptr->gamma.par,parptr->gamma.varflag,parptr->theta_g.par,
	   parptr->theta_g.varflag);
    *intptr++ = parptr->xl.varflag;
    *intptr++ = parptr->yl.varflag;
    *intptr++ = parptr->b.varflag;
    *intptr++ = parptr->gamma.varflag;
    *intptr++ = parptr->theta_g.varflag;
  }

  /*
   * Then the source paramters
   */

  srcptr = modinfo->fsource;
  for(i=0; i<nfsrc; i++,srcptr++) {
    printf("%8.4f %d %8.4f %d %8.4f %d %8.4f %d\n",
	   srcptr->x.par,srcptr->x.varflag,srcptr->y.par,srcptr->y.varflag,
	   srcptr->flux.par,srcptr->flux.varflag,srcptr->zf.par,
	   srcptr->zf.varflag);
    *intptr++ = srcptr->x.varflag;
    *intptr++ = srcptr->y.varflag;
    *intptr++ = srcptr->flux.varflag;
    *intptr++ = srcptr->zf.varflag;
  }

  modinfo->var = var;

  return modinfo;
}

/*.......................................................................
 *
 * Function check_nsrc
 *
 * Reads in the number of redshift factors in the observed image input
 *  file and calculates the number of sources implied (by calling the
 *  function n_sources).  This value is compared to the number of sources
 *  read in from the source input file and stored in modinfo->fsource.
 *  If the two don't match, an error is returned, otherwise the float
 *  array of redshift factors is returned.
 *
 * Input:  Modinfo *modinfo    model information
 *
 * Output: float *nzfact       array of redshift factors
 *
 */

float *check_nsrc(Modinfo *modinfo)
{
  int nz=1;               /* Number of redshift factors */
  float *nzfact=NULL;     /* Array of redshift factors */

  /*
   * Get the number of redshift factors from the image info
   */

  if(!(nzfact = n_sources(modinfo->posobs,modinfo->nobs,&nz))) {
    fprintf(stderr,"ERROR: check_nsrc\n");
    return NULL;
  }

  /*
   * Check to make sure number of source fluxes matches number of
   *  redshift factors
   */

  if(modinfo->nfsrc != nz) {
    fprintf(stderr,"ERROR: check_nsrc.  Number of source fluxes ");
    fprintf(stderr,"doesn't match number of redshift factors\n");
    return del_array(nzfact);
  }

  /*
   * Sort order of source fluxes in order to match order of redshift factor
   */

  if(sort_fsource(modinfo->fsource,nzfact,nz)) {
    fprintf(stderr,"ERROR: check_nsrc\n");
    return del_array(nzfact);
  }

  return nzfact;
}

/*.......................................................................
 *
 * Function n_sources
 *
 * Calculates the number of sources based on the number of redshift
 *   factors.
 *
 * Inputs:  float *zfact       list of redshift factors
 *          int nobs           number of observed image positions
 *          int *nz            number of sources
 *
 * Output:  float *nzfact      list of source redshift factors
 *
 */

float *n_sources(Posinfo *piobs, int nobs, int *nz)
{
  int i,j,k;           /* Looping variables */
  int flag;            /* Flag set to 1 if a new redshift factor is found */
  float *nzfact=NULL;  /* Array of redshift factors */
  float *tmpzf=NULL;   /* Temporary holder for redshift factors */

  /*
   * Allocate memory for temporary storage of redshift factors, and
   *  fill temporary array.
   */

  if(!(tmpzf = new_array(nobs,1))) {
    fprintf(stderr,"Error:  n_sources");
    return NULL;
  }

  for(i=0;i<nobs;i++)
    *(tmpzf+i) = piobs->zf;

  /*
   * Go through the image info and find the unique redshift factors.  
   *  Put these factors into the first nz positions in the temporay
   *  array.
   */

  for(i=0; i<nobs-1; i++)
    for(j=i+1; j<nobs; j++)
      if(piobs[i].zf != piobs[j].zf) {
	flag = 1;
	for(k=0; k<nobs; k++)
	  if(piobs[j].zf == *(tmpzf+k))
	    flag = 0;
	if(flag) {
	  *(tmpzf+*nz) = piobs[j].zf;
	  (*nz)++;
	}
      }

  /*
   * Allocate memory for the redshift factor array and fill it with the
   *  first nz members of the temporary array.
   */

  if(!(nzfact = new_array(*nz,1))) {
    fprintf(stderr,"Error:  n_sources");
    return NULL;
  }
  printf("\nThere are %d redshift factors.\n",*nz);
  for(i=0; i<*nz; i++) {
    *(nzfact+i) = *(tmpzf+i);
    printf("Redshift factor %d is %7.4f\n",i+1,*(nzfact+i));
  }
  printf("\n");

  /*
   * Clean up
   */

  tmpzf = del_array(tmpzf);

  return nzfact;
}

/*.......................................................................
 *
 * Function sort_fsource
 *
 * Sorts the order of the source fluxes to match that in the redshift
 *  factor list, based on the redshift factors associated with the
 *  source fluxes
 *
 * Inputs: Srcinfo *fsource    source fluxes (modified by this function)
 *         float *nzfact       redshift factors
 *         int nz              number of redshift factors (and source fluxes)
 *
 * Output: int (0 or 1)        0 ==> successful, 1 ==> error
 *
 */

int sort_fsource(Srcinfo *fsource, float *nzfact, int nz)
{
  int i,j;                /* Looping variables */
  float *fptr;            /* Pointer for navigating nzfact */
  Srcinfo *srcptr;        /* Pointer for navigating fsource */
  Srcinfo *tmpsrc=NULL;   /* Temporary container for sorted source info */
  Srcinfo *tmpptr;        /* Pointer for navigating tmpsrc */

  /*
   * Allocate memory for temporary storage
   */

  if(!(tmpsrc = new_srcinfo(nz,1))) {
    fprintf(stderr,"Error: sort_fsource\n");
    return 1;
  }

  /*
   * Put source info into temporary array in order that the sources are
   *  found in the nzfact array
   */

  tmpptr = tmpsrc;
  for(i=0,fptr=nzfact; i<nz; i++,fptr++) {
    for(j=0,srcptr=fsource; j<nz; j++,srcptr++)
      if(srcptr->zf.par == *fptr) {
	*tmpptr++ = *srcptr;
	break;
      }
  }

  /*
   * Transfer sorted source info into fsource
   */

  for(i=0,tmpptr=tmpsrc,srcptr=fsource; i<nz; i++,tmpptr++,srcptr++)
    *srcptr = *tmpptr;

  /*
   * Clean up
   */

  tmpsrc = del_srcinfo(tmpsrc);

  return 0;
}
