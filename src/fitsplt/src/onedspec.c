/*
 * onedspec.c
 *
 * A library of functions for handling and plotting one-dimensional
 *  spectra.
 *
 * 31Aug01 CDF   Split off from specplot.c
 * v13Sep01 CDF, Moved get_label_style, get_title, and get_clabel functions
 *                into setup_spec.c.
 * v17Feb03 CDF, Combined Pos spectrum array and Splotinfo structure into
 *                a new structure called Spectrum.
 * v18Feb03 CDF, Propogated the new Spectrum structure through into most
 *                of the functions.
 *               Added new read_onespec and read_multispec functions
 *                to take functionality previously in specplot.c.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cpgplot.h"
#include "structdef.h"
#include "plotfuncs.h"
#include "dataio.h"
#include "onedspec.h"

/*......................................................................
 *
 * Functions to handle memory allocation for arrays of Spectrum structures:
 *    new_spectrum  -   allocates the array
 *    del_spectrum  -   frees the array and the individual members
 *
 */

Spectrum *new_spectrum(int nspec) {
  Spectrum *newspec=NULL;  /* Array of spectra to be allocated */

  newspec = (Spectrum *) malloc(sizeof(Spectrum) * nspec);
  if(!newspec) {
    fprintf(stderr,"ERROR:  Insufficient memory for array of spectra.\n");
    return NULL;
  }

  /*
   * Initialize the spectrum parameters.
   */

  newspec->dolabel = 0;
  newspec->labpos = 'r';
  newspec->spec = NULL;

  return newspec;
}

Spectrum *del_spectrum(Spectrum *spec) {
  if(spec) 
    free(spec);

  return NULL;
}

/*......................................................................
 *
 * Functions to handle memory allocation for arrays of pointers for spectra:
 *    new_spec  -   allocates the array
 *    del_spec  -   frees the array and the individual members
 *
 */

Pos **new_spec(int nspec) {
  Pos **newspec={NULL};  /* Array of spectra to be allocated */

  newspec = (Pos **) malloc(sizeof(Pos *) * nspec);
  if(!newspec) {
    fprintf(stderr,"ERROR:  Insufficient memory for array of pointers\n");
    fprintf(stderr,"        to the spectra.");
    return NULL;
  }
  else
    return newspec;
}

Pos **del_spec(Pos **spec, int nspec) {
  int i;

  if(spec) {
    for(i=0; i<nspec; i++)
      del_pos(spec[i]);
    free(spec);
  }
  return NULL;
}

/*.......................................................................
 *
 * Function read_onespec
 *
 * Reads in one spectrum by calling read_spectrum..
 *
 * Inputs: char *filename      name of input file containing spectrum
 *
 * Output: Spectrum *spectrum  structure containing spectrum and associated
 *                              info
 *
 */

Spectrum *read_onespec(char *filename)
{
  int no_error=1;            /* Flag set to 0 on error */
  char line[MAXC];           /* General string variable */
  Spectrum *spectrum=NULL;   /* Array for spectrum */
  FILE *ifp=NULL;            /* Pointer to file containing spectrum */

  /*
   * Allocate memory for Spectrum structure
   */

  if(!(spectrum = new_spectrum(1))) {
    fprintf(stderr,"ERROR: read_onespec.\n");
    return NULL;
  }

  /*
   * Read in spectrum
   */

  if(read_spectrum(filename,spectrum))
    no_error = 0;

  /*
   * Return spectrum if no problems
   */

  if(no_error)
    return spectrum;
  else {
    fprintf(stderr,"ERROR: read_onespec.\n");
    return NULL;
  }
}

/*.......................................................................
 *
 * Function read_multispec
 *
 * Reads in multiple spectra by calling read_spectrum for each input
 *  spectrum.  In this case, the input filename passed to the function
 *  contains a list of files, each of which contains a spectrum.  This
 *  is in contrast to the read_onespec function, in which the passed
 *  filename is the name of the file containing the input spectrum.
 *
 * Inputs: char *filename      name of input file
 *         int *nspec          number of input spectra (set by this function)
 *
 * Output: Spectrum *spectrum  structure containing spectrum and associated
 *                              info
 *
 */

Spectrum *read_multispec(char *filename, int *nspec)
{
  int no_error=1;            /* Flag set to 0 on error */
  char line[MAXC];           /* General string variable */
  char specfile[MAXC];     /* Name of file containing an ASCII spectrum */
  Spectrum *spectrum=NULL;   /* Array for spectrum */
  Spectrum *sp;              /* Pointer to navigate spectrum */
  FILE *ifp=NULL;            /* Pointer to file containing list of spectra  */

  /*
   * Open file containing list of spectra and count number of spectra
   *  in the list.
   */

  if(!(ifp = open_readfile(filename))) {
    fprintf(stderr,"ERROR: read_multispec.\n");
    return NULL;
  }
  else {
    *nspec = n_lines(ifp,'#');
    rewind(ifp);
  }

  /*
   * Allocate memory for array of Spectrum structures
   */

  printf("\nread_multispec: Allocating memory for %d spectra.\n",*nspec);
  if(!(spectrum = new_spectrum(*nspec))) 
    no_error = 0;

  /*
   * Read in spectra
   */

  if(no_error) {
    sp = spectrum;
    while(fgets(line,MAXC,ifp) != NULL && no_error) {
      if(line[0] != '#') {
	if(sscanf(line,"%s",specfile) != 1)
	  no_error = 0;
	else {
	  if(read_spectrum(specfile,sp))
	    no_error = 0;
	  else
	    sp++;
	}
      }
    }
  }


  /*
   * Return spectrum if no problems
   */

  if(ifp)
    fclose(ifp);
  if(no_error)
    return spectrum;
  else {
    fprintf(stderr,"ERROR: read_multispec.\n");
    return NULL;
  }
}

/*.......................................................................
 *
 * Function read_spectrum
 *
 * Given an input filename, this function opens the file, counts the number
 *  of lines in the file, reads in the wavelength and flux data from the
 *  file, and stores those data in the appropriate array.
 *
 * **NB: The spectrum structure MUST be allocated before calling this
 *        function.
 *
 * Inputs: char *filename      name of input file
 *         Spectrum *spectrum  spectrum and associated info (read in by
 *                              this function).
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int read_spectrum(char *filename, Spectrum *spectrum)
{
  int no_error=1;            /* Flag set to 0 on error */
  char line[MAXC];           /* General string variable */
  FILE *ifp=NULL;            /* Pointer to file containing spectrum */

  /*
   * Open input file
   */

  if(!(ifp = open_readfile(filename))) {
    fprintf(stderr,"ERROR:  read_spectrum\n");
    return 1;
  }

  /*
   * Count number of lines in input file
   */

  if((spectrum->nlines = n_lines(ifp,'#')) < 0) {
    fprintf(stderr,"ERROR:  read_spectrum\n");
    return 1;
  }

  /*
   * Read in spectrum
   */

  if(no_error) {
    rewind(ifp);
    if(read_spec_list(spectrum,ifp))
      no_error = 0;
  }

  if(ifp)
    fclose(ifp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: read_spectrum.\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function read_spec_list
 *
 * Reads data from input file which contains two columns, wavelength and
 *  flux.  Puts data. as well asmin and max flux values into the Spectrum
 *  structure.
 *
 * Inputs: Spectrum *spectrum  wavelengths, fluxes, and associated info
 *         FILE *ifp           input file
 *
 * Output: int (0 or 1)        0 ==> success, 1==> error
 *
 */

int read_spec_list(Spectrum *spectrum, FILE *ifp)
{
 int i;            /* Counting variable */
 float min;        /* Minimum flux */
 float max;        /* Maximum flux */
 float fluxdiff;   /* Difference between min and max */
 char line[MAXC];  /* General string for reading input */
 Pos *sp;          /* Pointer to navigate spectrum->spec */

 /*
  * Allocate memory and initialize
  */

 if(!(spectrum->spec = new_pos(spectrum->nlines,1))) {
   fprintf(stderr,"ERROR: read_spec_list\n");
   return 1;
 }

 sp = spectrum->spec;
 i = 0;

 /*
  * Read in data
  */

 while(fgets(line,MAXC,ifp) != NULL) {
   if(line[0] != '#') {
     if(sscanf(line,"%lf %lf",&sp->x,&sp->y) != 2) {
       fprintf(stderr,"ERROR: read_spec_list\n");
       return 1;
     }

     /*
      * Initialize min and max
      */

     if(i == 0) {
       max = sp->y;
       min = sp->y;
       spectrum->minlambda = sp->x;
     }

     /*
      * Check if new min or max need to be set
      */

     if(sp->y > max)
       max = sp->y;
     if(sp->y < min)
       min = sp->y;
     sp++;
     i++;
   }
 }

 if(i != spectrum->nlines) {
   fprintf(stderr,"ERROR: read_spec_list.  Unexpected number of lines\n");
   return 1;
 }
 else {
   spectrum->maxlambda = (sp-1)->x;
   printf("read_spec_list: Read in %d lines.\n",spectrum->nlines);
   printf("read_spec_list: Min lambda = %7.1f, max lambda = %7.1f\n",
	  spectrum->minlambda,spectrum->maxlambda);
   printf("read_spec_list: Min flux = %e, max flux = %e\n",min,max);
 }

 fluxdiff = max - min;
 spectrum->maxflux = max + 0.35*fluxdiff;
 spectrum->minflux = min - 0.35*fluxdiff;

 return 0;
}

/*......................................................................
 *
 * Function read_line_list
 *
 * Reads list of spectral line names and rest wavelengths into a Specinfo
 *   structure.
 *
 * Inputs: int *goodlines      number of lines read into structure (set
 *                              by this function)    
 *
 * Output: Specinfo *specinfo  structure containing spectral info
 *
 * v25Apr99 CDF, Added second check for rest-wavelength file.  If it is
 *                not in the expected location, then check in the current
 *                working directory before querying the user for a new
 *                filename.
 */

Specinfo *read_line_list(int *goodlines)
{
  int ngood=0;
  int nlines=0;
  char specfile[MAXC];        /* Filename for rest-wavelength file */
  char line[MAXC];            /* General string variable */
  Specinfo *specinfo;         /* Structure containing spectral info */
  Specinfo *sptr;             /* General specinfo pointer */
  FILE *linfp;                /* File containing line list */

  /*
   * Open rest-wavelength file.  First check for home-workstation location
   */

  sprintf(specfile,"/Users/cdf/Spectroscopy/speclines.list");

  if(!(linfp = fopen(specfile,"r"))) {

    fprintf(stderr,"\n*** Warning: read_line_list. ***\n");
    fprintf(stderr,"*** %s not found. ***\n",specfile);
    fprintf(stderr,"  Checking current directory for speclines.list....\n");

    /*
     * Now check current directory for speclines.list file
     */

    sprintf(specfile,"speclines.list");
    if(!(linfp = open_readfile(specfile))) {
      fprintf(stderr,"\nERROR: read_line_list.\n");
      return NULL;
    }
  }
  else {
    printf("\nread_line_list: %s now open.\n",specfile);
  }

  /*
   * Count number of spectral lines in file
   */

  if((nlines = n_lines(linfp,'#')) < 0) {
    fprintf(stderr,"ERROR:  read_line_list\n");
    return NULL;
  }

  /*
   * Allocate specinfo container
   */

  if(!(specinfo = new_specinfo(nlines))) {
    fprintf(stderr,"\nError: read_line_list\n");
    return NULL;
  }

  /*
   * Read lines into structure and check line format
   */

  rewind(linfp);
  sptr = specinfo;
  while(fgets(line,MAXC,linfp) != NULL) {
    if(line[0] != '#') {
      switch(read_line(line,sptr)) {
      case 2: case 3:
	ngood++;
	sptr++;
	break;
      case -1:
	fprintf(stderr,"ERROR: read_line_list\n");
	return NULL;
	break;
      default:
	break;
      }
    }
  }

  *goodlines = ngood;
  if(*goodlines == 0) {
    fprintf(stderr,"\nERROR: No lines read from spectral line file\n");
    return NULL;
  }

  if(linfp)
    fclose(linfp);

  return specinfo;
}

/*......................................................................
 *
 * Function read_line
 *
 * Reads one line from the spectral line list and checks to see if it
 *   is in the right format.
 *
 * Inputs: char *line,           Input line from the file
 *         Specinfo *specinfo,   Pointer to spectral line structure
 *
 * Output: int nhash             Number of hash marks found
 *
 */

int read_line(char *line, Specinfo *specinfo)
{
  int nhash=0;                /* Number of hash marks on the input line */
  char *cptr1,*cptr2;         /* General char pointers */
  char *str1,*str2;
  char *str3,*str4;

  cptr1 = line;
  if((cptr2 = strchr(line,'#'))) {
    str2 = cptr2;
    nhash++;
    str1 = cptr1;
    while(cptr1 < cptr2)
      cptr1++;
    *cptr1='\0';
  }
  else
    return nhash;

  cptr1 = str2+1;
  if((cptr2 = strchr(cptr1,'#'))) {
    str3 = cptr2+1;
    nhash++;
    str2 = cptr1;
    while(cptr1 < cptr2)
      cptr1++;
    *cptr1='\0';
  }
  else
    return nhash;

  cptr1 = str3+1;
  if((cptr2 = strchr(cptr1,'#'))) {
    str4 = cptr2+1;
    nhash++;
    str3 = cptr1;
    while(cptr1 < cptr2)
      cptr1++;
    *cptr1='\0';
  }

  switch(nhash) {
  case 2:
    strcpy(specinfo->ionname,str1);
    sscanf(str2,"%f",&specinfo->lambda);
    strcpy(specinfo->label,str3);
    break;
  case 3:
    strcpy(specinfo->ionname,str1);
    sscanf(str2,"%f",&specinfo->lambda);
    strcpy(specinfo->label,str3);
    strcpy(specinfo->comment,str4);
    break;
  default:
    break;
  }

  return nhash;
}

/*.......................................................................
 *
 * Function display_spectrum
 *
 * Displays the spectrum on the desired plotting device, which is
 *  passed as an argument to the function.
 *
 * Inputs: Spectrum *spec      spectrum/spectra to be plotted
 *         int nspec           number of spectra contained in spectrum
 *                              structure
 *         int specnum         id number of spectrum to be plotted (ignored
 *                              for OVERPLOT option)
 *         Spsetup *setup      container for plotting parameters
 *         int option          method of plotting spectra
 *         int dobox           flag defining axis labeling
 *                              0 ==> no labels
 *                              1 ==> fully-labeled
 *                              2 ==> y label only
 *                              3 ==> x labels only
 *                              4 ==> no bottom axis, y labels only
 *                              5 ==> no top axis, y labels only
 *                              6 ==> no top axis, x labels only
 *
 * Output: int (0 or 1)         0 ==> success, 1 ==> error
 *
 */

int display_spectrum(Spectrum *spec, int nspec, int specnum, 
		     Spsetup *setup, int option, int dobox)
{
  int i;                  /* Looping variable */
  int no_error = 1;       /* Flag set to 0 on error */
  int basespec;           /* Base spectrum number */
  int ltype=1;            /* Line type for plot */
  int pdobox=1;           /* Local version of dobox */
  float pltxmin,pltxmax;  /* Min and max wavelengths to plot on x axis */
  float pltymin,pltymax;  /* Min and max fluxes to plot on y axis */
  Spectrum *sp;           /* Pointer to navigate spec array */

  /*
   * Check passed value of specnum
   */

  if(specnum < 0 || specnum >= nspec) {
    fprintf(stderr,"ERROR: display_spectrum.\n");
    fprintf(stderr," Bad choice of spectrum to plot.\n");
    fprintf(stderr," Spectrum chosen was %d.  Valid spectra are 0 -> %d\n\n",
	    specnum,nspec);
    return 1;
  }

  /*
   * Just in case the user has passed a bad specnum for the ONESPEC
   *  or OVERPLOT options, correct the error.
   */

  if(option == ONESPEC || option == OVERPLOT)
    basespec = 0;
  else
    basespec = specnum;

  /*
   * Define flux range, based on whether fluxes are in magnitudes or not.
   */

  sp = spec + basespec;
  set_pltlims_y(sp->minflux,sp->maxflux,setup->flchoice,&pltymin,
		&pltymax);

  /*
   * Define the wavelength range to be that occupied by first spectrum
   */

  pltxmin = setup->pltxmin;
  pltxmax = setup->pltxmax;

  /*
   * Set the character font.
   */

  cpgscf(2);

  /*
   * Overplot all the spectra if the OVERPLOT option has been chosen.
   */

  if(option == OVERPLOT) {
    for(i=0,sp=spec; i<nspec; i++,sp++) {
      if(i==0) {
	ltype = 1;
	pdobox = 1;
	cpgsci(1);
      }
      else {
	ltype = 2;
	pdobox = 0;
	cpgsci(2);
      }
      if(plot_spec(sp->spec,sp->nlines,pltymin,pltymax,pltxmin,pltxmax,
		   ltype,pdobox,setup->cheight))
	no_error = 0;
    }
  }

  /*
   * Otherwise plot only the chosen spectrum.
   */

  else {
    cpgsci(1);

    pdobox = dobox;
    sp = spec + basespec;
    if(plot_spec(sp->spec,sp->nlines,pltymin,pltymax,pltxmin,pltxmax,
		 1,pdobox,setup->cheight))
	no_error = 0;
  }

  /*
   * Add a label in lower right or left corner of plot, if desired
   */

  if(no_error) {
    sp = spec + basespec; 
    if(sp->dolabel == 0)
      sp->dolabel = get_clabel(sp->labtext,&sp->labpos);
    if(sp->dolabel == 1)
      labcorner(sp->labtext,sp->labpos,option);
  }

  /*
   * Return
   */

  cpgsci(1);
  if(no_error)
    return 0;
  else {
    fprintf(stderr,"Error: display_spectrum.\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function set_pltlims_y
 *
 * Sets y limits for plotting.
 *
 * Inputs: float minflux       passed minimum flux
 *         float maxflux       passed maximum flux
 *         int flchoice        flag set to 1 if fluxes in magnitudes
 *         float *pltymin      proper minimum y value (set by this function)
 *         float *pltymax      proper maximum y value (set by this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

void set_pltlims_y(float minflux, float maxflux, int flchoice, 
		   float *pltymin, float *pltymax)
{
  /*
   * Define reverse scale if user has input magnitudes, otherwise use
   *  a standard scale for the y axis.
   */

  if(flchoice == 1 || flchoice == 2) {
    *pltymin = maxflux;
    *pltymax = minflux;
  }
  else {
    *pltymin = minflux;
    *pltymax = maxflux;
  }
}

/*.......................................................................
 *
 * Function adjust_spec
 *
 * Adjusts the display limits through calls to get_vscale.  This allows
 *  the user to fiddle with the display limits and see the effects on
 *  the screen before plotting the image to a file.
 *
 * Inputs: Spectrum *spectrum  spectrum/spectra to be plotted
 *         int nspec           number of spectra contained in spectrum
 *                              structure
 *         int specnum         id number of spectrum to be plotted (ignored
 *                              for OVERPLOT option)
 *         Spsetup *setup      container for plotting parameters
 *         int option          method of plotting spectra
 *         int dobox           flag defining axis labeling
 *                              0 ==> no labels
 *                              1 ==> fully-labeled
 *                              2 ==> y label only
 *
 * Output: 0 or 1              0 ==> success, 1 ==> error
 *
 */

int adjust_spec(Spectrum *spectrum, int nspec, int specnum, 
		Spsetup *setup, int option, int dobox)
{
  int no_error=1;      /* Flag set to 0 on error */
  int do_adjust=1;     /* Flag set to 0 when no more adjustment is needed */
  char line[MAXC];     /* String variable for reading input */
  Spectrum *sp;        /* Pointer for navigating spectrum */

  /*
   * First see if the user wants to change the current values.
   */

  while(do_adjust && no_error) {
    printf("\n--------------------------------------------------\n\n");
    printf("\nDo you want to change the vertical scale");
    printf(" of the spectrum number %d? [y] ",specnum);
    fgets(line,MAXC,stdin);
    if(line[0] == 'n' || line[0] == 'N')
      do_adjust = 0;

    /*
     * If so, get the new values.
     */

    else {
      sp = spectrum + specnum;
      switch(get_vscale(&sp->minflux,&sp->maxflux)) {
      case -1:
	no_error = 0;
	break;
      case 1:
	cpgeras();
	if(display_spectrum(spectrum,nspec,specnum,setup,option,dobox))
	  no_error = 0;
	break;
      default:
	do_adjust = 0;
      }
    }
  }

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: adjust_spec\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function get_vscale
 *
 * Gets min and max values for vertical axis of plot.
 *
 * Inputs: float *minflux      current value of min flux (may be modified 
 *                              by this function)
 *         float *maxflux      current value of max flux (may be modified 
 *                              by this function)
 *
 * Output: int  (-1 to 1)      0 ==> no change, 1 ==> minflux and maxflux
 *                              changed, -1 ==> error
 *
 */

int get_vscale(float *minflux, float *maxflux)
{
  float min,max;    /* New values for min and max flux */
  char line[MAXC];  /* General string variable */

  printf(" Enter desired min and max flux or magnitude values or just hit\n");
  printf(" return to keep current scale (%g %g):  ",*minflux,*maxflux);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%f %f",&min,&max) != 2) {
      fprintf(stderr,"ERROR: Bad input format.  Enter values again: ");
      fgets(line,MAXC,stdin);
    }
    *minflux = min;
    *maxflux = max;
    return 1;
  }
  else
    return 0;
}

/*.......................................................................
 *
 * Function label_lines
 *
 * Chooses spectral lines to label
 *
 * Inputs: Specinfo *specinfo  spectral line information
 *         int goodlines       number of spectral lines
 *         float z             redshift of lines
 *         Spectrum *spectrum  spectrum and associated info
 *         Spsetup *setup      container for plotting parameters
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 * v13Sep01 CDF, Moved the passed labstyle and flchoice parameters into
 *                the Spsetup structure.
 * v18Sep01 CDF, Fixed a bug in the plotting of a source line that is not
 *                on the line list.
 */

int label_lines(Specinfo *specinfo, int goodlines, float z, Spectrum *spectrum,
		Spsetup *setup)
{
  int i=1;
  int no_error=1;        /* Flag set to 0 on error. */
  int flag;              /* If error on input */
  int choice;            /* ID number of chosen line */
  float marklamb;        /* x position of line marker */
  float dzl;             /* Optional label offset in Angstroms */
  char line[MAXC];       /* General string for reading input */
  char *cptr;            /* Character pointer for parsing input */
  Labinfo label;         /* Contains x position and text for label only */
  Specinfo *sptr;
  Specinfo tmpspec;
  FILE *ofp=NULL;        /* Output file pointer */

  /*
   * Open output file.
   */

  printf("\n*** Opening file to hold marker info. ***\n");
  if(!(ofp=open_writefile("foo.LIST"))) {
    fprintf(stderr,"ERROR: label_lines.\n");
    return 1;
  }
  printf("\n");

  /*
   * Query user for line IDs and plot lines.
   */

  printf("Choose line #%d. (Enter -1 to redisplay line list) [quit]:  ",i);
  fgets(line,MAXC,stdin);

  while(line[0] != '\n' && line[0] != 'q' && line[0] != 'Q') {

    /*
     * Clear error flag and reset offset
     */

    flag = 0;
    dzl = 0.0;

    /*
     * Check that input is a number
     */

    cptr = line;
    while(*cptr != '\n' && *cptr != '\0') {
      if(*cptr == ' ' || *cptr == '\t')
	cptr++;
      else
	break;
    }
    if(*cptr == 'q' || *cptr == 'Q')
      break;
    else if(*cptr == '-' && *(cptr+1) == '1') {
      print_line_list(specinfo,goodlines);
      flag = 1;
    }
    else if(*cptr < '0' || *cptr > '9') {
      fprintf(stderr,"Invalid input.  ");
      fprintf(stderr,"Enter a number from the list above.\n\n");
      flag = 1;
    }
    else if (sscanf(line,"%d %f",&choice,&dzl) > 0) {

      /*
       * Check to see if value is in acceptable range
       */
      
      if(choice < 0 || choice > goodlines+1) {
	fprintf(stderr,"ERROR: Choice is outside valid range\n\n");
	flag = 1;
      }
    }
    else
      flag = 1;

    if(!flag) {
      if(choice == 0) {
	if(get_sky_label(&marklamb,&label))
	  fprintf(stderr,"No label added.\n");
      }
      else if (choice == goodlines+1) {
	printf("Enter name of line as you want it to appear:  ");
	gets(line);
	strcpy(label.text,line);
	printf("Enter wavelength of feature:  ");
	fgets(line,MAXC,stdin);
	if(sscanf(line,"%f",&marklamb) == 1) {
	  label.x = (1+z) * marklamb;
	}
      }
      else {
	i++;
	sptr = specinfo + choice - 1;
	printf(" You have chosen %s %7.2f.\n",sptr->ionname,sptr->lambda);
	printf(" At z = %7.4f, the line is at %7.2f Angstroms\n\n",
	       z,(1+z)*sptr->lambda);
	marklamb = (1+z) * sptr->lambda;
	label.x = marklamb + dzl;
	strcpy(label.text,sptr->label);
      }

      /*
       * Label the line.
       */

      if(speclab(marklamb,label,spectrum,setup)) {
	fprintf(stderr,"ERROR: label_lines\n");
	fprintf(stderr,"ERROR: Line is probably out of range\n");
	i--;
	return 1;
      }

      /*
       * Add the line to the output list.
       */

      fprintf(ofp,"%8.3f %8.3f %s\n",marklamb,label.x,label.text);
    }

    printf("Choose line #%d (Enter -1 to redisplay line list) [quit]:  ",i);
    fgets(line,MAXC,stdin);
  }

  /*
   * Clean up and exit.
   */

  if(ofp)
    fclose(ofp);

  return 0;
}

/*.......................................................................
 *
 * Function label_lines_file
 *
 * Reads list of lines which were chosen.
 *
 * Inputs: Spectrum *spectrum  spectrum and associated info
 *         Spsetup *setup      container for plotting parameters
 *			       
 * Output: int (0 or 1)        0 ==> success.
 *
 * v13Sep01 CDF, Moved the passed labstyle and flchoice parameters into
 *                the Spsetup structure.
 */

int label_lines_file(Spectrum *spectrum, Spsetup *setup)
{
  float marklamb;       /* x position of marker */
  char foostring[MAXC];
  char *cptr;           /* Pointer for navigating input string */
  char line[MAXC];      /* General string for reading input. */
  Labinfo label;        /* Line label */
  FILE *ifp=NULL;       /* Input file pointer */

  /*
   * Open file -- for now hardwired to foo.LIST.
   */

  if(!(ifp = open_readfile("foo.LIST"))) {
    fprintf(stderr,"ERROR: label_lines_file.\n");
    return 1;
  }

  /*
   * Read input.
   */

  while(fgets(line,MAXC,ifp) != NULL) {
    if(sscanf(line,"%f %f %s",&marklamb,&label.x,foostring) != 3)
      fprintf(stderr,"ERROR: bad input on line %s",line);
    else {
      cptr = strstr(line,foostring);
      strcpy(label.text,cptr);
      label.text[strlen(label.text)-1] = '\0';
      if(speclab(marklamb,label,spectrum,setup)) {
	fprintf(stderr,"\nERROR: label_lines_file.  ");
	fprintf(stderr,"Problem with line at %8.2f Ang.\n",marklamb);
      }
    }
  }

  /*
   * Clean up and exit.
   */

  if(ifp)
    fclose(ifp);

  return 0;

}

/*......................................................................
 *
 * Function print_line_list
 *
 * Prints list of spectral lines and wavelengths
 *
 * Inputs:  Specinfo *specinfo,   structure array containing line info
 *          int goodlines,        number of valid data lines in array
 *
 * Output:  none
 *
 * v25Apr99 CDF, Added the previously "hidden" choices for a source line not
 *                 on the list and for sky lines.
 */

void print_line_list(Specinfo *specinfo, int goodlines)
{
  int i;                /* Looping variable */
  int half;             /* Half of goodlines */
  int isodd=0;          /* Flag set to 1 if goodlines is odd */
  Specinfo *sptr;       /* Pointer to navigate specinfo */

  half = goodlines / 2;
  if((goodlines % 2) == 1)
    isodd = 1;

  printf("\n\n");
  printf("                      SOURCE LINE SELECTION\n");
  printf("----------------------------------------------------------\n\n");
  printf("Choose the lines to plot from the following list.\n");
  printf("Enter one number at a time.  Just hit return to exit from ");
  printf("the line selection.\n\n");
  
  for(i=0,sptr=specinfo; i<half; i++,sptr++) {
    printf("%3d. %s %7.2f       %3d. %s %7.2f\n",
	   i+1,sptr->ionname,sptr->lambda,i+half+1+isodd,
	   (sptr+half+isodd)->ionname,(sptr+half+isodd)->lambda);
  }
  if(isodd)
    printf("%3d. %s %7.2f\n",i+1,sptr->ionname,sptr->lambda);

  printf("\n*** Special choices ***\n");
  printf("  0. Sky line\n");
  printf("%3d. Source line not on list\n",goodlines+1);
  
  printf("\n\n");

}

/*......................................................................
 *
 * Function speclab
 *
 * Labels lines in a spectrum
 *
 * Inputs: float marklamb      line mark position
 *         Labinfo label       label position and text
 *         Spectrum *spectrum  spectrum and associated info
 *         Spsetup *setup      container for plotting parameters
 *
 * Output:  int (0 or 1)       0 ==> success, 1 ==> failure
 *
 * v13Sep01 CDF, Moved the passed labstyle and flchoice parameters into
 *                the Spsetup structure.  Added flexibility in choosing
 *                the label character height (previously hardwired to 1.0)
 *                and vertical location of the labels for labstyle==2
 *                (previously hardwired to 0.06) via parameters contained
 *                in the Spsetup container.
 */

int speclab(float marklamb, Labinfo label, Spectrum *spectrum, Spsetup *setup)
{
  int i;
  int shortcount=0;
  int longcount=0;
  int plotchoice=2;            /* Sets direction of label */
  float zlamb;                 /* Redshifted wavelength of line */
  float shortave=0.0;          /* To get flux of line */
  float longave=0.0;           /* To approximate continuum */
  float fmax,fmin;             /* To get min and max fluxes in shortave */
  float fluxdiff0;             /* Difference btwn minflux and maxflux */
  float fluxdiff;              /* Scaled version of fluxdiff0 */
  Pos *spptr;                  /* Pointer to navigate spectrum->spec */

  /*
   * Set default plot parameters.
   */

  cpgsci(2);
  cpgsch(setup->labheight);
  fluxdiff0 = spectrum->maxflux - spectrum->minflux;
  fluxdiff = 0.625*fluxdiff0;

  /*
   * Set the redshifted wavelength
   */

  zlamb = marklamb;

  /*
   * First option is followed if tick marks are desired (labstyle == 1)
   */

  if(setup->labstyle == 1) {

    /*
     * Find the line and continuum average values
     */

    spptr = spectrum->spec;
    fmax = 0.0;
    fmin = 1.0e6;
    for(i=0; i<spectrum->nlines; i++,spptr++) {
      if(abs(spptr->x - zlamb) > LONGIN && 
	 abs(spptr->x - zlamb) < LONGOUT) {
	longcount++;
	longave += spptr->y;
      }
      else if(abs(spptr->x - zlamb) < SHORT) {
	shortcount++;
	shortave += spptr->y;
	if(spptr->y < fmin)
	  fmin = spptr->y;
	if(spptr->y > fmax)
	  fmax = spptr->y;
      }
    }

    if(shortcount == 0 || longcount == 0) {
      fprintf(stderr,"\nERROR: speclab.\n");
      fprintf(stderr,"  No points found to match averaging criteria.\n");
      return 1;
    }

    shortave = shortave/shortcount;
    longave = longave/longcount;

    /*
     * Now determine which way the tick mark should go, based on the
     *  values for the line and the continuum.  Also take into account
     *  the borders of the plot and the units being used (magnitudes or
     *  fluxes).
     */

    if((shortave < longave && (fmin - spectrum->minflux) > (0.14*fluxdiff))
       || (spectrum->maxflux - fmax) < (0.14*fluxdiff) ) {
      if(setup->flchoice == 1)
	plotchoice = 3;
      else
	plotchoice = 1;
    }
    else {
      if(setup->flchoice == 1)
	plotchoice = 4;
      else
	plotchoice = 2;
    }

    /*
     * Do the labeling
     */

    switch(plotchoice) {
    case 1:
      cpgmove(zlamb,fmin-0.02*fluxdiff);
      cpgdraw(zlamb,fmin-0.10*fluxdiff);
      cpgptxt(label.x,fmin-0.12*fluxdiff,90.0,1.0,label.text);
      break;
    case 2:
      cpgmove(zlamb,fmax+0.02*fluxdiff);
      cpgdraw(zlamb,fmax+0.10*fluxdiff);
      cpgptxt(label.x,fmax+0.12*fluxdiff,90.0,0.0,label.text);
      break;
    case 3:
      cpgmove(zlamb,fmin-0.02*fluxdiff);
      cpgdraw(zlamb,fmin-0.10*fluxdiff);
      cpgptxt(label.x,fmin-0.12*fluxdiff,90.0,0.0,label.text);
      break;
    case 4:
      cpgmove(zlamb,fmax+0.02*fluxdiff);
      cpgdraw(zlamb,fmax+0.10*fluxdiff);
      cpgptxt(label.x,fmax+0.12*fluxdiff,90.0,1.0,label.text);
      break;
    default:
      cpgmove(zlamb,fmax+0.02*fluxdiff);
      cpgdraw(zlamb,fmax+0.10*fluxdiff);
      cpgptxt(label.x,fmax+0.12*fluxdiff,90.0,0.0,label.text);
    }
  }

  /*
   * Draw dotted lines running full height of spectrum if labstyle == 2.
   * Put labels at bottom of page.
   */

  else if(setup->labstyle == 2) {
    cpgsls(4);
    cpgmove(zlamb,spectrum->minflux);
    cpgdraw(zlamb,spectrum->maxflux);
    cpgptxt(label.x,spectrum->minflux+(setup->labyfrac * fluxdiff0),90.0,0.0,
	    label.text);
    cpgsls(1);
  }

  /*
   * If a bad labstyle slipped through the error checking, just return.
   */

  else {
    fprintf(stderr,"ERROR: speclab -- Bad choice for label style.\n");
  }

  cpgsci(1);
  cpgsch(1.0);

  return 0;
}

/*.......................................................................
 *
 * Function get_sky_label
 *
 * Gets label position for a sky line.
 *
 * Inputs:  float *marklamb    position of marker (set by this function).
 *          Labinfo *label     position and text for label (also set here).
 *
 * Output:int (0 or 1)       0 ==> success, 1 ==> failure
 *
 */

int get_sky_label(float *marklamb, Labinfo *label)
{
  char line[MAXC];       /* General string container for input */

  /*
   * Get wavelength of sky line.
   */

  printf("Enter wavelength of sky feature:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f",marklamb) != 1) {
    fprintf(stderr,"\nError in input.  Enter wavelength again: ");
    fgets(line,MAXC,stdin);
  }

  /*
   * Set label text to be the Earth symbol and the label position to
   *  center the symbol above the line.
   */

  strcpy(label->text,"\\(2284)");
  label->x = *marklamb + 30.0;

  return 0;
}

/*......................................................................
 *
 * Function label_plot
 *
 * Puts x- and y-axis labels and title on the plot.
 *
 * Inputs:  char *title        title for the plot
 *          int titleopt       tells where to put the title
 *          int flchoice       flag set to 1 for fluxes in magnitudes
 *          float cheight      character height for labels
 *
 * Output:  int (0 or 1)       0 ==> success, 1 ==> failure
 *
 */

int label_plot(char *title, int titleopt, int flchoice, float cheight)
{
  /*
   * Set plot parameters
   */


  cpgsci(1);
  cpgsch(cheight);
  cpgvstd();
  cpgslw(2);

  /*
   * If the plot title is desired inside the plot (titleopt == 2) then
   *  place it inside with a pgmtxt command and then clear the title
   *  string.
   */

  if(titleopt == 2) {
    cpgmtxt("T",-1.0,0.5,0.5,title);
    title[0] = '\0';
  }

#if 0
  cpglab("Observed Wavelength (\\A)",
	 "F\\d \\gl \\u erg cm\\u -2\\d sec\\u -1 \\d \\A\\u -1\\d",title);
#endif
  switch(flchoice) {
  case 1:
    cpglab("Observed Wavelength (\\A)","Mag",title);
    break;
  case 2:
    cpglab("Observed Wavelength (\\A)","AB Mag",title);
    break;
  default:
    cpglab("Observed Wavelength (\\A)","Relative flux",title);
    break;
  }
#if 0
  cpgsci(0);
  cpglab("","",title);
  cpgsci(1);
  cpglab("","",title);
#endif
  cpgslw(1);

  return 0;
}

/*.......................................................................
 *
 * Function labcorner
 *
 * Puts a label in the lower left right corner of the plot
 *
 * Inputs: char *labtext       label text
 *         char side           side of plot in which to place label
 *         int option          method of plotting spectra
 *
 * Output: none
 *
 */

void labcorner(char *labtext, char side, int option)
{
  float vshift;  /* Amount to shift label from bottom axis in char height */

  if(option == MULTISPEC) {
    cpgsch(1.0);
    vshift = -0.5;
  }
  else {
    cpgsch(1.5);
    vshift = -1.0;
  }

  cpgscf(2);
  if(side == 'r')
    cpgmtxt("B",vshift,0.95,1.0,labtext);
  else
    cpgmtxt("B",vshift,0.03,0.0,labtext);
  cpgsch(1.0);
}

/*.......................................................................
 *
 * Function calc_d4000
 *
 * Calculates the strength of the 4000 Angstrom break (D_4000), which
 *  is defined as the ratio of average fluxes above and below the
 *  break, as:
 *
 *            S_ave (4050 - 4250)
 *   D_4000 = -------------------
 *            S_ave (3750 - 3950)
 *
 * (Bruzual 1983)
 *
 * Inputs: Pos *spectrum       array of wavelengths and fluxes
 *         int npoints         number of sampled point in the spectrum
 *
 * Output:  int (0 or 1)       0 ==> success, 1 ==> error
 *
 */

int calc_d4000(Pos *spectrum, int npoints)
{
  int i;             /* Looping variable */
  int nlsum=0;       /* Temporary variable used in calculating averages */
  int nhsum=0;       /* Temporary variable used in calculating averages */
  float z;           /* Redshift of object in question */
  float rlamb;       /* Restframe wavelength */
  float lsum=0.0;    /* Temporary variable used in calculating averages */
  float hsum=0.0;    /* Temporary variable used in calculating averages */
  float loave;       /* Average of the lower portion */
  float hiave;       /* Average of the higher portion */
  char line[MAXC];   /* General string variable for reading input lines */
  Pos *spptr;        /* Pointer for navigating spectrum */

  /*
   * Get redshift of system so that location of observed wavelengths can
   *  be calculated.
   */

  printf("\nEnter the redshift of the system with the 4000 Ang. break: ");
  gets(line);
  while(sscanf(line,"%f",&z) != 1 || z < 0.0) {
    fprintf(stderr,"ERROR: bad input.  Enter redshift again:  ");
    gets(line);
  }

  /*
   * Print out wavelengh ranges as an error check
   */

  printf("The Balmer break strength will be calculated using average fluxes");
  printf("\n in the ranges %5.0f -- %5.0f and %5.0f -- %5.0f for a ",
	 3750.0*(1+z),3950.0*(1+z),4050*(1+z),4250*(1+z));
  printf("redshift of %6.4f\n\n",z);

  /*
   * Sum up fluxes in the two wavelength ranges
   */

  for(i=0,spptr=spectrum; i<npoints; i++,spptr++) {
    rlamb = spptr->x / (1+z);
    if(rlamb > 3750.0 && rlamb < 3950.0) {
      lsum += spptr->y;
      nlsum++;
    }
    else if(rlamb > 4050.0 && rlamb < 4250.0) {
      hsum += spptr->y;
      nhsum++;
    }
  }

  /*
   * Calculate average fluxes and print out results
   */

  loave = lsum / nlsum;
  hiave = hsum / nhsum;
  printf("There were %d flux points in the upper range\n",nhsum);
  printf("and %d points in the lower range\n\n",nlsum);
  printf("Average flux above the break is %9.2e\n",hiave);
  printf("Average flux below the break is %9.2e\n",loave);
  printf("D_4000 = %6.2f\n\n",hiave/loave);

  return 0;
}

/*.......................................................................
 *
 * Function flambda_to_fnu
 *
 * Converts a spectrum from f_lambda units to f_nu units.
 *
 * Inputs: Pos *spectrum       input spectrum -- gets modified by this 
 *                              function.
 *         int npoints         number of points in the spectrum
 *         float *minflux      minimum flux in rescaled spectrum
 *         float *maxflux      maximum flux in rescaled spectrum
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int flambda_to_fnu(Pos *spectrum, int npoints, float *minflux,
		   float *maxflux)
{
  int i;               /* Looping variable */
  float dlambda;       /* Spacing of input array (assumed constant) */
  float dnu;           /* Width of bin in frequency */
  float newmin,newmax; /* Min and max in rescaled spectrum */
  float fluxdiff;      /* Final difference between newmax and newmin */
  Pos *sptr;           /* Pointer to navigate spectrum */

  /*
   * First calculate bin width in lambda (assumed to be constant)
   */

  printf("\nflambda_to_fnu: Converting to f_nu\n");
  dlambda = (spectrum+1)->x - spectrum->x;
  printf("flambda_to_fnu: dlambda = %f Angstroms.\n",dlambda);

  /*
   * Now step through spectrum converting to f_nu by using 
   *  f_nu dnu = f_lambda dlambda.  Note that f_nu has to be calculated at each
   *  point in the spectrum via:  dnu = c/lambda - c/(lambda+dlambda).
   *  Also note that to do this conversion, the speed of light is expressed in
   *  angstroms/sec.
   */

  newmin = 1.0e6;
  newmax = -1.0e6;
  for(i=0,sptr=spectrum; i<npoints; i++,sptr++) {
    dnu = (3.0e18/sptr->x) - (3.0e18/(sptr->x + dlambda));
    sptr->y *= (dlambda / dnu);
    if(sptr->y > newmax)
      newmax = sptr->y;
    if(sptr->y < newmin)
      newmin = sptr->y;
  }

  /*
   * Reset minflux and maxflux
   */

  fluxdiff = newmax - newmin;
  *minflux = newmin - 0.35*fluxdiff;
  *maxflux = newmax + 0.35*fluxdiff;

  return 0;
}

