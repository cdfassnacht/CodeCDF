/*
 * modfuncs.c
 *
 * A collection of functions for use in scaling difmap models
 *
 * 20Mar97 CDF
 * v28Aug97 CDF, Modified print_compmod and modfile_08 to do amplitude
 *                selfcal on 15 min timescales.
 * v29Aug97 CDF, Took out amplitude selfcal.
 * v30Nov97 CDF, Modified print_compmod to do a modelfit on central
 *                component in 1635.
 * v10Feb98 CDF, Updated location of model files to be consistent with
 *                new directory structure.
 * v23May98 CDF, Put model directory name into a #define constant
 *               Changed difmap_headers to average *.uvf data files to 10 sec
 *                integrations.
 *               Changed modfile_08 to write out both phasecal and gscale 
 *                models RMS values.
 * v25May98 CDF, Added a modtype parameter to source_flux function
 * v26May98 CDF, Added uvstat(chisq) function to modfile_08 to print out
 *                chisq of final model from model fitting.
 * v01Jun98 CDF, Modified source_flux to read in 1608 chisq value.
 *               Changed print_compmod to do uvstat on 1634 and 1635 too.
 * v03Jun98 CDF, Modified difmap_headers and print_compmod to read in
 *                the window around the central component in 1635, delete
 *                the clean components in it, and do a clean again in
 *                that window.  This is to try to compensate for possible
 *                variations in the central component.
 * v20Jun98 CDF, Took out separate fitting of 1635 central component, for now.
 * v22Mar00 CDF, Made source_flux more flexible so that it can handle
 *                calibrator outputs that are produced by difmap model
 *                fitting as well as the 1608 outputs.
 * v15Sep00 CDF, Changed source_flux to take an additional passed parameter,
 *                which is the date of observation in the form yyyymmdd.
 *               Added the set_log_obsdate, which takes the observing date
 *                in JD that has been read from the model file and converts
 *                it to MJD and then rounds to the nearest half-day.
 *               Cleaned up get_sinfo and get_clean_flux.
 *               Added new functions print_lensflux and print_calflux to
 *                print out data from the new combination of getflux.c and
 *                calflux.c (called getflux.c). 
 * v04May01 CDF, Modified scale_data and related functions to match updated 
 *                modscal.c and automapping scripts.
 * v09May01 CDF, Modified cal_flux, read_modfit, and get_best_cal functions
 *                to match updated best_cmod.c.
 * v21May01 CDF, Fixed bug in finding 3 lowest points in best_model.
 * v02Apr02 CDF, Took out rounding to nearest half day in set_log_obsdate
 * v14Jul03 CDF, Made source_flux more general to handle lenses other than
 *                1608.
 * v07Oct05 CDF, Made source_flux even more general.
 *               Adjusted print_lensflux and print_calflux to take account
 *                of the changes resulting from the updated source_flux
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "libfits.h"
#include "slalib.h"
#include "structdef.h"
#include "coords.h"
#include "dataio.h"
#include "modfuncs.h"

#define SHIFTDATE 50528.0   /* Date at which pointing for 1635 changed */
#define DOMODEL 1

/*.......................................................................
 *
 * Function new_cinfo
 *
 * Allocates dynamic memory for a pointer array of Compinfo structures
 *
 * Input:  int size            size of the array
 *
 * Output: Compinfo *newinfo   pointer to the new array.  NULL if error
 *
 */


Compinfo *new_cinfo(int size)
{
  Compinfo *newinfo;
 
  newinfo = (Compinfo *) malloc(sizeof(Compinfo) * size);
  if(!newinfo) {
    fprintf(stderr,"Insufficient memory for data array.\n");
    return NULL;
  }

  return newinfo;
}

/*.......................................................................
 *
 * Function del_cinfo
 *
 * Frees up memory allocated to float array
 *
 * Input:  float *cinfo        array to be freed
 *
 * Output: NULL
 */

Compinfo *del_cinfo(Compinfo *cinfo)
{
  if(cinfo)
    free(cinfo);
 
  return NULL;
}

/*.......................................................................
 *
 * Function get_ext
 *
 * Finds the sourcefile extension if the program is in interactive mode
 *
 * Inputs: char *ext           extension (set by this function)
 *
 * Output: void
 *
 */

void get_ext(char *ext)
{
  int extchoice;
  char line[MAXC];

  printf("Choose the filename extension\n");
  printf("  1. edt\n");
  printf("  2. incom\n");
  printf("Enter choice: [1] ");
  fgets(line,MAXC,stdin);
  if(line[0] != '\n')
    extchoice = 1;
  else if((sscanf(line,"%d",&extchoice)) != 1) {
    printf("Error in format -- assuming that choice is edt\n");
    extchoice = 1;
  }
  switch(extchoice) {
  case 1:
    strcpy(ext,"edt");
    break;
  case 2:
    strcpy(ext,"incom");
    break;
  default:
    printf("Not a valid choice -- setting choice to edt\n");
    strcpy(ext,"edt");
  }
}

/*.......................................................................
 *
 * Function default_caltimes
 *
 * Sets up the default phase self-cal time scales for each of the three
 *  sources.
 *
 * Inputs: char *root          source root name
 *         float *caltime      self-cal timescale (set by this function)
 *
 * Output: NONE
 *
 */

void default_caltimes(char *root, float *caltime)
{
  if(strcmp(root,"1634") == 0)
    *caltime = 1.0;
  else if(strcmp(root,"1635") == 0)
    *caltime = 1.0;
  else if(strcmp(root,"1608") == 0)
    *caltime = 1.0;
  else
    *caltime = 5.0;
}

/*.......................................................................
 *
 * Function get_modname
 *
 * Uses the source root name, filename extension, and date of observation
 *  to determine the name of the model file (in the model directory MODDIR)
 *  that will be used.  The filename is returned.
 *
 * Inputs: char *modname       model filename (set by function)
 *         char *root          root name
 *         char *ext           filename extention
 *         float obsdate       MJD observation date
 *
 * Output: None
 *
 */

void get_modname(char *modname, char *root, char *ext, float obsdate)
{
  if(strcmp(ext,"edt") == 0 || strcmp(root,"1608") == 0)
    sprintf(modname,"%s%s_a_and_b",MODDIR,root);
  else
    sprintf(modname,"%s%s_%s",MODDIR,root,ext);

  if(strcmp(root,"1635") == 0 && obsdate >= SHIFTDATE)
    strcat(modname,"_shift.mod");
  else
    strcat(modname,".mod");
}


/*.......................................................................
 *
 * Function open_modfile
 *
 * Checks for existence of a model file and, if the file exists opens it.
 *  Returns the file pointer.
 *
 * Inputs: char *modname       model file name
 *
 * Output: FILE *ifp           pointer to model file, NULL on error
 *
 */

FILE *open_modfile(char *modname)
{
  FILE *ifp;

  /*
   * Try to open file
   */

  if((ifp = fopen(modname,"r")) == NULL) {
    fprintf(stderr,"Error.  Cannot open %s\n",modname);
    fprintf(stderr,"Please check your root name and model file\n");
    return NULL;
  }
  else {
    printf("\nModel file %s now open.\n\n",modname);
  }

  return ifp;
}

/*.......................................................................
 *
 * Function read_cinfo
 *
 * Reads in data and puts it in a Compinfo array
 *
 * Inputs: Compinfo *cinfo     array of component info to be filled
 *         FILE *ifp           input file
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int read_cinfo(Compinfo *cinfo, FILE *ifp)
{
  Compinfo *ciptr;     /* Pointer to navigate cinfo array */
  char line[MAXC];     /* String variable for reading input */

  /*
   * First set temporary pointer to beginning of array
   */

  ciptr = cinfo;

  /*
   * Now read in the data
   */

  while(fgets(line,MAXC,ifp) != NULL) {

    if(line[0] != '!') {

      /*
       * Get the information out of the input line
       */

      switch(sscanf(line,"%f %f %f %f %f %f %d",&ciptr->flux,&ciptr->radius,
		    &ciptr->theta,&ciptr->maj,&ciptr->axrat,&ciptr->phi,
		    &ciptr->type)) {
      case 7:
	break;
      case 3:
	ciptr->maj = ciptr->axrat = ciptr->phi = 0.0;
	ciptr->type = 0;
	break;
      default:
	fprintf(stderr,"ERROR: read_cinfo.  Improper input format\n");
	return 1;
      }

      /*
       * Now increment the temporary pointer
       */

      ciptr++;
    }
  }

  return 0;
}

/*.......................................................................
 *
 * Function get_sinfo
 *
 * Gets minimum and maximum scaling and step size between them
 *
 * Inputs: float *min          minimum, set by this function
 *         float *max          maximum, set by this function
 *         float *step         step size, set by this function
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int get_sinfo(float *min, float *max, float *step)
{
  char line[MAXC];    /* General variable for reading input */

  /*
   * Initialize
   */

  *min = 0.9;
  *max = 1.1;
  *step = 0.02;

  /*
   * Read in values
   */

  printf("Enter minimum scale value [%4.1f]:  ",*min);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%f",min) != 1 || *min < 0.0) {
      fprintf(stderr,"get_sinfo: Bad input value.  Enter new value:  ");
      fgets(line,MAXC,stdin);
    }
  }

  printf("Enter maximum scale value [%4.1f]:  ",*max);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%f",max) != 1 || *max < 0.0 || *max <= *min) {
      fprintf(stderr,"get_sinfo: Bad input value.  Enter new value:  ");
      fgets(line,MAXC,stdin);
    }
  }

  printf("Enter step size [%4.1f]:  ",*step);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%f",step) != 1 || *step < 0.0) {
      fprintf(stderr,"get_sinfo: Bad input value.  Enter new value:  ");
      fgets(line,MAXC,stdin);
    }
  }

  return 0;
}

/*.......................................................................
 *
 * Function scale_data
 *
 * Loops through scale factors, multiplying the flux of each clean component
 *  by the current scale factor.  Creates model files containing the scaled
 *  clean components.  Also creates a difmap script file to do model fitting
 *  using the scaled models.
 *
 * Inputs: float min           minimum scale value
 *         float max           maximum scale value
 *         float step          step size for scales
 *         char *basemod       file name of "base" model, i.e. the one that
 *                              gets scaled by this function.
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int scale_data(float min, float max, float step, char *basemod)
{
  int i;                /* Looping variables */
  int no_error=1;       /* Flag set to zero on error */
  int nsteps;           /* Number of steps in scale factors */
  int nlines=0;         /* Number of lines in model input file */
  float totflux;        /* Total flux in input components */
  float scale;          /* Scale factor */
  float obsdate;        /* Date of observation */
  char outname[MAXC];   /* Model output file name */
  char scrname[MAXC];   /* Name of difmap output file */
  Compinfo *cinfo=NULL; /* Model component info */
  FILE *mfp=NULL;       /* Input base model file */
  FILE *scrfp=NULL;     /* Difmap script output file */

  /*
   * Calculate number of steps
   */

  nsteps = floor((max-min)/step) + 1;
  printf("\nThere will be %d model files produced.\n\n",nsteps);

  /*
   * Open base model file
   */

  if(!(mfp = open_readfile(basemod))) {
    fprintf(stderr,"ERROR: scale_data");
    return 1;
  }

  /*
   * Count number of lines in input file 
   */

  nlines = n_lines(mfp,'!');
  if(nlines == 0) {
    fprintf(stderr,"ERROR: scale_data.  No valid data lines in model file\n");
    return 1;
  }
  else {
    printf("\n%d lines in model file\n\n",nlines);
    rewind(mfp);
  }

  /*
   * Fill cinfo with data
   */

  if(!(cinfo = fill_cinfo(mfp,nlines))) {
    fprintf(stderr,"ERROR: scale_data\n");
    return 1;
  }

  /*
   * Calculate total flux in components, to be used in the scaling
   */

  if(calc_totflux(cinfo,nlines,&totflux)) {
    fprintf(stderr,"ERROR: scale_data\n");
    no_error = 0;
  }
#if 0
  /*
   * Open the script file for writing and put in header info
   */

  if(no_error) {
    sprintf(scrname,"%s.scr",root);
    printf("\nWriting difmap script file %s...\n\n",scrname);
    if((scrfp = fopen(scrname,"w")) == NULL) {
      fprintf(stderr,"ERROR: scale_data.  Cannot open %s for output\n",scrname);
      no_error = 0;
    }
  }

  if(no_error)
    difmap_headers(scrfp,root,ext,obsdate,caltime,interact);
#endif
  /*
   * Loop through the scale factors
   */

  for(i=0; i<nsteps; i++) {
    if(no_error) {

      /*
       * Define scale and set up file name for scaled model file
       */

      scale = min + i*step;
      sprintf(outname,"mod_%d",i+1);

      /*
       * Print scaled model components to model file
       */

      if(print_compmod(i+1,scale,totflux,cinfo,nlines,outname)) {
	fprintf(stderr,"ERROR:  scale_data\n");
	no_error = 0;
      }
    }
  }

  /*
   * Clean up
   */

  if(mfp)
    fclose(mfp);
  if(scrfp)
    fclose(scrfp);

  if(no_error)
    return 0;
  else
    return 1;
}

/*.......................................................................
 *
 * Function modfile_08
 *
 * Creates the model fitting file for the 1608 models
 *
 * Inputs: char *ext           UVF filename extension
 *         float caltime       phase self-cal timescale
 *         int interact        0 ==> batch mode
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int modfile_08(char *root, char *ext, float caltime, int interact)
{
  float obsdate;       /* Date of observation */
  char comment[MAXC];  /* Set to "!" if in batch mode */
  char outmod[MAXC];   /* Name of output model file */
  FILE *modfp=NULL;    /* Model file pointer */

  /*
   * First check that appropriate UVF file exists, and if it does, get
   *  the MJD of the observation.
   */

  if(check_uvf(root,ext,&obsdate))
    return 1;

  /*
   * Create the model file
   */

  printf("\nWriting difmap script file 1608.scr ....");
  if((modfp = fopen("1608.scr","w")) == NULL) {
    fprintf(stderr,"\n\nERROR: modfile_08.  ");
    fprintf(stderr,"Cannot open 1608.scr for output\n");
    return 1;
  }

  /*
   * Print lines that are in common for all sources and set interact
   *  string.
   */

  difmap_headers(modfp,root,ext,obsdate,caltime,interact);
  if(interact)
    strcpy(comment,"");
  else
    strcpy(comment,"!");

  /*
   * Print lines for model fitting or cleaning and phase selfcal
   */

  if(DOMODEL)
    fprintf(modfp,"modelfit 7\n");
  else {
    fprintf(modfp,"clrmod true\n");
    fprintf(modfp,"rmod %s1608_a_and_b.no_abcd.mod\n",MODDIR);
    fprintf(modfp,"rwin %s1608_a_and_b.peak.win\n",MODDIR);
    fprintf(modfp,"clean 200\n");
    fprintf(modfp,"keep\n");
    fprintf(modfp,"clean\n");
  }
  fprintf(modfp,"%smapl\n",comment);
  if(caltime == 0.0)
    fprintf(modfp,"selfcal\n");
  else
    fprintf(modfp,"selfcal false,false,%f\n",caltime);
  fprintf(modfp,"%smapl\n",comment);
  if(DOMODEL)
    fprintf(modfp,"modelfit 6\n");
  else
    fprintf(modfp,"clean 200\n");
  fprintf(modfp,"print \"Image rms phasecal = \", imstat(rms)\n");
  fprintf(modfp,"print \"UV plane chisq phasecal = \",uvstat(chisq)\n");

  /*
   * Print lines to save the pre-gscale model
   */

  if(strcmp(ext,"edt") == 0) 
    sprintf(outmod,"1608_phasecal.mod");
  else
    sprintf(outmod,"1608_%s_phasecal.mod\n",ext);
  fprintf(modfp,"wmod %s\n",outmod);
  if(!DOMODEL) {
    fprintf(modfp,"delwin\n");
    fprintf(modfp,"rwin %s1608_abcd.win\n",MODDIR);
    fprintf(modfp,"winmod\n");
    fprintf(modfp,"wmod 1608_abcd_phasecal.mod\n");
    fprintf(modfp,"rmod %s\n",outmod);
    fprintf(modfp,"rwin %s1608_a_and_b.peak.win\n",MODDIR);
  }


  /*
   * Print lines for gscale and model fitting
   */

  fprintf(modfp,"gscale\n");
  if(DOMODEL)
    fprintf(modfp,"modelfit 5\n");
  else
    fprintf(modfp,"clean 200\n");
  if(caltime == 0.0)
    fprintf(modfp,"selfcal\n");
  else
    fprintf(modfp,"selfcal false,false,%f\n",caltime);
  if(DOMODEL)
    fprintf(modfp,"modelfit 5\n");
  else
    fprintf(modfp,"clean 200\n");
  fprintf(modfp,"%smapl\n",comment);
  fprintf(modfp,"print \"Image rms gscale = \", imstat(rms)\n");
  fprintf(modfp,"print \"UV plane chisq gscale = \",uvstat(chisq)\n");

  /*
   * Print lines to save the post-gscale model
   */

  if(strcmp(ext,"edt") == 0)
    fprintf(modfp,"wmod 1608_gscale.mod\n");
  else
    fprintf(modfp,"wmod 1608_%s_gscale.mod\n",ext);
  if(!DOMODEL) {
    fprintf(modfp,"delwin\n");
    fprintf(modfp,"rwin %s1608_abcd.win\n",MODDIR);
    fprintf(modfp,"winmod\n");
    fprintf(modfp,"wmod 1608_abcd_gscale.mod\n");
  }

  printf("\n\n");

  if(modfp)
    fclose(modfp);

  return 0;
}

/*.......................................................................
 *
 * Function calc_totflux
 *
 * Sums up the total flux in all the components in the Compinfo array.
 *
 * Inputs: Compinfo *cinfo     component information
 *         int nlines          number of components
 *         float *totflux      total flux -- set by this function
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int calc_totflux(Compinfo *cinfo, int nlines, float *totflux)
{
  int i;           /* Looping variable */
  float sum;       /* Used to calculate summed flux */
  Compinfo *cptr;  /* Pointer to navigate cinfo */
  
  sum = 0.0;

  for(i=0,cptr=cinfo; i<nlines; i++,cptr++) {
    sum += cptr->flux;
  }

  /*
   * Convert flux densities from Janskys to millijanskys
   */

  *totflux = sum * 1000.0;
  printf("\nTotal flux in clean components = %9.3f mJy\n",*totflux);

  return 0;
}

/*.......................................................................
 *
 * Function print_compmod
 *
 * Given a scale factor, total flux and output name, creates an output
 *  model file for one of the comparison sources (1634 or 1635)
 *
 * Inputs: char *root          root name (1634 or 1635)
 *         int modnum          model number
 *         float scale         model scaling
 *         float totflux       total flux in model clean components
 *         Compinfo *cinfo     array containing model component info
 *         int nlines          number of model components
 *         char *outname       output file name
 *         FILE *scrfp         script (that runs the model fitting) file
 *                              pointer
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int print_compmod(int modnum, float scale, float totflux, 
		  Compinfo *cinfo, int nlines, char *outname)
{
  int i;           /* Looping variable */
  Compinfo *cptr;  /* Pointer for navigating cinfo */
  FILE *ofp=NULL;  /* Output file pointer */

  /*
   * Open output file
   */

  if(!(ofp = open_writefile(outname))) {
    fprintf(stderr,"ERROR: print_compmod.  Cannot open %s for output\n",
	    outname);
    return 1;
  }

  printf("Creating file %s with scaling %6.4f ==> flux = %9.3f mJy\n",
	 outname,scale,scale*totflux);
#if 0
  fprintf(scrfp,"rmod %s\n",outname);
  fprintf(scrfp,"print \"Model %d -- Flux = %9.3f mJy, Scale = %6.4f\"\n",
	  modnum,scale*totflux,scale);
  fprintf(scrfp,"uncal false,true\n");
  fprintf(scrfp,"gscale\n");
  fprintf(scrfp,"modelfit 0\n");
  if(strcmp(root,"1635") == 0) {
    fprintf(scrfp,"winmod true\n");
    fprintf(scrfp,"clean 100\n");
    fprintf(scrfp,"keep\n");
  }
  fprintf(scrfp,"print \"Chisq = \", uvstat(chisq)\n");
  fprintf(scrfp,"print \"Image rms = \", imstat(rms)\n");
#endif

  /*
   * Print scaled model to output file
   */

  for(i=0,cptr=cinfo; i<nlines; i++,cptr++) {
    if(i == 0)
      fprintf(ofp,"%#12.6gv %11.2f %8.3f",scale * cptr->flux, 
	      cptr->radius,cptr->theta);
    else
      fprintf(ofp,"%#12.6g  %11.2f %8.3f",scale * cptr->flux, 
	      cptr->radius,cptr->theta);
    if(cptr->type == 0)
      fprintf(ofp,"\n");
    else
      fprintf(ofp," %8.4f %8.6f %9.4f %d\n",
	      cptr->maj,cptr->axrat,cptr->phi,cptr->type);
  }

  if(ofp)
    fclose(ofp);

  return 0;
}

/*.......................................................................
 *
 * Function difmap_headers
 *
 * Prints out header information that is common to all scripts.
 *
 * Input:  FILE *scrfp         output file pointer
 *         char *root          source root name
 *         char *ext           filename extension
 *         float obsdate       date of observation
 *         float caltime       phase self-cal time interval
 *         int interact        0 ==> batch mode
 *
 * Output: void
 *
 */

void difmap_headers(FILE *scrfp, char *root, char *ext, float obsdate, 
		    float caltime, int interact)
{
  int doshift=0;       /* Flag set to 1 to shift phase center of 1635 map */
  int mapsize;         /* Size of map in pixels */
  float pixsize;       /* Size of pixels in arcsec */
  char modname[MAXC];  /* Model file name */
  char shortext[MAXC]; /* "Short" file extension */
  char comment[MAXC];  /* Used to comment out lines in batch mode */

  if(!interact)
    strcpy(comment,"!");
  else
    strcpy(comment,"");
  setvars(root,ext,obsdate,shortext,&mapsize,&pixsize,&doshift);
  fprintf(scrfp,"logfile %s%s.log\n",root,shortext);
  fprintf(scrfp,"obs %s_%s.uvf,10\n",root,ext);
  fprintf(scrfp,"select I\n");
  fprintf(scrfp,"rflags=\"m3\"\n");
  fprintf(scrfp,"mapunits arcsec\n");
  fprintf(scrfp,"uvwei 0,-1\n");
  fprintf(scrfp,"wtscal 11053.0\n");
  fprintf(scrfp,"float scale\n");
  fprintf(scrfp,"header\n");
  fprintf(scrfp,"%sdev /xw\n",comment);
  fprintf(scrfp,"%sradpl\n",comment);
  fprintf(scrfp,"%stpl\n",comment);
  if(doshift)
    fprintf(scrfp,"shift -7.3,0.1\n");
  fprintf(scrfp,"mapsize %d,%5.3f\n",mapsize,pixsize);
  fprintf(scrfp,"%smapl\n",comment);
  get_modname(modname,root,ext,obsdate);
  fprintf(scrfp,"rmod %s\n",modname);
  if(strcmp(root,"1608") == 0) {
    if(caltime == 0.0) {
      fprintf(scrfp,"selfcal\n");
    }
    else {
      fprintf(scrfp,"selfcal false,false,%3.1f\n",caltime);
    }
  }
  else
    fprintf(scrfp,"selfcal\n");
  fprintf(scrfp,"%sradpl\n",comment);
  fprintf(scrfp,"%smapl\n",comment);

}

/*.......................................................................
 *
 * Function setvars
 *
 * Fixes the values of the difmap parameters that are dependent on the
 *  source and/or model.
 *
 * Inputs: char *root          source root name
 *         char *ext           filename extension
 *         float obsdate       date of observation
 *         char *shortext      "short" filename extension
 *         int *mapsize        size of map in pixels
 *         float *pixsize      size of pixels in arcsec
 *         int *doshift        if true, source=1635 and date is before 21Mar97
 *
 * Output: none
 *
 */

void setvars(char *root, char *ext, float obsdate, char *shortext, 
	     int *mapsize, float *pixsize, int *doshift)
{
  float pixscal=1.0;    /* Pixel scaling factor */

  /*
   * Determine source-specific parameters
   */

  *doshift = 0;
  if(strcmp(root,"1635") == 0) {
    *mapsize = 512;
    if(obsdate < SHIFTDATE)
      *doshift = 1;
  }
  else {
    *mapsize = 256;
  }

  /*
   * Determine short filename extension
   */

  if(strcmp(ext,"edt") == 0) {
    strcpy(shortext,"");
  }
  else if(strcmp(ext,"incom") == 0) {
    strcpy(shortext,"_incom");
    pixscal = 4.0;
    *mapsize /= 2;
  }
  else
    sprintf(shortext,"_%s",ext);

  *pixsize = 0.05 * pixscal;

}
    
/*.......................................................................
 *
 * Function check_uvf
 *
 * Constructs filename of UVF file associated with the model files,
 *  and then calls get_obsdate to check existence of UVF file and
 *  get the observing date from the FITS header.
 *
 * Inputs: char *root          source root name
 *         char *ext           filename extension
 *         float *obsdate      set by function
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int check_uvf(char *root, char *ext, float *obsdate)
{
  char filename[MAXC];    /* UVF filename */

  sprintf(filename,"%s_%s.uvf",root,ext);

  if(get_obsdate(filename,obsdate)) {
    fprintf(stderr,"ERROR: check_uvf\n");
    return 1;
  }

  return 0;
}

/*.......................................................................
 *
 * Function get_obsdate
 *
 * Given a UVF filename, checks the existence of the file and if it
 *  exists, gets the observation date (in Julian days) from
 *  the FITS file header.
 *
 * Inputs: char *filename      UVF file name
 *         float *obsdate      observation date (set by this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int get_obsdate(char *filename, float *obsdate)
{
  int no_error=1;   /* Flag set to 0 on error */
  int yr,month,day; /* Date of observation converted to integers */
  int calflag;      /* Return flag of caldj */
  double mjd;       /* Modified Julian Date of observation */
  Fits *uvfits;     /* UVF file */
  Phdu *phdu;       /* Primary HDU of FITS file */

  if(!(uvfits = new_Fits(filename,1,1,0,0))) {
    fprintf(stderr,"ERROR: get_obsdate. %s does not exist\n",filename);
    return 1;
  }

  phdu = (Phdu *) uvfits->hdu;

  printf("\nThe date is %s\n",phdu->date_obs);

  /*
   * Convert to MJD
   */

  if(sscanf(phdu->date_obs,"%d/%d/%d",&day,&month,&yr) != 3 &&
     sscanf(phdu->date_obs,"%d-%d-%d",&yr,&month,&day) != 3) {
    fprintf(stderr,"ERROR: get_obsdate. DATE-OBS not in correct format\n");
    return 1;
  }

  if(yr < 1900)
    yr += 1900;
  slaCldj(yr,month,day,&mjd,&calflag);
  if(calflag) {
    fprintf(stderr,"ERROR: get_obsdate\n");
    no_error = 0;
  }

  if(no_error) {
    printf("MJD = %f\n",mjd);
    *obsdate = mjd;
  }

  uvfits = del_Fits(uvfits);

  if(no_error)
    return 0;
  else
    return 1;
}

/*.......................................................................
 *
 * Function set_log_obsdate
 *
 * Takes an observing date (in JD) as read from a difmap log and converts
 *  it to MJD.  
 *
 * Inputs: double *obsdate     observing date from log -- modified by
 *                              this function.
 *
 * Output: (none)
 *
 * v02Apr02 CDF, Took out rounding to nearest half day.
 */

void set_log_obsdate(double *obsdate)
{
  printf("\nset_log_obsdate: Obsdate from log file is JD = %9.3f\n",
	 *obsdate);

  /*
   * Convert obsdate to MJD
   */

  *obsdate -= 2400000.0;
  printf("set_log_obsdate: Convert obsdate to MJD = %8.3f\n",*obsdate);
#if 0
  *obsdate = 0.5 * (floor(2.0 * *obsdate + 0.5));
  printf("set_log_obsdate: *** Round to nearest half day.  Obsdate = ");
  printf("%7.2f ***\n",*obsdate);
#endif

  /*
   * Finally convert obsdate to 1608 monitoring format.
   */

  *obsdate -= 50000.0;
}

/*.......................................................................
 *
 * Function read_in_params
 *
 * If the optional input file exists, read in paramters from it
 *
 * Inputs: FILE *ifp           input file
 *         char *ext           model filename extension (set by this function)
 *         float *min          min. scale value (set by this function)
 *         float *max          max. scale value (set by this function)
 *         float *step         scale step value (set by this function)
 *         float *caltime      phasecal timescale (set by this function)
 *         int *interact       interactive flag (set by this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int read_in_params(FILE *ifp, char *root, char *ext, float *min, float *max, 
		   float *step, float *caltime, int *interact)
{
  int no_error=1;     /* Flag set to 0 on error */
  char line[MAXC];    /* General variable for reading in line */
  char keywd[MAXC];   /* Key word on line */

  /*
   * Initialize values to defaults
   */

  strcpy(ext,"edt");
  *min = 0.8;
  *max = 1.1;
  *step = 0.02;
  *interact = 1;
  default_caltimes(root,caltime);

  /*
   * Read in data
   */

  while(fgets(line,MAXC,ifp) != NULL && no_error) {
    if(strcmp(line,"") != 0) {
      if(sscanf(line,"%s",keywd) != 1)
	no_error = 0;
    }
    if(strcmp(keywd,"ext") == 0 || strcmp(keywd,"EXT") == 0) {
      if(sscanf(line,"%s %s",keywd,ext) != 2)
	no_error = 0;
    }
    else if(strcmp(keywd,"min") == 0 || strcmp(keywd,"MIN") == 0) {
      if(sscanf(line,"%s %f",keywd,min) != 2)
	no_error = 0;
    }
    else if(strcmp(keywd,"max") == 0 || strcmp(keywd,"MAX") == 0) {
      if(sscanf(line,"%s %f",keywd,max) != 2)
	no_error = 0;
    }
    else if(strcmp(keywd,"step") == 0 || strcmp(keywd,"STEP") == 0) {
      if(sscanf(line,"%s %f",keywd,step) != 2)
	no_error = 0;
    }
    else if(strcmp(keywd,"interact") == 0 || strcmp(keywd,"INTERACT") == 0) {
      if(sscanf(line,"%s %d",keywd,interact) != 2)
	no_error = 0;
    }
    else if(strcmp(keywd,"caltime") == 0 || strcmp(keywd,"CALTIME") == 0) {
      if(sscanf(line,"%s %f",keywd,caltime) != 2)
	no_error = 0;
    }
  }

  if(no_error) {
    printf("\nParameters were partially/fully set by an input file.\n");
    printf("   min  = %f\n",*min);
    printf("   max  = %f\n",*max);
    printf("   step = %f\n",*step);
    printf("   caltime = %f\n",*caltime);
    printf("   interact = %d\n\n",*interact);
    return 0;
  }
  else {
    fprintf(stderr,"ERROR: read_in_params.  Problem with input file\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function fill_cinfo
 *
 * Allocates memory for a Compinfo array and then fills the array
 *  by reading model data file (with the function read_cinfo)
 *
 * Input:  FILE *ifp           file containing model component info
 *         int nlines          number of lines in input file
 *
 * Output: Compinfo *cinfo     filled Compinfo array -- NULL on error
 *
 */

Compinfo *fill_cinfo(FILE *ifp, int nlines)
{
  Compinfo *cinfo=NULL;  /* Compinfo array to be filled */

  /*
   * Allocate array memory
   */

  if(!(cinfo = new_cinfo(nlines))) {
    fprintf(stderr,"ERROR:  fill_cinfo.\n");
    return NULL;
  }

  /*
   * Read in input
   */

  printf("Reading model file...");
  if(read_cinfo(cinfo,ifp)) {
    printf("\n\n");
    fprintf(stderr,"ERROR:  fill_cinfo.\n");
    return del_cinfo(cinfo);
  }
  else
    printf("  Finished\n\n");

  return cinfo;
}

/*.......................................................................
 *
 * Function new_modfit
 *
 * Allocates memory for a Modfit array of size "size".
 *
 * Input:  int size            size of the array
 *
 * Output: Modfit *new_modfit  pointer to the new array.  NULL on error
 *
 */

Modfit *new_modfit(int size)
{
  Modfit *new_modfit;

  new_modfit = (Modfit *) malloc(sizeof(Modfit) * size);
  if(!new_modfit) {
    fprintf(stderr,"Insufficient memory for Modfit array.\n");
    return NULL;
  }

  return new_modfit;
}

/*.......................................................................
 *
 * Function del_modfit
 *
 * Frees up memory associated with the modfit structure array.
 *
 * Input: Modfit *modfit       array to be freed
 *
 * Output: NULL
 *
 */

Modfit *del_modfit(Modfit *modfit)
{
  if(modfit)
    free(modfit);

  return NULL;
}

/*.......................................................................
 *
 * Function open_file
 *
 * Opens files expected by the program.  Returns NULL on error, since the
 *  file names are hard-wired at this point.
 *
 * Input:  char *name          file name of input file
 *
 * Output: FILE *ifp           pointer to successfully opened file.
 *                              NULL on error.
 *
 */

FILE *open_file(char *name)
{
  FILE *ifp;

  if((ifp = fopen(name,"r")) == NULL) {
    fprintf(stderr,"ERROR: open_file.  Cannot open %s\n",name);
    return NULL;
  }
  else
    return ifp;
}

/*.......................................................................
 *
 * Function cal_flux
 *
 * Finds the scale factor associated with the best-fit model for one 
 *  of the comparison sources (1634 or 1635).
 *
 * Inputs: char *filename      input filename
 *         float *scale        best-fit scale factor (set by this function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int cal_flux(char *filename, float *scale)
{
  int no_error=1;       /* Flag set to 0 on error */
  int nlines;           /* Number of lines in input file */
  Modfit *modfit=NULL;  /* Model fitting info */

  /*
   * Get model-fit info
   */

  if(!(modfit = read_modfit(filename,&nlines))) {
    fprintf(stderr,"ERROR: cal_flux\n");
    return 0;
  }

  /*
   * Get best-fit scale factor
   */

  if(no_error)
    if((*scale = best_model(modfit,nlines)) < 0.0)
      no_error = 0;

  /*
   * Clean up
   */

  modfit = del_modfit(modfit);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR:  cal_flux\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function get_best_cal
 *
 * Finds the flux and error associated with the best-fit model for one 
 *  of the comparison sources (1634 or 1635).
 * NB: This should only be used in programs after the program
 *  best_cmod has been run.
 *
 * Function called by:
 *   main in getflux.c
 *
 * Inputs: char *filename      input filename
 *         float *flux         best-fit flux (set by this function)
 *         float *err          rms error on best-fit model (set by function)
 *         float *chisq        chisq value for best-fit model (set by function)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int get_best_cal(char *filename, float *flux, float *err, float *chisq)
{
  int i;             /* Looping variable */
  int nlines;        /* Number of lines in input file */
  float bestchi;     /* Lowest value of chisq */
  Modfit *modfit;    /* Model fitting info */
  Modfit *mptr;      /* Pointer used for navigating modfit array */

  /*
   * Get model-fit info
   */

  if(!(modfit = read_modfit(filename,&nlines))) {
    fprintf(stderr,"ERROR: get_best_cal\n");
    return 0;
  }

  /*
   * Initialize
   */

  *flux = modfit->flux;
  *err = modfit->rms;
  bestchi = modfit->chisq;

  /*
   * Get best-fit flux and error
   */

  for(i=0,mptr=modfit; i<nlines; i++,mptr++)
    if(mptr->chisq < bestchi) {
      *flux = mptr->flux;
      *err = mptr->rms;
      bestchi = mptr->chisq;
    }

  *err *= 1000.0;
  *chisq = bestchi;

  /*
   * Clean up
   */

  modfit = del_modfit(modfit);

  return 0;
}

/*.......................................................................
 *
 * Function read_modfit
 *
 * Opens file with modfit information and reads the data into a Modfit
 *  array.
 *
 * Inputs: char *filename      input filename
 *         int *nlines         number of lines in input file. Set by this
 *                              function.
 *
 * Output: Modfit *modfit      filled array, NULL on error
 *
 */

Modfit *read_modfit(char *filename, int *nlines)
{
  int no_error=1;         /* Flag set to 0 on error */
  Modfit *modfit=NULL;    /* Model fitting info array */
  FILE *ifp=NULL;         /* Pointer to modfit info file */

  /*
   * Check for existence of input file and open it if it exists
   */

  if(!(ifp = open_readfile(filename))) {
    fprintf(stderr,"ERROR: read_modfit\n");
    return NULL;
  }

  /*
   * Find out the size of the input file
   */

  if((*nlines = n_lines(ifp,'#')) <= 0) {
    fprintf(stderr,"ERROR: read_modfit. No data in %s\n",filename);
    no_error = 0;
  }
  else
    printf("\nThere are %d lines in %s\n",*nlines,filename);

  /*
   * Allocate memory for input data
   */

  if(no_error)
    if(!(modfit = new_modfit(*nlines)))
      no_error = 0;

  /*
   * Read data into array
   */

  if(no_error)
    if(fill_modfit(modfit,*nlines,ifp) == 1)
      no_error = 0;

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

  if(no_error)
    return modfit;
  else {
    fprintf(stderr,"ERROR: read_modfit\n");
    return del_modfit(modfit);
  }
}

/*.......................................................................
 *
 * Function fill_modfit
 *
 * Reads input file and fills Modfit array with model fitting info
 *
 * Inputs: Modfit *modfit      model-fitting info array
 *         int nlines          number of lines in input file
 *         FILE *ifp           input file pointer
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int fill_modfit(Modfit *modfit, int nlines, FILE *ifp)
{
  int count=0;     /* Keeps track of number of input lines read */
  char line[MAXC]; /* General char string for reading input */
  Modfit *mptr;    /* Pointer to navigate modfit array */

  /*
   * Initialize
   */

  rewind(ifp);
  mptr = modfit;

  /*
   * Read in data
   */

  while(fgets(line,MAXC,ifp) != NULL) {
    if(line[0] != '#') {
      if(count >= nlines) {
	fprintf(stderr,"ERROR: fill_modfit.  Data reading error\n");
	return 1;
      }
      else if(sscanf(line,"%d %f %f %f %f",&mptr->model,&mptr->scale,
		  &mptr->flux,&mptr->chisq,&mptr->rms) != 5) {
	fprintf(stderr,"ERROR: fill_modfit.  Data not in correct format\n");
	return 1;
      }
      count++;
      mptr++;
    }
  }

  if(count != nlines) {
    fprintf(stderr,"ERROR: fill_modfit.\n");
    fprintf(stderr,"  Number of lines read does not match");
    fprintf(stderr," size of array (%d read, %d expected)\n",count,nlines);
    return 1;
  }

  return 0;
}

/*.......................................................................
 *
 * Function best_model
 *
 * Finds the best fitting model based on the model-fitting chisq by fitting
 *  a parabola to the three points with the lowest chisq (function
 *  fit_chi)
 *
 * Inputs: Modfit *modfit      model-fitting info array
 *         int nlines          number of lines in the array
 *
 * Output: float bestscal      scale factor of best fitting model,
 *                              neg. on error
 *
 * v21May01 CDF, Fixed error in finding 3 lowest points in modfit.
 */

float best_model(Modfit *modfit, int nlines)
{
  int i;                      /* Looping variable */
  int lowesti,loweri,lowi;    /* Indices associated with three lowest chisq */
  float slow,slower,slowest;  /* Three lowest values of the scaling */
  float clow,clower,clowest;  /* Three lowest values of chisq */
  float bestscale;            /* Scale factor of best fit model */
  Modfit *mptr;               /* Pointer used to navigate modfit array */

  if(nlines < 3) {
    fprintf(stderr,"ERROR: best_model.  ");
    fprintf(stderr,"There are fewer than 3 chisq values\n");
    return -999.9;
  }

  /*
   * Initialize
   */

  lowesti = 0;

  /*
   * Loop through the array to find the index corresponding
   *  to the lowest value of chisq.
   */

  for(i=0,mptr=modfit; i<nlines; i++,mptr++)
    if(mptr->chisq < (modfit+lowesti)->chisq)
      lowesti = i;

  /*
   * Now find the second-lowest value of chisq.
   */

  i = 0;
  while(i == lowesti)
    i++;
  loweri = i;

  for(i=0,mptr=modfit; i<nlines; i++,mptr++)
    if(mptr->chisq < (modfit+loweri)->chisq && i != lowesti)
      loweri = i;
  /*
   * Find the third-lowest value of chisq.
   */

  i = 0;
  while(i == lowesti || i == loweri)
    i++;
  lowi = i;

  for(i=0,mptr=modfit; i<nlines; i++,mptr++)
    if(mptr->chisq < (modfit+lowi)->chisq && i != lowesti && i != loweri)
      lowi = i;

  slow = modfit[lowi].scale;
  clow = modfit[lowi].chisq;
  slower = modfit[loweri].scale;
  clower = modfit[loweri].chisq;
  slowest = modfit[lowesti].scale;
  clowest = modfit[lowesti].chisq;

  if((bestscale = fit_parab(slow,clow,slower,clower,slowest,clowest)) < 0.0) {
    fprintf(stderr,"ERROR: best_model\n");
    return -999.0;
  }

  return bestscale;
}

/*.......................................................................
 *
 * Function source_flux
 *
 * Gets the component fluxes from a difmap *.mod file.  This now works for
 *  general lens fluxes as well as secondary calibrator source fluxes.
 *
 * Function called by:
 *   main in  getflux.c
 *
 * Inputs: char *date          observing date
 *         char *root          source rootname (e.g. 1608)
 *         char *modtype       model type (phasecal or gscale)
 *         int *ncomp          number of components -- set by this function
 *         float *rms          rms noise -- set by this function
 *         float *chisq        uv-plane chisq of fit -- set by this function
 *         float *obsdate      JD of observation -- set by this function
 *
 * Output: Secat *newdata      component information
 *
 * v22Mar00 CDF, Made more flexible by passing the source root name to
 *                the function, rather than hardwiring it to 1608.  This
 *                allows this function to be used for the secondary flux
 *                calibrators that can be fit with difmap model components.
 * v14Jul03 CDF, Made still more flexible in being able to handle lenses
 *                other than 1608.
 * v07Oct05 CDF, Made still more general by taking out completely
 *                any hardwired names.  Now returns the allocated
 *                flux array.
 *               Moved the file reading from this function into the new
 *                read_difmap function in dataio.c
 */

Secat *source_flux(char *date, char *root, char *modtype, int *ncomp,
		   float *rms, float *chisq, double *obsdate)
{
  int i;              /* Looping variable */
  int no_error=1;     /* Flag set to 0 on error */
  float *fptr;        /* Pointer to navigate flux array */
  char line[MAXC];    /* General char string for reading input */
  char modname[MAXC]; /* Name of model file */
  char keywd[MAXC];   /* Keyword value after a comment character */
  char junk[MAXC];    /* Junk string on comment line */
  Secat *comps=NULL;  /* Model component information */
  Secat *cptr;        /* Pointer for navigating comps */
  Skypos modcent;     /* RA, Dec of field center */

  /*
   * Get name of model file
   */

  sprintf(modname,"%s_%s_%s.mod",date,root,modtype);

  /*
   * Read in the components
   */

  if(!(comps = read_difmap(modname,'!',ncomp,&modcent))) {
    fprintf(stderr,"ERROR: source_flux\n");
    return 1;
  }
  printf("source_flux: Obtained data on %d components from %s\n",
	 ncomp,modname);

  /*
   * Initialize
   */

  *rms = 0.0;

  /*
   * Reopen model file for reading additional info
   */

  if(!(mfp = open_readfile(modname))) {
    fprintf(stderr,"ERROR: source_flux\n");
    return 1;
  }

  /*
   * Get additional info from model file
   */

  while((fgets(line,MAXC,mfp) != NULL) && no_error) {
    if(line[0] == '!') {
      if(sscanf(line,"! %s",keywd) == 1) {
	if(strcmp(keywd,"RMS") == 0) {
	  if(sscanf(line,"! RMS %s %f",junk,rms) != 2) {
	    fprintf(stderr,"ERROR: source_flux. Data file has incorrect ");
	    fprintf(stderr,"format\n");
	    no_error = 0;
	  }
	}
	else if(strcmp(keywd,"chisq") == 0) {
	  if(sscanf(line,"! chisq %s %f",junk,chisq) != 2) {
	    fprintf(stderr,"ERROR: source_flux. Data file has incorrect ");
	    fprintf(stderr,"format\n");
	    no_error = 0;
	  }
	}
	else if(strcmp(keywd,"Epoch") == 0) {
	  if(sscanf(line,"! Epoch %lf",obsdate) != 1) {
	    fprintf(stderr,"ERROR: source_flux. Data file has incorrect ");
	    fprintf(stderr,"format\n");
	    no_error = 0;
	  }
	}
      }
    }
  }

#if 0
  /*
   * If we have used CLEAN rather than modelfitting, open up new
   *  model file and get A, B, C, and D flux densities from that.
   */

  if(!DOMODEL && no_error) {
    sprintf(modname,"%s_abcd_%s.mod",root,modtype);
    if(get_clean_flux(flux,modname))
      no_error = 0;
  }
#endif
  /*
   * Convert fluxes to mJy
   */

  *rms *= 1000.0;

  /*
   * Clean up
   */

  if(mfp)
    fclose(mfp);

  if(no_error)
    return 0;
  else
    return 1;
}

/*.......................................................................
 *
 * Function get_clean_flux
 *
 * Takes a model file containing clean components that have come from
 *  clean windows surrounding the four components of 1608.  Reads in the
 *  window file to determine window boundaries, then sums up the flux
 *  densities of all components in each window.
 *
 * Inputs: float *flux         flux array to be filled by this function
 *         char *modname       name of file containing clean components
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int get_clean_flux(float *flux, char *modname)
{
  int i=0;             /* Looping variable */
  int no_error=1;      /* Flag set to 0 on error */
  int nwin;            /* Number of windows */
  float r,theta;       /* Component positions in r and theta */
  float compflux;      /* Component flux density */
  char winfile[MAXC];  /* File containing window information */
  char line[MAXC];     /* General string for reading input */
  Pos winmin[4];       /* Bottom left corners of windows */
  Pos winmax[4];       /* Top right corners of windows */
  Pos comppos;         /* Component position in x and y */
  FILE *mfp=NULL;      /* Model file pointer */
  FILE *wfp=NULL;      /* Window file pointer */

  /*
   * Open window file
   */

  sprintf(winfile,"%s1608_abcd.win",MODDIR);
  if(!(wfp = open_readfile(winfile))) {
    fprintf(stderr,"ERROR: get_clean_flux.\n");
    return 1;
  }

  /*
   * Read in windows -- there should be four of them
   */

  if((nwin = n_lines(wfp,'!')) != 4) {
    fprintf(stderr,"ERROR: get_clean_flux.  Expected 4 windows, found %d\n",
	    nwin);
    return 1;
  }
  else {
    rewind(wfp);
    while(fgets(line,MAXC,wfp) != NULL && no_error) {
      if(line[0] != '!') {
	if(sscanf(line,"%lf %lf %lf %lf",&winmin[i].x,&winmax[i].x,
		  &winmin[i].y,&winmax[i].y) != 4) {
	  fprintf(stderr,"ERROR: get_clean_flux.  Bad window file format.\n");
	  no_error = 0;
	}
	else {
	  i++;
	}
      }
    }
  }

  /*
   * Now open up model file
   */

  if(!(mfp = open_readfile(modname)))
    no_error = 0;

  /*
   * Step through model file and find fluxes that are inside the windows
   */

  while(fgets(line,MAXC,mfp) != NULL && no_error) {
    if(line[0] != '!') {
      if(sscanf(line,"%f %f %f",&compflux,&r,&theta) != 3) {
	fprintf(stderr,"ERROR: get_clean_flux.  Bad window file format.\n");
	no_error = 0;
      }
      else {

	/*
	 * Convert to r,theta to x,y
	 */

	rth2xy(r,theta,&comppos,1);

	/*
	 * See if this lies in any of the boxes
	 */

	for(i=0; i<4; i++) {
	  if(comppos.x > winmin[i].x && comppos.x < winmax[i].x &&
	     comppos.y > winmin[i].y && comppos.y < winmax[i].y) {
	    flux[i] += compflux;
	    break;
	  }
	}
      }
    }
  }

  /*
   * Clean up and exit
   */

  if(mfp)
    fclose(mfp);
  if(wfp)
    fclose(wfp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: get_clean_flux\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function print_lensflux
 *
 * Prints out flux densities a the lens system.
 * Called by getflux.c.
 *
 * Inputs: char *root          source root name (e.g., 1608)
 *         float *fluxp        phasecal flux densities
 *         float *fluxg        gscale flux densities
 *         int ncomp           number of components
 *         float rmsp          rms noise from phasecal map
 *         float rmsg          rms noise from gscale map
 *         float chisqp        chisq of fit to data - phasecal
 *         float chisqg        chisq of fit to data - gscale
 *         double obsdate      MJD - 50,000 of observing date
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 * v20Apr01 CDF, Changed output format to have day, the four main flux
 *                densities, followed by four columns of the RMS noise
 *                from the map, followed by the chisq of the fit.  The
 *                reason for repeating the RMS four times is to make the
 *                output format:
 *                  1) Compatible with the "mc" format written out by the
 *                      fluxplot.sm macros.
 *                  2) Future-compatible with runs of the automapping
 *                      procedure which will produce estimates of the errors
 *                      on each of the fitted values of the flux density.
 *               Added header information for the output files.
 * v15Jul03 CDF, Made more general so that can print out results for
 *                two-image lenses, too (actually lenses with an arbitrary
 *                number of images).
 */

int print_lensflux(char *root, float *fluxp, float *fluxg, int ncomp, 
		   float rmsp, float rmsg, float chisqp, float chisqg, 
		   double obsdate)
{
  int i;               /* Looping variable */
  int no_error=1;      /* Flag set to 0 on error */  
  int add_phead=0;     /* Flag set to 1 to print out header info */
  int add_ghead=0;     /* Flag set to 1 to print out header info */
  int add_dhead=0;     /* Flag set to 1 to print out header info */
  float *pptr,*gptr;   /* Pointers for navigating fluxp and fluxg */
  char pname[MAXC];    /* Output filename for phasecal fluxes */
  char gname[MAXC];    /* Output filename for gscale fluxes */
  char dname[MAXC];    /* Output filename for gscale - phasecal fluxes */
  FILE *pfp=NULL;      /* Output file pointer for phasecal fluxes */
  FILE *gfp=NULL;      /* Output file pointer for gscale fluxes */
  FILE *dfp=NULL;      /* Output file pointer for difference fluxes */

  /*
   * Set the output filenames
   */

  sprintf(pname,"%s%s_p.edt",REDDIR,root);
  sprintf(gname,"%s%s_g.edt",REDDIR,root);
  sprintf(dname,"%s%s_diff.edt",REDDIR,root);

  /*
   * Check for the existence of the output file.  If it doesn't 
   *  exist, then set flag for printing out header information.
   */

  if(file_exists(pname) == 0)
    add_phead = 1;
  if(file_exists(gname) == 0)
    add_ghead = 1;
  if(file_exists(dname) == 0)
    add_dhead = 1;

  /*
   * Open the output files for appending
   */

  if(!(pfp = open_appendfile(pname,1)))
    no_error = 0;
  if(!(gfp = open_appendfile(gname,1)))
    no_error = 0;
  if(!(dfp = open_appendfile(dname,1)))
    no_error = 0;

  /*
   * Print headers if required.
   */

  if(add_phead && no_error) {
    fprintf(pfp,"# Day    ");
    for(i=1; i<=ncomp; i++)
      fprintf(pfp,"  S_%d   ",i);
    for(i=1; i<=ncomp; i++)
      fprintf(pfp," rms   ");
    fprintf(pfp," chisq\n");
    fprintf(pfp,"#------  ");
    for(i=1; i<=ncomp; i++)
      fprintf(pfp,"------- ");
    for(i=1; i<=ncomp; i++)
      fprintf(pfp,"------ ");
    fprintf(pfp,"--------\n");
  }
  if(add_ghead && no_error) {
    fprintf(gfp,"# Day    ");
    for(i=1; i<=ncomp; i++)
      fprintf(gfp,"  S_%d   ",i);
    for(i=1; i<=ncomp; i++)
      fprintf(gfp," rms   ");
    fprintf(gfp," chisq\n");
    fprintf(gfp,"#------  ");
    for(i=1; i<=ncomp; i++)
      fprintf(gfp,"------- ");
    for(i=1; i<=ncomp; i++)
      fprintf(gfp,"------ ");
    fprintf(gfp,"--------\n");
  }
  if(add_dhead && no_error) {
    fprintf(dfp,"# Day    ");
    for(i=1; i<=ncomp; i++)
      fprintf(dfp,"(g-p)_%d ",i);
    fprintf(dfp,"\n");
    fprintf(dfp,"#------- ");
    for(i=1; i<=ncomp; i++)
      fprintf(dfp,"------- ");
    fprintf(dfp,"\n");
  }

  /*
   * Print data to output files
   */

  if(no_error) {

    /*
     * Print observing date
     */
    fprintf(pfp,"%7.2f  ",obsdate);
    fprintf(gfp,"%7.2f  ",obsdate);
    fprintf(dfp,"%7.2f  ",obsdate);

    /*
     * Print flux densities
     */

    for(i=0,pptr=fluxp,gptr=fluxg; i<ncomp; i++,pptr++,gptr++) {
      fprintf(pfp,"%7.3f ",*pptr);
      fprintf(gfp,"%7.3f ",*gptr);
      fprintf(dfp,"%7.3f ",(*gptr - *pptr));
    }

    /*
     * Print rms and chisq values
     */

    for(i=0; i<ncomp; i++) {
      fprintf(pfp,"%6.3f ",rmsp);
      fprintf(gfp,"%6.3f ",rmsg);
    }

    fprintf(pfp,"%8.4f\n",chisqp);
    fprintf(gfp,"%8.4f\n",chisqg);
    fprintf(dfp,"\n");
  }

  /*
   * Clean up and exit
   */

  if(pfp)
    fclose(pfp);
  if(gfp)
    fclose(gfp);
  if(dfp)
    fclose(dfp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: print_lensflux\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function print_calflux
 *
 * Prints out flux densities for a secondary flux density calibrator.
 * Called by getflux.c.
 *
 * Inputs: char *root          source root name (e.g., 1400)
 *         float *fluxp        phasecal flux densities
 *         float *fluxg        gscale flux densities
 *         int ncomp           number of components
 *         float rmsp          rms noise from phasecal map
 *         float rmsg          rms noise from gscale map
 *         float chisqp        chisq of fit to data - phasecal
 *         float chisqg        chisq of fit to data - gscale
 *         double obsdate      MJD - 50,000 of observing date
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int print_calflux(char *root, float *fluxp, float *fluxg, int ncomp, 
		  float rmsp, float rmsg, float chisqp, float chisqg, 
		  double obsdate)
{
  int i;               /* Looping variable */
  int no_error=1;      /* Flag set to 0 on error */  
  int add_header=0;    /* Flag set to 1 to print out header info */
  float *pptr,*gptr;   /* Pointers for navigating fluxp and fluxg */
  char outname[MAXC];  /* Output filename for phasecal fluxes */
  FILE *ofp=NULL;      /* Output file pointer for phasecal fluxes */

  /*
   * Set the output filename
   */

  sprintf(outname,"%s%s.edt",REDDIR,root);

  /*
   * Check for the existence of the output file.  If it doesn't 
   *  exist, then set flag for printing out header information.
   */

  if(file_exists(outname) == 0)
    add_header = 1;

  /*
   * Open the output files for appending
   */

  printf("\n");
  if(!(ofp = open_appendfile(outname,1)))
    no_error = 0;

  /*
   * Print data to output files
   */

  if(no_error) {

    /*
     * Print header if required.
     */

    if(add_header) {
      fprintf(ofp,"# Day      S_p      S_g     Diff     rms_p  chisq_p   ");
      fprintf(ofp," rms_g  chisq_g\n");
      fprintf(ofp,"#------  -------- -------- -------  ------- --------  ");
      fprintf(ofp,"------- --------\n");
    }

    /*
     * Print the observing date
     */

    fprintf(ofp,"%7.2f  ",obsdate);

    /*
     * Print flux densities
     */

    for(i=0,pptr=fluxp,gptr=fluxg; i<ncomp; i++,pptr++,gptr++) {
      fprintf(ofp,"%8.3f ",*pptr);
      fprintf(ofp,"%8.3f ",*gptr);
      fprintf(ofp,"%7.3f ",(*gptr - *pptr));
    }

    /*
     * Print rms and chisq values
     */

    fprintf(ofp," %7.4f %8.4f ",rmsp,chisqp);
    fprintf(ofp," %7.4f %8.4f\n",rmsg,chisqg);
  }

  /*
   * Clean up and exit
   */

  if(ofp)
    fclose(ofp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: print_lensflux\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function fit_parab
 *
 * Fits a parabola to 3 points (x1,y1), (x2,y2) and (x3,y3).  This can
 *  be done analytically because the equation for a parabola is parameterized
 *  by 3 values (a, b and c).  Thus you have 3 equations in 3 unknowns:
 *
 *     a x1^2 + b x1 + c = 0
 *     a x2^2 + b x2 + c = 0
 *     a x3^2 + b x3 + c = 0.
 *
 *  which gives:
 *
 *     a = ((y1-y3)(x2-x3) - (y2-y3)(x1-x3))/(x1-x2)(x2-x3)(x1-x3)
 *     b = (y2-y3)/(x2-x3) - a(x2+x3)
 *     c = y3 - a*x3^2 - b*x3
 *
 * Inputs: float x1
 *         float y1
 *         float x2
 *         float y2
 *         float x3
 *         float y3
 *
 * Output: float min           value of x at the minimum y
 *
 */

float fit_parab(float x1, float y1, float x2, float y2, float x3, float y3)
{
  float a,b,c;     /* Parabola parameters */

  a = ((y1-y3)*(x2-x3) - (y2-y3)*(x1-x3))/((x1-x2)*(x2-x3)*(x1-x3));
  b = (y2-y3)/(x2-x3) - a*(x2+x3);
  c = y3 - a*x3*x3 - b*x3;

  printf("\nfit_parab: For the 3 points:\n");
  printf("fit_parab:   (%f,%f)\n",x1,y1);
  printf("fit_parab:   (%f,%f)\n",x2,y2);
  printf("fit_parab:   (%f,%f)\n",x3,y3);
  printf("fit_parab: the parabola has parameters a=%g, b=%f, c=%f\n",a,b,c);
  printf("fit_parab: The minimum value occurs at x = %f\n\n",-b/(2*a));

  return -b/(2*a);
}
