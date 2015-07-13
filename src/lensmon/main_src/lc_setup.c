/*
 * lc_setup.c
 *
 * A library of functions to perform operations on the Setup structure
 *  used by the light-curve processing programs.
 *
 * 13Jun99 CDF,  Split out from lc_funcs.c
 * v16Sep99 CDF, Made separate interactive setup functions and summary
 *                functions for the separate programs (interpolation
 *                and delay-finding) by splitting up setup_interactive
 *                and setup_summary.
 *               Also split off all of the interactive functions in
 *                setup_interp to get_* to make things easier to read.
 * v12Sep00 CDF, Got rid of the SMTHENINT and SMIPINT options since these
 *                can now be emulated by running interp.c twice -- the
 *                first time w/SMONLY or SMINPLACE and the second time
 *                with INTONLY.
 * v13Sep00 CDF, Added the set_grid_params function, which sets the 
 *                parameters for the grid of points which have been 
 *                interpolated onto a regularly sampled grid.
 * v19Feb02 CDF, Changed setup_monte to also request time delay values
 *                if they haven't been set by setup_file.
 *               Changed setup_monte_summary to print out values contained
 *                in setup->tau0 rather than hardwired values (ALAG, etc.).
 * v04Sep02 CDF, Added a "tauset" variable to setup to indicate whether
 *                the tau0 values have been set in the input file and a
 *                "dtau" variable to set the stepsize used in the tau
 *                grid search.
 * v18Aug2005 CDF, Added a dmu variable to setup in order to control
 *                  the stepsize in the magnification grid.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "lc_funcs.h"
#include "lc_setup.h"


/*.......................................................................
 *
 * Function new_setup
 *
 * Allocates memory for a Setup structure array.
 *
 * Input:  int size            size of array
 *
 * Output: Setup *newsetup     new setup array
 *
 */

Setup *new_setup(int size)
{
  Setup *newsetup;

  newsetup = (Setup *) malloc((sizeof(Setup) * size));
  if(!newsetup) {
    fprintf(stderr,"\nnew_setup: Insufficient memory for Setup array.\n");
    return NULL;
  }

  /*
   * Initialize the setup parameters
   */

  newsetup->dochi = UNSET;
  newsetup->doxcorr = UNSET;
  newsetup->doacorr = NO;
  newsetup->dodisp = UNSET;
  newsetup->dodcf = NO;
  newsetup->docurvefit = UNSET;
  newsetup->dispchoice = UNSET;
  newsetup->d2delta = -1.0;
  newsetup->dosmooth = SMUNSET;
  newsetup->smtype = -1;
  newsetup->smwidth = 0.0;
  newsetup->ninterp = 0;
  newsetup->intstep = 0.0;
  newsetup->intstart = -1.0;
  newsetup->askstart = YES;
  newsetup->nvar = 0;
  newsetup->flagbad = UNSET;
  newsetup->meanchoice = UNSET;
  newsetup->mu0[0] = 0.0;
  newsetup->mu0[1] = 0.0;
  newsetup->mu0[2] = 0.0;
  newsetup->mu0[3] = 0.0;
  newsetup->dmu = 0.0005;
  newsetup->nmu = 0;
  newsetup->tauset = UNSET;
  newsetup->tau0[0] = 0.0;
  newsetup->tau0[1] = 0.0;
  newsetup->tau0[2] = 0.0;
  newsetup->tau0[3] = 0.0;
  newsetup->dtau = 0.0;
  newsetup->ntau = 0;
  newsetup->outfile = NULL;
  sprintf(newsetup->achifile,"chiba.dat");
  sprintf(newsetup->cchifile,"chibc.dat");
  sprintf(newsetup->dchifile,"chibd.dat");
  sprintf(newsetup->chilog,"stdout");
  sprintf(newsetup->xclog,"stdout");
  sprintf(newsetup->root,"mc_g");

  return newsetup;
}

/*.......................................................................
 *
 * Function del_setup
 *
 * Frees memory associated with Setup array
 *
 * Input:  Setup *setup        array to be freed
 *
 * Output: NULL
 *
 */

Setup *del_setup(Setup *setup)
{
  if(setup)
    free(setup);

  return NULL;
}

/*.......................................................................
 *
 * Function setup_file
 *
 * Fills in setup container with information from a file.
 *
 * Inputs: char *inname        name of file containing info
 *         Setup *setup        setup structure to be filled
 *
 * Output:int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int setup_file(Setup *setup, char *inname)
{
  int no_error=1;     /* Flag set to 0 on error */
  char keyword[MAXC]; /* Keyword at the beginning of each input line */
  char line[MAXC];    /* General string for reading input */
  char tmp[MAXC];     /* Temporary string variable */
  FILE *ifp=NULL;     /* Input file pointer */

  /*
   * Open input file
   */

  while(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: setup_file.\n");
    return 1;
  }

  /*
   * Read input from setup file.
   */

  printf("\nReading setup info from file %s\n\n",inname);
  while(fgets(line,MAXC,ifp) != NULL && no_error) {
    if(line[0] >= 32 && line[0] != '#') {
      switch(read_setup_line(line,keyword)) {
      case SETUPERR:
	no_error = 0;
	break;
      case DOCHI:
	if(sscanf(line,"%s %d",keyword,&setup->dochi) != 2 &&
	   setup->dochi < 0) {
	  fprintf(stderr,"ERROR: setup_file.  Bad value for dochi\n");
	  fprintf(stderr," Setting dochi = NO (0)\n");
	  setup->dochi = NO;
	}
	break;
      case DOXCORR:
	if(sscanf(line,"%s %d",keyword,&setup->doxcorr) != 2 &&
	   setup->doxcorr < 0) {
	  fprintf(stderr,"ERROR: setup_file.  Bad value for doxcorr\n");
	  fprintf(stderr," Setting doxcorr = NO\n");
	  setup->doxcorr = NO;
	}
	break;
      case DOACORR:
	if(sscanf(line,"%s %d",keyword,&setup->doacorr) != 2 &&
	   setup->doacorr < 0) {
	  fprintf(stderr,"ERROR: setup_file.  Bad value for doacorr\n");
	  fprintf(stderr," Setting doacorr = NO (0)\n");
	  setup->doacorr = NO;
	}
	break;
      case DODISP:
	if(sscanf(line,"%s %d",keyword,&setup->dodisp) != 2 &&
	   setup->dodisp < 0) {
	  fprintf(stderr,"ERROR: setup_file.  Bad value for dodisp\n");
	  fprintf(stderr," Setting dodisp = NO (0)\n");
	  setup->dodisp = NO;
	}
	break;
      case DODCF:
	if(sscanf(line,"%s %d",keyword,&setup->dodcf) != 2 &&
	   setup->dodcf < 0) {
	  fprintf(stderr,"ERROR: setup_file.  Bad value for dodcf\n");
	  fprintf(stderr," Setting dodcf = NO (0)\n");
	  setup->dodcf = NO;
	}
	break;
      case DOCURVEFIT:
	if(sscanf(line,"%s %d",keyword,&setup->docurvefit) != 2 &&
	   setup->docurvefit < 0) {
	  fprintf(stderr,"ERROR: setup_file.  Bad value for docurvefit\n");
	  fprintf(stderr," Setting docurvefit = NO (0)\n");
	  setup->docurvefit = NO;
	}
	break;
      case DISPCHOICE:
	if(sscanf(line,"%s %d",keyword,&setup->dispchoice) != 2 &&
	   (setup->dispchoice < D21 || setup->dispchoice > DLOVELL)) {
	  fprintf(stderr,"ERROR: setup_file.  Bad value for dispchoice\n");
	  fprintf(stderr," Setting dispchoice = 1\n");
	  setup->dispchoice = 1;
	}
	break;
      case D2DELTA:
	if(sscanf(line,"%s %f",keyword,&setup->d2delta) != 2 ||
	   setup->d2delta < 0.0) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for d2delta\n");
	  fprintf(stderr," Setting d2delta = 5.0\n");
	  setup->d2delta = 5.0;
	}
	break;
      case DOOVERLAP:
	if(sscanf(line,"%s %d",keyword,&setup->dooverlap) != 2 ||
	   setup->dooverlap < 0) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for dooverlap\n");
	  fprintf(stderr," Setting dooverlap = NO (0)\n");
	  setup->dooverlap = NO;
	}
	break;
      case OUTFILE:
	if(sscanf(line,"%s %s",keyword,&tmp) != 2) {
	  fprintf(stderr,
		  "ERROR: setup_file.  Bad value for outfile.\n");
	  fprintf(stderr,"NO output file will be set\n");
	  setup->outfile = NULL;
	}
	else {
	  if(!(setup->outfile = new_string(MAXC)))
	    setup->outfile = NULL;
	  else
	    strcpy(setup->outfile,tmp);
	}
	break;
      case ACHIFILE:
	if(sscanf(line,"%s %s",keyword,setup->achifile) != 2) {
	  fprintf(stderr,
		  "ERROR: setup_file.  Bad value for achifile.\n");
	  fprintf(stderr,"Setting achifile to chiba.dat.\n");
	  sprintf(setup->achifile,"chiba.dat");
	}
	break;
      case CCHIFILE:
	if(sscanf(line,"%s %s",keyword,setup->cchifile) != 2) {
	  fprintf(stderr,
		  "ERROR: setup_file.  Bad value for cchifile.\n");
	  fprintf(stderr,"Setting cchifile to chibc.dat.\n");
	  sprintf(setup->cchifile,"chibc.dat");
	}
	break;
      case DCHIFILE:
	if(sscanf(line,"%s %s",keyword,setup->dchifile) != 2) {
	  fprintf(stderr,
		  "ERROR: setup_file.  Bad value for dchifile.\n");
	  fprintf(stderr,"Setting dchifile to chibd.dat.\n");
	  sprintf(setup->dchifile,"chibd.dat");
	}
	break;
      case MU0:
	if(sscanf(line,"%s %f %f %f %f",keyword,&setup->mu0[0],&setup->mu0[1],
		  &setup->mu0[2],&setup->mu0[3]) != 5) {
	  fprintf(stderr,"ERROR: setup_file.  Bad inputs for m0\n");
	  fprintf(stderr," Setting mu0 = {0 0 0 0}\n");
	  setup->mu0[0] = 0.0;
	  setup->mu0[1] = 0.0;
	  setup->mu0[2] = 0.0;
	  setup->mu0[3] = 0.0;
	}
	break;
      case NMU:
	if(sscanf(line,"%s %d",keyword,&setup->nmu) != 2 ||
	   setup->nmu < 0) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for nmu\n");
	  fprintf(stderr," Setting nmu = 0\n");
	  setup->nmu = 0;
	}
	break;
      case TAU0:
	if(sscanf(line,"%s %lf %lf %lf %lf",keyword,&setup->tau0[0],
		  &setup->tau0[1],&setup->tau0[2],&setup->tau0[3]) != 5) {
	  fprintf(stderr,"ERROR: setup_file.  Bad inputs for tau0\n");
	  fprintf(stderr," Setting tau0 = {0 0 0 0}\n");
	  setup->tau0[0] = 0.0;
	  setup->tau0[1] = 0.0;
	  setup->tau0[2] = 0.0;
	  setup->tau0[3] = 0.0;
	}
	else
	  setup->tauset = YES;
	break;
      case DTAU:
	if(sscanf(line,"%s %lf",keyword,&setup->dtau) != 2 ||
	   setup->dtau < 0) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for dtau\n");
	  fprintf(stderr," Setting dtau = 0.0\n");
	  setup->dtau = 0.0;
	}
	break;
      case NTAU:
	if(sscanf(line,"%s %d",keyword,&setup->ntau) != 2 ||
	   setup->ntau < 0) {
	  fprintf(stderr,"ERROR: setup_file.  Bad input for ntau\n");
	  fprintf(stderr," Setting ntau = 0\n");
	  setup->ntau = 0;
	}
	break;
      case DOSMOOTH:
	if(sscanf(line,"%s %d",keyword,&setup->dosmooth) != 2 &&
	   setup->dosmooth > SMINPLACE) {
	  fprintf(stderr,"ERROR: setup_file.  Bad value for dosmooth\n");
	  setup->dosmooth = SMUNSET;
	}
	break;
      case BOXCAR:
	if(sscanf(line,"%s %f",keyword,&setup->smwidth) != 2 &&
	   (setup->smwidth <= 0.0 || setup->smwidth > 100.0)) {
	  fprintf(stderr,
		  "ERROR: setup_file.  Bad value for smoothing width.\n");
	  fprintf(stderr,"Setting smoothing width = 21 days.\n");
	  setup->smwidth = -1.0;
	}
	setup->smtype = SMBOXCAR;
	break;
      case MEDIAN:
	if(sscanf(line,"%s %f",keyword,&setup->smwidth) != 2 &&
	   (setup->smwidth <= 0.0 || setup->smwidth > 100.0)) {
	  fprintf(stderr,
		  "ERROR: setup_file.  Bad value for smoothing width.\n");
	  fprintf(stderr,"Setting smoothing width = 21 days.\n");
	  setup->smwidth = -1.0;
	}
	setup->smtype = SMMEDIAN;
	break;
      case VARBOX:
	if(sscanf(line,"%s %d",keyword,&setup->nvar) != 2 &&
	   (setup->nvar <= 0 || setup->nvar > 30)) {
	  fprintf(stderr,
		  "ERROR: setup_file.  Bad value for smoothing width.\n");
	  fprintf(stderr,"Setting number of smooth points to 5.\n");
	  setup->nvar = -1;
	}
	setup->smtype = SMVARBOX;
	break;
      case VARTRI:
	if(sscanf(line,"%s %d",keyword,&setup->nvar) != 2 &&
	   (setup->nvar <= 0 || setup->nvar > 30)) {
	  fprintf(stderr,
		  "ERROR: setup_file.  Bad value for smoothing width.\n");
	  fprintf(stderr,"Setting number of smooth points to 5.\n");
	  setup->nvar = -1;
	}
	setup->smtype = SMVARTRI;
	break;
      case TRIANGLE:
	if(sscanf(line,"%s %f",keyword,&setup->smwidth) != 2 &&
	   (setup->smwidth <= 0.0 || setup->smwidth > 100.0)) {
	  fprintf(stderr,
		  "ERROR: setup_file.  Bad value for smoothing width.\n");
	  fprintf(stderr,"Setting smoothing width = 21 days.\n");
	  setup->smwidth = -1.0;
	}
	setup->smtype = SMTRIANGLE;
	break;
      case GAUSS:
	if(sscanf(line,"%s %f",keyword,&setup->smwidth) != 2 &&
	   (setup->smwidth <= 0.0 || setup->smwidth > 100.0)) {
	  fprintf(stderr,
		  "ERROR: setup_file.  Bad value for smoothing width.\n");
	  fprintf(stderr,"Setting gaussian 1/e width = 21 days.\n");
	  setup->smwidth = -1.0;
	}
	setup->smtype = SMGAUSS;
	break;
      case NINTERP:
	if(sscanf(line,"%s %d",keyword,&setup->ninterp) != 2 &&
	   setup->ninterp < 0) {
	  fprintf(stderr,"ERROR: setup_file.  Bad value for ninterp.\n");
	  fprintf(stderr,"Setting ninterp = 0.\n\n");
	}
	break;
      case INTSTEP:
	if(sscanf(line,"%s %f",keyword,&setup->intstep) != 2 &&
	   (setup->intstep <= 0.0 || setup->intstep > 30.0)) {
	  fprintf(stderr,
		  "ERROR: setup_file.  Bad value for intstep.\n");
	  fprintf(stderr,"Setting interpolation step size = -1 days.\n");
	  setup->intstep = -1.0;
	}
	break;
      case INTSTART:
	if(sscanf(line,"%s %f",keyword,&setup->intstart) != 2 &&
	   setup->intstep <= 0.0) {
	  fprintf(stderr,
		  "ERROR: setup_file.  Bad value for intstart.\n");
	  fprintf(stderr,"Setting start day = -1.0.\n");
	  setup->intstart = -1.0;
	}
	break;
      case ASKSTART:
	if(sscanf(line,"%s %d",keyword,&setup->askstart) != 2) {
	  fprintf(stderr,
		  "ERROR: setup_file.  Bad value for askstart.\n");
	  fprintf(stderr,"Will ask for starting day for interpolation.\n");
	  setup->askstart = YES;
	}
	break;
      case CHILOG:
	if(sscanf(line,"%s %s",keyword,setup->chilog) != 2) {
	  fprintf(stderr,
		  "ERROR: setup_file.  Bad value for chisq logfile.\n");
	  fprintf(stderr,"Sending output to stdout.\n");
	}
	break;
      case XCLOG:
	if(sscanf(line,"%s %s",keyword,setup->xclog) != 2) {
	  fprintf(stderr,
		  "ERROR: setup_file.  Bad value for cross-corr. logfile.\n");
	  fprintf(stderr,"Sending output to stdout.\n");
	}
	break;
      case FLAGBAD:
	if(sscanf(line,"%s %d",keyword,&setup->flagbad) != 2) {
	  fprintf(stderr,"ERROR: setup_file.  Bad format for flagbad\n");
	  setup->flagbad = -1;
	}
	break;
      case MEANCHOICE:
	if(sscanf(line,"%s %d",keyword,&setup->meanchoice) != 2 ||
	   setup->meanchoice < 0 || setup->meanchoice > 2) {
	  fprintf(stderr,"ERROR: setup_file.  Bad format for meanchoice\n");
	  setup->meanchoice = UNSET;
	}
	break;
      case ROOT:
	if(sscanf(line,"%s %s",keyword,setup->root) != 2) {
	  fprintf(stderr,
		  "ERROR: setup_file.  Bad value for root.\n");
	  fprintf(stderr,"Setting root to mc_g.\n");
	  sprintf(setup->root,"mc_g");
	}
	break;
      default:
	printf("***WARNING: Not yet taking file info for keyword %s.\n",
	       keyword);
      }
    }
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: setup_file\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function read_setup_line
 *
 * Reads one line of the setup file and checks the keyword on the line.
 * Returns an integer value identifying the keyword value.
 *
 * Input:  char *line          input line
 *         char *keyword       keyword read from input file (modified by this
 *                              function)
 *
 * Output: int keyval          keyword value (see enumeration in fitsim.h)
 *                             SETUPERR for error
 *
 */

int read_setup_line(char *line, char *keyword)
{
  /*
   * Check to make sure that there is a keyword on the line.
   */

  if(sscanf(line,"%s",keyword) != 1) {
    fprintf(stderr,"ERROR: read_setup_line\n");
    keyword = NULL;
    return SETUPERR;
  }

  /*
   * If there is a keyword, read it and return the appropriate value
   */

  if(strcmp(keyword,"dochi") == 0 || strcmp(keyword,"DOCHI") == 0)
    return DOCHI;
  if(strcmp(keyword,"doxcorr") == 0 || strcmp(keyword,"DOXCORR") == 0)
    return DOXCORR;
  if(strcmp(keyword,"doacorr") == 0 || strcmp(keyword,"DOACORR") == 0)
    return DOACORR;
  if(strcmp(keyword,"dodisp") == 0 || strcmp(keyword,"DODISP") == 0)
    return DODISP;
  if(strcmp(keyword,"dodcf") == 0  || strcmp(keyword,"DODCF") == 0)
    return DODCF;
  if(strcmp(keyword,"docurvefit") == 0 || strcmp(keyword,"DOCURVEFIT") == 0)
    return DOCURVEFIT;
  if(strcmp(keyword,"dispchoice") == 0 || strcmp(keyword,"DISPCHOICE") == 0)
    return DISPCHOICE;
  if(strcmp(keyword,"d2delta") == 0 || strcmp(keyword,"D2DELTA") == 0)
    return D2DELTA;
  if(strcmp(keyword,"outfile") == 0 || strcmp(keyword,"OUTFILE") == 0)
    return OUTFILE;
  if(strcmp(keyword,"achifile") == 0 || strcmp(keyword,"ACHIFILE") == 0)
    return ACHIFILE;
  if(strcmp(keyword,"cchifile") == 0 || strcmp(keyword,"CCHIFILE") == 0)
    return CCHIFILE;
  if(strcmp(keyword,"dchifile") == 0 || strcmp(keyword,"DCHIFILE") == 0)
    return DCHIFILE;
  if(strcmp(keyword,"boxcar") == 0 || strcmp(keyword,"BOXCAR") == 0)
    return BOXCAR;
  if(strcmp(keyword,"median") == 0 || strcmp(keyword,"MEDIAN") == 0)
    return MEDIAN;
  if(strcmp(keyword,"varbox") == 0 || strcmp(keyword,"VARBOX") == 0)
    return VARBOX;
  if(strcmp(keyword,"vartri") == 0 || strcmp(keyword,"VARTRI") == 0)
    return VARTRI;
  if(strcmp(keyword,"triangle") == 0 || strcmp(keyword,"TRIANGLE") == 0)
    return TRIANGLE;
  if(strcmp(keyword,"gauss") == 0 || strcmp(keyword,"GAUSS") == 0)
    return GAUSS;
  if(strcmp(keyword,"chilog") == 0 || strcmp(keyword,"CHILOG") == 0)
    return CHILOG;
  if(strcmp(keyword,"xclog") == 0 || strcmp(keyword,"XCLOG") == 0)
    return XCLOG;
  if(strcmp(keyword,"dooverlap") == 0 || strcmp(keyword,"DOOVERLAP") == 0)
    return DOOVERLAP;
  if(strcmp(keyword,"mu0") == 0 || strcmp(keyword,"MU0") == 0)
    return MU0;
  if(strcmp(keyword,"nmu") == 0 || strcmp(keyword,"NMU") == 0)
    return NMU;
  if(strcmp(keyword,"tau0") == 0 || strcmp(keyword,"TAU0") == 0)
    return TAU0;
  if(strcmp(keyword,"dtau") == 0 || strcmp(keyword,"DTAU") == 0)
    return DTAU;
  if(strcmp(keyword,"ntau") == 0 || strcmp(keyword,"NTAU") == 0)
    return NTAU;
  if(strcmp(keyword,"ninterp") == 0 || strcmp(keyword,"NINTERP") == 0)
    return NINTERP;
  if(strcmp(keyword,"intstep") == 0 || strcmp(keyword,"INTSTEP") == 0)
    return INTSTEP;
  if(strcmp(keyword,"intstart") == 0 || strcmp(keyword,"INTSTART") == 0)
    return INTSTART;
  if(strcmp(keyword,"askstart") == 0 || strcmp(keyword,"ASKSTART") == 0)
    return ASKSTART;
  if(strcmp(keyword,"flagbad") == 0 || strcmp(keyword,"FLAGBAD") == 0)
    return FLAGBAD;
  if(strcmp(keyword,"meanchoice") == 0 || strcmp(keyword,"MEANCHOICE") == 0)
    return MEANCHOICE;
  if(strcmp(keyword,"dosmooth") == 0 || strcmp(keyword,"DOSMOOTH") == 0)
    return DOSMOOTH;
  if(strcmp(keyword,"root") == 0 || strcmp(keyword,"ROOT") == 0)
    return ROOT;

  /*
   * If none of the above checks have been satisfied, return the
   *  default value
   */

  return DEFAULT;
}

/*.......................................................................
 *
 * Function setup_interp
 *
 * Fills in parts of the setup structure necessary for _interpolation_
 *  functions that weren't filled in from the setup file.
 *
 * Inputs: Setup *setup        setup structure
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int setup_interp(Setup *setup)
{
  char line[MAXC];      /* General string for reading input */

#if 0
  /*
   * Decide whether or not to flag bad days
   */

  if(setup->flagbad < 0) {
    setup->flagbad = YES;
    printf("\nFlag bad days? [y] ");
    fgets(line,MAXC,stdin);
    if(line[0] == 'n' || line[0] == 'N')
      setup->flagbad = NO;
  }

  /*
   * Choose method for normalizing secondary flux calibrators
   */

  if(setup->meanchoice == UNSET)
    get_meanchoice(setup);
#endif
  /*
   * Decide whether to do smoothing, if an interpolative analysis method
   *  has been chosen.
   * NB: smoothing is only used for chisq, xcorr, and acorr methods.
   */

  if(setup->dosmooth == SMUNSET)
    get_interp_choice(setup);

  /*
   * Get step size(s) for smoothed, interpolated curves, if not set.
   * NB: This is not required for the SMINPLACE choice.
   */

  if(setup->dosmooth != SMINPLACE  && setup->dosmooth > NOSMOOTH)
    get_interp_step(setup);

  /*
   * If no interpolation is desired or if setup->dosmooth is set to
   *  linear interpolation only, nothing more needs to be set, so return.
   */
    
  if(setup->dosmooth <= INTONLY)
    return 0;

  /*
   * Get smoothing function if it is not set.
   */

  if(setup->smtype < 0)
    get_smooth_fn(setup);
    
  /*
   * Get smoothing window width if not set
   */

  if(setup->smwidth <= 0.0 && setup->smtype < SMVARBOX)
    get_smooth_width(setup);

  /*
   * Get number of points for variable-width smoothing functions, if
   *  not set
   */

  if(setup->nvar <= 0 && setup->smtype >= SMVARBOX)
    get_nvar(setup);

  return 0;
}

/*.......................................................................
 *
 * Function setup_delays
 *
 * Fills in parts of the setup structure relating to the delay-finding
 *  algorithms that weren't filled in from setup file.
 *
 * Inputs: Setup *setup        setup structure
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int setup_delays(Setup *setup)
{
  int dointerpolate=NO; /* Flag set to YES if interpolating */
  char line[MAXC];      /* General string for reading input */

  /*
   * Pick analysis methods if not set by setup file
   */

  if(setup->dochi == UNSET) {
    setup->dochi = NO;
    printf("\nDo chisq analysis? [n] ");
    fgets(line,MAXC,stdin);
    if(line[0] == 'y' || line[0] == 'Y')
      setup->dochi = YES;
  }

  if(setup->doxcorr == UNSET) {
    setup->doxcorr = NO;
    printf("\nDo cross-correlation analysis? [n] ");
    fgets(line,MAXC,stdin);
    if(line[0] == 'y' || line[0] == 'Y')
      setup->doxcorr = YES;
  }

  if(setup->doacorr == UNSET) {
    setup->doacorr = NO;
    printf("\nDo auto-correlation analysis? [n] ");
    fgets(line,MAXC,stdin);
    if(line[0] == 'y' || line[0] == 'Y')
      setup->doacorr = YES;
  }

  if(setup->dodisp == UNSET) {
    setup->dodisp = NO;
    printf("\nDo dispersion analysis? [n] ");
    fgets(line,MAXC,stdin);
    if(line[0] == 'y' || line[0] == 'Y')
      setup->dodisp = YES;
  }

  if(setup->dodcf == UNSET) {
    setup->dodcf = NO;
    printf("\nDo discrete correlation analysis? [n] ");
    fgets(line,MAXC,stdin);
    if(line[0] == 'y' || line[0] == 'Y')
      setup->dodcf = YES;
  }

  if(setup->docurvefit == UNSET) {
    setup->docurvefit = NO;
    printf("\nDo simultaneous curve fitting and chisq minimization? [n] ");
    fgets(line,MAXC,stdin);
    if(line[0] == 'y' || line[0] == 'Y')
      setup->docurvefit = YES;
  }

  /*
   * If one of the interpolative methods was chosen, then set
   *  the dointerpolate flag to YES, otherwise set setup->dosmooth 
   *  to NOSMOOTH
   */

  if(setup->dochi || setup->doxcorr || setup->doacorr)
    dointerpolate = YES;
  else
    setup->dosmooth = NOSMOOTH;

  /*
   * If dispersion analysis is desired, choose method
   */

  if(setup->dodisp && setup->dispchoice == UNSET) {
    setup->dispchoice = D21;
    printf("\nChoose dispersion analysis method\n");
    printf("%d. Pelt et al. D^2_1 method\n",D21);
    printf("%d. Pelt et al. D^2_1 method (> 2 light curves)\n",D21M);
    printf("%d. Pelt et al. D^2_2 method\n",D22);
    printf("%d. Lovell et al. D^2_2 method\n",DLOVELL);
    printf("------------------------------------------------------\n");
    printf("Enter choice: [%d] ",setup->dispchoice);
    fgets(line,MAXC,stdin);
    if(line[0] >= 32) {
      while(sscanf(line,"%d",&setup->dispchoice) != 1 || 
	    setup->dispchoice < D21 || setup->dispchoice > DLOVELL) {
	fprintf(stderr,"ERROR: setup_delays.  Bad choice. ");
	fprintf(stderr,"Enter choice again:  ");
	fgets(line,MAXC,stdin);
      }
    }
  }

  /*
   * Set d2delta parameter if required
   */

  if(setup->d2delta < 0.0 && 
     (setup->dispchoice == D22 || setup->dispchoice == DLOVELL)) {
    setup->d2delta = 5.0;
    printf("\nEnter value of delta for D^2 dispersion analysis: [%3.1f] ",
	   setup->d2delta);
    fgets(line,MAXC,stdin);
    if(line[0] >= 32) {
      while(sscanf(line,"%f",&setup->d2delta) != 1 || 
	    setup->d2delta < 0.0) {
	fprintf(stderr,"ERROR: setup_delays.  Bad choice. ");
	fprintf(stderr,"Enter choice again:  ");
	fgets(line,MAXC,stdin);
      }
    }
  }

  /*
   * Get names of logfiles if not set
   */

  if(setup->dochi && strcmp(setup->chilog,"stdout") == 0) {
    printf("\nEnter name of chisq logfile: [%s] ",setup->chilog);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%s",setup->chilog) != 1) {
	fprintf(stderr,"ERROR: setup_delays.  Enter logfile name again: ");
	fgets(line,MAXC,stdin);
      }
    }
  }

  if(setup->doxcorr && strcmp(setup->xclog,"stdout") == 0) {
    printf("Enter name of cross-correlation logfile: [%s] ",setup->xclog);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%s",setup->xclog) != 1) {
	fprintf(stderr,"ERROR: setup_delays.  Enter logfile name again: ");
	fgets(line,MAXC,stdin);
      }
    }
  }

  return 0;
}

/*.......................................................................
 *
 * Function setup_interp_summary
 *
 * Summarizes the information in the setup container relating to the
 *  interpolation algorithms.
 *
 * Inputs: Setup *setup        setup container
 *
 * Outputs: (none)
 *
 */

void setup_interp_summary(Setup *setup)
{
  printf("\n------------------------------------------------------------\n");
  printf("\nSetup parameters:\n");
  printf("   Method for normalizing secondary flux cal. light curves:\n");
  if(setup->meanchoice == 0) 
    printf("      Division by total flux from modscal model.\n");
  else if(setup->meanchoice == 1) 
    printf("      Division by mean over total length of observations.\n");
  else
    printf("      Separate normalization for different array configs.\n");
  if(setup->flagbad <= 0)
    printf("   Bad day flagging:    no\n");
  else
    printf("   Bad day flagging:    yes\n");
  printf("\n------------------------------------------------------------\n");
  printf("\nInterpolation parameters:\n");
  if(setup->dosmooth > INTONLY) {
    switch(setup->smtype) {
    case SMBOXCAR:
      printf("   Smoothing function:  Boxcar\n");
      printf("   Smoothing width:     %-4.1f\n",setup->smwidth);
      break;
    case SMMEDIAN:
      printf("   Smoothing function:  Median\n");
      printf("   Smoothing width:     %-4.1f\n",setup->smwidth);
      break;
    case SMTRIANGLE:
      printf("   Smoothing function:  Triangle\n");
      printf("   Smoothing width:     %-4.1f\n",setup->smwidth);
      break;
    case SMGAUSS:
      printf("   Smoothing function:  Gaussian\n");
      printf("   Smoothing width:     %-4.1f\n",setup->smwidth);
      break;
    case SMVARBOX:
      printf("   Smoothing function:  Variable-width boxcar\n");
      printf("   Number of points:    %d\n",setup->nvar);
      break;
    case SMVARTRI:
      printf("   Smoothing function:  Variable-width triangle\n");
      printf("   Number of points:    %d\n",setup->nvar);
      break;
    default:
      fprintf(stderr,"ERROR: Invalid smoothing function.\n");
    }
  }

  switch(setup->dosmooth) {
  case NOSMOOTH:
    break;
  case SMINPLACE:
    printf("   Smooth in place.  No interpolation.\n");
  default:
    printf("   Interpolation step:   %-5.2f\n",setup->intstep);
  }

  printf("\n------------------------------------------------------------\n");
}

/*.......................................................................
 *
 * Function setup_delays_summary
 *
 * Summarizes the information in the setup container relating to the
 *  delay-finding algorithms.
 *
 * Inputs: Setup *setup        setup container
 *
 * Outputs: (none)
 *
 */

void setup_delays_summary(Setup *setup)
{
  printf("\n------------------------------------------------------------\n");
  if(setup->infile1)
    printf("Input file 1: %s\n",setup->infile1);
  if(setup->infile2)
    printf("Input file 1: %s\n",setup->infile2);
  printf("\nAnalysis method             Type (optional) Size (optional)\n");
  printf("--------------------------- --------------- ---------------\n");
  if(setup->dochi)
    printf(" Chisq minimization\n");
  if(setup->doxcorr)
    printf(" Cross-correlation analysis\n");
  if(setup->doacorr)
    printf(" Auto-correlation analysis\n");
  if(setup->dodisp) {
    printf(" Dispersion analysis:");
    switch(setup->dispchoice) {
    case D21:
      printf("       Pelt D^2_1\n");
      break;
    case D21M:
      printf("       Pelt D^2_1 (>2 curves)\n");
      break;
    case D22:
      printf("       Pelt D^2_2");
      printf("       delta: %-6.2f\n",setup->d2delta);
      break;
    case DLOVELL:
      printf("       Lovell D^2_2");
      printf("       delta: %-6.2f\n",setup->d2delta);
      break;
    default:
      printf("       Pelt D^2_1\n");
    }
  }
  if(setup->dodcf)
    printf(" Discrete correlation function\n");
  if(setup->docurvefit)
    printf(" Simultaneous curve fitting\n");
  if(setup->dochi)
    printf("\n Chisq Logfile:       %s\n",setup->chilog);
  if(setup->doxcorr)
    printf(" Cross-corr logfile:  %s\n",setup->xclog);
  printf("------------------------------------------------------------\n");
}

/*.......................................................................
 *
 * Function get_meanchoice
 *
 * Gets method for normalizing the secondary flux calibrator light curves
 *  using the following system:
 *
 *    0 ==> Normalize using the total flux from the modscal model
 *    1 ==> Normalize using the mean flux for the entire observation
 *    2 ==> Normalize using separate means for each array configuration
 *           (A, BnA, or B)
 *
 * Inputs: Setup *setup        container holding meanchoice variable
 *
 * Output: none
 *
 */

void get_meanchoice(Setup *setup)
{
  char line[MAXC];      /* General string for reading input */

  setup->meanchoice = 2;
  printf("\nMethod for normalizing secondary flux cal light curves:\n");
  printf("  0. Normalize by total flux from modscal model.\n");
  printf("  1. Normalize by mean flux over entire time of observations.\n");
  printf("  2. Normalize each array configuration separately.\n");
  printf("------------------------------------------------------------\n");
  printf("  Choice? [%d] ",setup->meanchoice);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%d",&setup->meanchoice) != 1 || setup->meanchoice < 0
	  || setup->meanchoice > 2) {
      fprintf(stderr,"ERROR: get_meanchoice.  Bad format for meanchoice.");
      fprintf(stderr,"  Enter new choice:  ");
      fgets(line,MAXC,stdin);
    }
  }
}

/*.......................................................................
 *
 * Function get_interp_choice
 *
 * Gets the method to be used for smoothing and/or interpolating the
 *  input light curve.
 *
 * Inputs: Setup *setup        container holding setup info
 *
 * Output: none
 *
 */

void get_interp_choice(Setup *setup)
{
  char line[MAXC];      /* General string for reading input */

  /*
   * Set default value
   */

  setup->dosmooth = SMONLY;

  /*
   * Get choice
   */

  printf("\nSmoothing/interpolation choice :\n");
  printf("  %d. Linear interpolation\n",INTONLY);
  printf("  %d. Smoothing and interpolate\n",SMONLY);
  printf("  %d. Smooth in place with NO interpolation\n",SMINPLACE);
  printf("Enter choice: [%d] ",setup->dosmooth);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%d",&setup->dosmooth) != 1 
	  || setup->dosmooth < INTONLY || setup->dosmooth > SMINPLACE) {
      fprintf(stderr,"ERROR: get_interp_choice.  Bad smoothing choice.\n");
      fprintf(stderr,"  Enter new choice:  ");
      fgets(line,MAXC,stdin);
    }
  }
}

/*.......................................................................
 *
 * Function get_interp_step
 *
 * Gets the step size to be used in smoothing and/or interpolating the
 *  input light curve.
 *
 * Inputs: Setup *setup        container holding setup info
 *
 * Output: none
 *
 */

void get_interp_step(Setup *setup)
{
  char line[MAXC];      /* General string for reading input */

  /*
   * Only get step size if it hasn't already been set (setup->intstep < 0)
   */

  if(setup->intstep <= 0.0) {

    /*
     * Get step size, after first setting a default value
     */

    setup->intstep = 1.0;
    printf("\nEnter step size for interpolation: [%4.2f] ",setup->intstep);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%f",&setup->intstep) != 1 || 
	    setup->intstep < 0.0 || setup->intstep > 30.0) {
	fprintf(stderr,"ERROR: setup_interact.  Bad input.\n");
	fprintf(stderr,"Enter step size again:  ");
	fgets(line,MAXC,stdin);
      }
    }
  }
}

/*.......................................................................
 *
 * Function get_smooth_fn
 *
 * Gets the function to be used to smooth the input light curve.
 *
 * Inputs: Setup *setup        container holding setup info
 *
 * Output: none
 *
 */

void get_smooth_fn(Setup *setup)
{
  char line[MAXC];      /* General string for reading input */

  /*
   * Set default value
   */

  setup->smtype = SMBOXCAR;

  /*
   * Get user value
   */

  printf("\nSmoothing function:\n");
  printf("  %d. Boxcar\n",SMBOXCAR);
  printf("  %d. Median\n",SMMEDIAN);
  printf("  %d. Triangle\n",SMTRIANGLE);
  printf("  %d. Gaussian\n",SMGAUSS);
  printf("  %d. Variable width boxcar\n",SMVARBOX);
  printf("  %d. Variable width triangle\n",SMVARTRI);
  printf("Enter choice: [%d] ",setup->smtype);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%d",&setup->smtype) != 1 || setup->smtype < 0
	  || setup->smtype > SMVARTRI) {
      fprintf(stderr,"ERROR: get_smooth_fn. Invalid smoothing type.\n");
      fprintf(stderr,"  Enter value again: ");
      fgets(line,MAXC,stdin);
    }
  }
}

/*.......................................................................
 *
 * Function get_smooth_width
 *
 * Gets the "width" of the function to be used to smooth the input light curve.
 *
 * Inputs: Setup *setup        container holding setup info
 *
 * Output: none
 *
 */

void get_smooth_width(Setup *setup)
{
  char line[MAXC];      /* General string for reading input */

  /*
   * Choose a default width of 10 days
   */

  setup->smwidth = 10.0;

  /*
   * Get smoothing width
   */

  printf("\nEnter width of window for smoothing (in days): [%4.1f] ",
	 setup->smwidth);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%f",&setup->smwidth) != 1 || setup->smwidth < 1.0 || 
	  setup->smwidth > 100.0) {
      fprintf(stderr,"ERROR: get_smooth_width. ");
      fprintf(stderr,"Bad or poorly chosen input (too narrow or too wide).\n");
      fprintf(stderr,"Enter width again:  ");
      fgets(line,MAXC,stdin);
    }
  }
}

/*.......................................................................
 *
 * Function get_nvar
 *
 * Gets the number of points used to set the width of the variable-width
 *  smoothing boxes
 *
 * Inputs: Setup *setup        container holding setup info
 *
 * Output: none
 *
 */

void get_nvar(Setup *setup)
{
  char line[MAXC];      /* General string for reading input */

  /*
   * Choose a default value of 5 days
   */

  setup->nvar = 5;

  /*
   * Get number of points
   */

  printf("\nEnter number of points for smoothing: [%d] ",
	 setup->nvar);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%d",&setup->nvar) != 1 || setup->nvar <= 0 || 
	  setup->nvar > 30) {
      fprintf(stderr,"ERROR: setup_interact.  Bad input.  ");
      fprintf(stderr,"Enter number again:  ");
      fgets(line,MAXC,stdin);
    }
  }

}

/*.......................................................................
 *
 * Function set_grid_params
 *
 * Sets the starting date and number of points for the light curves which
 *  result from interpolating the raw light curves onto a regularly 
 *  sampled grid.
 *
 */

void set_grid_params(Setup *setup, Fluxrec *raw, int nraw)
{
  char line[MAXC];       /* General string for reading input */

  /*
   * Set interpolation start date, if not already set.
   */

  printf("\n");
  if(setup->intstart < 0) {
    setup->intstart = raw->day;
    if(setup->askstart == YES && setup->dosmooth != SMINPLACE) {
      printf("set_grid_params: Enter start date for interpolation: [%7.2f] ",
	     setup->intstart);
      fgets(line,MAXC,stdin);
      if(line[0] != '\n') {
	while(sscanf(line,"%f",&setup->intstart) != 1 || 
	      setup->intstart < 0) {
	  fprintf(stderr,"  ERROR: Bad value.  Enter start date again:  ");
	  fgets(line,MAXC,stdin);
	}
      }
    }
  }

  /*
   * Calculate the size of the final interpolated or smoothed array that 
   *  will come out of all of the smoothing functions except SMINPLACE.
   *  For SMINPLACE option, just set size to nraw.
   */

  if(setup->dosmooth == SMINPLACE)
    setup->ninterp = nraw;
  else if(setup->intstep > 0) {
    setup->ninterp = 1 + (int) (((raw+nraw-1)->day - setup->intstart) / 
      setup->intstep);
  }

  /*
   * Print out results.
   */

  printf("set_grid_params: Grid start date = %7.2f.\n",setup->intstart);
  printf("set_grid_params: N_points = %d\n",setup->ninterp);
  printf("set_grid_params: Grid spacing = %5.2f\n",setup->intstep);

}
