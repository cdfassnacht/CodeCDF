/*
 * getflux.c
 *
 * This program takes the output from the difmap modelfitting scripts for
 *  the lens monitoring program and finds the fluxes associated with the
 *  best fit models.  These fluxes are appended, along with the MMJD, 
 *  to a file in the Analysis directory containing all the
 *  daily fluxes.
 *
 * Usage: getflux [date] [source_root_name]
 *
 * 15Jan97 CDF
 * v19Mar97 CDF, Added parabola fitting to the chisq curves to give a better
 *                estimate of the scale factor associated with the min. chisq
 * v15Jun97 CDF, Changed output file to match file extension of input file.
 * v25May98 CDF, Included fluxes of other sources in 1608 field.
 *               Made separate output files for 1608-phasecal, 1608-gscale,
 *                and the comparison sources (1634 and 1635).
 * v01Jun98 CDF, Added processing of 1608 modelfit chisq values.
 * v20Jul99 CDF, Split comparison source output file (compfluxes.[ext]) into
 *                two files (1634.[ext] and 1635.[ext])
 * v16Sep99 CDF, Modified to read observing date from 1608 log file (in JD
 *                format) rather than calculating MJD from the DATE-OBS
 *                header card, which only gives integer days.
 * v26Feb00 CDF, Only run on 1608 files for now.
 *               Improved output format for STDOUT.
 *               Added output file for difference between gscale and
 *                phasecal flux densities.
 * v14Sep00 CDF, Added rounding of observed MJD to nearest half day.
 * v15Sep00 CDF, Changed input to be date of observation in form yyyymmdd
 *                in order to cut down on useless copying and re-naming
 *                in automapping procedures.
 *               Moved obsdate conversion into new set_log_obsdate function
 *                in modfuncs.c.
 *               Moved printing functions into the new print_lensflux
 *                function in modfuncs.c.
 *               *** Made general enough to also take over functionality
 *                ** of calflux.c, so calflux.c has been deleted. ***
 * v21May01 CDF, Added initialization of some variables that needed it.
 * v02Apr02 CDF, Took out rounding to nearest half day.
 * v14Jul03 CDF, Modified to handle lenses other than 1608.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structdef.h"
#include "dataio.h"
#include "modfuncs.h"

#define N1608 4 /* Number of sources in 1608 field, including 4 components */

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;                  /* Looping variable */
  int no_error=1;         /* Flag set to 0 on error */
  int ncomp=1;            /* Number of components */
  int srctype=CAL;        /* Type of source (lens or cal) */
  float rmsp=0.0;         /* RMS noise from phasecal map */
  float rmsg=0.0;         /* RMS noise from gscale map */
  float chisqp=0.0;       /* Chisq value from phasecal modelfit */
  float chisqg=0.0;       /* Chisq value from gscale modelfit */
  float diff;             /* Difference between gscale and phasecal fluxes */
  double obsdate;         /* Date of observation */
  char root[MAXC];        /* Source root name */
  char datename[MAXC];    /* Observation date as used in input filenames */
  char flag[MAXC];        /* Either 'lens' or 'cal' */
  char baddayname[MAXC];  /* Name for bad day file */
  Secat *fluxp=NULL;      /* Component fluxes from phasecal model */
  Secat *fluxg=NULL;      /* Component fluxes from gscale model */
  Secat *pptr,*gptr;      /* Pointers to navigate flux arrays */
  FILE *dayfp=NULL;       /* Output file pointer for temporary day file */
  FILE *badfp=NULL;       /* Output file pointer for bad day file */

  /*
   * Check command line format
   */

  if(argc != 4) {
    fprintf(stderr,"\n*** Usage: getflux [date] [source_name] [flag] ");
    fprintf(stderr,"***\n\n");
    fprintf(stderr,"The date is the observing date in the form yyyymmdd\n");
    fprintf(stderr,"The source name is how the source is identified in the\n");
    fprintf(stderr,"    input files.  It is often the first part of the IAU");
    fprintf(stderr," name.\n");
    fprintf(stderr,"    e.g., 1608 for 1608+656");
    fprintf(stderr,"The flag should be either lens or cal\n");
    fprintf(stderr,"For example:\n\n");
    fprintf(stderr,"             getflux 20000214 1608 lens\n\n");
    return 1;
  }

  /*
   * Set datename, root, and flag variables from command line
   */

  strcpy(datename,argv[1]);
  strcpy(root,argv[2]);
  strcpy(flag,argv[3]);

  /*
   * Check that the flag variable is correct
   */

  if(strcmp(flag,"lens") == 0)
    srctype = LENS;
  else if(strcmp(flag,"cal") == 0)
    srctype = CAL;
  else {
    fprintf(stderr,"\n\nERROR: Flag must be either 'lens' or 'cal'\n");
    fprintf(stderr,"  You entered %s\n",flag);
    return 1;
  }

  /*
   * Process data, once for phasecal model and once for gscale model
   */

  printf("\n");
  if(no_error)
    if(!(fluxp = source_flux(datename,root,"p",&ncomp,&rmsp,&chisqp,
		    &obsdate)) == 1)
      no_error = 0;

  if(no_error)
    if(!(fluxg = source_flux(datename,root,"g",&ncomp,&rmsg,&chisqg,
		    &obsdate)) == 1)
      no_error = 0;

  /*
   * Convert JD to MJD 
   */

  if(no_error)
    set_log_obsdate(&obsdate);

  /*
   * Print out results to STDOUT
   */

  if(no_error) {
    printf("\nBest-fit flux densities for %s (mJy):\n",root);
    printf("---------------------------------------\n\n");
    printf("   Component  S(phasecal)  S(gscale)   dS(mJy)   dS/S(\%)\n");
    printf("   ---------  -----------  ---------  ---------  -------\n");
    for(i=0,pptr=fluxp,gptr=fluxg; i<ncomp; i++,pptr++,gptr++) {
      diff = gptr->fauto - pptr->fauto;
      printf("   %4d         %7.3f     %7.3f     %6.3f  %7.3f\n",
	     i+1,pptr->fauto,gptr->fauto,diff,100.0*diff / gptr->fauto);
    }
    printf("\n  Phasecal fit parameters: rms = %8.5f mJy/beam, ",rmsp);
    printf("chisq = %8.4f\n",chisqp);
    printf("  Gscale fit parameters:   rms = %8.5f mJy/beam, chisq = %8.4f\n",
	   rmsg,chisqg);
    printf("\n");
  }

  /*
   * Open tmpday and badday output files
   */

  if(srctype == LENS) {
    if(!(dayfp = open_writefile("tmpday")))
      no_error = 0;
    sprintf(baddayname,"%s%s_badday.edt",REDDIR,root);
    if(!(badfp = open_appendfile(baddayname,1)))
      no_error = 0;
  }

  /*
   * Print data to output files
   */

  if(no_error) {
    if(srctype == LENS) {
      if(print_lensflux(root,fluxp,fluxg,ncomp,rmsp,rmsg,chisqp,chisqg,
			obsdate))
	no_error = 0;
      fprintf(dayfp,"DAY %7.2f\n",obsdate);
      fprintf(badfp,"%7.2f 0 %s\n",obsdate,datename);
    }
    else {
      if(print_calflux(root,fluxp,fluxg,ncomp,rmsp,rmsg,chisqp,chisqg,
		       obsdate))
	no_error = 0;
    }
  }

  /*
   * Clean up
   */

  fluxp = del_array(fluxp);
  fluxg = del_array(fluxg);
  if(dayfp)
    fclose(dayfp);
  if(badfp)
    fclose(badfp);

#if 0
  /***********************************************************************
   *
   * Old 1634 and 1635 stuff -- keep just in case we go back to AF310 
   * analysis.
   */

  float flux34;           /* Flux for 1634 best-fit model */
  float flux35;           /* Flux for 1635 best-fit model */
  float rms34,rms35;      /* RMS errors on best-fit models */
  float chisq34,chisq35;  /* chisq values for best-fit models */
  char out1634[MAXC];     /* Output filename for 1634  */
  char out1635[MAXC];     /* Output filename for 1635 */
  FILE *ofp1634=NULL;     /* Output file pointer for 1634  */
  FILE *ofp1635=NULL;     /* Output file pointer for 1635 */

  /*
   * Process 1634 data
   */

  if((get_best_cal("1634",ext,&flux34,&rms34,&chisq34)) == 1)
    no_error = 0;

  /*
   * Process 1635 data
   */

  if(no_error)
    if((get_best_cal("1635",ext,&flux35,&rms35,&chisq35)) == 1)
      no_error = 0;

  /*
   * Print out results to STDOUT
   */

  if(no_error) {
    printf("\nBest fit 1634 flux = %9.5f, rms = %10.3e, chisq = %8.3f\n",
	   flux34,rms34,chisq34);
    printf("Best fit 1635 flux = %9.5f, rms = %10.3e, chisq = %8.3f\n",
	   flux35,rms35,chisq35);
  }

  /*
   * Set the output filenames
   */

  sprintf(out1634,"%s1634.edt",REDDIR);
  sprintf(out1635,"%s1635.edt",REDDIR);

  /*
   * Open the output files
   */

  if((ofp1634 = fopen(out1634,"a")) == NULL) {
    fprintf(stderr,"ERROR.  Cannot open output file %s\n",out1634);
    no_error = 0;
  }
  if((ofp1635 = fopen(out1635,"a")) == NULL) {
    fprintf(stderr,"ERROR.  Cannot open output file %s\n",out1635);
    no_error = 0;
  }

  /*
   * Print data to output files
   */

  if(no_error) {
    fprintf(ofp1634,"%6.1f %8.3f %8.5f %8.3f\n",
	    obsdate,flux34,rms34,chisq34);
    fprintf(ofp1635,"%6.1f %8.3f %8.5f %8.3f\n",
	    obsdate,flux35,rms35,chisq35);
  }

  if(ofp1634)
    fclose(ofp1634);
  if(ofp1635)
    fclose(ofp1635);

  /*
   ***********************************************************************
   */
#endif

  if(no_error) {
    printf("\nFinished with program getflux.c.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"ERROR: Exiting program getflux.c\n\n");
    return 1;
  }
}

