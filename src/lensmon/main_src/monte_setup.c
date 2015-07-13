/*
 * monte_setup.c
 *
 * A library of functions to perform operations on the Setup structure
 *  used by the Monte-Carlo-related programs.
 *
 * 2002Apr01 CDF,  A revision of lc_setup.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "lc_funcs.h"
#include "monte_setup.h"


/*.......................................................................
 *
 * Function new_MC_Setup
 *
 * Allocates memory for a Setup structure array.
 *
 * Input:  int size            size of array
 *
 * Output: MC_Setup *newsetup  new setup array
 *
 */

MC_Setup *new_MC_Setup(int size)
{
  MC_Setup *newsetup;

  newsetup = (MC_Setup *) malloc((sizeof(MC_Setup) * size));
  if(!newsetup) {
    fprintf(stderr,"\nnew_setup: Insufficient memory for Setup array.\n");
    return NULL;
  }

  /*
   * Initialize the setup parameters
   */

  newsetup->mu0[0] = 0.0;
  newsetup->mu0[1] = 0.0;
  newsetup->mu0[2] = 0.0;
  newsetup->mu0[3] = 0.0;
  newsetup->nmu = 0;
  newsetup->tau0[0] = 0.0;
  newsetup->tau0[1] = 0.0;
  newsetup->tau0[2] = 0.0;
  newsetup->tau0[3] = 0.0;
  newsetup->ntau = 0;
  newsetup->nideal = 0;
  newsetup->nobs = 0;
  newsetup->mcsig[0] = 0.0;
  newsetup->mcsig[1] = 0.0;
  newsetup->mcsig[2] = 0.0;
  newsetup->mcsig[3] = 0.0;
  newsetup->ftol = 0.0;

  return newsetup;
}

/*.......................................................................
 *
 * Function del_MC_Setup
 *
 * Frees memory associated with Setup array
 *
 * Input:  MC_Setup *MC_Setup    array to be freed
 *
 * Output: NULL
 *
 */

MC_Setup *del_MC_Setup(MC_Setup *mcsetup)
{
  if(mcsetup)
    free(mcsetup);

  return NULL;
}

/*.......................................................................
 *
 * Function MC_Setup_file
 *
 * Fills in MC_Setup container with information from a file.
 *
 * Inputs: char *inname        name of file containing info
 *         MC_Setup *setup     setup structure to be filled
 *
 * Output:int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int MC_Setup_file(MC_Setup *setup, char *inname)
{
  int no_error=1;     /* Flag set to 0 on error */
  char keyword[MAXC]; /* Keyword at the beginning of each input line */
  char line[MAXC];    /* General string for reading input */
  FILE *ifp=NULL;     /* Input file pointer */

  /*
   * Open input file
   */

  while(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: MC_Setup_file.\n");
    return 1;
  }

  /*
   * Read input from setup file.
   */

  printf("\nReading MC setup info from file %s\n\n",inname);
  while(fgets(line,MAXC,ifp) != NULL && no_error) {
    if(line[0] >= 32 && line[0] != '#') {
      switch(read_mcsetup_line(line,keyword)) {
      case MCSETUPERR:
	no_error = 0;
	break;
      case MCMU0:
	if(sscanf(line,"%s %f %f %f %f",keyword,&setup->mu0[0],&setup->mu0[1],
		  &setup->mu0[2],&setup->mu0[3]) != 5) {
	  fprintf(stderr,"ERROR: MC_Setup_file.  Bad inputs for m0\n");
	  fprintf(stderr," Setting mu0 = {0 0 0 0}\n");
	  setup->mu0[0] = 0.0;
	  setup->mu0[1] = 0.0;
	  setup->mu0[2] = 0.0;
	  setup->mu0[3] = 0.0;
	}
	break;
      case MCNMU:
	if(sscanf(line,"%s %d",keyword,&setup->nmu) != 2 ||
	   setup->nmu < 0) {
	  fprintf(stderr,"ERROR: MC_Setup_file.  Bad input for nmu\n");
	  fprintf(stderr," Setting nmu = 0\n");
	  setup->nmu = 0;
	}
	break;
      case MCTAU0:
	if(sscanf(line,"%s %f %f %f %f",keyword,&setup->tau0[0],
		  &setup->tau0[1],&setup->tau0[2],&setup->tau0[3]) != 5) {
	  fprintf(stderr,"ERROR: MC_Setup_file.  Bad inputs for tau0\n");
	  fprintf(stderr," Setting tau0 = {0 0 0 0}\n");
	  setup->tau0[0] = 0.0;
	  setup->tau0[1] = 0.0;
	  setup->tau0[2] = 0.0;
	  setup->tau0[3] = 0.0;
	}
	break;
      case MCNTAU:
	if(sscanf(line,"%s %d",keyword,&setup->ntau) != 2 ||
	   setup->ntau < 0) {
	  fprintf(stderr,"ERROR: MC_Setup_file.  Bad input for ntau\n");
	  fprintf(stderr," Setting ntau = 0\n");
	  setup->ntau = 0;
	}
	break;
      case MCSIG:
	if(sscanf(line,"%s %f %f %f %f",keyword,&setup->mcsig[0],
		  &setup->mcsig[1],&setup->mcsig[2],&setup->mcsig[3]) != 5) {
	  fprintf(stderr,"ERROR: MC_Setup_file.  Bad inputs for mcsig\n");
	  fprintf(stderr," Setting mcsig = {0 0 0 0}\n");
	  setup->mcsig[0] = 0.0;
	  setup->mcsig[1] = 0.0;
	  setup->mcsig[2] = 0.0;
	  setup->mcsig[3] = 0.0;
	}
	break;
      case MCROOT:
	if(sscanf(line,"%s %s",keyword,setup->root) != 2) {
	  fprintf(stderr,
		  "ERROR: MC_Setup_file.  Bad value for root.\n");
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
    fprintf(stderr,"ERROR: MC_Setup_file\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function read_mcsetup_line
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

int read_mcsetup_line(char *line, char *keyword)
{
  /*
   * Check to make sure that there is a keyword on the line.
   */

  if(sscanf(line,"%s",keyword) != 1) {
    fprintf(stderr,"ERROR: read_mcsetup_line\n");
    keyword = NULL;
    return SETUPERR;
  }

  /*
   * If there is a keyword, read it and return the appropriate value
   */

  if(strcmp(keyword,"mu0") == 0 || strcmp(keyword,"MU0") == 0)
    return MCMU0;
  if(strcmp(keyword,"nmu") == 0 || strcmp(keyword,"NMU") == 0)
    return MCNMU;
  if(strcmp(keyword,"tau0") == 0 || strcmp(keyword,"TAU0") == 0)
    return MCTAU0;
  if(strcmp(keyword,"ntau") == 0 || strcmp(keyword,"NTAU") == 0)
    return MCNTAU;
  if(strcmp(keyword,"mcsig") == 0 || strcmp(keyword,"MCSIG") == 0)
    return MCSIG;
  if(strcmp(keyword,"root") == 0 || strcmp(keyword,"ROOT") == 0)
    return MCROOT;

  /*
   * If none of the above checks have been satisfied, return the
   *  default value
   */

  return DEFAULT;
}

/*.......................................................................
 *
 * Function setup_monte
 *
 * Fills in the unfilled (by MC_Setup_file) parts of the setup structure 
 *  relating to the algorithms that generate the Monte Carlo curves.
 *
 * Inputs: MC_Setup *setup     setup structure
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int setup_monte(MC_Setup *setup)
{
  if(setup->mu0[0] == 0.0) {
    if(get_mu0_mc(setup)) {
      fprintf(stderr,"ERROR: setup_monte.\n");
      return 1;
    }
  }

  if(setup->tau0[0] == 0.0) {
    if(get_tau0_mc(setup)) {
      fprintf(stderr,"ERROR: setup_monte.\n");
      return 1;
    }
  }

  if(setup->mcsig[0] == 0.0) {
    if(get_mcsig(setup)) {
      fprintf(stderr,"ERROR: setup_monte.\n");
      return 1;
    }
  }

  return 0;
}

/*.......................................................................
 *
 * Function setup_monte_summary
 *
 * Summarizes the information in the setup container relating to the
 *  Monte Carlo algorithms.
 *
 * Inputs: MC_Setup *setup     setup container
 *
 * Outputs: (none)
 *
 */

void setup_monte_summary(MC_Setup *setup)
{
  int i; /* Looping varaible */

  printf("\n------------------------------------------------------------\n");
  printf("\nMonte Carlo information:\n");
  printf(" Delays:         %6.2f %6.2f %6.2f %6.2f\n",
	 setup->tau0[0],setup->tau0[1],setup->tau0[2],setup->tau0[3]);
  printf(" Magnifications: ");
  for(i=0; i<4; i++)
    printf("%6.4f ",setup->mu0[i]);
  printf("\n Sigmas:         %6.4f %6.4f %6.4f %6.4f\n",
	 setup->mcsig[0],setup->mcsig[1],setup->mcsig[2],setup->mcsig[3]);
  printf(" Output file root name: %s\n",setup->root);
  printf("\n------------------------------------------------------------\n");
}

/*.......................................................................
 *
 * Function get_mu0_mc
 *
 * Gets information about the magnifications used.
 *
 * Inputs: MC_Setup *setup     setup structure
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int get_mu0_mc(MC_Setup *setup)
{
  char line[MAXC];  /* General string variable for reading input */

  printf("\nget_mu0_mc: The four image magnifications haven't been set.\n");
  printf("get_mu0_mc: Enter values all on this line: ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f %f %f %f",&setup->mu0[0],&setup->mu0[1],
	       &setup->mu0[2],&setup->mu0[3]) != 4) {
    fprintf(stderr,"ERROR: get_mu0_mc.  Bad inputs for mu0\n");
    fprintf(stderr,"Enter values again:  ");
    fgets(line,MAXC,stdin);
  }

  return 0;
}


/*.......................................................................
 *
 * Function get_tau0_mc
 *
 * Gets information about the delays used.
 *
 * Inputs: MC_Setup *setup     setup structure
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int get_tau0_mc(MC_Setup *setup)
{
  char line[MAXC];  /* General string variable for reading input */

  printf("\nget_tau0_mc: The four time delays haven't been set.\n");
  printf("get_tau0_mc: Enter values all on this line: ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f %f %f %f",&setup->tau0[0],&setup->tau0[1],
	       &setup->tau0[2],&setup->tau0[3]) != 4) {
    fprintf(stderr,"ERROR: get_tau0_mc.  Bad inputs for tau0\n");
    fprintf(stderr,"Enter values again:  ");
    fgets(line,MAXC,stdin);
  }

  return 0;
}

/*.......................................................................
 *
 * Function get_mcsig
 *
 * Gets information about the width of the random gaussian distribution
 *  used to generated the realizations.
 *
 * Inputs: MC_Setup *setup     setup structure
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int get_mcsig(MC_Setup *setup)
{
  char line[MAXC];  /* General string variable for reading input */

  printf("\nget_mcsig: The four values of sigma used to create the\n");
  printf("get_mcsig:  haven't been set.\n");
  printf("get_mcsig: Enter values all on this line: ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f %f %f %f",&setup->mcsig[0],&setup->mcsig[1],
	       &setup->mcsig[2],&setup->mcsig[3]) != 4) {
    fprintf(stderr,"ERROR: get_mcsig.  Bad inputs for mcsig\n");
    fprintf(stderr,"Enter values again:  ");
    fgets(line,MAXC,stdin);
  }

  return 0;
}

/*.......................................................................
 *
 * Function get_nout
 *
 * Gets the number of realizations to create.
 *
 * Inputs: MC_Setup *setup     setup structure
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> failure
 *
 */

int get_nout(MC_Setup *setup)
{
  char line[MAXC];  /* General string variable for reading input */

  printf("Enter number of Monte Carlo realizations to create: [%d] ",
	 setup->nout);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%d",&setup->nout) != 1 && setup->nout < 0) {
      fprintf(stderr,"ERROR: Bad input value.  Enter number again:  ");
      fgets(line,MAXC,stdin);
    }
  }

  return 0;
}


