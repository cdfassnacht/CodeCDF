/*
 * setup_spec.c
 *
 * A library of functions for obtaining and storing plotting parameters
 *  used in plotting one-dimensional spectra (e.g., with specplot.c).
 *
 * 13Sep01 CDF   Split off from onedspec.c
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

/*.......................................................................
 *
 * Function new_spsetup
 *
 * Allocates memory for a Spsetup structure array.
 *
 * Output: Spsetup *newspsetup  new spsetup array
 *
 */

Spsetup *new_spsetup()
{
  Spsetup *newspsetup;

  newspsetup = (Spsetup *) malloc(sizeof(Spsetup));
  if(!newspsetup) {
    fprintf(stderr,"Insufficient memory for Spsetup array.\n");
    return NULL;
  }

  /*
   * Initialize the spsetup parameters.
   */

  newspsetup->option = ONESPEC;
  newspsetup->flchoice = 3;
  newspsetup->labstyle = 1;
  newspsetup->labyfrac = 0.02;
  newspsetup->labheight = 1.0;
  newspsetup->cheight = 1.2;
  newspsetup->pltxmin = 0.0;
  newspsetup->pltxmax = 0.0;

  return newspsetup;
}

/*.......................................................................
 *
 * Function del_spsetup
 *
 * Frees memory associated with Spsetup array
 *
 * Input:  Spsetup *spsetup        array to be freed
 *
 * Output: NULL
 *
 */

Spsetup *del_spsetup(Spsetup *spsetup)
{
  if(spsetup)
    free(spsetup);

  return NULL;
}

/*.......................................................................
 *
 * Function get_flchoice
 *
 * Sets flag if fluxes are in magnitudes.
 *
 * Input:  Spsetup *setup     container for plotting parameters.
 *
 * Output: (none)
 *
 */

void get_flchoice(Spsetup *setup)
{
  char line[MAXC];  /* General string for reading input */

  printf("\nFlux units for the spectrum:\n");
  printf("  1. Magnitudes\n");
  printf("  2. AB Magnitudes\n");
  printf("  3. Anything else (relative flux, ergs/cm^2/sec/Ang, etc.)\n");
  printf("Enter choice:  [%d] ",setup->flchoice);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%d",&setup->flchoice) != 1 || 
	  setup->flchoice < 1 || setup->flchoice > 3) {
      fprintf(stderr,"ERROR: invalid choice.  Enter choice again:  ");
      fgets(line,MAXC,stdin);
    }
  }
}

/*.......................................................................
 *
 * Function get_cheight
 *
 * Gets height of characters used in external labesl.
 *
 * Input:  Spsetup *setup     container for plotting parameters.
 *
 * Output: (none)
 *
 */

void get_cheight(Spsetup *setup)
{
  char line[MAXC];  /* General string for reading input */

  printf("\nEnter character height for external plot labels: [%3.1f]  ",
	 setup->cheight);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%f",&setup->cheight) != 1 || 
	  setup->cheight < 0.0) {
      fprintf(stderr,"ERROR: Bad input format.  Enter character height: ");
      fgets(line,MAXC,stdin);
    }
  }
}

/*.......................................................................
 *
 * Function get_labheight
 *
 * Gets height o labels for spectral lines.
 *
 * Input:  Spsetup *setup     container for plotting parameters.
 *
 * Output: (none)
 *
 */

void get_labheight(Spsetup *setup)
{
  char line[MAXC];  /* General string for reading input */

  printf("\nEnter character height for spectral line labels: [%3.1f]  ",
	 setup->labheight);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%f",&setup->labheight) != 1 || 
	  setup->labheight < 0.0) {
      fprintf(stderr,"ERROR: Bad input format.  Enter character height: ");
      fgets(line,MAXC,stdin);
    }
  }
}

/*.......................................................................
 *
 * Function get_label_style
 *
 * Gets style of labeling spectral lines.
 *
 * Inputs: Spsetup *setup     container for plotting parameters.
 *         int option          plotting option (e.g. MULTISPEC)
 *
 * Output:  int (0 or 1)       0 ==> success, 1 ==> failure
 *
 */

int get_label_style(Spsetup *setup, int option)
{
  char line[MAXC];     /* General string variable */

  /*
   * If plotting in the MULTISPEC style, automatically use the vertical
   *  dotted line label style.
   */

  if(option == MULTISPEC) {
    setup->labstyle = 2;
  }

  /*
   * Otherwise, query the user about the line label style.
   */

  else {
    printf("\nLine label style options:\n");
    printf("  1. Short tick marks\n");
    printf("  2. Vertical dotted lines covering height of plot\n");
    printf("--------------------------------------------------\n");
    printf(" Enter choice: [%d] ",setup->labstyle);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%d",&setup->labstyle) != 1 || 
	    setup->labstyle < 1 || setup->labstyle > 2) {
	fprintf(stderr,"ERROR: Invalid input.  Enter choice:  ");
	fgets(line,MAXC,stdin);
      }
    }
  }

  /*
   * If using the vertical dotted line-style, find out where to put the
   *  line labels vertically.
   */

  if(setup->labstyle == 2) {
    printf("\nEnter the position where the line labels should appear\n");
    printf(" vertically, as a fraction of the plot height.  For example,\n");
    printf(" 0.0 would be at the very bottom of the plot, 0.5 would be\n");
    printf(" half-way up, and 1.0 would be the very top of the plot. ");
    printf(" [%5.3f]  ",setup->labyfrac);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%f",&setup->labyfrac) != 1 || 
	    setup->labyfrac < 0.0 || setup->labyfrac > 1.0) {
	fprintf(stderr,"ERROR: Enter value between 0 and 1:  ");
	fgets(line,MAXC,stdin);
      }
    }
  }

  return 0;
}

/*.......................................................................
 *
 * Function get_xrange
 *
 * Sets desired display range for wavelengths.
 *
 * Input:  Spsetup *setup     container for plotting parameters.
 *         Spectrum *spec     spectra
 *
 * Output: (none)
 *
 */

void get_xrange(Spsetup *setup, Spectrum *spec)
{
  char line[MAXC];     /* General string for reading input */

  /*
   * Set default wavelength range to be that occupied by the first spectrum.
   */

  setup->pltxmin = spec->minlambda;
  setup->pltxmax = spec->maxlambda;

  /*
   * See if a different range is desired.
   */

  printf("\nEnter wavelength range to be displayed: [%7.1f %7.1f]  ",
	 setup->pltxmin,setup->pltxmax);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%lf %lf",&setup->pltxmin,&setup->pltxmax) != 2) {
      fprintf(stderr,"ERROR: Bad input format.  Enter wavelength range: ");
      fgets(line,MAXC,stdin);
    }
  }
}

/*.......................................................................
 *
 * Function get_title
 *
 * Gets optional title for plot and location for title.
 *
 * Inputs: char *title         title for plot (set by this function)
 *         int *titleopt       title location (set by this function)
 *
 * Output:  int (0 or 1)       0 ==> success, 1 ==> failure
 *
 */

int get_title(char *title, int *titleopt)
{
  char line[MAXC];    /* General string variable */

  /*
   * Get title location
   */

  printf("\nPlot title options:\n");
  printf("  0. No title\n");
  printf("  1. Title on top of plot OUTSIDE box\n");
  printf("  2. Title on top of plot INSIDE box\n");
  printf("-------------------------------------\n");
  printf(" Enter choice [%d]:  ",*titleopt);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%d",titleopt) != 1 || *titleopt < 0 || 
	  *titleopt > 2) {
      fprintf(stderr,"ERROR: invalid input.  Enter new value:  ");
      fgets(line,MAXC,stdin);
    }
  }

  /*
   * Now get the title string
   */

  if(*titleopt == 0)
    title[0] = '\0';
  else {
    printf("\nEnter title for plot:  ");
    fgets(title,MAXC,stdin);

    /*
     * Strip the newline out of the title array, replacing it by a \0
     */
    
    title[strlen(title)-1] = '\0';
  }

  return 0;
}

/*.......................................................................
 *
 * Function get_clabel
 *
 * Gets label text for optional corner label, and also gets label
 *  position (left or right lower corner).
 *
 * Inputs: char *labtext       text of label (set by this function)
 *         char *labside       position of label (set by this function)
 *
 * Outputs: int dolabel       1 ==> label desired. -1 ==> no label
 *                              desired or error.
 *
 */

int get_clabel(char *labtext, char *labside)
{
  char line[MAXC];  /* Line of text to put in corner */

  /*
   * See if label is desired and, if so, get position.
   */

  printf("\nLabel in lower right or left corner?\n");
  printf("Enter l for left corner, r for right corner, \n");
  printf(" or just hit return for no label: [no label] ");
  fgets(line,MAXC,stdin);
  if(line[0] == 'r' || line[0] == 'R')
    *labside = 'r';
  else if(line[0] == 'l' || line[0] == 'L')
    *labside = 'l';
  else
    return -1;

  /*
   * If label is desired, get text.
   */

  printf(" Enter label for corner:  ");
  fgets(line,MAXC,stdin);
  if(line[0] == '\n')
    return -1;
  else {
    line[strlen(line)-1] = '\0';
    strcpy(labtext,line);
  }

  return 1;
}

