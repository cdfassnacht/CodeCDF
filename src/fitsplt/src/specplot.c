/*
 * specplot.c
 *
 * Plots a spectrum from an ascii file containing two columns, the wavelength
 *  and the flux.  The user enters a redshift or redshifts and then labels
 *  lines in the spectrum.  The label position is computed by redshifting the
 *  rest wavelength of the line by the appropriate ammount
 *
 * Usage: specplot [input_file]
 *
 * 20Jul96 CDF
 * v05Jul97 CDF, Put more processing into functions.  Better documentation
 *               Mark lines with vertical dashed lines
 * v07Jul97 CDF, Change line marking back to short line segments.  Fixed
 *               bug in displaying line list choices.
 * v30Jul97 CDF, Better error checking. 
 *               Added option to reprint line list.
 *               Added option to shift labels to avoid overlaps.
 * v15Dec97 CDF, Added option to plot fluxes which are on a magnitude scale.
 *               Combined functionality of lambda and flux arrays into a
 *                Pos array since the Pos (x,y) structure is equally 
 *                appropriate for (lambda,flux) pairs.
 *               Added option to plot a second spectrum, e.g. an error
 *                spectrum or a model.
 * v17Dec97 CDF, User now can choose limits for flux axis on plot.
 * v23Feb98 CDF, Interactive setting of character height for plot labels.
 * v25Apr99 CDF, Modified several functions to make "official" the previously 
 *                "hidden" options for plotting sky lines and source lines not
 *                 on the list.
 *               Added label_sky_lines function.
 * v30Apr99 CDF, Added interactive adjustment of vertical scale.
 * v01May99 CDF, Modified scale-adjustment routine.
 * v03Aug99 CDF, Allowed more flexibility in placement of title.
 *               Added choice of short line segments or vertical dotted lines
 *                to label spectral lines.
 *               Added choice of converting spectrum from f_lambda units on y
 *                axis to f_nu.
 *               Added choice of placing internal label in lower left or lower
 *                right corners.
 * v06Aug99 CDF, Fixed speclab to allow shifting of user-entered lines.
 * v31Aug01 CDF, Split off many functions into onedspec.c.
 * v14Sep01 CDF, Major rewrite, including:
 *                1) Addition of multispec and overplotting options.
 *                2) Introduction of Spsetup container to hold many of the
 *                   variables that control the plotting.
 * v14Feb03 CDF, Moved new_spec and del_spec functions into onedspec.c
 * v18Feb03 CDF, Major modifications resulting from new Spectrum structure
 *                (a combination of the old Pos spectrum and the
 *                Splotinfo structure).  Moved all file I/O into new
 *                read_onespec and read_multispec functions in onedspec.c
 * v08Sep03 CDF, Change to have different calls to the plot_open_* functions
 *                in the hardcopy section, depending on whether a
 *                single spectrum or multiple spectra are being plotted.
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
 * Function declarations
 *
 */

void specplot_help();
int specplot_parse(char *optchoice);
int display_multispec(Spectrum *spectrum, int nspec, int maxspec, 
		      Spsetup *setup, int option);
int adjust_multispec(Spectrum *spectrum, int nspec, 
		     Spsetup *setup, int option);

/*.......................................................................
 *
 * Main Program
 */

int main(int argc, char *argv[])
{
  int i;                   /* Looping variable */
  int option;              /* Command line run option */
  int no_error=1;          /* Flag set to 0 on error */
  int nspec=0;             /* Number of input spectra */
  int nlines=0;            /* Number of lines in input file */
  int goodlines=0;         /* Number of lines in spectral info file */
  int dolabel=0;           /* Flag set to 1 if internal label is desired */
  int titleopt=0;          /* Variable determining location of plot title */
  int fileflag=0;          /* Flag set to 1 if output ps file is desired */
  int nz=1;                /* Number of redshifts used */
  float z;                 /* Redshifts of lines */
  float offset=0.0;        /* Flux offset to apply to second spectrum */
  char labpos;             /* Position (l or r) of optional intern. label */
  char title[MAXC];        /* Optional title for plot */
  char labtext[MAXC];      /* Optional internal label */
  char plotfile[MAXC];     /* Name of optional output postscript file */
  char plotdev[MAXC];      /* PGPLOT device */
  char infile[MAXC];       /* Generic input file name */
  char line[MAXC];         /* General string variable */
  Spectrum *spectrum=NULL; /* Spectrum/spectra to be plotted */
  Spectrum *sptr;          /* Pointer for navigating spectra */
  Specinfo *specinfo=NULL; /* Structure containing spectral line info */
  Spsetup *setup=NULL;     /* Structure with general plotting info */

  /*
   * Check command line format
   */

  if(argc < 3) {
    specplot_help();
    return 1;
  }
  else
    option = specplot_parse(argv[1]);

  if(option == BADOPT){
    specplot_help();
    return 1;
  }

  /*
   * Get input filename
   */

  strcpy(infile,argv[2]);

  /*
   * Read in spectra
   */

  if(no_error) {
    switch(option) {
    case ONESPEC:
      nspec = 1;
      if(!(spectrum = read_onespec(infile)))
	no_error = 0;
      break;
    case MULTISPEC: case OVERPLOT:
      if(!(spectrum = read_multispec(infile,&nspec)))
	no_error = 0;
      break;
    default:
      no_error = 0;
      break;
    }
  }

#if 0
  if(argc == 3) {
    if(!(spectrum2 = read_spectrum(argv[2],&nlines2,minflux2,maxflux2))) {
      no_error = 0;
    }
    else {
      printf("Enter offset to be applied to second spectrum [0.0]: ");
      fgets(line,MAXC,stdin);
      if(line[0] != '\n') {
	while(sscanf(line,"%f",&offset) != 1) {
	  fprintf(stderr,"ERROR: Bad input format.  Enter offset  ");
	  fgets(line,MAXC,stdin);
	}
      }

      /*
       * Add offset
       */

      if(offset != 0.0) {
	printf("Adding offset to second spectrum..\n");
	for(i=0,sptr=spectrum2; i<nlines2; i++,sptr++)
	  sptr->y += offset;
      }
    }
  }
#endif

  /*
   * Allocate memory for spsetup structure
   */

  if(no_error)
    if(!(setup = new_spsetup()))
      no_error = 0;
  /*
   * Get units for fluxes
   */

  if(no_error)
    get_flchoice(setup);

  /*
   * As long as vertical axis is not in magnitudes, give user the
   *  choice of converting to f_nu
   */
#if 0
  if(no_error && setup->flchoice == 3) {
    printf("\nChange spectrum from f_lambda to f_nu? [n] ");
    fgets(line,MAXC,stdin);
    if(line[0] == 'y' || line[0] == 'Y')
      if(flambda_to_fnu(spectrum[0],nlines,minflux,maxflux))
	no_error = 0;
  }
#endif
  /*
   * Get the character height for external labels and line labels
   */

  if(no_error) {
    get_cheight(setup);
    get_labheight(setup);
  }

  /*
   * Choose line-labeling style
   */

  if(no_error)
    if(get_label_style(setup,option))
      no_error = 0;

  /*
   * Get wavelength range to display
   */

  if(no_error)
    get_xrange(setup,spectrum);

  /*
   * Plot the spectrum/spectra, and give the user a chance to adjust
   *  the scale.
   */

  if(no_error) {
    printf("\n");

    /*
     * Call for one spectrum or overplotted spectra
     */

    if(option == ONESPEC || option == OVERPLOT) {
      strcpy(plotdev,getenv("PGPLOT_DEV"));
      if(plot_open_spec(setup->cheight,plotdev))
	no_error = 0;
      if(no_error)
	if(display_spectrum(spectrum,nspec,0,setup,option,1))
	  no_error = 0;
      if(no_error)
	if(adjust_spec(spectrum,nspec,0,setup,option,1))
	  no_error = 0;
    }

    /*
     * Call for multispec -- complicated use of viewports
     */

    if(option == MULTISPEC) {
      if(plot_open_gen(setup->cheight,"/aqt"))
	no_error = 0;
      if(no_error)
	if(adjust_multispec(spectrum,nspec,setup,option))
	  no_error = 0;
    }

  }

  /*
   * Get the number of redshifts
   */

  if(no_error) {
    printf("\nNumber of redshifts used in line ID's [%d]:  ",nz);
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%d",&nz) != 1 || nz < 0) {
	fprintf(stderr,"ERROR: Bad input format for number of redshifts.\n");
	fprintf(stderr,"Enter number of redshifts:  ");
	fgets(line,MAXC,stdin);
      }
    }
  }

  /*
   * Read in list of rest wavelengths
   */

  if(no_error)
    if(!(specinfo = read_line_list(&goodlines)))
      no_error = 0;

  /*
   * Label the lines
   */

  if(no_error && nz>0) {
    for(i=0; i<nz; i++) {
      printf("\nEnter redshift %d:  ",i+1);
      fgets(line,MAXC,stdin);
      while(sscanf(line,"%f",&z) != 1 || z < 0.0) {
	fprintf(stderr,"ERROR: Bad input format\n");
	printf("Enter redshift %d:  ",i+1);
	fgets(line,MAXC,stdin);
      }
      print_line_list(specinfo,goodlines);
      if(label_lines(specinfo,goodlines,z,spectrum,setup))
	no_error = 0;
    }
  }

  /*
   * Label the plot
   */

  if(no_error) {
    if(get_title(title,&titleopt))
      no_error = 0;
    else
      label_plot(title,titleopt,setup->flchoice,setup->cheight);
  }

  /*
   * See if output postscript file is desired and get filename
   *  if it is.
   */

  printf("\nPlot image to a postscript file? [y] ");
  fgets(line,MAXC,stdin);
  if(line[0] == 'n' || line[0] == 'N') {
    fileflag = 0;
  }
  else {
    printf("Enter output filename:  ");
    fgets(line,MAXC,stdin);
    if(line[0] == '\n') {
      fprintf(stderr,"No filename -- no output file will be produced.\n");
      fileflag = 0;
    }
    else {
      while(sscanf(line,"%s",plotfile) != 1) {
	fprintf(stderr,"ERROR: Bad input.  Enter filename again:  ");
	fgets(line,MAXC,stdin);
      }
      fileflag = 1;
    }
  }

  /*
   * If output file is desired, close xwindow and make plot.
   */

  if(fileflag && no_error) {
    plot_close();
    printf("\nWriting to file %s\n",plotfile);
    sprintf(plotdev,"%s/vps",plotfile);
    if(option == ONESPEC || option == OVERPLOT) {
      plot_open_spec(setup->cheight,plotdev);
      if(display_spectrum(spectrum,nspec,0,setup,option,1))
	no_error = 0;
    }
    if(option == MULTISPEC) {
      plot_open_gen(setup->cheight,plotdev);
      if(display_multispec(spectrum,nspec,nspec,setup,option))
	no_error = 0;
    }
    if(nz>0)
      if(label_lines_file(spectrum,setup))
	no_error = 0;
    label_plot(title,titleopt,setup->flchoice,setup->cheight);
  }

  /*
   * See if output gif file is desired and get filename
   *  if it is.
   */

  printf("\nPlot image to a GIF file? [y] ");
  fgets(line,MAXC,stdin);
  if(line[0] == 'n' || line[0] == 'N') {
    fileflag = 0;
  }
  else {
    printf("Enter output filename:  ");
    fgets(line,MAXC,stdin);
    if(line[0] == '\n') {
      fprintf(stderr,"No filename -- no output file will be produced.\n");
      fileflag = 0;
    }
    else {
      while(sscanf(line,"%s",plotfile) != 1) {
	fprintf(stderr,"ERROR: Bad input.  Enter filename again:  ");
	fgets(line,MAXC,stdin);
      }
      fileflag = 1;
    }
  }

  /*
   * If output file is desired, close xwindow and make plot.
   */

  if(fileflag && no_error) {
    plot_close();
    printf("\nWriting to file %s\n",plotfile);
    sprintf(plotdev,"%s/gif",plotfile);
    if(option == ONESPEC || option == OVERPLOT) {
      plot_open_spec(setup->cheight,plotdev);
      if(display_spectrum(spectrum,nspec,0,setup,option,1))
	no_error = 0;
    }
    if(option == MULTISPEC) {
      plot_open_gen(setup->cheight,plotdev);
      if(display_multispec(spectrum,nspec,nspec,setup,option))
	no_error = 0;
    }
    if(nz>0)
      if(label_lines_file(spectrum,setup))
	no_error = 0;
    label_plot(title,titleopt,setup->flchoice,setup->cheight);
  }

  /*
   * Clean up
   */

  plot_close();
  setup = del_spsetup(setup);
  specinfo = del_specinfo(specinfo);
  spectrum = del_spectrum(spectrum);

  if(no_error) {
    printf("\nFinished with program specplot.c\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR: Exiting program specplot.c\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function specplot_help
 *
 * Prints out helpful information when the program is called with no 
 *  arguments.
 *
 * Inputs: none
 *
 * Output: none
 *
 */

void specplot_help()
{
  printf("\n");
  printf("*****************************************************************");
  printf("\n\n");
  printf("\nProgram: specplot -- ");
  printf("Plots spectra and labels lines.\n\n");
  printf("Usage: specplot -1 input_file\n");
  printf("       specplot -o input_list_file\n");
  printf("       specplot -m input_list_file\n\n");
  printf("OPTIONS:\n");
  printf(" -1   Plot one spectrum only.  Input file is the ASCII spectrum.");
  printf("\n\n");
  printf(" -o   Overplot two spectra, where the second one is plotted in\n");
  printf("      a dashed line.  For this option, the input file consists of\n");
  printf("      a two-element list containing the two file names of the\n");
  printf("      ASCII spectra, one per line.\n\n");
  printf(" -m  Plot multiple spectra, each in its own box.\n");
  printf("      For this option, the input file consists of a list \n");
  printf("      containing the filenames of the ASCII spectra, one per \n");
  printf("      line.  All spectra will be marked with the same chosen\n");
  printf("      set of spectral lines.\n\n");
  printf("*****************************************************************");
  printf("\n\n");
}

/*.......................................................................
 *
 * Function specplot_parse
 *
 * Sets the plotting option depending on the second element of the command
 *  line (e.g., '-m')
 *
 * Inputs: char *optchoice     option entered on command line.
 *
 * Output: int option          one of option enumeration values.
 *
 */

int specplot_parse(char *optchoice) {
  if(strcmp(optchoice,"-1") == 0)
    return ONESPEC;
  else if(strcmp(optchoice,"-o") == 0)
    return OVERPLOT;
  else if(strcmp(optchoice,"-m") == 0)
    return MULTISPEC;
  else
    return BADOPT;
}

/*.......................................................................
 *
 * Function: display_multispec
 *
 * Plots multiple spectra, each in its own window.  Not all of the spectra
 *  need to be displayed.  The total number displayed is set by the passed
 *  maxspec variable.
 *
 * Inputs: Spectrum *spectrum  Structure containing spectra to be plotted
 *                              and associated info
 *         int nspec           number of spectra contained in spectrum
 *                              structure
 *         int maxspec         total number of spectra actually displayed.
 *         Spsetup *setup      container for plotting parameters
 *         int option          method of plotting spectra
 *
 * Output: int (0 or 1)         0 ==> success, 1 ==> error
 *
 */

int display_multispec(Spectrum *spectrum, int nspec, int maxspec, 
		      Spsetup *setup, int option)
{
  int i;              /* Looping variable */
  int dobox=2;        /* Flag defining plot axis labeling */
  float offset;       /* Offset from edges, in fractional terms */
  float vpmin,vpmax;  /* Min and max y coords for the viewport */
  float xmin;         /* Left-hand edge of viewport in fractional coords */
  float xmax;         /* Right-hand edge of viewport in fractional coords */
  float vpabsmin;     /* Absolute minimum viewport coordinate */
  float vpabsmax;     /* Absolute maximum viewport coordinate */
  float vpabsdiff;    /* Fractional height of viewport */
  float specmin;      /* Bottom edge of region for plotting spectra */
  float specdiff;     /* Fractional height of viewport occupied by spectra */

  /*
   * Define size of viewport containing all of the plots.
   *
   * ** NB: Make the viewport slightly wider than it is tall (in a fractional
   *    sense) so that macros like the AASTeX \plotone command, which scale
   *    figures such that the width fills the space alloted and then scale
   *    the height to preserve the aspect ratio, don't overrun the page.
   */

#if 0
  xmin=0.12;      
  xmax=0.88;      
  vpabsmin=0.12;  
  vpabsmax=0.88;  
#endif
  offset = setup->cheight / 10.0;
  xmin = offset;
  xmax = 1.0 - (offset / 2.0);
  vpabsmin = 0.75 * offset;
  vpabsmax = 1.0 - 2.0 * offset;
  specmin = vpabsmin + 0.5 * offset;
  vpabsdiff = vpabsmax - vpabsmin;
  specdiff = vpabsmax - specmin;

  /*
   * Loop through the spectra, defining different viewports for each
   *  one.  The first spectrum is plotted on top.
   */

  for(i=0; i<maxspec; i++) {

    /*
     * Define y min and y max of viewport and then set the viewport.
     */

    vpmax = vpabsmax - specdiff*i/nspec;
    vpmin = vpabsmax - specdiff*(i+1)/nspec;
    cpgsvp(xmin,xmax,vpmin,vpmax);

    /*
     * Plot the spectrum
     */
  
    if(i==(nspec-1))
      dobox = 4;
    else
      dobox = 2;
    if(display_spectrum(spectrum,nspec,i,setup,option,dobox)) {
      fprintf(stderr,"ERROR: display_multispec.\n");
      return 1;
    }
  }

  /*
   *  If all spectral have been plotted, add space at bottom for line
   *   labels.
   */

  if(i == nspec) {
    cpgsvp(xmin,xmax,vpabsmin,specmin);
    plot_specbox(setup->pltxmin,setup->pltxmax,0.0,1.0,6,setup->cheight);
  }

  /*
   * Define the viewport to cover the entire multiplot region.
   */

  cpgsvp(xmin,xmax,vpabsmin,vpabsmax);

  /*
   * Re-define the window to match that of the first spectrum, for
   *  compatability with ONESPEC and OVERPLOT options.
   */

  cpgswin(setup->pltxmin,setup->pltxmax,spectrum->minflux,spectrum->maxflux);

  return 0;
}

/*.......................................................................
 *
 * Function adjust_multispec
 *
 * Adjusts the display limits for multiple spectra, each plotted in its
 *  own window.  The main work is done through calls to get_vspec and
 *  display_multispec.
 *
 * Inputs: Spectrum *spectrum  Structure containing spectra to be plotted
 *                              and associated info
 *         int nspec           number of spectra contained in spectrum
 *                              structure
 *         Spsetup *setup      container for plotting parameters
 *         int option          method of plotting spectra
 *
 * Output: 0 or 1              0 ==> success, 1 ==> error
 *
 */

int adjust_multispec(Spectrum *spectrum, int nspec, Spsetup *setup, 
		     int option)
{
  int i;               /* Loping variable */
  int no_error=1;      /* Flag set to 0 on error */
  int do_adjust=1;     /* Flag set to 0 when no more adjustment is needed */
  char line[MAXC];     /* String variable for reading input */
  Spectrum *sp;        /* Pointer to navigate spectrum */

  /*
   * Loop through the spectra, displaying them and inquiring if any 
   *  adjustment is needed.
   */

  for(i=0; i<nspec; i++) {

    /*
     * Display all the spectra up to and including the one being adjusted.
     */

    cpgeras();
    if(display_multispec(spectrum,nspec,i+1,setup,option))
	no_error = 0;

    /*
     * Query if adjustment is desired.
     */

    while(do_adjust && no_error) {
      printf("\n--------------------------------------------------\n\n");
      printf("\nDo you want to change the vertical scale");
      printf(" of spectrum number %d? [y] ",i+1);
      fgets(line,MAXC,stdin);
      if(line[0] == 'n' || line[0] == 'N')
	do_adjust = 0;

      /*
       * If so, get the new values.
       */

      else {
	sp = spectrum + i;
	switch(get_vscale(&sp->minflux,&sp->maxflux)) {
	case -1:
	  no_error = 0;
	  break;
	case 1:
	  cpgeras();
	  if(display_multispec(spectrum,nspec,i+1,setup,option))
	    no_error = 0;
	  break;
	default:
	  do_adjust = 0;
	}
      }
    }

    /*
     * Reset do_adjust.
     */

    do_adjust = 1;
  }

  /*
   * Return
   */

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: adjust_multispec\n");
    return 1;
  }
}

