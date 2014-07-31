/*
 * sext2reg.c
 *
 * Usage: sext2reg [catfile] [informat] [outfile] [outformat] [dolabel]
 *
 * Description: This program reads in a SExtractor catalog and outputs
 *  a ds9 *.reg file to show the sources.  This reg file will show the
 *  Kron-like apertures that SExtractor shows as ellipses in its output
 *  "APERTURES" fits file.  However, being a text file rather than a fits
 *  file, the reg file is a more efficient way to keep the information.
 *  If the input file is just a list of positions, the reg file will show
 *  just circles rather than ellipses.    
 *
 * Inputs:
 *   catfile   -  Name of file containing SExtractor catalog
 *   informat  -  An integer representing the format of the
 *                 SExtractor catalog.
 *   outfile   -  Name of output region file
 *   outformat -  Either "radec", or "xy", representing the format
 *                 of the output file.
 *   dolabel   -  A flag set to 0 for no labels and 1 for labels
 *
 * Revision history:
 *  24Jul2006: First working version. Chris Fassnacht (CDF)
 *  21Nov2006: (CDF) More general revision that will work with more input 
 *              formats.
 *  30May2006: (CDF) Make the inclusion of labels an option accessed from the
 *              command line. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "catlib.h"


/*.......................................................................
 *
 * Function declarations
 */

void sext2reg_help();
int write_reg(Secat *secat, int ncat, int informat, char *outname, 
	      char *outformat, int dolabel);


/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;                   /* Looping variable */
  int no_error=1;          /* Flag set to 0 on error */
  int ncat;                /* Number of lines in the final catalog */
  int informat;            /* Format of input file */
  int dolabel;             /* Include labels? */
  char catfile[MAXC];      /* Filename for input catalog */
  char outfile[MAXC];      /* Filename for input catalog */
  char outformat[MAXC];    /* Format of the output file */
  Secat *secat=NULL;       /* Data array from catalog, after purging */
  Secat *sptr;             /* Pointer to navigate secat */

  /*
   * Check the command line invocation
   */

  if(argc < 6) {
    sext2reg_help();
    return 1;
  }
  printf("\n");

  /*
   * Get the inputs from the command line.
   */

  strcpy(catfile,argv[1]);
  if(sscanf(argv[2],"%d",&informat) != 1)
    no_error = 0;
  strcpy(outfile,argv[3]);
  strcpy(outformat,argv[4]);
  if(sscanf(argv[5],"%d",&dolabel) != 1)
    no_error = 0;

  /*
   * Fill the data structure
   */

  printf("informat = %d\n",informat);
  if(no_error)
    if(!(secat = read_secat(catfile,'#',&ncat,informat)))
      no_error = 0;

  /*
   * Write output file
   */

  if(no_error)
     if(write_reg(secat,ncat,informat,outfile,outformat,dolabel))
       no_error = 0;

  /*
   * Clean up and exit
   */

  secat = del_secat(secat);

  if(no_error) {
    printf("\nProgram sext2reg finished.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"\nERROR.  Exiting sext2reg.\n\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function write_reg
 *
 * Description: 
 *  Writes an output region file, given an input SExtractor file.
 *
 * Inputs:
 *  Secat *secat               input SExtractor catalog structure
 *  int ncat                   number of objects in catalog
 *  char *outname              output region filename
 *  char *outformat            output format (radec or xy)
 *  int informat               input file format
 *
 * Output:
 *  int (0 or 1)               0 ==> success, 1 ==> error
 */

int write_reg(Secat *secat, int ncat, int informat, char *outname, 
	      char *outformat, int dolabel)
{
  int i;               /* Looping variable */
  int no_error=1;      /* Flag set to 0 on error */
  float rfac;          /* Factor by which to multiply a and b */ 
  double textx,texty;  /* Location for text in label */
  Secat *sptr;         /* Pointer to navigate secat */
  FILE *ofp=NULL;      /* Output file pointer */


  /*
   * Open output file
   */

  if(!(ofp = open_writefile(outname))) {
    fprintf(stderr,"ERROR: write_reg.\n");
    return 1;
  }

  /*
   * Print header
   */

  fprintf(ofp,"global color=green\n");

  /*
   * Print ellipse positions, either in (x,y) coordinates 
   *  or in (ra,dec) format 
   */

  if(strcmp(outformat,"radec") == 0) {
    for(i=0,sptr=secat; i<ncat; i++,sptr++) {
      switch(informat) {
#if 0
      case 3: case 6: case 7:
	rfac = sptr->r_kron / 3600.0 / 20.0;
	fprintf(ofp,"fk5;ellipse(%f,%f,%f,%f,%.0f)\n",
		sptr->alpha,sptr->delta,rfac * sptr->a_im,rfac * sptr->b_im,
		sptr->theta);
	break;
#endif
      default:
	fprintf(ofp,"fk5;circle(%f,%f,0.0007)\n",sptr->alpha,sptr->delta);
	break;
      }
      if(dolabel) {
	textx = sptr->alpha - 0.0007;
	texty = sptr->delta + 0.0007;
	switch(informat) {
	case 1: case 2: case 9: case 10:
	  break;
	default:
	  sprintf(sptr->name,"%d",sptr->id);
	  break;
	}
	fprintf(ofp,"fk5;text(%f,%f) # text={%s}\n",textx,texty,sptr->name);
      }
    }
  }
  else {
    for(i=0,sptr=secat; i<ncat; i++,sptr++) {
      switch(informat) {
      case 6: case 7: case 18:
	rfac = sptr->r_kron;
	fprintf(ofp,"image;ellipse(%.0f,%.0f,%.0f,%.0f,%.0f)\n",
		sptr->x, sptr->y, rfac * sptr->a_im, rfac * sptr->b_im,
		sptr->theta);
	break;
      default:
	fprintf(ofp,"image;circle(%.0f,%.0f,20)\n",
		sptr->x, sptr->y);
	break;
      }
      if(dolabel) {
	textx = sptr->x + 5.0;
	texty = sptr->y + 5.0;
	switch(informat) {
	case 1: case 2: case 9: case 10:
	  break;
	default:
	  sprintf(sptr->name,"%d",sptr->id);
	  break;
	}
	fprintf(ofp,"image;text(%.1f,%.1f) # text={%s}\n",textx,texty,sptr->name);
      }
    }
  }

  /*
   * Exit
   */

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: write_reg\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function sext2reg_help
 *
 * Description:
 *  Describes program and inputs
 *  
 * Inputs:
 *  (none)
 *
 * Output:
 *  (none)
 */

void sext2reg_help() {
  char line[MAXC];   /* String for getting input */

  fprintf(stderr,"\nUsage: sext2reg [catfile] [informat] [outfile] ");
  fprintf(stderr,"[outformat] [dolabel]\n\n");
  fprintf(stderr,"Inputs:\n");
  fprintf(stderr,"  catfile   - file containing the SExtractor catalog\n");
  fprintf(stderr,"  informat  - indicates the format of the input file");
  fprintf(stderr," (see below)\n");
  fprintf(stderr,"  outfile   - output ds9 region file\n");
  fprintf(stderr,"  outformat - either \"radec\", or \"xy\"\n");
  fprintf(stderr,"              (sets the format of the output region file)\n");
  fprintf(stderr,"  dolabel   - include labels? either 1 (yes) or 0 (no)\n\n");
  fprintf(stderr,"Hit return for possible input formats:  ");
  fgets(line,MAXC,stdin);
  secat_format();
}
