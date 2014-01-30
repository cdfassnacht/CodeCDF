/*
 * cosmocalc.c
 *
 * This program computes various cosmological quantities based on an input
 *  redshift and world model.  H_0 is assumed to be unknown.
 *
 * First working version: 07Nov1999 by Chris Fassnacht (CDF)
 * Revision history:
 *  v28Jan2002 CDF,  Added distance modulus.
 *  v11Jul2006 CDF,  Modified to use new get_cosmo and calc_cosdist functions
 *                    in cosmo.c
 *  v02Feb2007 CDF,  Incorporated batch mode which formerly was in 
 *                    cosmo_multiz.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "cosmo.h"

/*.......................................................................
 *
 * Function declarations
 *
 */

int cosmocalc_interactive();
int cosmocalc_batch(char *inname, char *outname);
void cosmocalc_help();

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int no_error=1;  /* Flag set to 0 on error */

  /*
   * Check input line format
   */

  switch(argc) {
  case 1:
    cosmocalc_help();
    return 1;
    break;
  case 2:
    if(strcmp(argv[1],"-i") == 0) {
      /* Interactive input */
      if(cosmocalc_interactive())
	no_error = 0;
    }
    else {
      cosmocalc_help();
    }
    break;
  case 3:
    if(strcmp(argv[1],"-b") == 0) {
      fprintf(stderr,"\nERROR: batch mode requires both an input and an ");
      fprintf(stderr,"output filename.\n\n");
      no_error = 0;
    }
    else if(strcmp(argv[1],"-i") == 0) {
      if(cosmocalc_interactive())
	no_error = 0;
    }
    else {
      cosmocalc_help();
    }
    break;
  case 4:
    /* File-based input with a designated output filename */
    if(strcmp(argv[1],"-b") == 0) {
      if(cosmocalc_batch(argv[2],argv[3]))
	no_error = 0;
    }
    else if(strcmp(argv[1],"-i") == 0) {
      if(cosmocalc_interactive())
	no_error = 0;
    }
    else {
      cosmocalc_help();
    }
    break;
  default:
    fprintf(stderr,"\n***WARNING: Too many arguments.***\n\n");
    cosmocalc_help();
    return 1;
    break;
  }

  return 0;
}

/*.......................................................................
 *
 * Function cosmocalc_interactive
 *
 * Takes redshift and cosmology entered by the user and prints output.
 *
 * Inputs: none
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 */



int cosmocalc_interactive()
{
  int contin=1;       /* Flag set to 0 to break loop */
  double z;           /* Lens redshift */
  double dz=0.0;      /* Error on lens redshift */
  double tmpdl;       /* Temporary storage in D_L^2 calculation */
  Cosmo cosmo;        /* Cosmological world model */
  Cosdist cosdist;    /* Structure for all of the distance measures */
  char line[MAXC];     /* General string for reading input */

  /*
   * Get redshift
   */

  if(get_valerr(&z,&dz,"source redshift") == 0) {
    fprintf(stderr,"\nERROR.  Exiting cosmocalc.\n\n");
    return 1;
  }

  /*
   * Initialize world model
   */

  cosmo.omega_m = 0.3;
  cosmo.omega_de = 0.7;
  cosmo.w = -1.0;

  /*
   * Loop through world models, if desired
   */

  while(contin) {

    /*
     * Get the cosmological world model
     */

    get_cosmo(&cosmo);

    /*
     * Calculate the distance measures
     */

    cosdist = calc_cosdist(0,z,cosmo);
    tmpdl = 4.0 * PI * cosdist.d_l * cosdist.d_l;

    /*
     * Print out quantities of interest
     */

    printf("\n--------------------------------------------------\n\n");
    printf("Inputs:\n");
    printf("  z             = %f\n",z);
    printf("  Omega_m       = %f\n",cosmo.omega_m);
    printf("  Omega_DE      = %f\n",cosmo.omega_de);
    if(cosmo.omega_de != 0.0)
      printf("  w             = %f\n",cosmo.w);
    printf("\nCalculated quantities:\n");
    printf("   Distance modulus          = %5.2f - 5 log h mag\n",cosdist.DM);
    printf("   Angular diameter distance = %4.0f h^{-1} Mpc\n",
	   cosdist.d_a/MPC2CM);
    printf("   Luminosity distance       = %4.0f h^{-1} Mpc\n",
	   cosdist.d_l/MPC2CM);
    printf("   4*PI*D_L^2                =  %8.3e h^{-2} cm^2\n",
	   tmpdl);
    printf("   1 arcsec                  = %6.3f h^{-1} kpc\n",
	   cosdist.d_a * 1000.0 /(MPC2CM * RAD2ASEC));
    printf("   1.0 h^{-1} Mpc            = %5.2f arcmin\n",
	   RAD2ASEC * MPC2CM /(cosdist.d_a * 60.0));
    printf("   1.0 h^{-1} _comoving_ Mpc = %5.2f arcmin\n",
	   RAD2ASEC * MPC2CM/(cosdist.d_m * 60.0));
    printf("   Lookback time             =  %8.3e h^{-1} yr\n",
	   cosdist.t_l/YR2SEC);
    printf("   H(z)                      = %6.1f h km/s/Mpc\n",
	   cosdist.hz);
    printf("\n--------------------------------------------------\n\n");

    /*
     * Continue?
     */

    printf("Another cosmological world model? (y/n) [n] ");
    fgets(line,MAXC,stdin);
    if(line[0] == 'Y' || line[0] == 'y')
      contin = 1;
    else
      contin = 0;
  }

  printf("\n");
  return 0;
}

/*.......................................................................
 *
 * Function cosmocalc_batch
 *
 * Takes redshifts and cosmology included in a file and prints output.
 *
 * Inputs: 
 *  char *inname        input filename
 *  char *outname       output filename
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 */

int cosmocalc_batch(char *inname, char *outname)
{
  int i;              /* Looping variable */
  int contin=1;       /* Flag set to 0 to break loop */
  int no_error=1;     /* Flag set to 0 on error */
  int nz=0;           /* Number of redshifts */
  double *z=NULL;     /* Lens redshift */
  double *zptr;       /* Pointer to navigate z */
  double dz=0.0;      /* Error on lens redshift */
  double tmpdl;       /* Temporary storage in D_L^2 calculation */
  char comment='#';   /* Comment character in input file */
  char line[MAXC];     /* General string for reading input */
  Cosmo cosmo;        /* Cosmological world model */
  Cosdist cosdist;    /* Structure for all of the distance measures */
  FILE *ifp=NULL;     /* Input file pointer */
  FILE *ofp=NULL;     /* Output file pointer */

  /*
   * Open files
   */

  if(!(ifp = open_readfile(inname)))
    no_error = 0;
  if(!(ofp = open_writefile(outname)))
    no_error = 0;

  /*
   * Get number of lines in input file
   */

  if(no_error) {
    if((nz = n_lines(ifp,'#')) == 0) {
      fprintf(stderr,"ERROR: No valid data in %s.\n",inname);
      no_error = 0;
    }
    else {
      rewind(ifp);
    }
  }

  /*
   * Allocate memory for redshift array.
   */

  if(no_error) {
    if(!(z = new_doubarray(nz)))
      no_error = 0;
    else
      zptr = z;
  }

  /*
   * Read in redshifts
   */

  while(no_error && fgets(line,MAXC,ifp) != NULL) {
    if(line[0] != comment) {
      if(sscanf(line,"%lf",zptr) != 1) {
	fprintf(stderr,"ERROR: Invalid data format in %s\n",inname);
	no_error = 0;
      }
      else
	zptr++;
    }
  }

  /*
   * Print out header info
   */

  if(no_error) {
    fprintf(ofp,"# Column info\n");
    fprintf(ofp,"#------------\n");
    fprintf(ofp,"#   1. Redshift\n");
    fprintf(ofp,"#   2. Omega_m\n");
    fprintf(ofp,"#   3. Omega_DE\n");
    fprintf(ofp,"#   4. w\n");
    fprintf(ofp,"#   5. H(z) in h km/s/Mpc\n");
    fprintf(ofp,"#   6. t_l (lookback time), in h^{-1} years\n");
    fprintf(ofp,"#   7. DM (distance modulus), in mags - 5 log h\n");
    fprintf(ofp,"#   8. D_A (ang. diam. distance, in h^{-1} Mpc\n");
    fprintf(ofp,"#   9. D_L (luminosity distance, in h^{-1} Mpc\n");
    fprintf(ofp,"#  10. 4*PI*D_L^2, in h^{-2} cm^2\n");
    fprintf(ofp,"#  11. Physical distance subtending 1 arcsec, ");
    fprintf(ofp,"in h^{-1} kpc\n");
    fprintf(ofp,"#  12. Angle subtended by 1 h^{-1} Mpc, in arcmin\n");
    fprintf(ofp,"#  13. Angle subtended by 1 h^{-1} COMOVING Mpc, in arcmin\n");
    fprintf(ofp,"#\n#\n");
    fprintf(ofp,"#  z   Om_m Om_L  w   H(z)    t_l     DM   D_A  D_L  ");
    fprintf(ofp,"4*PI*DL^2 1asec  1Mpc  1Mpc\n");
    fprintf(ofp,"#----- ---- ---- ---- ---- --------- ----- ---- ---- ");
    fprintf(ofp,"--------- ----- ----- -----\n");
  }

  /*
   * Initialize world model
   */

  cosmo.omega_m = 0.3;
  cosmo.omega_de = 0.7;
  cosmo.w = -1.0;

  /*
   * Loop through world models, if desired
   */

  while(no_error && contin) {

    /*
     * Get the cosmological world model
     */

    get_cosmo(&cosmo);

    /*
     * Loop through all the redshifts
     */

    for(i=0,zptr=z; i<nz; i++,zptr++) {

      /*
       * Calculate the distance measures
       */

      cosdist = calc_cosdist(0,*zptr,cosmo);
      tmpdl = 4.0 * PI * cosdist.d_l * cosdist.d_l;

      /*
       * Print out quantities of interest to the output file
       */

      fprintf(ofp,"%6.3f %4.2f %4.2f %4.1f %4.0f %9.3e %5.2f %4.0f %4.0f ",
	      *zptr,cosmo.omega_m,cosmo.omega_de,cosmo.w,
	      cosdist.hz,cosdist.t_l/YR2SEC,
	      cosdist.DM,cosdist.d_a/MPC2CM,cosdist.d_l/MPC2CM);
      fprintf(ofp,"%9.3e %5.2f %5.2f %5.2f\n",
	      tmpdl,
	      cosdist.d_a*1000.0/(MPC2CM*RAD2ASEC),
	      RAD2ASEC*MPC2CM/(cosdist.d_a*60.0),
	      RAD2ASEC*MPC2CM/(cosdist.d_m*60.0));
    }

    /*
     * Continue?
     */

    printf("\nCalculated parameters for (%4.2f,%4.2f,%4.2f) and wrote to %s\n",
	   cosmo.omega_m,cosmo.omega_de,cosmo.w,outname);
    printf("Another cosmological world model? (y/n) [n] ");
    fgets(line,MAXC,stdin);
    if(line[0] == 'Y' || line[0] == 'y')
      contin = 1;
    else
      contin = 0;
  }

  /*
   * Clean up and exit
   */

  z = del_doubarray(z);
  if(ifp)
    fclose(ifp);
  if(ofp)
    fclose(ofp);

  printf("\n");
  if(no_error) {
    printf("Finished with cosmo_multiz.c.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"ERROR: cosmo_multiz.c\n\n");
    return 1;
  }

  return 0;
}

/*.......................................................................
 *
 * Function cosmocalc_help
 *
 * Prints out helpful information when the program is called with no 
 *  arguments.
 *
 * Inputs: none
 *
 * Output: none
 *
 */

void cosmocalc_help()
{
  printf("\n*******************************************************\n\n");
  printf("\nProgram: cosmocalc -- ");
  printf("Calculates cosmological parameters given redshifts\n\n");
  printf("Usage: cosmocalc -i\n");
  printf("       cosmocalc -b input_file output_file\n\n");
  printf("OPTIONS:\n");
  printf(" -i   Interactive mode.\n\n");
  printf(" -b   Convert positions in a batch mode from an input file.\n");
  printf("      For this option, each line of the input file contains a\n");
  printf("      redshift\n\n");
  printf("*******************************************************\n\n");
}
