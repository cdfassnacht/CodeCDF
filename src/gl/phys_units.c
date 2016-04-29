/*
 * phys_units.c
 *
 * This program takes the redshifts of a lens and background source and
 *  computes various quantities including, angular diameter distances,
 *  and masses inside the Einstein ring of the lens.
 *
 * First working version: 05 Jul 1997 by Chris Fassnacht (CDF)
 *
 * Revision history:
 *  16Oct1997 CDF, Calculate errors on angular diameter distances given the
 *                  errors on the redshifts.
 *  18Oct1997 CDF, Better I/O
 *  05Aug1998 CDF, Added the option to have non-zero Lambda cosmologies.
 *  08Aug1998 CDF, Added calculations of velocity dispersion and circular
 *                  velocities.
 *  09Aug1998 CDF, Correct circular velocities for effect of finding
 *                  total mass inside cylinder of radius r_E rather than
 *                  in a sphere of radius r_E, which is what you want.
 *                 Moved input of cosmological parameters into get_omega in
 *                  cosmo.c.
 *  29Sep1998 CDF, Moved calculations into calc_params and calc_errs.
 *                 Added loop over cosmological world models.
 *  12Jul2006 CDF, Cleaned up input and transferred many of the calculations
 *                  to the new calc_cosdist function in cosmo.c
 *                 **** NB: calc_errs has not been fixed yet ****
 *  28Apr2016 CDF, Changed the output (when printing to an output file)
 *                  to handle the case when one of the input redshifts is
 *                  unknown (indicated by a value less than 0).
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "cosmo.h"
#include "dataio.h"

#define DTH0 2.845397e-6
#define PI 3.141592653589793

/*.......................................................................
 *
 * Function declarations
 *
 */

void phys_units_help();
int calc_interactive();
int calc_batch(char *inname, char *outname);
void calc_params(double zl, double zs, double dzl, double dzs, double theta, 
		 double dtheta, Cosmo cosmo, FILE *ofp);
#if 0
void calc_errs(double zl, double zs, double dzl, double dzs, double theta, 
	       double dtheta, double omega0, double lambda, double dl,
	       double ds, double dls, double d, double me);
#endif

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int no_error = 1;     /* Flag set to 0 on error */

  /*
   * Check input line format
   */

  switch(argc) {
  case 1:
    phys_units_help();
    return 1;
    break;
  case 2:      /* Interactive input */
    if(strcmp(argv[1],"-i") == 0) {
      if(calc_interactive())
	no_error = 0;
    }
    else {
      phys_units_help();
    }
    break;
  case 3:   /* File-based input with NO designated output filename */
    if(strcmp(argv[1],"-b") == 0) {
      if(calc_batch(argv[2],NULL))
	no_error = 0;
    }
    else if(strcmp(argv[1],"-i") == 0) {
      if(calc_interactive())
	no_error = 0;
    }
    else {
      phys_units_help();
    }
    break;
  case 4:   /* File-based input with a designated output filename */
    if(strcmp(argv[1],"-b") == 0) {
      if(calc_batch(argv[2],argv[3]))
	no_error = 0;
    }
    else if(strcmp(argv[1],"-i") == 0) {
      if(calc_interactive())
	no_error = 0;
    }
    else {
      phys_units_help();
    }
    break;
  default:
    fprintf(stderr,"\n***WARNING: Too many arguments.***\n\n");
    phys_units_help();
    return 1;
    break;
  }

  /*
   * Clean up and exit
   */


  if(no_error) {
    printf("\nCompleted program phys_units.\n\n");
    return 0;
  }
  else {
    fprintf(stderr,"ERROR.  Exiting program\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function phys_units_help
 *
 * Prints out helpful information when the program is called with no 
 *  arguments.
 *
 * Inputs: none
 *
 * Output: none
 *
 */

void phys_units_help()
{
  printf("\n**********************************************************\n\n");
  printf("\nProgram: phys_units -- \n");
  printf(" Calculates gravitational lens parameters for given cosmologies.\n\n");
  printf("Usage: phys_units -i\n");
  printf("       phys_units -b input_file [output_file]\n");
  printf("Note that the output filename is optional for the -b option.\n\n");
  printf("OPTIONS:\n");
  printf(" -i   Interactive mode.\n\n");
  printf(" -b   Batch mode:\n");
  printf("      Each line of the input file contains redshift and \n");
  printf("      cosmological information.\n");
  printf("      Input file format: z_l z_s theta omega_m omega_de w\n");
  printf("          z_l      = lens redshift\n");
  printf("          z_s      = source redshift\n");
  printf("          theta    = image separation\n");
  printf("          omega_m  = Omega_matter\n");
  printf("          omega_de = Omega_dark-energy (Omega_Lambda)\n");
  printf("          w        = Dark energy equation of state parameter\n");
  printf("                     (Set to -1 for a cosmological constant)\n");
  printf("      Use a # at the beginning of a line to comment it out.\n\n");
  printf("***********************************************************\n\n");
}

/*.......................................................................
 *
 * Function calc_interactive
 *
 * Interactive call(s) to calc_params
 *
 * Inputs: none
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 */

int calc_interactive()
{
  int contin=1;        /* Flag set to 0 to break loop */
  double zl;           /* Lens redshift */
  double zs;           /* Source redshift */
  double dzl=0.0;      /* Error on lens redshift */
  double dzs=0.0;      /* Error on source redshift */
  double theta;        /* Angular separation in arcsec */
  double dtheta=0.0;   /* Error in angular separation in arcsec */
  Cosmo cosmo;         /* Cosmological world model */
  char line[MAXC];     /* General string for reading input */

  /*
   * Initialize world model
   */

  cosmo.omega_m = 0.3;
  cosmo.omega_de = 0.7;
  cosmo.w = -1.0;

  /*
   * Get redshifts and image separation
   */

  if(get_valerr(&zl,&dzl,"lens redshift") == 0) {
    fprintf(stderr,"\nERROR.  Exiting phys_units.\n\n");
    return 1;
  }

  if(get_valerr(&zs,&dzs,"source redshift") == 0) {
    fprintf(stderr,"\nERROR.  Exiting phys_units.\n\n");
    return 1;
  }

  if(get_valerr(&theta,&dtheta,"image separation (arcsec)") == 0) {
    fprintf(stderr,"\nERROR.  Exiting phys_units.\n\n");
    return 1;
  }

  /*
   * Loop through world models, if desired
   */

  while(contin) {

    /*
     * Get the cosmological world model
     */

    printf("\n");
    get_cosmo(&cosmo);

    /*
     * Calculate parameters
     */

    calc_params(zl,zs,dzl,dzs,theta,dtheta,cosmo,NULL);

    /*
     * Continue?
     */

    printf("Continue? (y/n) [y] ");
    fgets(line,MAXC,stdin);
    if(strcmp(line,"") != 0) {
      if(line[0] == 'N' || line[0] == 'n')
	contin = 0;
    }
  }

  return 0;
}

/*.......................................................................
 *
 * Function calc_batch
 *
 * Reads in lines from a file as inputs to calc_params.
 *
 * Inputs: char *inname        input filename
 *         char *outname       output filename (NULL if none designated on
 *                              command line)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 */

int calc_batch(char *inname, char *outname)
{
  int no_error = 1;     /* Flag set to 0 on error */
  double zl;           /* Lens redshift */
  double zs;           /* Source redshift */
  double dzl=0.0;      /* Error on lens redshift */
  double dzs=0.0;      /* Error on source redshift */
  double theta;        /* Angular separation in arcsec */
  double dtheta=0.0;    /* Error in angular separation in arcsec */
  Cosmo cosmo;        /* Cosmological world model */
  char newoutn[MAXC];    /* New output name on error */
  char line[MAXC];       /* General string for reading input */
  FILE *ifp=NULL;       /* Input file pointer */
  FILE *ofp=NULL;       /* Output file pointer */
 
  /*
   * Open the input file
   */

  if(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: calc_batch.\n");
    return 1;
  }

  /*
   * Open the output file if one is desired
   */

  if(outname == NULL) {
    printf("\nEnter name of output file if one is desired or\n");
    printf("just hit return for no output file: [no output file] ");
    fgets(line,MAXC,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%s",newoutn) != 1) {
	fprintf(stderr,"ERROR: Bad input.  Enter filename again:  ");
	fgets(line,MAXC,stdin);
      }
      if(!(ofp = open_writefile(newoutn))) {
	fprintf(stderr,"ERROR: calc_batch\n");
	return 1;
      }
    }
  }
  else {
    if(!(ofp = open_writefile(outname))) {
      fprintf(stderr,"ERROR: calc_batch\n");
      return 1;
    }
    fprintf(ofp,"# z_l    z_s  theta  O_m   O_de   w    D_l  D_s D_ls   D\n");
  }

  /*
   * Read in the list of lens system redshifts and cosmologies
   */

  while(fgets(line, MAXC, ifp) != NULL && no_error) {
    if(line[0] != '#') {
      if(sscanf(line,"%lf %lf %lf %lf %lf %lf",
		&zl,&zs,&theta,&cosmo.omega_m,&cosmo.omega_de,&cosmo.w) != 6) {
	  fprintf(stderr,"ERROR.  Bad input file format.\n");
	  fprintf(stderr,"This program requires the following format:\n");
	  fprintf(stderr,"  z_l z_s theta omega_m omega_de w.\n");
	  no_error = 0;
      }
      else {
	calc_params(zl,zs,dzl,dzs,theta,dtheta,cosmo,ofp);
      }
    }
  }


  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);
  if(ofp)
    fclose(ofp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: calc_batch\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function calc_params
 *
 * Calculates the physical parameters (angular diameter distances,
 *  Einstein masses, and velocity dispersions) for a lens system
 *  in a given cosmology.  If the redshift errors are non-zero,
 *  call calc_errs to calculate uncertainties.
 *
 * Inputs: double zl            lens redshift
 *         double zs            source redshift
 *         double dzl           error on zl
 *         double dzs           error on zs
 *         double theta         image separation in arcsec
 *         double dtheta        uncertainty in theta
 *         Cosmo cosmo          cosmological world model
 *         FILE *ofp            output file pointer
 *
 * Output: none
 *
 */

void calc_params(double zl, double zs, double dzl, double dzs, double theta, 
		 double dtheta, Cosmo cosmo, FILE *ofp)
{
  double d;            /* D = D_l D_s / D_ls */
  double theta2;       /* 0.5 * theta (estimate of Einstein ring radius) */
  double me;           /* Mass enclosed in Einstein ring */
  double r_phys;       /* Physical radius of Einstein ring */
  double v_circ;       /* Implied circular velocity, given M_E */
  double sigbn;          /* Velocity dispersion from Blandford & Narayan eqn. */
  double sigma_crit;     /* Sigma_critical */
  Cosdist cdl,cds,cdls;  /* Distance measures */

  /*
   * Calculate the distances with calc_cosdist
   */

  cdl = calc_cosdist(0.0,zl,cosmo);
  cds = calc_cosdist(0.0,zs,cosmo);
  cdls = calc_cosdist(zl,zs,cosmo);

  /*
   * Compute D and M_E -- put D in Gpc
   */

  d = cdl.d_a * cds.d_a / (cdls.d_a * 1000.0 * MPC2CM);
  theta2 = theta / 2.0;
  me = C*C*MPC2CM*1000.0*d*theta2*theta2 / (4.0*G*RAD2ASEC*RAD2ASEC*MSUN);

  /*
   * Compute Sigma_crit
   */

  sigma_crit = C*C*cds.d_a / (4.0 * PI * G * cdl.d_a * cdls.d_a);

  /*
   * Compute v_circ and sigma.
   *  NB: The factor of sqrt (2/PI) in the circular velocity calculation
   *       comes from the difference between the total mass enclosed in
   *       the cylinder of radius r_phys (what we have) and the mass inside
   *       the sphere of radius r_phys (what we want).  See Shu, "The
   *       Physical Universe", pp. 299-301.
   */

  r_phys = theta2 * cdl.d_a / RAD2ASEC;
  v_circ = sqrt(2 * G * me * MSUN / (PI * r_phys))/KM2CM;
  sigbn = 300.0 * sqrt(theta2 * cds.d_a / (2.6 * cdls.d_a));

  /*
   * Print out results
   */

  if(ofp) {
    fprintf(ofp,"%6.3f %6.3f %5.2f %5.2f %5.2f %5.2f ",
	    zl,zs,theta,cosmo.omega_m,cosmo.omega_de,cosmo.w);
    if(cdl.d_a < 0.) {
      fprintf(ofp," -99 ");
      d = -99.;
    }
    else
      fprintf(ofp,"%4.0f ",cdl.d_a/MPC2CM);

    if(cds.d_a < 0.) {
      fprintf(ofp," -99 ");
      d = -99.;
    }
    else
      fprintf(ofp,"%4.0f ",cds.d_a/MPC2CM);

    if(cdls.d_a < 0.)
      fprintf(ofp," -99 -99\n");
    else
      fprintf(ofp,"%4.0f %7.4f\n",cdls.d_a/MPC2CM,d);
  }
  else {
    fprintf(stdout,
	    "\nFor this system with (Omega_M = %5.3f, Omega_Lambda = %5.3f)\n",
	    cosmo.omega_m,cosmo.omega_de);
    fprintf(stdout,"   D_l     = %4.0f h^{-1} Mpc\n",cdl.d_a/MPC2CM);
    fprintf(stdout,"   D_s     = %4.0f h^{-1} Mpc\n",cds.d_a/MPC2CM);
    fprintf(stdout,"   D_ls    = %4.0f h^{-1} Mpc\n",cdls.d_a/MPC2CM);
    fprintf(stdout,"   D       = %7.4f h^{-1} Gpc\n",d);
    fprintf(stdout,"   Sigma_c = %5.2f h g/cm^2\n",sigma_crit);
    fprintf(stdout,"   M_E     = %9.3e h^{-1} M_sun\n",me);
    fprintf(stdout,"   R_E     = %8.3f h^{-1} kpc\n",r_phys*1000.0/MPC2CM);
    fprintf(stdout,"   v_circ  = %5.0f km/sec\n",v_circ);
    fprintf(stdout,"   sigma   = %5.0f km/sec\n\n",sigbn);
  }

  /*
   * Compute errors on computed values, if redshift errors entered
   */
#if 0
  if(dzl > 0.0 || dzs > 0.0)
    calc_errs(zl,zs,dzl,dzs,theta,dtheta,cosmo,cdl,cds,cdls,me);
#endif

}

#if 0
/*.......................................................................
 *
 * Function calc_errs
 *
 * Calculates the uncertainties in physical parameters derived in
 *  calc_params.
 *
 * Inputs: double zl            lens redshift
 *         double zs            source redshift
 *         double dzl           error on zl
 *         double dzs           error on zs
 *         double theta         image separation in arcsec
 *         double dtheta        uncertainty in theta
 *         Cosmo cosmo          cosmological world model
 *         Cosdist cdl          container for cosmological distances to zl
 *         Cosdist cds          container for cosmological distances to zs
 *         Cosdist cdls         container for cosmological distances between 
 *                               zl and zs
 *         double me            Einstein ring mass
 *
 * Output: none
 *
 */

void calc_errs(double zl, double zs, double dzl, double dzs, double theta, 
	       double dtheta, Cosmo cosmo, Cosdist cdl, Cosdist cds, 
	       Cosdist cdls,  double me)
{
  double ddl,dds,ddls; /* Errors on D_l, D_s and D_ls */
  double dd;           /* Error on D */
  double dme;          /* Error on M_E */

  /*
   * Error on D -- answer in Gpc
   */

  d = cdl.d_a * cds.d_a / (cdls.d_a * 1000.0 * MPC2CM);
  dd = d * sqrt(ddl*ddl/(dl*dl) + dds*dds/(ds*ds) + ddls*ddls/(dls*dls));

  /*
   * Error on M_E, answer in solar masses
   */

  dme = me * sqrt(dd*dd/(d*d) + 4.0*dtheta*dtheta/(theta*theta));

  /*
   * Print out answers
   */

  printf("\nUncertainties on the above quantities are:\n");
  printf("   dD_l  = %7.4f h^{-1} Mpc\n",ddl);
  printf("   dD_s  = %7.4f h^{-1} Mpc\n",dds);
  printf("   dD_ls = %7.4f h^{-1} Mpc\n",ddls);
  printf("   dD    = %7.4f h^{-1} Gpc\n",dd);
  printf("   dM_E  = %9.3e h^{-1} M_sun\n\n",dme);

}
#endif
