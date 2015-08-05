/*
 * tmp_er.c
 *
 * Usage: tmp_er <zl> <zs> <vdisp> [output_file_name]
 *
 *   where the output filename is optional.
 *
 * This program takes the redshifts of a lens and background source and
 *  some representation of the mass of the lens (e.g. its central velocity
 *  dispersion) and computes the Einstein ring radius.
 * NB: For now the only input allowed for the mass is the velocity dispersion
 *  for a singular isothermal sphere.
 *
 * A temporary version of calc_er.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cosmo.h"
#include "dataio.h"
#include "structdef.h"

#define DTH0 2.845397e-6

/*.......................................................................
 *
 * Function declarations
 *
 */

void sis_radius(double zl, double zs, double dzl, double dzs, double sigv, 
		double dsigv, double omega0, double lambda, FILE *ofp);

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int i;              /* Looping variable */
  int doflatloop=0;   /* Flag set to 1 to do a flat-cosmology loop */
  int contin=1;       /* Flag set to 0 to break loop */
  double zl;          /* Lens redshift */
  double zs;          /* Source redshift */
  double dzl=0.0;     /* Error on lens redshift */
  double dzs=0.0;     /* Error on source redshift */
  double omega0=1.0;  /* Omega_0 and default starting value */
  double lambda=0.0;  /* Omega_Lambda and default starting value */
  double sigv;        /* Velocity dipsersion in km/sec */
  double dsigv=0.0;   /* Error on velocity dispersion */
  char line[MAXC];    /* General string for reading input */
  char outname[MAXC]; /* Output filename, if desired */
  FILE *ofp=NULL;     /* Pointer for optional output file */

  /*
   * Open output file, if one is desired.
   */

  if(argc == 5) {
    strcpy(outname,argv[4]);
    printf("\n");
    if(!(ofp = open_appendfile(outname,1))) {
      fprintf(stderr,"\nERROR: Exiting program calc_er.c\n\n");
      return 1;
    }
  }

  /*
   * Get redshifts and velocity dispersion.
   */

  if(sscanf(argv[1],"%lf",&zl) != 1) {
    fprintf(stderr,"ERROR: exiting program.\n");
    return 1;
  }
  if(sscanf(argv[2],"%lf",&zs) != 1) {
    fprintf(stderr,"ERROR: exiting program.\n");
    return 1;
  }
  if(sscanf(argv[3],"%lf",&sigv) != 1) {
    fprintf(stderr,"ERROR: exiting program.\n");
    return 1;
  }

  /*
   * Do flat-cosmology loop
   */

  for(i=0; i<10; i++) {
    omega0 = 1.0 - i*0.1;
    lambda = 0.0 + i*0.1;
    sis_radius(zl,zs,dzl,dzs,sigv,dsigv,omega0,lambda,ofp);
  }
  omega0 = 0.01;
  lambda = 0.99;
  sis_radius(zl,zs,dzl,dzs,sigv,dsigv,omega0,lambda,ofp);
  printf("\n");

  /*
   * Clean up and exit
   */

  if(ofp)
    fclose(ofp);

  return 0;
}

/*.......................................................................
 *
 * Function sis_radius
 *
 * Calculates the Einstein ring for a lens system
 *  in a given cosmology.  If the redshift errors are non-zero,
 *  call calc_errs to calculate uncertainties.
 *
 * For now, this just works with a singular isothermal sphere, for which
 *  it is very straightforward to calculate the Einstein ring radius for
 *  a given velocity dispersion.  The appropriate formula is:
 *
 *  theta_E = (4 * PI * sigma^2 / c^2) * (D_ls / D_s)
 *
 * Inputs: double zl           lens redshift
 *         double zs           source redshift
 *         double dzl          error on zl
 *         double dzs          error on zs
 *         double sigv         velocity dispersion in km/sec
 *         double dsigv        uncertainty in sigv
 *         double omega0       Omega_M
 *         double lambda       Omega_Lambda
 *         FILE *ofp           optional output file pointer
 *
 * Output: none
 *
 */

void sis_radius(double zl, double zs, double dzl, double dzs, double sigv, 
		double dsigv, double omega0, double lambda, FILE *ofp)
{
  double dl;           /* Ang. diam. dist. to lens */
  double ds;           /* Ang. diam. dist. to source */
  double dls;          /* Ang. diam. dist. btwn source and lens */
  double d;            /* D_l D_s / D_ls */
  double me;           /* Mass enclosed in Einstein ring */
  double r_phys;       /* Physical radius of Einstein ring */
  double v_circ;       /* Implied circular velocity, given M_E */
  double sigbn;        /* Velocity dispersion from Blandford & Narayan eqn. */
  double theta_E;      /* Einstein ring radius */

  /*
   * Compute the angular diameter distances in Mpc
   */

  dl = ang_dist(0.0,zl,omega0,lambda) / MPC2CM;
  ds = ang_dist(0.0,zs,omega0,lambda) / MPC2CM;
  dls = ang_dist(zl,zs,omega0,lambda) / MPC2CM;

  /*
   * Compute theta_E.  The factor of 1 x 10^10 is included because
   *  c is defined (in cosmo.h) in cm/sec whereas we have defined
   *  sigv in km/sec.
   */

  theta_E = 4.0e10 * PI * sigv * sigv * dls / (C * C * ds);

  /*
   * Convert theta_E to arcsec
   */

  theta_E *= 206265.0;

  /*
   * Print out results
   */

  printf("\nFor this system\n");
  printf("   Omega_M = %6.3f\n",omega0);
  printf("   Omega_L = %6.3f\n",lambda);
  printf("   D_l     = %4.0f h^{-1} Mpc\n",dl);
  printf("   D_s     = %4.0f h^{-1} Mpc\n",ds);
  printf("   D_ls    = %4.0f h^{-1} Mpc\n",dls);
  printf("   theta_E = %8.5f arcsec\n",theta_E);

  /*
   * Print out results to file if file pointer has been set.
   */

  if(ofp) {
    fprintf(ofp,"%6.3f %6.3f %6.1f  %5.2f %5.2f  %4.0f %4.0f %4.0f %9.5f\n",
	    zl,zs,sigv,omega0,lambda,dl,ds,dls,theta_E);
  }

  /*
   * Compute errors on computed values, if redshift errors entered
   */
#if 0
  if(dzl > 0.0 || dzs > 0.0)
    calc_errs(zl,zs,dzl,dzs,theta,dtheta,omega0,lambda,dl,ds,dls,d,me);
#endif
}

