/*
 * sis_vdisp.c
 *
 * This program computes quantities associated with a singular isothermal
 *  sphere (e.g., Einstein radius, shear) given inputs.  The mass of the
 *  SIS that is given to the program will be in terms of the 1D velocity
 *  dispersion.
 *
 * 18Sep01 CDF,  A modification of cosmocalc.c and phys_units.c.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "cosmo.h"

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int contin=1;        /* Flag set to 0 to break loop */
  double offset;       /* Offset of the SIS center from the lens, in arcsec */
  double sigv;         /* Velocity dispersion of the SIS in km/sec */
  double dsigv=0.0;    /* Error on velocity dispersion */
  double v2;           /* Square of velocity dispersion in cm^2/sec^2 */
  double zl,zs;        /* SIS and source redshift */
  double dzl,dzs;      /* Error on redshifts */
  double omega0=1.0;   /* Omega_0 (default starting value is EdS) */
  double lambda=0.0;   /* Omega_Lambda (default starting value is EdS) */
  double dal,das,dals; /* Angular diameter distances */
  double theta_E;      /* Einstein ring radius, in arcsec */
  double gamma;        /* Shear */
  char line[MAXC];     /* General string for reading input */

  /*
   * Initialize
   */

  dzl = dzs = 0.0;

  /*
   * Get SIS velocity dispersion and calculate its square (converting
   *  to cm/sec from km/sec.
   */

  if(get_valerr(&sigv,&dsigv,"SIS velocity dispersion (km/sec)") == 0) {
    fprintf(stderr,"\nERROR.  Exiting calc_er.c\n\n");
    return 1;
  }
  v2 = 1.0e10 * sigv * sigv;

  /*
   * Get redshifts
   */

  if(get_valerr(&zl,&dzl,"redshift of SIS") == 0) {
    fprintf(stderr,"\nERROR.  Exiting sis_sigv.\n\n");
    return 1;
  }
  if(get_valerr(&zs,&dzs,"redshift of background source") == 0) {
    fprintf(stderr,"\nERROR.  Exiting sis_sigv.\n\n");
    return 1;
  }


  /*
   * Get SIS velocity dispersion
   */

  printf("\nEnter offset of SIS center from lens system (arcsec): ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%lf",&offset) != 1 || offset < 0.0) {
    fprintf(stderr,"ERROR.  Bad input.  Enter offst again.  ");
    fgets(line,MAXC,stdin);
  }

  /*
   * Loop through world models, if desired
   */

  while(contin) {

    /*
     * Get the cosmological parameters
     */

    printf("\n");
    get_omega(&omega0,&lambda);
    printf("  Omega_m       = %f\n",omega0);
    printf("  Omega_Lambda  = %f\n\n",lambda);
 
    /*
     * Compute the angular diameter distances in Mpc
     */

    dal = ang_dist(0.0,zl,omega0,lambda) / MPC2CM;
    das = ang_dist(0.0,zs,omega0,lambda) / MPC2CM;
    dals = ang_dist(zl,zs,omega0,lambda) / MPC2CM;

    /*
     * Compute the Einstein ring radius, in arcsec
     */

    theta_E = RAD2ASEC * 4 * PI * v2 * dals / (C * C * das);

    /*
     * Compute the shear if offset is larger than theta_E
     */

    if(offset > theta_E)
      gamma = theta_E / (2.0 * offset);

    /*
     * Print out quantities of interest
     */

    printf("\n--------------------------------------------------\n\n");
    printf("Inputs:\n");
    printf("  z_l           = %f\n",zl);
    printf("  z_s           = %f\n",zs);
    printf("  Omega_m       = %f\n",omega0);
    printf("  Omega_Lambda  = %f\n",lambda);
    printf("  sigma_V       = %5.0f km/sec\n\n",sigv);
    printf("Calculated quantities:\n");
    printf("   D_l           = %4.0f h^{-1} Mpc\n",dal);
    printf("   D_s           = %4.0f h^{-1} Mpc\n",das);
    printf("   D_ls          = %4.0f h^{-1} Mpc\n",dals);
    printf("   D_ls / D_s    = %6.3f\n",dals/das);
    printf("   theta_E       = %6.3f arcsec\n",theta_E);
    if(offset > theta_E) 
      printf("   shear         = %5.3f\n",gamma);
    printf("\n--------------------------------------------------\n\n");

    /*
     * Continue?
     */

    printf("Another (Omega_m, Omega_Lambda)? (y/n) [y] ");
    fgets(line,MAXC,stdin);
    if(line[0] == 'N' || line[0] == 'n')
      contin = 0;
  }

  printf("\n");
  return 0;
}

