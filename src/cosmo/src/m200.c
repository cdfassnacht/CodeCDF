/*
 * m200.c
 *
 * This program computes M_200 for an isothermal halo of a given velocity
 *  dispersion.
 *
 * First working version: 25Jan2007 by Chris Fassnacht (CDF)
 * Revision history:
 *  02Feb2007 CDF, Changed to reflect the new, more general function mn_iso
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cosmo.h"

#define MAX 80
#define PI 3.141592653589793

/*.......................................................................
 *
 * Main program
 *
 */

int main(int argc, char *argv[])
{
  int contin=1;       /* Flag set to 0 to break loop */
  double z;           /* Lens redshift */
  double dz=0.0;      /* Error on lens redshift */
  double omega0=0.3;  /* Omega_0 (default starting value is approx WMAP) */
  double lambda=0.7;  /* Omega_Lambda (default starting value is approx WMAP) */
  double dl;          /* Luminosity distance */
  double dm;          /* Distance modulus */
  double t_l;         /* Lookback time */
  double sigma;       /* 1 dimensional velocity dispersion for the halo */
  double dsig;        /* Optional error on sigma */
  double m200;        /* M_200 */
  Cosmo cosmo;        /* Cosmological world model */
  Cosdist cosdist;    /* Structure for all of the distance measures */
  char line[MAX];     /* General string for reading input */

  /*
   * Get redshift and velocity dispersion
   */

  if(get_valerr(&z,&dz,"source redshift") == 0) {
    fprintf(stderr,"\nERROR.  Exiting m200.\n\n");
    return 1;
  }

  if(get_valerr(&sigma,&dsig,"velocity dispersion") == 0) {
    fprintf(stderr,"\nERROR.  Exiting m200.\n\n");
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
     * Calculate M_200 for a SIS
     */

    m200 = mn_sis(z,cosmo,sigma,200.0);

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
    printf("   M_200                     = %g h M_solar\n",m200);
    printf("\n--------------------------------------------------\n\n");

    /*
     * Continue?
     */

    printf("Another cosmological world model? (y/n) [y] ");
    fgets(line,MAX,stdin);
    if(line[0] == 'N' || line[0] == 'n')
      contin = 0;
  }

  printf("\n");
  return 0;
}

