/*
 * optmags.c
 *
 * This program calculates the absolute magnitude of an object given
 *  its apparent magnitude, a redshift and a k-correction.
 *
 * 17Oct1997 CDF,  First working version
 * v09Aug1998 CDF, Added capacity for making calculations in non-zero
 *                  lambda cosmologies by incorporating the revised
 *                  ang_dist function and the new get_omega function from
 *                  cosmo.c
 * v10Aug1998 CDF, Added inclusion of a color term and then calculation of
 *                  a luminosity in the resulting band.
 * v02Feb2007 CDF, Replaced old calculation of D_L with cosmo library
 *                  function calc_cosdist
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "dataio.h"
#include "cosmo.h"

int main(int argc, char *argv[])
{
  float z;          /* Redshift of source */
  float k;          /* K-correction */
  float appmag;     /* Apparent magnitude of source */
  float absmag;     /* Absolute magnitude of source */
  float color;      /* Color term used to convert to new band */
  float lum;        /* Luminosity in new band */
  char line[MAXC];  /* General string for reading input */
  Cosmo cosmo;      /* Cosmological world model */
  Cosdist cosdist;  /* Structure for all of the distance measures */

  /*
   * Get apparent magnitude of source
   */

  printf("Enter the apparent magnitude of the source:  ");
  gets(line);
  while(sscanf(line,"%f",&appmag) != 1) {
    fprintf(stderr,"ERROR.  Bad input.  Enter new value:  ");
    gets(line);
  }

  /*
   * Get redshift of source
   */

  printf("Enter the redshift of the source:  ");
  gets(line);
  while(sscanf(line,"%f",&z) != 1 || z < 0.0) {
    fprintf(stderr,"ERROR.  Bad input.  Enter new value:  ");
    gets(line);
  }

  /*
   * Get k-correction
   */

  printf("Enter the k-correction for the source:  ");
  gets(line);
  while(sscanf(line,"%f",&k) != 1) {
    fprintf(stderr,"ERROR.  Bad input.  Enter new value:  ");
    gets(line);
  }

  /*
   * Initialize world model
   */

  cosmo.omega_m = 0.3;
  cosmo.omega_de = 0.7;
  cosmo.w = -1.0;

  /*
   * Get the cosmological world model
   */

  get_cosmo(&cosmo);

  /*
   * Calculate the distance measures
   */

  cosdist = calc_cosdist(0,z,cosmo);

  /*
   * Convert the luminosity distance to Mpc
   */

  cosdist.d_l /= MPC2CM;
  printf("\nWith the entered values, the luminosity distance is ");
  printf("%6.1f h^{-1} Mpc\n\n",cosdist.d_l);

  /*
   * Calculate the absolute magnitude, M = m - 5 log (D_L / 10pc) - K
   */

  absmag = appmag - 5*log10(cosdist.d_l * 1.0e5) - k;

  printf("For this source we have:\n");
  printf("  Apparent magnitude = %7.3f\n",appmag);
  printf("  Distance modulus   = %7.3f\n",5*log10(cosdist.d_l*1.0e5));
  printf("  k-correction       = %7.3f\n",k);
  printf("  Absolute magnitude = %7.3f + 5 log(h)\n\n",absmag);

  /*
   * Get the color term
   */

  printf("Enter color term C in sense that M_new = M_old + C:  ");
  gets(line);
  while(sscanf(line,"%f",&color) != 1) {
    fprintf(stderr,"ERROR.  Bad input.  Enter color again:  ");
    gets(line);
  }

  /*
   * Calculate absolute magnitude in new band and get luminosity
   */

  if((lum = mag_to_lum(absmag+color)) >= 0.0) {
    printf("\n  Absolute magnitude in new band = %7.3f + 5 log(h)\n",
	   absmag+color);
    printf("  Luminosity in new band         = %7.2e h^{-2} L_sun\n",lum);
  }
  else {
    printf("\n  Absolute magnitude in new band = %7.3f + 5 log(h)\n",
	   absmag+color);
    printf("  Cannot calculate luminosity in that band yet.\n");
  }

  return 0;
}
