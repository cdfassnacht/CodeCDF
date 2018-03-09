/*
 * cosmo.c
 *
 * A collection of functions used to calculate cosmological quantities.
 *
 *    get_valerr
 *    get_cosmo
 *    calc_cosdist
 *    peebles_E
 *    mn_sis
 *    mag_to_lum
 *
 * Revision history:
 *  05Jul1997  CDF, First split out of other GL programs.
 *  v16Oct1997 CDF, Added error calculation for angular diameter distances.
 *  v05Aug1998 CDF, Moved calculation for Lambda = 0 case into 
 *                   ang_dist_no_lambda
 *                  Added ang_dist_lambda for non-zero Lambda cosmologies.
 *  v09Aug1998 CDF, Moved the cosmology inputs from phys_units.c (and optmags.c)
 *                   into the new get_omega function.
 *  v10Aug1998 CDF, Added mag_to_lum to convert absolute magnitudes into
 *                   luminosities in solar units.
 *  v09Aug2005 CDF, Added new peebles_E function to calculate E(z)
 *                  Added new hz function to calculate H(z)
 *                  Added notes on the end for age of Universe calculations
 *  v11Jul2006 CDF, Added new function get_cosmo that takes the place of
 *                   the old get_omega function.
 *                  Added new function calc_cosdist that takes the place of
 *                   the old ang_dist, lookback_time, and hz functions while
 *                   adding more capability.
 *                   **NB: This new function does not include errors.  See
 *                     GL notebook #5, p. 72 for formula.
 *                  Modified peebles_E to use dark energy with constant w
 *                   but where w is not necessarily -1.
 *  v25Jan2007 CDF, Added the new mn_sis function
 *  v02Feb2007 CDF, Removed the old get_omega and hz functions after fixing all
 *                   programs that called them.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cosmo.h"

#define TINY 1.0e-18       /* To avoid division by zero */
#define MAXLINE 100
#define NSTEP 100          /* Number of steps in numerical integrations */
#define DZ 0.001


/*.......................................................................
 *
 * Function get_valerr
 *
 * Gets a value and an optional error on the value.
 *
 * Inputs: double *val         value to be set by this function
 *         double *dval        error on the value (optional)
 *         char *valdef        string defining the value (e.g. "redshift" 
 *                               or "velocity dispersion")
 *
 * Output: int nget            number of quantities set by this function
 *                              (should be 1 or 2).  0 ==> error.
 *
 */

int get_valerr(double *val, double *dval, char *valdef)
{
  int readok=0;       /* Flag set to 1 when data successfully read */
  int nget=0;         /* Number of quantities set by this function */
  char line[MAXLINE]; /* General string for reading variables */

  readok = 0;
  printf("\nEnter value and optional error for the %s:  ",valdef);
  while(readok == 0) {
    fgets(line,MAXLINE,stdin);
    switch(sscanf(line,"%lf %lf",val,dval)) {
    case 1:
      if(*val < 0.0) {
	fprintf(stderr,"ERROR: improper input.  Enter value(s) again:  ");
	break;
      }
      else {
	readok = 1;
	nget = 1;
	break;
      }
    case 2:
      if(*val < 0.0 || *dval < 0.0) {
	fprintf(stderr,"ERROR: improper input.  Enter value(s) again:  ");
	break;
      }
      else {
	readok = 1;
	nget = 2;
	break;
      }
    default:
      fprintf(stderr,"ERROR: improper input.  Enter value(s) again:  ");
    }
  }

  return nget;
}

/*.......................................................................
 *
 * Function get_cosmo
 *
 * Gets the cosmological parameters for the world model - namely
 *  Omega_m, Omega_DE, and w.
 *
 * Inputs: Cosmo *cosmo        Structure containing world model
 *
 * Output: (none)
 *
 */

void get_cosmo(Cosmo *cosmo)
{
  char line[MAXLINE];    /* General string for input */

  printf("\nSet cosmological world model:\n");

  /*
   * Update values
   */

  printf("  Enter Omega_m [%4.2f]:  ",cosmo->omega_m);
  fgets(line,MAXLINE,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%lf",&cosmo->omega_m) != 1 || cosmo->omega_m < 0.0) {
      fprintf(stderr,"  ERROR: improper input.  Enter value again:  ");
      fgets(line,MAXLINE,stdin);
    }
  }

  printf("  Enter Omega_DE [%4.2f]:  ",cosmo->omega_de);
  fgets(line,MAXLINE,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%lf",&cosmo->omega_de) != 1) {
      fprintf(stderr,"  ERROR: improper input.  Enter value again:  ");
      fgets(line,MAXLINE,stdin);
    }
  }

  if(cosmo->omega_de != 0.0) {
    printf("  Enter w [%4.2f]:  ",cosmo->w);
    fgets(line,MAXLINE,stdin);
    if(line[0] != '\n') {
      while(sscanf(line,"%lf",&cosmo->w) != 1) {
	fprintf(stderr,"  ERROR: improper input.  Enter value again:  ");
	fgets(line,MAXLINE,stdin);
      }
    }
  }
}

/*.......................................................................
 *
 * Function calc_cosdist
 *
 * Given a cosmology, calculates various distance measures between two
 *  redshifts (z1 and z2).  
 *  Formulas are from Hogg's cosmological distance measures writeup:
 *   astro-ph/9905116
 *  **NB: This new function does not include errors.  See
 *    GL notebook #5, p. 72 for formula.
 *
 * Given the redshift in the input redshifts, this function will
 *  calculate the following quantities (all between z1 and z2):
 *    D_C  - Comoving distance (line of sight)
 *    D_M  - Comoving distance (for calculating transverse separations)
 *    D_A  - Angular diameter distance
 *    D_L  - Luminosity distance
 *    DM   - Distance modulus
 *    t_L  - Lookback time
 *    H(z) - Hubble constant at z2
 *
 * NOTE: if either of the input redshifts is less than 0.0, then set all
 *  values = -99.0
 *
 *
 * Inputs: double z1           lower redshift
 *         double z2           lower redshift
 *         Cosmo cosmo         cosmological world model
 *
 * Output: Cosdist *cosdist    structure containing distance measures
 *
 */

Cosdist calc_cosdist(double z1, double z2, Cosmo cosmo)
{
  int i=0;              /* Looping variable */
  double ztmp;          /* Stepped value of redshift */
  double dsum=0.0;      /* Running total for distance numerical integration */
  double tsum=0.0;      /* Running total for time numerical integration */
  double omega_k;       /* Omega_k = 1 - Omega_m - Omega_Lambda */
  double ez;            /* Peebles' E(z) */
  double sqrtomk;       /* Square root of absolute value of Omega_k */
  double zp1;           /* 1 + z */
  Cosdist cosdist;      /* Structure to be filled with distance measures */

  /*
   * First check that both input redshifts are valid.
   * If either one is less than 0, then set all the cosdist values to -99
   *  and return.
   */

  if(z1<0. || z2<0.) {
    cosdist.ez  = -99.;
    cosdist.hz  = -99.;
    cosdist.t_l = -99.;
    cosdist.d_c = -99.;
    cosdist.d_m = -99.;
    cosdist.d_a = -99.;
    cosdist.d_l = -99.;
    cosdist.DM  = -99.;
    return cosdist;
  }

  /*
   * Calculate Omega_k (curvature density) for this cosmology
   */

  omega_k = 1.0 - cosmo.omega_m - cosmo.omega_de;
  sqrtomk = sqrt(fabs(omega_k));

  /*
   * Calculate quantities that do not need numerical integration, 
   *  namely E(z) and H(z)
   */

  cosdist.ez = peebles_E(z2,cosmo);
  cosdist.hz = 100.0 * cosdist.ez;

  /*
   * Set initial step in numerical integration.  When we do the
   *  numerical integration, we will step through the redshift
   *  range in steps of dz, centering each step halfway between
   *  integral multiples of dz.
   */

  ztmp = z1 + 0.5*DZ;

  /*
   * Numerically integrate.  The only quantities that require
   *  numerical integration are D_C and t_L.  The other quantities
   *  can be derived from these.  The formulas are:
   *
   *          c   (z2    dz'
   *   D_C = ---  |    ------
   *         H_0  )z1   E(z')
   *
   * and
   *
   *          1   (z2       dz'
   *   t_l = ---  |    --------------
   *         H_0  )z1  (1 + z') E(z')
   */

  while(ztmp < z2) {
    ztmp += DZ;
    ez = peebles_E(ztmp,cosmo);
    dsum += DZ / ez;
    tsum += DZ / ((1.0 + ztmp) * ez);
    i++;
  }

  cosdist.d_c = C * dsum / H0;
  cosdist.t_l = tsum / H0;

  /*
   * Calculate the transverse comoving distance (D_M) from D_C.
   * This will depend on the curvature.
   */

  if(omega_k > 0.0)
    cosdist.d_m = C * sinh(sqrtomk * cosdist.d_c * H0 / C) / (H0 * sqrtomk);
  else if(omega_k < 0.0)
    cosdist.d_m = C * sin(sqrtomk * cosdist.d_c * H0 / C) / (H0 * sqrtomk);
  else
    cosdist.d_m = cosdist.d_c;

  /*
   * Calculate the remaining distances from D_M
   */

  zp1 = 1.0 + z2;
  cosdist.d_a = cosdist.d_m / zp1;
  cosdist.d_l = cosdist.d_m * zp1;

  /*
   * Calculate the distance modulus from D_L
   *  The factor of 3.1e19 comes from converting the 10pc in the
   *  denominator into cm.
   */

  cosdist.DM = 5.0 * log10(cosdist.d_l / 3.1e19);

  /*
   * Return to calling function
   */

  return cosdist;
}


/*.......................................................................
 *
 * Function peebles_E
 *
 * Calculates E(z) as defined in Peebles "Physical Cosmology"
 *  This function is proportional to (da/dt)/a, and thus is
 *  useful in calculating H(z) and all sorts of distance measures.
 *
 * E(z) can be written several different ways.  However, following
 *  Peebles, take
 *
 *   E(z) = sqrt[omega_m*(1+z)^3 + omega_k*(1+z)^2 + omega_DE*(1+z)^{3*(1+w)}]
 *
 *     where omega_k = 1 - (omega_m + omega_lam)
 *
 * Inputs: double z              redshift
 *         Cosmo cosmo           cosmological world model
 *
 * Output: double E              E(z)
 *
 */

double peebles_E(double z, Cosmo cosmo)
{
  double g;             /* Alias for 1 + z */
  double g2;            /* Square of g */
  double omega_k;       /* Omega_k = 1 - (Omega_m - Omega_DE) */
  double de_term;       /* Dark energy multiplier */
  double E;             /* E(z) */

  omega_k = 1.0 - cosmo.omega_m - cosmo.omega_de;
  g = 1.0 + z;
  g2 = g * g;
  if(cosmo.omega_de != 0.0)
    de_term = pow(g,3.0 * (1.0 + cosmo.w));
  else
    de_term = 1.0;
  E = sqrt(cosmo.omega_m * g * g2 + omega_k * g2 + cosmo.omega_de * de_term);

  return E;
}

/*.......................................................................
 *
 * Function rn_generic
 *
 * Calculates r_n (e.g., r_200) for a halo given a halo mass (M_n) in solar
 * masses.  This calculation uses the following formula:
 *
 *   r_n = (2 G M_n / n H(z)^2)^(1/3)
 *
 * Inputs: double z             redshift of halo
 *         Cosmo cosmo          world model
 *         double m_n           mass of the halo, in solar masses
 *         double n             factor by which the average density within
 *                               r_n is greater than the critical density
 *                               of the Universe at redshift z.
 *
 * Output: r_n                  radius enclosing the overdensity.
 *
 */

double rn_generic(double z, Cosmo cosmo, double m_n, double n)
{
  double hz;       /* H(z) in km/s/Mpc */
  double r_n;      /* Radius enclosing the overdensity */

  /*
   * Calculate H(z) and r_n
   * Note that M_n has to be converted to grams by multiplying by MSUN
   */

  hz = 100.0 * peebles_E(z,cosmo);
  r_n = 2.0 * G * m_n * MSUN / (n * hz * hz);
  r_n = pow(r_n,(1.0/3.0));

  return r_n;
}

/*.......................................................................
 *
 * Function rn_sis
 *
 * Calculates r_n (e.g., r_200) for a SIS model.  This calculation uses
 *  the value of M_n derived by a direct integration of rho_sis(r) from
 *  Narayan and Bartelmann.
 *
 * Inputs: double z             redshift of halo
 *         Cosmo cosmo          world model
 *         double sigma         velocity dispersion of the halo, in km/s
 *         double n             factor by which the average density within
 *                               r_n is greater than the critical density
 *                               of the Universe at redshift z.
 *
 * Output: r_n                  radius enclosing the overdensity.
 *
 */

double rn_sis(double z, Cosmo cosmo, double sigma, double n)
{
  double r_n;   /* Radius enclosing the overdensity */

  r_n = rn_isothermal(z,cosmo,sigma,n,2.0);

  return r_n;
}

/*.......................................................................
 *
 * Function rn_carlberg
 *
 * Calculates r_n (e.g., r_200) for a SIS model.  This calculation uses
 *  the value of M_n taken as the virial mass from Carlberg, Yee and
 *  Ellingson (1997)
 *
 * Inputs: double z             redshift of halo
 *         Cosmo cosmo          world model
 *         double sigma         velocity dispersion of the halo, in km/s
 *         double n             factor by which the average density within
 *                               r_n is greater than the critical density
 *                               of the Universe at redshift z.
 *
 * Output: r_n                  radius enclosing the overdensity.
 *
 */

double rn_carlberg(double z, Cosmo cosmo, double sigma, double n)
{
  double r_n;   /* Radius enclosing the overdensity */

  r_n = rn_isothermal(z,cosmo,sigma,n,3.0);

  return r_n;
}

/*.......................................................................
 *
 * Function rn_isothermal
 *
 * Calculates r_n (e.g., r_200) for a halo given a halo mass (M_n) in solar
 *  masses.  This calculation assumes that the halo has a singular isothermal
 *  sphere profile with a 1D velocity dispersion sigma.  
 * Generically, this can be expressed as
 *
 *   r_n = sqrt(2*K/n) sigma / H(z)
 *
 * where the K represents the confusion between a straight integration of
 * the SIS profile from Narayan and Bartelmann (K = 2) or the "virial mass"
 * value from Carlberg, Yee, & Ellingson (1997; K = 3).
 *
 * Inputs: double z             redshift of halo
 *         Cosmo cosmo          world model
 *         double sigma         velocity dispersion of the halo, in km/s
 *         double n             factor by which the average density within
 *                               r_n is greater than the critical density
 *                               of the Universe at redshift z.
 *         double K             factor to differentiate between SIS and 
 *                               Carlberg et al.
 *
 * Output: r_n                  radius enclosing the overdensity.
 *
 */

double rn_isothermal(double z, Cosmo cosmo, double sigma, double n, double K)
{
  double hz;       /* H(z) in km/s/Mpc */
  double r_n;      /* Radius enclosing the overdensity */

  /*
   * Calculate H(z) and r_n
   * Note that sigma has to be converted to cm/s by multiplying by 10^5
   */

  hz = 100.0 * peebles_E(z,cosmo);
  r_n = sqrt(2.0 * K / n) * sigma * 1.0e5 / hz;

  return r_n;
}

/*.......................................................................
 *
 * Function mn_sis
 *
 * Calculates M_n (e.g., M_200) for a SIS mass distribution.
 * Uses M_n(sis) = 4 sigma^3 / G * H(z) * sqrt(n).
 *
 * Inputs: double z             redshift of halo
 *         Cosmo cosmo          world model
 *         double sigma         SIS 1-D velocity dispersion in km/s
 *         double n             factor by which the average density within
 *                               r_n is greater than the critical density
 *                               of the Universe at redshift z.
 *    
 * Output: double m_n           mass of the halo within r_n
 *
 */

double mn_sis(double z, Cosmo cosmo, double sigma, double n)
{
  double sigcm;    /* Velocity dispersion in cm/s */
  double hz;       /* H(z) */
  double m_n;     /* M_n */

  /*
   * Do the calculations
   */

  sigcm = sigma * 1.0e5;
  hz = 100.0 * peebles_E(z,cosmo);
  m_n = 4.0 * sigcm * sigcm * sigcm / (G * hz * KMSMPC2S * sqrt(n));

  /*
   * Convert M_n to solar masses
   */

  m_n /= MSUN;

  return m_n;
}

/*.......................................................................
 *
 * Function mag_to_lum
 *
 * Converts an absolute magnitude into a luminosity in terms of solar
 *  luminosities.
 *
 * Inputs: double absmag        absolute magnitude
 *
 * Output: double lum           luminosity (<0.0 on error)
 *
 */

double mag_to_lum(double absmag)
{
  int band=OTHER;     /* Band in which absmag is measured */
  double lum;          /* Luminosity in solar luminosities */
  double msol;         /* Solar magnitude in band of choice */
  char line[MAXLINE]; /* String variable for reading input */

  /*
   * Get the band in which absmag is calculated
   */

  printf("Band in which absolute magnitude is given:\n");
  printf("  %d. U\n",U);
  printf("  %d. B\n",B);
  printf("  %d. V\n",V);
  printf("  %d. R\n",R);
  printf("  %d. I\n",I);
  printf("  %d. Other\n",OTHER);
  printf(" Enter choice:  [%d] ",band);
  fgets(line,MAXLINE,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%d",&band) != 1 || band < OTHER || band > I) {
      fprintf(stderr,"ERROR: Bad value for input band.  Enter new value:  ");
      fgets(line,MAXLINE,stdin);
    }
  }

  /*
   * Set the appropriate solar luminosity
   */

  switch(band) {
  case U:
    msol = MSUN_U;
    break;
  case B:
    msol = MSUN_B;
    break;
  case V:
    msol = MSUN_V;
    break;
  case R:
    msol = MSUN_R;
    break;
  case I:
    msol = MSUN_I;
    break;
  case OTHER:
    printf("Other bands are not yet supported.\n");
    return -1.0;
    break;
  default:
    printf("Other bands are not yet supported.\n");
    return -1.0;
  }

  /*
   * Calculate the luminosity
   */

  lum = pow(10.0,0.4*(msol-absmag));

  return lum;
}

/*.......................................................................
 *
 * Age of Universe
 *
 * Use analytic equations from Kolb & Turner if appropriate, otherwise
 *  use approximation from Carroll et al. Lambda paper.
 *
 */

