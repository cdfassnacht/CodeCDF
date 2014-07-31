/* wfc3_texp.c
 *
 * Calculates an exposure time for WFC3 given an input count rate and desired
 *  SNR.
 *
 */

#include <stdio.h>
#include <math.h>

#define DIAM 250     /* Diameter of the telescope, in cm */
#define SCAL 0.61    /* Ratio of Destiny A-Omega to HST/WFC3 A-Omega */
#define QE 0.449     /* Quantum efficiency at ~1.1 microns */
#define RN 16.0      /* Read noise per pixel, in electrons */
/* #define NPIX 16.0  */  /* Number of pixels  */
#define PIX 0.13     /* Pixel scale in arcsec */
#define SKY 0.023    /* Sky noise in photons/sec/Angstrom/(sq. arcsec) */
#define BW  5060.0   /* Approximate F110W bandwidth */
#define DC 0.12      /* Dark current in electrons/s/pix */
#define THERM 0.1    /* Thermal background in phot/s/pix */
#define MAXC 100     /* Max number of characters in a line */

int main(int argc, char *argv[])
{
  float texp;   /* Exposure time */
  float npix;   /* Number of pixels covered by source */
  float rsrc;   /* Count rate from source */
  float rsky;   /* Count rate from sky */
  float rtherm; /* Count rate from thermal background */
  float rdc;    /* Count rate from dark current */
  float Nsrc;   /* Number of counts from source */
  float Nsky;   /* Number of counts from sky */
  float Nbkgd;  /* Number of counts from thermal background */
  float Ndc;    /* Number of counts from dark current */
  float sn;     /* Desired SNR */
  float sn2;    /* SNR*SNR */
  float noise;  /* Total noise in observation */
  float tmp;
  float tmpa;   /* Temporary holder for solving quadratic equation */
  float tmpb;   /* Temporary holder for solving quadratic equation */
  float tmpc;   /* Temporary holder for solving quadratic equation */
  char line[MAXC]; /* General string varaible */

  /*
   * Get inputs 
   */

  printf("\nCount rate from source in counts/sec\n");
  printf("  Get this value from wfc3tp.pro in IDL/destiny directory.\n");
  printf("  Enter value here:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f",&rsrc) != 1) {
    fprintf(stderr,"ERROR: bad input\n");
    fprintf(stderr,"Enter count rate again:  ");
    fgets(line,MAXC,stdin);
  }

  printf("Enter number of pixels covered by source:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f",&npix) != 1) {
    fprintf(stderr,"ERROR: bad input\n");
    fprintf(stderr,"Enter number of pixels again:  ");
    fgets(line,MAXC,stdin);
  }

  printf("Enter desired SNR:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f",&sn) != 1) {
    fprintf(stderr,"ERROR: bad input\n");
    fprintf(stderr,"Enter SNR again:  ");
    fgets(line,MAXC,stdin);
  }

  /*
   * Convert inputs and detector/telescope parameters into count rates
   */

  rsky = SCAL * QE * SKY * BW * npix * PIX * PIX;
  rdc = DC * npix;
  rtherm = QE * THERM * npix;
  tmp = rsrc + rsky + rdc + rtherm;
  sn2 = sn*sn;

  /*
   * Put counts into standard variables for solving the quadratic equation
   */

  tmpa = rsrc * rsrc;
  tmpb = -1.0 * sn2 * tmp;
  tmpc = -1.0 * sn2 * npix * RN * RN;

  /*
   * Calculate exposure time needed to achieve the requested SNR
   */

  texp = (-1.0*tmpb + sqrt(tmpb*tmpb - 4.0 * tmpa * tmpc)) / (2.0 * tmpa);
  noise = sqrt(tmp*texp + npix * RN * RN);

  /*
   * Print out results
   */

  printf("\nFor this observation the input parameters are:\n");
  printf(" Source count rate     = %g e-/sec\n",rsrc);
  printf(" Sky count rate        = %g e-/sec\n",rsky);
  printf("    Sky                = %f phot/sec/Ang/arcsec^2\n",SKY);
  printf("    QE                 = %f\n",QE);
  printf("    Filter bandwidth   = %f Ang\n",BW);
  printf("    Pixel size         = %f arcsec\n",PIX);
  printf("    N_pix in source    = %f\n",npix);
  printf(" Thermal count rate    = %f e-/sec\n",rtherm);
  printf(" Dark current rate     = %f e-/sec\n",rdc);
  printf(" Read noise            = %f e-\n",RN);
  printf("Giving:\n");
  printf("  Source counts     = %g cts\n",rsrc*texp);
  printf("  Sky counts        = %g cts\n",rsky*texp);
  printf("  Thermal counts    = %g cts\n",rtherm*texp);
  printf("  Dark counts       = %g cts\n",rdc*texp);
  printf("Associated noise\n");
  printf("  source            = %g cts\n",sqrt(rsrc*texp));
  printf("  sky               = %g cts\n",sqrt(rsky*texp));
  printf("  thermal           = %g cts\n",sqrt(rtherm*texp));
  printf("  dark current      = %g cts\n",sqrt(rdc*texp));
  printf("  read noise        = %g cts\n",RN*sqrt(npix));
  printf(" Total noise        = %g cts\n\n",noise);
  printf("Which means you need an exposure time of %6.0f sec to get\n",texp);
  printf(" a SNR of %6.2f\n",sn);

  return 0;
}
