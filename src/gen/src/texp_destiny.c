/* texp_destiny.c
 *
 * Calculates an exposure time for WFC3 given an input count rate and desired
 *  SNR.
 *
 */

#include <stdio.h>
#include <math.h>

#define DIAM 250     /* Diameter of the telescope, in cm */
#define QE 0.449     /* Quantum efficiency at ~1.1 microns */
#define RN 16.0      /* Read noise per pixel, in electrons */
#define PIX 0.13     /* Pixel scale in arcsec */
#define SKYP 0.023    /* Sky noise in photons/sec/Angstrom/arcsec^ */
#define FSKY 1.0e-18 /* Sky flux at F110W in erg/s/cm^/Ang/arcsec^2 */
#define BW  5060.0   /* Approximate F110W bandwidth */
#define DC 0.12      /* Dark current in electrons/s/pix */
#define THERM 0.1    /* Thermal background in phot/s/pix */
#define MAXC 100     /* Max number of characters in a line */
#define PI 3.141592653589793

int main(int argc, char *argv[])
{
  int instrument; /* Observing instrument */
  int flformat;   /* Format of input flux */
  float c;        /* Speed of light in cm/s */
  float h;        /* Planck's constant in erg s */
  float lambda_e; /* Effective wavelength of F110W filter */
  float influx;   /* Input source flux or surface brightness */
  float sbsrc;    /* Source surface brightness in flux units */
  float npix;     /* Number of detector pixels covered by the source */
  float srcarea;  /* Area of source, in arcsec^2 */
  float srcflux;  /* Flux density of source, in erg/s/cm^2/Ang */
  float srcmag;   /* Observed magnitude of the source, in AB mags */
  float srcmu;    /* Observed surface brightness, in AB mags/arcsec^2 */
  float mag2flux; /* Exponent in converting mags to flux units */
  float psf;      /* Sigma of PSF, in arcsec */
  float gain;     /* Detector gain in ADU/e- */
  float qe;       /* Detector qe in e-/photon */
  float eta;      /* Throughput of telescope/camera/filter assembly */
  float tptot;    /* Total throughput of telescope + detector qe */
  float dark;     /* Dark current, in e-/s/pix */
  float therm;    /* Thermal background, in e-/s/pix */
  float readnoise;     /* Read noise, in e- */
  float pixsize;  /* Detector pixel size */
  float atel;     /* Telescope primary area in cm^2 */
  float dtel;     /* Telescope primary area in cm */
  float bw;       /* Filter bandwidth */
  float texp;     /* Exposure time */
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
   * Initialize
   */

  c = 3.0e10;
  h = 6.6e-27;
  lambda_e = 1.17e-4;  /* Effective wavelength of F110W filter, in cm */

  /*
   * Get instrument
   */

  printf("Choose instrument:\n");
  printf("   1. NICMOS/NIC3\n");
  printf("   2. WFC3\n");
  printf("   3. Destiny\n");
  printf(" Enter choice:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%d",&instrument) != 1) {
    fprintf(stderr,"ERROR: bad input\n");
    fprintf(stderr,"Enter instrument again:  ");
    fgets(line,MAXC,stdin);
  }

  /*
   * Set instrument parameters
   */

  switch(instrument) {
  case 1: /* NICMOS/NIC3 */
    gain = 6.5;          /* gain for NICMOS/NIC3 */
    qe = 0.3;            /* QE for NICMOS/NIC3 at F110W */
    eta = 0.744 * 0.84;  /* Telescope throughput times filter throughput */
    tptot = qe * eta;    /* Total throughput - photons to e- */
    readnoise = 22.0;         /* Read noise, in e- */
    dark = 0.3;          /* Dark current, in e-/s/pix */
    therm = 9.821e-5;    /* Thermal background, in e-/s/pix */
    pixsize = 0.200;     /* Pixel size in arcsec */
    atel = 45238.9;      /* HST area in cm^2 */
    dtel = 240;          /* HST diameter in cm */
    bw = 6000.0;         /* Assumed filter bandwidth in Ang */
    break;
  case 2: /* WFC3 */
    qe = 0.42;           /* Average QE */
    tptot = 0.24;        /* Average total throughput from IDL task */
    readnoise = 15.6;         /* Read noise, in e- */
    dark = 0.12;         /* Dark current, in e-/s/pix */
    therm = qe * 0.1;    /* Thermal background, in e-/s/pix */
    pixsize = 0.130;     /* Pixel size in arcsec */
    atel = 45238.9;      /* HST area in cm^2 */
    dtel = 240;          /* HST diameter in cm */
    bw = 5060.0;         /* Assumed filter bandwidth in Ang */
    break;
  case 3: /* Destiny - assume WFC3 except for telescope area and pixel size */
    qe = 0.42;           /* Average QE */
    tptot = 0.24;        /* Average total throughput from IDL task */
    readnoise = 15.6;         /* Read noise, in e- */
    dark = 0.12;         /* Dark current, in e-/s/pix */
    therm = qe * 0.1;    /* Thermal background, in e-/s/pix */
    pixsize = 0.150;     /* Pixel size in arcsec */
    atel = 21382.0;      /* HST area in cm^2 */
    dtel = 165;          /* HST diameter in cm */
    bw = 5060.0;         /* Assumed filter bandwidth in Ang */
    break;
  default:
    fprintf(stderr,"\nERROR: Bad choice of instrument.\n");
    fprintf(stderr,"Exiting texp_destiny.\n\n");
    return 1;
  }
  psf = 1.22 * 1.2e-4 * 206265 / dtel; /* PSF sigma, in arcsec, for F110W */
  printf("\nPSF = %f arcsec = %f pix\n",psf,psf/pixsize);
  printf(" so 1 pix = %f psf\n",pixsize/psf);

  /*
   * Get input flux
   */


  printf("\nChoose units for source flux or surface brightness:\n");
  printf("   1. Surface brightness in erg/s/cm^2/Ang/arcsec^2\n");
  printf("   2. Surface brightness in AB mag/arcsec^2\n");
  printf("   3. Flux density in erg/s/cm^2/Ang\n");
  printf("   4. Flux density in AB mag\n");
  printf(" Enter choice:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%d",&flformat) != 1 || flformat<1 || flformat>4) {
    fprintf(stderr,"ERROR: bad input\n");
    fprintf(stderr,"Enter unit format again:  ");
    fgets(line,MAXC,stdin);
  }

  printf("\nEnter area on the detector covered by the source, in pixels:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f",&npix) != 1) {
    fprintf(stderr,"ERROR: bad input\n");
    fprintf(stderr,"Enter value again:  ");
    fgets(line,MAXC,stdin);
  }
  srcarea = npix * pixsize * pixsize;
  printf("Solid angle covered by source = %f arcsec^2.\n",srcarea);


  printf("Now enter value of source flux/surface brightness:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f",&influx) != 1) {
    fprintf(stderr,"ERROR: bad input\n");
    fprintf(stderr,"Enter value again:  ");
    fgets(line,MAXC,stdin);
  }

#if 0
  printf("\nCount rate from source in counts/sec\n");
  printf("  Get this value from wfc3tp.pro in IDL/destiny directory.\n");
  printf("  Enter value here:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f",&rsrc) != 1) {
    fprintf(stderr,"ERROR: bad input\n");
    fprintf(stderr,"Enter count rate again:  ");
    fgets(line,MAXC,stdin);
  }
#endif
  /*
   * Convert input flux or surface brightness into surface brightness
   *  in flux units.
   */

  switch(flformat) {
  case 1:
    srcflux = influx * srcarea;
    srcmag = -2.5 * log10(srcflux) - 21.1;
    srcmu = -2.5 * log10(influx) - 21.1;
    break;
  case 2:
    mag2flux = (influx + 21.0)/-2.5;
    srcflux = pow(10.0,mag2flux) * srcarea;
    srcmag = -2.5 * log10(srcflux) - 21.1;
    srcmu = influx;
    break;
  case 3:
    srcflux = influx;
    srcmag = -2.5 * log10(srcflux) - 21.1;
    srcmu =  -2.5 * log10(srcflux/srcarea) - 21.1;
    break;
  case 4:
    srcmag = influx;
    mag2flux = (influx + 21.0)/-2.5;
    srcflux = pow(10.0,mag2flux);
    srcmu =  -2.5 * log10(srcflux/srcarea) - 21.1;
    break;
  default:
    fprintf(stderr,"\nERROR. Not a valid choice for source units.\n");
    fprintf(stderr,"Exiting program\n");
    return 0;
  }

  /*
   * Get desired SNR
   */

  printf("Enter desired SNR:  ");
  fgets(line,MAXC,stdin);
  while(sscanf(line,"%f",&sn) != 1) {
    fprintf(stderr,"ERROR: bad input\n");
    fprintf(stderr,"Enter SNR again:  ");
    fgets(line,MAXC,stdin);
  }

  /*
   * Convert inputs and detector/telescope parameters into count rates
   *  in e-/sec
   */

  rsrc = srcflux * tptot * atel * bw * lambda_e / (h * c);
  rsky = FSKY * srcarea * tptot * atel * bw * lambda_e / (h * c);
  rdc = dark * npix;
  rtherm = therm * npix;
  tmp = rsrc + rsky + rdc + rtherm;
  sn2 = sn*sn;

  /*
   * Put counts into standard variables for solving the quadratic equation
   */

  tmpa = rsrc * rsrc;
  tmpb = -1.0 * sn2 * tmp;
  tmpc = -1.0 * sn2 * npix * readnoise * readnoise;

  /*
   * Calculate exposure time needed to achieve the requested SNR
   */

  texp = (-1.0*tmpb + sqrt(tmpb*tmpb - 4.0 * tmpa * tmpc)) / (2.0 * tmpa);
  noise = sqrt(tmp*texp + npix * readnoise * readnoise);

  /*
   * Print out results
   */

  printf("\nFor this observation the input parameters are:\n");
  printf(" Input flux            = %g erg/s/cm^2/Ang\n",srcflux);
  printf("                       = %g AB mags\n",srcmag);
  printf(" Input spec. intensity = %g erg/s/cm^2/Ang/arcsec^2\n",
	 srcflux/srcarea);
  printf("                       = %g AB mags/arcsec^2\n",
	 srcmu);
  printf(" Source count rate     = %g e-/sec\n",rsrc);
  printf(" Sky count rate        = %g e-/sec\n",rsky);
  printf("    Sky                = %g erg/sec/cm^2/Ang/arcsec^2\n",FSKY);
  printf("    QE                 = %f\n",QE);
  printf("    Filter bandwidth   = %f Ang\n",BW);
  printf("    Pixel size         = %f arcsec\n",pixsize);
  printf("    N_pix in source    = %f\n",npix);
  printf(" Thermal count rate    = %f phot/sec\n",rtherm);
  printf(" Dark current rate     = %f e-/sec\n",rdc);
  printf(" Read noise            = %f e-\n",readnoise);
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
  printf("  read noise        = %g cts\n",readnoise*sqrt(npix));
  printf(" Total noise        = %g cts\n\n",noise);
  printf("Which means you need an exposure time of %6.0f sec to get\n",texp);
  printf(" a SNR of %6.2f\n",sn);

  return 0;
}
