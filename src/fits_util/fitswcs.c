/*.......................................................................
 *
 * fitswcs.c
 *
 * A series of functions dealing with conversions between WCS info
 *  in fits headers.
 *
 * Functions included in this library:
 * -----------------------------------
 *  new_Imaxis
 *  del_Imaxis
 *  print_wcs_info
 *  read_wcs_info
 *  update_pixel_scale
 *  update_rotation
 *  update_ref_pixel
 *  wcs_init_guess
 *  rscale_to_cdmatrix
 *  cdmatrix_to_rscale
 *  ccdxy_to_deg
 *  radec_to_ccdxy
 *  darcsec_to_ccdxy
 *
 *
 *-----------------------------------------------------------------------
 * Revision history:
 * -----------------
 *  2007Jun27 - First split off from fitsxy2wcs.c. Chris Fassnacht (CDF)
 *  2007Jun29 - CDF - Split cdmatrix_to_rscale function from read_wcs_info.
 *  2007Jul08 - CDF - Added rscale_to_cdmatrix and invert_cdmatrix functions.
 *  2007Jul09 - CDF - Added update_pixel_scale and update_rotation functions.
 *  2007Jul11 - CDF - Added update_ref_pixel to update CRPIX and CRVAL.
 *   Improved handling of non-WCS info in read_wcs_info.
 *  2007Jul17 - CDF - Fixed a bug in update_ref_pixel
 *  2008Jul08 - CDF - Added new darcsec_to_ccdxy and radec_to_ccdxy functions
 *  2008Jul10 - CDF - Moved new_Imaxis and del_Imaxis functions from fitsim.h
 *   Moved functionality of old axis_info function in fitsim.h into 
 *    read_wcs_info function.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "coords.h"
#include "libfits.h"
#include "fitswcs.h"

/*.......................................................................
 *
 * Function new_Imaxis
 *
 * Allocates memory for a Imaxis structure array.
 *
 * Input:  int size            size of array
 *
 * Output: Imaxis *newimaxis   new setup array
 *
 */

Imaxis *new_Imaxis(int size)
{
  Imaxis *newimaxis;

  newimaxis = (Imaxis *) malloc((sizeof(Imaxis) * size));
  if(!newimaxis) {
    fprintf(stderr,"Insufficient memory for Imaxis array.\n");
    return NULL;
  }

  return newimaxis;
}

/*.......................................................................
 *
 * Function del_Imaxis
 *
 * Frees memory associated with Imaxis array
 *
 * Input:  Imaxis *imaxis      array to be freed
 *
 * Output: NULL
 *
 */

Imaxis *del_Imaxis(Imaxis *imaxis)
{
  if(imaxis)
    free(imaxis);

  return NULL;
}

/*.......................................................................
 *
 * Function print_wcs_info
 *
 * Prints out information contained in the WCS keywords
 *
 * Inputs:
 *   WCSinfo wcs               WCS information
 *
 * Outputs: (none)
 *
 */

void print_wcs_info(WCSinfo wcs)
{
  int i;   /* Looping variable */

  printf("WCS keywords\n");
  printf("------------\n");
  printf("# n      CDn_1        CDn_2        CDELTn   CRPIXn    CRVALn\n");
  printf("#--- ------------- ------------- --------- -------- ------------\n");
  for(i=0; i<2; i++) {
    printf("# %d  %+12e %+12e %+8f %8.2f  %f\n",
	   i+1,wcs.cdcol1[i],wcs.cdcol2[i],3600.0*wcs.cdelt[i],
	   wcs.crpix[i],wcs.crval[i]);
  }
}

/*.......................................................................
 *
 * Function print_rscale
 *
 * Prints rotation and pixel scale from image.
 *
 * Inputs:
 *  WCSinfo wcs                WCS information
 *
 * Output: (none)
 *  
 */

void print_rscale(WCSinfo wcs)
{
  printf("Pixel scales (arcsec/pix):    %f %f\n",
	 wcs.pixscale[0]*3600.0,wcs.pixscale[1]*3600.0);
  printf("Image rotations (N through E): %f %f\n",
	 180.0*wcs.rot[0]/PI,180.0*wcs.rot[1]/PI);
  printf("Average image rotation (N through E): %f\n",
	 180.0*wcs.rotation/PI);
}

/*.......................................................................,
 *
 * Function read_wcs_info
 *
 * Reads WCS information from the header keywords of a fits file.
 *
 * Inputs: Fits *fits          Input fits file
 *         WCSinfo wcs         WCS information
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int read_wcs_info(Fits *fits, WCSinfo *wcs)
{
  enum {
    CD1_1,
    CD1_2,
    CD2_1,
    CD2_2,
    CRVAL1,
    CRVAL2,
    CRPIX1,
    CRPIX2,
    RA,
    DEC
  };

  int i;              /* Looping variable */
  int count=0;        /* Keeps track of valid WCS inside a loop */
  int cdmatx=0;       /* Flag set to 1 when CD matrix is found */
  int validaxes=0;    /* Flag set to 1 if axes are successfully read */
  int raaxis=0;       /* Axis containing RA info */
  int decaxis=0;      /* Axis containing Dec info */
  Skypos crvalpos;    /* Temporary container for CRVAL in hms format */
  Phdu *phdu;         /* Temporary PHDU container for phdu */
  Hdu *hdu;           /* Temporary HDU container for phdu */
  Imaxis *axis;       /* Temporary image axis container */
  Imaxis *axptr;      /* Pointer for navigating wcs->axisinfo */
  Fitkey key;         /* Used for reading FITS header */
  Fitkey keys[]={
    {"CD1_1", 0, CD1_1, DAT_DBL, NULL, NULL},
    {"CD1_2", 0, CD1_2, DAT_DBL, NULL, NULL},
    {"CD2_1", 0, CD2_1, DAT_DBL, NULL, NULL},
    {"CD2_2", 0, CD2_2, DAT_DBL, NULL, NULL},
    {"CRVAL1", 0, CRVAL1, DAT_DBL, NULL, NULL},
    {"CRVAL2", 0, CRVAL2, DAT_DBL, NULL, NULL},
    {"CRPIX1", 0, CRPIX1, DAT_DBL, NULL, NULL},
    {"CRPIX2", 0, CRPIX2, DAT_DBL, NULL, NULL},
    /* {"RA", 0, RA, DAT_STR, NULL, NULL},
       {"DEC", 0, RA, DAT_STR, NULL, NULL} */
  };

  /*
   * Initialize the validwcs flag
   */

  wcs->validwcs = 0;

  /*
   * Set the appropriate pointer locations
   */

  phdu = (Phdu *) fits->hdu;
  hdu = fits->hdu;
  new_hline(hdu,0);

  /*
   * Make sure that the image is at least two-dimensional
   */

  if(phdu->naxis < 2) {
    fprintf(stderr,"ERROR: read_wcs_info.  Fits file has dimension < 2\n");
    return 1;
  }

  /*
   * Allocate the axisinfo array
   */

  if(!(wcs->axisinfo = new_Imaxis(phdu->naxis))) {
    fprintf(stderr,"ERROR: read_wcs_info.  \n");
    return 1;
  }

  /*
   * Initialize the CD matrix
   */

  for(i=0; i<2; i++) {
    wcs->cdcol1[i] = 0.0;
    wcs->cdcol2[i] = 0.0;
  }

  /*
   * Read the CD matrix
   */

  while(next_key(fits,hdu,keys,sizeof(keys)/sizeof(Fitkey),EOH_SEEK,&key) 
	== KEY_FOUND) {
    switch(key.keyid) {
    case CD1_1:
      wcs->cdcol1[0] = KEYDBL(key);
      cdmatx = 1;
      break;
    case CD1_2:
      wcs->cdcol2[0] = KEYDBL(key);
      cdmatx = 1;
      break;
    case CD2_1:
      wcs->cdcol1[1] = KEYDBL(key);
      cdmatx = 1;
      break;
    case CD2_2:
      wcs->cdcol2[1] = KEYDBL(key);
      cdmatx = 1;
      break;
    case CRVAL1:
      wcs->crval[0] = KEYDBL(key);
      break;
    case CRVAL2:
      wcs->crval[1] = KEYDBL(key);
      break;
    case CRPIX1:
      wcs->crpix[0] = KEYDBL(key);
      break;
    case CRPIX2:
      wcs->crpix[1] = KEYDBL(key);
      break;
    default:
      break;
    }
  }

  /*
   * Do calculations if CDm_n headers are found
   */

  if(cdmatx) {

    /*
     * Find rotation and pixel scale implied by CDm_n values
     */

    cdmatrix_to_rscale(wcs);
  }

  /*
   * Now try to read axis information from header - the reason for
   *  doing this is to check that the axis types are RA* and DEC*.
   *  If not, treat this as pixel values.
   * Also, set the cdelt values if the cd matrix has not been set
   *  and if the axis types are RA* and DEC*
   */

  for(i=0,axptr=wcs->axisinfo; i< phdu->naxis; i++,axptr++) {

    /*
     * Check first to see if axis information can even be read
     */

    if(!(axis=get_axis(phdu,i+1))) {
      printf("  read_wcs_info:  ");
      printf("Could not read axis %d info from header.\n",i+1);
      validaxes = 0;
      cdmatx = 0;
    }

    /*
     * If it can, then check to see if the ctype can be read.
     */

    else if(!axis->ctype) {
      printf("  read_wcs_info:  No CTYPE for axis %d\n",i+1);
      validaxes = 0;
      cdmatx = 0;
    }

    /*
     * If ctype can be read, put into the wcs->axisinfo container
     */

    else {
      validaxes = 1;
      *axptr = *axis;
    }
  }

  /*
   * Now load CDELT info if the CD matrix has not already been read
   */

  if(validaxes && cdmatx == 0) {
    if((raaxis = find_axis(phdu,"RA",2,1))>0) {
      axptr = wcs->axisinfo + raaxis - 1;
      if((wcs->cdelt[0] = axptr->cdelt)>0)
	cdmatx = 1;
      else {
	fprintf(stderr,"ERROR. No CD matrix and no RA CDELT info\n");
	cdmatx = 0;
      }
    }
    if((decaxis = find_axis(phdu,"DEC",3,1))>0) {
      axptr = wcs->axisinfo + decaxis - 1;
      if((wcs->cdelt[1] = axptr->cdelt)>0)
	cdmatx = 1;
      else {
	fprintf(stderr,"ERROR. No CD matrix and no DEC CDELT info\n");
	cdmatx = 0;
      }
    }
  }

  /*
   * If cdmatx is still 0, it means that neither the CD matrix nor the
   * CDELT vectors are present.
   */

  if(!cdmatx) {
    fprintf(stderr," read_wcs_info: No WCS info in header\n");
    return 1;
  }
  else {
    wcs->validwcs = 1;
    wcs->pixscale[0] = -wcs->cdelt[0];
    wcs->pixscale[1] = wcs->cdelt[1];
    deg2spos(wcs->crval[0],wcs->crval[1],&crvalpos);
    sprintf(wcs->fitsra,"%02d:%02d:%07.4f",
	    crvalpos.hr,crvalpos.min,crvalpos.sec);
    sprintf(wcs->fitsdec,"%+03d:%02d:%07.4f",
	    crvalpos.deg,crvalpos.amin,crvalpos.asec);
  }

  /*
   * Clean up and return
   */

  return 0;
}

/*.......................................................................
 *
 * Function update_pixel_scale
 *
 * Updates the pixel scale and CD matrix
 *
 * Inputs:
 *  WCSinfo *wcs               WCS information
 *
 * Output: (none)
 *  
 */

void update_pixel_scale(WCSinfo *wcs)
{
  char line[MAXC];    /* Generic string for reading input */

  printf("Updating pixel scale\n");
  printf("--------------------\n");
  printf(" Current values:\n");
  printf("   x pixel scale: %8.4f arcsec/pix\n",wcs->pixscale[0]*3600.0);
  printf("   y pixel scale: %8.4f arcsec/pix\n",wcs->pixscale[1]*3600.0);

  /*
   * Do x pixel scale first
   */

  printf("\nEnter the new x pixel scale in arcsec/pix: [%8.4f] ",
	 wcs->pixscale[0]*3600.0);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%lf",&wcs->pixscale[0]) !=1 && 
	  wcs->pixscale[0]<0.0) {
      fprintf(stderr,"ERROR: Bad input format for pixel scale. ");
      fprintf(stderr,"Enter pixel scale again: ");
      fgets(line,MAXC,stdin);
    }
    wcs->pixscale[0] /= 3600.0;
  }

  /*
   * Now do y pixel scale, assuming as a default that it is the
   *  same as the (possibly updated) x pixel scale
   */

  wcs->pixscale[1] = wcs->pixscale[0];
  printf("\nEnter the new y pixel scale in arcsec/pix: [%8.4f] ",
	 wcs->pixscale[1]*3600.0);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%lf",&wcs->pixscale[1]) !=1 && 
	  wcs->pixscale[1]<0.0) {
      fprintf(stderr,"ERROR: Bad input format for pixel scale. ");
      fprintf(stderr,"Enter pixel scale again: ");
      fgets(line,MAXC,stdin);
    }
    wcs->pixscale[1] /= 3600.0;
  }

  /*
   * Finally, update the CD matrix to reflect new pixel scale
   */

  wcs->cdelt[0] = -wcs->pixscale[0];
  wcs->cdelt[1] = wcs->pixscale[1];
  rscale_to_cdmatrix(wcs);
}

/*.......................................................................
 *
 * Function update_rotation
 *
 * Inputs:
 *  WCSinfo *wcs               WCS information
 *
 * Output: (none)
 *  
 */

void update_rotation(WCSinfo *wcs)
{
  char line[MAXC];    /* Generic string for reading input */

  printf("Updating image position angle\n");
  printf("-----------------------------\n");
  printf(" Current value: %+7.2f\n",wcs->rotation * 180.0 / PI);

  printf("\nEnter the new position angle in degrees: [%+7.2f] ",
	 wcs->rotation*180.0/PI);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%lf",&wcs->rotation) !=1) {
      fprintf(stderr,"ERROR: Bad input format for position angle. ");
      fprintf(stderr,"Enter position again: ");
      fgets(line,MAXC,stdin);
    }
    wcs->rotation *= (PI / 180.0);
  }

  /*
   * Finally, update the CD matrix to reflect new rotation angle
   */

  wcs->rot[0] = wcs->rotation;
  wcs->rot[1] = wcs->rotation;
  rscale_to_cdmatrix(wcs);
  print_wcs_info(*wcs);
  printf("\n");
  print_rscale(*wcs);
}

/*.......................................................................
 *
 * Function update_ref_pixel
 *
 * Updates the reference pixel -- crpix*, crval*, or both
 *
 * Inputs:
 *  WCSinfo *wcs               WCS information
 *
 * Output: (none)
 *  
 */

void update_ref_pixel(WCSinfo *wcs)
{
  int i;              /* Looping variable */
  char line[MAXC];    /* Generic string for reading input */
  Pos tmppos;         /* Temporary storage for updated CRVAL */
  Skypos skypos;      /* CRVAL in hms format*/

  printf("Updating reference pixel\n");
  printf("------------------------\n");
  printf(" Current values:\n");
  printf("      CRPIX               CRVAL\n");
  printf("  -------------  ---------------------------\n");
  deg2spos(wcs->crval[0],wcs->crval[1],&skypos);
  printf("  %6.2f %6.2f  ",wcs->crpix[0],wcs->crpix[1]);
  print_skypos(stdout,skypos);
  printf("\n");

  /*
   * Do crpix first
   */

  printf("\nEnter the new x and y values for the reference pixel (CRPIX)\n");
  printf("  [%6.2f %6.2f] ",wcs->crpix[0],wcs->crpix[1]);
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    while(sscanf(line,"%lf %lf",&wcs->crpix[0],&wcs->crpix[1]) !=2) {
      fprintf(stderr,"ERROR: Bad input format for CRPIX ");
      fprintf(stderr,"Enter both x and y, separated by a space: ");
      fgets(line,MAXC,stdin);
    }
  }

  /*
   * Now do crval
   */

  wcs->pixscale[1] = wcs->pixscale[0];
  printf("\nEnter the new CRVAL as RA and Dec: [");
  print_skypos(stdout,skypos);
  printf("] ");
  fgets(line,MAXC,stdin);
  if(line[0] != '\n') {
    tmppos = read_skypos();
    wcs->crval[0] = tmppos.x;
    wcs->crval[1] = tmppos.y;
  }
}

/*.......................................................................
 *
 * Function wcs_init_guess
 *
 * Sets the initial guess for the WCS info of a fits file.  In the simplest
 *  case, the fits file will have a valid and fairly accurate CD matrix
 *  that can be accessed through read_wcs_info.  If, however, there is
 *  no WCS info in the fits header, the user can enter his or her own
 *  initial guess for the WCS info.
 *
 */

WCSinfo wcs_init_guess(Fits *fits)
{
  int i;            /* Looping variable */
  int trytargpos=0; /* Flag set to 1 to try search for target RA,Dec keywords */
  int getpos=0;     /* Flag set to 1 to get RA,Dec interactively */
  char line[MAXC];  /* General string for reading input */
  Pos tmppos;       /* Temporary holder for RA and Dec as decimal degrees */
  Phdu *phdu;       /* Temporary PHDU container for phdu */
  WCSinfo wcsinit;  /* Initial guess for WCS information */

  /*
   * First try to read in WCS info from the fits file CD matrix, etc.
   * If this is found and approved, return the values from the fits header
   */

  printf("\nChecking for WCS information in fits header.\n");
  if(read_wcs_info(fits,&wcsinit)) {
    printf("\nNo WCS information found in standard place\n");
  } else {
    printf("\nFound WCS information in fits header:\n\n");
    print_wcs_info(wcsinit);
    printf("\n");
    print_rscale(wcsinit);
#if 0
    printf("\nUse this WCS information in fits header as initial guess? [y] ");
    fgets(line,MAXC,stdin);
    if(line[0] !=  'n' && line[0] != 'N')
      return wcsinit;
#endif
    return wcsinit;
  }

  /*
   * If WCS information is not found in fits header or is bad, then have
   *  to get the information another way.  Start by setting the CRPIXn
   *  values as roughly the center of the image
   */

  phdu = (Phdu *) fits->hdu;
  for(i=0; i<2; i++)
    wcsinit.crpix[i] = phdu->dims[i] / 2.0;

  /*
   * At this point, program should check the fits header for RA and DEC
   *  keywords, which usually reflect the rough telescope pointing.
   * Add this in later
   */

  getpos = 1;

  /* 
   * Get RA and Dec if necessary
   */

  if(getpos) {
    printf("\n\nEnter the rough RA and Dec for the center of the image.\n");
    printf("-------------------------------------------------------\n");
    tmppos = read_skypos();
    wcsinit.crval[0] = tmppos.x;
    wcsinit.crval[1] = tmppos.y;
  }

  /*
   * Set default values for PA and pixel scale
   */

  wcsinit.pixscale[0] = wcsinit.pixscale[1] = 0.1 / 3600.0;
  wcsinit.rotation = 0.0;
  printf("\nSetting default values for pixel scale and PA:\n");
  printf(" Pixel scale (x and y): 0.1 arcsec/pix\n");
  printf(" Position angle (N->E): 0.0 degrees\n");

  /*
   * Get rough pixel scale
   */

  printf("\nSetting the rough pixel scale for the data in arcsec/pix\n");
  update_pixel_scale(&wcsinit);

  /*
   * Get rough PA
   */

  printf("\nEnter the rough PA for the image in degrees east of north\n");
  update_rotation(&wcsinit);

  /*
   * Return values
   */

  printf("\nConverted information into a CD matrix:\n");
  print_wcs_info(wcsinit);
  print_rscale(wcsinit);
  return(wcsinit);
}

/*.......................................................................
 *
 * Function rscale_to_cdmatrix
 *
 * Converts a pixel scale and rotation into a CD matrix.  Does this
 *  via the following formulas (taken from 
 *   http://stsdas.stsci.edu/documents/SUG/UG_21.html):
 *
 * CD1_1 = cdcol1[0] = -pixscale_1 * cos(theta_1)
 * CD1_2 = cdcol2[0] = -pixscale_2 * sin(theta_2)
 * CD2_1 = cdcol1[1] = -pixscale_1 * sin(theta_1)
 * CD2_2 = cdcol2[1] =  pixscale_2 * cos(theta_2)
 *
 * Inputs:
 *  WCSinfo *wcs               Container for pixel scale, rotation, and
 *                              output CD matrix.  This structure gets
 *                              modified by the function.
 *
 * Output: (none)
 */

void rscale_to_cdmatrix(WCSinfo *wcs)
{
  double rot0,rot1;   /* Rotation angles in radians */

  rot0 = wcs->rot[0];
  rot1 = wcs->rot[1];

  wcs->cdcol1[0] = -1.0 * wcs->pixscale[0] * cos(rot0);
  wcs->cdcol1[1] = -1.0 * wcs->pixscale[0] * sin(rot0);
  wcs->cdcol2[0] = -1.0 * wcs->pixscale[1] * sin(rot1);
  wcs->cdcol2[1] = wcs->pixscale[1] * cos(rot1);
}

/*.......................................................................
 *
 * Function cdmatrix_to_rscale
 *
 * Converts a CD matrix into apixel scale and rotation.
 *
 * Inputs:
 *  WCSinfo *wcs               Container for pixel scale, rotation, and
 *                              CD matrix.  The pixel scale and rotation get
 *                              modified by the function.
 *
 * Output: (none)
 */

void cdmatrix_to_rscale(WCSinfo *wcs)
{
  int cdsgn;          /* Sign of CDm_n determinant */
  double det;         /* Determinant of CDm_n matrix */

  det = wcs->cdcol1[0] * wcs->cdcol2[1] - wcs->cdcol2[0] * wcs->cdcol1[1];
  if(det<0) {
    cdsgn = -1;
  }
  else {
    cdsgn = 1;
    fprintf(stderr,"  WARNING: read_wcs_info - Astrometry is for a ");
    fprintf(stderr,"right-handed coordinate system.\n");
  }
  if(wcs->cdcol2[0] == 0 && wcs->cdcol1[1] == 0) {
    wcs->rot[0] = 0.0;
    wcs->rot[1] = 0.0;
  }
  else {
    wcs->rot[0] = atan2(cdsgn * wcs->cdcol1[1],cdsgn * wcs->cdcol1[0]);
    wcs->rot[1] = atan2(-1.0 * wcs->cdcol2[0], wcs->cdcol2[1]);
  }
  if(fabs(wcs->rot[0] - wcs->rot[1]) * 180.0/PI < 2.0)
    wcs->rotation = (wcs->rot[0] + wcs->rot[1])/2.0;
  else
    wcs->rotation = wcs->rot[0];

  /*
   * Now convert CDn_m headers into CDELTn values.  For now assume
   *  CDELT1 = sqrt(CD1_1^2 + CD2_1^2) and CDELT2 = sqrt(CD1_2^2 + CD2_2^2).
   */

  wcs->cdelt[0] = -sqrt(wcs->cdcol1[0]*wcs->cdcol1[0] +
			wcs->cdcol1[1]*wcs->cdcol1[1]);
  wcs->cdelt[1] = sqrt(wcs->cdcol2[0]*wcs->cdcol2[0] + 
		       wcs->cdcol2[1]*wcs->cdcol2[1]);
}

/*.......................................................................
 *
 * Function invert_cdmatrix
 *
 * Inverts the CD matrix, assumed to be 2x2 for now.
 *
 * Inputs:
 *  WCSinfo wcs                structure containing input CD matrix
 *  WCSinfo *wcsinv            output structure containing inverted matrix
 *
 * Output:
 *  int (0 or 1)               0 ==> success, 1 ==> error
 *
 */

int invert_cdmatrix(WCSinfo wcs, WCSinfo *wcsinv)
{
  double det;    /* Determinant of input CD matrix */

  /*
   * Compute determinant
   */

  det = wcs.cdcol1[0] * wcs.cdcol2[1] - wcs.cdcol2[0] * wcs.cdcol1[1];
  if(det == 0) {
    fprintf(stderr,"ERROR: invert_cdmatrix.  CD matrix is singular");
    return 1;
  }

  /*
   * Invert matrix
   */

  wcsinv->cdcol1[0] = wcs.cdcol2[1] / det;
  wcsinv->cdcol2[0] = -1.0 * wcs.cdcol2[0] / det;
  wcsinv->cdcol1[1] = -1.0 * wcs.cdcol1[1] / det;
  wcsinv->cdcol2[1] = wcs.cdcol1[0] / det;

  return 0;
}

/*.......................................................................
 *
 * Function ccdxy_to_deg
 *
 * Given a WCS solution (in the form of a CD matrix, CRPIX, CRVAL),
 *  converts a set of (x,y) positions on a CCD into (alpha,delta) sky positions.
 *
 * The proper formula should be:
 *
 *  ( alpha )   (alpha_0)   (dalpha/dx  dalpha/dy)   (x - x_0)
 *  (       ) = (       ) + (                    ) * (       )
 *  ( delta )   (delta_0)   (ddelta/dx  ddelta/dy)   (y - y_0)
 *
 *  where the vectors and matrix map as:
 *     (alpha_0 delta_0)  are the CRVAL* values (in decimal degrees)
 *     (x_0 y_0)          are the CRPIX* values (in pixels)
 *     the matrix         contains the CDm_n values (in decimal degrees per pix)
 *
 * Inputs:
 *  Pos *xypos                 (x,y) positions on the CCD
 *  int npos                   number of (x,y) positions to convert
 *  WCSinfo                    WCS information
 *
 * Output:
 *  Pos *adpos                 (alpha,delta) positions on the sky in decimal
 *                              degrees
 */

Pos *ccdxy_to_deg(Pos *xypos, int npos, WCSinfo wcs)
{
  int i;             /* Looping variable */
  int no_error=1;    /* Flag set to 0 on error */
  Pos cent;          /* CRVAL values, expressed as a position structure */
  Pos *dpos=NULL;    /* Offset from (x_0,y_0) in arcsec */
  Pos *adpos=NULL;   /* Output array of sky positions */
  Pos *xptr,*aptr;   /* Pointers to the position arrays */
  Pos xydiff;        /* Offset of (x,y) from (x_0,y_0) */

  /*
   * Allocate memory for the output array 
   */

  if(!(dpos = new_pos(npos,1))) {
    fprintf(stderr,"ERROR: ccdxy_to_deg.\n");
    return NULL;
  }

  /*
   * Loop through list, converting (x,y) to (d_alpha,d_delta) in arcsec
   */

  for(i=0,xptr=xypos,aptr=dpos; i<npos; i++,xptr++,aptr++) {
    xydiff.x = xptr->x - wcs.crpix[0];
    xydiff.y = xptr->y - wcs.crpix[1];
    aptr->x = wcs.cdcol1[0] * xydiff.x + wcs.cdcol2[0] * xydiff.y;
    aptr->y = wcs.cdcol1[1] * xydiff.x + wcs.cdcol2[1] * xydiff.y;
  }

  /*
   * Convert the offsets into (alpha,delta) in decimal degrees
   */

  cent.x = wcs.crval[0];
  cent.y = wcs.crval[1];
  if(!(adpos = ddeg2deg(cent,dpos,npos)))
    no_error = 0;

  /*
   * Clean up and exit
   */

  dpos = del_pos(dpos);
  if(no_error)
    return adpos;
  else {
    fprintf(stderr,"ERROR: ccdxy_to_deg.\n");
    return NULL;
  }
}

/*.......................................................................
 *
 * Function radec_to_ccdxy
 *
 * Given a list of (RA,Dec) positions in a Secat array and the WCS information
 *  from a fits header, converts the (RA,Dec) positions into (x,y) positions
 *  on the CCD.
 *
 * Inputs:
 *  Fits *fits                 fits file container
 *  Secat *secat               container for (RA,Dec) position list. NB: this
 *                              list is modified to contain the calculated
 *                              (x,y) positions in its x and y members.
 *  int ncat                   number of positions in list
 *
 * Output:
 *  int (SUCCESS or ERROR)
 *
 */

int radec_to_ccdxy(Fits *fits, Secat *secat, int ncat)
{
  int no_error=1;    /* Flag set to 0 on error */
  Secat *sptr;       /* Pointer to navigate secat */
  Secat crvalpos;    /* Secat container for CRVALn values */
  WCSinfo wcs;       /* Container for WCS information */

  /*
   * Read the WCS information from the fits header
   */

  if(read_wcs_info(fits,&wcs))
    no_error = 0;
  if(no_error) {
    printf("\n");
    print_wcs_info(wcs);
    printf("\n");
    crvalpos.alpha = wcs.crval[0];
    crvalpos.delta = wcs.crval[1];
  }

  /*
   * Compute offsets in arcsec between each of the catalog positions and the
   *  position designated by the CRVAL values
   */

  if(no_error)
    if(secat2offset(secat,ncat,DEGS,crvalpos,DEGS) == ERROR)
      no_error = 0;
  
  /*
   * Now convert the offsets in arcsec into pixel values and print out results
   */

  if(no_error)
    if(darcsec_to_ccdxy(secat,ncat,wcs) == ERROR)
      no_error = 0;

  /*
   * Return
   */

  if(no_error)
    return SUCCESS;
  else {
    fprintf(stderr,"ERROR: radec_to_ccdxy\n");
    return ERROR;
  }
}

/*.......................................................................
 *
 * Function darcsec_to_ccdxy
 *
 * Converts a list of offsets in arcsec to (x,y) positions on a CCD
 *  given the CCD WCS information.
 *
 * Inputs:
 *  Secat *secat               container for offsets (stored in dx and dy)
 *                              Note that the computed (x,y) positions will
 *                              be stored in the x and y members of the
 *                              Secat structure
 *  int ncat                   number of offsets in catalog
 *  WCSinfo wcs                CCD WCS information
 *
 * Output:
 *  SUCCESS or ERROR
 */

int darcsec_to_ccdxy(Secat *secat, int ncat, WCSinfo wcs)
{
  int i;           /* Looping variable */
  int no_error=1;  /* Flag set to 0 on error */
  Secat *sptr;     /* Pointer to navigate secat */
  WCSinfo wcsinv;  /* Container for inverse of CD matrix */

  /*
   * Invert the CD matrix in preparation for converting offsets to (x,y) value
   */

  if(invert_cdmatrix(wcs,&wcsinv))
    no_error = 0;

  /*
   * Loop through secat, calculating (x,y) positions and storing them in
   *  secat->x and secat->y
   */

  if(no_error) {
    for(i=0,sptr=secat; i<ncat; i++,sptr++) {
      sptr->x = wcs.crpix[0] + (sptr->dx * wcsinv.cdcol1[0] + 
				sptr->dy * wcsinv.cdcol2[0])/3600.0;
      sptr->y = wcs.crpix[1] + (sptr->dx * wcsinv.cdcol1[1] + 
				sptr->dy * wcsinv.cdcol2[1])/3600.0;
    }
  }

  /*
   * Return
   */

  if(no_error)
    return SUCCESS;
  else {
    fprintf(stderr,"ERROR: darcsec_to_ccdxy\n");
    return ERROR;
  }
}
