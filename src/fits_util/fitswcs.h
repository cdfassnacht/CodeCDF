#ifndef fitswcs_h
#define fitswcs_h

#include "structdef.h"
#include "coords.h"
#include "libfits.h"

/*.......................................................................
 *
 * Structure definitions
 *
 */

typedef struct {
  int validwcs;        /* Flag set to 1 if the file has valid WCS info */
  Imaxis *axisinfo;    /* Key values for each axis */
  double cdcol1[2];    /* 1st column of CD matrix: (CD1_1 and CD2_1) */
  double cdcol2[2];    /* 2nd column of CD matrix: (CD1_2 and CD2_2) */
  double cdelt[2];     /* CDELT vector */
  double crpix[2];     /* CRPIX vector */
  double crval[2];     /* CRVAL vector */
  double pixscale[2];  /* Pixel scale along the two axes */
  double rot[2];       /* WCS axis rotation (individual) */
  double rotation;     /* WCS axis rotation (N->E) */
  double targ_ra;      /* Target RA (what telescope was pointing at) */
  double targ_dec;     /* Target Dec (what telescope was pointing at) */
  char fitsra[80];     /* RA of crval in hh:mm:ss.ss format */
  char fitsdec[80];    /* Dec of crval in dd:mm:ss.ss format */
} WCSinfo;

/*.......................................................................
 *
 * Function declarations
 *
 */

Imaxis *new_Imaxis(int size);
Imaxis *del_Imaxis(Imaxis *imaxis);
void print_wcs_info(WCSinfo wcs);
void print_rscale(WCSinfo wcs);
int read_wcs_info(Fits *fits, WCSinfo *wcs);
void update_pixel_scale(WCSinfo *wcs);
void update_rotation(WCSinfo *wcs);
void update_ref_pixel(WCSinfo *wcs);
WCSinfo wcs_init_guess(Fits *fits);
void rscale_to_cdmatrix(WCSinfo *wcs);
void cdmatrix_to_rscale(WCSinfo *wcs);
int invert_cdmatrix(WCSinfo wcs, WCSinfo *wcsinv);
Pos *ccdxy_to_deg(Pos *xypos, int npos, WCSinfo wcs);
int darcsec_to_ccdxy(Secat *secat, int ncat, WCSinfo wcs);
int radec_to_ccdxy(Fits *fits, Secat *secat, int ncat);

#endif
