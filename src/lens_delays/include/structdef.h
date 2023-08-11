#ifndef structdef_h
#define structdef_h

#define MAXC 1000
#define PI 3.141592653589793

/*.......................................................................
 *
 * Enumerations
 *
 */

enum {
  SUCCESS,
  ERROR
}; /* Enumeration for returns of integer functions */

enum {
  DEGS,
  HMS
}; /* Enumeration for format of sky position */


/*.......................................................................
 *
 * Structure definitions
 *
 */

typedef struct {
  double x;
  double y;
  double xerr;
  double yerr;
  int flag;
  char label[MAXC];
} Pos;

typedef struct {
  int hr;
  int min;
  double sec;
  int deg;
  int amin;
  double asec;
  char label[MAXC];
} Skypos;

typedef struct {
  double x;          /* First coord */
  double y;          /* Second coord */
  double z;          /* Third coord */
  int dataflag;      /* Flag */
  char text[MAXC];   /* Label */
} Datastruct;        /* General structure for 3D coords, a flag and a label */

typedef struct {
  int nval;          /* Number of values along the axis */
  float minval;      /* Minimum axis value */
  float maxval;      /* Maximum axis value */
} Axisinfo;          /* Structure for information about an axis in a grid */

typedef struct {
  float x;           /* x coordinate */
  float y;           /* y coordinate */
  char text[MAXC];   /* Label text */
  float size;        /* Size of marker or label text */
} Labinfo;

typedef struct {
  char name[MAXC];   /* Source name */
  int id;            /* ID number */
  double dx;         /* x displacement from lens system */
  double dy;         /* y displacement from lens system */
  double dpos;       /* Total displacement from lens system */
  double dpostheta;  /* For expressing dpos in (r,theta) format */
  Skypos skypos;     /* RA  and Dec of object, in HMS format */
  double alpha;      /* RA of object, in decimal degrees */
  double delta;      /* Dec of object, in decimal degrees */
  float u;           /* u-band magnitude */
  float g;           /* g-band magnitude */
  float r;           /* r-band magnitude */
  float i;           /* i-band magnitude */
  float z;           /* z-band magnitude */
  float zspec;       /* Redshift */
  float zspecerr;    /* Error on redshift */
  float class;       /* Stellarity class */
} SDSScat;

typedef struct {
  char name[MAXC];   /* Source name */
  int id;            /* ID number */
  double x;          /* First coord */
  double y;          /* Second coord */
  double dx;         /* x displacement (in arcsec) from lens system */
  double dy;         /* y displacement (in arcsec) from lens system */
  double dpos;       /* Total displacement from lens system */
  double dpostheta;  /* For expressing dpos in (r,theta) format */
  Skypos skypos;     /* RA  and Dec of object, in HMS format */
  double alpha;      /* RA of object, in decimal degrees */
  double delta;      /* Dec of object, in decimal degrees */
  float miso;        /* Isophotal magnitude */
  float misoerr;     /* Error on isophotal magnitude */
  float fiso;        /* Flux in isophotal area */
  float fisoerr;     /* Error on isophotal flux */
  float mtot;        /* Total magnitude by "auto" method */
  float merr;        /* Error on mtot */
  float fauto;       /* Flux from "automatic" aperture */
  float fautoerr;    /* Error on fauto */
  float r_kron;      /* Kron radius */
  float bkgd;        /* Background level at centroid */
  float thresh;      /* Threshold in counts */
  float muthresh;    /* Threshold surface brightness */
  float isoarea;     /* Area enclosed by the isophotal region */
  int naper;         /* Number of apertures (must be <= 100) */
  float maper[100];  /* Aperture magnitude(s) */
  float mapererr[100]; /* Error(s) on aperture magnitude(s) */
  float faper[100];  /* Aperture flux(es) */
  float fapererr[100]; /* Error(s) on aperture flux(es) */
  float ma1,ma2,ma3; /* Old-style aperture magnitudes. DELETE SOON */
  float ma1err,ma2err,ma3err; /* Old-style aperture magnitude errors */
  float fa1,fa2,fa3; /* Old-style aperture fluxes. DELETE SOON */
  float fa1err,fa2err,fa3err; /* Old-style aperture flux errors */
  float a_im;        /* Profile RMS along major axis */
  float b_im;        /* Profile RMS along major axis */
  float theta;       /* PA of major axis (math definition) */
  float fwhm;        /* FWHM of image */
  float r2;          /* Radius enclosing 20% of detected flux */
  float r5;          /* Radius enclosing half of detected flux */
  float r8;          /* Radius enclosing 80% of detected flux */
  float class;       /* Star/gal classifier (0.0 --> 1.0 (most starlike)) */
  int fitflag;       /* Flag returned by SExtractor */
  int nmatch;        /* Number of matches found */
  int matchcat[100]; /* Catalog in which the match was found */
  int matchid[1000]; /* Flag indicating a match has been found */
  int matchflag[10]; /* NOT USED ANYMORE*/
  float sep[1000];   /* Positional separation between catalogs */
  float zspec;       /* Redshift of source */
  float zspecerr;    /* Error on redshift */
} Secat;             /* Structure for SExtractor output info */

typedef struct {
  char colname[MAXC]; /* Column name */
  int colnum;         /* Column position */
  char nelem;         /* Number of elements (e.g., for mag_aper) */
  char datatype[3];   /* Data type of column (e.g. %f) */
} SHeader; /* SExtractor catalog header structure */

typedef struct {
  float par;       /* Model parameter */
  int varflag;     /* Flag set to 1 for a free parameter */
} Modparam;

typedef struct {
  int pottyp;
  Modparam xl;
  Modparam yl;
  Modparam b;
  Modparam gamma;
  Modparam theta_g;
} Lensparams;

typedef struct {
  float ph;
  float dphdx;
  float dphdy;
  float dphdxx;
  float dphdyy;
  float dphdxy;
} Poten;

typedef struct {
  float x;
  float y;
  float zf;
  float flux;
  Pos spos;
  float sflux;
} Posinfo;

typedef struct {
  Modparam x;
  Modparam y;
  Modparam flux;
  Modparam zf;
} Srcinfo;

typedef struct {
  Lensparams *lp;
  int nmod;
  Srcinfo *fsource;
  int nfsrc;
  int *var;
  int nvar;
  int nlinks;
  Posinfo *posobs;
  int nobs;
  Posinfo *posimmod;
  int nimmod;
  float pixscale;
  float poserr;
  float ferr;
  int chisqflag;
} Modinfo;

typedef struct {
  int idnum;
  char ionname[MAXC];
  float lambda;
  char label[MAXC];
  char comment[MAXC];
} Specinfo;

typedef struct {
  float day;
  float flux;
  float err;
  int match;
} Fluxrec; 

/*.......................................................................
 *
 * Function declarations
 *
 */

char *new_string(int size);
char *del_string(char *string);

float *new_array(int size1,int size2);
float *del_array(float *array);

float **new_float2(int n2);
float **del_float2(float **float2, int n2);

int *new_intarray(int size1,int size2);
int *del_intarray(int *array);

double *new_doubarray(int size);
double *del_doubarray(double *array);

double **new_doub2(int ndoub2);
double **del_doub2(double **doub2, int ndoub2);

Datastruct *new_datastruct(int size);
Datastruct *del_datastruct(Datastruct *datastruct);

Axisinfo *new_axisinfo(int size);
Axisinfo *del_axisinfo(Axisinfo *axisinfo);

Secat *new_secat(int size);
Secat *del_secat(Secat *secat);

SDSScat *new_sdsscat(int size);
SDSScat *del_sdsscat(SDSScat *sdsscat);

Lensparams *new_lp(int n);
Lensparams *del_lp(Lensparams *lp);

Modparam *new_modparam(int size);
Modparam *del_modparam(Modparam *modparam);

Pos *new_pos(int n1,int n2);
Pos *del_pos(Pos *pos);
float **pos2xy(Pos *pos, int npts, int *ngood);

Poten *new_poten(int n1, int n2);
void reinit_poten(Poten *potptr,int n1,int n2);
Poten *del_poten(Poten *potptr);

Posinfo *new_posinfo(int n1, int n2);
Posinfo *del_posinfo(Posinfo *posinfo);

Srcinfo *new_srcinfo(int n1, int n2);
Srcinfo *del_srcinfo(Srcinfo *srcinfo);

Modinfo *new_modinfo();
int modpars_to_modinfo(float *modpars, Modinfo *modinfo);
Modinfo *del_modinfo(Modinfo *info);

int modpars_to_struct(float *valptr, int *varptr, int *maxvar, float *mpptr,
		      float *linkvals);

Specinfo *new_specinfo(int n);
Specinfo *del_specinfo(Specinfo *specinfo);

Labinfo *new_labinfo(int n);
Labinfo *del_labcinfo(Labinfo *labinfo);

#endif
