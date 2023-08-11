#ifndef cosmo_h
#define cosmo_h

#define C 3.0e10          /* Speed of light in cm/sec */
#define G 6.67e-8         /* Gravitational constant in cgs units */
#define H0 3.23e-18       /* Because H_0=100 h km/sec/Mpc = 3.2e-18 h sec^-1 */
#define KMSMPC2S 3.23e-20 /* Conversion from km/sec/Mpc to sec^{-1} */
#define MPC2CM 3.1e24     /* Conversion from Mpc to cm */
#define KM2CM 1.0e5       /* Conversion from km to cm */
#define RAD2ASEC 206265.0 /* To convert radians to arcsec */
#define YR2SEC 3.156e7    /* To convert years to seconds */
#define MSUN 1.99e33      /* Solar mass in g */
#define LSUN 3.90e33      /* Solar luminosity in erg/sec */
#define MSUN_U 5.51       /* Solar magnitude in the U band */
#define MSUN_B 5.41       /* Solar magnitude in the B band */
#define MSUN_V 4.79       /* Solar magnitude in the V band */
#define MSUN_R 4.27       /* Solar magnitude in the R band */
#define MSUN_I 4.01       /* Solar magnitude in the I band */

/*.......................................................................
 *
 * Enumeration for band
 *
 */

enum {
  OTHER, 
  U, 
  B,
  V,
  R,
  I
};

/*.......................................................................
 *
 * Structure definitions
 *
 */

typedef struct {
  double omega_m;    /* Omega_matter */
  double omega_de;   /* Omega_dark-energy */
  double w;          /* Dark energy "equation of state parameter */
  double h;          /* Hubble Constant */
} Cosmo;

typedef struct {
  double z;       /* Redshift  */
  double zerr;    /* Error on redshift */
  double a;       /* The scale factor at z, assuming a(t_0) = 1 */
  double ez;      /* Peeble's E(z) */
  double d_c;     /* Comoving distance to z (line-of sight) */
  double d_m;     /* Transverse comoving (or proper motion) distance to z */
  double d_a;     /* Angular diameter distance to z */
  double d_l;     /* Luminosity distance to z */
  double DM;      /* Distance modulus to z */
  double t_l;     /* Lookback time to z */
  double age;     /* Age of Universe at z */
  double hz;      /* Hubble Constant at z */
} Cosdist;

/*.......................................................................
 *
 * Function declarations
 *
 */

int get_valerr(double *val, double *dval, char *valdef);
void get_cosmo(Cosmo *cosmo);
Cosdist calc_cosdist(double z1, double z2, Cosmo cosmo);
double rn_generic(double z, Cosmo cosmo, double m_n, double n);
double rn_sis(double z, Cosmo cosmo, double sigma, double n);
double rn_carlberg(double z, Cosmo cosmo, double sigma, double n);
double rn_isothermal(double z, Cosmo cosmo, double sigma, double n, double K);
double mn_sis(double z, Cosmo cosmo, double sigma, double n);
double peebles_E(double z, Cosmo cosmo);
double mag_to_lum(double absmag);

#endif
