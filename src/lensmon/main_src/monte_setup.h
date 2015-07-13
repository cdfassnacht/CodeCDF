#ifndef monte_setup_h
#define monte_setup_h

/*.......................................................................
 *
 * Enumeration for setup file parameters
 *
 */

enum {
  MCMU0,
  MCNMU,
  MCTAU0,
  MCNTAU,
  MCSIG,
  NIDEAL,
  NOBS,
  NOUT,
  FTOL,
  MCROOT,
  MCDEFAULT,
  MCSETUPERR
};

/*.......................................................................
 *
 * Structure definitions
 *
 */

typedef struct {
  float mu0[4];         /* Initial guesses for magnifications */
  int nmu;              /* Number of mag steps on either side of mu0 */
  float tau0[4];        /* Initial guess for delays */
  int ntau;             /* Number of delay steps on either side of tau0 */
  int nideal;           /* Number of points in ideal light curve */
  int nobs;             /* Number of points in observed light curves */
  float mcsig[4];       /* Sigma for width of gaussian deviate distribution */
  int nout;             /* Number of output curves */
  float ftol;           /* Tolerance for comparision of float variables */
  char root[MAXC];      /* Root name for files */
} MC_Setup;

/*.......................................................................
 *
 * Function declarations
 *
 */

MC_Setup *new_MC_Setup(int size);
MC_Setup *del_MC_Setup(MC_Setup *mcsetup);
int MC_Setup_file(MC_Setup *mcsetup, char *inname);
int read_mcsetup_line(char *line, char *keyword);
int setup_monte(MC_Setup *mcsetup);
void setup_monte_summary(MC_Setup *mcsetup);
int get_mu0_mc(MC_Setup *setup);
int get_tau0_mc(MC_Setup *setup);
int get_mcsig(MC_Setup *setup);
int get_nout(MC_Setup *setup);

#endif
