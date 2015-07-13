#ifndef modfuncs_h
#define modfuncs_h

#define MMAX 100
#define REDDIR "../Analysis/"   /* Output directory */
#define MODDIR "../../Models/"  /* Directory containing models */


/*.......................................................................
 *
 * Enumeration for calibrator/lens separation.
 *
 */

enum {
  CAL,
  LENS
};

/*.......................................................................
 *
 * Structure definitions
 *
 */

typedef struct {
  float flux;    /* Component flux */
  float radius;  /* Distance from phase center */
  float theta;   /* PA w/r.t. phase center */
  float maj;     /* Major axis (for gaussian components) */
  float axrat;   /* Axial ratio ( "    "         "     ) */
  float phi;     /* Comp. PA   (  "    "         "     ) */
  int type;      /* Model type (1 if gaussian, 0 if point) */
} Compinfo;

typedef struct {
  int model;
  float scale;
  float flux;
  float chisq;
  float rms;
} Modfit;

/*.......................................................................
 *
 * Function declarations here
 *
 */

void get_ext(char *ext);
void default_caltimes(char *root, float *caltime);
Compinfo *new_cinfo(int size);
Compinfo *del_cinfo(Compinfo *cinfo);
void get_modname(char *modname, char *root, char *ext, float obsdate);
FILE *open_modfile(char *modname);
Compinfo *fill_cinfo(FILE *ifp, int nlines);
int read_cinfo(Compinfo *cinfo, FILE *ifp);
int get_sinfo(float *min, float *max, float *step);
int scale_data(float min, float max, float step, char *basemod);
int modfile_08(char *root, char *ext, float caltime, int interact);
int calc_totflux(Compinfo *cinfo, int nlines, float *totflux);
int print_compmod(int modnum, float scale, float totflux, 
		  Compinfo *cinfo, int nlines, char *outname);
void difmap_headers(FILE *scrfp, char *root, char *ext, float obsdate, 
		    float caltime, int interact);
void setvars(char *root, char *ext, float obsdate, char *shortext, 
	     int *mapsize, float *pixsize, int *doshift);
int check_uvf(char *root, char *ext, float *obsdate);
int get_obsdate(char *filename, float *obsdate);
void set_log_obsdate(double *obsdate);
int read_in_params(FILE *ifp, char *root, char *ext, float *min, float *max, 
		   float *step, float *caltime, int *interact);
Modfit *new_modfit(int size);
Modfit *del_modfit(Modfit *modfit);
int cal_flux(char *filename, float *scale);
int get_best_cal(char *filename, float *flux, float *err, float *chisq);
Modfit *read_modfit(char *filename, int *nlines);
FILE *open_file(char *name);
int fill_modfit(Modfit *modfit, int nlines, FILE *ifp);
float best_model(Modfit *modfit, int nlines);
Secat *source_flux(char *date, char *root, char *modtype, int *ncomp,
		   float *rms, float *chisq, double *obsdate);
int get_clean_flux(float *flux, char *modname);
int print_lensflux(char *root, float *fluxp, float *fluxg, int ncomp, 
		   float rmsp, float rmsg, float chisqp, float chisqg, 
		   double obsdate);
int print_calflux(char *root, float *fluxp, float *fluxg, int ncomp, 
		  float rmsp, float rmsg, float chisqp, float chisqg, 
		  double obsdate);
float fit_parab(float x1, float y1, float x2, float y2, float x3, float y3);

#endif
