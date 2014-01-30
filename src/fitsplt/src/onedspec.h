#ifndef onedspec_h
#define onedspec_h

#include "structdef.h"

/*.......................................................................
 *
 * Definitions
 *
 */

#define MAXFLUX 99999.99
#define MINFLUX -99999.99
#define SHORT 10.0
#define LONGIN 50.0
#define LONGOUT 100.0

/*.......................................................................
 *
 * Enumerations
 *
 */

enum {
  BADOPT = -1,
  ONESPEC,
  OVERPLOT,
  MULTISPEC
};

/*.......................................................................
 *
 * Structure definitions
 *
 */

typedef struct {
  int option;          /* Type of plotting to be done (see enumeration) */
  int flchoice;        /* Choice of flux units */
  int labstyle;        /* Choice of style for labeling spectral lines */
  float labyfrac;      /* y location of line labels as fraction of plot */
  float labheight;     /* Character height for line labels */
  float cheight;       /* Character height for plot labels */
  double pltxmin;      /* Minimum wavelength to be plotted */
  double pltxmax;      /* Maximum wavelength to be plotted */
} Spsetup;

typedef struct {
  Pos *spec;           /* Actual spectrum */
  int nlines;          /* Number of points in the spectrum */
  float minlambda;     /* Wavelength of first point in spectrum */
  float maxlambda;     /* Wavelength of last point in spectrum */
  float minflux;       /* Minimum flux level to plot */
  float maxflux;       /* Maximum flux level to plot */
  int dolabel;         /* Flag set to 1 if internal label is desired */
  char labpos;         /* Position for optional label ('l' or 'r') */
  char labtext[MAXC];  /* Optional label text */
} Spectrum;

/*.......................................................................
 *
 * Function declarations
 */

Spectrum *new_spectrum(int nspec);
Spectrum *del_spectrum(Spectrum *spec);
Pos **new_spec(int nspec);
Pos **del_spec(Pos **spec, int nspec);
Spectrum *read_onespec(char *filename);
Spectrum *read_multispec(char *filename, int *nspec);
int read_spectrum(char *filename, Spectrum *spectrum);
int read_spec_list(Spectrum *spectrum, FILE *ifp);
Specinfo *read_line_list(int *goodlines);
int read_line(char *line, Specinfo *specinfo);
int display_spectrum(Spectrum *spec, int nspec, int specnum, 
		     Spsetup *setup, int option, int dobox);
int adjust_spec(Spectrum *spectrum, int nspec, int specnum, 
		Spsetup *setup, int option, int dobox);
void set_pltlims_y(float minflux, float maxflux, int flchoice, float *pltymin, 
		   float *pltymax);
int get_vscale(float *minflux, float *maxflux);
void print_line_list(Specinfo *specinfo, int goodlines);
int label_lines(Specinfo *specinfo, int goodlines, float z, Spectrum *spectrum,
		Spsetup *setup);
int label_lines_file(Spectrum *spectrum, Spsetup *setup);
int speclab(float marklamb, Labinfo label, Spectrum *spectrum, Spsetup *setup);
int get_sky_label(float *marklamb, Labinfo *label);
int label_plot(char *title, int titleopt, int flchoice, float cheight);
void labcorner(char *text, char side, int option);
int calc_d4000(Pos *spectrum, int npoints);
int flambda_to_fnu(Pos *spectrum, int npoints, float *minflux, 
		   float *maxflux);

/*
 * Functions in setup_spec.c
 */

Spsetup *new_spsetup();
Spsetup *del_spsetup(Spsetup *spsetup);
void get_flchoice(Spsetup *setup);
void get_cheight(Spsetup *setup);
void get_labheight(Spsetup *setup);
int get_label_style(Spsetup *setup, int option);
void get_xrange(Spsetup *setup, Spectrum *spec);
int get_title(char *title, int *titleopt);
int get_clabel(char *labtext, char *labside);

#endif
