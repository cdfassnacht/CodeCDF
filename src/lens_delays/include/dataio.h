#ifndef dataio_h
#define dataio_h

#include "structdef.h"
#include "coords.h"

#define MAX 1000

int file_exists(char *filename);
FILE *open_readfile(char *filename);
FILE *open_writefile(char *filename);
FILE *open_appendfile(char *filename, int verbose);
int n_lines(FILE *ifp, char comchar);
int n_cols(char *line, char comchar, int verbose);
Datastruct *read_datastruct(char *inname, char comment, int *nlines,
			    int format);
Secat *read_difmap(char *inname, char comment, int *nlines, Skypos *pos0);
Secat *read_secat(char *inname, char comment, int *nlines, int format);
Secat *read_secat2(char *inname, char comment, int *nlines);
int write_secat(Secat *secat, int ncat, char *outname, int format);
void secat_format();
Secat *read_distcalc(char *inname, char comment, int *nlines, int format);
SDSScat *read_sdss(char *inname, char comment, int *nlines, int format);
int write_sdss(SDSScat *scat, int ncat, char *outname, int format);
int print_offsets(Skypos cent, Skypos *skypos, Pos *offsets, int noffsets,
		  FILE *ofp);

#endif
