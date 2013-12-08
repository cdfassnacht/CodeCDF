#ifndef catlib_h
#define catlib_h

#include "structdef.h"

#define MASTERLIM 64
#define COMPLIM 64

int find_lens(Secat *cat, int ncat, int *lensindex);
Secat find_closest(Pos cpos, Secat *secat, int ncat, int verbose);
Secat *purge_cat(Secat *in_cat, int nincat, char *catname, int *npurged, 
		 int purgeflag);
int dposcmp(const void *v1, const void *v2);
int dposcmp_sdss(const void *v1, const void *v2);
int dcmp(const void *v1, const void *v2);

#endif
