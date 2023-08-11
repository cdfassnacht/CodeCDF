#ifndef coords_h
#define coords_h

#include "structdef.h"
#define CMAXC 200

/*.......................................................................
 *
 * Structure definitions
 *
 */

/*.......................................................................
 *
 * Function declarations
 *
 */

Skypos *new_skypos(int size);
Skypos *del_skypos(Skypos *skypos);
void print_skypos(FILE *fp, Skypos pos);
Pos read_skypos();
void rth2xy(double r, double theta, Pos *pos, int isastro);
void xy2rth(Pos pos, double *r, double *theta, int isastro);
void deg2rad(double dalpha, double ddelta, double *ralpha, double *rdelta);
void spos2rad(Skypos spos, double *alpha, double *delta);
void rad2spos(double alpha, double delta, Skypos *spos);
void spos2deg(Skypos spos, double *alphadeg, double *deltadeg);
void deg2spos(double alphdeg, double deltdeg, Skypos *spos);
Pos *dspos2xy(Skypos cent, Skypos *spos, int npos);
Pos ddeg2xy(double alphadeg0, double deltadeg0, double alphadeg,
	    double deltadeg);
Pos rad2offset(double alpha1, double delta1, double alpha2, double delta2);
int offset2rad(double alpha0, double delta0, Pos offset, double *alpha,
	       double *delta);
Skypos *dspos2spos(Skypos cent, Pos *pos, int npos);
Pos *ddeg2deg(Pos cent, Pos *pos, int npos);
int mod_center(FILE *ifp, Skypos *center);
int secat2offset(Secat *secat, int ncat, int posformat, Secat centpos,
		 int centformat);

#endif
