/*.......................................................................
 *
 * structdef.c - a library to define specialized structures and their
 *                pointer arrays.
 *
 * v12Oct01 CDF, Added pos2xy function to split x and y components out
 *                of a Pos array.
 * v12Jul03 CDF, Moved new_secat and del_secat functions from matchcat.c
 *               Moved new_skypos and del_skypos functions from coords.c
 */

#include <stdlib.h>
#include <stdio.h>
#include "structdef.h"

/*.......................................................................
 *
 * Strings
 *
 * Functions:
 *  new_string - allocates memory for the string
 *  del_string - frees the character array
 */

char *new_string(int size)
{
  int i;
  char *string;

  string = (char *) malloc(sizeof(char) * size);
  if(!string) {
    fprintf(stderr,"Insufficient memory for string of length %d.\n",
	    size);
    return NULL;
  }
  else
    return string;
}

char *del_string(char *string)
{
  if(string)
    free(string);

  return NULL;
}


/*......................................................................
 *
 * Two-dimensional float arrays
 *
 * Functions:
 *    new_array  -   allocates the array
 *    del_array  -   frees the array
 *
 */

float *new_array(int size1,int size2)
{
  int i;
  float *array,*ptr;

  array = (float *) malloc(sizeof(float) * size1 * size2);
  if(!array) {
    fprintf(stderr,"Insufficient memory for %dx%d data array.\n",
	    size1,size2);
    return NULL;
  }
  ptr = array;
  for (i=0;i<size1*size2;i++,ptr++)
    *ptr = 0.0;
  return array;
}

float *del_array(float *array)
{
  if(array)
    free(array);

  return NULL;
}

/*......................................................................
 *
 * Function to allocate memory for a 2D array of pointers for float:
 *    new_float2  -   allocates the array
 *    del_float2  -   frees the array and the individual members
 *
 */

float **new_float2(int n2) {
  float **newfloat2={NULL};  /* Array of floats to be allocated */

  newfloat2 = (float **) malloc(sizeof(float *) * n2);
  if(!newfloat2) {
    fprintf(stderr,"ERROR:  Insufficient memory for array of pointers\n");
    fprintf(stderr,"        to the spectra.");
    return NULL;
  }
  else
    return newfloat2;
}

float **del_float2(float **float2, int n2) {
  int i;

  if(float2) {
    for(i=0; i<n2; i++)
      float2[i] = del_array(float2[i]);
    free(float2);
  }

  return NULL;
}

/*......................................................................
 *
 * Two-dimensional integer arrays
 *
 * Functions:
 *    new_intarray  -   allocates the array
 *    del_intarray  -   frees the array
 *
 */

int *new_intarray(int size1,int size2)
{
  int i;
  int *array,*ptr;

  array = (int *) malloc(sizeof(int) * size1 * size2);
  if(!array) {
    fprintf(stderr,"Insufficient memory for %dx%d data array.\n",
	    size1,size2);
    return NULL;
  }
  ptr = array;
  for (i=0;i<size1*size2;i++,ptr++)
    *ptr = 0;
  return array;
}

int *del_intarray(int *array)
{
  if(array)
    free(array);

  return NULL;
}

/*......................................................................
 *
 * Double arrays
 *
 * Functions:
 *    new_doubarray  -   allocates the array
 *    del_doubarray  -   frees the array
 *
 */

double *new_doubarray(int size)
{
  int i;
  double *array,*ptr;

  array = (double *) malloc(sizeof(double) * size);
  if(!array) {
    fprintf(stderr,"Insufficient memory for %d data array.\n",size);
    return NULL;
  }
  ptr = array;
  for (i=0; i<size; i++,ptr++)
    *ptr = 0.0;
  return array;
}

double *del_doubarray(double *array)
{
  if(array)
    free(array);

  return NULL;
}

/*......................................................................
 *
 * Function to allocate memory for a 2D array of pointers for doubles:
 *    new_doub2  -   allocates the array
 *    del_doub2  -   frees the array and the individual members
 *
 */

double **new_doub2(int ndoub2) {
  double **newdoub2={NULL};  /* Array of doubles to be allocated */

  newdoub2 = (double **) malloc(sizeof(double *) * ndoub2);
  if(!newdoub2) {
    fprintf(stderr,"ERROR:  Insufficient memory for array of pointers\n");
    fprintf(stderr,"        to the spectra.");
    return NULL;
  }
  else
    return newdoub2;
}

double **del_doub2(double **doub2, int ndoub2) {
  int i;

  if(doub2) {
    for(i=0; i<ndoub2; i++)
      doub2[i] = del_doubarray(doub2[i]);
    free(doub2);
  }

  return NULL;
}

/*.......................................................................
 *
 * Function new_datastruct
 *
 * Allocates dynamic memory for a pointer array of Datastruct structures
 *
 * Input:  int size            size of the array
 *
 * Output: Datastruct *newinfo    pointer to the new array.  NULL if error
 *
 */

Datastruct *new_datastruct(int size)
{
  int i;
  Datastruct *newinfo;
  Datastruct *dptr;
 
  newinfo = (Datastruct *) malloc(sizeof(Datastruct) * size);
  if(!newinfo) {
    fprintf(stderr,"new_datastruct: \n");
    fprintf(stderr,"Insufficient memory for Datastruct array.\n");
    return NULL;
  }

  /*
   * Initialize values in the array
   */

  for(i=0,dptr=newinfo; i<size; i++,dptr++) {
    dptr->x = dptr->y = dptr->z = 0.0;
    dptr->dataflag = 0;
  }

  return newinfo;
}

/*.......................................................................
 *
 * Function del_datastruct
 *
 * Frees up memory allocated to Datastruct array
 *
 * Input:  Datastruct *datastruct    array to be freed
 *
 * Output: NULL
 */

Datastruct *del_datastruct(Datastruct *datastruct)
{
  if(datastruct)
    free(datastruct);
 
  return NULL;
}

/*......................................................................
 *
 * Functions to allocate and free up pointer arrays of the Pos
 *   structure, which contains positions in the form of two doubles,
 *   x and y, the errors on the positions (xerr and yerr), and a
 *   flag (flag).
 *
 * Functions:
 *    new_pos  -   allocates the array
 *    del_pos  -   frees the array
 *    pos2xy   -   extracts the x and y components from a Pos array
 *
 */

Pos *new_pos(int n1,int n2)
{
  int i;           /* Looping variable */
  Pos *pos=NULL;   /* New array to be allocated */
  Pos *pptr;       /* Pointer to navigate array */

  /*
   * Allocate memory for array
   */

  pos = (Pos *) malloc(sizeof(Pos) * n1 * n2);
  if(!pos) {
    fprintf(stderr,"Insufficient memory for a %dx%d position array.\n",n1,n2);
    return NULL;
  }

  /*
   * Initialize the array
   */
  
  for(i=0,pptr=pos; i < n1*n2; i++,pptr++) {
    pptr->x = pptr->y = pptr->xerr = pptr->yerr = 0.0;
    pptr->flag = 0;
  }

  return pos;
}

Pos *del_pos(Pos *pos)
{
  if (pos)
    free (pos);

  return NULL;
}

float **pos2xy(Pos *pos, int npts, int *ngood)
{
  int i;                 /* Looping variable */
  float **xy={NULL};     /* x and y components of the Pos array */
  float *xptr;           /* Pointer to navigate the x array */
  float *yptr;           /* Pointer to navigate the y array */
  Pos *pptr;             /* Pointer to navigate the Pos array */

  /*
   * Allocate memory
   */

  if(!(xy = new_float2(2))) {
    fprintf(stderr,"ERROR: pos2xy.\n");
    return NULL;
  }

  if(!(xy[0] = new_array(npts,1))) {
    fprintf(stderr,"ERROR: pos2xy\n");
    return del_float2(xy,2);
  }

  if(!(xy[1] = new_array(npts,1))) {
    fprintf(stderr,"ERROR: pos2xy\n");
    return del_float2(xy,2);
  }

  /*
   * Loop through Pos array transferring x and y to the xy array.
   */

  *ngood = 0;
  for(i=0,pptr=pos,xptr=xy[0],yptr=xy[1]; i<npts; 
      i++,pptr++,xptr++,yptr++) {
    if(pptr->flag == 0) {
      *xptr = (float) pptr->x;
      *yptr = (float) pptr->y;
      (*ngood)++;
    }
  }

  return xy;
}

/*.......................................................................
 *
 * Function new_skypos
 *
 * Allocates dynamic memory for a pointer array of Skypos structures
 *
 * Input:  int size            size of the array
 *
 * Output: Skypos *newinfo     pointer to the new array.  NULL if error
 *
 */
 
 
Skypos *new_skypos(int size)
{
  Skypos *newinfo;
 
  newinfo = (Skypos *) malloc(sizeof(Skypos) * size);
  if(!newinfo) {
    fprintf(stderr,"Insufficient memory for data array.\n");
    return NULL;
  }
 
  return newinfo;
}
 
/*.......................................................................
 *
 * Function del_skypos
 *
 * Frees up memory allocated to Skypos array
 *
 * Input:  Skypos *skypos      array to be freed
 *
 * Output: NULL
 */
 
Skypos *del_skypos(Skypos *skypos)
{
  if(skypos)
    free(skypos);
 
  return NULL;
}

/*.......................................................................
 *
 * Function new_axisinfo
 *
 * Allocates dynamic memory for a pointer array of Axisinfo structures
 *
 * Input:  int size            size of the array
 *
 * Output: Axisinfo *newinfo    pointer to the new array.  NULL if error
 *
 */

Axisinfo *new_axisinfo(int size)
{
  int i;
  Axisinfo *newinfo;
  Axisinfo *aptr;
 
  newinfo = (Axisinfo *) malloc(sizeof(Axisinfo) * size);
  if(!newinfo) {
    fprintf(stderr,"new_axisinfo: \n");
    fprintf(stderr,"Insufficient memory for Axisinfo array.\n");
    return NULL;
  }

  /*
   * Initialize values in the array
   */

  for(i=0,aptr=newinfo; i<size; i++,aptr++) {
    aptr->nval = 0;
    aptr->minval = aptr->maxval = 0.0;
  }

  return newinfo;
}

/*.......................................................................
 *
 * Function del_axisinfo
 *
 * Frees up memory allocated to Axisinfo array
 *
 * Input:  Axisinfo *axisinfo    array to be freed
 *
 * Output: NULL
 */

Axisinfo *del_axisinfo(Axisinfo *axisinfo)
{
  if(axisinfo)
    free(axisinfo);
 
  return NULL;
}

/*.......................................................................
 *
 * Function new_secat
 *
 * Allocates dynamic memory for a pointer array of Secat structures
 *
 * Input:  int size            size of the array
 *
 * Output: Secat *newinfo    pointer to the new array.  NULL if error
 *
 */

Secat *new_secat(int size)
{
  int i;
  Secat *newinfo;
  Secat *sptr;
 
  newinfo = (Secat *) malloc(sizeof(Secat) * size);
  if(!newinfo) {
    fprintf(stderr,"new_secat: \n");
    fprintf(stderr,"Insufficient memory for Secat array.\n");
    return NULL;
  }

  /*
   * Initialize values in the array
   */

  for(i=0,sptr=newinfo; i<size; i++,sptr++) {
    sptr->x = sptr->y = 0.0;
    sptr->ma1 = sptr->ma2 = sptr->ma3 = sptr->mtot = sptr->miso = -99.0;
    sptr->matchflag[0] = sptr->nmatch = 0;
    sptr->matchid[0] = 0;
  }

  return newinfo;
}

/*.......................................................................
 *
 * Function del_secat
 *
 * Frees up memory allocated to Secat array
 *
 * Input:  Secat *secat    array to be freed
 *
 * Output: NULL
 */

Secat *del_secat(Secat *secat)
{
  if(secat)
    free(secat);
 
  return NULL;
}

/*.......................................................................
 *
 * Function new_sdsscat
 *
 * Allocates dynamic memory for a pointer array of SDSScat structures
 *
 * Input:  int size            size of the array
 *
 * Output: SDSScat *newinfo    pointer to the new array.  NULL if error
 *
 */

SDSScat *new_sdsscat(int size)
{
  int i;
  SDSScat *newinfo;
  SDSScat *sptr;
 
  newinfo = (SDSScat *) malloc(sizeof(SDSScat) * size);
  if(!newinfo) {
    fprintf(stderr,"new_sdsscat: \n");
    fprintf(stderr,"Insufficient memory for SDSScat array.\n");
    return NULL;
  }

  /*
   * Initialize values in the array
   */

  for(i=0,sptr=newinfo; i<size; i++,sptr++) {
    sptr->alpha = sptr->delta = 0.0;
    sptr->u = sptr->g = sptr->r = sptr->i = sptr->z = -99.0;
  }

  return newinfo;
}

/*.......................................................................
 *
 * Function del_sdsscat
 *
 * Frees up memory allocated to SDSScat array
 *
 * Input:  SDSScat *sdsscat    array to be freed
 *
 * Output: NULL
 */

SDSScat *del_sdsscat(SDSScat *sdsscat)
{
  if(sdsscat)
    free(sdsscat);
 
  return NULL;
}

/*......................................................................
 *
 * Functions to allocate and free up pointer arrays of the Lensparam
 *   structure, which contains information on lens model parameters.
 *
 * Functions:
 *    new_lp  -   allocates the array
 *    del_lp  -   frees the array
 *
 */

Lensparams *new_lp(int n)
{
  Lensparams *lp;

  lp = (Lensparams *) malloc(sizeof(Lensparams)*n);
  if(!lp) {
    fprintf(stderr,"Insufficient memory for lens parameters.\n");
    return NULL;
  }
  return lp;
}

Lensparams *del_lp(Lensparams *lp)
{
  if(lp)
    free(lp);

  return NULL;
}

/*......................................................................
 *
 * Functions to allocate and free up pointer arrays of the Modparam
 *   structure, which contains a parameter value (a float) and a
 *   flag (an int) which tells whether or not the parameter varies in
 *   the model fitting
 *
 * Functions:
 *    new_modparam  -   allocates the array
 *    del_modparam  -   frees the array
 *
 */

Modparam *new_modparam(int size)
{
  Modparam *modparam;

  modparam = (Modparam *) malloc(sizeof(Modparam) * size);
  if(!modparam) {
    fprintf(stderr,"Insufficient memory for a %d Modparam array.\n",size);
    return NULL;
  }
  return modparam;
}

Modparam *del_modparam(Modparam *modparam)
{
  if (modparam)
    free (modparam);

  return NULL;
}

/*......................................................................
 *
 * Functions to allocate and free up pointer arrays of the Poten
 *   structure, which contains information on a two-dimensional 
 *   gravitational potential and its first and second derivatives.
 *
 * Functions:
 *    new_poten    -   allocates the array
 *    reinit_poten -   sets all members of an existing structure array to
 *                       zero
 *    del_poten    -   frees the array
 *
 */

Poten *new_poten(int n1,int n2)
{
  int i;
  Poten *ptmp,*tmpptr;

  ptmp = (Poten *) malloc(sizeof(Poten) * n1 * n2);
  if(!ptmp) {
    fprintf(stderr,"Insufficient memory for a %dx%d potential array.\n",n1,n2);
    return NULL;
  }
  tmpptr = ptmp;
  for(i=0;i<n1*n2;i++,tmpptr++) {
    tmpptr->ph = 0.0;
    tmpptr->dphdx = 0.0;
    tmpptr->dphdy = 0.0;
    tmpptr->dphdxx = 0.0;
    tmpptr->dphdyy = 0.0;
    tmpptr->dphdxy = 0.0;
  }

  return ptmp;
}

void reinit_poten(Poten *potptr,int n1,int n2)
{
  int i;
  Poten *tmpptr;

  tmpptr = potptr;
  for(i=0;i<n1*n2;i++,tmpptr++) {
    tmpptr->ph = 0.0;
    tmpptr->dphdx = 0.0;
    tmpptr->dphdy = 0.0;
    tmpptr->dphdxx = 0.0;
    tmpptr->dphdyy = 0.0;
    tmpptr->dphdxy = 0.0;
  }
}

Poten *del_poten(Poten *potptr)
{
  if(potptr)
    free(potptr);

  return NULL;
}

/*......................................................................
 *
 * Functions to allocate and free up pointer arrays of the Posinfo
 *   structure, which contains information on model image positions
 *   and fluxes, as well as the position and flux for the background
 *   source associated with the image.
 *
 * Functions:
 *    new_posinfo  -   allocates the array
 *    del_posinfo  -   frees the array
 *
 */

Posinfo *new_posinfo(int n1, int n2)
{
  Posinfo *newposinfo;

  newposinfo = (Posinfo *) malloc(sizeof(Posinfo)*n1*n2);
  if(!(newposinfo)) {
    fprintf(stderr,"Insufficient memory for posinfo.\n");
    return NULL;
  }

  return newposinfo;
}

Posinfo *del_posinfo(Posinfo *posinfo)
{
  if(posinfo)
    free(posinfo);

  return NULL;
}

/*......................................................................
 *
 * Functions to allocate and free up pointer arrays of the Srcinfo
 *   structure, which contains information on background source
 *   positions and fluxes.
 *
 * Functions:
 *    new_srcinfo  -   allocates the array
 *    del_srcinfo  -   frees the array
 *
 */

Srcinfo *new_srcinfo(int n1, int n2)
{
  Srcinfo *newsrcinfo;

  newsrcinfo = (Srcinfo *) malloc(sizeof(Srcinfo)*n1*n2);
  if(!(newsrcinfo)) {
    fprintf(stderr,"Insufficient memory for srcinfo.\n");
    return NULL;
  }

  return newsrcinfo;
}

Srcinfo *del_srcinfo(Srcinfo *srcinfo)
{
  if(srcinfo)
    free(srcinfo);

  return NULL;
}

/*......................................................................
 *
 * Functions to allocate and free up pointer arrays of the Modinfo
 *   structure, which contains all information on a lens model.
 *
 * Functions:
 *    new_modinfo         -   allocates the array
 *    modpars_to_modinfo  -   converts model-fitting output to a
 *                               modinfo structure
 *    del_modinfo         -   frees the array
 *
 */

Modinfo *new_modinfo()
{
  Modinfo *newinfo;

  newinfo = (Modinfo *) malloc(sizeof(Modinfo));
  if(!newinfo) {
    fprintf(stderr,"Insufficient memory for lens parameters.\n");
    return NULL;
  }

  newinfo->lp = NULL;
  newinfo->nmod = 0;
  newinfo->fsource = NULL;
  newinfo->nfsrc = 0;
  newinfo->var = NULL;
  newinfo->nvar = 0;
  newinfo->nlinks = 0;
  newinfo->posobs = NULL;
  newinfo->nobs = 0;
  newinfo->posimmod = NULL;
  newinfo->nimmod = 0;
  newinfo->pixscale = 0.0;
  newinfo->poserr = 0.0;
  newinfo->ferr = 0.0;
  newinfo->chisqflag = 0;
  newinfo->nlinks = 0;

  return newinfo;
}

int modpars_to_modinfo(float *modpars, Modinfo *modinfo)
{
  int i;
  int maxvar;
  int *intptr;
  float *linkvals;
  float *fptr;
  Srcinfo *srcptr;
  Lensparams *parptr;

  if (!(linkvals = new_array(modinfo->nlinks,1))) {
    fprintf(stderr,"ERROR: modpars_to_modinfo\n");
    return 1;
  }

  intptr = modinfo->var;
  fptr = modpars;
  maxvar = 1;
  parptr = modinfo->lp;
  for(i=0; i<modinfo->nmod; i++,parptr++) {
    switch(modpars_to_struct(&parptr->xl.par,intptr++,&maxvar,fptr,linkvals)) {
    case 0:
      break;
    case 1:
      fptr++;
      break;
    default:
      fprintf(stderr,"Error: modinfo_to_info\n");
      return 1;
    }
    switch(modpars_to_struct(&parptr->yl.par,intptr++,&maxvar,fptr,
			     linkvals)) {
    case 0:
      break;
    case 1:
      fptr++;
      break;
    default:
      fprintf(stderr,"Error: modinfo_to_info\n");
      return 1;
    }
    switch(modpars_to_struct(&parptr->b.par,intptr++,&maxvar,fptr,linkvals)) {
    case 0:
      break;
    case 1:
      fptr++;
      break;
    default:
      fprintf(stderr,"Error: modinfo_to_info\n");
      return 1;
    }
    switch(modpars_to_struct(&parptr->gamma.par,intptr++,&maxvar,fptr,linkvals)) {
    case 0:
      break;
    case 1:
      fptr++;
      break;
    default:
      fprintf(stderr,"Error: modinfo_to_info\n");
      return 1;
    }
    switch(modpars_to_struct(&parptr->theta_g.par,intptr++,&maxvar,
			     fptr,linkvals)) {
    case 0:
      break;
    case 1:
      fptr++;
      break;
    default:
      fprintf(stderr,"Error: modinfo_to_info\n");
      return 1;
    }
  }

  srcptr = modinfo->fsource;
  for(i=0; i<modinfo->nfsrc; i++,srcptr++) {
    switch(modpars_to_struct(&srcptr->x.par,intptr++,&maxvar,fptr,linkvals)) {
    case 0:
      break;
    case 1:
      fptr++;
      break;
    default:
      fprintf(stderr,"Error: modinfo_to_info\n");
      return 1;
    }
    switch(modpars_to_struct(&srcptr->y.par,intptr++,&maxvar,fptr,linkvals)) {
    case 0:
      break;
    case 1:
      fptr++;
      break;
    default:
      fprintf(stderr,"Error: modinfo_to_info\n");
      return 1;
    }
    switch(modpars_to_struct(&srcptr->flux.par,intptr++,&maxvar,fptr,
			     linkvals)) {
    case 0:
      break;
    case 1:
      fptr++;
      break;
    default:
      fprintf(stderr,"Error: modinfo_to_info\n");
      return 1;
    }
    switch(modpars_to_struct(&srcptr->zf.par,intptr++,&maxvar,fptr,linkvals)) {
    case 0:
      break;
    case 1:
      fptr++;
      break;
    default:
      fprintf(stderr,"Error: modinfo_to_info\n");
      return 1;
    }
  }

  linkvals = del_array(linkvals);

  return 0;
}

Modinfo *del_modinfo(Modinfo *info)
{
  info->lp = del_lp(info->lp);
  info->fsource = del_srcinfo(info->fsource);
  info->var = del_intarray(info->var);
  info->posobs = del_posinfo(info->posobs);
  info->posimmod = del_posinfo(info->posimmod);

  if(info)
    free(info);

  return NULL;
}

/*.......................................................................
 *
 * Function modpars_to_struct
 *
 * Converts model parameter output from fitting routine into float 
 *   components of a structure.
 *
 * Inputs:  float *valptr      pointer to the structure component
 *          int *varptr        tells whether the component was modified
 *                               by the model fitting
 *          int *maxvar        keeps track of varying parameters
 *          float *mpptr       pointer to the model parameters
 *          float *linkvals    values of linked varying parameters
 *
 * Output:  int (0, 1 or -1)   0 ==> don't increment the pointers
 *                             1 ==> increment the pointers
 *                             -1 ==> error
 *
 */

int modpars_to_struct(float *valptr, int *varptr, int *maxvar, float *mpptr,
		      float *linkvals)
{
  if (*varptr == 1) {
    *valptr = *mpptr;
    return 1;
  }
  else if (*varptr > *maxvar) {
    *valptr = *mpptr;
    *maxvar = *varptr;
    *(linkvals+(*varptr)-2) = *mpptr;
    return 1;
  }
  else if (*varptr > 1) {
    *valptr = *(linkvals+(*varptr)-2);
    return 0;
  }
  else
    return 0;

  return 0;
}


/*......................................................................
 *
 * Functions to allocate and free up pointer arrays of the Specinfo
 *   structure, which contains information on spectral lines
 *
 * Functions:
 *    new_specinfo  -   allocates the array
 *    del_specinfo  -   frees the array
 *
 */

Specinfo *new_specinfo(int n)
{
  Specinfo *newspecinfo;

  newspecinfo = (Specinfo *) malloc(sizeof(Specinfo)*n);
  if(!(newspecinfo)) {
    fprintf(stderr,"Insufficient memory for specinfo.\n");
    return NULL;
  }

  return newspecinfo;
}

Specinfo *del_specinfo(Specinfo *specinfo)
{
  if(specinfo)
    free(specinfo);

  return NULL;
}

/*......................................................................
 *
 * Functions to allocate and free up pointer arrays of the Labinfo
 *   structure, which contains positions and text of labels for images
 *
 * Functions:
 *    new_labinfo  -   allocates the array
 *    del_labinfo  -   frees the array
 *
 */

Labinfo *new_labinfo(int n)
{
  Labinfo *newlabinfo;

  newlabinfo = (Labinfo *) malloc(sizeof(Labinfo)*n);
  if(!(newlabinfo)) {
    fprintf(stderr,"Insufficient memory for labinfo.\n");
    return NULL;
  }

  return newlabinfo;
}

Labinfo *del_labinfo(Labinfo *labinfo)
{
  if(labinfo)
    free(labinfo);

  return NULL;
}

