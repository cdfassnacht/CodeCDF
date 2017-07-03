/*
 * dataio.c
 *
 * This is a library of useful data I/O functions:
 *
 * Functions included in this library:
 * -----------------------------------
 *  file_exists     - checks for the existence of a file
 *  open_readfile   - opens a read-only file
 *  open_writefile  - opens a writeable file
 *  open_appendfile - opens a file for appending
 *  n_lines         - counts the number of data lines in a file
 *  n_cols          - counts the number of columns in a line
 *  read_datastruct - creates a Datastruct array and fills it from a file.
 *  read_difmap     - reads a difmap model file
 *  read_secat      - reads a SExtractor catalog
 *  read_secat2     - a more generalized version of read_secat that will
 *                    not (hopefully) need a format flag.  Not yet
 *                    functional.
 *  write_secat     - writes a Secat array to an output file
 *  secat_format    - describes the formats acceptable by the secat functions
 *  read_distcalc   - reads the output file produced by distcalc.c
 *  read_sdss       - reads a SDSS-style file, in which magnitudes in multiple
 *                     bands are contained in a single input file
 *  write_sdss      - writes out a SDSS-style file, with more positional info
 *  print_offsets   - prints offsets such as those produced by distcalc.c
 *
 *-----------------------------------------------------------------------
 * Revision history:
 * -----------------
 * v15Sep2000 CDF, Added new file_exists and open_appendfile functions.
 * v12Jul2003 CDF, Moved read_secat function from matchcat.c
 *                 Added read_distcalc function
 * v13Jul2003 CDF, Added calculation of astronomical PA to read_distcalc
 * v20Jul2003 CDF, Added processing of alpha, delta pairs in read_secat
 * v07Aug2003 CDF, Added a format flag to read_distcalc to enable reading
 *                  of either input to or output of distcalc.c.
 * v11Mar2004 CDF, Added another format flag to read_secat to read output
 *                  of 2MASS server, among others.
 * v05Aug2004 CDF, Added a distcalc format flag to read_secat
 * v20Jul2005 CDF, Added a USNO star format flag to read_secat
 * v08Aug2005 CDF, Added a format flag to read_datastruct
 * v2005Sep08 CDF, Added error checking on format code in write_secat
 *                 Added handling of 2MASS (processed by catsort) to write_secat
 * v2005Oct07 CDF, Added a new read_difmap function (basically the read_comps
 *                  function from pos.c in the class directory)
 * v2006Feb08 CDF, Fixed a bug in the read_difmap function
 * v2006Jun07 CDF, Added a new format to read_datastruct
 * v2007Jan06 CDF, Added a new format to read_secat
 * v2007Jan09 CDF, Added a new function, secat_format, to print out help
 *                  on possible formats for a sextractor file.
 *                 Added new functions read_sdss and write_sdss to deal with
 *                  SDSS-style files, which contain magnitudes in multiple 
 *                  bands all within a single input file.
 * v2007Jan16 CDF, Replaced an old hard-wired SExtractor input format in
 *                  read_secat with a more generic format.  Also impacts
 *                  write_secat.
 * v2007Jul09 CDF, Modifications to read_secat, write_secat, and secat_format
 * v2007Jul29 CDF, Added format=14 to write_secat and secat_format
 * v2007Aug02 CDF, Added format=15 to read_secat, write_secat, and secat_format
 * v2008Jun26 CDF, Modified format 15 to match the new redshift catalog format
 * v2009Jan24 CDF, Changed format 0 from id,x,y to id,alpha,delta
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structdef.h"
#include "coords.h"
#include "dataio.h"

/*.......................................................................
 *
 * Function file_exists
 *
 * Checks for the existence of a file by trying to open it in read-only
 *  mode.  If the file exists, return a 1 (YES), otherwise return a
 *  0 (NO).
 *
 * Inputs: char *filename      filename to check
 *
 * Output: int (0 or 1)        *** NB: 1 ==> file exists, 0 ==> no file.
 *
 */

int file_exists(char *filename)
{
  int exists=0;    /* Set to 1 if file exits */
  FILE *fp=NULL;   /* File pointer */

  /*
   * Try to open the output file in read-only mode.
   */

  if((fp = fopen(filename,"r")) != NULL)
    exists = 1;

  /*
   * Close the output file, if it was opened.
   */

  if(fp)
    fclose(fp);

  return exists;
}

/*.......................................................................
 *
 * Function open_readfile
 *
 * Checks for existence of an input file and, if the file exists opens it.
 *  Returns the file pointer.
 *
 * Inputs: char *filename      input file name
 *
 * Output: FILE *ifp           pointer to opened input file
 *
 */

FILE *open_readfile(char *filename)
{
  char line[MAX];  /* General string for reading input */
  FILE *ifp=NULL;  /* File pointer for opened file */

  /*
   * Try to open file
   */

  while((ifp = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"ERROR: open_readfile.  Cannot open %s\n",filename);
    fprintf(stderr," Enter new file name:  ");
    fgets(line,MAX,stdin);
    while(sscanf(line,"%s",filename) != 1 || 
	  (ifp = fopen(filename,"r")) == NULL) {
      fprintf(stderr,"ERROR: bad input.  Enter filename again: ");
      fgets(line,MAX,stdin);
    }
  }

  printf("open_readfile: %s now open\n",filename);
  return ifp;
}

/*.......................................................................
 *
 * Function open_writefile
 *
 * Opens an output file for writing.
 *  Returns the file pointer.
 *
 * Inputs: char *filename      input file name
 *
 * Output: FILE *ifp           pointer to opened input file
 *
 */

FILE *open_writefile(char *filename)
{
  char line[MAX];    /* General string for reading input */
  char newname[MAX]; /* New name needed if error opening file */
  FILE *ofp=NULL;    /* File pointer for opened file */

  /*
   * Try to open file
   */

  if((ofp = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"ERROR: open_writefile.  Cannot open %s.\n",filename);
    fprintf(stderr," Enter new file name:  ");
    fgets(line,MAX,stdin);
    while(sscanf(line,"%s",newname) != 1 || 
	  (ofp = fopen(newname,"w")) == NULL) {
      fprintf(stderr,"ERROR: bad input.  Enter filename again: ");
      fgets(line,MAX,stdin);
    }
  }
  
  printf("open_writefile: %s now open\n",filename);
  return ofp;
}

/*.......................................................................
 *
 * Function open_appendfile
 *
 * Opens an output file for appending.
 *  Returns the file pointer.
 *
 * Inputs: char *filename      input file name
 *         int verbose         set to 1 for verbose output
 *
 * Output: FILE *ifp           pointer to opened input file
 *
 */

FILE *open_appendfile(char *filename, int verbose)
{
  char line[MAX];  /* General string for reading input */
  FILE *ofp=NULL;  /* File pointer for opened file */

  /*
   * Try to open file
   */

  while((ofp = fopen(filename,"a")) == NULL) {
    fprintf(stderr,"ERROR: open_appendfile.  Cannot open %s.\n",filename);
    fprintf(stderr," Enter new file name:  ");
    fgets(line,MAX,stdin);
    while(sscanf(line,"%s",filename) != 1 || 
	  (ofp = fopen(filename,"w")) == NULL) {
      fprintf(stderr,"ERROR: bad input.  Enter filename again: ");
      fgets(line,MAX,stdin);
    }
  }

  if(verbose)
    printf("open_appendfile: %s now open\n",filename);
  return ofp;
}

/*.......................................................................
 *
 * Function n_lines
 *
 * Reads input file and counts number of data lines
 *
 * Input:  FILE *ifp           input file
 *         char comchar        comment character
 *
 * Output: int nlines          number of lines
 *
 */

int n_lines(FILE *ifp, char comchar)
{
  int nlines=0;
  char line[MAX];

  /*
   * Check input
   */

  if(ifp == NULL) {
    fprintf(stderr,"ERROR: n_lines.  Null input file.\n");
    return 0;
  }

  while(fgets(line,MAX,ifp) != NULL) {

    /*
     * First check to see that the line doesn't have the characteristic
     *  of a comment line (which begins with a comment character)
     */

    if(line[0] != comchar) {

      /*
       * Increment the nlines counter
       */

      nlines++;
    }
  }

  if(nlines == 0)
    fprintf(stderr,"** Warning.  n_lines:  No valid lines in input file\n");
  else
    printf("n_lines: %d lines in input file\n",nlines);

  return nlines;
}

/*.......................................................................
 *
 * Function n_cols
 *
 * Reads the number of columns in a line, which is passed to the
 *  function as a string
 *
 * Input:  char *line          input line
 *         char comchar        comment character
 *         int verbose         print out information if set to 1
 *
 * Output: int ncols           number of columns in the line
 *
 */

int n_cols(char *line, char comchar, int verbose)
{
  int ncols=0;      /* Number of columns in line */
  char *cptr;       /* Pointer to navigate line */

  /*
   * Check input
   */

  if(line == NULL) {
    fprintf(stderr,"ERROR: n_cols.  Null line passed to function.\n");
    return 0;
  }

  /*
   * Initialize pointer
   */

  cptr = line;

  /*
   * Step through line until end
   */

  while(*cptr != '\n' && *cptr != '\0') {

  /*
   * Step through any white space
   *  NB: Does not include tab character for now
   */

    while(*cptr == ' ')
      cptr++;

    if(*cptr != '\n' && *cptr != '\0') {
      ncols++;
      while(*cptr != ' ' && *cptr != '\n' && *cptr != '\0')
	cptr++;
    }
  }

  /*
   * Return the number of columns
   */

  if(verbose)
    printf("%d columns found.\n",ncols);
  if(ncols == 0)
    fprintf(stderr,"** Warning.  n_cols:  No valid cols in input file\n");

  return ncols;
}

/*.......................................................................,
 *
 * Function read_datastruct
 *
 * Reads a 2 or 3 column input file and puts results into a Datastruct array,
 *  which has its memory allocation performed in the function.
 * NB:  This function requires input in the form %f %f %f %d %s
 *
 * Inputs: char *inname        name of input file
 *         char comment        comment character
 *         int *nlines         number of lines in input file -- set by
 *         int format          format of input file
 *                              1: %f %f %f %d %s
 *                              2: %d %f %f %f %s
 *                              3: %f %f %f
 *                              this function.
 *         
 * Output: Datastruct *newdata filled array
 *
 * v07Jun2006 CDF, Added a new format for reading files of form x y z.
 */

Datastruct *read_datastruct(char *inname, char comment, int *nlines,
			    int format)
{
  int no_error=1;           /* Flag set to 0 on error */
  int ncols;                /* Number of columns in input file */
  char line[MAXC];          /* General input string */
  Datastruct *newdata=NULL; /* Filled datastruct array */
  Datastruct *dptr;         /* Pointer to navigate datastruct */
  FILE *ifp=NULL;           /* Input file pointer */

  /*
   * Check the input format
   */

  if(format<1 || format>3) {
    fprintf(stderr,"ERROR: read_datastruct.  Format %d is not supported.\n",
	    format);
    return NULL;
  }

  /*
   * Open input file
   */

  if(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: read_datastruct.\n");
    return NULL;
  }

  /*
   * Get number of lines in input file
   */

  if((*nlines = n_lines(ifp,comment)) == 0) {
    fprintf(stderr,"ERROR: read_datastruct.  No valid data in input file.\n");
    no_error = 0;
  }
  else {
    rewind(ifp);
  }

  /*
   * Allocate memory for datastruct and point dptr at the beginning of
   *  the array.
   */

  if(no_error) {
    if(!(newdata = new_datastruct(*nlines)))
      no_error = 0;
    else
      dptr = newdata;
  }

  /*
   * Read in data
   */

  while(no_error && fgets(line,MAXC,ifp) != NULL) {
    if(line[0] != comment) {
      switch(format) {
      case 1:
	ncols = sscanf(line,"%lf %lf %lf %d %s",&dptr->x,&dptr->y,&dptr->z,
		       &dptr->dataflag,dptr->text);
	if(ncols<1 || ncols>5) {
	  fprintf(stderr,"ERROR: read_datastruct.  Bad input format in %s.\n",
		  inname);
	  fprintf(stderr," For input format = %d,",format);
	  fprintf(stderr," data must be in format %%f %%f %%f %%d %%s\n");
	  fprintf(stderr," i.e., float, float, float, integer, string\n");
	  no_error = 0;
	}
	else
	  dptr++;
	break;
      case 2:
	ncols = sscanf(line,"%d %lf %lf %lf %s",
		       &dptr->dataflag,&dptr->x,&dptr->y,&dptr->z,
		       dptr->text);
 	if(ncols<1 || ncols>5) {
	  fprintf(stderr,"ERROR: read_datastruct.  Bad input format in %s.\n",
		  inname);
	  fprintf(stderr," For input format = %d,",format);
	  fprintf(stderr," data must be in format %%f %%f %%f %%d %%s\n");
	  fprintf(stderr," i.e., float, float, float, integer, string\n");
	  no_error = 0;
	}
	else
	  dptr++;
	break;
      case 3:
	ncols = sscanf(line,"%lf %lf %lf",&dptr->x,&dptr->y,&dptr->z);
 	if(ncols<1 || ncols>3) {
	  fprintf(stderr,"ERROR: read_datastruct.  Bad input format in %s.\n",
		  inname);
	  fprintf(stderr," For input format = %d,",format);
	  fprintf(stderr," data must be in format %%f %%f %%f %%d %%s\n");
	  fprintf(stderr," i.e., float, float, float, integer, string\n");
	  no_error = 0;
	}
	else
	  dptr++;
	break;
      default:
	fprintf(stderr,"ERROR: read_datastruct.  Bad format\n");
	no_error = 0;
      }
    }
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

  if(no_error) {
    printf("\nread_datastruct: %s has %d columns and %d lines\n",inname,
	   ncols,*nlines);
    return newdata;
  }
  else {
    fprintf(stderr,"ERROR: read_datastruct.\n");
    return del_datastruct(newdata);
  }
}

/*.......................................................................,
 *
 * Function read_difmap
 *
 * Reads a the output file from a run of difmap model fitting and puts
 *  the results into a double array, which has its memory allocation 
 *  performed in the function.
 *
 * Inputs: char *inname        name of input file
 *         char comment        comment character
 *         int *nlines         number of lines in input file -- set by
 *                              this function.
 *         Skypos *pos0        field center -- set by this function
 *         
 * Output: Secat *newdata      filled array
 *
 */

Secat *read_difmap(char *inname, char comment, int *nlines, Skypos *pos0)
{
  int i;
  int no_error = 1;         /* Flag set to 0 on error */
  int ngauss = 0;           /* Number of Gaussian components */
  int npoint = 0;           /* Number of point-source components */
  int ncomps = 0;           /* Total number of components */
  int nmodval = 0;          /* Number of parameters describing component */
  float ratio;              /* Axis ratio for gaussian components */
  char line[MAXC];          /* General string variable for reading input */
  char *cptr;               /* Pointer to navigate line */
  Secat *newdata=NULL;      /* Filled secat array */
  Secat *sptr;              /* Pointer to navigate secat */
  FILE *ifp=NULL;           /* Input file pointer */

  /*
   * Open input file
   */

  if(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: read_distcalc.\n");
    return NULL;
  }

  /*
   * Get position of map phase center (RA and Dec)
   */

  if(mod_center(ifp,pos0)) {
    return NULL;
  }
  else {
    printf(" Central position: ");
    print_skypos(stdout,*pos0);
    printf("\n");
  }

  /*
   * Get number of lines in input file
   */

  if((*nlines = n_lines(ifp,comment)) == 0) {
    fprintf(stderr,"ERROR: read_distcalc.  No valid data in input file.\n");
    no_error = 0;
  }
  else {
    printf("read_distcalc: Found %d data lines in %s\n",*nlines,inname);
    rewind(ifp);
  }

  /*
   * Allocate memory for secat and point sptr at the beginning of
   *  the array.
   */

  if(no_error) {
    if(!(newdata = new_secat(*nlines)))
      no_error = 0;
    else
      sptr = newdata;
  }

  /*
   * Loop through the input file
   */

  while(fgets(line,MAXC,ifp) != NULL && no_error) {

    if(line[0] != comment) {

      /*
       * Change all "v"s to spaces
       */

      cptr = line;
      for(i=0; i<=MAXC && *cptr != '\n'; i++,cptr++) {
	if (*cptr == 'v')
	  *cptr = ' ';
      }

      /*
       * Read in line, allowing for possibility of point-source components.
       */

      switch((nmodval = 
	      sscanf(line,"%f %lf %lf %f %f %f",&sptr->fauto,
		     &sptr->dpos,&sptr->dpostheta,
		     &sptr->a_im,&ratio,&sptr->theta))) {
      case 6:
	sptr->b_im = sptr->a_im * ratio; /* May not be right */
	ngauss++;
	break;
      case 3:
	npoint++;
	break;
      default:
	fprintf(stderr," ERROR: Bad input format. nmodval = %d\n",nmodval);
	no_error = 0;
      }
  
      /*
       * Convert flux to mJy
       */
      
      sptr->fauto *= 1000.0;
  
      /*
       * Convert r to arcsec
       */

      sptr->dpos /= 1000.0;

      sptr++;
    }
  }
  
  /*
   * Clean up and exit
   */

  ncomps = ngauss + npoint;
  if(ncomps != *nlines) {
    fprintf(stderr,"ERROR: read_difmap.  Number of components does not ");
    fprintf(stderr,"equal number of data lines (%d vs. %d)\n",ncomps,*nlines);
    no_error = 0;
  }
  else {
    printf(" Read %d gaussian and %d point-source components.\n",ngauss,
	   npoint);
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

 if(no_error)
    return newdata;
  else {
    fprintf(stderr,"ERROR: read_difmap\n");
    return NULL;
  }
}

/*.......................................................................,
 *
 * Function read_secat
 *
 * Reads the output catalog from a SExtractor run and puts the results
 *  into a Secat array, which has its memory allocation 
 *  performed in the function.
 *
 * Inputs: char *inname        name of input file
 *         char comment        comment character
 *         int *nlines         number of lines in input file -- set by
 *                              this function.
 *         int format          flag describing format of input file
 *                              (for up-to-date listing, see the
 *                               secat_format function)
 *         
 * Output: Secat *newdata filled array
 *
 * v20Jul2003 CDF, Added read-in of alpha, delta pairs and the
 *                  conversion into RA, Dec pairs.
 * v20Jul2005 CDF, Added a version for the output from the USNO star catalog
 * v06Jan2007 CDF, Added a another input format
 * v19Jun2007 CDF, Added a another input format
 * v09Jul2007 CDF, Changed format 3
 * v2009Jan24 CDF, Modified format 0 to be ID, alpha, delta
 * v2010Dec19 CDF, Added format 17 for photometry
 */

Secat *read_secat(char *inname, char comment, int *nlines, int format)
{
  int no_error=1;           /* Flag set to 0 on error */
  int lc=0;                 /* Running counter for line number */
  int nskip;                /* Nubmer of header lines to skip when reading in */
  int ncols;                /* Number of columns in input file */
  int nexp;                 /* Number of columns expected in input file */
  int count=0;              /* Used to set ID number for catalogs with no IDs */
  int rahr,ramin,decamin;   /* Components of hms format */
  double rasec,decasec;     /* Components of hms format */
  double djunk;             /* Variable used to read in unrequired data */
  char line[MAXC];          /* General input string */
  char outformat[MAXC];     /* Output format string */
  char decdeg[3];           /* Degrees of declination, in hms format */
  Pos tmppos;               /* Temporary position holder */
  Secat *newdata=NULL;      /* Filled secat array */
  Secat *sptr;              /* Pointer to navigate secat */
  FILE *ifp=NULL;           /* Input file pointer */

  printf("read_secat: Input file: %s. Input format = %d\n",inname,format);

  /*
   * If file is in distcalc format, call read_distcalc
   */

  if(format == 5) {
    if(!(newdata = (read_distcalc(inname,'#',nlines,0)))) {
      fprintf(stderr,"ERROR: read_secat");
      return NULL;
    }
    else
      return newdata;
  }

  /*
   * Set up number of expected columns and number of lines to
   *  skip when reading in data
   */

  nskip = 0;
  switch(format) {
  case 0:
    nexp = 3;
    break;
  case 1:
    nexp = 2;
    break;
  case 2:
    nexp = 6;
    break;
  case 3:
    nexp = 8;
    break;
  case 4:
    nskip = 14;
    nexp = 2;
    break;
  case 6:
    nexp = 39;
    break;
  case 7:
    nexp = 48;
    break;
  case 8:
    nskip = 3;
    nexp = 7;
    break;
  case 9:
    nexp = 7;
    break;
  case 10:
    nexp = 3;
    break;
  case 13:
    nexp = 20;
    break;
  case 15:
    nexp = 5;
    break;
  case 16:
    nexp = 8;
    break;
  case 17:
    nexp = 47;
    break;
  case 18:
    nexp = 6;
    break;
  default:
    fprintf(stderr,"ERROR: read_secat. Not a valid format\n");
    return NULL;
  }

  /*
   * Open input file
   */

  if(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: read_secat.\n");
    return NULL;
  }

  /*
   * Get number of lines in input file
   */

  if((*nlines = n_lines(ifp,comment)) == 0) {
    fprintf(stderr,"ERROR: read_secat.  No valid data in input file.\n");
    no_error = 0;
  }
  else {
    rewind(ifp);
    *nlines = *nlines - nskip;
  }


  /*
   * Allocate memory for secat and point sptr at the beginning of
   *  the array.
   */

  if(no_error) {
    if(!(newdata = new_secat(*nlines)))
      no_error = 0;
    else
      sptr = newdata;
  }

  /*
   * Read in data
   */

  while(no_error && fgets(line,MAXC,ifp) != NULL) {
    /* Have to skip the first nskip lines */
    while(lc < nskip) {
      lc++;
      if(fgets(line,MAXC,ifp) == NULL) {
	fprintf(stderr,"ERROR: read_secat. Bad file format for %s\n",
		inname);
	return(del_secat(newdata));
      }
    }
    lc++;
    if(line[0] != comment) {
      switch(format) {
      case 0:
	ncols = sscanf(line,"%lf %lf %lf",&djunk,&sptr->alpha,&sptr->delta);
	sptr->id = (int) djunk;
	sprintf(sptr->name,"%d",sptr->id);
	break;
      case 1:
	ncols = sscanf(line,"%lf %lf",&sptr->alpha,&sptr->delta);
	count++;
	sptr->id = count; 
	sprintf(sptr->name,"%d",count);
	break;
      case 2:
	ncols = 
	  sscanf(line,"%d %d %lf %d %d %lf",
		 &sptr->skypos.hr,&sptr->skypos.min,&sptr->skypos.sec,
		 &sptr->skypos.deg,&sptr->skypos.amin,&sptr->skypos.asec);
	count++;
	sptr->id = count; 
	sprintf(sptr->name,"%d",count);

	/*
	 * Convert RA,Dec to alpha,delta 
	 */

	spos2deg(sptr->skypos,&sptr->alpha,&sptr->delta);
	break;
      case 3:
	ncols = 
	  sscanf(line,"%lf %lf %d %d %lf %d %d %lf",
		 &sptr->x,&sptr->y,
		 &sptr->skypos.hr,&sptr->skypos.min,&sptr->skypos.sec,
		 &sptr->skypos.deg,&sptr->skypos.amin,&sptr->skypos.asec);

	/*
	 * Convert RA,Dec to alpha,delta 
	 */

	spos2deg(sptr->skypos,&sptr->alpha,&sptr->delta);
	break;
      case 4:
	ncols = sscanf(line,"%lf %lf",&sptr->alpha,&sptr->delta);
	count++;
	sptr->id = count; 
	sprintf(sptr->name,"%d",sptr->id);

	/*
	 * Convert alpha,delta to RA,Dec
	 */

	deg2spos(sptr->alpha,sptr->delta,&sptr->skypos);
	break;
      case 6:
	ncols = 
	  sscanf(line,"%d %lf %lf %d %f %f %f %g %g %f %f %g %g %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %lf %lf",
		 &sptr->id,&sptr->x,&sptr->y,&sptr->fitflag,&sptr->class,
		 &sptr->miso,&sptr->misoerr,&sptr->fiso,&sptr->fisoerr,
		 &sptr->mtot,&sptr->merr,&sptr->fauto,&sptr->fautoerr,
		 &sptr->r_kron,&sptr->bkgd,&sptr->thresh,&sptr->muthresh,
	         &sptr->isoarea,&sptr->a_im,
		 &sptr->b_im,&sptr->theta,&sptr->fwhm,
		 &sptr->maper[0],&sptr->maper[1],&sptr->maper[2],
		 &sptr->mapererr[0],&sptr->mapererr[1],&sptr->mapererr[2],
		 &sptr->faper[0],&sptr->faper[1],&sptr->faper[2],
		 &sptr->fapererr[0],&sptr->fapererr[1],&sptr->fapererr[2],
		 &sptr->r2,&sptr->r5,&sptr->r8,&sptr->alpha,&sptr->delta);
	sprintf(sptr->name,"%d",sptr->id);

	/*
	 * Convert alpha,delta to RA,Dec
	 */

	deg2spos(sptr->alpha,sptr->delta,&sptr->skypos);
	break;
      case 7:
	ncols = 
	  sscanf(line,"%d %lf %lf %d %f %f %f %g %g %f %f %g %g %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %lf %lf %lf %lf %lf %d %d %lf %d %d %lf",
		 &sptr->id,&sptr->x,&sptr->y,&sptr->fitflag,&sptr->class,
		 &sptr->miso,&sptr->misoerr,&sptr->fiso,&sptr->fisoerr,
		 &sptr->mtot,&sptr->merr,&sptr->fauto,&sptr->fautoerr,
		 &sptr->r_kron,&sptr->bkgd,&sptr->thresh,&sptr->muthresh,
	         &sptr->isoarea,&sptr->a_im,
		 &sptr->b_im,&sptr->theta,&sptr->fwhm,
		 &sptr->ma1,&sptr->ma2,&sptr->ma3,
		 &sptr->ma1err,&sptr->ma2err,&sptr->ma3err,
		 &sptr->fa1,&sptr->fa2,&sptr->fa3,
		 &sptr->fa1err,&sptr->fa2err,&sptr->fa3err,
		 &sptr->r2,&sptr->r5,&sptr->r8,&sptr->alpha,&sptr->delta,
		 &sptr->dx,&sptr->dy,&sptr->dpos,
		 &sptr->skypos.hr,&sptr->skypos.min,&sptr->skypos.sec,
		 &sptr->skypos.deg,&sptr->skypos.amin,&sptr->skypos.asec);
	sprintf(sptr->name,"%d",sptr->id);

	/*
	 * Convert dx,dy to r,theta
	 */

	tmppos.x = -sptr->dx;
	tmppos.y = sptr->dy;
	xy2rth(tmppos,&djunk,&sptr->dpostheta,1);
	break;
      case 8:
	ncols = 
	  sscanf(line,"%d %d %lf %d %d %lf %f",
		 &sptr->skypos.hr,&sptr->skypos.min,&sptr->skypos.sec,
		 &sptr->skypos.deg,&sptr->skypos.amin,&sptr->skypos.asec,
		 &sptr->mtot);
	sptr->id = count + 1; 
	sprintf(sptr->name,"%d",sptr->id);
	count++;
	/*
	 * Convert RA,Dec to alpha,delta 
	 */

	spos2deg(sptr->skypos,&sptr->alpha,&sptr->delta);
	break;
      case 9:
	ncols = 
	  sscanf(line,"%s %d %d %lf %d %d %lf %f %f",
		 sptr->name,
		 &sptr->skypos.hr,&sptr->skypos.min,&sptr->skypos.sec,
		 &sptr->skypos.deg,&sptr->skypos.amin,&sptr->skypos.asec,
		 &sptr->mtot,&sptr->merr);
	sptr->id = count + 1; 
	count++;
	/*
	 * Convert RA,Dec to alpha,delta 
	 */

	spos2deg(sptr->skypos,&sptr->alpha,&sptr->delta);
	break;
      case 10:
	ncols = 
	  sscanf(line,"%s %lf %lf %f %f",
		 sptr->name,&sptr->alpha,&sptr->delta,&sptr->mtot,&sptr->merr);
	sptr->id = count + 1; 
	count++;
	/*
	 * Convert alpha,delta to RA,Dec
	 */
	deg2spos(sptr->alpha,sptr->delta,&sptr->skypos);
	break;
      case 13:
	ncols = 
	  sscanf(line,
		 "%d %lf %lf %lf %lf %d %f %f %f %g %g %f %f %f %f %f %f %f %f %f",
		 &sptr->id,&sptr->alpha,&sptr->delta,
		 &sptr->x,&sptr->y,&sptr->fitflag,&sptr->class,
		 &sptr->mtot,&sptr->merr,&sptr->fauto,&sptr->fautoerr,
		 &sptr->r_kron,&sptr->a_im,&sptr->b_im,&sptr->theta,&sptr->fwhm,
		 &sptr->bkgd,&sptr->thresh,&sptr->muthresh,&sptr->isoarea);
	sprintf(sptr->name,"%d",sptr->id);

	/*
	 * Convert alpha,delta to RA,Dec
	 */

	deg2spos(sptr->alpha,sptr->delta,&sptr->skypos);
	break;
      case 15:
	ncols = 
	  sscanf(line,"%s %lf %lf %f %f",
		 sptr->name,&sptr->alpha,&sptr->delta,&sptr->zspec,
		 &sptr->zspecerr);
	sptr->id = count + 1; 
	count++;
	/*
	 * Convert alpha,delta to RA,Dec
	 */

	deg2spos(sptr->alpha,sptr->delta,&sptr->skypos);
	break;
      case 16:
	ncols = 
	  sscanf(line,"%s %d %d %lf %d %d %lf %f",
		 sptr->name,
		 &sptr->skypos.hr,&sptr->skypos.min,&sptr->skypos.sec,
		 &sptr->skypos.deg,&sptr->skypos.amin,&sptr->skypos.asec,
		 &sptr->zspec);
	sptr->id = count + 1; 
	count++;
	/*
	 * Convert RA,Dec to alpha,delta 
	 */

	spos2deg(sptr->skypos,&sptr->alpha,&sptr->delta);
	break;
      case 17:
	ncols = 
	  sscanf(line,"%d %lf %lf %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %f",
		 &sptr->id,&sptr->x,&sptr->y,&sptr->mtot,&sptr->merr,
		 &sptr->maper[0],&sptr->maper[1],&sptr->maper[2],
		 &sptr->maper[3],&sptr->maper[4],&sptr->maper[5],
		 &sptr->maper[6],&sptr->maper[7],&sptr->maper[8],
		 &sptr->maper[9],&sptr->maper[10],&sptr->maper[11],
		 &sptr->maper[12],&sptr->maper[13],&sptr->maper[14],
		 &sptr->maper[15],&sptr->maper[16],&sptr->maper[17],
		 &sptr->maper[18],&sptr->maper[19],
		 &sptr->mapererr[0],&sptr->mapererr[1],&sptr->mapererr[2],
		 &sptr->mapererr[3],&sptr->mapererr[4],&sptr->mapererr[5],
		 &sptr->mapererr[6],&sptr->mapererr[7],&sptr->mapererr[8],
		 &sptr->mapererr[9],&sptr->mapererr[10],&sptr->mapererr[11],
		 &sptr->mapererr[12],&sptr->mapererr[13],&sptr->mapererr[14],
		 &sptr->mapererr[15],&sptr->mapererr[16],&sptr->mapererr[17],
		 &sptr->mapererr[18],&sptr->mapererr[19],
		 &sptr->fitflag,&sptr->class);
	sprintf(sptr->name,"%d",sptr->id);
	break;
      case 18:
	ncols = sscanf(line,"%lf %lf %f %f %f %f",
		       &sptr->x,&sptr->y,&sptr->a_im,&sptr->b_im,&sptr->theta,
		       &sptr->r_kron);
	break;
      }
      if(ncols < nexp) {
	fprintf(stderr,
		"ERROR: read_secat.  Bad input format in %s. (line = %d)\n",
		inname,lc);
	fprintf(stderr," File must contain at least %d columns.",nexp);
	fprintf(stderr," --  it contained %d columns.\n",ncols);
	no_error = 0;
      }
      else
	sptr++;
    }
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

  if(no_error) {
    printf("read_secat: %s has %d columns and %d lines\n",inname,
	   ncols,*nlines);
    return newdata;
  }
  else {
    fprintf(stderr,"ERROR: read_secat.\n");
    return del_secat(newdata);
  }
}

/*.......................................................................,
 *
 * Function read_secat2
 *
 * Reads the output catalog from a SExtractor run and puts the results
 *  into a Secat array, which has its memory allocation 
 *  performed in the function.
 *
 * Inputs: char *inname        name of input file
 *         char comment        comment character
 *         int *nlines         number of lines in input file -- set by
 *                              this function.
 *         
 * Output: Secat *newdata filled array
 *
 * v23Apr2006 CDF, First attempt to generalize read_secat
 */

Secat *read_secat2(char *inname, char comment, int *nlines)
{
  int no_error=1;           /* Flag set to 0 on error */
  int lc=0;                 /* Running counter for line number */
  int nskip;                /* Nubmer of header lines to skip when reading in */
  int ncols;                /* Number of columns in input file */
  int nexp;                 /* Number of columns expected in input file */
  int count=0;              /* Used to set ID number for catalogs with no IDs */
  double djunk;             /* Variable used to read in unrequired data */
  char line[MAXC];          /* General input string */
  char outformat[MAXC];     /* Output format string */
  Pos tmppos;               /* Temporary position holder */
  Secat *newdata=NULL;      /* Filled secat array */
  Secat *sptr;              /* Pointer to navigate secat */
  FILE *ifp=NULL;           /* Input file pointer */

  /*
   * Open input file
   */

  if(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: read_secat.\n");
    return NULL;
  }

  /*
   * Get number of lines in input file
   */

  if((*nlines = n_lines(ifp,comment)) == 0) {
    fprintf(stderr,"ERROR: read_secat.  No valid data in input file.\n");
    no_error = 0;
  }
  else {
    rewind(ifp);
  }
#if 0
  /*
   * Read in header information
   */

  while(no_error && fgets(line,MAXC,ifp) != NULL) {
    if(line[0] == comment) {
      if(sscanf(line,"%s %d %s",junk,&col,keyword) != 3) {
	fprintf(stderr,"ERROR: read_secat.  Bad format in header for %s\n",
		inname);
	no_error = 0;
      }
    }
  }


  /*
   * Allocate memory for secat and point sptr at the beginning of
   *  the array.
   */

  if(no_error) {
    if(!(newdata = new_secat(*nlines)))
      no_error = 0;
    else
      sptr = newdata;
  }
#endif
  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

  if(no_error) {
    printf("read_secat2: %s has %d columns and %d lines\n",inname,
	   ncols,*nlines);
    return newdata;
  }
  else {
    fprintf(stderr,"ERROR: read_secat2.\n");
    return del_secat(newdata);
  }
}

/*.......................................................................
 *
 * Function read_secat_header
 *
 * Extracts the column definitions from the header information in a
 *  SExtractor catalog
 */

SHeader *read_secat_header(FILE *ifp, int *nhead)
{
  SHeader *sheader=NULL;

  return sheader;
}

/*.......................................................................,
 *
 * Function write_secat
 *
 * Writes an array of Secat structures to an output file.
 *
 * Inputs: Secat *secat        secat array
 *         int ncat            number of catalog members
 *         char *outname       name of output file
 *         int format          flag describing format of input file
 *                              (see secat_format for format descriptions)
 *         
 * Output: int (0 or 1)       0 ==> OK.  1 ==> error.
 *
 * v2005Sep08 CDF, Added error checking on format code.
 *                 Added output format for 2MASS format processed by catsort
 *                  -- write out in USNO format (format 8)
 * v2007Jul09 CDF, Added format 3
 * v2007Jul29 CDF, Added format 14
 * v2007Aug02 CDF, Added format 15
 * v2009Jan24 CDF, Modified format 0 to be ID, alpha, delta
 *
 */

int write_secat(Secat *secat, int ncat, char *outname, int format)
{
  int i;            /* Looping variable */
  int no_error=1;   /* Flag set to 0 on error */
  Secat *sptr;      /* Pointer to navigate secat */
  FILE *ofp=NULL;   /* Output file pointer */


  /*
   * Open output file
   */

  if(!(ofp = open_writefile(outname))) {
    fprintf(stderr,"ERROR: write_secat.\n");
    return 1;
  }

  /*
   * Add header for format 4 --> format 8
   */

  if(format == 4)
    fprintf(ofp,".\n.\n.\n");

  /*
   * Write out data
   */

  for(i=0,sptr=secat; i<ncat; i++,sptr++) {
    if(no_error) {
      switch(format) {
      case 0:
	fprintf(ofp,"%19s %11.7f %+11.7f %02d %02d %08.5f %+03d %02d %07.4f\n",
		sptr->name,sptr->alpha,sptr->delta,
		sptr->skypos.hr,sptr->skypos.min,sptr->skypos.sec,
		sptr->skypos.deg,sptr->skypos.amin,sptr->skypos.asec);
	break;
      case 1:
	fprintf(ofp,
		"%-10s %02d %02d %08.5f %+03d %02d %07.4f ",
		sptr->name,sptr->skypos.hr,sptr->skypos.min,sptr->skypos.sec,
		sptr->skypos.deg,sptr->skypos.amin,sptr->skypos.asec);
	fprintf(ofp,"%12.8f %+12.8f %8.2f %8.2f %8.2f\n",
		sptr->alpha,sptr->delta,sptr->dx,sptr->dy,sptr->dpos);
	break;
      case 2:
	fprintf(ofp,"%04d %8.2f %8.2f %2d %4.2f %7.3f %7.3f %7.3f %7.3f ",
		sptr->id,sptr->x,sptr->y,sptr->fitflag,sptr->class,
		sptr->miso,sptr->misoerr,sptr->mtot,sptr->merr);
	fprintf(ofp,"%6.2f %8.2f %9.4f %7.3f %5.0f %8.3f %8.3f %6.1f %8.3f ",
		sptr->r_kron,sptr->bkgd,sptr->thresh,sptr->muthresh,
		sptr->isoarea,sptr->a_im,sptr->b_im,sptr->theta,sptr->fwhm);
	fprintf(ofp,"%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f ",
		sptr->ma1,sptr->ma2,sptr->ma3,
		sptr->ma1err,sptr->ma2err,sptr->ma3err);
	fprintf(ofp,"%9.2f %9.2f %9.2f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f ",
		sptr->fa1,sptr->fa2,sptr->fa3,
		sptr->fa1err,sptr->fa2err,sptr->fa3err,
		sptr->r2,sptr->r5,sptr->r8);
	fprintf(ofp,"%11.7f %+11.7f %8.2f %8.2f %8.2f ",
		sptr->alpha,sptr->delta,sptr->dx,sptr->dy,sptr->dpos);
	fprintf(ofp,"%02d %02d %08.5f %+03d %02d %07.4f\n",
		sptr->skypos.hr,sptr->skypos.min,sptr->skypos.sec,
		sptr->skypos.deg,sptr->skypos.amin,sptr->skypos.asec);
	break;
      case 3:
	fprintf(ofp,
		"%7.2f %7.2f %02d %02d %08.5f %+03d %02d %07.4f %11.7f %+11.7f",
		sptr->x,sptr->y,
		sptr->skypos.hr,sptr->skypos.min,sptr->skypos.sec,
		sptr->skypos.deg,sptr->skypos.amin,sptr->skypos.asec,
		sptr->alpha,sptr->delta);
	fprintf(ofp,"%s %05d  %6.2f\n",sptr->name,sptr->id,sptr->fwhm);
	break;
      case 4:
	fprintf(ofp,"%02d %02d %08.5f %+03d %02d %07.4f 0.00 %11.7f %+11.7f\n",
		sptr->skypos.hr,sptr->skypos.min,sptr->skypos.sec,
		sptr->skypos.deg,sptr->skypos.amin,sptr->skypos.asec,
		sptr->alpha,sptr->delta);
	break;
      case 6:
	fprintf(ofp,"%04d %8.2f %8.2f %2d %4.2f %7.3f %7.3f %10g %10g ",
		sptr->id,sptr->x,sptr->y,sptr->fitflag,sptr->class,
		sptr->miso,sptr->misoerr,sptr->fiso,sptr->fisoerr);
	fprintf(ofp,"%7.3f %7.3f %10g %10g ",sptr->mtot,sptr->merr,
		sptr->fauto,sptr->fautoerr);
	fprintf(ofp,"%6.2f %8.2f %9.4f %7.3f %5.0f %8.3f %8.3f %6.1f %8.3f ",
		sptr->r_kron,sptr->bkgd,sptr->thresh,sptr->muthresh,
		sptr->isoarea,sptr->a_im,sptr->b_im,sptr->theta,sptr->fwhm);
	fprintf(ofp,"%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f ",
		sptr->ma1,sptr->ma2,sptr->ma3,
		sptr->ma1err,sptr->ma2err,sptr->ma3err);
	fprintf(ofp,"%9.2f %9.2f %9.2f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f ",
		sptr->fa1,sptr->fa2,sptr->fa3,
		sptr->fa1err,sptr->fa2err,sptr->fa3err,
		sptr->r2,sptr->r5,sptr->r8);
	fprintf(ofp,"%11.7f %+11.7f %8.2f %8.2f %8.2f ",
		sptr->alpha,sptr->delta,sptr->dx,sptr->dy,sptr->dpos);
	/*
	fprintf(ofp,"%02d %02d %08.5f %+03d %02d %07.4f",
		sptr->skypos.hr,sptr->skypos.min,sptr->skypos.sec,
		sptr->skypos.deg,sptr->skypos.amin,sptr->skypos.asec);
	*/
	fprintf(ofp,"\n");
	break;
      case 9: case 10:
	fprintf(ofp,"%s %11.7f %+11.7f  %5.2f %5.2f  ",
		sptr->name,sptr->alpha,sptr->delta,sptr->mtot,sptr->merr);
	fprintf(ofp,"%02d %02d %08.5f %+03d %02d %07.4f  ",
		sptr->skypos.hr,sptr->skypos.min,sptr->skypos.sec,
		sptr->skypos.deg,sptr->skypos.amin,sptr->skypos.asec);
	fprintf(ofp,"%8.2f %8.2f %8.2f\n",
		sptr->dx,sptr->dy,sptr->dpos);
	break;
      case 14:
	fprintf(ofp,
		"%7.2f %7.2f %02d:%02d:%07.4f %+03d:%02d:%06.3f %11.7f %+11.7f",
		sptr->x,sptr->y,
		sptr->skypos.hr,sptr->skypos.min,sptr->skypos.sec,
		sptr->skypos.deg,sptr->skypos.amin,sptr->skypos.asec,
		sptr->alpha,sptr->delta);
	fprintf(ofp," %s %6.2f\n",sptr->name,sptr->fwhm);
	break;
#if 0
      case 15:
	fprintf(ofp,"05d %11.7f %+11.7f %7.4f %6.4f\n",
		sptr->id,sptr->alpha,sptr->delta,sptr->zspec,sptr->zspecerr);
	break;
#endif
      default:
	fprintf(stderr,"ERROR: write_secat.  ");
	fprintf(stderr,"Format %d is not valid for this function.\n",format);
	no_error = 0;
      }
    }
  }

  /*
   * Clean up and exit
   */

  if(ofp)
    fclose(ofp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: write_secat.\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function secat_format
 *
 * Prints out the acceptable formats for a secat file.
 *
 * Inputs: (none)
 *
 * Output: (none)
 *
 */

void secat_format()
{
  fprintf(stderr,"     0 ==> ID, alpha, delta\n");
  fprintf(stderr,"     1 ==> alpha, delta\n");
  fprintf(stderr,"     2 ==> hr, min, sec, deg, amin, asec\n");
  fprintf(stderr,"     3 ==> x, y, hr, min, sec, deg, amin, asec\n");
  fprintf(stderr,"     4 ==> Format returned by 2MASS server\n");
  fprintf(stderr,"        alpha, delta, ...many other cols...\n");
  fprintf(stderr,"        (Skip the first 14 lines)\n");
  fprintf(stderr,"     5 ==> distcalc.c output format\n");
  fprintf(stderr,"     6 ==> Format produced by run_sext*sh in Lenses ");
  fprintf(stderr,"directory.\n");
  fprintf(stderr,"     7 ==> Format produced by running catsort on format=6\n");
  fprintf(stderr,"     8 ==> Format returned by USNO A2 star web page\n");
  fprintf(stderr,"        hr, min, sec, deg, amin, asec, mag\n");
  fprintf(stderr,"        (Skip the first 3 lines)\n");
  fprintf(stderr,"     9 ==> ");
  fprintf(stderr,"name, hr, min, sec, deg, amin, asec, mag, magerr\n");
  fprintf(stderr,"        (where mag and magerr are optional columns)\n");
  fprintf(stderr,"        e.g., format returned by USNO B1 star web page\n");
  fprintf(stderr,"    10 ==> name, alpha, delta, mag, magerr\n");
  fprintf(stderr,"        (where mag and magerr are optional columns)\n");
  fprintf(stderr,"        e.g., from SDSS object-crossid page\n");
  fprintf(stderr,"    11 ==> SDSS format from Imaging Query page:\n");
  fprintf(stderr,"        alpha, delta, u, g, r, i, z, ... perhaps more\n");
  fprintf(stderr,"    12 ==> modified SDSS format:\n");
  fprintf(stderr,"        id, alpha, delta, u, g, r, i, z, dalpha, ddelta, ");
  fprintf(stderr,"dtot, RA(hms) Dec(hms)\n");
  fprintf(stderr,"    13 ==> id, alpha, delta, x, y, flag, class, ... \n");
  fprintf(stderr,"    14 ==> x, y, hr:min:sec, deg:amin:asec\n");
  fprintf(stderr,"    15 ==> master tab format:\n");
  fprintf(stderr,"        name, alpha, delta, z, zerr\n");
  fprintf(stderr,"    16 ==> name, hr, min, sec, deg, amin, asec, z\n");
  fprintf(stderr,"    17 ==> standard star photometry format:\n");
  fprintf(stderr,"        id, x, y, mauto, errauto, (20 mag_apers), ");
  fprintf(stderr,"(20 mag_aper errors)\n");
  fprintf(stderr,"\n");
}

/*.......................................................................,
 *
 * Function read_distcalc
 *
 * Reads a the output file from a distcalc run and puts the results
 *  into a Secat array, which has its memory allocation 
 *  performed in the function.
*
 * Inputs: char *inname        name of input file
 *         char comment        comment character
 *         int *nlines         number of lines in input file -- set by
 *                              this function.
 *         int format          format of input file:
 *                              0 ==> format of input to distcalc.c
 *                              1 ==> format of output of distcalc.c
 *         
 * Output: Secat *newdata filled array
 *
 * v07Aug03 CDF, Added a format flag to enable reading
 *                of either input to or output of distcalc.c.
 */

Secat *read_distcalc(char *inname, char comment, int *nlines, int format)
{
  int no_error=1;           /* Flag set to 0 on error */
  int lc=0;                 /* Running counter for line number */
  int ncols;                /* Number of columns in input file */
  int nexp;                 /* Number of columns expected in input file */
  double djunk;             /* Variable used to read in unrequired data */
  char line[MAXC];          /* General input string */
  char ew,ns;               /* Characters indicating offset direction */
  Pos tmppos;               /* Temporary position holder */
  Secat *newdata=NULL;      /* Filled secat array */
  Secat *sptr;              /* Pointer to navigate secat */
  FILE *ifp=NULL;           /* Input file pointer */

  /*
   * Open input file
   */

  if(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: read_distcalc.\n");
    return NULL;
  }

  /*
   * Get number of lines in input file
   */

  if((*nlines = n_lines(ifp,comment)) == 0) {
    fprintf(stderr,"ERROR: read_distcalc.  No valid data in input file.\n");
    no_error = 0;
  }
  else {
    printf("read_distcalc: Found %d data lines in %s\n",*nlines,inname);
    rewind(ifp);
  }

  /*
   * Allocate memory for secat and point sptr at the beginning of
   *  the array.
   */

  if(no_error) {
    if(!(newdata = new_secat(*nlines)))
      no_error = 0;
    else
      sptr = newdata;
  }

  /*
   * Set up number of expected columns
   */

  if(format)
    nexp = 12;
  else
    nexp = 7;

  /*
   * Read in data
   */

  while(no_error && fgets(line,MAXC,ifp) != NULL) {
    lc++;
    if(line[0] != comment) {
      if(format) {
	ncols = 
	  sscanf(line,"%s %d %d %lf %d %d %lf %lf %c %lf %c %lf",sptr->name,
		 &sptr->skypos.hr,&sptr->skypos.min,&sptr->skypos.sec,
		 &sptr->skypos.deg,&sptr->skypos.amin,&sptr->skypos.asec,
		 &sptr->dx,&ew,&sptr->dy,&ns,&sptr->dpos);
      }
      else {
	ncols = 
	  sscanf(line,"%s %d %d %lf %d %d %lf",sptr->name,
		 &sptr->skypos.hr,&sptr->skypos.min,&sptr->skypos.sec,
		 &sptr->skypos.deg,&sptr->skypos.amin,&sptr->skypos.asec);
      }
      if(ncols != nexp) {
	fprintf(stderr,
		"ERROR: read_distcalc.  Bad input format in %s. (line = %d)\n",
		inname,lc);
	fprintf(stderr," File must contain at least %d columns.",nexp);
	fprintf(stderr," --  it contained %d columns.\n",ncols);
	no_error = 0;
      }
      else {
	/*
	 * Convert RA,Dec to alpha,delta 
	 */

	spos2deg(sptr->skypos,&sptr->alpha,&sptr->delta);

	if(format) {
	  /* Astronomy convention */
	  if(ew == 'W') 
	    sptr->dx = -sptr->dx;
	  if(ns == 'S') 
	    sptr->dy = -sptr->dy;
	  tmppos.x = -sptr->dx;
	  tmppos.y = sptr->dy;
	  xy2rth(tmppos,&djunk,&sptr->dpostheta,1);
	}
	sptr++;
      }
    }
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

  if(no_error) {
    printf("read_distcalc: %s has %d columns and %d lines\n",inname,
	   ncols,*nlines);
    return newdata;
  }
  else {
    fprintf(stderr,"ERROR: read_distcalc.\n");
    return del_secat(newdata);
  }
}

/*.......................................................................
 *
 * Function read_sdss
 *
 * Reads a catalog generated by the SDSS Imaging Query website
 *  into an array of sdsscat containers.
 *
 * The expected format of the input file is
 *  alpha, delta, u, g, r, i, z, class*, ..spectroscopic parameters
 *  *NB: The class column is not actually part of the format
 *   returned by the SDSS web site, but is added by an awk script
 *   before this function is called.
 *
 * Inputs: char *inname        name of input file
 *         char comment        comment character
 *         int *nlines         number of lines in input file -- set by
 *                              this function.
 *         
 * Output: SDSScat *newdata    filled array
 */

SDSScat *read_sdss(char *inname, char comment, int *nlines, int format)
{
  int no_error=1;           /* Flag set to 0 on error */
  int lc=0;                 /* Running counter for line number */
  int nskip=0;              /* Nubmer of header lines to skip when reading in */
  int ncols;                /* Number of columns in input file */
  int nexp;                 /* Number of columns expected in input file */
  int count=0;              /* Used to set ID number for catalogs with no IDs */
  double djunk;             /* Variable used to read in unrequired data */
  char line[MAXC];          /* General input string */
  char outformat[MAXC];     /* Output format string */
  Pos tmppos;               /* Temporary position holder */
  SDSScat *newdata=NULL;    /* Filled secat array */
  SDSScat *sptr;            /* Pointer to navigate secat */
  FILE *ifp=NULL;           /* Input file pointer */

  /*
   * Initialize the number of expected columns
   */

  nexp = 8;

  /*
   * Open input file
   */

  if(!(ifp = open_readfile(inname))) {
    fprintf(stderr,"ERROR: read_sdss.\n");
    return NULL;
  }

  /*
   * Get number of lines in input file
   */

  if((*nlines = n_lines(ifp,comment)) == 0) {
    fprintf(stderr,"ERROR: read_sdss.  No valid data in input file.\n");
    no_error = 0;
  }
  else {
    rewind(ifp);
  }

  /*
   * Allocate memory for SDSScat and point sptr at the beginning of
   *  the array.
   */

  if(no_error) {
    if(!(newdata = new_sdsscat(*nlines)))
      no_error = 0;
    else
      sptr = newdata;
  }

  /*
   * Read in data
   */

  while(no_error && fgets(line,MAXC,ifp) != NULL) {
    if(line[0] != comment) {
      ncols = sscanf(line,"%lf %lf %f %f %f %f %f %f",
		     &sptr->alpha,&sptr->delta,
		     &sptr->u,&sptr->g,&sptr->r,&sptr->i,&sptr->z,
		     &sptr->class);
      if(ncols != nexp) {
	fprintf(stderr,
		"ERROR: read_sdss.  Bad input format in %s. (line = %d)\n",
		inname,lc);
	fprintf(stderr," File must contain at least %d columns.",nexp);
	fprintf(stderr," --  it contained %d columns.\n",ncols);
	no_error = 0;
      }

      /*
       * Convert alpha,delta to RA,Dec
       */

      deg2spos(sptr->alpha,sptr->delta,&sptr->skypos);
      sptr++;
    }
  }

  /*
   * Clean up and exit
   */

  if(ifp)
    fclose(ifp);

  if(no_error) {
    printf("read_sdss: %s has %d columns and %d lines\n",inname,
	   ncols,*nlines);
    return newdata;
  }
  else {
    fprintf(stderr,"ERROR: read_secat.\n");
    return del_sdsscat(newdata);
  }
}

/*.......................................................................,
 *
 * Function write_sdss
 *
 * Writes an array of SDSScat structures to an output file.
 *
 * Inputs: SDSScat *scat       SDSScat array
 *         int ncat            number of catalog members
 *         char *outname       name of output file
 *         int format          flag describing format of input file
 *                              0 ==> ID, x, y
 *                              1 ==> Format produced by run_sext.sh in
 *                                 Lenses/SExtractor directory
 *                              2 ==> Temporary run_sext.sh format during
 *                                 format change
 *                              3 ==> Format produced by matchcat.c
 *                              4 ==> Specialized format from running
 *                                 catsort on a file from the 2MASS server
 *                              6 ==> Another hardwired SExtractor format
 *         
 * Output: int (0 or 1)       0 ==> OK.  1 ==> error.
 *
 *
 */

int write_sdss(SDSScat *scat, int ncat, char *outname, int format)
{
  int i;            /* Looping variable */
  int no_error=1;   /* Flag set to 0 on error */
  SDSScat *sptr;    /* Pointer to navigate scat */
  FILE *ofp=NULL;   /* Output file pointer */


  /*
   * Open output file
   */

  if(!(ofp = open_writefile(outname))) {
    fprintf(stderr,"ERROR: write_sdss.\n");
    return 1;
  }

  /*
   * Write out data -- right now there are no differences between
   *  the two formats
   */

  for(i=0,sptr=scat; i<ncat; i++,sptr++) {
    if(no_error) {
      fprintf(ofp,"%04d %12.8f %+12.8f   %6.3f %6.3f %6.3f %6.3f %6.3f   ",
	      sptr->id,sptr->alpha,sptr->delta,
	      sptr->u,sptr->g,sptr->r,sptr->i,sptr->z);
      fprintf(ofp,"%4.2f %8.2f %8.2f %8.2f    ",
	      sptr->class,sptr->dx,sptr->dy,sptr->dpos);
      fprintf(ofp,"%02d %02d %08.5f %+03d %02d %07.4f",
	      sptr->skypos.hr,sptr->skypos.min,sptr->skypos.sec,
	      sptr->skypos.deg,sptr->skypos.amin,sptr->skypos.asec);
      /* At a later time add spectroscopic info */
      fprintf(ofp,"\n");
    }
  }

  /*
   * Clean up and exit
   */

  if(ofp)
    fclose(ofp);

  if(no_error)
    return 0;
  else {
    fprintf(stderr,"ERROR: write_sdss.\n");
    return 1;
  }
}

/*.......................................................................
 *
 * Function print_offsets
 *
 * Prints out the offsets between a central position and a series of
 *  secondary positions
 *
 * Inputs: Skypos cent         central position (RA, Dec)
 *         Skypos *skypos      secondary positions (RA, Dec)
 *         Pos *offsets        (x,y) offsets
 *         int noffsets        number of offsets
 *         FILE *ofp           output file pointer (NULL if no output file)
 *
 * Output: int (0 or 1)        0 ==> success, 1 ==> error
 *
 */

int print_offsets(Skypos cent, Skypos *skypos, Pos *offsets, int noffsets,
		  FILE *ofp)
{
  int i;           /* Looping variable */
  double r,theta;  /* Temporary polar coordinate holders */
  char ew[5];      /* Direction of delta_RA */
  char ns[5];      /* Direction of delta_Dec */
  char line[MAXC];      /* General string for reading input */
  Skypos *sptr;    /* Pointer used to navigate skypos */
  Pos *optr;       /* Pointer used to navigate offsets */

  /*
   * Print out the central position
   */

  if(ofp) {
    fprintf(ofp,"#\n");
    fprintf(ofp,"#\n");
    fprintf(ofp,"#                CENTRAL SOURCE:  %s\n",cent.label);
    fprintf(ofp,"#\n");
    fprintf(ofp,"#Central source:  %02d %02d %07.4f  %+03d %02d %07.4f\n",
	    cent.hr,cent.min,cent.sec,cent.deg,cent.amin,cent.asec);
    fprintf(ofp,"#\n");
    fprintf(ofp,"#                                         ");
    fprintf(ofp,"   Shift FROM central source\n");
    fprintf(ofp,"#                                         ");
    fprintf(ofp,"   TO listed source (arcsec)\n");
    fprintf(ofp,"#Source          RA              Dec      ");
    fprintf(ofp,"    d_RA       d_Dec    d_tot  PA(N->E)\n");
    fprintf(ofp,"#----------  -------------  -------------  ");
    fprintf(ofp,"---------  ---------  ------- -------\n");
  }
  else {
    printf("#\n");
    printf("#\n");
    if(strcmp(cent.label,"") != 0) {
      printf("#                CENTRAL SOURCE:  %s\n",cent.label);
    }
    printf("#\n");
    printf("#\n");
    printf("#Central source:  %02d %02d %07.4f  %+03d %02d %07.4f\n",
	   cent.hr,cent.min,cent.sec,cent.deg,cent.amin,cent.asec);
    printf("#\n");
    printf("#                                         ");
    printf("   Shift FROM central source\n");
    printf("#                                         ");
    printf("   TO listed source (arcsec)\n");
    printf("#Source          RA              Dec      ");
    printf("    d_RA       d_Dec    d_tot  PA(N->E)\n");
    printf("#----------  -------------  -------------  ");
    printf("---------  ---------  ------- -------\n");
  }

  for(i=0,sptr=skypos,optr=offsets; i<noffsets; i++,sptr++,optr++) {

    /*
     * Determine the direction of the shift
     */

    if(optr->x >= 0.0)
      sprintf(ew,"E");
    else
      sprintf(ew,"W");
    if(optr->y >= 0.0)
      sprintf(ns,"N");
    else
      sprintf(ns,"S");

    /*
     * Convert the x,y offsets into r,theta offsets
     */

    xy2rth(*optr,&r,&theta,1);

    /*
     * Print output to file or to STDOUT
     */

    if(ofp) {
      fprintf(ofp,"%11s  %02d %02d %07.4f  %+03d %02d %06.3f ",
	      sptr->label,sptr->hr,sptr->min,sptr->sec,
	      sptr->deg,sptr->amin,sptr->asec);
      fprintf(ofp,"%8.2f %1s %8.2f %1s %8.2f  %+6.1f\n",
	      fabs(optr->x),ew,fabs(optr->y),ns,r,theta);
    }
    else {
      printf("%11s  %02d %02d %07.4f  %+03d %02d %06.3f ",
	      sptr->label,sptr->hr,sptr->min,sptr->sec,
	      sptr->deg,sptr->amin,sptr->asec);
      printf("%8.2f %1s %8.2f %1s %8.2f  %+6.1f\n",
	      fabs(optr->x),ew,fabs(optr->y),ns,r,theta);
    }
  }

  if(ofp) {
    fprintf(ofp,"#\n");
    fprintf(ofp,
	    "#-----------------------------------------------------------");
    fprintf(ofp,"-----------------\n");
  }
  else {
    printf("#\n");
    printf("#-----------------------------------------------------------");
    printf("-----------------\n");
  }
    
  return 0;
}
 
/*.......................................................................
 *
 * Function read_floats
 *
 * Reads data from an input file that contains a specified number
 *  of columns.  All data are read in as floats -- the function calling
 *  read_floats can later change the data type if needed.
 *
 * Inputs: char *inname        name of input file
 *         int ncols           number of columns to read
 *         char comment        comment character
 *         int *nlines         number of lines in input file -- set by
 *                              this function.
 *
 * Output: float *newdata[]    arrays of floats read from input file
 *
 */

/*
 * float *read_floats[](char *inname, int ncols, char comment, int *nlines)
 */
