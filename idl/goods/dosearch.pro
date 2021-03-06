pro dosearch, catfile, format=format, root=root, datadir=datadir, band=band, $
              mskfac=mskfac, usegauss=usegauss, dofitell=dofitell, $
              dodisplay=dodisplay
;
; Procedure dosearch
;
; Description: A driver to read in a catalog of image properties,
;  generated from SExtractor.  Based on the object IDs and image
;  properties, this procedure creates the proper inputs and then
;  calls the medsub procedure (formerly known as gopost) for each object 
;  in the catalog.
;  ***NB: The postage stamp cutouts of the objects must already have
;    been created by some other process.
;
; Inputs: catfile     (string)     name of catalog file.  The file is
;                                   assumed to have (at least) six
;                                   columns in the ff. order:
;                                    [0]: id number
;                                    [1]: x position of object
;                                    [2]: y position of object
;                                    [3]: PA of semimajor axis 
;                                      (astronomical convention)
;                                    [4]: ellipticity (1 - b/a)
;                                    [5]: "size" - for now, SExtractor's FWHM
;         [format=]   (int)        format of input catalog file to use:
;                                   1 ==> use individual SExtractor catalogs
;                                         produced by stampinfo.sh
;                                         (in this case catfile contains only
;                                          id numbers)
;                                   2 ==> use modified GOODS pipeline
;                                         SExtractor catalogs, as used for
;                                         the initial v1.7 lens search
;                                   Default = 1
;         [root=]     (string)     root name for input files --
;                                   default is 'goods'
;         [band=]     (string)     observing band (default is z)
;         [mskfac=]   (float)      factor by which galaxy "size" is 
;                                   multiplied to define region of
;                                   "good" pixels for the mask in the
;                                   medsub procedure.
;         [/dofitell]              use fitell procedures to determine
;                                   galaxy center, ellipticity, and
;                                   PA.  Otherwise, just use
;                                   SExtractor values, which this
;                                   procedure obtains from the input catfile.
;         [/dodisplay]             display images if set
;
; Revision history:
;  2003Feb25 Chris Fassnacht -- A simple modification of the old dosearch.pro
;  2003Feb28 Chris Fassnacht -- Added the optional "band" passed parameter.
;  2003Apr02 Chris Fassnacht -- Corrected for math orientation here rather
;                                than in medsub.pro.
;  2003Apr16 Chris Fassnacht -- Added option to use gauss2dfit PA by setting
;                                usegauss flag
;  2003Jun22 Chris Fassnacht -- Modified to use the new secat structure.
;

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, 'This procedure acts as an interface to call medsub'
    print, ''
    print, $
      'syntax: dosearch, catfile [, format=format,root=root,datadir=datadir, '
    print, '            band=band,mskfac=mskfac,/usegauss,/dofitell,/dodisplay]'
    print, ''
    print, '  Default values:'
    print, '    format = 1'
    print, "    root = 'goods'"
    print, "    datadir = 'z-band'"
    print, "    band = 'z'"
    print, '    mskfac  = 1.0'
    print, ''
    return
endif

; If optional input variables have not been set, use defaults

if n_elements(format) eq 0 then format = 1
if not keyword_set(root) then root = 'goods'
if not keyword_set(datadir) then datadir = 'z-band'
if not keyword_set(band) then band = 'z'
if not keyword_set(mskfac) then mskfac = 1.0
if not keyword_set(dofitell) then dofitell = 0
if not keyword_set(dodisplay) then dodisplay = 0
if not keyword_set(usegauss) then usegauss = 0

; Read in catalog

if format eq 2 then begin
  read_lscat, catfile, secat
endif else begin
  read_secat, catfile, secat
endelse
ncat = n_elements(secat)

; Loop on input list, calling medsub for each object

print, ''
if ncat gt 0 then begin
    print, 'dosearch: Starting image loop.'
    print, '---------------------------------------------------------------'
    for i=1,ncat do begin
        basestr = strtrim(secat[i-1].id,1)
        case strlen(basestr) of
            1: idstr = '0000'+basestr
            2: idstr = '000'+basestr
            3: idstr = '00'+basestr
            4: idstr = '0'+basestr
            else: idstr = basestr
        endcase
        filename = root+'_'+idstr+'_'+band+'.fits'
        outbase = root+'_'+idstr+'_'+band
        print, ''
        print, ' *** STARTING OBJECT #',strtrim(i,1),' (',filename,'). ***'
        print, ''
        gshape = [secat[i-1].fwhm, secat[i-1].ellip, secat[i-1].pa]
        medsub, filename, outbase, secat[i-1], datad=datadir, maskfac=mskfac, $
           usegauss=usegauss, dofitell=dofitell, dodisplay=dodisplay
        print, ''
        print, ' *** FINISHED WITH OBJECT #',strtrim(i,1),' (',filename,'). ***'
        print, ''
        print, '---------------------------------------------------------------'
    endfor
endif else begin
    print, ''
    print, '*** ERROR: dosearch.  No objects found in catalog list. ***'
    print, ''
endelse



end
