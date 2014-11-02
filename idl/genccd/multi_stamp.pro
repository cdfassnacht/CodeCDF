pro multi_stamp, fitsfile, catfile, fitsdir=fitsdir,catdir=catdir, $
   outroot=outroot,outdir=outdir,stsize=stsize,suffix=suffix
;
; Procedure multi_stamp
;
; Description: A driver to read in an image from a fits file of a
;  field and to read in a catalog of objects detected in that field.
;  Then call the poststamp procedure for each object in the catalog.
;
; Inputs: fitsfile    (string)     name of input FITS file
;         catfile     (string)     name of catalog file.  The file is
;                                   assumed to have (at least) three
;                                   columns in the ff. order:
;                                    [0]: id of object
;                                    [1]: x position of object
;                                    [2]: y position of object
;         [fitsdir=]  (string)     directory of input fits file --
;                                   default is '.'
;         [catdir=]   (string)     directory of input catalog file --
;                                   default is '.'
;         [outroot=]  (string)     root name for output file --
;                                   default is 'cutout'
;         [outdir=]   (string)     directory for output file --
;                                   default is '.'
;         [stsize=]   (int)        size for postage stamp image, in
;                                   pixels -- default = 201
;         [suffix=]   (string)     suffix for output file --
;                                   default is no suffix
;
; Revision history:
;  2003Mar31 Chris Fassnacht -- First working version.  Loosely based on the
;                                old dosearch.pro.
;  2003Jul08 Chris Fassnacht -- Fixed to use the new form of poststamp,
;                                 which is much more efficient in terms of
;                                 memory.
;

; Check input format

if n_params() lt 2 then begin
    print, ''
    print, $
      'syntax: multi_stamp, fitsfile, catfile [,fitsdir=fitsdir,catdir=catdir,'
    print, '       outroot=outroot,outdir=outdir,stsize=stsize,suffix=suffix]'
    print, ''
    print, '  Default values:'
    print, "    fitsdir = '.'"
    print, "    catdir = '.'"
    print, "    outroot = 'cutout'"
    print, '    stsize  = 201'
    print, ''
    return
endif

; If optional input variables have not been set, use defaults

if not keyword_set(fitsdir) then fitsdir = '.'
if not keyword_set(catdir) then catdir = '.'
if not keyword_set(outdir) then outdir = '.'
if not keyword_set(outroot) then outroot = 'cutout'
if not keyword_set(stsize) then stsize = 201

; Read in catalog with readcol

incat=catdir+'/'+catfile
print, 'multi_stamp: Reading catalog from ',incat
readcol, catfile, id, x, y, format='A,F,F', comment='#'

; Read in fits file with mrdfits
; ** NB: For now assume no errors

infits=fitsdir+'/'+fitsfile

; Loop on input list, calling poststamp for each object

print, ''
if n_elements(id) gt 0 then begin
   print, 'multi_stamp: Starting image loop.'
   print, '---------------------------------------------------------------'
   for i=0,n_elements(id)-1 do begin
      basestr = id[i]
      case strlen(basestr) of
         1: idstr = '0000'+basestr
         2: idstr = '000'+basestr
         3: idstr = '00'+basestr
         4: idstr = '0'+basestr
         else: idstr = basestr
       endcase
       if not keyword_set(suffix) then outname = outroot+'_'+idstr $
          else outname = outroot+'_'+idstr+'_'+suffix
       print, ''
       print, ' *** STARTING OBJECT #',strtrim(i+1,1),' (',outname,'). ***'
       print, ''
       outim=poststamp(infits,x[i],y[i],stsize,stsize,outhead=outhead)
       outfits = outdir+'/'+outname+'.fits'
       writefits, outfits, outim, outhead
       print, ''
       print, ' *** FINISHED WITH OBJECT #',strtrim(i+1,1),' (',outname,'). ***'
       print, ''
       print, '---------------------------------------------------------------'
   endfor
endif else begin
   print, ''
   print, '*** ERROR: multi_stamp.  No objects found in catalog list. ***'
   print, ''
endelse

end
