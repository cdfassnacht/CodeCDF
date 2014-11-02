pro multi_er, catfile, prefix=prefix, fitsdir=fitsdir,catdir=catdir, $
   outfile=outfile,outdir=outdir,stsize=stsize,suffix=suffix,_extra=_extra
;
; Procedure multi_er
;
; Description: Calls goods_er for each object in a list, which is passed
;               to the procedure as catfile.
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
;         [outfile=]  (string)     output file name --
;                                   default is 'er.dat'
;         [outdir=]   (string)     directory for output file --
;                                   default is '.'
;
; Revision history:
;  2003Apr29 Chris Fassnacht -- First working version.
;

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, $
      'syntax: multi_er, catfile [,fitsdir=fitsdir,catdir=catdir,'
    print, '       outfile=outfile,outdir=outdir]'
    print, ''
    print, '  Default values:'
    print, "    fitsdir = '.'"
    print, "    catdir = '.'"
    print, "    outfile = 'er.dat'"
    print, ''
    return
endif

; If optional input variables have not been set, use defaults

if not keyword_set(fitsdir) then fitsdir = '.'
if not keyword_set(catdir) then catdir = '.'
if not keyword_set(outdir) then outdir = '.'
if not keyword_set(outfile) then outfile = 'er.dat'

; Read in catalog with readcol

incat=catdir+'/'+catfile
print, 'multi_er: Reading catalog from ',incat
readcol, incat, id, format='A', comment='#'
nid = n_elements(id)

; Create holders for xlens, ylens, and r_e

obj = strarr(nid)
xlens = fltarr(nid)
ylens = fltarr(nid)
r_e = dblarr(nid)

; Loop on input list, calling goods_er for each object

print, ''
if n_elements(id) gt 0 then begin
   print, 'multi_er: Starting loop.'
   print, '---------------------------------------------------------------'
   for i=0,nid-1 do begin
      if keyword_set(prefix) then basestr = prefix+'_'+id[i] else $
         basestr = id[i]
      goods_er, basestr, outinfo
      xlens[i] = outinfo.xlens
      ylens[i] = outinfo.ylens
      r_e[i] = outinfo.r_e
   endfor
endif else begin
   print, ''
   print, '*** ERROR: multi_er.  No objects found in catalog list. ***'
   print, ''
endelse

forprint, id, xlens, ylens, r_e
forprint, id, xlens, ylens, r_e, text=outfile

end
