pro goods_reval, catfile, prefix=prefix, rap=rap, catdir=catdir, $
   outfile=outfile,outdir=outdir,stsize=stsize,suffix=suffix,_extra=_extra
;
; Procedure goods_reval
;
; Description: Calls eval_resid for each object in a list, which is passed
;               to the procedure as catfile.
;
; Inputs: fitsfile    (string)     name of input FITS file
;         catfile     (string)     name of catalog file.  The file is
;                                   assumed to have (at least) three
;                                   columns in the ff. order:
;                                    [0]: id of object
;                                    [1]: x position of object
;                                    [2]: y position of object
;         [indir=]    (string)     directory of input fits files --
;                                   default is 'z-band'
;         [rsddir=]   (string)     directory of input residualt files --
;                                   default is 'Resid'
;         [catdir=]   (string)     directory of input catalog file --
;                                   default is '.'
;         [outfile=]  (string)     output file name --
;                                   default is 'rsdval.dat'
;         [outdir=]   (string)     directory for output file --
;                                   default is '.'
;
; Revision history:
;  2003Jun14 Chris Fassnacht -- First working version.  Based loosely
;                                on multi_phot.pro
;

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, $
      'syntax: goods_reval, catfile [,prefix=prefix,rap=rap,catdir=catdir,'
    print, '       outfile=outfile,outdir=outdir]'
    print, ''
    print, '  Default values:'
    print, "    catdir = '.'"
    print, '    rap = 30'
    print, "    outdir = '.'"
    print, "    outfile = 'rsdval.dat'"
    print, ''
    return
endif

; If optional input variables have not been set, use defaults

if not keyword_set(catdir) then catdir = '.'
if not keyword_set(outdir) then outdir = '.'
if not keyword_set(outfile) then outfile = 'rsdval.dat'
if not keyword_set(rap) then rap=30

; Read in catalog with readcol

incat=catdir+'/'+catfile
print, 'goods_reval: Reading catalog from ',incat
readcol, catfile, id,xlens,ylens, format='A,D,D', comment='#'
nid = n_elements(id)

; Create holder for rsdrat

rsdrat = fltarr(nid)

; Loop on input list, calling eval_resid for each object

print, ''
if n_elements(id) gt 0 then begin
   print, 'goods_reval: Starting loop.'
   print, '---------------------------------------------------------------'
   xcent = 101
   ycent = 101
   rbkgd = 60
   for i=0,nid-1 do begin
      if keyword_set(prefix) then begin
         inname = prefix+'_'+id[i]+'_z.fits' 
         rsdname = prefix+'_'+id[i]+'_z_resid.fits' 
      endif else begin
         inname = id[i]+'_z.fits'
         rsdname = id[i]+'_z_resid.fits'
      endelse
      eval_resid, inname,rsdname,xcent,ycent,rap,rbkgd,outstr
      rsdrat[i] = outstr.rsdrat
      print, ''
      print, '---------------------------------------------------------------'
      print, ''
   endfor
endif else begin
   print, ''
   print, '*** ERROR: goods_reval.  No objects found in catalog list. ***'
   print, ''
endelse

print,''
print,''
print,'goods_reval: Final info'
print,'----------------------------'
fmt='(A15,F10.3)'
;forprint, id, rsdrat, format=fmt
forprint, id, rsdrat, text=outfile, format=fmt

end
