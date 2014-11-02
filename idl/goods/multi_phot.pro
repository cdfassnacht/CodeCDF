pro multi_phot, catfile, prefix=prefix, catdir=catdir, $
   outfile=outfile,outdir=outdir,stsize=stsize,suffix=suffix,_extra=_extra
;
; Procedure multi_phot
;
; Description: Calls er_phot for each object in a list, which is passed
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
      'syntax: multi_phot, catfile [,catdir=catdir,'
    print, '       outfile=outfile,outdir=outdir]'
    print, ''
    print, '  Default values:'
    print, "    catdir = '.'"
    print, "    outfile = 'phot.dat'"
    print, ''
    return
endif

; If optional input variables have not been set, use defaults

if not keyword_set(catdir) then catdir = '.'
if not keyword_set(outdir) then outdir = '.'
if not keyword_set(outfile) then outfile = 'phot.dat'

; Read in catalog with readcol

incat=catdir+'/'+catfile
print, 'multi_phot: Reading catalog from ',incat
readcol, catfile, id,xlens,ylens,er, format='A,D,D,D', comment='#'
nid = n_elements(id)

; Create holders for b_e, v_e, i_e and z_e

b_e = fltarr(nid)
v_e = fltarr(nid)
i_e = fltarr(nid)
z_e = fltarr(nid)

; Loop on input list, calling goods_er for each object

print, ''
if n_elements(id) gt 0 then begin
   print, 'multi_phot: Starting loop.'
   print, '---------------------------------------------------------------'
   for i=0,nid-1 do begin
      if keyword_set(prefix) then basestr = prefix+'_'+id[i] else $
         basestr = id[i]
      er_phot, basestr, xlens[i],ylens[i],er[i],outinfo
      b_e[i] = outinfo.b_e[0]
      v_e[i] = outinfo.v_e[0]
      i_e[i] = outinfo.i_e[0]
      z_e[i] = outinfo.z_e[0]
   endfor
endif else begin
   print, ''
   print, '*** ERROR: multi_phot.  No objects found in catalog list. ***'
   print, ''
endelse

print,''
print,''
print,'multi_phot: Final photometry'
print,'----------------------------'
print,'     Object      ER(pix) ER(asec)  B_E   V_E   i_E   z_E'
er_arcsec = er * 0.05
fmt='(A15,F7.2,F8.2,F9.2,F6.2,F6.2,F6.2)'
forprint, id, er, er_arcsec, b_e, v_e, i_e, z_e, format=fmt
forprint, id, er, er_arcsec, b_e, v_e, i_e, z_e, text=outfile, format=fmt

end
