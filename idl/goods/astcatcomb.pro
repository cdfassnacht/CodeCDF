pro astcatcomb, cat1, cat2, outfile, prefix=prefix, fitsdir=fitsdir, $
   catdir=catdir,outdir=outdir,stsize=stsize,suffix=suffix,_extra=_extra
;
; Procedure astcatcomb
;
; Description: Combines two astrometric catalogs into a format appropriate
;               for an input file for distcalc.c in its -p option.
;
; Inputs: cat1        (string)     name of first catalog file.  The
;                                   format is: id hr min sec deg amin asec
;         cat2        (string)     name of second catalog file.  The
;                                   format is the same as for cat1.
;         outfile     (string)     output file name
;         [fitsdir=]  (string)     directory of input fits file --
;                                   default is '.'
;         [catdir=]   (string)     directory of input catalog file --
;                                   default is '.'
;         [outdir=]   (string)     directory for output file --
;                                   default is '.'
;
; Revision history:
;  2003Apr29 Chris Fassnacht -- First working version.
;

; Check input format

if n_params() lt 3 then begin
    print, ''
    print, $
      'syntax: astcatcomb, cat1, cat2, outfile'
    print, ''
    print, ''
    return
endif

; If optional input variables have not been set, use defaults

if not keyword_set(catdir) then catdir = '.'
if not keyword_set(outdir) then outdir = '.'

; Read in catalogs with readcol

incat1=catdir+'/'+cat1
print, 'astcatcomb: Reading catalog from ',incat1
readcol, incat1, id1,hr1,min1,sec1,deg1,amin1,asec1, $
   format='a,i,i,d,i,i,d', comment='#'
n1 = n_elements(id1)

incat2=catdir+'/'+cat2
print, 'astcatcomb: Reading catalog from ',incat2
readcol, incat2, id2,hr2,min2,sec2,deg2,amin2,asec2, $
   format='a,i,i,d,i,i,d', comment='#'
n2 = n_elements(id2)

; Check that catalogs have the same size

if n1 ne n2 then begin
   print, ''
   print, '*** ERROR: Catalogs have different sizes. ***'
   print, '   Catalog 1: '+incat1+' has n_lines = ',n1
   print, '   Catalog 2: '+incat2+' has n_lines = ',n2
   print, ''
   return
endif

; Print output to screen and file

fmt='(A,2(I4,I3,F7.3,I4,I3,F7.3))'

forprint, id1,hr1,min1,sec1,deg1,amin1,asec1,hr2,min2,sec2,deg2,amin2,asec2, $
   format=fmt
forprint, id1,hr1,min1,sec1,deg1,amin1,asec1,hr2,min2,sec2,deg2,amin2,asec2, $
   format=fmt,text=outfile,/nocomment

end
