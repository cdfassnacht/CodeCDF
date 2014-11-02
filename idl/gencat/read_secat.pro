pro read_secat, catfile, outinfo, catdir=catdir, format=format
;
; Procedure read_secat
;
; Description: Reads in data from a catalog (produced by SExtractor) and
;  returns the information in a structure array.
;
; Inputs: catfile    (string)      base name of input FITS file (root+id number)
;         outinfo    (structarr)   output structure array containing catalog
;                                   info
;         [catdir=]  (string)      directory containing the catalog -- default
;                                   is 'Stampcat'
;         [format=]  (int)         Format of input catalog:
;                                   0 = id number only
;                                   1 = stampinfo.sh 
;                                   2 = GOODS pipeline
;                                   (default value is 1)
;
;
; Revision history:
;  2003Jul01 Chris Fassnacht -- First working version.
;  2003Jul02 Chris Fassnacht -- Added id-only format
;

; Check input format

if n_params() lt 2 then begin
    print, ''
    print, 'syntax: read_secat, catfile, outinfo'
    print, '         [,catdir=catdir,format=format]'
    print, '   Format values:'
    print, '    0 ==> id number only'
    print, '    1 ==> stampinfo.sh'
    print, '    2 ==> GOODS pipeline'
    print, ''
    print, 'Default values:'
    print, "  catdir:  'Stampcat'"
    print, '  format:   1'
    print, ''
    return
endif

; Initialize some variables and set optional variables to defaults if
;  not set by the function call

if not keyword_set(catdir) then catdir='Stampcat'
if n_elements(format) eq 0 then format=1

; Set input filename

infile = catdir+'/'+catfile

; Read in catalog with readcol, depending on format
; *** Right now, only use default format 1 *** FIX LATER ***

case format of
   0:    begin
      print, 'read_secat: Reading catalog from ',infile
      readcol, infile, id,comment='#',format=''
    end
   1:    begin
      print, 'read_secat: Reading catalog from ',infile
      readcol, infile, x,y,flg,clss,miso,eiso,jnk,jnk,rk,bkg,thrsh,muthr,aiso, $
                  a,b,pa,comment='#',format=''
    end
   else: begin
      print, ''
      print, 'ERROR: read_secat.  Not a valid format'
      print, ''
      return
    end
endcase

print, ''

; Get number of elements in arrays

ncat = 0
ncat=n_elements(x)
if ncat eq 0 then begin
   print, ''
   print, '*** ERROR: read_secat.  No objects found in catalog list. ***'
   print, ''
endif

; Set up structure array

def_secat, template
outinfo = replicate(template,ncat)

; Read info into structure array

for i=0,ncat-1 do begin
   outinfo[i].x = x[i]
   outinfo[i].y = y[i]
   outinfo[i].class = clss[i]
   outinfo[i].flag = flg[i]
   outinfo[i].m_iso = miso[i]
   outinfo[i].merr_iso = eiso[i]
   outinfo[i].a = a[i]
   outinfo[i].b = b[i]
   outinfo[i].ellip = (1 - b[i]/a[i])
   outinfo[i].pa = pa[i]
endfor

end
