pro multi_stamp_color, bfits, gfits, rfits, catfile, $
                       bwtfile=bwtfile, gwtfile=gwtfile, rwtfile=rwtfile, $
                       fitsdir=fitsdir,catdir=catdir,stsize=stsize, $
                       lsig=lsig,hsig=hsig,jpegroot=jpegroot
;
; Procedure multi_stamp_color
;
; Description: A driver to read in an image from 3 fits files of a
;  field and to read in a catalog of objects detected in that field.
;  Then call the poststamp_color procedure for each object in the catalog.
;
; Inputs: 
;   bfits       (string)     name of shortest wavelength fits file
;   gfits       (string)     name of middle wavelength fits file
;   rfits       (string)     name of longest wavelength fits file
;   catfile     (string)     name of catalog file.  The file is
;                             assumed to have (at least) three
;                             columns in the ff. order:
;                              [0]: id of object
;                              [1]: x position of object
;                              [2]: y position of object
;   [bwtfile=] (string)      optional weight file for bfits
;   [gwtfile=] (string)      optional weight file for gfits
;   [rwtfile=] (string)      optional weight file for rfits
;   [fitsdir=]  (string)     directory of input fits file --
;                             default is '.'
;   [catdir=]   (string)     directory of input catalog file --
;                             default is '.'
;   [stsize=]   (int)        size for postage stamp image, in
;                             pixels -- default = 201
;   [lsig=]     (float)      number of sigma below mean to set lower
;                             display limit (default = 1.0)
;   [hsig=]     (float)      number of sigma above mean to set upper
;                             display limit (default = 30.0)
;   [jpegroot=] (string)     root name for jpeg output files.
;                             Procedure will not create jpegs unless
;                             this parameter is set.
;
; Revision history:
;  2006Feb02 Chris Fassnacht -- First working version.
;  2007Aug13 Chris Fassnacht -- Added the possibility of including
;                               weight images for the input fits
;                               files.
;                               Added jpeg capability through new
;                               jpegroot passed parameter
;

; Check input format

if n_params() lt 4 then begin
    print, ''
    print, $
      'syntax: multi_stamp_color, bfits, gfits, rfits, catfile '
    print, '       [bwtwfits=bwtfits,gwtfits=gwtfits,rwtfits=rwtfits,'
    print, '       ,fitsdir=fitsdir,catdir=catdir,lsig=lsig,hsig=hsig,'
    print, '       stsize=stsize]'
    print, ''
    print, 'Optional parameters:'
    print, '  stsize:        Size of image axis in pixels. Default=201'
    print, '  bwtfits, etc.: Weight files associated with the input fits files'
    print, '  fitsdir:       Directory containing the fits files.  Default is'
    print, '                  the current working directory'
    print, '  catdir:        Directory containing the catalog file.  Default is'
    print, '                  the current working directory'
    print, '  lsig:          Lower limit of display range, in units of'
    print, '                  sigma below the clipped mean.  Default=1.0'
    print, '  hsig:          Upper limit of display range, in units of'
    print, '                  sigma above the clipped mean.  Default=30.0'
    print, '  jpegroot:      Root name for optional output jpeg files.  There'
    print, '                  is no default.  If this keyword is not set,'
    print, '                  then no output jpeg files will be created.'
    print, ''
    return
endif

; Set up control variable

convar = 0
if n_elements(jpegroot) gt 0 then convar = convar + 1

; If optional input variables have not been set, use defaults

if not keyword_set(fitsdir) then fitsdir = '.'
if not keyword_set(catdir) then catdir = '.'
if not keyword_set(stsize) then stsize = 201
if not keyword_set(lsig) then lsig=1.0
if not keyword_set(hsig) then hsig=30.0

; Read in catalog with readcol

incat=catdir+'/'+catfile
print, 'multi_stamp_color: Reading catalog from ',incat
readcol, catfile, id, x, y, format='A,F,F', comment='#'
neid = n_elements(id)
if(neid eq 0) then begin
   print, ''
   print, '*** ERROR: multi_stamp_color.  No objects found in catalog list. ***'
   print, ''
   return
endif

; Set up input fits files

inb=fitsdir+'/'+bfits
ing=fitsdir+'/'+gfits
inr=fitsdir+'/'+rfits
if n_elements(bwtfile) gt 0 then begin 
    inbwt=fitsdir+'/'+bwtfile
    convar = convar + 2
endif
if n_elements(gwtfile) gt 0 then begin
    ingwt=fitsdir+'/'+gwtfile
    convar = convar + 4
endif
if n_elements(rwtfile) gt 0 then begin
    inrwt=fitsdir+'/'+rwtfile
    convar = convar + 8
endif

; Loop on input list, calling poststamp_color for each object

; Start the loop

userin = ''
i=0
mousebutt = 0
print, ''
print, 'multi_stamp_color: Starting image loop.'
print, '---------------------------------------------------------------'
for i=0,neid-1 do begin
    print, ''
    print, ' *** STARTING OBJECT #',strtrim(i+1,1),'  ***'
    print, ''
    ; Set up name for jpeg file
    if n_elements(jpegroot) gt 0 then $
      jpfile=jpegroot+'_'+strtrim(x[i],2)+'_'+strtrim(y[i],2)+'.jpg'
    ; Do call to poststamp_color depending on value of convar
    ; Note that the case statements get called in order, so put most likely
    ;  values of convar earliest
    case convar of
        0: poststamp_color,inb,ing,inr,x[i],y[i],stsize,stsize,lsig=lsig,$
          hsig=hsig
        1: poststamp_color,inb,ing,inr,x[i],y[i],stsize,stsize,lsig=lsig,$
          hsig=hsig,jpfile=jpfile
        14: poststamp_color,inb,ing,inr,x[i],y[i],stsize,stsize,lsig=lsig,$
          hsig=hsig,bwtfile=inbwt,gwtfile=ingwt,rwtfile=inrwt
        15: poststamp_color,inb,ing,inr,x[i],y[i],stsize,stsize,lsig=lsig,$
          hsig=hsig,jpfile=jpfile,bwtfile=inbwt,gwtfile=ingwt,rwtfile=inrwt
        2: poststamp_color,inb,ing,inr,x[i],y[i],stsize,stsize,lsig=lsig,$
          hsig=hsig,bwtfile=inbwt
        3: poststamp_color,inb,ing,inr,x[i],y[i],stsize,stsize,lsig=lsig,$
          hsig=hsig,jpfile=jpfile,bwtfile=inbwt
        4: poststamp_color,inb,ing,inr,x[i],y[i],stsize,stsize,lsig=lsig,$
          hsig=hsig,gwtfile=ingwt
        5: poststamp_color,inb,ing,inr,x[i],y[i],stsize,stsize,lsig=lsig,$
          hsig=hsig,jpfile=jpfile,gwtfile=ingwt
        6: poststamp_color,inb,ing,inr,x[i],y[i],stsize,stsize,lsig=lsig,$
          hsig=hsig,bwtfile=inbwt,gwtfile=ingwt
        7: poststamp_color,inb,ing,inr,x[i],y[i],stsize,stsize,lsig=lsig,$
          hsig=hsig,jpfile=jpfile,bwtfile=inbwt,gwtfile=ingwt
        8: poststamp_color,inb,ing,inr,x[i],y[i],stsize,stsize,lsig=lsig,$
          hsig=hsig,rwtfile=inrwt
        9: poststamp_color,inb,ing,inr,x[i],y[i],stsize,stsize,lsig=lsig,$
          hsig=hsig,jpfile=jpfile,rwtfile=inrwt
        10: poststamp_color,inb,ing,inr,x[i],y[i],stsize,stsize,lsig=lsig,$
          hsig=hsig,bwtfile=inbwt,rwtfile=inrwt
        11: poststamp_color,inb,ing,inr,x[i],y[i],stsize,stsize,lsig=lsig,$
          hsig=hsig,jpfile=jpfile,bwtfile=inbwt,rwtfile=inrwt
        12: poststamp_color,inb,ing,inr,x[i],y[i],stsize,stsize,lsig=lsig,$
          hsig=hsig,gwtfile=ingwt,rwtfile=inrwt
        13: poststamp_color,inb,ing,inr,x[i],y[i],stsize,stsize,lsig=lsig,$
          hsig=hsig,jpfile=jpfile,gwtfile=ingwt,rwtfile=inrwt
    endcase
    print, ''
    print, ' *** FINISHED WITH OBJECT #',strtrim(i+1,1),' (ID = ',id[i],')  ***'
    print, ''
    read, userin, prompt='Enter return to go to next object: '
    print, ''
    print, '---------------------------------------------------------------'
endfor

end
