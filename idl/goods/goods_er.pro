pro goods_er, infile, outinfo, band=band
;
; Procedure goods_er
;
; Description: Provides interactive measurements of distances from selected
;               objects in an image to the central object.
;
; Inputs: infile     (string)      base name of input FITS file (root+id number)
;         outinfo    (structure)   output structure containing info on the
;                                   image
;         [band=]    (string)      band (b,v,i,z) for input image -- default
;                                   is 'v'
;
;
; Revision history:
;  2003Apr28 Chris Fassnacht -- First working version.
;  2003Jul02 Chris Fassnacht -- Moved catalog reading to read_secat.pro
;                               Moved central source matching to
;                                find_centsrc.pro
;

; Check input format

if n_params() lt 2 then begin
    print, ''
    print, 'syntax: goods_er, infile, outinfo, [,band=band]'
    print, ''
    return
endif

; Initialize some variables and set optional variables to defaults if
;  not set by the function call

if not keyword_set(band) then band='v'

; Set directory names

datadir = band+'-band/'
catdir = 'Stampcat'
residdir = 'Resid/'
fitdir = 'Fit/'
maskdir = 'Mask/'

; Read in catalog with read_secat

catname = infile+'.cat'
read_secat, catname, secat, catdir=catdir, format=1
ncat = n_elements(secat)


; Read postage-stamp fits file into a float array called stamp
; ** NB: For now assume no errors

datafile = datadir+infile+'_'+band+'.fits'
print, ''
print, 'goods_er: Reading image from ',datafile
stamp = mrdfits(datafile, 0, inhead)

; Display input image

dispim, stamp

; Go through input catalog and find object closest to the image center.

find_centsrc, stamp, secat, best

; Display best object parameters

tvellipse, best.a, best.b, best.x, best.y, best.pa, /data, color=255

; Get image positions by clicking on them

print, ''
mousebutt = 0
rbi = 0
j = 0
while mousebutt ne 2 do begin
   print, ''
   print, 'Hit left mouse button to compute and save distance to central object'
   print, 'Hit right mouse button on two objects to compute their distances'
   print, 'Hit middle mouse button to quit'
   print, ''
   cursor,xx,yy,/down
   mousebutt = !mouse.button
   print, xx, yy
   if mousebutt eq 1 then begin
      dx = xx - best.x
      dy = yy - best.y
      offset = sqrt(dx*dx + dy*dy)
      print, 'Offset = ',offset,' pixels or ',offset*0.05,' arcsec'
      if j eq 0 then erarr = offset else erarr = [erarr, offset]
   endif
   if mousebutt eq 4 then begin
      if rbi eq 0 then begin
         x0 = xx
         y0 = yy
         print,''
         print,'*** Right mouse ==> first position marked. ***'
         print,'*** Hit right mouse again to compute distance between objects'
      endif else begin
         dx = xx - x0
         dy = yy - y0
         offset = sqrt(dx*dx + dy*dy)
         print, 'Offset = ',offset,' pixels or ',offset*0.05,' arcsec'
         if rbi eq 1 then offarr = offset else offarr = [offarr, offset]
      endelse
      rbi = rbi + 1
   endif
   j = j + 1
endwhile

; Fix empty arrays if needed

if n_elements(erarr) eq 0 then erarr = 0.0
if n_elements(offarr) eq 0 then offarr = 0.0

; Transfer info to output structure

outinfo = create_struct('object',infile,$
                        'band',band,$
                        'xlens',best.x,$
                        'ylens',best.y,$
                        'r_e',erarr,$
                        'r_2pt',offarr)

end
