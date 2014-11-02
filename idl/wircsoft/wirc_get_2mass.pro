; automatically retrieve a 2mass dataset for the
; wirc photometric calibration routines
;
; reads in the headers of the input fits filelist
; and figures out an appropriate radius and center
;
; NOTE - the input file list must have at least approximate astrometry (i.e.
;        the "ff_*" files!
;
pro wirc_get_2mass,filelist,outfile,$
                   VERBOSE=VERBOSE


; read in the input file list
readcol,filelist,infile,format="A"


nfiles=n_elements(infile)
ra=fltarr(nfiles)
dec=fltarr(nfiles)

; first, read in all the ra and dec values
for i=0,nfiles-1 do begin

   hd=headfits(infile[i])
   ra[i]=sxpar(hd,'CRVAL1')
   dec[i]=sxpar(hd,'CRVAL2')

endfor

; get the pixel scale
; cdelt will come back in degrees
extast,hd,hd_struct
getrot,hd_struct,crota,cdelt
naxis1=sxpar(hd,'NAXIS1')
naxis2=sxpar(hd,'NAXIS2')

; find the center of the image set
ra_center=mean(ra)
dec_center=mean(dec)


; now find the width in x and y
; we assume crpix is in the center of the array
;
ra_width=max(ra)-min(ra)+(cdelt[0]*naxis1)
dec_width=max(dec)-min(dec)+(cdelt[1]*naxis2)


; compute radius in arcminutes
; the root-2 is eeded for a square geometry, plus 5%
radius=1.05*1.414*60.*(max([ra_width,dec_width])/2)

; print some informational messages
if keyword_set(VERBOSE) then begin
   print,'Using RA,DEC center of ',ra_center,dec_center,' degrees.'
   print,'Radius is ',radius,' arcminutes.'
endif



; get the 2mass catalog
tmass=queryvizier('2MASS-PSC',[ra_center,dec_center],radius,/CANADA)


; the catalog itself is a structure
;RAJ2000 DEJ2000 _2MASS JMAG E_JMAG HMAG E_HMAG KMAG E_KMAG 
;  QFLG RFLG BFLG CFLG XFLG AFLG

nstars=n_elements(tmass._2MASS)

; write the basic information out

get_lun,unit
openw,unit,outfile

for i=0,nstars-1 do begin

printf,unit,tmass[i].RAJ2000,tmass[i].DEJ2000,tmass[i].JMAG,$
       tmass[i].HMAG,tmass[i].KMAG


endfor


free_lun,unit

end