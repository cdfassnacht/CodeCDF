; used for checking the image final calibration
;
pro wirc_finalcalib,infile,$
                     TMASSFILE=TMASSFILE,$ optional name of input 2mass table
                     MATCHRADIUS=MATCHRADIUS,$ arcseconds
                     FILTER=FILTER,$ filter override, JHK=123
                     PLOTROOT=PLOTROOT,$ ; rootname for storing plotfiles
                     VERBOSE=VERBOSE ; /verbose, print messages
                     

; set up the display window
 window,xsize=640,ysize=480
 
 
if keyword_set(MATCHRADIUS) NE 1 then MATCHRADIUS=4. 
 

; first, read the input image
if keyword_set(VERBOSE) then print,'Reading ',infile
image=readfits(infile,hd,/silent)
naxis1=sxpar(hd,'NAXIS1')
naxis2=sxpar(hd,'NAXIS2')
magzpt=sxpar(hd,'MAGZPT')
filter2=sxpar(hd,'AFT')
crval1=sxpar(hd,'crval1')
crval2=sxpar(hd,'crval2')
getrot,hd,roty,cdelt

; condition the image for sextractor
;
; set nans > saturated
sex_image=image
badpix=where(finite(sex_image,/infinity),count)
if count GT 0 then sex_image(badpix)=999999
badpix=where(finite(sex_image,/nan),count)
if count GT 0 then sex_image(badpix)=999999


; write temporary image for sextractor, and run
;  to get the comparison source list
if keyword_set(VERBOSE) then print,'Running sextractor.'
writefits,'sex_temp.fits',sex_image,hd
spawn,'/usr/local/bin/sex sex_temp.fits -c sex_files/finalcalib.sex > sex.log'

; delete the temporary file
spawn,'rm sex_temp.fits'


; this checks for validity of the sextractor solution
sex_detect=file_test('finalcalib.txt',/zero_length)

; read the file if non-zero length
if (sex_detect NE 1) then readcol,'finalcalib.txt',cntr,x,y,mag_iso,mag_auto,mag_best,isoarea,ra,dec,theta,flags,fwhm,elong,ellip,star,/silent
print,'Number of sextractor detected objects: '+string(n_elements(cntr))


if keyword_set(TMASSFILE) then begin

         ; read the input 2mass table
         readcol,TMASSFILE,mra,mdec,j,h,k

endif else begin

print,'Fetching 2mass data directly.'

; get the 2mass data
;   first get the center RA/DEC
xyad,hd,round(naxis1/2),round(naxis2/2),ra_center,dec_center
radius=60.*1.05*0.5*sqrt(naxis1^2+naxis2^2)

; get the 2mass catalog
if keyword_set(VERBOSE) then print,'Fetching 2MASS catalog at ',ra_center,dec_center,' radius',radius

tmass=queryvizier('2MASS-PSC',[ra_center,dec_center],radius,/CANADA)

; the catalog itself is a structure
;RAJ2000 DEJ2000 _2MASS JMAG E_JMAG HMAG E_HMAG KMAG E_KMAG 
;  QFLG RFLG BFLG CFLG XFLG AFLG
print,'2mass data acquired.'
nstars=n_elements(tmass._2MASS)

mra=tmass.RAJ2000
mdec=tmass.DECJ2000
j=tmass.JMAG
h=tmass.HMAG
k=tmass.KMAG


endelse

if keyword_set(FILTER) then begin
  
case 1 of
   (filter EQ 3):begin 
             compmag=k
             print,'Comparing vs. 2MASS Ks.'
             end
   (filter EQ 2):begin
             compmag=h
             print,'Comparing vs. 2MASS H.'
             end
   (filter EQ 1): begin
             compmag=j
             print,'Comparing vs. 2MASS J.'
             end
   else: begin
             compmag=k
             print,'No 2MASS filter found, defaulting to K.'
         end  
  endcase
  
endif else begin

case 1 of
   (filter2 EQ 'Ks__(2.15)'):begin 
             compmag=k
             print,'Comparing vs. 2MASS Ks.'
             end
   (filter2 EQ 'H__(1.64)'):begin
             compmag=h
             print,'Comparing vs. 2MASS H.'
             end
   (filter2 EQ 'J__(1.25)'): begin
             compmag=j
             print,'Comparing vs. 2MASS J.'
             end
   else: begin
             compmag=k
             print,'No 2MASS filter found, defaulting to K.'
         end
endcase

endelse
   


n=n_elements(ra)
dist=mra
sep_x=ra
sep_y=ra
sep_xy=ra
fwhm_match=ra
dist[*]=0
match=ra
match[*]=-999
nan=sqrt(-1)
sep_x[*]=nan
sep_y[*]=nan
sep_xy[*]=nan
fwhm_match[*]=nan

; set some objects as bad
valid=ra
valid[*]=0
;good=where((x GT 625)and(x lt 2341)and(y GT 699)and(y lt 2375))
;valid[good]=1
valid[*]=1



; do the matching
print,'Matching 2mass stars....'
for i=0,n-1 do begin

if (valid[i] EQ 1)then begin

dist=sphdist(ra[i],dec[i],mra,mdec,/degrees)

; separation must be less than 4"
sep=min(dist,ind)
if (sep LT (MATCHRADIUS/3600.)) then begin
    match[i]=ind
    sep_xy[i]=sep*3600.
    sep_y[i]=3600.*(dec[i]-mdec[ind])
    dec_mult=cos(dec[i]*3.14159/180.)
    sep_x[i]=dec_mult*3600.*(ra[i]-mra[ind])
    fwhm_match[i]=fwhm[ind]

    print,i,ind,sep*3600.
    endif
    
    
endif    
    
endfor



hist=histogram(sep_y,binsize=0.1,min=-MATCHRADIUS,max=MATCHRADIUS,/nan,locations=xaxis)
plot,xaxis,hist,xtitle='Dec Offset Arcseconds',ytitle='Number'

if keyword_set (PLOTROOT) then begin
     outname=strcompress(PLOTROOT+'_final_decoffset.gif',/REMOVE_ALL)
     write_gif,outname,tvrd()
    endif

print,'Median DEC offset ',median(sep_y),' arcseconds'
crval2=crval2-(median(sep_y)/3600.)
sxaddpar,hd,'CRVAL2',crval2


  print,'Press any key to continue...'
  result=get_kbrd()  
  
hist=histogram(sep_x,binsize=0.1,min=-MATCHRADIUS,max=MATCHRADIUS,/nan,locations=xaxis)
plot,xaxis,hist,xtitle='RA Offset Arcseconds',ytitle='Number'

if keyword_set (PLOTROOT) then begin
     outname=strcompress(PLOTROOT+'_final_raoffset.gif',/REMOVE_ALL)
     write_gif,outname,tvrd()
    endif


print,'Median RA offset ',median(sep_x),' arcseconds'
dec_mult=cos(crval2*3.14159/180.)
crval1=crval1+(median(sep_y)/(3600.*dec_mult))
sxaddpar,hd,'CRVAL1',crval1


  print,'Press any key to continue...'
  result=get_kbrd()  
  
; generate a skyview take file for diagnostic purposes
; print,'Making skyview take file: final_calib.take'
; xr=x(where(match GT -1))
; yr=y(where(match GT -1))
; maketake,'final_calib.take',xr,yr


; match the magnitudes
matchmag=ra
matchmag[*]=nan

; set the output values
for i=0,n-1 do begin
if (match[i] GT -1) then begin
    matchmag[i]=compmag[match[i]]
   endif
endfor



print,'Matched ',n_elements(where(match GT -1)),' stars.'

offset=matchmag-mag_best
k=moment(offset,/nan)
lo=median(offset)-5*(sqrt(k[1]))
hi=median(offset)+5*(sqrt(k[1]))

plot,matchmag,offset,psym=2,yrange=[lo,hi],xtitle='2mass',ytitle='Zeropoint'

good=where(((offset-median(offset)) LT 0.5)and ((offset-median(offset)) GT -0.5)and(matchmag GT 12)and(matchmag LT 15))

k=moment(offset[good])
print,'Mean Zeropoint is ',k[0],' +/- ',sqrt(k[1])
print,'Median Zeropoint is ',median(offset[good])

; draw the mean +3-sigma lines in red
linex=[0,20]
liney=[k(0),k(0)]
oplot,linex,liney,color='0000ff'x
liney=[k(0)+3*sqrt(k[1]),k(0)+3*sqrt(k[1])]
oplot,linex,liney,LINESTYLE=2,color='0000ff'x
liney=[k(0)-3*sqrt(k[1]),k(0)-3*sqrt(k[1])]
oplot,linex,liney,LINESTYLE=2,color='0000ff'x

; draw the median line in green
liney=[median(offset[good]),median(offset[good])]
oplot,linex,liney,color='00ff00'x

sxaddpar,hd,'MEANZPT',k[0],' mean magnitude zeropoint'
sxaddpar,hd,'MEDZPT',median(offset[good]),' median magnitude zeropoint'
sxaddpar,hd,'SIGMAMAG',sqrt(k[1]),' 1-sigma scatter in zeropoint'


if keyword_set (PLOTROOT) then begin
     outname=strcompress(PLOTROOT+'_final_magzpt.gif',/REMOVE_ALL)
     write_gif,outname,tvrd()
    endif
    
    
    
  print,'Press any key to continue...'
  result=get_kbrd()  
  
hist=histogram(fwhm*abs(cdelt[0]*3600.),binsize=0.09,min=0,max=3,/nan,locations=xaxis)
plot,xaxis,hist,xtitle='Measured FWHM (arcseconds)'
  
if keyword_set (PLOTROOT) then begin
     outname=strcompress(PLOTROOT+'_final_seeing.gif',/REMOVE_ALL)
     write_gif,outname,tvrd()
    endif  
  
  
q=sort(fwhm)
q1=fwhm[q[ceil(n_elements(fwhm)*0.25)]]*abs(cdelt[0]*3600.)
print,'Median FWHM: ',median(fwhm)*abs(cdelt[0]*3600.),' arcec'
sxaddpar,hd,'MEDSEE',median(fwhm)*abs(cdelt[0]*3600.),' median observed FWHM (arcsec)'
print,'First quartile observed FWHM: ',q1,' arcsec'
sxaddpar,hd,'Q1SEE',q1,' first quartile observed FWHM (arcsec)'




; get the name for the output file
outbase=strmid(infile,0,strlen(infile)-5)
outname=strcompress(outbase+'_finalcalib.fits',/remove_all)
print,'Writing ',outname
writefits,outname,image,hd



end
