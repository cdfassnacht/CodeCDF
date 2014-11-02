; adjust the wirc wcs from runningskysub
;
;
;
pro wirc_fixwcs,inlist,catalog,$
                OFFSETS=OFFSETS,$   ; 2-element vector of pixel offsets
                VERBOSE=VERBOSE,$   ; print diagnostic messages
                NOFOLLOW=NOFOLLOW,$ ; adjust the offsets frame-to-frame
                COMPRESS=COMPRESS   ; compress the output

; handle the presets
if (keyword_set(OFFSETS) NE 1) then OFFSETS=[0.,0.]

; read in the input file list
readcol,inlist,infile,format="A"

; read the 2MASS catalog
readcol,catalog,mra,mdec,j,h,k


; determine which magnitude to compare against
hd=headfits(infile[0])
filter2=sxpar(hd,'AFT')

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



nfiles=n_elements(infile)
nstars=n_elements(mra)
match=intarr(nstars)
matchdist=fltarr(nstars)
match_raoffset=fltarr(nstars)
match_decoffset=fltarr(nstars)
matchmag=fltarr(nstars)
onimage=intarr(nstars)
nan=sqrt(-1)

for j=0,nfiles-1 do begin

print,'--------------------------------------------------------------------'
print,infile[j]
print,' '
image=readfits(infile[j],hd,/silent)
exptime=sxpar(hd,'EXPTIME')
filter2=sxpar(hd,'AFT')
coadds=sxpar(hd,'COADDS')
lst=sxpar(hd,'LST')
tmpra= sxpar(hd,'CRVAL1')
tmpdec= sxpar(hd,'CRVAL2')
get_coords,timetemp,instring=lst+' 0:0:0'

; print some info about the input image        
print,filter2,exptime,coadds
print,'Local Sidereal Time: ',timetemp[0],' hours'
print,'Initial RA, Dec of file (from telescope pointing)',tmpra,tmpdec

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
writefits,'sex_temp.fits',sex_image,hd
spawn,'/usr/local/bin/sex sex_temp.fits -c sex_files/fixwcs.sex > sex.log'


; this checks for validity of the sextractor solution
sex_detect=file_test('fixwcs.txt',/zero_length)

; read the file if non-zero length
if (sex_detect NE 1) then readcol,'fixwcs.txt',cntr,x,y,magbest,magerr,threshold,flags,fwhm,/silent
print,'Number of sextractor detected objects: '+string(n_elements(cntr))
if (n_elements(cntr) LT 10) then sex_detect=1 ; if < 10 declare failure


; start the big loop
if (sex_detect NE 1) then begin


;let's get the fwhm
; must do this before applying offsets
goodfwhm=where((x GT 200) and (x LT 1800) and $
               (y GT 200) and (x LT 1800) and $
               (flags EQ 0),numpix)
if numpix GT 0 then medianfwhm=0.248*median(fwhm(goodfwhm)) else medianfwhm=-99.
print,'Median seeing: ',medianfwhm,' arcseconds.'



; note that offsets are the values to move in X,Y to the real star
;   position from the current position
; the "1" is needed to fix difference between sextractor and xyad
; 
; note that some offsets priobably still remain
;
x=x+OFFSETS[0]-1.
y=y+OFFSETS[1]-1.

; convert sextractor list to RA/DEC
xyad,hd,x,y,ra,dec


; figure out which catalog stars are on the image
adxy,hd,mra,mdec,mx,my
onimage[*]=0
onimage(where((mx LT 2000)and(mx GT 48)and(my LT 2000)and(my GT 48)))=1
print,'Number of potential match stars: '+string(n_elements(where(onimage EQ 1)))

matchmag[*]=nan
match[*]=-99
matchdist[*]=nan
match_raoffset[*]=nan
match_decoffset[*]=nan

for i=0,nstars-1 do begin


 if (onimage[i] EQ 1)then begin

     dist=sphdist(mra[i],mdec[i],ra,dec,/degrees)
     sep=min(dist,ind)
     match[i]=ind
     matchdist[i]=sep
     match_raoffset[i]=mra[i]-ra[ind]
     match_decoffset[i]=mdec[i]-dec[ind]
     matchmag[i]=compmag[i]-magbest[ind]

     endif


endfor

k=moment(match_raoffset,/nan)*3600.
print,'RA Median, Mean, Sigma offset (arcsec) of local fit:'
print,median(match_raoffset)*3600.,k[0],sqrt(k[1])
ravar=k[1]
print,'DEC Median, Mean, Sigma offset (arcsec) of local fit:'
k=moment(match_decoffset,/nan)*3600.
print,median(match_decoffset)*3600.,k[0],sqrt(k[1])
decvar=k[1]
print,'Magnitude zeropoint: '+string(median(matchmag))

; write mag zeropoint data
matchloc=where(match GT -1)
magnum=n_elements(matchloc)
if (matchloc[0] EQ -1) then magnum=0
sxaddpar,hd,'MAGNUM',magnum

if (magnum LT 3) then begin
          magzpt=-99 
          print,'Match to 2MASS FAILED!!!'
     endif else begin
          magzpt=median(matchmag)
          sxaddpar,hd,'MAGZPT',magzpt,' Magnitude zeropoint relative to 2MASS'
     endelse




crval1=sxpar(hd,'CRVAL1')
crval2=sxpar(hd,'CRVAL2')
;cdelt1=sxpar(hd,'CDELT1')
;cdelt2=sxpar(hd,'CDELT2')
cdelt1=-1.0*sxpar(hd,'PIXSCAL1')
cdelt2=sxpar(hd,'PIXSCAL2')


if (magnum GT 3) then begin
  crval2=crval2+(OFFSETS[1]*cdelt2)+median(match_decoffset)
  dec_mult=cos(crval2*3.14159/180.)
  crval1=crval1+(OFFSETS[0]*cdelt1/dec_mult)+median(match_raoffset)
  print,'Total offsets ',((OFFSETS[0]*cdelt1/dec_mult)+median(match_raoffset))*3600.,$
                   ((OFFSETS[1]*cdelt2)+median(match_decoffset))*3600.,$
                   ' arcseconds.'
             
  sxaddpar,hd,'CRVAL1',crval1
  sxaddpar,hd,'CRVAL2',crval2
  sxaddpar,hd,'PNTREF',1,' Pointing Refinement Applied'
  sxaddpar,hd,'PNTSTARS',magnum,' Number of stars used'
  sxaddpar,hd,'PNTERR',sqrt(ravar+decvar),' Pointing error in arcseconds'
sxaddpar,hd,'RAOFF',((OFFSETS[0]*cdelt1/dec_mult)+median(match_raoffset))*3600.,' RA Offset in arcseconds' 
  sxaddpar,hd,'DECOFF',((OFFSETS[1]*cdelt2)+median(match_decoffset))*3600.,' DEC Offset in arcseconds'
  sxaddpar,hd,'SEEING',medianfwhm,' Median FWHM in arcseconds'
endif

if (magnum LE 3) then sxaddpar,hd,'PNTREF',0,' Pointing Refinement NOT Applied'


; derive the output name
; get the output file basename
if (strlowcase(strmid(infile[j],2,3,/reverse)) EQ '.gz' ) then begin
       outbase=strmid(infile[j],0,strlen(infile[j])-3)
    endif else begin
       outbase=infile[j]
    endelse
    
outname=strcompress('wcs_'+outbase,/remove_all)

; write the output
if keyword_set(COMPRESS) then begin
          writefits,outname,image,hd,/COMPRESS
      endif else begin
          writefits,outname,image,hd
      endelse


; adjust the offsets
;  so that the code follows it's own solutions, applying deltas
;
if (keyword_set(NOFOLLOW) NE 1)then begin
       print,'Before offsets (pixels) ',offsets
       offsets[1]=OFFSETS[1]+median(match_decoffset)/cdelt2
       offsets[0]=OFFSETS[0]+(median(match_raoffset)*dec_mult/cdelt1)
       print,'After offsets (pixels) ',offsets
    endif


endif


endfor


end
