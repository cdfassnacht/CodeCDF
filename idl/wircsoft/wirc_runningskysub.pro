; note - timeradius is equal to the half-width
;  went for IO-intensive solution
;  otherwise memory demands might be too high to be accomodated.
;
; cognoscenti will recognize descent from dim-sum
;
; 4/29/08  cal file handling re-written  JAS
; 5/01/08  added file compression option
;          stopped saving sky frames by default
; 5/15/08  added some logic for multiple coadds JAS
; 5/29/08  fixed a problem in the where statement for sextractor
;
;
pro wirc_runningskysub,inlist,$
                       timeradius,$ ; radius (in images) forward and back
                       RUN=RUN,$ ; should be 1 or 2
                       FLATFILE=FLATFILE,$ ; alternative flatfield
                       DARKFILE=DARKFILE,$ ; alternative dark
                       NOQUADCLEAN=NOQUADCLEAN,$ ; don't run quadrant scrubber
                       SAVESKY=SAVESKY,$ ; save the sky files
                       COMPRESS=COMPRESS ; gzip all output files

; the RUN keyword indicates if this is the first pass for mask generation
if (keyword_set(RUN) NE 1) then RUN=1


; get the time
start_time=systime(/seconds)

; convolution kernel for growing the masks
kernel=fltarr(5,5)
kernel[*]=1

; read in the input list
readcol,inlist,infile,format="a"
nfiles=n_elements(infile)
; this list will be manipulated to decide which files get included 
list=fltarr(nfiles)
list[*]=0

; get the calibration frames
hd=headfits(infile[0])
exptime=sxpar(hd,'EXPTIME')
filter2=sxpar(hd,'AFT')
coadds=sxpar(hd,'COADDS')
print,filter2,exptime


; determine which flatfield to use
if (keyword_set(FLATFILE) NE 1) then begin
case 1 of
   (filter2 EQ 'Ks__(2.15)'):begin 
             flat=readfits('kflat.fits')
             print,'Kflat selected.'
             end
   (filter2 EQ 'H__(1.64)'):begin
             flat=readfits('hflat.fits')
             print,'Hflat selected.'
             end
   (filter2 EQ 'J__(1.25)'): begin
             flat=readfits('jflat.fits')
             print,'Jflat selected.'
             end
   else: begin
         flat=fltarr(2048,2048)
         flat[*,*]=1.
         print,'No flat selected.'
         end
   endcase
endif else begin
   flat=readfits(FLATFILE)
   endelse
 
 
 
; find the right dark calibration file
if (keyword_set(DARKFILE) NE 1) then begin 

; get the automatic dark files and test for existence
darkfile=strcompress('dark_'+string(floor(exptime))+'_'+string(coadds)+'.fits',/remove_all)
filetest=file_test(darkfile)
filetest2=file_test('darkmodel.fits')

if (filetest EQ 1) then begin
      dark=readfits(darkfile)
      print,'Reading ',darkfile
      endif
      
if ((filetest NE 1) and(filetest2 EQ 1)) then begin
      print,'Can not find matching exptime dark!'
      print,'Falling back to dark current model.'
      darkmodel=readfits('darkmodel.fits')  
      dark=darkmodel[*,*,0]+exptime*coadds*darkmodel[*,*,1] 
    endif
    
if ((filetest NE 1) and (filetest2 NE 1)) then begin
      print,' No dark info found!'
      print,' No dark will be used.'
      dark=fltarr(2048,2048)
      dark[*]=0
   endif
      
endif else begin
   dark=readfits(DARKFILE)
   print,'Using manually specified dark.'
   print,'Reading ',darkfile
   endelse

; make the running sky sub frame
for i=0,nfiles-1 do begin

list[*]=0

; figure out infile subrange
case 1 of
  (i LT timeradius): list[0:(timeradius*2)]=1
  (i GT (nfiles-timeradius-1)): list[(nfiles-1-2*timeradius):nfiles-1]=1
   else: list[(i-timeradius):(i+timeradius)]=1
endcase
      
;this eliminates the current image from being included in the stack
list[i]=0
      
print,'--------------'      
print,i,'    ',infile[i] 
print,'--------------'  
print,infile[where(list EQ 1)]
      
;load the background images
numimages=n_elements(where(list EQ 1))
images=fltarr(2048,2048,numimages)
masks=fltarr(2048,2048,numimages)
masks[*,*,*]=0

count=0
for j=0,nfiles-1 do begin

    if (list[j] EQ 1) then begin
           images[*,*,count]=readfits(infile[j])
           filename=strcompress('mask_'+infile[j],/remove_all)
           if (RUN EQ 2) then masks[*,*,count]=convol(readfits(filename),kernel)
           count=count+1
           endif
           
endfor

; now read in the actual data
image2=readfits(infile[i],hd)-dark
data_background=median(image2)


; condition the input sky stack a little by
;    first subtracting the dark, then resetting the medians
;    equal to that in the input data frame
for j=0,numimages-1 do begin
        images[*,*,j]=images[*,*,j]-dark
        images[*,*,j]=images[*,*,j]*(data_background/median(images[*,*,j]))
endfor

; zap the input images based on the masks
if (RUN EQ 2) then images(where(masks NE 0))=sqrt(-1)

; generate the output array
out=fltarr(2048,2048)
out[*]=0


; generate the median sky - note that median rejects all the nans
;for j=0,2047 do begin
;  for k=0,2047 do begin

;     out[j,k]=median(images[j,k,*])

; endfor
;endfor
out=median(images,DIMENSION=3)



; get the output file basename
if (strlowcase(strmid(infile[i],2,3,/reverse)) EQ '.gz' ) then begin
       outbase=strmid(infile[i],0,strlen(infile[i])-3)
    endif else begin
       outbase=infile[i]
    endelse

; finally, write the output sky image
if keyword_set(SAVESKY) then begin
outname=strcompress('sky_'+outbase,/remove_all)
if keyword_set(COMPRESS) then begin
             writefits,outname,out,hd,/COMPRESS
         endif else begin
             writefits,outname,out,hd
         endelse
endif


; add the basic WCS
ra=sxpar(hd,'RA')   ; from the TCS
dec=sxpar(hd,'DEC')
get_coords,coords,instring=(ra+' '+dec),/quiet
sxaddpar,hd,'CRVAL1',coords[0]*15
sxaddpar,hd,'CRVAL2',coords[1]
sxaddpar,hd,'CRPIX1',1024.0
sxaddpar,hd,'CRPIX2',1024.0
sxaddpar,hd,'CROTA2',1.4 ; from Thompson
sxaddpar,hd,'CDELT1',(-0.2487/3600.) ; from Thompson
sxaddpar,hd,'CDELT2',(0.2487/3600.)
sxaddpar,hd,'CTYPE1','RA---TAN'
sxaddpar,hd,'CTYPE2','DEC--TAN'


; write the background to the header
rawmed=median(image2)
sxaddpar,hd,'RAWMED',rawmed,' median value of raw image'


; now make the actual data image
out=(image2-out)/flat

; scrub it with the quadrant scrubber
out2=fltarr(2048,2048)
out2[*]=0
if (keyword_set(NOQUADCLEAN) NE 1) then begin
              print,'Cleaning quadrants.'
              wirc_cleanquadrants,out,out2
          endif else begin
              out2=out
          endelse

; and write the output
outname=strcompress('ff_'+outbase,/remove_all)
sxaddpar,hd,'WIRCTRAD',timeradius
sxaddpar,hd,'WIRCRUN',RUN
if keyword_set(COMPRESS) then begin
          writefits,outname,out2,hd,/COMPRESS
  endif else begin
          writefits,outname,out2,hd
  endelse


; if this were the mask generating run
if (RUN EQ 1) then begin    
    ; run sextractor
      badpix=where(finite(out2,/infinity),badpixcount)
         if badpixcount GT 0 then out2[badpix]=999999
      badpix=where(finite(out2,/nan),badpixcount)
         if badpixcount GT 0 then out2[badpix]=999999
      writefits,'sex_temp.fits',out2,hd
      print,'Running SExtractor'
      spawn,'/usr/local/bin/sex sex_temp.fits -c sex_files/preproc_mask.sex '
      mask=readfits('check.fits',mask_hd)



; do a little more processing
      mask(where(flat LE 0.6))=0
      maskpixels=where(mask GT 0)
      if (maskpixels[0] NE -1) then  mask(maskpixels)=1 
      sxaddpar,hd,'WIRCTRAD',timeradius
      sxaddpar,mask_hd,'WIRCMSK',1
      outname=strcompress('mask_'+outbase,/remove_all)
      if keyword_set(COMPRESS) then begin
            writefits,outname,mask,mask_hd,/COMPRESS
         endif else begin
            writefits,outname,mask,mask_hd
         endelse

endif


endfor

; print the total time
stop_time=systime(/seconds)
print,(stop_time-start_time)/60.,' minutes elapsed.'


end
