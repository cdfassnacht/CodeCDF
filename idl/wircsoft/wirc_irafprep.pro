; prepare the images for processing by IRAF
;
;
;  060108 Added reprojection option JAS
;  061008 Added fixzp option for combining multiple days
;
pro wirc_irafprep,inlist,$
                  MASKLIST=MASKLIST,$ ; optional list of the object mask files
                  COVMAP=COVMAP,$ ; produce a flag image for a coverage map
                  REPROJECT=REPROJECT,$ ; use SWARP to reproject images
                  FLATCUT=FLATCUT,$ ; level in flat below which to nan data
                  SWARPFILE=SWARPFILE,$ ; name of swarp coadd to keep
                  FIXZP=FIXZP ; force scaling factor to be magzpt=25


; let's start by processing the command-line defaults

; set FLATCUT equal to a value if not set
if (keyword_set(FLATCUT) NE 1) then FLATCUT=0.7

; set the initial bad pixel value
if keyword_set(REPROJECT) then badpixval=sqrt(-1) else badpixval=-99999

; read in the input file list
readcol,inlist,infile,format="A"
nfiles=n_elements(infile)

; read some header items from the first image
hd=headfits(infile[0])
exptime=sxpar(hd,'EXPTIME')
filter2=sxpar(hd,'AFT')
print,'Filter and exposure time are ',filter2,exptime

; if a list of mask files is present, read them too
if keyword_set(MASKLIST) then readcol,MASKLIST,maskfile,format="A"


; read in the appropriate flat
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
 
 
 
 
; load the zeropoint array
magzpt=fltarr(nfiles)
print,'Loading magnitude zeropoints.'
for i=0,nfiles-1 do begin
     hd=headfits(infile[i])
     magzpt[i]=sxpar(hd,'MAGZPT')
endfor
 
; and print some stats about it
medianzpt=median(magzpt)
print,'Median magnitude zeropoint is ',medianzpt
k=moment(magzpt)
print,'Average is ',k[0],'+/-',sqrt(k[1])


; create the bad data mask array
;  0=good, 1=bad
bad=fltarr(2048,2048)
bad[*]=0

  
; now to kill pixels based on the flat
; first, kill all nan'd or infinite points
flatbad=where(finite(flat,/nan) or finite(flat,/infinity),count)
if (count GT 0) then bad(flatbad)=1
; then kill any data in low sensitivity regions
flatbad=where(flat LT FLATCUT,count)
if (count GT 0) then bad(flatbad)=1
 
 
; load the zeropoint array
print,'Loading magnitude zeropoints.'
for i=0,nfiles-1 do begin
     hd=headfits(infile[i])
     magzpt[i]=sxpar(hd,'MAGZPT')
endfor
 
; and print some stats about it
medianzpt=median(magzpt)
print,'Median magnitude zeropoint is ',medianzpt
k=moment(magzpt)
print,'Average is ',k[0],'+/-',sqrt(k[1])




; if reprojection option is selected, do it
if keyword_set(REPROJECT) then begin

; check to see if the temporary directory exists
; then make it if it doesn't
test=file_test('prep_temp')
if (test EQ 1) then spawn,'rm -r prep_temp'
spawn,'mkdir prep_temp'


endif



 
 
; read in the images and start killing pixels
;  by setting them to an absurdly low level
for i=0,nfiles-1 do begin

;read in the image
image=readfits(infile[i],hd,/silent)
    
; optionally read in the mask file
; this will improve the image statistics
if keyword_set(MASKLIST) then begin
    mask=readfits(maskfile[i],mask_hd,/silent)
   endif
       
; zap anything bad in the flat
;  set up above
badpix=where((bad NE 0),count)
if (count GT 0) then image(badpix)=badpixval
    
; kill non-numbers in the processed data images
badpix=where(finite(image,/infinity),count)
if (count GT 0) then image(badpix)=badpixval
badpix=where(finite(image,/nan),count)
if (count GT 0) then image(badpix)=badpixval
      

    
; compute irafscl and iraf weight

if keyword_set(FIXZP) then begin
        irafscl=10^((magzpt[i]-25.0)/(-2.5))
endif else begin
        irafscl=10^((magzpt[i]-medianzpt)/(-2.5))
        endelse
        
print,'Frame level MagZpt is ',magzpt[i]
print,'IRAFSCL is ',irafscl
sxaddpar,hd,'IRAFSCL',irafscl
sxaddpar,hd,'IRAFWGT',(1./irafscl)
    
; get the background pixels
if keyword_set(MASKLIST) then begin
      background=image(where((mask EQ 0) and (image GT -90000)))
   endif else begin
      background=image(where(image GT -90000))
   endelse
      
print,n_elements(background),' background pixels.'    
k=moment(background,/nan)
print,'Median background ',median(background)
sxaddpar,hd,'BACKMED',median(background)
sxaddpar,hd,'BACKMEAN',k[0]
sxaddpar,hd,'BACKSIG',sqrt(k[1])
    

; this converts the data to 0 or the total number of images
; use this option to make a set of images that can be averaged together to
;    make a coverage map (in image units)
if (keyword_set(COVMAP) and (keyword_set(REPROJECT) NE 1))then begin
     image(where(image GT -900))=1.
     image(where(image LT -900))=0
   endif
   
; this version is for using sextractor   
if (keyword_set(COVMAP) and keyword_set(REPROJECT))then begin
     badpix=where(image NE image)
     image[*]=1.
     image[badpix]=sqrt(-1)
   endif
   
   

; get the output filename basename
if (strlowcase(strmid(infile[i],2,3,/reverse)) EQ '.gz' ) then begin
       outbase=strmid(infile[i],0,strlen(infile[i])-3)
   endif else begin
       outbase=infile[i]
   endelse

    
; write the prepped file  
if keyword_set(REPROJECT) then begin
              outname=strcompress('prep_temp/prep_'+outbase,/remove_all)
          endif else begin
              outname=strcompress('prep_'+outbase,/remove_all)
          endelse
print,'Writing ',outname
writefits,outname,image,hd

    
endfor



; now add a bunch more code to swarp the files

if keyword_set(REPROJECT) then begin

print,'Reprojecting images'

; now we must run swarp

; first make a list of the files in the temp directory
spawn,'ls prep_temp/prep*fits > prep_temp.list'

; now swarp them
if keyword_set(COVMAP) then begin
    spawn,"swarp prep_temp/*fits -c sex_files/wirc_irafprep_covmap.swarp"
endif else begin
    spawn,"swarp prep_temp/*fits -c sex_files/wirc_irafprep.swarp"
endelse

; make new filelists
spawn,"ls prep_temp/*resamp.fits > prep_temp2.list"
spawn,"ls prep_temp/*resamp.weight.fits > prep_temp3.list"


readcol,'prep_temp2.list',resamp,format="A"
readcol,'prep_temp3.list',weight,format="A"

for i=0,n_elements(resamp)-1 do begin

img=readfits(resamp[i],hd,/silent)
img2=readfits(weight[i],/silent)

badpix=where(img2 EQ 0,count)


if ((count GT 0) and keyword_set(COVMAP)) then img[badpix]=0

if ((count GT 0) and (keyword_set(COVMAP) NE 1)) then img[badpix]=-9999

outbase=strmid(resamp[i],10,strlen(resamp[i])-22)
outname=strcompress(outbase+'.fits',/remove_all)
print,'Writing ',outname
writefits,outname,img,hd


endfor


; if specified, save the swarp file
if keyword_set(SWARPFILE) then begin
   csh_command='cp prep_temp/coadd.fits '+SWARPFILE
   spawn,csh_command
endif


; delete temp files
spawn,'rm prep_temp.list'
spawn,'rm prep_temp2.list'
spawn,'rm prep_temp3.list'
spawn,'rm -r prep_temp'

endif





end
