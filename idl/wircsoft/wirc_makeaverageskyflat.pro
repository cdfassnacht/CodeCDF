; make a skyflat
pro wirc_makeaverageskyflat,inlist,$
                     masklist,$
                     DARKFILE=DARKFILE,$
                     VERBOSE=VERBOSE

; read the input filelist 
readcol,inlist,infile,format="A"
readcol,masklist,maskfile,format="A"
nfiles=n_elements(infile)


; read in the first image header
hd=headfits(infile[0],/silent)
naxis1=sxpar(hd,'NAXIS1')
naxis2=sxpar(hd,'NAXIS2')
filter2=sxpar(hd,'AFT')
exptime=sxpar(hd,'EXPTIME')
coadds=sxpar(hd,'COADDS')


accum=fltarr(naxis1,naxis2)
n=fltarr(naxis1,naxis2)
accum[*]=0.
n[*]=0.



; automatic filenaming
if (keyword_set(OUTNAME) NE 1) then begin
case 1 of
   (filter2 EQ 'Ks__(2.15)'):begin 
             outname='kskyflat.fits'
             end
   (filter2 EQ 'H__(1.64)'):begin
             outname='hskyflat.fits'
             end
   (filter2 EQ 'J__(1.25)'): begin
             outname='jskyflat.fits'
             end
   else: begin
             print,'I do no know what to call the output!'
             print,'Defaulting to skyflat.fits.'
             outname='skyflat.fits'
         end
   endcase
endif 

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

; load the temp array
for k=0,nfiles-1 do begin

img=readfits(infile[i])
mask=readfits(maskfile[i])
badpix=where(mask NE 0,count)
if (count GT 0) then img[badpix]=0.
accum=accum+img

goodpix=where(mask EQ 0,count)
if (count GT 0) then n[goodpix]=n[goodpix]+1.

endfor



; take the average
out=accum/n


out=out-dark
out=out/median(out)

if keyword_set(VERBOSE) then print,'Writing ',outname
writefits,outname,out,hd



end