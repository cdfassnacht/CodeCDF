;  Generate a dark current model for WIRC
;
;  plane 0 = bias
;  plane 1 = dark current in DN/sec
;
;  FULLAUTO scans the input files for darks
;
;  last modified 4/29/08 JAS
;
pro wirc_makedarkmodel,filelist,$        ; input list of files
                       OUTNAME=OUTNAME,$ ; specify alternative output filename
                       VERBOSE=VERBOSE,$ ; print diagnostic messages
                       FULLAUTO=FULLAUTO ; find and use darks automatically


readcol,filelist,infile,format="A"
nfiles=n_elements(infile)
head=headfits(infile[0])
naxis1=sxpar(head,'NAXIS1')
naxis2=sxpar(head,'NAXIS2')
darks=fltarr(nfiles)
darks[*]=0


; if running in automatic mode
if keyword_set(FULLAUTO) then begin

    for i=0,nfiles-1 do begin ; figure out which frames are darks and coadds=1
   
          head=headfits(infile[i])
          filter1=sxpar(head,'FORE')
          filter2=sxpar(head,'AFT')
          coadds=sxpar(head,'COADDS')
     if ((filter1 EQ 'BrGamma__(2.17)')and(filter2 EQ 'J__(1.25)')and(coadds eq 1)) $
          then darks[i]=1
   
   
   endfor
endif else begin
     darks[*]=1
endelse

; find how many input files were darks
count=where(darks EQ 1,ndarks)
if keyword_set(VERBOSE) then print,'Found ',ndarks,' out of ',nfiles


; initialize some storage arrays
temp=fltarr(naxis1,naxis2,ndarks)
out=fltarr(naxis1,naxis2,2)
times=fltarr(ndarks)

; make the image header
mkhdr,out_hd,out
sxaddhist,'Processed by wirc_makedarkmodel',out_hd
sxaddhist,'Plane 0 = bias',out_hd
sxaddhist,'Plane 1 = dark current (dn/sec)',out_hd

; load the images
counter=0
for i=0,nfiles-1 do begin

   if (darks[i] EQ 1) then begin
      if keyword_set(VERBOSE) then print,'Reading dark ',infile[i]
      temp[*,*,counter]=readfits(infile[i],hd,/silent)
      coadds=sxpar(hd,'COADDS')
      if (COADDS NE 1) then temp[*,*,counter]=temp[*,*,counter]/coadds
      times[counter]=sxpar(hd,'EXPTIME')
      sxaddhist,infile[i],out_hd
      counter=counter+1
   endif
   
endfor

;stop

; now collapse and fit the data
if keyword_set(VERBOSE) then print,'Fitting dark model.'
for i=0,naxis1-1 do begin
   for j=0,naxis2-1 do begin
   
   result=linfit(times,temp[i,j,*])
   out[i,j,*]=result
  
   
   endfor
  endfor


; write the result
if (keyword_set(outname) NE 1) then outname='darkmodel.fits'
if keyword_set(VERBOSE) then print,'Writing ',outname

writefits,outname,out,out_hd

end
