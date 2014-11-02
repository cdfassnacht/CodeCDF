; simple median stack and difference engine
;  for making dome-flats
;
;
;  5/1/08 totally rewritten to ease dataflow  JAS
;
pro wirc_makeflat,inlist,$
                   dome_on,$ ; 2-element vector
                   dome_off,$; 2-element vector
                   OUTNAME=OUTNAME,$
                   VERBOSE=VERBOSE

; read the input filelist (from wirc_logsheet)
readcol,inlist,infile,format="A"



; read in the first image header
hd=headfits(infile[dome_on[0]],/silent)
naxis1=sxpar(hd,'NAXIS1')
naxis2=sxpar(hd,'NAXIS2')
filter2=sxpar(hd,'AFT')

tempon=fltarr(naxis1,naxis2,dome_on[1]-dome_on[0]+1)
tempoff=fltarr(naxis1,naxis2,dome_off[1]-dome_off[0]+1)
out=fltarr(naxis1,naxis2)


; automatic filenaming
if (keyword_set(OUTNAME) NE 1) then begin
case 1 of
   (filter2 EQ 'Ks__(2.15)'):begin 
             outname='kflat.fits'
             end
   (filter2 EQ 'H__(1.64)'):begin
             outname='hflat.fits'
             end
   (filter2 EQ 'J__(1.25)'): begin
             outname='jflat.fits'
             end
   else: begin
             print,'I do no know what to call the output!'
             print,'Defaulting to flat.fits.'
             outname='flat.fits'
         end
   endcase
endif 

; add some header keywords
history='Flat made by wirc_makeflat'
sxaddhist,history,hd
history=string(2*(dome_on[1]-dome_on[0]+1))+' files used.'
sxaddhist,history,hd



; load the arrays

; load the dome-on array
if keyword_set(VERBOSE) then print,'Dome On:'
sxaddhist,'Dome on:',hd
counter=0
for i=dome_on[0],dome_on[1] do begin
   tempon[*,*,counter]=readfits(infile[i],/silent)
   sxaddhist,infile[i],hd
   if keyword_set(VERBOSE) then begin
           med=median(tempon[*,*,counter])
           print,infile[i],med
           endif
   counter=counter+1
   endfor


; now for the off array
if keyword_set(VERBOSE) then print,'Dome Off:'
sxaddhist,'Dome off:',hd
counter=0
for i=dome_off[0],dome_off[1] do begin
   tempoff[*,*,counter]=readfits(infile[i],/silent)
   sxaddhist,infile[i],hd
   if keyword_set(VERBOSE) then begin
        med=median(tempoff[*,*,counter])
        print,infile[i],med
        endif
   counter=counter+1
   endfor


if keyword_set(VERBOSE) then print,'Computing Flat.'
for i=0,naxis1-1 do begin
  for j=0,naxis2-1 do begin
  
  out[i,j]=median(tempon[i,j,*])-median(tempoff[i,j,*])
  
  endfor
endfor

; normalize to a median of 1
out=out/median(out)

if keyword_set(VERBOSE) then print,'Writing ',outname
writefits,outname,out,hd



end