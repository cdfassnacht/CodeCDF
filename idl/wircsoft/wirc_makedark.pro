; median stack a dark image 
;
;
; 5/15/08 added # of coadds to the filenames  JAS 
;
pro wirc_makedark,inlist,$ ; list of input fits files
                  start_index,$ ; index number from logsheet.txt
                  stop_index,$ ; index number from logsheet.txt
                  OUTNAME=OUTNAME,$ ; override the default output name
                  VERBOSE=VERBOSE ; print diagnostic messages

readcol,inlist,infile,format="A"

head=headfits(infile[start_index])
naxis1=sxpar(head,'NAXIS1')
naxis2=sxpar(head,'NAXIS2')
exptime=sxpar(head,'EXPTIME')
coadds=sxpar(head,'COADDS')

; if no output specified, derive a name
if (keyword_set(OUTNAME) NE 1) then $  
        outname=strcompress('dark_'+string(floor(exptime))+'_'+string(coadds)+'.fits',/remove_all)
       


nfiles=stop_index-start_index+1
temp=fltarr(naxis1,naxis2,nfiles)
out=fltarr(naxis1,naxis2)
temp[*,*,*]=0

; add some header keywords
history='Dark made by wirc_makedark'
sxaddhist,history,head
history=string(nfiles)+' files used.'
sxaddhist,history,head


; read in the files
for i=0,nfiles-1 do begin
     temp[*,*,i]=readfits(infile[start_index+i],/silent)
     sxaddhist,infile[start_index+i],head
     if keyword_set(VERBOSE) then print,infile[start_index+i]
endfor

for i=0,naxis1-1 do begin
   for j=0,naxis2-1 do begin
       out[i,j]=median(temp[i,j,*])
       endfor
       endfor
       
       
if keyword_set(VERBOSE) then print,'Writing ',outname   
writefits,outname,out,head



end