; simple median stack and difference engine
;  for making dome-flats
;
pro wirc_makeflat_fromlists,onlist,offlist,outflat

readcol,onlist,on_file,format="A"
readcol,offlist,off_file,format="A"


bcd=readfits(on_file[0],hd)
naxis1=sxpar(hd,'NAXIS1')
naxis2=sxpar(hd,'NAXIS2')
tempon=fltarr(naxis1,naxis2,n_elements(on_file))
tempoff=fltarr(naxis1,naxis2,n_elements(off_file))
out=fltarr(naxis1,naxis2)


; load the arrays

for i=0,n_elements(on_file)-1 do begin
   tempon[*,*,i]=readfits(on_file[i])
   endfor



for i=0,n_elements(off_file)-1 do begin
   tempoff[*,*,i]=readfits(off_file[i])
   endfor


for i=0,naxis1-1 do begin
  for j=0,naxis2-1 do begin
  
  out[i,j]=median(tempon[i,j,*])-median(tempoff[i,j,*])
  
  endfor
endfor


out=out/median(out)

writefits,outflat,out



end