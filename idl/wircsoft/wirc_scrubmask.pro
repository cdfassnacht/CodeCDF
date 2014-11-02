; takes the bad pixels _out_ of the mask
pro wirc_scrubmask,masklist


readcol,masklist,maskfile,format="a"
nfiles=n_elements(maskfile)

count=fltarr(2048,2048)
count[*]=0


for i=0,nfiles-1 do begin

    print,'Reading ',maskfile[i]
    mask=readfits(maskfile[i])
    count=count+mask

endfor

turnoff=where(count GT (nfiles/2))
turnoffnum=n_elements(turnoff)
print,turnoffnum

if (turnoff[0] NE -1) then begin

print,'Fixing masks.'

; now fix the masks
for i=0,nfiles-1 do begin


   ; turn pixels off
   mask=readfits(maskfile[i],hd)
   mask(turnoff)=0
 
   sxaddhist,'wirc_scrubmask run',hd
   writefits,maskfile[i],mask,hd


endfor

endif


end