pro wirc_plotmedian,inlist

readcol,inlist,infile,format="A"

nfiles=n_elements(infile)
med=fltarr(nfiles)

for i=0,nfiles-1 do begin


    im=readfits(infile[i],/silent)
    med[i]=median(im)
    print,med[i]

endfor

stop

end
