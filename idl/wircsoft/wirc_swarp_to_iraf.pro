; 
; this script reads in a set of swarp reprojected images
;  and then writes them with pixel formatting suitable for IRAF
;
pro wirc_swarp_to_iraf,swarplist

readcol,swarplist,swarp,format="A"
nfiles=n_elements(swarp)




for i=0,nfiles-1 do begin

outbase=strmid(swarp[i],0,strlen(infile[i])-4)

img=readfits(swarp[i],hd)
weight=readfits(outbase+'.weight.fits')


badpix=where(weight EQ 0)
img[badpix]=-9999.


outname=strcompress('iraf_swarp_'+outbase,/remove_all)

writefits,outname,img,hd



endfor





end