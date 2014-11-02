; prepare the images for processing by SWARP
pro wirc_swarpprep,inlist

nan=sqrt(-1)

readcol,inlist,infile,format="A"
  nfiles=n_elements(infile)
  hd=headfits(infile[0])
  exptime=sxpar(hd,'EXPTIME')
  filter2=sxpar(hd,'AFT')
  print,filter2,exptime

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
 
 bad=fltarr(2048,2048)
 bad[*]=0
 flatbad=fltarr(2048,2048)
 flatbad[*]=0
 magzpt=fltarr(nfiles)
 
 
 ; set bad pixels based on the flat
 if (n_elements(where(flat NE flat)) GT 1) then flatbad(where(flat ne flat))=1
 flatbad(where(flat LT 0.7))=1
 
 bad=bad+flatbad
 
 ; load the zeropoint array
 print,'Loading magnitude zeropoints.'
 for i=0,nfiles-1 do begin
     hd=headfits(infile[i])
     magzpt[i]=sxpar(hd,'MAGZPT')
 endfor
 
 medianzpt=median(magzpt)
 print,'Median magnitude zeropoint is ',medianzpt

 
weight=flat
weight(where(flatbad NE 0))=0

 
 
; read in the images and massage before passing to IRAF 
for i=0,nfiles-1 do begin

    image=readfits(infile[i],hd)
    image(where(flatbad NE 0))=nan
    
    ; compute scales and weights
    irafscl=10^((magzpt[i]-medianzpt)/(-2.5))
    print,'FLXSCALE is ',irafscl
    sxaddpar,hd,'FLXSCALE',irafscl
    sxaddpar,hd,'IRAFWGT',(1./irafscl)
    
    
    outname=strcompress('prep_'+infile[i],/remove_all)
    writefits,outname,image,hd
    
    
    parts=strsplit(outname,'.',/extract)
    outname=parts[0]+'.weight.'+parts[1]
    outname=strcompress(outname,/remove_all)
    writefits,outname,weight*(1./irafscl)
    
    
endfor


end