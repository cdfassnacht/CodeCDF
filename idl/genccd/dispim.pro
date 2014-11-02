pro dispim, img, inhead=inhead, lsig=lsig, hsig=hsig, absmin=absmin, $
   absmax=absmax, imgsz=imgsz,  $
   psfile=psfile, jpfile=jpfile, _extra=_extra

;
; Procedure dispim
;
; Description: Given an input image, displays the image, which also can be 
;               exported into either a postscript file, a jpeg file, or both.
;              Just like mydisp.pro, but with the input being an image
;               rather than a fits file.
;
; Inputs: image       (floatarr)   input image
;         [inhead=]   (stringarr)  optional FITS header associated with
;                                   input image
;         [lsig=]     (float)      number of sigma below mean to set lower
;                                   display limit (default = 1.0)
;         [hsig=]     (float)      number of sigma above mean to set upper
;                                   display limit (default = 30.0)
;         [absmin=]   (float)      value of lower display limit -- if set then
;                                   overrides lsig
;         [absmax=]   (float)      value of upper display limit -- if set then
;                                   overrides hsig
;         [imgsz=]    (float)      image size in arcsec.  If chosen, will
;                                   produce an output image of imgsz x imgsz
;         [psfile=]   (string)     name for output postscript file (if not set,
;                                   no file is created)
;         [jpfile=]   (string)     name for output jpeg file (if not set,
;                                   no file is created)
;
; Revision history:
;  2003Mar31 Chris Fassnacht -- First working version.
;  2003Apr02 Chris Fassnacht -- Fixed several small bugs.
;

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, $
      'syntax: dispim, img [,inhead=inhead,lsig=lsig,hsig=hsig,absmin=absmin,'
    print, '  absmax=absmax,imgsz=imgsz,psfile=psfile,jpfile=jpfile]'
    print, ''
    print, ' Default value for lsig is 1.0'
    print, ' Default value for hsig is 30.0'
    print, ' *** If absmin is set, it overrides lsig. ***'
    print, ' *** If absmax is set, it overrides hsig. ***'
    print, ''
    return
endif

; Set directory names and defaults

if not keyword_set(lsig) then lsig=1.0
if not keyword_set(hsig) then hsig=30.0
pixscale = 0.05

; Get image information 

sz=size(img,/dimen)
print, ''
if keyword_set(inhead) then begin
   cd1 = fxpar(inhead, 'cd1_*', count=count1)
   cd2 = fxpar(inhead, 'cd2_*', count=count2)
   if (count1 eq 0) then begin
      pixscale=0.05 
   endif else begin
      pixscale=sqrt(cd1[0]*cd1[0] + cd1[1]*cd1[1])*3600.0
   endelse
   axtitle='arcsec'
endif else begin
   pixscale=1
   axtitle='pixels'
endelse
print, 'pixscale = ',pixscale
print, ''

; Set up scaled versions of the input file

imscal, img, scalimg, imgstat, nsig=3.0, lsig=lsig, hsig=hsig

; Check for absmin, absmax overrides

rescal = 0
if n_elements(absmin) gt 0 then begin
   llim = absmin
   rescal = 1
endif else llim = imgstat.min
if n_elements(absmax) gt 0 then begin
   hlim = absmax
   rescal = 1
endif else hlim = imgstat.max

if rescal gt 0 then scalimg = bytscl(img,min=llim,max=hlim)

; Set displayed image size

angxsz = sz[0]*pixscale/2.
angysz = sz[1]*pixscale/2.
print, ''
print, 'Native image size is ',angxsz*2.,axtitle,' x',angysz*2.,axtitle
if n_elements(imgsz) gt 0 then begin
   imgxsz = imgsz / 2.
   imgysz = imgsz / 2.
endif else begin
   imgxsz = angxsz
   imgysz = angysz
endelse
print, 'Displayed image size will be:',imgxsz*2.,axtitle,' x',imgysz*2.,axtitle
if keyword_set(inhead) then begin
   xr = [-1,1]*angxsz
   yr = [-1,1]*angysz
   ixr = [-1,1]*imgxsz
   iyr = [-1,1]*imgysz
endif else begin
   ixr = [0,sz[0]]
   iyr = [0,sz[1]]
   xr = [0,sz[0]]
   yr = [0,sz[1]]
endelse

; Set up viewport window

subcellarray, [1],[1], newpan, newsubpan

; Plot images

print, ''
print, 'Plotting image -- display limits are:'
print, '   lower = clipped mean - ',lsig,' sigma.'
print, '   upper = clipped mean + ',hsig,' sigma.'
print, ''

plotimage,scalimg,/preserve,_extra=_extra,panel=newpan[0,0,*], $
   subpan=newsubpan[0,0,*],imgxr=ixr,imgyr=iyr, $
   xrange=xr,yrange=yr, $
   xtit=axtitle,ytit=axtitle
;  ,title=fitsfile

; Plot to a postscript file

if n_elements(psfile) gt 0 then begin
   print,' Saving image in postscript file '+psfile+'.ps.'
   ps_open,psfile,/color,/ps_fonts
   device,xsize=6,ysize=6,/inches,/portrait,/times
   plotimage,scalimg,/preserve,_extra=_extra,panel=newpan[0,0,*], $
      subpan=newsubpan[0,0,*],imgxr=ixr,imgyr=iyr, $
      xrange=xr,yrange=yr, $
      xtit=axtitle,ytit=axtitle
;     ,title=fitsfile
   ps_close
endif

; Plot to a jpeg file (no border or labels)

if n_elements(jpfile) gt 0 then begin
   write_jpeg,jpfile,scalimg,true=1
   print,' Image saved in JPEG file '+jpfile+'.'
endif

end 
