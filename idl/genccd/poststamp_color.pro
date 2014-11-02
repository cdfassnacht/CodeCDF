pro poststamp_color, bfits, gfits, rfits, centx, centy, nx, ny, $
                     bwtfile=bwtfile, gwtfile=gwtfile, rwtfile=rwtfile, $
                     lsig=lsig, hsig=hsig, psfile=psfile, jpfile=jpfile
;
; Function poststamp_color
;
; Description: Takes 3 input registered fits files that contain data
; in 3 bands and displays a "truecolor" version of a smaller "postage
;  stamp" image of size nx by ny centered at (centx, centy).  
;
; Inputs:
;   bfits      (string)      name of shortest wavelength fits file
;   gfits      (string)      name of middle wavelength fits file
;   rfits      (string)      name of longest wavelength fits file
;   centx      (float)       x position (in input image) to be
;                             used as the center of the postage stamp
;   centy      (float)       y position (in input image) to be
;                             used as the center of the postage stamp
;   nx         (int)         output size in x dimension
;   ny         (int)         output size in y dimension
;   [bwtfile=] (string)      optional weight file for bfits
;   [gwtfile=] (string)      optional weight file for gfits
;   [rwtfile=] (string)      optional weight file for rfits
;   [lsig=]    (float)       number of sigma below mean to set lower
;                             display limit (default = 1.0)
;   [hsig=]    (float)       number of sigma above mean to set upper
;                             display limit (default = 15.0)
;
; Output:
;   (none)
;
; Revision history:
;  2006Feb02 Chris Fassnacht -- First rough version.
;  2007Aug09 Chris Fassnacht -- Added the optional input weight files.
;

; Check input format

if n_params() lt 5 then begin
    print, ''
    print, 'poststamp_color: Given 3 input fits files, an (x,y) location,'
    print, ' and an image size, produce and display a 3-color image.'
    print, ''
    print, 'syntax:  poststamp_color, bfits, gfits, rfits, centx, centy, nx, ny'
    print, '          [bwtwfits=bwtfits,gwtfits=gwtfits,rwtfits=rwtfits,'
    print, '          lsig=lsig,hsig=hsig,psfile=psfile,jpfile=jpfile])'
    print, ''
    print, 'Optional parameters:'
    print, '  bwtfits, etc.: Weight files associated with the input fits files'
    print, '  lsig:          Lower limit of display range, in units of'
    print, '                  sigma below the clipped mean.  Default=1.0'
    print, '  hsig:          Upper limit of display range, in units of'
    print, '                  sigma above the clipped mean.  Default=15.0'
    print, '  psfile:        Output postscript file name, if desired.  This'
    print, '                  procedure automatically adds the ".ps" at end'
    print, '  jpfile:        Output jpeg file name, if desired' 
    nullim = fltarr(1)
    return
endif

; Set defaults

if not keyword_set(lsig) then lsig=1.0
if not keyword_set(hsig) then hsig=15.0

; Create blank images of the appropriate sizes to store the cutouts.

bim = fltarr(nx,ny)
gim = fltarr(nx,ny)
rim = fltarr(nx,ny)
bwtim = fltarr(nx,ny)
gwtim = fltarr(nx,ny)
rwtim = fltarr(nx,ny)

; Get dimensions of input images

inhead=headfits(bfits)
sz=fxpar(inhead,'NAXIS*')
nx0 = sz[0]
ny0 = sz[1]
print, ''
print, 'poststamp_color: Input image dimensions: ', nx0, ny0

; Make sure that the postage stamp is smaller than the input image

if ((nx gt nx0) or (ny gt ny0)) then begin
    print, 'poststamp_color: One or more output dimensions are larger than '
    print, 'poststamp_color:  input dimensions.'
    return
endif

; Make everything an integer

cx = round(centx)
cy = round(centy)
inx = round(nx)
iny = round(ny)
halfx = round(inx/2)
halfy = round(iny/2)
print, 'poststamp_color: Centering postage stamp at ', $
    strtrim(cx,1),',',strtrim(cy,1)
print, 'poststamp_color: Requested postage stamp dimensions are ', $
    strtrim(inx,1),',',strtrim(iny,1)

; Set limits for copying, checking for edges of input image.

x10 = cx - halfx
x20 = x10 + inx - 1
y10 = cy - halfy
y20 = y10 + iny - 1
x1 = cx - halfx > 0
x2 = x10 + inx - 1 < (nx0 - 1)
y1 = cy - halfy > 0
y2 = y10 + iny - 1 < (ny0 - 1)
xsize = x2 - x1 + 1
ysize = y2 - y1 + 1
if x10 lt 0 then ox1 = inx - xsize else ox1 = 0
if y10 lt 0 then oy1 = iny - ysize else oy1 = 0
if x20 gt (nx0-1) then ox2 = xsize - 1 else ox2 = inx - 1
if y20 gt (ny0-1) then oy2 = ysize - 1 else oy2 = iny - 1

; Read in selected portion of input images and copy into appropriate
;  range of output image
; Copy input image to output image

print, 'poststamp_color: Using ranges ',strtrim(x1,1),'-',strtrim(x2,1), $
    ' and ',strtrim(y1,1),'-',strtrim(y2,1)

fxread,bfits,imgtmp,hd,x1,x2,y1,y2
bim[ox1:ox2,oy1:oy2] = imgtmp
fxread,gfits,imgtmp,hd,x1,x2,y1,y2
gim[ox1:ox2,oy1:oy2] = imgtmp
fxread,rfits,imgtmp,hd,x1,x2,y1,y2
rim[ox1:ox2,oy1:oy2] = imgtmp

; Read in weight files, if such files have been designated.

if n_elements(bwtfile) gt 0 then begin
   fxread,bwtfile,imgtmp,hd,x1,x2,y1,y2
   bwtim[ox1:ox2,oy1:oy2] = imgtmp
endif
if n_elements(gwtfile) gt 0 then begin
   fxread,gwtfile,imgtmp,hd,x1,x2,y1,y2
   gwtim[ox1:ox2,oy1:oy2] = imgtmp
endif
if n_elements(rwtfile) gt 0 then begin
   fxread,rwtfile,imgtmp,hd,x1,x2,y1,y2
   rwtim[ox1:ox2,oy1:oy2] = imgtmp
endif

; Get information from first file's header

cd1 = fxpar(hd, 'cd1_*', count=count1)
cd2 = fxpar(hd, 'cd2_*', count=count2)

if (count1 lt 2) then begin
   pixscale=0.05 
endif else begin
    pixscale=sqrt(cd1[0]*cd1[0] + cd1[1]*cd1[1])*3600.0
endelse
print, 'pixscale = ',pixscale

; Get image scaling

hssz = size(hsig,/dimen)
if(hssz ne 3) then imhs=[hsig, hsig, hsig] else imhs = hsig

; Set up image masks based on weight files, if selected, and then
;  create scaled versions of the input files, either using the masks
;  (if created) or without any masks

if n_elements(bwtfile) gt 0 then begin
   bmask = intarr(nx,ny)
   bwbad = where(bwtim eq 0, nbad)
   if nbad gt 0 then bmask(bwbad) = 1
   imscal, bim, bimsc, bstat, hsig=imhs[0], nsig=3.0, goodmask=bmask
endif else begin
   imscal, bim, bimsc, bstat, hsig=imhs[0], nsig=3.0
endelse

if n_elements(gwtfile) gt 0 then begin
   gmask = intarr(nx,ny)
   gwbad = where(gwtim eq 0, nbad)
   if nbad gt 0 then gmask(gwbad) = 1
   imscal, gim, gimsc, gstat, hsig=imhs[1], nsig=3.0, goodmask=gmask
endif else begin
   imscal, gim, gimsc, gstat, hsig=imhs[1], nsig=3.0
endelse

if n_elements(rwtfile) gt 0 then begin
   rmask = intarr(nx,ny)
   rwbad = where(rwtim eq 0, nbad)
   if nbad gt 0 then rmask(rwbad) = 1
   imscal, rim, rimsc, rstat, hsig=imhs[2], nsig=3.0, goodmask=rmask
endif else begin
   imscal, rim, rimsc, rstat, hsig=imhs[2], nsig=3.0
endelse

; Set up 3-plane array for inputs to true-color images

rgb=bytarr(3,nx,ny)
rgb[0,*,*] = rimsc
rgb[1,*,*] = gimsc
rgb[2,*,*] = bimsc

; Set displayed image size

angxsz = nx*pixscale/2.
angysz = ny*pixscale/2.
print, ''
print, 'Native image size is ',angxsz*2.,' arcsec x',angysz*2.,' arcsec'
if n_elements(imgsz) gt 0 then begin
   imgxsz = imgsz / 2.
   imgysz = imgsz / 2.
endif else begin
   imgxsz = angxsz
   imgysz = angysz
endelse
print, 'Displayed image size will be:',imgxsz*2.,' arcsec x',imgysz*2.,' arcsec'
print, 'Upper value of display ranges are:'
print, '   ',bfits,': ',imhs[0],' sigma above the clipped mean'
print, '   ',gfits,': ',imhs[1],' sigma above the clipped mean'
print, '   ',rfits,': ',imhs[2],' sigma above the clipped mean'

; Set up viewport window

subcellarray, [1],[1], newpan, newsubpan

; Plot images

print, ''
print, 'Plotting 3-color image'
print, ''

plotimage,rgb,/preserve,_extra=_extra,panel=newpan[0,0,*], $
   subpan=newsubpan[0,0,*],imgxr=[-1,1]*angxsz,imgyr=[-1,1]*angysz, $
   xrange=[-1,1]*imgxsz,yrange=[-1,1]*imgysz, $
   xtit='arcsec',ytit='arcsec',title=bigroot
;xyouts,-0.9*sz[0]*pixscale/2.,0.8*sz[1]*pixscale/2.,'BVI',/data, $
;   charsize=3.0, font=1

; Plot to a postscript file

if n_elements(psfile) gt 0 then begin
   print,' Saving image in postscript file '+psfile+'.ps.'
   ps_open,psfile,/color,/ps_fonts
   device,xsize=6,ysize=6,/inches,/portrait,/times
   plotimage,rgb,/preserve,_extra=_extra,panel=newpan[0,0,*], $
      subpan=newsubpan[0,0,*],imgxr=[-1,1]*angxsz,imgyr=[-1,1]*angysz, $
      xrange=[-1,1]*imgxsz,yrange=[-1,1]*imgysz, $
      xtit='arcsec',ytit='arcsec',title=bigroot
;   xyouts,-0.9*sz[0]*pixscale/2.,0.8*sz[1]*pixscale/2.,'BVI',/data, $
;      charsize=3.0, font=0, color=255
   ps_close
endif

; Plot to a jpeg file (no border or labels)

if n_elements(jpfile) gt 0 then begin
   write_jpeg,jpfile,rgb,true=1
   print,' Image saved in JPEG file '+jpfile+'.'
endif

end
