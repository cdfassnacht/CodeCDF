;+
;
; NAME: xybinsub
;
; PURPOSE: bin up oversampled images.
;
; USAGE:
;   binimage = xybinsub(image, xbinfact, ybinfact)
;
; ARGUMENTS:
;   image: an oversampled image of some sort.
;   xbinfact, ybinfact: the binning factors to apply in x and y.
;
; RETURNS: a binned image.
;
; WRITTEN: 2006may abolton@cfa
;
;-

function xybinsub, image, xbinfact, ybinfact

xbf = fix(abs(xbinfact)) > 1
ybf = fix(abs(ybinfact)) > 1
sz = size(image)
nxb = sz[1] / xbf
nyb = sz[2] / ybf
imstack = fltarr(nxb, nyb, xbf*ybf)
xind = (xbf * lindgen(nxb)) # replicate(1L, nyb)
yind = replicate(1L, nxb) # (ybf * lindgen(nyb))

for i = 0, xbf-1 do for j = 0, ybf-1 do imstack[*,*,xbf*j+i] = image[xind+i,yind+j]

oimage = total(imstack, 3)

return, oimage
end
