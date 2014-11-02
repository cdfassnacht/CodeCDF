pro imcrop, image, goodmask, cropim, cropmask, xmin, xmax, ymin, ymax

;
; Procedure imcrop
;
; Description: Given an input image and a mask of good values, crops
;  the image to the minimum rectangular area that contains all the
;  good pixels.
;
; Inputs: image      (floatarray)  input image
;         goodmask   (floatarray)  mask of good values (set to 1 for
;                                   good pixels and 0 for bad pixels).
;
; Output accessed through passed parameters:
;         cropim     (floatarray)  cropped image
;         cropmask   (floatarray)  cropped version of mask
;         xmin       (float)       minimum x value
;         xmax       (float)       maximum x value
;         ymin       (float)       minimum x value
;         ymax       (float)       maximum x value
;
; Revision history:
;  2002Aug06 Chris Fassnacht -- Moved from gopost.pro
;

; Check input format

if n_params() lt 4 then begin
    print, ''
    print, 'syntax: imcrop, image, goodmask, cropim, cropmask, xmin, xmax'
    print, ' ymin, ymax'
    print, ''
    return
endif

; Make sure that there are good pixels

wgood = where(goodmask gt 0, ngood)
if ngood eq 0 then begin
    print, ''
    print, '*** ERROR: imcrop.  No good pixels in area mask! ***'
    print, '           Cannot continue!'
    return
endif

; Find the size of the image.

nx = (size(image,/dimen))[0]
ny = (size(image,/dimen))[1]

; Crop the image to its minimum size in x and y.  
; The x and y arrays are set up through the matrix multiplication operator (#)

x = findgen(nx) # replicate(1.0,ny)
y = replicate(1.0,nx) # findgen(ny)
goodx = x[wgood]
goody = y[wgood]
xmax = max(goodx, min=xmin)
ymax = max(goody, min=ymin)
print, 'imcrop: Cropping input image to range [',strtrim(xmin,1),':', $
  strtrim(xmax,1),',',strtrim(ymin,1),':',strtrim(ymax,1),']'
cropim = image[xmin:xmax,ymin:ymax]
cropmask = goodmask[xmin:xmax,ymin:ymax]

end
