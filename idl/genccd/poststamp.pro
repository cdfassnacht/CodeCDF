Function poststamp, fitsfile, centx, centy, nx, ny, outhead=outhead, $
                    display=display
;
; Function poststamp
;
; Description: Takes an input fitsfile and cuts out a smaller "postage
;  stamp" image of size nx by ny centered at (centx, centy).  Returns
;  the postage stamp image.
;
; Inputs: fitsfile   (string)      name of input fits file
;         centx      (float)       x position (in input image) to be
;                                   used as the center of the postage stamp
;         centy      (float)       y position (in input image) to be
;                                   used as the center of the postage stamp
;         nx         (int)         output size in x dimension
;         ny         (int)         output size in y dimension
;         [/display]               set this flag to display the extracted image
;
; Output accessed through passed parameters:
;         [outhead=] (strarray)    FITS header of output image
;
; Output: cutout     (floatarray)  cutout image
;
; Revision history:
;  2002Apr05 Chris Fassnacht -- First rough version.
;  2002Apr09 Chris Fassnacht -- Changed order of inputs to be more
;                               logical.
;  2002Apr12 Chris Fassnacht -- Added optional inhead and outhead
;                               keywords to allow the creation of a
;                               fits header for the postage stamp image. 
;  2002Apr18 Chris Fassnacht -- Error checking against non-integer values. 
;  2003May04 Chris Fassnacht -- Modified to take an input fits file and
;                                not an "image" (data part of a fits file).
;                                This is now made efficient by only reading
;                                in the appropriate subarray of the input
;                                file using the optional range parameters of
;                                fxread.
;  2006Jan23 Chris Fassnacht -- Added an optional display parameter to 
;                                display the cutout.
;

; Check input format

if n_params() lt 5 then begin
    print, ''
    print, 'syntax: outim = poststamp(fitsfile, centx, centy, nx, ny'
    print, '                         [, outhead=outhead,/display])'
    print, ''
    print, '  Set the display flag to display the cutout'
    nullim = fltarr(1)
    return, nullim
endif

; Create a blank image of the appropriate size to store the cutout.

cutout = fltarr(nx,ny)

; Get dimensions of input image

inhead=headfits(fitsfile)
sz=fxpar(inhead,'NAXIS*')
nx0 = sz[0]
ny0 = sz[1]
print, ''
print, 'poststamp: Input image dimensions: ', nx0, ny0

; Make sure that the postage stamp is smaller than the input image

if ((nx gt nx0) or (ny gt ny0)) then begin
    print, 'poststamp: One or more output dimensions are larger than '
    print, 'poststamp:  input dimensions.  Returning the input image.'
    cutout = mrdfits(fitsfile,0,inhead)
    return, cutout
endif

; Make everything an integer

cx = round(centx)
cy = round(centy)
inx = round(nx)
iny = round(ny)
halfx = round(inx/2)
halfy = round(iny/2)
print, 'poststamp: Centering postage stamp at ', $
    strtrim(cx,1),',',strtrim(cy,1)
print, 'poststamp: Requested postage stamp dimensions are ', $
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

; Read in selected portion of input image and copy into appropriate
;  range of output image
; Copy input image to output image

print, 'poststamp: Using ranges ',strtrim(x1,1),'-',strtrim(x2,1), $
    ' and ',strtrim(y1,1),'-',strtrim(y2,1)
fxread,fitsfile,imgtmp,hd,x1,x2,y1,y2
cutout[ox1:ox2,oy1:oy2] = imgtmp

; Create output header 

print, 'poststamp: Creating FITS header for postage stamp image.'
outhead = inhead
fxaddpar, outhead, 'naxis1', inx
fxaddpar, outhead, 'naxis2', iny
histstr = '  Postage stamp cutout centered at (' + strtrim(cx,1)
histstr = histstr + ',' + strtrim(cy,1) + ')'
fxaddpar, outhead, 'history', histstr

print, ''

; Display the image using dispim if display flag is set

if keyword_set(display) then dispim,cutout,inhead=outhead

return, cutout

end
