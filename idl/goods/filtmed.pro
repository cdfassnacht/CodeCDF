Pro filtmed, yin, yout, width, mask=mask
;
; Procedure filtmed
;
; Description: Does a median filtering on yin by running a window
;  of width=width along yin and replacing each point by the median
;  value in the box centered on that point.  The first and last
;  (width/2) values in yin are not altered.  This is exactly the
;  same as IDL's median command used with a width, with one exception --
;  the possibility of including a mask.  If a mask is passed to the
;  function, then bad pixels (indicated by the places where the mask
;  is set to 0) are not included in the median.
;
; Inputs: yin        (floatarray)  input vector
;         width      (floatarray)  width of smoothing window, in array
;                                   entries (see help for IDL median)
;         [mask=]    (floatarray)  input bad pixel mask (bad pixels
;                                   set to 0)
;
; Output accessed through passed parameters:
;         yout       (floatarray)  output median-smoothed vector
;                                   set to 0)
;
; Revision history:
;  2002Apr23 Chris Fassnacht -- First working version.
;

; Check input format

if n_params() lt 3 then begin
    print, ''
    print, 'syntax: filtmed, yin, yout, width [, mask=mask]'
    print, ''
    return
endif

; Check number of elements in array

ny = n_elements(yin)
if ny le width then begin
    print, ''
    print, 'ERROR: filtmed.  Input array too small compared to width.'
    yout = yin
    return
endif

; Initialize yout to be the same as yin

yout = yin

; Start running the window

halfwidth = long(width / 2)

for i=halfwidth,ny-halfwidth-1 do begin
    tmpy = yin[i-halfwidth:i+halfwidth]
    if keyword_set(mask) then begin
        tmpmask = mask[i-halfwidth:i+halfwidth]
        wgood = where(tmpmask gt 0.0, ngood)
        if ngood gt 0 then goody = tmpy[wgood]
    endif else begin
        goody = tmpy
    endelse
    yout[i] = median(goody)
endfor

end
