Pro medfit_interp, xin, yin, yout
;
; Procedure medfit_interp
;
; Description: Takes the input x and y vectors and interpolates
;               onto a regular grid in x by averaging all points
;               in annuli in x.  The procedure then creates the output
;               y vector by replacing each input y value by its
;               interpolated value.
;    *** NB: This procedure assumes that x has been sorted. ***
;
; Inputs: xin        (floatarray)  input irregularly spaced x vector
;         yin        (floatarray)  input y values
;
; Output accessed through passed parameters:
;         yout       (floatarray)  output interpolated y values
;
; Revision history:
;  2002Apr16 Chris Fassnacht -- First working version.
;

; Check input format

if n_params() lt 3 then begin
    print, ''
    print, 'syntax: medfit_interp, xin, yin, yout
    print, ''
    return
endif

; Find the limits for the x interpolation

xmax = max(xin,min=xmin)
xinterpmin = floor(xmin)
xinterpmax = ceil(xmax)
ninterp = 2 * (1L + (xinterpmax - xinterpmin))

; Interpolate onto half-pixel steps

xinterp = findgen(ninterp)/2.0 + xinterpmin

; Set up vectors to be filled

yinterp = fltarr(ninterp - 1)
yout = 0.0 * xin
nbint = 0.0 * yin

; Calculate the interpolated y value at each xinterp by calculating
;  the mean of all the input y points that fall inside the x binwidth.

for i=0L,ninterp-2 do begin
    wbin = where (xin ge xinterp[i] and xin lt xinterp[i+1],nbin)
    case nbin of
        0: yinterp[i] = 0.0
        1: yinterp[i] = yin[wbin]
        else: yinterp[i] = mean(yin[wbin])
    endcase
    if nbin gt 0 then yout[wbin] = yinterp[i]
    nbint[i] = nbin
endfor
; xinterp2 = findgen(ninterp-1)/2.0 + xinterpmin + 0.25

end
