pro run_fitell, image, inmask, inrms, impar, doellip, fitstring, $
   dodisplay=dodisplay

;
; Procedure run_fitell
;
; Description: A procedure to call fitell twice for an input image.
;               The first time the call uses the r^(1/4)-law fit and
;               the second uses an exponential disk fit.
;
; Inputs: image      (floatarray)  input image
;         inmask     (floatarray)  mask showing good pixels (1) and
;                                   bad pixels (0)
;         inrms      (float)       RMS noise in input image
;         [/dodisplay]             display images if set
;
; Output accessed through passed parameters:
;         impar      (floatarray)  parameters of fit
;         doellip    (byte)        flag set to 1 if r^(1/4) fit is
;                                   better
;         fitstring  (string)      string indicating which fit is better
;
; Revision history:
;  2002Aug06 Chris Fassnacht -- Moved from gopost.pro
;
;

; Check input format

if n_params() lt 6 then begin
    print, ''
    print, 'syntax: run_fitell, image, inmask, inrms, impar, doellip,'
    print, 'fitstring, [, /dodisplay]'
    print, ''
    return
endif

; Get indices of good pixels

wmask = where(inmask gt 0)

; Fit both an elliptical (r^1/4 law) and exponential disk model to the
;  input image.

print, ''
print, 'run_fitell: Fitting initial galaxy models (takes a while...)'
print, 'run_fitell: -------------------------------------------------------'
print, ''
print, 'ELLIPTICAL GALAXY FIT:'
fit_e = fitell(image,par_e,err_e,noise=inrms,mask=inmask)
resid_e = image - fit_e
goodresid = resid_e[wmask]
rms_e =  stddev(goodresid,/double)
print, 'RMS in residual image: ',rms_e
print, ''
print, 'EXPONENTIAL DISK GALAXY FIT:'
fit_s = fitell(image,par_s,err_s,noise=inrms,mask=inmask,/exponential)
resid_s = image - fit_s
goodresid = resid_s[wmask]
rms_s =  stddev(goodresid,/double)
print, 'RMS in residual image: ',rms_s
print, ''

; Choose best-fitting model as the one which produces the lowest rms
;  in the residual image.

if rms_e lt rms_s then begin
    doellip = 1B
    impar = par_e
    print, 'run_fitell: Elliptical galaxy fit is better.'
    fitstring = 'ELLIPTICAL GALAXY FIT:'
endif else begin
    doellip = 0B
    impar = par_s
    print, 'run_fitell: Exponential disk galaxy fit is better.'
    fitstring = 'EXPONENTIAL DISK GALAXY FIT:'
endelse

; Display residual images if requested

print, ''
if keyword_set(dodisplay) then begin
    print, '*** Displaying residual image (elliptical fit). ***'
    display, resid_e, dmin=-10*inrms, dmax=20*inrms
    print, '*** Displaying residual image (exponential fit). ***'
    display, resid_s, dmin=-10*inrms, dmax=20*inrms
    print, ''
endif

end
