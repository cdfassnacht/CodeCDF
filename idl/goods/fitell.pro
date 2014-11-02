Function fitell, image, parameters, errors, noise=noise, keeppeak=keeppeak, $
   mask=mask, exponential=exponential, subpixfactor=subpixfactor, psf=psf, $
   circular=circular, silent=silent, tol=tol
;
; Accept an image containing an elliptical galaxy; find the best-fitting
; parameters; return parameters, errors
;
; Accept noise image (assumes noise = 1. otherwise)
; Accept mask image (=0 where image is to be excluded)
; Blanks the highest pixel unless keeppeak is set
; Returns best-fit image
;
; Parameters:
; IMAGE    (input)      The pixel values; assumes background has been subtracted
; PARAMETERS (output)   fitted background, x0, y0, counts, halflight, 
;                        axisratio, positionangle
; ERRORS (output)       estimated errors in fitted parameters
; NOISE  (input)        RMS noise image, same size as input image
; KEEPPEAK (input)      Set to keep the peak pixel; default is to
;                        exclude the peak from the fit
; MASK (input)          Mask image, same size as input; set to 0 for
;                        pixels to be excluded
; EXPONENTIAL (input)   Set to change the fitted function from an
;                        elliptical-galaxy-like fit (r^(1/4) law) to
;                        an exponential disk fit.
; SUBPIXFACTOR (input)  Subpixellization factor to be used internally; 
;                        default 5, or 1 if PSF is specified
; PSF (input)           PSF for convolution; assumed normalized; 
;                        must be given in subpixels
; CIRCULAR (input)      Set to limit the fit to axially symmetric
;                        galaxy; PARAMETERS (5:6) and ERRORS(5:6)
;                        will be set to [1,0] and [0,0] on output
; SILENT (input)        Suppresses explicit printing of results
;
; Revision history:
;   2002Apr04 Stefano Casertano -- First working version.
;   2002Apr08 Chris Fassnacht --  Modified to also fit an exponential disk if
;       requested, through the use of the optional "exponential" keyword.  
;       Default fit is still an elliptical (r^(1/4)) fit.
;
;

; Information passed to subroutine via common:
; Size of image (nx, ny)
; Mask image (pixelmask)
; Subpixellation factor (subpix)
; Use PSF and its value (psf, usepsf)
;
common galimage, inputimage, pixelmask, nx, ny, subpix, psf0, usepsf
inputimage = image
imagesize = size(inputimage)
nx = imagesize(1)
ny = imagesize(2)
;
; Convolution requested?  If not, subpix defaults to 5
if (keyword_set(psf)) then begin
   usepsf = 1B
   psf0 = psf
   if (keyword_set(subpixfactor) eq 0) then subpixfactor = 1L
endif else begin
   usepsf = 0B
   if (keyword_set(subpixfactor) eq 0) then subpixfactor = 5L
endelse
subpix = subpixfactor
;
; Determine initial guesses for xcen, ycen, counts, halflight
; Use IDL function GAUSS2DFIT, restricted to fit along x and y
gaussianfit = GAUSS2DFIT (image, par0)
; par0 is the result vector:
; [0] constant term
; [1] scale factor (peak value)
; [2] x width (sigma)
; [3] y width (sigma)
; [4] x center
; [5] y center
; [6] rotation (defaults to zero)
; counts = 2 pi sigmax sigmay * peak
if (keyword_set(silent) eq 0) then print, 'Initial parameters ', par0
;
if (keyword_set(circular) eq 0) then circular=0B else circular=1B
width = sqrt (par0(2)*par0(3))
counts = par0(1) * 2*3.14159*width^2
;
if (circular) then begin
   parameters = fltarr(5)
   parameters(0) = par0(0) & parameters(1) = par0(4) & parameters(2) = par0(5) 
   parameters(3) = counts & parameters(4) = width
endif else begin
   parameters = fltarr(7)
   parameters(0) = par0(0) & parameters(1) = par0(4) & parameters(2) = par0(5)
   parameters(3) = counts & parameters(4) = width & parameters(5) = 0.8
endelse
;
; Now prepare the actual vector with weights
if (keyword_set(noise) eq 0) then noise=1.
weight = image*0. + 1. / noise^2
; which will make a vector of the right form whether noise is a scalar or a vector
;
; set to zero the weight of masked pixels
pixelmask = weight gt 0
if (keyword_set(mask)) then pixelmask = mask and pixelmask
;
; set to zero the weight of the brightest pixel
if (keyword_set(keeppeak) eq 0) then begin
   ispeak = max (image, max_subscript)
   pixelmask(max_subscript) = 0B
endif
w0 = where (pixelmask eq 0, nzero)
if (nzero gt 0) then weight(w0) = 0.
;
; Now make the 1-d vector curvefit wants...
acc = where (weight gt 0., nacc)
if (nacc le 0) then begin
   print, 'No valid points, exiting'
   fitted = image * 0.
   return, fitted
endif
;
y = image(acc)
w = weight (acc)
x = y*0.
; x is a dummy vector here needed by curvefit
;
if (keyword_set(tol) eq 0) then tol=1.e-3
;
; Choose the function used by curvefit depending on the exponential
;  keyword.
;
if (keyword_set(exponential)) then begin
   yfit = curvefit (x,y,w,parameters,errors,func='exp_fitell', $
                    /noder,chisq=chisq,tol=tol)
endif else begin
   yfit = curvefit (x,y,w,parameters,errors,func='ell_fitell', $
                    /noder,chisq=chisq,tol=tol)
endelse
;
; Output results, unless silent is set
if (keyword_set(silent) eq 0) then begin
   print, format='(a,2f10.4,a,2f10.4)','center x, y = ', parameters(1:2), ' +/- ', errors(1:2)
   print, ' background   = ', parameters(0), ' +/- ', errors(0), ' counts/pixel'
   print, ' total signal = ', parameters(3), ' +/- ', errors(3), ' counts to infinity'
   print, ' half-light   = ', parameters(4), ' +/- ', errors(4), ' pixels'
   if (circular eq 0B) then begin
      print, ' axis ratio   = ', parameters(5), ' +/- ', errors(5)
      print, ' pos angle    = ', parameters(6), ' +/- ', errors(6), ' deg'
   endif
endif
;
if (circular) then begin
   parameters = [parameters, 1., 0.]
   errors = [errors, 0., 0.]
endif
;
if (usepsf) then begin
   fitted = mkgal (parameters(1), parameters(2), parameters(3), parameters(4), $
      axisratio=parameters(5), posang=parameters(6), nx=nx, ny=ny, $
      subpixfactor=subpixpfactor, dmin=0.35, psf=psf0) + parameters(0)
endif else begin
   fitted = mkgal (parameters(1), parameters(2), parameters(3), parameters(4), $
      axisratio=parameters(5), posang=parameters(6), nx=nx, ny=ny, $
      subpixfactor=subpixpfactor, dmin=0.35) + parameters(0)
endelse
;
return, fitted
end
