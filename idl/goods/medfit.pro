Pro medfit, image, fitim, modpar, nsig=nsig, inmask=inmask, $
            outmask=outmask, rkeep=rkeep, dodisplay=dodisplay
;
; Procedure medfit
;
; Description: Takes an input image and parameters output from the
;  fitell procedure and creates a "fit" image.  The "fit" is created
;  by doing a median smoothing of the counts as a function of
;  elliptical radius.  The parameters for computing the elliptical
;  radius are set by the modpar parameters, which have come from
;  fitell or SExtractor or something similar.  The program outputs the
;  median-smoothed image and, if requested, a mask of "bad" pixels
;  (i.e., those more than nsig-sigma displaced from the median at
;  each elliptical radius).
;
; Inputs: image      (floatarray)  input image
;         modpar     (floatarray)  vector containing fit parameters:
;                                   modpar[0] = background level
;                                   modpar[1] = x location of center
;                                   modpar[2] = y location of center
;                                   modpar[3] = total counts
;                                   modpar[4] = half-light radius
;                                   modpar[5] = axis ratio
;                                   modpar[6] = PA of ellipse
;         [nsig=]    (float)       number of sigma used in sigma
;                                   clipping
;         [inmask=]  (floatarray)  input bad pixel mask (bad pixels
;                                   set to 0)
;         [rkeep=]   (float)       elliptical radius within which
;                                   all mask values are set to 1
;                                   (i.e., good) even if they have
;                                   been set to 0 by sigclip.  This
;                                   allows the center of the galaxy
;                                   profile to be kept for future
;                                   fitting rounds.
;         [/dodisplay]             display radial profiles if set
;
; Output accessed through passed parameters:
;         fitim      (floatarray)  output median-smoothed image
;         [outmask=] (floatarray)  output bad pixel mask (bad pixels
;                                   set to 0)
;
; Revision history:
;  2002Apr15 Chris Fassnacht -- First rough version.
;  2002Apr16 Chris Fassnacht -- Added interpolation of median-filtered
;      values onto a regular half-pixel grid in elliptical radius.
;      Tried input mask, but made things worse in the center of the
;      galaxy, where too many pixels had been masked out.
;  2002Apr17 Chris Fassnacht -- Added the optional rkeep keyword
;      which, if set, restores all pixels in the mask interior to 
;      an elliptical radius of rkeep back to the "good" value, i.e. 1.
;      This may improve the situation in the center of the galaxy.
;      Added one iteration of masking and re-doing the median
;      filtering.
;  2002Apr23 Chris Fassnacht -- Changed the procedure which does the
;      median filtering.  Now instead of using IDL's canned "median"
;      routine, call the newly-written "filtmed.pro", which does
;      exactly the same thing as the IDL median routine, but also
;      allows a bad pixel mask to be passed and does the filtering
;      correctly with the mask.
;

; Check input format

if n_params() lt 3 then begin
    print, ''
    print, $
      'syntax: medfit, image, fitim, modpar [, nsig=, inmask=, outmask=, rkeep=]'
    print, ''
    return
endif

; Check dimensions of input image and set size of elliptical radius
;  vector equal to the larger of the two.

nx0 = (size(image,/dimen))[0]
ny0 = (size(image,/dimen))[1]
print, ''
print, 'medfit: Input image dimensions: ', nx0, ny0

ner = nx0 > ny0

; Generate the elliptical radius array using the input parameters.
; Note that the axis ratio required by dist_ellipse is (a/b) and not
;  the more usual (b/a).

print, 'medfit: Computing elliptical distance using:'
print, 'medfit:  (b/a) = ',modpar[5]
print, 'medfit:  PA =    ',modpar[6] 
dist_ellipse, er, ner, modpar[1], modpar[2], 1.0/modpar[5], modpar[6]

; Sort by elliptical radius

isort = sort(er)
xsrt = er(isort)
ysrt = image(isort)

; Mask the data if given an input mask

; ygood = ysrt
; if keyword_set(inmask) then begin
;     imasksrt = inmask(isort)
;     wgood = where(imasksrt gt 0,ngood)
;     if ngood gt 0 then ygood = ysrt(wgood)
; endif

; Do a median filtering on the data set, including the input mask,
;  if set.

msrt = 0.0 * ysrt
if keyword_set(inmask) then begin
    imasksrt = inmask(isort)
    filtmed, ysrt, msrt, 21, mask=imasksrt
endif else filtmed, ysrt, msrt, 21

; msrt = median(ygood,21)
mresid = ysrt - msrt

; Interpolate the median-filtered values and calculate the
;  residuals between the fit and the data.

; medfit_interp, xsrt, msrt, mint
; mresid = ysrt - mint

; Create the bad-pixel mask by doing a nsig-sigma clipping on
;  fit residuals.  If nsig has not been set, use the default
;  value of 3.0.
; NB: The RMS value used will be a local (in the radial sense)
;  RMS computed in regularly-spaced bins in elliptical radius
;  of size xbin.

if not keyword_set(nsig) then nsig = 3.0
sigmask = mresid
xbin = 3.0

radsigclip, xsrt, mresid, nsig, sigmask, xbin=xbin
if keyword_set(rkeep) then begin
    print, 'medfit: Setting central pixels in mask to 1 (r<', $
      strtrim(rkeep,1),')'
    wkeep = where(xsrt lt rkeep,nkeep)
    if nkeep gt 0 then sigmask[wkeep] = 1
endif

; Re-do the median-filtering without the newly-masked pixels

if keyword_set(inmask) then combmask = imasksrt < sigmask $
else combmask = sigmask
; wgood = where(combmask gt 0,ngood)
; if ngood gt 0 then ygood = ysrt(wgood)
; msrt = median(ygood,21)

filtmed, ysrt, msrt, 21, mask=combmask
mresid = ysrt - msrt

; ; Re-do the interpolation

; medfit_interp, xsrt, msrt, mint
; mresid = ysrt - mint

; Re-do the sigma clipping

radsigclip, xsrt, mresid, nsig, sigmask, xbin=xbin
if keyword_set(rkeep) then begin
    print, 'medfit: Setting central pixels in mask to 1 (r<', $
      strtrim(rkeep,1),')'
    wkeep = where(xsrt lt rkeep,nkeep)
    if nkeep gt 0 then sigmask[wkeep] = 1
endif

; Create the median-sorted image

print, ''
print, 'medfit: Creating final image.'
; mtmp = 0 * mint
; mtmp[isort] = mint
mtmp = 0 * msrt
mtmp[isort] = msrt
fitim = reform(mtmp,nx0,ny0)

if keyword_set(outmask) then begin
    msktmp = 0 * sigmask
    msktmp[isort] = sigmask
    outmask = reform(msktmp,nx0,ny0)
endif

end
