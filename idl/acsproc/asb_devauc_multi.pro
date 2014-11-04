;+
;
; NAME: asb_devauc_multi
;
; PURPOSE: Generate a 2d deVaucouleurs with
;   specified parameters.
;
; USAGE:
;   image = asb_devauc_multi(x, y, p [, basis=basis, rmin=rmin, $
;    psf=psf, logpars=logpars, rpow=rpow, subs=subs, _Extra=junk])
;
; INPUTS:
;   x and y give the 2d points at which to evaluate.
;   p is the parameter vector:
;     p[0] = Surface brightness at r=0 (but see ("rpow" option)
;     p[1] = intermediate-axis half-light radius
;     p[2] = x-center
;     p[3] = y-center
;     p[4] = axis ratio
;     p[5] = counterclockwise major-axis rotation w.r.t x-axis
;      another 6 (or 12 or 18, etc.) parameter-stacks may be
;       concatenated onto "p" to generate a multiple deVaucouleur image.
;   basis: keyword to set in order to return not a summend image
;    of components but rather a basis image stack.
;   rmin: minimum radius, implemented as a core radius.
;   psf: psf with which to convolve the image
;   rpow: multiply surface-brightness parameter by this power of
;    half-light radius to give the desired brightness parameter.
;    Why?  Because we may suppress parameter covariance this way.
;   logpars: set this keyword to indicate that the *natural log* of
;    the brightness and radius parameters is supplied, rather than
;    the parameter values themselves.
;   subs: integer factor by which the coordinate images subsample the
;    desired output image.  Output image will be binned by this factor
;    before psf convolution and output.
;    Obviously will only work if the coordinate images are contiguous.
;
;
; Written Jan 2005 by A.. Bolton, MIT.
; Minor changes and renaming: abolton@cfa 2006apr
; Added exponent variation and "basis" functionality,
;  renamed again: abolton@cfa 2006jun
; Made special devauc case 2006sep abolton@cfa.
;
;-

function asb_devauc_multi, x, y, p, basis=basis, rmin=rmin, psf=psf, $
 logpars=logpars, rpow=rpow, subs=subs, _Extra=junk
common tmp, zmod

deg2rad = !pi / 180.
k_dev = 7.66925001

ss_loc = keyword_set(subs) ? fix(subs) : 1

if (n_elements(x) ne n_elements(y)) then begin
  splog, ' x and y dimensions mis-matched!'
  return, 0
end

npar = 6L
ndev = n_elements(p) / npar
if (ndev lt 1) then begin
  splog, ' insufficient number of parameters!'
  return, 0
endif

rml = keyword_set(rmin) ? rmin : 0.

nx = (size(x))[1] / ss_loc
ny = (size(x))[2] / ss_loc

; Loop over all specified components:
zmod = fltarr(nx, ny, ndev)
for i = 0L, ndev-1 do begin
  ioff = i * npar
  a = keyword_set(logpars) ? exp(p[0+ioff]) : p[0+ioff]
  sig = keyword_set(logpars) ? exp(p[1+ioff]) : p[1+ioff]
  if keyword_set(rpow) then a = a / sig^rpow
  xcen = p[2+ioff]
  ycen = p[3+ioff]
  q = p[4+ioff]
  phi = p[5+ioff] * deg2rad
  xp = (x-xcen) * cos(phi) + (y-ycen) * sin(phi)
  yp = (y-ycen) * cos(phi) - (x-xcen) * sin(phi)
  r_ell = sqrt((xp^2)*q + (yp^2)/q + rml^2)/abs(sig)
;  if (min(r_ell) le 0.) then stop
  z_this = a * exp(-k_dev*r_ell^0.25)
  if (ss_loc gt 1) then z_this = xybinsub(z_this, ss_loc, ss_loc) / float(ss_loc^2)
  if keyword_set(psf) then z_this = convolve(z_this, psf)
  if keyword_set(basis) then zmod[*,*,i] = z_this $
   else zmod = zmod + z_this
endfor

return, zmod
end
