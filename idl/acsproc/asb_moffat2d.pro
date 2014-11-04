;+
;
; NAME: asb_moffat2d
;
; PURPOSE: Generate a 2d Moffat with
;   specified parameters.  a.k.a. a cored
;   powerlaw, I think...
;
; USAGE:
;   moffat_im = asb_moffat2d(x, y, p)
;
; Functional form adopted is:
;  I(r) = I_0 / (1. + (r / rscale)^2)^gamma
;
; INPUTS:
;   x and y give the 2d points at which to evaluate.
;   p is the parameter vector:
;     p[0] = scale (I_0)
;     p[1] = intermediate-axis rscale
;     p[2] = "gamma" (?) exponent of moffat profile
;     p[3] = x-center
;     p[4] = y-center
;     p[5] = axis ratio
;     p[6] = counterclockwise major-axis rotation w.r.t x-axis
;   another 7 (or 14 or 21, etc.) parameter stacks may be
;   concatenated onto "p" to generate a multiple-Moffat image.
;
; Written Feb 2005 by A.. Bolton, MIT, following g2dfunc
; Minor revisions and renaming, abolton@cfa 2006apr
; Adapted from asb_sersic2d, abolton@cfa 2006apr
;
;-

function asb_moffat2d, x, y, p, _Extra=junk

deg2rad = !pi / 180.

if (n_elements(x) ne n_elements(y)) then begin
  splog, ' x and y dimensions mis-matched!'
  return, 0
end

npar = 7L
nmoffat = n_elements(p) / npar
if (nmoffat lt 1) then begin
  splog, ' insufficient number of parameters!'
  return, 0
endif

; Loop over all specified Moffats:
zmod = 0. * x
for i = 0L, nmoffat-1 do begin
  ioff = i * npar
  a = p[0+ioff]
  rscale = p[1+ioff]
  gamma = p[2+ioff]
  xcen = p[3+ioff]
  ycen = p[4+ioff]
  q = p[5+ioff]
  phi = p[6+ioff] * deg2rad
  xp = (x-xcen) * cos(phi) + (y-ycen) * sin(phi)
  yp = (y-ycen) * cos(phi) - (x-xcen) * sin(phi)
  rr = sqrt((xp^2)*q + (yp^2)/q)
  zmod = zmod + a / (1. + (rr/rscale)^2)^gamma
endfor

return, zmod
end
