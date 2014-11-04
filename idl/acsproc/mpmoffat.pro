;+
;
; NAME: mpmoffat
;
; PUPOSE: asb_moffat2d wrapper for MPFIT
;
; abolton@cfa 2007apr
;
;-

function mpmoffat, par, x=x, y=y, data=data, ivar=ivar, $
 deviates=deviates, psf=psf

modim = asb_moffat2d(x, y, par)

if keyword_set(psf) then modim = convolve(modim, psf)

if keyword_set(deviates) then modim = ((data-modim) * sqrt(ivar))[*]

return, modim
end

