;+
;
; NAME: mpmoffatstack
;
; PUPOSE: stack of moffat images for MPFIT
;
; abolton@cfa 2007apr
;
;-

function mpmoffatstack, par, x=x, y=y, data=data, ivar=ivar, $
 deviates=deviates

    ni = (size(x))[1]
    nj = (size(x))[2]
    nf = n_elements(x) / (ni * nj)

    modim = 0. * x
    if (size(par))[0] gt 0 then begin
        moffatpar = reform(par, 7, nf)
        for i = 0, nf-1 do modim[*,*,i] = asb_moffat2d( $
         x[*,*,i], y[*,*,i], moffatpar[*,i])

        if keyword_set(deviates) then modim = ((data-modim) * sqrt(ivar))[*]
    endif else message, "MPMOFFATSTACK ERROR par = NAN"

    return, modim
end

