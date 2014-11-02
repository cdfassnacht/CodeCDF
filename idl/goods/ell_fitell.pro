Pro ell_fitell, x, apar, y
;
; Procedure to call mkgal to create a galaxy image using a set of
;  input paramters.  The galaxy is set to be an elliptical-type, i.e.,
;  r^(1/4) law.  This procedure is used in fitting to an input image.
;
; Revision history:
;   2002Apr04 Stefano Casertano -- First working version.
;   2002Apr08 Chris Fassnacht -- Slightly improved documentation.
;
common galimage, inputimage, pixelmask, nx, ny, subpix, psf, usepsf
;
; print,apar
background = apar(0)
x0 = apar(1)
y0 = apar(2)
counts = apar(3)
halflight = apar(4)
circular = 1B
if (n_elements(apar) gt 5) then begin
   circular = 0B
   axrat = apar(5)
   posang = apar(6)
endif else begin
   axrat = 1.
   posang = 0.
endelse
;
acc = where (pixelmask ne 0, nacc)
;
if (nacc le 0) then begin
; no valid point!!!  Should never happen
   y=0.
   return
endif
;
if (axrat gt 1 or axrat le 0 or halflight le 0. or counts lt 0.) then begin
   y = replicate (0.0, nacc)
endif else begin
;
   if (circular) then begin
      if (usepsf) then begin
         imfit = mkgal (x0, y0, counts, halflight, nx=nx, ny=ny, dmin=0.35, subpix=subpix, psf=psf)
      endif else begin
         imfit = mkgal (x0, y0, counts, halflight, nx=nx, ny=ny, dmin=0.35, subpix=subpix)
      endelse
   endif else begin
      if (usepsf) then begin
         imfit = mkgal (x0, y0, counts, halflight, axisratio=axrat, posang=posang, $
            nx=nx, ny=ny, dmin=0.35, subpix=subpix, psf=psf)
      endif else begin
         imfit = mkgal (x0, y0, counts, halflight, axisratio=axrat, posang=posang, $
            nx=nx, ny=ny, dmin=0.35, subpix=subpix)
      endelse
   endelse
   y = imfit (acc) + background
endelse
;
return
end
