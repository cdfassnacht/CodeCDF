Function mkgal, x0, y0, counts, halflight, axisratio=axisratio, posang=posang, $
   nx=nx, ny=ny, exponential=exponential, subpixfactor=subpixfactor, dmin=dmin, psf=psf
;
; Make the image of a galaxy centered at x0, y0 with given total 
; counts and half-light radius.  Make galaxy round unless 
; axis ratio is explicitly set to non-unity.  Position angle in degrees,
; counterclockwise from +y direction.  Use de Vaucouleurs law
; unless exponential is set.  Default size of the image is 128x128.
; Size in pixels; axis ratio must be <= 1.
; Allows convolution with PSF, which is assumed to be subpixellated.
;
if (keyword_set(nx) eq 0) then nx = 128
if (keyword_set(ny) eq 0) then ny = nx
;
circular = 1B
if (keyword_set(axisratio) ne 0) then begin
   circular = 0B
endif else begin
   circular = 1B
   axisratio = 1.
endelse
;
if (keyword_set(posang) eq 0) then posang = 0.
isrquarter = 1B
if (keyword_set(exponential) ne 0) then isrquarter = 0B
subpixellate = 0B
if (keyword_set(subpixfactor) ne 0) then begin
   subfactor=long(subpixfactor)
   subpixellate = 1B
   effnx = nx * subfactor
   effny = ny * subfactor
   ; Note: with origin 0, if the pixel position is interpreted as "centered",
   ; it transforms correctly with subpixellation.  However, "rebin" effectively
   ; shifts pixels by an amount equal to (subfactor-1)/(2.*subfactor).  Example:
   ; subfactor=10; original pixel 0 corresponds to x=0, original pixel 1 to x=1.0
   ; Subpixellated pixel 0 corresponds to x=0, subpixellated 1 to x=0.1, 10 to x=1.0.
   ; In rebinning, 0 is [0:9], which means effectively position 0.45 (=9/20).
   ; Thus the effective positions must be shifted by (subfactor-1.)/(2.*subfactor)
   ; in original coordinates, or (subfactor-1.)/2. in subpixel coordinates.
   shift = (subfactor-1.)/2.
   effx0 = x0 * subfactor + shift
   effy0 = y0 * subfactor + shift
endif else begin
   subfactor = 1L
   effnx = nx & effny = ny & effx0 = x0 & effy0 = y0
endelse
;
; Generate the (elliptical) distance mask
; Note: this part can be made *much* more efficient if the routine
; is called repeatedly.
;
nmax = effnx > effny
if (circular) then begin
   dist_circle, dd, nmax, effx0, effy0
endif else begin
   invratio = 1./axisratio
   dist_ellipse, dd, nmax, effx0, effy0, invratio,  posang
endelse
distmask = dd (0:effnx-1,0:effny-1) / float(subfactor)
image = fltarr (effnx,effny)
;
; Generate the counts
; Note: for the r^1/4 law, the intensity at the very center is extremely 
; high.  It is recommended to set a "minimum" distance of about 0.4 subpixels,
; so that the central pixel is not inordinately high.  (0.4 subpixels gives 
; approximately the correct answer if the pixel is perfectly centered and 
; the half-light radius is approx. 5 pixels.  The "correct" approach would 
; be to compute an approximate integral of the flux, possibly on a
; subpixellated basis.)
;
if (keyword_set(dmin) eq 0) then begin
   dmin = 0.
   if (isrquarter) then dmin = 0.4
endif
;
ellradius = distmask > dmin
;
; Now compute actual flux per (sub)pixel, scaled to unit total flux.
; Set a maximum (elliptical) distance to avoid "bad" calculations:
;     exponential: 10 exponential scalelengths
;     r^1/4:       100 half-light radii
;
; Flux is given by:
;     exponential: peak * exp (- ell_r / scalelength)
;     r^1/4:       hlb * exp(-7.67 ((ell_r/halflight)^(1/4)-1) )
;
; The total flux is:
;     exponential: pi * peak * scalelength^2 * axisratio
;     r^1/4:       22.64 * hlb * halflight^2 * axisratio
;
if (isrquarter) then begin
   hlb = counts / (22.64*axisratio*halflight^2)
   rmax = halflight * 100.
   acc = where (ellradius le rmax, nacc)
   if (nacc gt 0) then image (acc) = hlb * exp (-7.67*( (ellradius(acc)/halflight)^0.25 -1.))
endif else begin
   scalelength = halflight / 1.68
   rmax = scalelength * 10.
   peak = counts / (2.*3.14159 * scalelength^2 * axisratio)
   acc = where (ellradius le rmax, nacc)
   if (nacc gt 0) then image (acc) = peak * exp (- ellradius(acc) / scalelength)
endelse
; Convolve with PSF - edge-truncate sets points beyond image edge to zero
; PSF is assumed to be normalized to sum 1
if (keyword_set(psf) ne 0) then image = convol (image, psf, /edge_truncate)
;
; print,dummy
if (subpixellate) then image = rebin (image, nx, ny)
return, image
end
