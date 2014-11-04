;+
;
; NAME: fft_register
;
; PURPOSE: register pixel-shifted images by maximizing
;   FFT cross correlation.  Limits to a specific window range.
;   Will only work for the special case of no rotations
;   and no change in point-to-point differential astrometry.
;
; USAGE:
;   fft_register, ims, shwin, xshift, yshift [, nofrac=nofrac]
;
; ARGUMENTS:
;   ims: stack of nx X ny X nf images
;   shwin: maximum shift to consider.  Also defines masking window.
;   nofrac: set this keyword to suppress fractional shifting.
;
; RETURNS:
;   xshift, yshift:  nf X nf arrays of relative x and y
;     pixel shifts.  the [i,j]th element of each is the
;     active shift to apply to image j to bring it into
;     the frame of image i.  i/j asymmetry is due to the
;     mask that is applied to frame i but not j under the
;     FFT.  This also explains the non-zero diagonal entries.
;
; WRITTEN abolton@cfa 2007apr
;
;-

pro fft_register, ims, shwin, xshift, yshift, nofrac=nofrac

nx = (size(ims))[1]
ny = (size(ims))[2]
nf = (size(ims))[3]

; The mask:
cmask = 0B * bytarr(nx, ny)
cmask[shwin:nx-shwin-1,shwin:ny-shwin-1] = 1B
;cmask[*] = 1B

; The FFTs:
imft = dcomplexarr(nx, ny, nf)
mkft = dcomplexarr(nx, ny, nf)
for i = 0, nf-1 do imft[*,*,i] = fft(ims[*,*,i], /double)
for i = 0, nf-1 do mkft[*,*,i] = fft(cmask * ims[*,*,i], /double)

; Loop over images in i and j
xshift = fltarr(nf,nf)
yshift = fltarr(nf,nf)
shbase = findgen(2*shwin+1) - shwin
hw = 1
for i = 0L, nf-1 do begin
  for j = 0L, nf-1 do begin
    ctest = real_part(fft(mkft[*,*,i] * conj(imft[*,*,j]), /inverse))
    ctest = (shift(ctest, shwin, shwin))[0:2*shwin,0:2*shwin]
    maxval = max(ctest, wmax)
    imax = wmax mod (2*shwin+1)
    jmax = wmax / (2*shwin+1)
    csub = ctest[imax-hw:imax+hw,jmax-hw:jmax+hw]
    peakfrac = keyword_set(nofrac) ? 0. : sfit_fpeak(csub) - float(hw)
    peakwhole = [shbase[imax], shbase[jmax]]
    tshift = peakwhole + peakfrac
    xshift[i,j] = tshift[0]
    yshift[i,j] = tshift[1]
  endfor
endfor

return
end
