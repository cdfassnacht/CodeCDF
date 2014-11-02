pro imscal, image, out, stat, lsig=lsig, hsig=hsig, nsig=nsig, $
  goodmask=goodmask

;
; Function imscal
;
; Description: Runs sigma clipping on the input image, sets the pixel
;               range to be mean-lsig*rms to mean+hsig*rms, where the
;               mean and rms are the values computed from the sigma-clipped
;               image and not the input image.
;
; Inputs: image       (floatarr)   input image
;         out         (bytearr)    output scaled image
;         stat        (struct)     output structure containing image info
;         [lsig=]     (float)      number of sigma below mean to set lower
;                                   display limit (default = 1.0)
;         [hsig=]     (float)      number of sigma above mean to set upper
;                                   display limit (default = 15.0)
;         [nsig=]     (float)      level of sigma clipping.  e.g., for
;                                   3-sigma clipping, set nsig = 3.0
;         [goodmask=] (intarr)     mask indicating which regions to
;                                   use in the sigma clipping.  Bad pixels
;                                   set to 1.
;
; Revision history:
;  2003Mar18 Chris Fassnacht -- Transferred functionality from truecolor.pro
;                               First working version.
;  2004Sep05 Chris Fassnacht -- Added a goodmask passed parameter to
;   restrict the area of the image used for sigma-clipping if necessary.
;  2007Aug13 Chris Fassnacht -- Fixed a bug that occurs if there is an
;                               input mask but there are no good
;                               pixels in the mask.

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, $
      'syntax: imscal, image, out [,lsig=lsig,hsig=hsig,nsig=nsig,stat=stat]'
    print, ''
    print, ' Default value for lsig is 1.0'
    print, ' Default value for hsig is 15.0'
    print, ' Default value for nsig is 3.0'
    print, ''
    return
endif

; Set values to defaults if not passed on command line

if n_elements(lsig) eq 0 then lsig = 1.0
if n_elements(hsig) eq 0 then hsig = 15.0
if n_elements(nsig) eq 0 then nsig = 3.0

; Use mask if passed as a parameter

if n_elements(goodmask) gt 0 then begin
   wgood = where(goodmask eq 0, ngood)
   if ngood gt 0 then begin
       tmpimg = image(wgood)
       print,'imscal: Using mask'
   endif else begin
       print, ''
       print, 'imscal: WARNING. No good pixels in input mask for requested area'
       print, 'imscal: Setting output scaled image to 0'
       insz = size(image,/dimension)
       out = bytarr(insz[0],insz[1])
       return
   endelse
endif else begin
   tmpimg = image
endelse

; Run sigclip to get image properties

sigclip, tmpimg, outstat, nsig

; Set up scaled versions of the input files

stat = outstat
stat.min = outstat.mean - lsig*outstat.rms
stat.max = outstat.mean + hsig*outstat.rms

out = bytscl(image,min=stat.min,max=stat.max)


end 
