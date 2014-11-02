pro lrislog, fitsfile, inhead=inhead, lsig=lsig, hsig=hsig, absmin=absmin, $
   absmax=absmax, imgsz=imgsz,  $
   psfile=psfile, jpfile=jpfile, _extra=_extra

;
; Procedure lrislog
;
; Description: Given an input LRIS fits file, generates an output log
;  file containing important keyword values.
;
; Inputs: fitsfile      (string)     input image name
;
; Revision history:
;  2007Feb13 Chris Fassnacht -- First working version.
;

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, $
      'syntax: lrislog, fitsfile'
    print, ''
    return
endif

; Get head information of input image

inhead=headfits(fitsfile)

; Get image information - For each keyword, check to see if the
;  keyword was read in.  If not, replace it by a dummy value.

sz=size(img,/dimen)
frame = fxpar(inhead, 'frameno')
object = fxpar(inhead, 'object', count=count)
if (count eq 0) then object='---'
texp = fxpar(inhead, 'ttime', count=count)
if (count eq 0) then texp=-999.0
lst = fxpar(inhead, 'st', count=count)
if (count eq 0) then lst='---'
airmass = fxpar(inhead, 'airmass', count=count)
if (count eq 0) then airmass='---'
slitname = fxpar(inhead, 'slitname', count=count)
if (count eq 0) then slitname='---'
dichroic = fxpar(inhead, 'dichname', count=count)
if (count eq 0) then dichroic='---'
bgris = fxpar(inhead, 'grisname', count=count)
if (count eq 0) then bgris='---'
bfilt = fxpar(inhead, 'blufilt', count=count)
if (count eq 0) then bfilt='---'
rgrat = fxpar(inhead, 'graname', count=count)
if (count eq 0) then rfilt='---'
rfilt = fxpar(inhead, 'redfilt', count=count)
if (count eq 0) then rfilt='---'

print, frame, object, texp, lst, airmass, slitname, dichroic, bgris, bfilt, rgrat, rfilt

end 
