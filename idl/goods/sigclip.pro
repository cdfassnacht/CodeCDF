pro sigclip, image, outstat, nsig, ftol=ftol, mask=mask, verbose=verbose
;
; Procedure sigclip
;
; Description: Takes an input image does an iterative f*sigma clipping
;               until the RMS of the clipped image converges.  The
;               factor f is passed as the nsig parameter.
;
; Inputs: image (floatarray)      input image
;         outstat (strcut)        output structure containing image statistics
;         nsig (float)            level of sigma clipping.  e.g., for
;                                  3-sigma clipping, set nsig = 3.0
;         [ftol=] (float)         minimum fractional change in rms to 
;                                  continue clipping process.
;         [/verbose]              set this flag for verbose output
;
; Outputs accessed through passed parameters:
;         [mask=] (floatarray)    mask of bad pixels (bad pixels set
;                                  to 0)
;
; Revision history:
;  2002Apr05 Chris Fassnacht -- First rough version.
;  2002Apr09 Chris Fassnacht -- Made most of the inputs optional for
;      simplest calling format.  Improved error checking slightly.
;  2003Mar18 Chris Fassnacht -- Moved all of the formerly optional passed
;                                info variables such as immean, imrms, etc.
;                                into a new output structure called outstat
;

; Check input format

if n_params() lt 2 then begin
    print, ''
    print, $
      'syntax: sigclip, image, outstat, nsig [, ftol=, mask=, /verbose]'
    print, ''
    return
endif

; Set up output structure

outstat =  create_struct('ngood',0l,$
                         'mean',0.0d,$
                         'median',0.0d,$
                         'rms',0.0d,$
                         'mode',0.0d,$
                         'min',0.0d,$
                         'max',0.0d)
   
; Compute mean and rms of input image to use as starting values

imstat = moment(image,/double)
clipmean = imstat[0]
cliprms = sqrt(imstat[1])
lastrms = 1000.0 * cliprms

; Set float tolerance

if not keyword_set(ftol) then ftol = 0.1

; Loop on clipping

print, format='(a,f5.2,a)', 'sigclip: Clipping at ',nsig,' sigma.'
print, 'sigclip: ------------------------'

repeat begin
   lastrms = cliprms
   tmpgood = 0
   clipsum = 0.0d
   wgood = where (abs(image - clipmean) lt (nsig * cliprms), tmpgood)
   if tmpgood gt 0 then begin
      goodim = image(wgood)
      goodstat = moment(goodim,/double)
      clipmean = goodstat[0]
      cliprms = sqrt(goodstat[1])
      print, format='(a,f11.5,a,f11.5,a,i)', 'sigclip: Mean = ',clipmean, $
             ' RMS = ',cliprms,' ngood = ',tmpgood
   endif else begin
      print, 'ERROR: sigclip. No valid points in clipped array'
      tmpgood = -1
      return
   endelse

endrep until (abs(lastrms - cliprms)/lastrms lt ftol)

; Create mask of bad pixels, where bad pixels have value 0 and all others
;  have value 1.

if keyword_set(mask) then begin
    print, 'sigclip: Creating bad pixel mask.'
    mask = 0.0 * image + 1.0
    wbad = where (abs(image - clipmean) gt (nsig * cliprms), nbad)
    if nbad gt 0 then mask(wbad) = 0.
endif

; Transfer values into the output structure
outstat.ngood = tmpgood
outstat.mean = clipmean
outstat.rms = cliprms

if keyword_set(verbose) then help,/structure,outstat

end

