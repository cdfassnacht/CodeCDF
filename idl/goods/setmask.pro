pro setmask, nx, ny, secat, maskfac, goodmask, retval

;
; Procedure imcrop
;
; Description: Given a set of shape parameters, creates a mask of "good"
;  values, i.e, set goodmask=1 for pixels inside the ellipse defined by
;  the shape parameters and goodmask=0 outside the ellipse.
;
; Inputs: nx         (float)       x dimension of mask image
;         ny         (float)       y dimension of mask image
;         secat      (structure)   Structure containing ellipse shape parameters
;                                   "size"
;                                   ellipticity (1 - b/a)
;                                   PA (N->E in astronomical sense)
;         maskfac    (float)       multipicative factor to define the 
;                                   semimajor axis of the good ellipse, i.e.,
;                                   a_good = a_input * maskfac
;
; Output accessed through passed parameters:
;         goodmask   (floatarray)  mask of good values (set to 1 for
;                                   good pixels and 0 for bad pixels).
;         retval     (float)       integer value signifying error (retval=1)
;                                   or success (retval=0)
;
; Revision history:
;  2003Feb25 Chris Fassnacht -- Moved from gopost.pro
;  2003Apr03 Chris Fassnacht -- Took out correction for PA
;

; Set default return value to error setting (retval = 1).  Only
;  set to success value (retval = 0) if procedure reaches the end without
;  error.

retval = 1

; Check input format

if n_params() lt 4 then begin
    print, ''
    print, 'syntax: setmask, image, secat, maskfac, goodmask'
    print, ''
    return
endif

; Set the mask pixels to bad where
;  their elliptical distance is more than maskfac times the
;  semimajor axis.

print, ''
print, 'setmask: Input shape parameters:'
print, 'setmask:  size = ',secat.a
print, 'setmask:  ellipticity =',secat.ellip
print, 'setmask:  PA =',secat.astropa
print, 'setmask: Therefore:'
print, 'setmask:  (b/a) = ',1.0 - secat.ellip
print, 'setmask:  a = ',secat.a / (1.0 - secat.ellip)
ner = [nx, ny]
astropa = secat.astropa
stampx0 = round(nx / 2)
stampy0 = round(ny / 2)
shaperatio = 1.0 / (1.0 - secat.ellip) ; a/b = 1 / (1 - ellipticity)
masksize = maskfac*secat.a*shaperatio
dist_ellipse, er, ner, stampx0, stampy0, shaperatio, astropa
if masksize gt 20 then emax = masksize else emax = 20
print, ''
print, 'setmask: Masking pixels outside ellipse defined by:'
print, 'setmask:   Center = ',stampx0,stampy0
print, 'setmask:   Semimajor axis = ',emax
print, 'setmask:   (a/b) = ',shaperatio
print, 'setmask:   PA = ',astropa
wbad = where(er gt emax, nbad)
if nbad gt 0 then goodmask[wbad]=0

; Check output

wgood = where(goodmask gt 0, ngood)
if ngood eq 0 then begin
    print, ''
    print, '*** ERROR: setmask.  No good pixels in area mask! ***'
    print, '           Cannot continue!'
    return
endif
print, 'setmask: Mask set with ',ngood,' good pixels and ',nbad,' bad pixels.'

retval = 0

end
