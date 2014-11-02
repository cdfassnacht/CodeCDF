pro showpost, image, xsearch, ysearch, stampsize
;
; Procedure showpost
;
; Description: A driver to show postage stamp cutouts of a large image
;  given an input list of coordinates.
;
; Inputs: image      (floatarray)  input image
;         xsearch    (floatarray)  vector of x coordinates of input objects
;         ysearch    (floatarray)  vector of y coordinates of input objects
;
; Revision history:
;  2002Apr18 Chris Fassnacht -- First working version.
;  2003Mar18 Chris Fassnacht -- Changed to reflect new version of sigclip
;

; Check input format

if n_params() lt 4 then begin
    print, ''
    print, $
      'syntax: showpost, image, xsearch, ysearch, stampsize'
    print, ''
    return
endif

; Loop on input list, calling poststamp, sigclip, and display for each object

stmean=1.0
strms = 1.0
for i=1,(size(xsearch,/dimen))[0] do begin
    tmpim = poststamp(image,xsearch[i-1],ysearch[i-1],stampsize,stampsize)
    sigclip, tmpim, ststat, 3.0
    strms = ststat.rms
    stmean = ststat.mean
    print, ''
    print, '*** Displaying object ',strtrim(i,1),'. ***'
    display, tmpim, dmin=stmean-1.5*strms, dmax=stmean+20*strms
endfor

end
