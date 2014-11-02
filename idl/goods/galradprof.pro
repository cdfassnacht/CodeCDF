pro galradprof, image, x0, y0, annulus, annct
;
; Procedure avect_rings
;
; Description: Takes an input image, a central position, and a vector 
;               indicating radii of concentric annuli and returns a vector
;               giving the average number of counts within each annulus.
;              ** NB: There is no background subtraction done.  This 
;                 procedure assumes that the input image has been 
;                 background-subtracted.
;
; Inputs: floatarray image        input image
;         float x0                x position of center
;         float y0                y position of center
;         floatarray annulus      radii of concentric annuli
;
; Outputs (accessed through passed parameters):
;         floatarray annct        average counts in each annulus
;
; Revision history:
;  2002Apr09 Chris Fassnacht -- First rough version.
;

; Check input format

if n_params() ne 5 then begin
    print, 'ERROR: galradprof'
    print, 'syntax: galradprof, image, x0, y0, annulus, annct'
    print, ''
endif

n = long(size(annulus, /n_dimen) - 1)
for i=0L, n, 1 do begin
    



rout = FLTARR( n + 1 )
iout = FLTARR( n + 1 )
j = LONG( 0 )
FOR r=0.0, rmax, dr DO BEGIN
   Sumann, image, xcen, ycen, r, r+dr, 0.0, area, psum, nsum, d1, d2, d3, d4
   rout[ j ] = r + dr / 2
   iout[ j ] = ( psum + nsum ) / area - skyback
   j = j + 1
ENDFOR


end

