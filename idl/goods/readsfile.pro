pro readsfile, fitsfile, outimage, outhead, catfile, outvec
;
; Procedure readsfile
;
; Description: Reads the files needed for the lens search algorithms
;  driven by dosearch.pro.
;
; Inputs: fitsfile    (string)     name of input FITS file
;         catfile     (string)     name of catalog file
;
; Outputs accessed through passed parameters:
;         outimage    (fltarray)   output image array
;         outhead     (strarray)   FITS header of outimage
;         outvec      (fltarray)   a 5 x N array, where each of the
;                                   N rows is associated with one
;                                   row in the catalog file.  The
;                                   five columns are:
;                                    [0]: x position of object
;                                    [1]: y position of object
;                                    [2]: semimajor axis of object
;                                    [3]: semiminor axis of object
;                                    [4]: PA of semimajor axis 
;                                      (astronomical convention)
;
; Revision history:
;  2002Apr25 Chris Fassnacht -- First working version.
;

; Check input format

if n_params() lt 5 then begin
    print, ''
    print, 'syntax: readsfile, fitsfile, outimage, outhead, catfile, outvec'
    print, ''
    return
endif

; Read in fits file with mrdfits
; ** NB: For now assume no errors

outimage = mrdfits(fitsfile, 0, outhead)

; Read in catalog with readcol

readcol, catfile, x, y, a, b, pa

; Put catalog parameters into output vector

if n_elements(x) gt 0 then begin
    outvec = fltarr(5,n_elements(x))
    outvec[0,*] = x
    outvec[1,*] = y
    outvec[2,*] = a
    outvec[3,*] = b
    outvec[4,*] = pa
endif else outvec = fltarr(5)

end
