pro read_lscat, catfile, secat, rotangle=rotangle
;
; Procedure read_lscat
;
; Description: A procedure to read in a SExtractor catalog
;               and extract from it the parameters needed for the strong-lens
;               search algorithms.
;
; Inputs: catfile     (string)     name of catalog file.  The file is
;                                   assumed to be in the GOODS order.
;                                    [0]: id number
;                                    [1]: x position of object
;                                    [2]: y position of object
;                                    [3]: PA of semimajor axis 
;                                      (astronomical convention)
;                                    [4]: ellipticity (1 - b/a)
;                                    [5]: "size" - for now, SExtractor's FWHM
;         secat       (struct)     output structure array containing
;                                   catalog info.
;         [rotangle=] (float)      additional angle to rotate PA.  Note that
;                                   this procedure will automatically convert
;                                   the SExtractor math-convention PA to
;                                   an astronomy-convention PA
;
; Revision history:
;  2003Jun22 Chris Fassnacht -- First version
;

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, 'This procedure reads in a SExtractor catlog'
    print, ''
    print, $
      'syntax: read_lscat, catfile, secat [,rotangle=rotangle]'
    print, ''
    print, ''
    return
endif

; Set up output structure

def_secat, temp

; Read in catalog with readcol

print, 'read_lscat: Reading catalog from ',catfile
format1 ='i,f,f,f,f,f,f'
readcol, catfile, id, x, y, theta, ellip, fwhm, miso, $
   format=format1, comment='#'

; Get number of catalog entries 

ncat = n_elements(id)

; Convert PA to astronomical convention, and add any additional rotation
;  if necessary
;
; *** NB: v1 files needed pa = theta - 90.0 - 22.44 to account
;       for rotation to north up.  v1.7 files don't have this
;       rotation added so pa = theta - 90.0

if keyword_set(rotangle) then pa = theta - 90.0 - rotangle else $
   pa = theta - 90.0

; Fill the structure array

print, ''
if ncat gt 0 then begin

    ; Create output structure array
    secat = replicate(temp,ncat)

    ; Load the array
    for i=0,ncat-1 do begin
       secat[i].id = id[i]
       secat[i].x = x[i]
       secat[i].y = y[i]
       secat[i].pa = pa[i]
       secat[i].ellip = ellip[i]
       secat[i].fwhm = fwhm[i]
       secat[i].m_iso = miso[i]
       secat[i].a = fwhm[i]
    endfor
endif else begin
    print, ''
    print, '*** ERROR: read_lscat.  No objects found in catalog list. ***'
    print, ''
endelse

end
