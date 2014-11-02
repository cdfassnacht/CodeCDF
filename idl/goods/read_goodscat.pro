pro read_goodscat, catfile, secat, rotangle=rotangle
;
; Procedure read_goodscat
;
; Description: A procedure to read in a SExtractor catalog produced by
;               the GOODS processing for the ACS data.
;
; Inputs: catfile     (string)     name of catalog file.  The file is
;                                   assumed to be in the GOODS order.
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
      'syntax: read_goodscat, catfile, secat [,rotangle=rotangle]'
    print, ''
    print, ''
    return
endif

; Set up output structure

def_secat, temp

; Read in catalog with readcol

print, 'read_goodscat: Reading catalog from ',catfile
format1 ='i,f,f,d,d,f,f,f,x,x,f,f,f,f,f,f,f'
readcol, catfile, id, x, y, alpha, delta, isoarea, theta, ellip, r_kron, $
   r_flux1, r_flux2, r_flux3, fwhm, class, flag, format=format1, comment='#'

; Get number of catalog entries 

ncat = n_elements(id)

; Convert PA to astronomical convention, and add any additional rotation
;  if necessary

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
       secat[i].alpha = alpha[i]
       secat[i].delta = delta[i]
       secat[i].isoarea = isoarea[i]
       secat[i].pa = pa[i]
       secat[i].ellip = ellip[i]
       secat[i].r_kron = r_kron[i]
       secat[i].r_flux1 = r_flux1[i]
       secat[i].r_flux2 = r_flux2[i]
       secat[i].r_flux3 = r_flux3[i]
       secat[i].fwhm = fwhm[i]
       secat[i].class = class[i]
       secat[i].flag = flag[i]
    endfor
endif else begin
    print, ''
    print, '*** ERROR: read_goodscat.  No objects found in catalog list. ***'
    print, ''
endelse

end
