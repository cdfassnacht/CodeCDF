pro def_secat, secat
;
; Procedure def_secat
;
; Description: Creates a structure that contains information from a
;               SExtractor catalog.
;
; Inputs: secat (struct)          output structure containing secat info
;
; Revision history:
;  2003Apr11 Chris Fassnacht -- First version.
;

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, 'syntax: def_secat, secat'
    print, ''
    return
endif

; Set up output structure

secat =  create_struct('id',0,$
                      'x',0.0,$
                      'y',0.0,$
                      'alpha',0.d0,$
                      'delta',0.d0,$
                      'isoarea',0.0,$
                      'pa',0.0,$
                      'astropa',0.0,$
                      'ellip',0.0,$
                      'r_kron',0.0,$
                      'r_flux1',0.0,$
                      'r_flux2',0.0,$
                      'r_flux3',0.0,$
                      'a',0.0,$
                      'b',0.0,$
                      'fwhm',0.0,$
                      'class',0.0,$
                      'flag',0.0,$
                      'm_iso',0.0,$
                      'merr_iso',0.0 $
                      )

end

