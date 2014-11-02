pro plotcat_aitoff, catfile, format=format, psym=psym, psfile=psfile, $
   jpfile=jpfile, _extra=_extra

;
; Function plotcat_aitoff
;
; Description: Takes an input catalog of (for now) two columns, the
;  first being the RA in decimal degrees and the second being the Dec
;  in decimal degrees, and plots the points on an Aitoff projection
;  of the sky
;
; Inputs: catfile    (string)      name of input catalog file
;         [format=]  (int)         format of input catalog file:
;                                   0 => ra, dec
;                                   1 => id, ra, dec
;                                   In both cases, ra and dec are in
;                                    decimal degrees
;                                   default = 0
;         [psym=]    (int)         plotting symbol (default = 6, which
;                                   gives open squares)
;         [psfile=]  (string)      name for output postscript file (if not set,
;                                   no file is created)
;         [jpfile=]  (string)      name for output jpeg file (if not set,
;                                   no file is created)
;         _extra     ----          any extra parameters to be passed
;                                   to procedures called by this procedure
;
; Output: [none for right now]
;
; Revision history:
;  2006Sep18 Chris Fassnacht -- First working version
;

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, 'syntax: plotcat_aitoff, catfile'
    print, '         [,format=format,psym=psym,psfile=psfile,jpfile=jpfile]'
    print, '   Format values:'
    print, '     0 ==> ra, dec'
    print, '     1 ==> id, ra, dec'
    print, '    NB: ra and dec are in decimal degrees'
    print, ''
    print, 'Default values:'
    print, '  format:   0'
    print, '  psym:     6 (open squares)'
    print, ''
    return
endif

; Set things to default values if not set on command line

if n_elements(format) eq 0 then format=0
if n_elements(psym) eq 0 then psym=6


; Read in the catalog from the input file

case format of
   0:    begin
      print, 'read_secat: Reading catalog from ',catfile
      readcol, catfile, ra, dec, comment='#',format=''
    end
   1:    begin
      print, 'read_secat: Reading catalog from ',catfile
      readcol, catfile, id, ra, dec, comment='#',format='a,f,f'
    end
   else: begin
      print, ''
      print, 'ERROR: read_secat.  Not a valid format'
      print, ''
      return
    end
endcase

; Set up the Aitoff projection grid

aitoff_grid, /label,/new

; Convert RA and Dec to x, y on the Aitoff grid

aitoff, ra, dec, x, y

; Actually plot the points

plots, x, y, psym=psym, color=255, _extra=_extra

; Plot to a postscript file if requested

if n_elements(psfile) gt 0 then begin
   print,' Saving image in postscript file '+psfile+'.ps.'
   ps_open,psfile,/color,/ps_fonts
   device,xsize=6,ysize=4,/inches,/portrait,/times
   aitoff_grid, /label, /new
   plots, x, y, psym=psym
   ps_close
endif


end
