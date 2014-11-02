pro multi_mydisp, catfile, root, istart, _extra=_extra

;
; Procedure multi_mydisp
;
; Description: Interactive interface for calling mydisp.pro.  The
;               mouse is used to control the calling sequence.
;                 * right mouse button to move forward in your list
;                 * left mouse button to move back in your list
;                 * middle mouse button to stop
;
; Inputs: catfile     (string)     name of catalog file.  The file is
;                                   assumed to have (at least) one
;                                   columns. with the first column being the
;                                   id number
;         root        (string)     root part of name for all input files.
;                                   Along with the id, defines the filenames.
;         istart      (int)        optional variable to start at somewhere
;                                   besides the first object in the
;                                   catalog.
;
; Revision history:
;  2003Mar31 Chris Fassnacht -- First working version.
;

; Check input format

if n_params() lt 2 then begin
    print, ''
    print, $
      'syntax: multi_mydisp, catfile, root, (istart)'
    print, '   [,_extra=_extra]'
    print, ''
    return
endif

; Set directory names and defaults

if not keyword_set(indir) then indir='.'

; Read in catalog with readcol

print, ''
print, 'Reading catalog from ',catfile
readcol, catfile, id, format='A'
neid=n_elements(id)

; If optional input variables have not been set, use defaults

if n_elements(istart) eq 0 then istart=0l
if istart lt 0 or istart gt (neid-1) then begin
   print, ''
   print,'ERROR: multi_mydisp.  Need to give start value in range 0 --'+$
      strcompress(neid-1)
   print, '*** Exiting ***'
   print, ''
   return
endif

; Start the loop

i=istart
mousebutt = 0
while mousebutt ne 2 do begin
   if(i ge 0 and i lt neid) then begin
      print,strcompress(i)+': '+id[i]
      basestr = id[i]
      case strlen(basestr) of
         1: idstr = '0000'+basestr
         2: idstr = '000'+basestr
         3: idstr = '00'+basestr
         4: idstr = '0'+basestr
         else: idstr = basestr
       endcase
       fitsfile = root+'_'+idstr+'_z.fits'
       mydisp, fitsfile, _extra=_extra
   endif else begin
      if i lt 0 then begin
         print, ''
         print, '*** Already at beginning of list ***'
         i = i + 1
        print, ''
      endif else begin
         print, ''
         print, '*** Already at end of list ***'
         i = i - 1
         print, ''
      endelse
   endelse
   print, 'Hit left mouse button to go to next frame'
   print, 'Hit right mouse button to go to previous frame'
   print, 'Hit middle mouse button to quit'
   print, ''
   cursor,x,y,/wait
   mousebutt = !mouse.button
   if mousebutt eq 4 then i=i+1
   if mousebutt eq 1 then i=i-1
endwhile
        
end 


