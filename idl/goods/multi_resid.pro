pro multi_resid, catfile, istart, prefix=prefix, resid1=resid1, resid2=resid2, $
   _extra=_extra

;
; Procedure multi_resid
;
; Description: Interactive interface for calling plot_resid.pro.  The
;               mouse is used to control the calling sequence.
;                 * right mouse button to move forward in your list
;                 * left mouse button to move back in your list
;                 * middle mouse button to stop
;
; Inputs: catfile     (string)     name of catalog file.  The file is
;                                   assumed to have (at least) one
;                                   columns. with the first column being the
;                                   id number or full source name
;         istart      (int)        optional variable to start at somewhere
;                                   besides the first object in the
;                                   catalog.
;         [prefix=]   (string)     optional source-name prefix, if candlist
;                                   contains only id numbers
;         [resid1=]   (string)     name of directory containing residuals
;                                   default is 'Resid'
;         [resid2=]   (string)     optional second directory of residuals,
;                                   used if side-by-side comparision is
;                                   desired.
;
; Revision history:
;  2003Feb27 Chris Fassnacht -- A modification of Lexi's eidolwrap.pro.
;                                (first working version).
;  2003Feb28 Chris Fassnacht -- Added some error checking for beginning
;                                and end of list.
;  2003Mar03 Chris Fassnacht -- Realized that the input catalog only had
;                                to contain the id number of the sources,
;                                and fixed code to reflect this.
;  2003Mar05 Chris Fassnacht -- Added option to plot a second set of
;                                residuals for side-by-side comparision
;  2003Mar13 Chris Fassnacht -- Changed button order so that (right-mouse)
;                                now advances in the list and (left-mouse)
;                                goes back
;  2003Jul02 Chris Fassnacht -- Changed passed "root" parameter into
;                                optional "prefix" parameter.
;

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, $
      'syntax: multi_resid, catfile, (istart)'
    print, '   [,prefix=prefix,resid1=resid1,resid2=resid2]'
    print, ''
    print, " Default value for resid1 is 'Resid'"
    print, ' resid2 is only used if a side-by-side comparision of residuals'
    print, '  is desired.  There is no default setting.'
    print, ''
    return
endif

; Read in catalog with readcol

print, ''
print, 'Reading catalog from ',catfile
readcol, catfile, id, format='A', comment='#'
neid=n_elements(id)

; If optional input variables have not been set, use defaults

if not keyword_set(resid1) then resid1 = 'Resid'
if n_elements(istart) eq 0 then istart=0l
if istart lt 0 or istart gt (neid-1) then begin
   print, ''
   print,'ERROR: multi_resid.  Need to give start value in range 0 --'+$
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
       if keyword_set(prefix) then idstr=prefix+'_'+idstr
       if keyword_set(resid2) then $
         plot_resid, idstr, resid1=resid1, resid2=resid2, _extra=_extra $
          else plot_resid, idstr, resid1=resid1, _extra=_extra
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
   cursor,x,y,/down
   mousebutt = !mouse.button
   if mousebutt eq 4 then i=i+1
   if mousebutt eq 1 then i=i-1
endwhile
        
end 


