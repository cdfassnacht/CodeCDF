pro plot_panel, im, column, row, setup, panel, subpanel, noerase=noerase, $
                _extra=_extra

;
; Procedure plot_panel
;
; Description: Plots an image in one panel of a multi-panel  plot.
;
; Inputs: im          (bytearr)    Image to be plotted in position 1
;         setup       (structure)  Container for all the plotting info.
;         panel       (?)          Specification of plot position
;         subpanel    (?)          Specification of plot position
;         [/noerase]               If set, don't erase before drawing first
;                                   position.
;
; Revision history:
;  2003May15 Chris Fassnacht -- Derived from plot_4x
;

; Check input format

if n_params() lt 6 then begin
    print, ''
    print, $
      'syntax: plot_panel, im, column, row, setup, panel, subpanel'
    print, '                [,/noerase,_extra=_extra]'
    print, ''
    return
endif

; Determine if x-axes should be labeled

notick = replicate(' ',20)
if (setup.doxlab) then begin
   xlab='arcsec' 
   ticknx=''
endif else begin
   xlab=''
   ticknx=notick
endelse
if (setup.doylab) then begin
   ylab='arcsec' 
   tickny=''
endif else begin
   ylab=''
   tickny=notick
endelse

; Plot and label 

plotimage,im,/preserve,_extra=_extra, $
   panel=panel[column,row,*],subpan=subpanel[column,row,*], $
   imgxr=[-1,1]*setup.fullxsz, imgyr=[-1,1]*setup.fullysz, $
   xrange=[-1,1]*setup.imgxsz,yrange=[-1,1]*setup.imgysz, $
   xtit=xlab,ytit=ylab,xtickn=ticknx,ytickn=tickny,$
   color=setup.axiscolor,noerase=noerase
xyouts,setup.labx,setup.laby,setup.lab1,/data,font=0, $
   color=setup.namecolor
xyouts,setup.namex,setup.namey,setup.name,/data,align=0.5,font=0, $
   color=setup.namecolor


end
