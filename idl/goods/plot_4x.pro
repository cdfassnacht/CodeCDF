pro plot_4x, im1, im2, im3, im4, setup4, panel, subpanel, noerase=noerase

;
; Procedure cand_plot
;
; Description: Plots 4 images in a 4-panel plot, with the images all in
;               one row.  Position 1 is the leftmost position, position 4
;               is the rightmost position.
;
; Inputs: im1         (bytearr)    Image to be plotted in position 1
;         im2         (bytearr)    Image to be plotted in position 2
;         im3         (bytearr)    Image to be plotted in position 3
;         im4         (bytearr)    Image to be plotted in position 4
;         setup4      (structure)  Container for all the plotting info.
;         panel       (?)          Specification of plot position
;         subpanel    (?)          Specification of plot position
;         [/noerase]               If set, don't erase before drawing first
;                                   position.
;
; Revision history:
;  2003Mar20 Chris Fassnacht -- Moved out of plot_cand.pro
;  2003Apr11 Chris Fassnacht -- Made the color of the text in the first
;                                window part of the passed setup4 structure.
;

; Check input format

if n_params() lt 7 then begin
    print, ''
    print, $
      'syntax: plot_4x, im1, im2, im3, im4, setup, panel, subpanel'
    print, ''
    return
endif

; Determine if x-axes should be labeled

notick = replicate(' ',20)
if (setup4.doxlab) then begin
   xlab='arcsec' 
   ticknx=''
endif else begin
   xlab=''
   ticknx=notick
endelse
if (setup4.doylab) then begin
   ylab='arcsec' 
   tickny=''
endif else begin
   ylab=''
   tickny=notick
endelse

; Plot and label position 1

plotimage,im1,/preserve,_extra=_extra, $
   panel=panel[0,setup4.row,*],subpan=subpanel[0,setup4.row,*], $
   imgxr=[-1,1]*setup4.fullxsz, imgyr=[-1,1]*setup4.fullysz, $
   xrange=[-1,1]*setup4.imgxsz,yrange=[-1,1]*setup4.imgysz, $
   xtit=xlab,ytit=ylab,xtickn=ticknx,ytickn=tickny,$
   color=setup4.axiscolor,noerase=noerase
xyouts,setup4.labx,setup4.laby,setup4.lab1,/data,font=0, $
   color=setup4.namecolor
xyouts,setup4.namex,setup4.namey,setup4.name,/data,align=0.5,font=0, $
   color=setup4.namecolor

; Plot and label position 2

plotimage,im2,/preserve,_extra=_extra, $
   panel=panel[1,setup4.row,*],subpan=subpanel[1,setup4.row,*], $
   imgxr=[-1,1]*setup4.fullxsz, imgyr=[-1,1]*setup4.fullysz, $
   xrange=[-1,1]*setup4.imgxsz,yrange=[-1,1]*setup4.imgysz, $
   xtit=xlab,xtickn=ticknx,ytickn=notick,color=setup4.axiscolor,/noerase
xyouts,setup4.labx,setup4.laby,setup4.lab2,/data,color=setup4.axiscolor,font=0

; Plot and label position 3

plotimage,im3,/preserve,_extra=_extra, $
   panel=panel[2,setup4.row,*],subpan=subpanel[2,setup4.row,*], $
   imgxr=[-1,1]*setup4.fullxsz, imgyr=[-1,1]*setup4.fullysz, $
   xrange=[-1,1]*setup4.imgxsz,yrange=[-1,1]*setup4.imgysz, $
   xtit=xlab,xtickn=ticknx,ytickn=notick,color=setup4.axiscolor,/noerase
xyouts,setup4.labx,setup4.laby,setup4.lab3,/data,color=setup4.axiscolor,font=0

; Plot and label position 4

plotimage,im4,/preserve,_extra=_extra, $
   panel=panel[3,setup4.row,*],subpan=subpanel[3,setup4.row,*], $
   imgxr=[-1,1]*setup4.fullxsz, imgyr=[-1,1]*setup4.fullysz, $
   xrange=[-1,1]*setup4.imgxsz,yrange=[-1,1]*setup4.imgysz, $
   xtit=xlab,xtickn=ticknx,ytickn=notick,color=setup4.axiscolor,/noerase
xyouts,setup4.labx,setup4.laby,setup4.lab4,/data,color=setup4.axiscolor,font=0

end
