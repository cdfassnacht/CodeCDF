pro plot_7x, im1, im2, im3, im4, im5, im6, im7, setup, panel, subpanel, $
             noerase=noerase

;
; Procedure cand_plot
;
; Description: Plots 7 images in a 7-panel plot, with the images all in
;               one row.  Position 1 is the leftmost position, position 7
;               is the rightmost position.
;
; Inputs: im1         (bytearr)    Image to be plotted in position 1
;         im2         (bytearr)    Image to be plotted in position 2
;         im3         (bytearr)    Image to be plotted in position 3
;         im4         (bytearr)    Image to be plotted in position 4
;         im5         (bytearr)    Image to be plotted in position 5
;         im6         (bytearr)    Image to be plotted in position 6
;         im7         (bytearr)    Image to be plotted in position 7
;         setup       (structure)  Container for all the plotting info.
;         panel       (?)          Specification of plot position
;         subpanel    (?)          Specification of plot position
;         [/noerase]               If set, don't erase before drawing first
;                                   position.
;
; Revision history:
;  2003Mar20 Chris Fassnacht -- Moved out of plot_cand.pro
;  2003Apr11 Chris Fassnacht -- Made the color of the text in the first
;                                window part of the passed setup structure.
;

; Check input format

if n_params() lt 10 then begin
    print, ''
    print, $
      'syntax: plot_7x, im1,im2,im3,im4,im5,im6,im7,setup,panel,subpanel
    print, ''
    return
endif

; Determine if x-axes should be labeled

if (setup.doxlab) then xlab='arcsec' else xlab=''

; Plot and label position 1

plotimage,im1,/preserve,_extra=_extra, $
   panel=panel[0,setup.row,*],subpan=subpanel[0,setup.row,*], $
   imgxr=[-1,1]*setup.fullxsz, imgyr=[-1,1]*setup.fullysz, $
   xrange=[-1,1]*setup.imgxsz,yrange=[-1,1]*setup.imgysz, $
   xtit=xlab,ytit='arcsec',color=setup.axiscolor,noerase=noerase
;xyouts,setup.labx,setup.laby,setup.lab1,/data,color=setup.axiscolor,font=0
;xyouts,setup.namex,setup.namey,setup.name,/data,color=setup.axiscolor,align=0.5,font=0
xyouts,setup.labx,setup.laby,setup.lab1,/data,font=0, $
   color=setup.namecolor
xyouts,setup.namex,setup.namey,setup.name,/data,align=0.5,font=0, $
   color=setup.namecolor

; Plot and label position 2

plotimage,im2,/preserve,_extra=_extra, $
   panel=panel[1,setup.row,*],subpan=subpanel[1,setup.row,*], $
   imgxr=[-1,1]*setup.fullxsz, imgyr=[-1,1]*setup.fullysz, $
   xrange=[-1,1]*setup.imgxsz,yrange=[-1,1]*setup.imgysz, $
   xtit=xlab,ytickn=replicate(' ',20),color=setup.axiscolor,/noerase
xyouts,setup.labx,setup.laby,setup.lab2,/data,color=setup.axiscolor,font=0

; Plot and label position 3

plotimage,im3,/preserve,_extra=_extra, $
   panel=panel[2,setup.row,*],subpan=subpanel[2,setup.row,*], $
   imgxr=[-1,1]*setup.fullxsz, imgyr=[-1,1]*setup.fullysz, $
   xrange=[-1,1]*setup.imgxsz,yrange=[-1,1]*setup.imgysz, $
   xtit=xlab,ytickn=replicate(' ',20),color=setup.axiscolor,/noerase
xyouts,setup.labx,setup.laby,setup.lab3,/data,color=setup.axiscolor,font=0

; Plot and label position 4

plotimage,im4,/preserve,_extra=_extra, $
   panel=panel[3,setup.row,*],subpan=subpanel[3,setup.row,*], $
   imgxr=[-1,1]*setup.fullxsz, imgyr=[-1,1]*setup.fullysz, $
   xrange=[-1,1]*setup.imgxsz,yrange=[-1,1]*setup.imgysz, $
   xtit=xlab,ytickn=replicate(' ',20),color=setup.axiscolor,/noerase
xyouts,setup.labx,setup.laby,setup.lab4,/data,color=setup.axiscolor,font=0

; Plot and label position 5

plotimage,im5,/preserve,_extra=_extra, $
   panel=panel[4,setup.row,*],subpan=subpanel[4,setup.row,*], $
   imgxr=[-1,1]*setup.fullxsz, imgyr=[-1,1]*setup.fullysz, $
   xrange=[-1,1]*setup.imgxsz,yrange=[-1,1]*setup.imgysz, $
   xtit=xlab,ytickn=replicate(' ',20),color=setup.axiscolor,/noerase
xyouts,setup.labx,setup.laby,setup.lab5,/data,color=setup.axiscolor,font=0

; Plot and label position 6

plotimage,im6,/preserve,_extra=_extra, $
   panel=panel[5,setup.row,*],subpan=subpanel[5,setup.row,*], $
   imgxr=[-1,1]*setup.fullxsz, imgyr=[-1,1]*setup.fullysz, $
   xrange=[-1,1]*setup.imgxsz,yrange=[-1,1]*setup.imgysz, $
   xtit=xlab,ytickn=replicate(' ',20),color=setup.axiscolor,/noerase
xyouts,setup.labx,setup.laby,setup.lab6,/data,color=setup.axiscolor,font=0

; Plot and label position 7

plotimage,im7,/preserve,_extra=_extra, $
   panel=panel[6,setup.row,*],subpan=subpanel[6,setup.row,*], $
   imgxr=[-1,1]*setup.fullxsz, imgyr=[-1,1]*setup.fullysz, $
   xrange=[-1,1]*setup.imgxsz,yrange=[-1,1]*setup.imgysz, $
   xtit=xlab,ytickn=replicate(' ',20),color=setup.axiscolor,/noerase
xyouts,setup.labx,setup.laby,setup.lab7,/data,color=setup.axiscolor,font=0

end
