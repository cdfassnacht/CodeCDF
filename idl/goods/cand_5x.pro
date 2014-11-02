pro cand_5x, candlist, prefix=prefix, rsddir=rsddir, fitdir=fitdir, $
   imgsz=imgsz, hsig=hsig, rhsig=rhsig, psfile=psfile, _extra=_extra

;
; Procedure cand_5x
;
; Description: Makes a multi-panel plot containing true-color images of
;               a galaxy.
;               Like cand_stack, but instead of each row containing
;               4 images associated with one target (BVi,z,fit,resid),
;               each row has 5 images, one for each target (BVi only).
;               So, for N input galaxies, the
;               output is a  5 x N/5 grid of panels.
;
; Inputs: candlist    (string)     file containing list of sources
;         [prefix=]   (string)     optional source-name prefix, if candlist
;                                   contains only id numbers
;         [rsddir=]   (string)     name of directory containing residuals
;                                   default is 'Resid'
;         [fitdir=]   (string)     name of directory containing model fits
;                                   default is 'Fit'
;         [imgsz=]    (floatarr)   image size in arcsec.  If chosen, will
;                                   produce an output image of imgsz x imgsz
;         [psfile=]   (string)     name for output postscript file (if not set,
;                                   no file is created)
;         [hsig=]     (float)      upper limit for display range in units of
;                                   the rms in the clipped image.  The
;                                   upper limit wil be hsig*rms + mean.
;                                   Default value is 30.0.
;         [rhsig=]    (float)      same as hsig, but for the residual image
;                                   only.  Default value is 10.0.
;
; Revision history:
;  2003May15 Chris Fassnacht -- A slight modification of cand_stack.pro
;                                Very kludgy with no error checking.


; Check input format

if n_params() lt 1 then begin
    print, ''
    print, 'Description: Prints out a vertical stack of N/2 '
    print, '  plot_cand-like plots.  The input list (candlist) contains'
    print, '  the ID numbers of the lens candidates to be plotted.'
    print, '  If N>35 the output plot will not fit on a standard letter-sized '
    print, '   page in postscript mode.'
    print, ''
    print, 'syntax: cand_5x, candlist [,prefix=prefix,,rsddir=rsddir,'
    print, '           fitdir=fitdir,imgsz=imgsz,hsig=hsig,rhsig=rhsig,'
    print, '           psfile=psfile]'
    print, ''
    print, " Default value for rsddir is 'Resid'"
    print, " Default value for fitdir is 'Fit'"
    print, ' hsig sets the upper limit (in units of the RMS) for all '
    print, '   plots except the residual plot.  Default value is 30.0'
    print, ' rhsig sets the upper limit for the residual plot.  Default '
    print, '   value is 10.0'
    print, ''
    return
endif

; Create setup structure

def_pltsetup, setup

; Set directory names and other defaults

bdir = 'b-band'
vdir = 'v-band'
idir = 'i-band'
zdir = 'z-band'
if not keyword_set(rsddir) then rsddir = 'Resid'
if not keyword_set(fitdir) then fitdir = 'Fit'
pixscale = 0.05
if n_elements(hsig) eq 0 then hsig = 30.0
if n_elements(rhsig) eq 0 then rhsig = 10.0
psfigsz = 1.44
psbordsz = 0.43
xsize = 5.0 * psfigsz + psbordsz

; Read in catalog with readcol and set number of sources.

print, ''
print, 'Reading candidate list from ',candlist
readcol, candlist, id, format='A'
neid=n_elements(id)

; Set up viewport window 
; ***NB: 7 is the max to fit on a standard letter page in postscript file.***

ycells = replicate(1,neid/5)
subcellarray, [1,1,1,1,1], ycells, newpan, newsubpan
ysize = 1.0 * (neid/5) * psfigsz + psbordsz


; Plot postscript if desired

if n_elements(psfile) gt 0 then begin
   print,'Saving image in postscript file '+psfile+'.ps.'
   if (xsize lt 8.5) then xoff=0.5*(8.5-xsize) else xoff=0.5
   if (ysize lt 11.0) then yoff=0.5*(11.0 - ysize) else yoff=0.5
   ps_open,psfile,/color,/ps_fonts
;   device,xsize=6,ysize=3,/inches,/portrait,/times
   device,bits=8,/portrait,/times,xsize=xsize,ysize=ysize,yoffset=yoff,$
      xoffset=xoff,/inches
   setup.namecolor = 0
   setup.axiscolor = 0
   loadct,0
   invct
endif

; Set row and column counters

column = -1
row = neid/5 - 1

; Start the loop

for i=0,(neid-1) do begin

   ; Increment row and column counters

   column = column + 1
   if column eq 5 then begin
      column = 0
      row = row - 1
   endif

   ; Set filenames

   if n_elements(prefix) gt 0 then bigroot1 = prefix+'_'+id[i] else $
      bigroot1 = id[i]
   bfits1=bdir+'/'+bigroot1+'_b.fits'
   vfits1=vdir+'/'+bigroot1+'_v.fits'
   ifits1=idir+'/'+bigroot1+'_i.fits'
   
   ; Read in files
   
   print, ''
   print, 'cand_5x: Reading input files:'
   print, 'cand_5x:   ',bfits1
   bim1 = mrdfits(bfits1, 0, bhead1)
   print, 'cand_5x:   ',vfits1
   vim1 = mrdfits(vfits1, 0, vhead1)
   print, 'cand_5x:   ',ifits1
   iim1 = mrdfits(ifits1, 0, ihead1)
   sz=size(bim1,/dimen)
   print, ''
   
   
   ; Set up scaled versions of the input files
   
   imscal, bim1, bimsc1, bstat, lsig=1.5, hsig=hsig, nsig=3.0
   imscal, vim1, vimsc1, vstat, lsig=1.5, hsig=hsig, nsig=3.0
   imscal, iim1, iimsc1, istat, lsig=1.5, hsig=hsig, nsig=3.0
   
   ; Set up 3-plane array for inputs to true-color images
   
   bvi1=bytarr(3,sz[0],sz[1])
   bvi1[0,*,*] = iimsc1
   bvi1[1,*,*] = vimsc1
   bvi1[2,*,*] = bimsc1
  
   ; Set displayed image size and label locations
   
   doxlab=0
   doylab=0
   setup.fullxsz = sz[0]*pixscale / 2.
   setup.fullysz = sz[1]*pixscale / 2.
   print, ''
   print, 'Native image size is ',setup.fullxsz*2.,' arcsec x',$
   	 setup.fullysz*2.,' arcsec'
   if n_elements(imgsz) gt 0 then begin
        if n_elements(imgsz) lt neid then begin
           setup.imgxsz = imgsz / 2.
           setup.imgysz = imgsz / 2.
        endif else begin
           setup.imgxsz = imgsz[i] / 2.
           setup.imgysz = imgsz[i] / 2.
           doxlab = 0
           doylab = 0
        endelse
   endif else begin
   	setup.imgxsz = setup.fullxsz
   	setup.imgysz = setup.fullysz
   endelse
   print, 'Displayed image size will be:',setup.imgxsz*2., $
   	 ' arcsec x',setup.imgysz*2.,' arcsec'
   print, 'For these plots, hsig = ',hsig,' and rhsig = ',rhsig
   
   ; Fill in other setup parameters
   
   setup.column = i
   setup.labx = -0.9*setup.imgxsz
   setup.laby = -0.8*setup.imgysz
   setup.namex = 0.0
   setup.namey = 0.8*setup.imgysz
   setup.name = bigroot1
   setup.lab1 = 'BVi'
   if (i eq 0) then setup.doylab = doylab else setup.doylab = 0
   setup.doxlab=doxlab
   
   ; Plot images
   
   print, ''
   print, 'Plotting data for ',bigroot1,'   --',i
   print, ''
   if (i eq 0) then begin
      plot_panel, bvi1, column, row, setup, newpan, newsubpan
   endif else begin
      plot_panel, bvi1, column, row, setup, newpan, newsubpan, /noerase
   endelse
endfor


; Plot postscript if desired

if n_elements(psfile) gt 0 then begin
   print, ''
   print, 'Final postscript image occupies ',xsize,' x',ysize,' inches'
   print,''
   ps_close
endif


end 
