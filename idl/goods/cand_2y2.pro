pro cand_2y2, candlist, prefix=prefix, rsddir=rsddir, fitdir=fitdir, $
   imgsz=imgsz, hsig=hsig, rhsig=rhsig, psfile=psfile, _extra=_extra

;
; Procedure cand_2y2
;
; Description: Makes a multi-panel plot containing true-color images of
;               a galaxy and residual images created by subtracting off
;               some form of fit to the galaxy light distribution.
;               Like cand_stack, but instead of each column containing
;               4 images associated with one target (BVi,z,fit,resid),
;               each column has 4 images, with 2 each for 2 targets
;               (BVi and resid).  So, for N input galaxies, the
;               output is a  N/2x4 grid of panels.
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
;  2003May14 Chris Fassnacht -- A slight modification of cand_stack.pro
;                                Very kludgy with no error checking.


; Check input format

if n_params() lt 1 then begin
    print, ''
    print, 'Description: Prints out a vertical stack of N/2 '
    print, '  plot_cand-like plots.  The input list (candlist) contains'
    print, '  the ID numbers of the lens candidates to be plotted.'
    print, '  If N>14 the output plot will not fit on a standard letter-sized '
    print, '   page in postscript mode.'
    print, ''
    print, 'syntax: cand_2y2, candlist [,prefix=prefix,,rsddir=rsddir,'
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
ysize = 4.0 * psfigsz + psbordsz

; Read in catalog with readcol and set number of sources.

print, ''
print, 'Reading candidate list from ',candlist
readcol, candlist, id, format='A'
neid=n_elements(id)

; Set up viewport window 
; ***NB: 7 is the max to fit on a standard letter page in postscript file.***

xcells = replicate(1,neid/2)
subcellarray, xcells, [1,1,1,1], newpan, newsubpan
xsize = 1.0 * (neid/2) * psfigsz + psbordsz


; Plot postscript if desired

if n_elements(psfile) gt 0 then begin
   print,'Saving image in postscript file '+psfile+'.ps.'
   if (xsize lt 8.5) then xoff=0.5*(8.5-ysize) else xoff=0.5
   yoff = 0.5*(11.0-ysize)
   ps_open,psfile,/color,/ps_fonts
;   device,xsize=6,ysize=3,/inches,/portrait,/times
   device,bits=8,/portrait,/times,xsize=xsize,ysize=ysize,yoffset=yoff,$
      xoffset=xoff,/inches
   setup.tcolor1 = 0
   setup.tcolor3 = 0
   loadct,0
   invct
endif

; Start the loop

for i=0,((neid/2)-1) do begin

   ; Set filenames

   i1 = 2 * i
   i2 = (2 * i) + 1
   if n_elements(prefix) gt 0 then bigroot1 = prefix+'_'+id[i1] else $
      bigroot1 = id[i1]
   if n_elements(prefix) gt 0 then bigroot2 = prefix+'_'+id[i2] else $
      bigroot2 = id[i2]
   bfits1=bdir+'/'+bigroot1+'_b.fits'
   vfits1=vdir+'/'+bigroot1+'_v.fits'
   ifits1=idir+'/'+bigroot1+'_i.fits'
   rsd1fits=rsddir+'/'+bigroot1+'_z_resid.fits'
   bfits2=bdir+'/'+bigroot2+'_b.fits'
   vfits2=vdir+'/'+bigroot2+'_v.fits'
   ifits2=idir+'/'+bigroot2+'_i.fits'
   rsd2fits=rsddir+'/'+bigroot2+'_z_resid.fits'
   
   ; Read in files
   
   print, ''
   print, 'cand_2y2: Reading input files:'
   print, 'cand_2y2:   ',bfits1
   bim1 = mrdfits(bfits1, 0, bhead1)
   print, 'cand_2y2:   ',vfits1
   vim1 = mrdfits(vfits1, 0, vhead1)
   print, 'cand_2y2:   ',ifits1
   iim1 = mrdfits(ifits1, 0, ihead1)
   print, 'cand_2y2:   ',rsd1fits
   rsd1im = mrdfits(rsd1fits, 0, rsd1head)
   print, 'cand_2y2:   ',bfits2
   bim2 = mrdfits(bfits2, 0, bhead2)
   print, 'cand_2y2:   ',vfits2
   vim2 = mrdfits(vfits2, 0, vhead2)
   print, 'cand_2y2:   ',ifits2
   iim2 = mrdfits(ifits2, 0, ihead2)
   print, 'cand_2y2:   ',rsd1fits
   rsd2im = mrdfits(rsd2fits, 0, rsd2head)
   sz=size(bim1,/dimen)
   print, ''
   
   
   ; Set up scaled versions of the input files
   
   imscal, bim1, bimsc1, bstat, lsig=1.5, hsig=hsig, nsig=3.0
   imscal, vim1, vimsc1, vstat, lsig=1.5, hsig=hsig, nsig=3.0
   imscal, iim1, iimsc1, istat, lsig=1.5, hsig=hsig, nsig=3.0
   imscal, rsd1im, reimsc1, restat, lsig=1.5, hsig=rhsig, nsig=3.0
   imscal, bim2, bimsc2, bstat, lsig=1.5, hsig=hsig, nsig=3.0
   imscal, vim2, vimsc2, vstat, lsig=1.5, hsig=hsig, nsig=3.0
   imscal, iim2, iimsc2, istat, lsig=1.5, hsig=hsig, nsig=3.0
   imscal, rsd2im, reimsc2, restat, lsig=1.5, hsig=rhsig, nsig=3.0
   
   ; Set up 3-plane array for inputs to true-color images
   
   bvi1=bytarr(3,sz[0],sz[1])
   bvi1[0,*,*] = iimsc1
   bvi1[1,*,*] = vimsc1
   bvi1[2,*,*] = bimsc1
   bvi2=bytarr(3,sz[0],sz[1])
   bvi2[0,*,*] = iimsc2
   bvi2[1,*,*] = vimsc2
   bvi2[2,*,*] = bimsc2
  
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
           setup.imgxsz2 = imgsz / 2.
           setup.imgysz2 = imgsz / 2.
        endif else begin
           setup.imgxsz = imgsz[i1] / 2.
           setup.imgysz = imgsz[i1] / 2.
           setup.imgxsz2 = imgsz[i2] / 2.
           setup.imgysz2 = imgsz[i2] / 2.
           doxlab = 0
           doylab = 0
        endelse
   endif else begin
   	setup.imgxsz = setup.fullxsz
   	setup.imgysz = setup.fullysz
   	setup.imgxsz2 = setup.fullxsz
   	setup.imgysz2 = setup.fullysz
   endelse
   print, 'Displayed image size (1) will be:',setup.imgxsz*2., $
   	 ' arcsec x',setup.imgysz*2.,' arcsec'
   print, 'Displayed image size (2) will be:',setup.imgxsz2*2., $
   	 ' arcsec x',setup.imgysz2*2.,' arcsec'
   print, 'For these plots, hsig = ',hsig,' and rhsig = ',rhsig
   
   ; Fill in other setup parameters
   
   setup.column = i
   setup.labx = -0.9*setup.imgxsz
   setup.laby = -0.8*setup.imgysz
   setup.labx2 = -0.9*setup.imgxsz2
   setup.laby2 = -0.8*setup.imgysz2
   setup.namex = 0.0
   setup.namey = 0.8*setup.imgysz
   setup.namex2 = 0.0
   setup.namey2 = 0.8*setup.imgysz2
   setup.name = bigroot1
   setup.name2 = bigroot1
   setup.name3 = bigroot2
   setup.name4 = bigroot2
   setup.lab1 = 'BVi'
   setup.lab2 = 'Resid'
   setup.lab3 = 'BVi'
   setup.lab4 = 'Resid'
   if (i eq 0) then setup.doylab = doylab else setup.doylab = 0
   setup.doxlab=doxlab
   
   ; Plot images
   
   print, ''
   print, 'Plotting data for ',bigroot1,bigroot2,'   --',i1,i2
   print, ''
   if (i eq 0) then begin
      plot_2y2, bvi1, reimsc1,bvi2, reimsc2, $
        setup, newpan, newsubpan
   endif else begin
      plot_2y2, bvi1, reimsc1,bvi2, reimsc2, $
        setup, newpan, newsubpan, /noerase
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
