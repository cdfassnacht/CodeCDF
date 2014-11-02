pro cand_plot, id, prefix=prefix, rsddir=rsddir, fitdir=fitdir, imgsz=imgsz, $
   hsig=hsig, rhsig=rhsig, psfile=psfile, _extra=_extra

;
; Procedure cand_plot
;
; Description: Makes a multi-panel plot containing true-color images of
;               a galaxy and residual images created by subtracting off
;               some form of fit to the galaxy light distribution.
;               Very similar to plot_resid, but with a different set of
;               images and layout.
;
; Inputs: id          (string)     root part of name for all input files.
;                                   May be modified by the optional prefix
;                                   parameter
;         [prefix=]   (string)     optional source-name prefix, if candlist
;                                   contains only id numbers
;         [rsddir=]   (string)     name of directory containing residuals
;                                   default is 'Resid'
;         [fitdir=]   (string)     name of directory containing model fits
;                                   default is 'Fit'
;         [imgsz=]    (float)      image size in arcsec.  If chosen, will
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
;  2003Mar18 Chris Fassnacht -- A revision of plot_resid
;  2003Mar19 Chris Fassnacht -- Added image-size optional passed parameter.
;                               Put a lot of plotting info into a new
;                                structure called setup4.
;  2003Mar20 Chris Fassnacht -- Added optional postscript file parameter.
;                               Moved actual plotting into new plot_4x.pro.
;  2003Mar27 Chris Fassnacht -- Added second row showing B, V, i, z plots.
;                               Added "_extra" parameter.
;  2003Apr01 Chris Fassnacht -- Added fitdir optional parameter.
;  2003Apr09 Chris Fassnacht -- Fixed bugs in setting plot size.
;  2003Apr11 Chris Fassnacht -- Added hsig and rhsig as passed parameters
;                                rather than hard-wiring the values.
;                               Moved definition of the setup structure into
;                                def_pltsetup
;  2003Apr27 Chris Fassnacht -- Changed the passed "root" parameter into an
;                                optional "prefix" parameter
;

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, $
      'syntax: cand_plot, root, [prefix=prefix,rsddir=rsddir,fitdir=fitdir,'
    print, '       imgsz=imgsz,hsig=hsig,rhsig=rhsig,psfile=psfile]'
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

def_pltsetup, setup4

; Set directory names and defaults

bdir = 'b-band'
vdir = 'v-band'
idir = 'i-band'
zdir = 'z-band'
if not keyword_set(rsddir) then rsddir = 'Resid'
if not keyword_set(fitdir) then fitdir = 'Fit'
pixscale = 0.05
if n_elements(hsig) eq 0 then hsig = 30.0
if n_elements(rhsig) eq 0 then rhsig = 10.0
zsig = hsig

; Set filenames

if n_elements(prefix) gt 0 then bigroot = prefix+'_'+id else bigroot = id
bfits=bdir+'/'+bigroot+'_b.fits'
vfits=vdir+'/'+bigroot+'_v.fits'
ifits=idir+'/'+bigroot+'_i.fits'
zfits=zdir+'/'+bigroot+'_z.fits'
rsd1fits=rsddir+'/'+bigroot+'_z_resid.fits'
fitfits=fitdir+'/'+bigroot+'_z_fit.fits'

; Read in files

print, ''
print, 'cand_plot: Reading input files:'
print, 'cand_plot:   ',bfits
bim = mrdfits(bfits, 0, bhead)
print, 'cand_plot:   ',vfits
vim = mrdfits(vfits, 0, vhead)
print, 'cand_plot:   ',ifits
iim = mrdfits(ifits, 0, ihead)
print, 'cand_plot:   ',zfits
zim = mrdfits(zfits, 0, zhead)
print, 'cand_plot:   ',rsd1fits
rsd1im = mrdfits(rsd1fits, 0, rsd1head)
print, 'cand_plot:   ',fitfits
fitim = mrdfits(fitfits, 0, fithead)
sz=size(zim,/dimen)
print, ''


; Set up scaled versions of the input files

imscal, bim, bimsc, bstat, lsig=1.5, hsig=hsig, nsig=3.0
imscal, vim, vimsc, vstat, lsig=1.5, hsig=hsig, nsig=3.0
imscal, iim, iimsc, istat, lsig=1.5, hsig=hsig, nsig=3.0
imscal, zim, zimsc, zstat, lsig=1.5, hsig=zsig, nsig=3.0
imscal, rsd1im, reimsc, restat, lsig=1.5, hsig=rhsig, nsig=3.0
fitsc = bytscl(fitim,min=zstat.min,max=zstat.max)

; Set up 3-plane array for inputs to true-color images

bvi=bytarr(3,sz[0],sz[1])
bvi[0,*,*] = iimsc
bvi[1,*,*] = vimsc
bvi[2,*,*] = bimsc

biz=bytarr(3,sz[0],sz[1])
biz[0,*,*] = zimsc
biz[1,*,*] = iimsc
biz[2,*,*] = bimsc

; bviz=bytarr(3,sz[0],sz[1])
; bviz[0,*,*] = zimsc
; bviz[1,*,*] = viimsc
; bviz[2,*,*] = bimsc

; Set up viewport window


subcellarray, [1,1,1,1],[1,1], newpan, newsubpan

; Get information from zband file's header

cd1 = fxpar(zhead, 'cd1_*', count=count1)
cd2 = fxpar(zhead, 'cd2_*', count=count2)

if (count1 eq 0) then begin
   pixscale=0.05 
endif else begin
    pixscale=sqrt(cd1[0]*cd1[0] + cd1[1]*cd1[1])*3600.0
endelse

; Set displayed image size and label locations

setup4.fullxsz = sz[0]*pixscale / 2.0
setup4.fullysz = sz[1]*pixscale / 2.0
print, ''
print, 'pixscale = ',pixscale
print, 'Native image size is ',setup4.fullxsz*2.,' arcsec x',$
   setup4.fullysz*2.,' arcsec'
if n_elements(imgsz) gt 0 then begin
  setup4.imgxsz = imgsz / 2.
  setup4.imgysz = imgsz / 2.
endif else begin
  setup4.imgxsz = setup4.fullxsz
  setup4.imgysz = setup4.fullysz
endelse
print, 'Displayed image size will be:',setup4.imgxsz*2., $
   ' arcsec x',setup4.imgysz*2.,' arcsec'

; Fill in other setup parameters

setup4.row = 1
setup4.labx = -0.9*setup4.imgxsz
setup4.laby = -0.8*setup4.imgysz
setup4.namex = 0.0
setup4.namey = 0.8*setup4.imgysz
setup4.name = bigroot
setup4.lab1 = 'BVi'
setup4.lab2 = 'z'
setup4.lab3 = 'Fit'
setup4.lab4 = 'Resid'

; Create setup for second row

setupb = setup4
setupb.row = 0
setupb.name = ''
setupb.lab1 = 'B'
setupb.lab2 = 'V'
setupb.lab3 = 'i'
setupb.lab4 = 'z'


; Plot images

print, ''
print, 'Plotting data for ',bigroot
print, ''
plot_4x, bvi, zimsc, fitsc, reimsc, setup4, newpan, newsubpan
plot_4x, bimsc, vimsc, iimsc, zimsc, setupb, newpan, newsubpan, /noerase

; Plot postscript if desired

if n_elements(psfile) gt 0 then begin
   print,'Saving image in postscript file '+psfile+'.ps.'
   ps_open,psfile,/color,/ps_fonts
;   device,xsize=6,ysize=3,/inches,/portrait,/times
   device,bits=8,/portrait,/times
   setup4.namecolor = 0
   loadct,0
   invct
   plot_4x, bvi, zimsc, fitsc, reimsc, setup4, newpan, newsubpan
   plot_4x, bimsc, vimsc, iimsc, zimsc, setupb, newpan, newsubpan, /noerase
   ps_close
endif


end 
