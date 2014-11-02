pro plot_resid, id, prefix=prefix, resid1=resid1, resid2=resid2, hsig=hsig, $
   rhsig=rhsig, secat=secat, _extra=_extra

;
; Procedure plot_resid
;
; Description: Makes a multi-panel plot containing true-color images of
;               a galaxy and residual images created by subtracting off
;               some form of fit to the galaxy light distribution.
;
; Inputs: id          (string)     source id, either id number only or
;                                   full source name
;         [prefix=]   (string)     optional source-name prefix, if id is
;                                   only an id number
;         [resid1=]   (string)     name of directory containing residuals
;                                   default is 'Resid'
;         [resid2=]   (string)     optional second directory of residuals,
;                                   used if side-by-side comparision is
;                                   desired.
;         [hsig=]     (float)      upper limit for display range in units of
;                                   the rms in the clipped image.  The
;                                   upper limit wil be hsig*rms + mean.
;                                   Default value is 30.0.
;         [rhsig=]    (float)      same as hsig, but for the residual image
;                                   only.  Default value is 10.0.
;         [secat=]    (struct)     structure containing SExtractor info.
;
; Revision history:
;  2003Feb26 Chris Fassnacht -- A modification of Lexi's eidol.pro.
;                                (first working version).
;  2003Feb27 Chris Fassnacht -- Changed from a two-panel to a 4-panel
;                                plot, with two copies (one full-sized,
;                                and one zoomed) each of a bvi true-color 
;                                plot and the residual plot
;  2003Feb27 Chris Fassnacht -- Added informational printout and some
;                                titles for plots.
;  2003Mar05 Chris Fassnacht -- Added option to plot a second set of
;                                residuals for side-by-side comparision
;  2003Mar18 Chris Fassnacht -- Changed to reflect new version of sigclip
;                               Moved image scaling into new imscal procedure
;                               Took out explicit creation of zoomed images --
;                                functionality now replaced with the 
;                                xrange and yrange parameters of plotimage.
;  2003Apr11 Chris Fassnacht -- Added hsig and rhsig as passed parameters
;                                rather than hard-wiring the values.
;  2003Jun22 Chris Fassnacht -- Added ability to plot an ellipse on the
;                                z-band image through the new optional
;                                secat parameter
;

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, 'syntax: plot_resid, id [,prefix=prefix, resid1=resid1, '
    print, '                    resid2=resid2, hsig=hsig, rhsig=rhsig, '
    print, '                    secat=secat, _extra=_extra]'
    print, ''
    print, " Default value for resid1 is 'Resid'"
    print, ' resid2 is only used if a side-by-side comparision of residuals'
    print, '  is desired.  There is no default setting.'
    print, ''
    return
endif

; Set directory names and defaults

bdir = 'b-band'
vdir = 'v-band'
idir = 'i-band'
zdir = 'z-band'
if not keyword_set(resid1) then resid1 = 'Resid'
if keyword_set(resid2) then tworesid=1B else tworesid=0B
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
rsd1fits=resid1+'/'+bigroot+'_z_resid.fits'
if(tworesid) then rsd2fits=resid2+'/'+bigroot+'_z_resid.fits'

; Read in files

print, ''
print, 'plot_resid: Reading input files:'
print, 'plot_resid:   ',bfits
bim = mrdfits(bfits, 0, bhead)
print, 'plot_resid:   ',vfits
vim = mrdfits(vfits, 0, vhead)
print, 'plot_resid:   ',ifits
iim = mrdfits(ifits, 0, ihead)
print, 'plot_resid:   ',zfits
zim = mrdfits(zfits, 0, zhead)
print, 'plot_resid:   ',rsd1fits
rsd1im = mrdfits(rsd1fits, 0, rsd1head)
sz=size(zim,/dimen)
if(tworesid) then begin
  print, 'plot_resid:   ',rsd2fits
  rsd2im = mrdfits(rsd2fits, 0, rsd1head)
  sz2 = size(rsd2im,/dimen)
endif
print, ''


; Set up scaled versions of the input files

imscal, bim, bimsc, bstat, lsig=1.5, hsig=hsig, nsig=3.0
imscal, vim, vimsc, vstat, lsig=1.5, hsig=hsig, nsig=3.0
imscal, iim, iimsc, istat, lsig=1.5, hsig=hsig, nsig=3.0
imscal, zim, zimsc, zstat, lsig=1.5, hsig=zsig, nsig=3.0
imscal, rsd1im, reimsc, restat, lsig=1.5, hsig=rhsig, nsig=3.0

if(tworesid) then begin
   re2imsc = bytscl(rsd2im,min=restat.min,max=restat.max)
endif

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

if(tworesid) then subcellarray, [1,1,1,1],[1,1], newpan, newsubpan else $
   subcellarray, [1,1,1],[1,1], newpan, newsubpan

; Set displayed image size

fullxsz = sz[0]*pixscale / 2.0
fullysz = sz[1]*pixscale / 2.0
if(tworesid) then begin
  re2xsz = sz2[0]*pixscale / 2.0
  re2ysz = sz2[1]*pixscale / 2.0
endif
zoomxsz = 2.0
zoomysz = 2.0
print, ''

; Plot images

print, ''
print, 'Plotting data for ',bigroot
print, ''

; Plot full-sized images

plotimage,bvi,/preserve,_extra=_extra,panel=newpan[0,1,*], $
   subpan=newsubpan[0,1,*],imgxr=[-1,1]*fullxsz, imgyr=[-1,1]*fullysz, $
   xrange=[-1,1]*fullxsz,yrange=[-1,1]*fullysz, $
   xtit='arcsec',ytit='arcsec',title=bigroot
xyouts,-0.9*fullxsz,0.8*fullysz,'BVI',/data

plotimage,zimsc,/preserve,_extra=_extra,panel=newpan[1,1,*], $
   subpan=newsubpan[1,1,*],imgxr=[-1,1]*fullxsz, imgyr=[-1,1]*fullysz, $
   xrange=[-1,1]*fullxsz,yrange=[-1,1]*fullysz, $
   xtit='arcsec',ytickn=replicate(' ',20),title="Full-sized",/noerase
xyouts,-0.9*fullxsz,0.8*fullysz,'z',/data

if keyword_set(secat) then $
   tvellipse, secat.a, secat.b, secat.x, secat.y, secat.pa, /data, color=255

plotimage,reimsc,/preserve,_extra=_extra,panel=newpan[2,1,*], $
   subpan=newsubpan[2,1,*],imgxr=[-1,1]*fullxsz, imgyr=[-1,1]*fullysz, $
   xrange=[-1,1]*fullxsz,yrange=[-1,1]*fullysz, $
   xtit='arcsec',ytickn=replicate(' ',20),title="Full-sized",/noerase
if(tworesid) then labstr1='Residuals: '+resid1 else labstr1='Residuals'
xyouts,-0.9*fullxsz,0.8*fullysz,labstr1,/data

if(tworesid) then begin
   if(sz2[0] gt sz[0]) then pltsz=sz else pltsz=sz2
   plotimage,re2imsc,/preserve,_extra=_extra,panel=newpan[3,1,*], $
      subpan=newsubpan[3,1,*],imgxr=[-1,1]*re2xsz,imgyr=[-1,1]*re2ysz, $
      xrange=[-1,1]*fullxsz,yrange=[-1,1]*fullysz, $
      xtit='arcsec',ytickn=replicate(' ',20),/noerase
   labstr2='Residuals: '+resid2
   xyouts,-0.9*fullxsz,0.8*fullysz,labstr2,/data
endif

; Plot zoomed images

plotimage,bvi,/preserve,_extra=_extra,panel=newpan[0,0,*], $
   subpan=newsubpan[0,0,*],imgxr=[-1,1]*fullxsz, imgyr=[-1,1]*fullysz, $
   xrange=[-1,1]*zoomxsz,yrange=[-1,1]*zoomysz, $
   xtit='arcsec',ytit='arcsec',/noerase,title='Zoom'
xyouts,-0.9*zoomxsz,0.8*zoomysz,'BVI',/data

plotimage,zimsc,/preserve,_extra=_extra,panel=newpan[1,0,*], $
   subpan=newsubpan[1,0,*],imgxr=[-1,1]*fullxsz, imgyr=[-1,1]*fullysz, $
   xrange=[-1,1]*zoomxsz,yrange=[-1,1]*zoomysz, $
   xtit='arcsec',ytickn=replicate(' ',20),/noerase
xyouts,-0.9*zoomxsz,0.8*zoomysz,'z',/data

plotimage,reimsc,/preserve,_extra=_extra,panel=newpan[2,0,*], $
   subpan=newsubpan[2,0,*],imgxr=[-1,1]*fullxsz, imgyr=[-1,1]*fullysz, $
   xrange=[-1,1]*zoomxsz,yrange=[-1,1]*zoomysz, $
   xtit='arcsec',ytickn=replicate(' ',20),/noerase
xyouts,-0.9*zoomxsz,0.8*zoomysz,labstr1,/data

if(tworesid) then begin
   plotimage,re2imsc,/preserve,_extra=_extra,panel=newpan[3,0,*], $
      subpan=newsubpan[3,0,*],imgxr=[-1,1]*re2xsz,imgyr=[-1,1]*re2ysz, $
      xrange=[-1,1]*zoomxsz,yrange=[-1,1]*zoomysz, $
      xtit='arcsec',ytickn=replicate(' ',20),/noerase
   xyouts,-0.9*zoomxsz,0.8*zoomysz,labstr2,/data
endif


end 
