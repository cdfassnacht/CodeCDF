pro goods_bvi, id, prefix=prefix, bdir=bdir, vdir=vdir, idir=idir, $
   imgsz=imgsz, hsig=hsig, psfile=psfile, _extra=_extra

;
; Procedure goods_bvi
;
; Description: Given a source name, plots a BVi true-color image 
;  showing the source.
;
; Inputs: id          (string)     root part of name for all input files.
;                                   May be modified by the optional prefix
;                                   parameter
;         [prefix=]   (string)     optional source-name prefix, if candlist
;                                   contains only id numbers
;         [bdir=]     (string)     name of directory containing b-band files
;                                   default is 'b-band'
;         [vdir=]     (string)     name of directory containing v-band files
;                                   default is 'v-band'
;         [idir=]     (string)     name of directory containing i-band files
;                                   default is 'i-band'
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
;
; Revision history:
;  2003Jul31 Chris Fassnacht -- A revision of goods_bviz
;

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, $
      'syntax: goods_bvi, root, [prefix=prefix,bdir=bdir,vdir=vdir,'
    print, '       idir=idir,imgsz=imgsz,hsig=hsig,psfile=psfile]'
    print, ''
    print, " Default value for bdir is 'b-band'"
    print, " Default value for vdir is 'v-band'"
    print, " Default value for idir is 'i-band'"
    print, ' hsig sets the upper limit (in units of the RMS) for all '
    print, '   plots.  Default value is 30.0'
    print, ''
    return
endif

; Create setup structure

def_pltsetup, setup

; Set directory names and defaults

if not keyword_set(bdir) then bdir = 'b-band'
if not keyword_set(vdir) then vdir = 'v-band'
if not keyword_set(idir) then idir = 'i-band'
pixscale = 0.05
if n_elements(hsig) eq 0 then hsig = 30.0
zsig = hsig

; Set filenames

if n_elements(prefix) gt 0 then bigroot = prefix+'_'+id else bigroot = id
bfits=bdir+'/'+bigroot+'_b.fits'
vfits=vdir+'/'+bigroot+'_v.fits'
ifits=idir+'/'+bigroot+'_i.fits'

; Read in files

print, ''
print, 'goods_bvi: Reading input files:'
print, 'goods_bvi:   ',bfits
bim = mrdfits(bfits, 0, bhead)
print, 'goods_bvi:   ',vfits
vim = mrdfits(vfits, 0, vhead)
print, 'goods_bvi:   ',ifits
iim = mrdfits(ifits, 0, ihead)
sz=size(vim,/dimen)
print, ''


; Set up scaled versions of the input files

imscal, bim, bimsc, bstat, lsig=1.5, hsig=hsig, nsig=3.0
imscal, vim, vimsc, vstat, lsig=1.5, hsig=hsig, nsig=3.0
imscal, iim, iimsc, istat, lsig=1.5, hsig=hsig, nsig=3.0

; Set up 3-plane array for inputs to true-color images

bvi=bytarr(3,sz[0],sz[1])
bvi[0,*,*] = iimsc
bvi[1,*,*] = vimsc
bvi[2,*,*] = bimsc

; Set up viewport window


subcellarray, [1],[1], newpan, newsubpan

; Get information from zband file's header

cd1 = fxpar(vhead, 'cd1_*', count=count1)
cd2 = fxpar(vhead, 'cd2_*', count=count2)

if (count1 eq 0) then begin
   pixscale=0.05 
endif else begin
    pixscale=sqrt(cd1[0]*cd1[0] + cd1[1]*cd1[1])*3600.0
endelse

; Set displayed image size and label locations

setup.fullxsz = sz[0]*pixscale / 2.0
setup.fullysz = sz[1]*pixscale / 2.0
print, ''
print, 'pixscale = ',pixscale
print, 'Native image size is ',setup.fullxsz*2.,' arcsec x',$
   setup.fullysz*2.,' arcsec'
if n_elements(imgsz) gt 0 then begin
  setup.imgxsz = imgsz / 2.
  setup.imgysz = imgsz / 2.
endif else begin
  setup.imgxsz = setup.fullxsz
  setup.imgysz = setup.fullysz
endelse
print, 'Displayed image size will be:',setup.imgxsz*2., $
   ' arcsec x',setup.imgysz*2.,' arcsec'

; Fill in other setup parameters

setup.row = 0
setup.labx = -0.9*setup.imgxsz
setup.laby = -0.8*setup.imgysz
setup.namex = 0.0
setup.namey = 0.8*setup.imgysz
setup.name = bigroot
setup.lab1 = 'BVi'

; Plot images

print, ''
print, 'Plotting data for ',bigroot
print, ''
setup.column = 0
plot_panel, bvi, 0, 0, setup, newpan, newsubpan

; Plot postscript if desired

if n_elements(psfile) gt 0 then begin
   print,'Saving image in postscript file '+psfile+'.ps.'
   ps_open,psfile,/color,/ps_fonts,/encap
;   device,xsize=6,ysize=3,/inches,/portrait,/times
   device,bits=8,/portrait,/times
   setup.namecolor = 0
   loadct,0
   invct
   setup.lab1 = 'BVi'
   plot_panel, bvi, 0, 0, setup, newpan, newsubpan
   ps_close
endif


end 
