pro goodstrue, bfits, gfits, rfits, lsig=lsig, hsig=hsig, psfile=psfile, $
   jpfile=jpfile, tffile=tffile, imgsz=imgsz

;
; Procedure goodstrue
;
; Description: Given three input fits files taken in three bands, makes
;               a true-color image which can be exported into either a
;               postscript file, a jpeg file, or both.
;
; Inputs: bfits       (string)     name of fits file to be the "b" image
;         gfits       (string)     name of fits file to be the "g" image
;         rfits       (string)     name of fits file to be the "r" image
;         [lsig=]     (float)      number of sigma below mean to set lower
;                                   display limit (default = 1.0)
;         [hsig=]     (float)      number of sigma above mean to set upper
;                                   display limit (default = 15.0)
;         [psfile=]   (string)     name for output postscript file (if not set,
;                                   no file is created)
;         [jpfile=]   (string)     name for output jpeg file (if not set,
;                                   no file is created)
;         [tffile=]   (string)     name for output tiff file (if not set,
;                                   no file is created)
;         [imgsz=]    (float)      image size in arcsec.  If chosen, will
;                                   produce an output image of imgsz x imgsz
;
; Revision history:
;  2003Mar09 Chris Fassnacht -- First working version.
;  2003Mar18 Chris Fassnacht -- Added option for postscript file.
;                               Also added optional imgsz parameter to set
;                                output image size in arcsec
;                               Changed to reflect new version of sigclip.
;                               Moved image scaling into new imscal procedure
;

; Check input format

if n_params() lt 2 then begin
    print, ''
    print, $
      'syntax: goodstrue, bfits,gfits,rfits[,lsig=lsig,hsig=hsig,'
    print, '  imgsz=imgsz,psfile=psfile,jpfile=jpfile,tffile=tffile]'
    print, ''
    print, ' Default value for lsig is 1.0'
    print, ' Default value for hsig is 15.0'
    print, ''
    return
endif

; Set directory names and defaults

bdir = 'b-band'
vdir = 'v-band'
idir = 'i-band'
zdir = 'z-band'
if not keyword_set(lsig) then lsig=1.0
if not keyword_set(hsig) then hsig=15.0

; Read in files and then scale them
;
;  **NB: the "temporary" command is a tricky way of freeing memory

print, ''
print, 'goodstrue: Reading input files:'
print, 'goodstrue:   ',bfits
bim = mrdfits(bfits, 0, bhead)
sz=size(bim,/dimen)
bimsc = bytscl(temporary(bim),min=-0.0005,max=0.0194)
print, 'goodstrue:   ',gfits
gim = mrdfits(gfits, 0, ghead)
gimsc = bytscl(temporary(gim),min=-0.0008,max=0.0364)
print, 'goodstrue:   ',rfits
rim = mrdfits(rfits, 0, rhead)
rimsc = bytscl(temporary(rim),min=-0.0007,max=0.0334)
print, ''

; Set up scaled versions of the input files

;imscal, rim, rimsc, rstat, hsig=hsig, nsig=3.0
;imscal, gim, gimsc, gstat, hsig=hsig, nsig=3.0
;imscal, bim, bimsc, bstat, hsig=hsig, nsig=3.0

; Set up 3-plane array for inputs to true-color images

rgb=bytarr(3,sz[0],sz[1])
rgb[0,*,*] = rimsc
rgb[1,*,*] = gimsc
rgb[2,*,*] = bimsc

; Get information from first file's header

cd1 = fxpar(rhead, 'cd1_*', count=count1)
cd2 = fxpar(rhead, 'cd2_*', count=count2)

if (count1 lt 2) then begin
   pixscale=0.05 
endif else begin
    pixscale=sqrt(cd1[0]*cd1[0] + cd1[1]*cd1[1])*3600.0
endelse
print, 'pixscale = ',pixscale

; Set displayed image size

angxsz = sz[0]*pixscale/2.
angysz = sz[1]*pixscale/2.
print, ''
print, 'Native image size is ',angxsz*2.,' arcsec x',angysz*2.,' arcsec'
if n_elements(imgsz) gt 0 then begin
   imgxsz = imgsz / 2.
   imgysz = imgsz / 2.
endif else begin
   imgxsz = angxsz
   imgysz = angysz
endelse
print, 'Displayed image size will be:',imgxsz*2.,' arcsec x',imgysz*2.,' arcsec'

; Set up viewport window

subcellarray, [1],[1], newpan, newsubpan

; Plot images

print, ''
print, 'Plotting goodstrue image
print, ''

plotimage,rgb,/preserve,_extra=_extra,panel=newpan[0,0,*], $
   subpan=newsubpan[0,0,*],imgxr=[-1,1]*angxsz,imgyr=[-1,1]*angysz, $
   xrange=[-1,1]*imgxsz,yrange=[-1,1]*imgysz, $
   xtit='arcsec',ytit='arcsec',title=bigroot
xyouts,-0.9*sz[0]*pixscale/2.,0.8*sz[1]*pixscale/2.,'BVI',/data, $
   charsize=3.0, font=1

; Plot to a postscript file

if n_elements(psfile) gt 0 then begin
   print,' Saving image in postscript file '+psfile+'.ps.'
   ps_open,psfile,/color,/ps_fonts
   device,xsize=6,ysize=6,/inches,/portrait,/times
   plotimage,rgb,/preserve,_extra=_extra,panel=newpan[0,0,*], $
      subpan=newsubpan[0,0,*],imgxr=[-1,1]*angxsz,imgyr=[-1,1]*angysz, $
      xrange=[-1,1]*imgxsz,yrange=[-1,1]*imgysz, $
      xtit='arcsec',ytit='arcsec',title=bigroot
   xyouts,-0.9*sz[0]*pixscale/2.,0.8*sz[1]*pixscale/2.,'BVI',/data, $
      charsize=3.0, font=0, color=255
   ps_close
endif

; Plot to a jpeg file (no border or labels)

if n_elements(jpfile) gt 0 then begin
   write_jpeg,jpfile,rgb,true=1
   print,' Image saved in JPEG file '+jpfile+'.'
endif

; Plot to a tiff file (no border or labels)

if n_elements(tffile) gt 0 then begin
   write_tiff,tffile,rgb
   print,' Image saved in TIFF file '+tffile+'.'
endif

end 
