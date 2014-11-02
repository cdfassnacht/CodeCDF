pro mydisp, fitsfile, indir=indir, lsig=lsig, hsig=hsig, psfile=psfile, $
   jpfile=jpfile, imgsz=imgsz, _extra=_extra

;
; Procedure mydisp
;
; Description: Given an input fits file, displays an image which can be 
;               exported into either a postscript file, a jpeg file, or both.
;
; Inputs: fitsfile    (string)     name of input fits file
;         [lsig=]     (float)      number of sigma below mean to set lower
;                                   display limit (default = 1.0)
;         [hsig=]     (float)      number of sigma above mean to set upper
;                                   display limit (default = 30.0)
;         [psfile=]   (string)     name for output postscript file (if not set,
;                                   no file is created)
;         [jpfile=]   (string)     name for output jpeg file (if not set,
;                                   no file is created)
;         [imgsz=]    (float)      image size in arcsec.  If chosen, will
;                                   produce an output image of imgsz x imgsz
;
; Revision history:
;  2003Mar31 Chris Fassnacht -- First working version.
;  2003Apr02 Chris Fassnacht -- Fixed several small bugs.
;

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, $
      'syntax: mydisp, fitsfile [indir=indir,lsig=lsig,hsig=hsig,'
    print, '  imgsz=imgsz,psfile=psfile,jpfile=jpfile]'
    print, ''
    print, " Default value for indir is '.'"
    print, ' Default value for lsig is 1.0'
    print, ' Default value for hsig is 30.0'
    print, ''
    return
endif

; Set directory names and defaults

bdir = 'b-band'
vdir = 'v-band'
idir = 'i-band'
zdir = 'z-band'
if not keyword_set(indir) then indir='.'
if not keyword_set(lsig) then lsig=1.0
if not keyword_set(hsig) then hsig=30.0
pixscale = 0.05

; Read in file

infits = indir+'/'+fitsfile
print, ''
print, 'mydisp: Reading input file:'
print, 'mydisp:   ',infits
img = mrdfits(infits, 0, inhead)
sz=size(img,/dimen)
print, ''

; Get information from first file's header

cd1 = fxpar(inhead, 'cd1_*', count=count1)
cd2 = fxpar(inhead, 'cd2_*', count=count2)

if (count1 eq 0) then begin
   pixscale=0.05 
endif else begin
    pixscale=sqrt(cd1[0]*cd1[0] + cd1[1]*cd1[1])*3600.0
endelse
print, 'pixscale = ',pixscale
print, ''

; Set up scaled versions of the input files

imscal, img, scalimg, imgstat, nsig=3.0, lsig=lsig, hsig=hsig

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
print, 'Plotting image -- display limits are:'
print, '   lower = clipped mean - ',lsig,' sigma.'
print, '   upper = clipped mean + ',hsig,' sigma.'
print, ''

plotimage,scalimg,/preserve,_extra=_extra,panel=newpan[0,0,*], $
   subpan=newsubpan[0,0,*],imgxr=[-1,1]*angxsz,imgyr=[-1,1]*angysz, $
   xrange=[-1,1]*imgxsz,yrange=[-1,1]*imgysz, $
   xtit='arcsec',ytit='arcsec',title=fitsfile

; Plot to a postscript file

if n_elements(psfile) gt 0 then begin
   print,' Saving image in postscript file '+psfile+'.ps.'
   ps_open,psfile,/color,/ps_fonts
   device,xsize=6,ysize=6,/inches,/portrait,/times
   plotimage,scalimg,/preserve,_extra=_extra,panel=newpan[0,0,*], $
      subpan=newsubpan[0,0,*],imgxr=[-1,1]*angxsz,imgyr=[-1,1]*angysz, $
      xrange=[-1,1]*imgxsz,yrange=[-1,1]*imgysz, $
      xtit='arcsec',ytit='arcsec',title=fitsfile
   ps_close
endif

; Plot to a jpeg file (no border or labels)

if n_elements(jpfile) gt 0 then begin
   write_jpeg,jpfile,scalimg,true=1
   print,' Image saved in JPEG file '+jpfile+'.'
endif

end 
