pro goods_stamp, field, centx, centy, size, outroot
;
; Function goods_stamp
;
; Description: A kludgy hardwired procedure that, for a given x and y
;               calls poststamp for each of the 4 goods bands (b,v,i,z).
;               The fits file names are hardwired, with the only variable
;               being whether the data are from the northern or southern
;               field.
;
; Inputs: field      (string)      GOODS field (n or s)
;         centx      (float)       x position (in input image) to be
;                                   used as the center of the postage stamp
;         centy      (float)       y position (in input image) to be
;                                   used as the center of the postage stamp
;         size       (intarr)      output size.  if a non-square image is
;                                   desired, this should be passed as a
;                                   2-dim array
;         outroot    (string)      root name for output files
;
; Revision history:
;  2002May06 Chris Fassnacht -- First rough version.
;

; Check input format

if n_params() lt 5 then begin
    print, ''
    print, 'syntax: goods_stamp, field, centx, centy, size, outroot'
    print, ''
    print, " field is either 'n' or 's'"
    return
endif

; Set file names

mosdir = '/home/cdf/GOODS/Mosaics/'
bfits = mosdir+field+'1b_mosV1.5_scale1_drz.fits'
vfits = mosdir+field+'123v_mosV2.0_scale1_drz.fits'
ifits = mosdir+field+'123i_mosV2.0_scale1_drz.fits'
zfits = mosdir+field+'123z_mosV1.8_scale1_drz.fits'

; Set sizes

if n_elements(size) eq 2 then begin
   nx = size[0]
   ny = size[1]
endif else begin
   nx = size[0]
   ny = size[0]
endelse

; Call poststamp for each of the four bands

bim = poststamp(bfits,centx,centy,nx,ny,outhead=bhead)
vim = poststamp(vfits,centx,centy,nx,ny,outhead=vhead)
iim = poststamp(ifits,centx,centy,nx,ny,outhead=ihead)
zim = poststamp(zfits,centx,centy,nx,ny,outhead=zhead)

; Write output fits files for each of the four bands

print, ''
print, 'goods_stamp: Writing out fits files for ',outroot
print, 'goods_stamp: --------------------------------------------------------'

bfile = outroot+'_b.fits'
vfile = outroot+'_v.fits'
ifile = outroot+'_i.fits'
zfile = outroot+'_z.fits'

print, 'goods_stamp:   Writing ',bfile
writefits, bfile, bim, bhead
print, 'goods_stamp:   Writing ',vfile
writefits, vfile, vim, vhead
print, 'goods_stamp:   Writing ',ifile
writefits, ifile, iim, ihead
print, 'goods_stamp:   Writing ',zfile
writefits, zfile, zim, zhead

print, 'goods_stamp: Finished for object ',outroot

end
