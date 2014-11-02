pro goods_stamp_ra, field, ra, dec, size, outroot
;
; Function goods_stamp_ra
;
; Description: A kludgy hardwired procedure that, for a given ra and dec
;               calls stamp_ra for each of the 4 goods 
;               bands (b,v,i,z). The fits file names are hardwired, with the 
;               only variable being whether the data are from the northern or 
;               southern field.
;
; Inputs: field      (string)      GOODS field (n or s)
;         ra         (double)      RA (in input image) to be
;                                   used as the center of the postage stamp
;         dec        (double)      Dec (in input image) to be
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
    print, 'syntax: goods_stamp_ra, field, ra, dec, size, outroot'
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

; Call stamp_ra for each of the four bands

bim = stamp_radec(bfits,ra,dec,size[0],outhead=bhead)
vim = stamp_radec(vfits,ra,dec,size[0],outhead=vhead)
iim = stamp_radec(ifits,ra,dec,size[0],outhead=ihead)
zim = stamp_radec(zfits,ra,dec,size[0],outhead=zhead)

; Write output fits files for each of the four bands

print, ''
print, 'goods_stamp_ra: Writing out fits files for ',outroot
print, 'goods_stamp_ra: --------------------------------------------------------'

bfile = outroot+'_b.fits'
vfile = outroot+'_v.fits'
ifile = outroot+'_i.fits'
zfile = outroot+'_z.fits'

print, 'goods_stamp_ra:   Writing ',bfile
writefits, bfile, bim, bhead
print, 'goods_stamp_ra:   Writing ',vfile
writefits, vfile, vim, vhead
print, 'goods_stamp_ra:   Writing ',ifile
writefits, ifile, iim, ihead
print, 'goods_stamp_ra:   Writing ',zfile
writefits, zfile, zim, zhead

print, 'goods_stamp_ra: Finished for object ',outroot

end
