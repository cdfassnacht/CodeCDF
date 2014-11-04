; master_phot_bells.pro
;
; master ACS photometry of SLACS and BELLS targets.
;
; Goals are to (1) reproduce SLACS deVauc photometry,
; then (2) apply to BELLS, then (3) try flat-weight fitting.
;
; does b-spline and deVauc/Sersic fitting in a variety
; of different ways.
;
; abolton@cfa 2007jun
; updated bolton@utah 2011jun
;

;metadir = '/Volumes/raid_b/ACS/metadata/'
;rawdir = '/Volumes/raid_b/ACS/rawdata/'
;maskdir = '/Volumes/raid_b/ACS/maskdata/'
;metafile = 'meta_comb_03.fits'
;ms = mrdfits(metadir + metafile, 1)
;nc = n_elements(ms)
;files = rawdir + string(ms.proposid, format='(i5.5)') $
; + '/' + strtrim(ms.root_new,2) + '.fits'

function master_phot_bells, imgfile=imgfile, maskfile=maskfile, bizfile=bizfile

;imgfile = '/Volumes/bloo/slacs/ACS/rawdata/10886/SLACSJ1706+3304_54102_814_4.fits'
;maskfile = '/Volumes/bloo/slacs/ACS/maskdata/SLACSJ1706+3304_mask.fits'
;;;;bfile = '/Volumes/bloo/slacs/ACS/rawdata/10886/SLACSJ1706+3304_54102_814_4_bsmods.fits'
;bizfile = '/Volumes/bloo/slacs/ACS/bizzle/02/YES/SLACSJ1706+3304_54102_814_4_biz02.fits'


; PSF half-width:
phw = 18
; pixel size:
dpix = 0.05

; Breakpoints for b-spline photometry:
rbkpt = dpix*[1., 2., 4., 6., 8., 12., 16., 22., 30., $
 50., 60., 100., 150., 200., 300., 400., 500., 600., 800.]

;for i = 0, nc-1 do begin
;for i = 1,1 do begin

; Read in files:
img = float(mrdfits(imgfile,1))
err = float(mrdfits(imgfile,2))
nim = mrdfits(imgfile,3)
psf = float(mrdfits(imgfile,4))
nx = (size(img))[1]
ny = (size(img))[2]

; Trim PSF to size:
    nx_psf = (size(psf))[1]
    junk = max(psf, wmax)
    pxc = wmax mod nx_psf
    pyc = wmax / nx_psf
    psf = psf[pxc-phw:pxc+phw,pyc-phw:pyc+phw]
    psf = psf / total(psf)

; Make inverse variance and get masks:
ivar = (nim gt 0) * (err gt 0.) / (err^2 + (nim le 0) + (err le 0.))
masks = mrdfits(maskfile, 1)
; Kluge to get something that works with Joel's format:
if (n_elements(masks.fmask) eq n_elements(img)) then $
   masks = {fmask: mrdfits(bizfile,6), jmask: mrdfits(bizfile,5), $
            xc: ((size(img))[1] - 1) / 2, yc: ((size(img))[2] - 1) / 2}

    xhw = (size(masks.fmask))[1] / 2
    yhw = (size(masks.fmask))[2] / 2
    xc = masks.xc
    yc = masks.yc

ivar[xc-xhw:xc+xhw,yc-yhw:yc+yhw] = $
 ivar[xc-xhw:xc+xhw,yc-yhw:yc+yhw] * masks.fmask * masks.jmask
ymask = 0B * byte(ivar)
ymask[xc-xhw:xc+xhw,yc-yhw:yc+yhw] = 1B
pmask = replicate(0B, nx, ny)
pmask[phw:nx-(phw+1),phw:ny-(phw+1)] = 1B
ivar = ivar * pmask
pmask = 0

; Get the right pixel statistics:
djs_iterstat, img * sqrt(ivar), sigma=errsig
ivar = ivar / errsig^2

; Make RA and Dec images:
aimg = -dpix * (findgen(nx) - xc) # replicate(1.0, ny)
dimg = dpix * replicate(1.0, nx) # (findgen(ny) - yc)

; Get centering parameters from the b-spline fits:
;     bfile = rawdir + string(ms[i].proposid, format='(i5.5)') $
;      + '/' + strtrim(ms[i].root_new,2) + '_bsmods.fits'
;     bstruc = mrdfits(bfile,1)
     bhdr = headfits(bizfile, exten=7)
     b_xc = sxpar(bhdr, 'NAXIS1') / 2
     b_yc = sxpar(bhdr, 'NAXIS2') / 2
     b_offset_dec = (sxpar(bhdr, 'YCENTER') - b_xc) * dpix
     b_offset_ra = - (sxpar(bhdr, 'XCENTER') - b_yc) * dpix
     b_axisrat = sxpar(bhdr, 'AXISRAT')
     b_axisangl = sxpar(bhdr, 'AXISANGL')

; Starting deVaucouleurs parameters:
dev_amp_zero = 10.
dev_rad_zero = 1.
dspars = float([dev_amp_zero, dev_rad_zero, b_offset_dec, $
                b_offset_ra, b_axisrat, b_axisangl])

; Do the initial deVaucouleurs fit:
;  dspars = [10., 1., bstruc.bspar]
  devpars = mpfit2dfun('asb_devauc_multi', dimg, aimg, $
   img, 0.*img, dspars, weights=ivar, functargs={psf: psf, rmin: 0.1*dpix}, $
   bestnorm=bestnorm, perror=perror, niter=niter, yfit=yfit, covar=covar)

; Make a mask for any outlying junk:
sigthresh = 4.0
newmask = ((img-yfit) * sqrt(ivar)) ge sigthresh
newmask = newmask * (ymask eq 0)
ivar = ivar * (newmask eq 0)

; Do the next deVaucouleurs fit:
;  devpars = mpfit2dfun('asb_devauc_multi', dimg, aimg, $
;   img, 0.*img, devpars, weights=ivar, functargs={psf: psf, rmin: 0.1*dpix}, $
;   bestnorm=bestnorm, perror=perror, niter=niter, yfit=yfit, covar=covar)

;;; If you dare...
; Subsample the grid and do one last MPFIT:
subs = 5
xbase = 0.5 * (2.*findgen(nx*subs) + 1. - float(subs))/float(subs)
ybase = 0.5 * (2.*findgen(ny*subs) + 1. - float(subs))/float(subs)
abig = bilinear(aimg, xbase, ybase)
dbig = bilinear(dimg, xbase, ybase)
xbase = 0
ybase = 0

  devpars = mpfit2dfun('asb_devauc_multi', dbig, abig, $
   img, 0.*img, devpars, weights=ivar, $
   functargs={psf: psf, rmin: 0.1*dpix, subs: subs}, $
   bestnorm=bestnorm, perror=perror, niter=niter, yfit=yfit, covar=covar)

; Keep it positive!
devpars[1] = abs(devpars[1])

; A fit with flat weighting:
med_ivar = median(ivar[where(ivar gt 0.)])
ivar_flat = med_ivar * (ivar gt 0.)

  devpars_flat = mpfit2dfun('asb_devauc_multi', dbig, abig, $
   img, 0.*img, devpars, weights=ivar_flat, $
   functargs={psf: psf, rmin: 0.1*dpix, subs: subs}, $
   bestnorm=bestnorm_flat, perror=perror_flat, niter=niter_flat, $
   yfit=yfit_flat, covar=covar_flat)

ndof = long(total(ivar gt 0.)) - n_elements(devpars)

; Now for the (ugh!) Sersic models:
serpars = mpfit2dfun('asb_sersic_multi', dimg, aimg, img, 0.*img, $
 [devpars, 4.], weights=ivar, functargs={psf: psf, rmin: 0.1*dpix}, $
 bestnorm=sersic_bestnorm, perror=sersic_perror, niter=sersic_niter, $
 yfit=sersic_fit, covar=sersic_covar)

  serpars = mpfit2dfun('asb_sersic_multi', dbig, abig, $
   img, 0.*img, serpars, weights=ivar, $
   functargs={psf: psf, rmin: 0.1*dpix, subs: subs}, $
   bestnorm=sersic_bestnorm, perror=sersic_perror, niter=sersic_niter, $
   yfit=sersic_fit, covar=sersic_covar)

; ... and with flat weighting:
  serpars_flat = mpfit2dfun('asb_sersic_multi', dbig, abig, $
   img, 0.*img, serpars, weights=ivar_flat, $
   functargs={psf: psf, rmin: 0.1*dpix, subs: subs})

; Turn this all into proper band photometry:

; Get the deVaucouleurs fit parameters and the photometric
; header keyword parameters:

; Get the band keywords from the header:
  hdr = headfits(imgfile, exten=1)
  photflam = sxpar(hdr, 'PHOTFLAM')
  photplam = sxpar(hdr, 'PHOTPLAM')

;  modstruc = mrdfits(modfiles[i],1)
;  devpars[*,i] = modstruc.pars & $
;  deverrs[*,i] = modstruc.perror & $
;endfor

; Convert deVaucouleurs model amplitude to counts, and to AB magnitude:
k_dev = 7.66925001
cps_dev = ((2. * !pi * devpars[0] * devpars[1]^2) $
   * 4. * gamma(8.) / k_dev^8.) / (dpix^2)
cps_dev_flat = ((2. * !pi * devpars_flat[0] * devpars_flat[1]^2) $
   * 4. * gamma(8.) / k_dev^8.) / (dpix^2)
abmag_zpt = -2.5 * alog10(photflam) - 21.10 - 5. * alog10(photplam) + 18.6921
ab_dev = -2.5 * alog10(cps_dev) + abmag_zpt
ab_dev_flat = -2.5 * alog10(cps_dev_flat) + abmag_zpt

; Get the rest of the deVaucouleurs parameters straightened out:
cps_deverr = cps_dev * perror[0] / devpars[0]
acs_rdev = devpars[1]
acs_rdeverr = perror[1]
acs_qdev = devpars[4]
acs_qdeverr = perror[4]
acs_padev = devpars[5]
acs_padeverr = perror[5]
acs_nser = serpars[6]
acs_nsererr = sersic_perror[6]

; The flat-weight parameters:
acs_rdev_flat = devpars_flat[1]
acs_qdev_flat = devpars_flat[4]
acs_padev_flat = devpars_flat[5]
acs_nser_flat = serpars_flat[6]

if (acs_qdev gt 1.) then begin & $
  acs_qdev = 1./acs_qdev & $
  acs_qdeverr = acs_qdeverr * acs_qdev^2 & $
  acs_padev = acs_padev + 90. & $
endif

if (acs_qdev_flat gt 1.) then begin & $
  acs_qdev_flat = 1./acs_qdev_flat & $
;  acs_qdeverr = acs_qdeverr * acs_qdev^2 & $
  acs_padev_flat = acs_padev_flat + 90. & $
endif

; Normalize PA range between 0 and 180:

; regular weights:
pa = acs_padev
old_pa = pa
pa = pa mod 180.
ptest = pa ge 0.
pa = pa * ptest + (180. - abs(pa)) * (1B - ptest)
acs_padev = pa

; flat weights:
pa = acs_padev_flat
old_pa = pa
pa = pa mod 180.
ptest = pa ge 0.
pa = pa * ptest + (180. - abs(pa)) * (1B - ptest)
acs_padev_flat = pa

; Do the b-spline photometry:
bfit = bspline_ellip(dspars[2:*], x=dimg, y=aimg, invvar=ivar, $
 data=img, psf=psf, ntheta=0, rbkpt=rbkpt, sset=sset)

; Store the bspline photometry of interest.
; Really we're just evaluating the ellipsoidal
; b-spline surface brightness on a baseline.
; Anything else we want to do can be done later.

bbase = (findgen(800) + 1.) / 100.
bvalu = bspline_valu(bbase, sset)

newstruc = $
 {photflam: 0.d0, photplam: 0.d0, abmag_zpt: 0.d0, abmag_dev: 0., $
 cps_dev: 0., cps_deverr: 0., acs_rdev: 0., acs_rdeverr: 0., acs_qdev: 0., acs_qdeverr: 0., $
 acs_padev: 0., acs_padeverr: 0., acs_nser: 0., acs_nsererr: 0., $
 abmag_dev_flat: 0., cps_dev_flat: 0., acs_rdev_flat: 0., acs_qdev_flat: 0., $
 acs_padev_flat: 0., acs_nser_flat: 0., $
 bsr_arcs: fltarr(n_elements(bbase)), bsflux: fltarr(n_elements(bbase))}

newstruc.photflam = photflam
newstruc.photplam = photplam
newstruc.abmag_zpt = abmag_zpt
newstruc.abmag_dev = ab_dev
newstruc.abmag_dev_flat = ab_dev_flat
newstruc.cps_dev = cps_dev
newstruc.cps_dev_flat = cps_dev_flat
newstruc.cps_deverr = cps_deverr
newstruc.acs_rdev = acs_rdev
newstruc.acs_rdev_flat = acs_rdev_flat
newstruc.acs_rdeverr = acs_rdeverr
newstruc.acs_qdev = acs_qdev
newstruc.acs_qdev_flat = acs_qdev_flat
newstruc.acs_qdeverr = acs_qdeverr
newstruc.acs_padev = acs_padev
newstruc.acs_padev_flat = acs_padev_flat
newstruc.acs_padeverr = acs_padeverr
newstruc.acs_nser = acs_nser
newstruc.acs_nser_flat = acs_nser_flat
newstruc.acs_nsererr = acs_nsererr
newstruc.bsr_arcs = bbase
newstruc.bsflux = bvalu

return, newstruc
end


;catfile = '/Volumes/bloo/slacs/catalog/slacs_bintable.fits'
;modfile = '/Volumes/bloo/slacs/ACS/rawdata/10886/SLACSJ1706+3304_54102_814_4_photmod.fits'

