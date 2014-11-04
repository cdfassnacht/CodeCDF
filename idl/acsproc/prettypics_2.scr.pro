;
; prettypics_2.scr.pro
;
; Script to generate even prettier pictures of the
; lenses and their models.
;
; abolton@cfa 2007jun
;

metadir = getenv('SLACS_DATAROOT')+'ACS/metadata/'
metafile = 'spAll_ACS_04b.fits'
bizdir = getenv('SLACS_DATAROOT')+'ACS/bizzle/02/'

sp = mrdfits(metadir+metafile,1)
sp = sp[where(sp.asbgrade gt 0, nsp)]


; Find the biz and model files associated with these:
bfiles = strarr(nsp)
mfiles = strarr(nsp)
cfiles = strarr(nsp)
for i = 0L, nsp-1 do begin & $
  ff = file_search(bizdir+'*/'+sp[i].root_new+'_biz02.fits') & $
  nf = n_elements(ff) & $
  if (nf ne 1) then print, 'UH-OH: ', i & $
  bfiles[i] = ff[0] & $
  ff = file_search(bizdir+'*/'+sp[i].root_new+'_sie02.fits') & $
  nf = n_elements(ff) & $
  if (nf ne 1) then print, 'UH-OH: ', i & $
  mfiles[i] = ff[0] & $
  ff = file_search(bizdir+'*/'+sp[i].root_new+'_cml02.fits') & $
  nf = n_elements(ff) & $
  if (nf ne 1) then print, 'UH-OH: ', i & $
  cfiles[i] = ff[0] & $
endfor

; Get the necessary images, etc:
hw = 50
imgstack = fltarr(2*hw+1,2*hw+1,nsp)
bspstack = fltarr(2*hw+1,2*hw+1,nsp)
modstack = fltarr(2*hw+1,2*hw+1,nsp)  ; for SIE models
cmodstack = fltarr(2*hw+1,2*hw+1,nsp) ; for CML models
srcstack = fltarr(2*hw+1,2*hw+1,nsp)  ; for SIE source plane
csrcstack = fltarr(2*hw+1,2*hw+1,nsp) ; for CML source plane
invstack = fltarr(2*hw+1,2*hw+1,nsp)
mskstack = bytarr(2*hw+1,2*hw+1,nsp)  ; junk-mask stack -- not used?

magnif = 2  ; magnification factor for source plane
for i = 0L, nsp-1 do begin & $
  im = mrdfits(bfiles[i],0) & $
  iv = mrdfits(bfiles[i],1) & $
  jm = mrdfits(bfiles[i],5) & $
  mo = mrdfits(bfiles[i],7) & $
  nx = (size(im))[1] & $
  ny = (size(im))[2] & $
  xc = nx/2 & $
  yc = ny/2 & $
  imgstack[*,*,i] = im[xc-hw:xc+hw,yc-hw:yc+hw] & $
  bspstack[*,*,i] = mo[xc-hw:xc+hw,yc-hw:yc+hw] & $
  mskstack[*,*,i] = jm[xc-hw:xc+hw,yc-hw:yc+hw] & $
  invstack[*,*,i] = iv[xc-hw:xc+hw,yc-hw:yc+hw] & $

  lmod = mrdfits(mfiles[i],1) & $
  lft = mrdfits(mfiles[i],2) & $
  nx2 = (size(lmod.modim))[1] & $
  ny2 = (size(lmod.modim))[2] & $
  xc2 = nx2/2 & $
  yc2 = ny2/2 & $
  modstack[*,*,i] = lmod.modim[xc2-hw:xc2+hw,yc2-hw:yc2+hw] & $
  testim = sie_lens_gauss([0.,lmod.fpar[1:*]], x=lft.x/float(magnif), $
   y=lft.y/float(magnif), subs=lft.subs, $
   psf=lft.psf, $
   rmin=lft.rmin) & $
  srcstack[*,*,i] = testim[xc2-hw:xc2+hw,yc2-hw:yc2+hw] & $

  clmod = mrdfits(cfiles[i],1) & $
  clft = mrdfits(cfiles[i],2) & $
  nx2 = (size(clmod.modim))[1] & $
  ny2 = (size(clmod.modim))[2] & $
  xc2 = nx2/2 & $
  yc2 = ny2/2 & $
  cmodstack[*,*,i] = clmod.modim[xc2-hw:xc2+hw,yc2-hw:yc2+hw] & $
  testim = cml_lens_gauss([0.,clmod.fpar[1:*]], x=clft.x/float(magnif), $
   y=clft.y/float(magnif), subs=clft.subs, $
   psf=clft.psf, $
   rmin=clft.rmin, defl_x=clft.defl_x, defl_y=clft.defl_y) & $
  csrcstack[*,*,i] = testim[xc2-hw:xc2+hw,yc2-hw:yc2+hw] & $
endfor

; Loop with which to attempt caustic/critical-curve location:

; X and Y super-sampled arrays, in pixels:
pixsamp = 100L
subsamp = 5L
rangevec = findgen(2*pixsamp*subsamp)/float(subsamp) $
 - float(pixsamp) + 0.5 / float(subsamp)
bigx = rangevec # replicate(1., 2*pixsamp*subsamp)
bigy = transpose(bigx)
del = 1./float(subsamp)

; Tiny circle for finding the radial cut:
n_rad = 1000
r_rad = .001
t_rad = 2 * !pi * findgen(n_rad+1) / float(n_rad)
xrad = r_rad * cos(t_rad)
yrad = r_rad * sin(t_rad)
for i = 0L, nsp-1 do begin & $
;for i = 0L, 5 do begin & $
  print, i & $
  lmod = mrdfits(mfiles[i],1) & $
  lpar = lmod.fpar[0:4] & $
  sie_grad, bigx+del, bigy, xg_xp, yg_xp, par=lpar & $
  sie_grad, bigx-del, bigy, xg_xm, yg_xm, par=lpar & $
  sie_grad, bigx, bigy+del, xg_yp, yg_yp, par=lpar & $
  sie_grad, bigx, bigy-del, xg_ym, yg_ym, par=lpar & $
  xg_x = (xg_xp - xg_xm) / (2. * del) & $
  xg_y = (xg_yp - xg_ym) / (2. * del) & $
  yg_x = (yg_xp - yg_xm) / (2. * del) & $
  yg_y = (yg_yp - yg_ym) / (2. * del) & $
  imagnif = (1.- xg_x) * (1. - yg_y) - xg_y * yg_x & $
;  contour, imagnif, bigx, bigy, levels=[0.], path_info=path_info, path_xy=path_xy
  contour, imagnif, bigx, bigy, levels=[0.], path_info=path_info, path_xy=path_xy, /path_data_coords & $
  pmax = max(path_info.n, wmax) & $
;  atv, abs(imagnif)
  poff = path_info[wmax].offset & $
  pnpt = path_info[wmax].n & $
  cr_x = reform(path_xy[0,poff:poff+pnpt-1]) & $
  cr_y = reform(path_xy[1,poff:poff+pnpt-1]) & $
;  splot, cr_x, cr_y, ps=3, color=2
;  atvplot, cr_x, cr_y, ps=3
  ; Map the critical curves to the source plane:
  sie_grad, cr_x, cr_y, xg_cau, yg_cau, par=lpar & $
  x_cau = cr_x - xg_cau & $
  y_cau = cr_y - yg_cau & $
;  soplot, x_cau, y_cau, ps=3, color=3
;  atvplot, x_cau, y_cau, ps=3
; The radial "cut"
  sie_grad, xrad+lpar[3], yrad+lpar[4], xg_cut, yg_cut, par=lpar & $
  x_cut = xrad+lpar[3] - xg_cut & $
  y_cut = yrad+lpar[4] - yg_cut & $
  tag = 'tag' + string(i, format='(i4.4)') & $
  this_cut_x = create_struct(tag, x_cut) & $
  this_cut_y = create_struct(tag, y_cut) & $
  this_crit_x = create_struct(tag, cr_x) & $
  this_crit_y = create_struct(tag, cr_y) & $
  this_cau_x = create_struct(tag, x_cau) & $
  this_cau_y = create_struct(tag, y_cau) & $
  cut_x = (i eq 0) ? this_cut_x : struct_addtags(cut_x, this_cut_x) & $
  cut_y = (i eq 0) ? this_cut_y : struct_addtags(cut_y, this_cut_y) & $
  crit_x = (i eq 0) ? this_crit_x : struct_addtags(crit_x, this_crit_x) & $
  crit_y = (i eq 0) ? this_crit_y : struct_addtags(crit_y, this_crit_y) & $
  cau_x = (i eq 0) ? this_cau_x : struct_addtags(cau_x, this_cau_x) & $
  cau_y = (i eq 0) ? this_cau_y : struct_addtags(cau_y, this_cau_y) & $
endfor

i=-1

i=i+1
splot, crit_x.(i), crit_y.(i), xrange=[-50,50], yrange=[-50,50], color=3
soplot, cut_x.(i), cut_y.(i), color=4
soplot, cau_x.(i), cau_y.(i), color=2


;;;;;;;;;here

; The pretty images:
p_img = imgstack * (invstack gt 0.) + (bspstack + modstack) * (invstack le 0.)

;i=-1
;i=i+1 & atv, [p_img[*,*,i]-bspstack[*,*,i],modstack[*,*,i],srcstack[*,*,i],cmodstack[*,*,i],csrcstack[*,*,i]]


; Get various percentiles within the images:
lev_d = fltarr(nsp)
lev_r = fltarr(nsp)

npix = (2L*hw+1L)^2
for i = 0L, nsp-1 do begin & $
  vec_d = (bspstack[*,*,i])[*] & $
  vec_r = (modstack[*,*,i])[*] & $
  vec_d = vec_d[sort(vec_d)] & $
  vec_r = vec_r[sort(vec_r)] & $
; For B/W:
  lev_d[i] = vec_d[round(0.97*npix)] & $
  lev_r[i] = vec_r[round(0.99*npix)] & $
; For color:
;  lev_d[i] = vec_d[round(0.999*npix)] & $
;  lev_r[i] = vec_r[round(0.999*npix)] & $
endfor


;i=-1
;i=i+1 & atv, [p_img[*,*,i], bspstack[*,*,i]], max=lev97d[i], min=-0.25*lev97d[i]

;i=-1
;i=i+1 & atv, [p_img[*,*,i]- bspstack[*,*,i], modstack[*,*,i]], max=lev_r[i], min=-0.25*lev_r[i]


; Looks like, at least tentatively, we want 97th percentile for the
; direct scaling and 99th percentile for the residual scaling.


snames = 'SDSS' + strmid(sp.uniqname, 5)

; For internal development of make_prettypic_2.pro:
;i = 2
;image = p_img[*,*,i]
;resid = p_img[*,*,i] - bspstack[*,*,i]
;modim =  modstack[*,*,i]
;modres = p_img[*,*,i] - bspstack[*,*,i] - modstack[*,*,i]
;srcim = srcstack[*,*,i]
;modim2 = cmodstack[*,*,i]
;modres2 = p_img[*,*,i] - bspstack[*,*,i] - cmodstack[*,*,i]
;srcim2 = csrcstack[*,*,i]
;d_max=lev_d[i]
;d_min=-0.25*lev_d[i]
;r_max=lev_r[i]
;r_min=-0.25*lev_r[i]
;sysname='junk'
;x_cau = cau_x.(i)
;y_cau = cau_y.(i)
;x_crit = crit_x.(i)
;y_crit = crit_y.(i)
;x_cut = cut_x.(i)
;y_cut = cut_y.(i)

; For B/W:
;.com make_prettypic_2
for i = 0L, nsp-1 do $
;for i = 1L, 1L do $
 make_prettypic_2, p_img[*,*,i], p_img[*,*,i] - bspstack[*,*,i], $
 modstack[*,*,i],  p_img[*,*,i] - bspstack[*,*,i] - modstack[*,*,i], srcstack[*,*,i], $
 cmodstack[*,*,i],  p_img[*,*,i] - bspstack[*,*,i] - cmodstack[*,*,i], csrcstack[*,*,i], $
 d_max=lev_d[i], d_min=-0.25*lev_d[i], r_max=lev_r[i], r_min=-0.25*lev_r[i], $
 x_cau=cau_x.(i)*float(magnif),y_cau=cau_y.(i)*float(magnif),x_crit=crit_x.(i),$
 y_crit=crit_y.(i),x_cut=cut_x.(i)*float(magnif),y_cut=cut_y.(i)*float(magnif), $
 sysname=snames[i] ;, /showplot

; For color:
;for i = 0L, nsp-1 do $
;for i = 0L, 0L do $
; make_prettypic_2, p_img[*,*,i], p_img[*,*,i] - bspstack[*,*,i], $
; modstack[*,*,i],  p_img[*,*,i] - bspstack[*,*,i] - modstack[*,*,i], srcstack[*,*,i], $
; cmodstack[*,*,i],  p_img[*,*,i] - bspstack[*,*,i] - cmodstack[*,*,i], csrcstack[*,*,i], $
; d_max=lev_d[i], d_min=-0.05*lev_d[i], r_max=lev_r[i], r_min=-0.05*lev_r[i], $
; sysname=snames[i]

; Try some straightforward TeX formatting:
;for i = 0L, nsp/2 - 1 do begin & $
;  print, '\centerline{\scalebox{0.4}{\includegraphics{ppic' + snames[2*i] + $
;   '.eps}\hspace*{2pt}\includegraphics{ppic' + snames[2*i+1] + '.eps}}}' & $
;endfor


for i = 0L, nsp - 1 do print, $
 '\centerline{\scalebox{0.7}{\includegraphics{vppic' + snames[i] + '.eps}}}'



