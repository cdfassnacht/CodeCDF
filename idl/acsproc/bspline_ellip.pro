;+
;
; NAME: bspline_ellip
;
; PURPOSE: Elliptical bspline image modeler.  Written to
;   work with MPFIT if keyword "deviates" is set.
;
; ARGUMENTS:
;  par: parameter vector:
;    p[0] = x-center
;    p[1] = y-center
;    p[2] = axis ratio
;    p[3] = major-axis PA (degrees countercl. from x axis)
;    p[4:7], p[8:11], etc.: parameters for more components.
;
;  x, y: x and y coordinate images.  If RA and Dec, note that
;    you should use x=dec and y=ra for a rhight-handed coordinate
;    system on the sky.
;  data: data image to fit.
;  invvar: inverse-variance of data.
;  psf: PSF image, if PSF convolution is to be done.
;  deviates: set this keyword to return vector of weighted
;    deviates as desired by MPFIT.
;  ntheta: vector of multipole identifiers, if multipole fitting
;    is to be done.
;  rbkpt: vector of radial breakpoint locations.
;
; WRITTEN: abolton@cfa 2006may.
;  Fixed "action" bug for multi-component case, abolton@cfa 2006jun.
;
;-

function bspline_ellip, par, x=x, y=y, data=data, invvar=invvar, $
 psf=psf, deviates=deviates, ntheta=ntheta, rbkpt=rbkpt, sset=sset

deg2rad = !pi / 180.
psf_loc = keyword_set(psf) ? psf : 0B

if (not keyword_set(ntheta)) then ntheta = 0
nmultip = n_elements(ntheta)

l_data = data[*]
l_invvar = invvar[*]

; Generate action matrices:
nimg = n_elements(par) / 4

nx = (size(data))[1]
ny = (size(data))[2]

for i = 0, nimg-1 do begin
  sset = 0
  xc = par[4*i]
  yc = par[4*i+1]
  q  = par[4*i+2]
  pa = par[4*i+3] * deg2rad

  xp = (x-xc) * cos(pa) + (y-yc) * sin(pa)
  yp = (y-yc) * cos(pa) - (x-xc) * sin(pa)
  r = sqrt((xp^2)*q + (yp^2)/q)
  theta = atan(yp, xp)

  sub_action = acsproc_bspline_psf_action(r, theta, psf=psf_loc, $
   rbkpt=rbkpt, ntheta=ntheta, sset=sset, /silent)
  nco = (size(sub_action))[3]
  sub_action = reform(temporary(sub_action), nx*ny, nco)

  action = (i eq 0) ? sub_action : [[action], [sub_action]]
endfor

iaction = action * (l_invvar # replicate(1., nco*nimg))
alpha = transpose(iaction) # action
beta = transpose(iaction) # l_data
ialpha = invert(alpha)
fcoeff = ialpha # beta
fmodel = action # fcoeff
fmodel = reform(fmodel, nx, ny)

sset.xmin = -!pi
sset.xmax = !pi
sset = replicate(sset, nimg)

icoeff = 1./ialpha[lindgen(nco*nimg),lindgen(nco*nimg)]

for i = 0, nimg-1 do begin
  sset[i].coeff = reform(fcoeff[i*nco:(i+1)*nco-1], nmultip, nco/nmultip)
  sset[i].icoeff = reform(icoeff[i*nco:(i+1)*nco-1], nmultip, nco/nmultip)
endfor

if keyword_set(deviates) then fmodel = ((data - fmodel) * sqrt(invvar))[*]

return, fmodel
end
