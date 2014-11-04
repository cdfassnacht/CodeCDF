function bspline_ellip_valu, x, y, par, sset, psf=psf, r=r, theta=theta

nimg = n_elements(sset) < (n_elements(par)/4)

deg2rad = !pi / 180.

for i = 0, nimg-1 do begin
  xc = par[4*i]
  yc = par[4*i+1]
  q  = par[4*i+2]
  pa = par[4*i+3] * deg2rad

  xp = (x-xc) * cos(pa) + (y-yc) * sin(pa)
  yp = (y-yc) * cos(pa) - (x-xc) * sin(pa)
  r = sqrt((xp^2)*q + (yp^2)/q)
  theta = atan(yp, xp)

  sub_image = bspline_radial_valu(r, theta, sset[i])
  image = (i eq 0) ? sub_image : image + sub_image
endfor

if keyword_set(psf) then image = convolve(image, psf)

return, image
end
