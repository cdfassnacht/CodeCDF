function f, b
  common tmp2, p2

  r_e=p2[1]
  nser=p2[6]
  r_b=p2[7]
  alpha=p2[8]
  gama=p2[9]
  return, 1.0+igamma(2.0*nser, b*(r_b/r_e)^(1.0/nser))-2.0*igamma(2.0*nser, b)
end

function df, b
  common tmp2, p2
  nser=p2[6]
  return, exp(-b)*b^(nser+1)
end


function core_sersic, x, y, p, basis=basis, rmin=rmin, psf=psf, logpars=logpars, rpow=rpow, subs=subs, _Extra=junk

common tmp, zmod
common tmp2, p2
!Except=0
p2=p
deg2rad=!pi/180.

ss_loc=keyword_set(subs) ? fix(subs) : 1

if (n_elements(x) ne n_elements(y)) then begin
  splog, ' x and y dimensions mis-matched!'
  return, 0
end

npar=10L
ndev=n_elements(p)/npar
if (ndev lt 1) then begin
  splog, ' insufficient number of parameters!'
  return, 0
endif


rml=keyword_set(rmin) ? rmin : 0.

nx=(size(x))[1]/ss_loc
ny=(size(x))[2]/ss_loc

; Loop over all specified components:
zmod=fltarr(nx, ny, ndev)
for i=0L, ndev-1 do begin
  ioff=i * npar
  a=p[0+ioff]
  sig=p[1+ioff]
  xcen=p[2+ioff]
  ycen=p[3+ioff]
  q=p[4+ioff]
  phi=p[5+ioff] * deg2rad
  nser=p[6+ioff]
  r_b=p[7+ioff]
  alpha=p[8+ioff]
  gama=p[9+ioff]

;  b_ser_0=(1.9992*nser-0.3271)
;  tag=0
;  while (tag eq 0) do begin &$
;     b_ser_1=b_ser_0-f(b_ser_0)/df(b_ser_0) &$
;     if (abs(b_ser_1-b_ser_0)/b_ser_0 gt 10.^(-1)) then tag=0 else tag=1 &$
;     b_ser_0=b_ser_1 &$
;  endwhile
;  b_ser=b_ser_1

  b_ser=(1.9992*nser-0.3271)
  if (b_ser le 0.0) then b_ser=1.0D

  xp=(x-xcen)*cos(phi)+(y-ycen)*sin(phi)
  yp=(y-ycen)*cos(phi)-(x-xcen)*sin(phi)
  r_ell=sqrt((xp^2)*q+(yp^2)/q+rml^2)

  if (alpha gt 50.0) then begin 
     z_this=0.0*r_ell
     wh_PL=where(x le r_b, c_PL)
     wh_Sersic=where(x gt r_b, c_Sersic)
     if (c_PL gt 0) then z_this[wh_PL]=a*exp(-b_ser*(r_b/sig)^(1.0D/nser))*(r_b/x[wh_PL])^gama
     if (c_Sersic gt 0) then z_this[wh_Sersic]=a*exp(-b_ser*(x[wh_Sersic]/sig)^(1.0D/nser))
  endif else begin 
     z_this=a*(1.0+(r_b/r_ell)^alpha)^(gama/alpha)*exp(-b_ser*((r_ell^alpha+r_b^alpha)/sig^alpha)^(1.0/(alpha*nser)))
;  endelse
  if (ss_loc gt 1) then z_this=xybinsub(z_this, ss_loc, ss_loc)/float(ss_loc^2)
  if keyword_set(psf) then z_this=convolve(z_this, psf)
  if keyword_set(basis) then zmod[*,*,i]=z_this $
   else zmod=zmod+z_this
endfor

return, zmod

end
