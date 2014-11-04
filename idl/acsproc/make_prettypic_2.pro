pro make_prettypic_2, image, resid, modim, modres, srcim, $
 modim2, modres2, srcim2, d_max=d_max, $
 r_max=r_max, d_min=d_min, r_min=r_min, sysname=sysname, showplot=showplot, $
 x_cau=in_x_cau, y_cau=in_y_cau, x_crit=in_x_crit, y_crit=in_y_crit, x_cut=in_x_cut, y_cut=in_y_cut

x_cau = [in_x_cau, in_x_cau[0]]
y_cau = [in_y_cau, in_y_cau[0]]
x_crit = [in_x_crit, in_x_crit[0]]
y_crit = [in_y_crit, in_y_crit[0]]
x_cut = [in_x_cut, in_x_cut[0]]
y_cut = [in_y_cut, in_y_cut[0]]


imx = (size(image))[1]
imy = (size(image))[2]

big_im = rebin(image, 2*imx, 2*imy, /sample)
big_res = rebin(resid, 2*imx, 2*imy, /sample)

bvalu = 255B
;bvalu = 0B
sc1 = padimage(bytscl(modim, min=r_min, max=r_max), 1, value = bvalu)
sc2 = padimage(bytscl(modres, min=r_min, max=r_max), 1, value = bvalu)
sc3 = padimage(bytscl(srcim, min=r_min, max=r_max), 1, value = bvalu)
sc1b = padimage(bytscl(modim2, min=r_min, max=r_max), 1, value = bvalu)
sc2b = padimage(bytscl(modres2, min=r_min, max=r_max), 1, value = bvalu)
sc3b = padimage(bytscl(srcim2, min=r_min, max=r_max), 1, value = bvalu)
sc4 = padimage([[sc1b, sc2b, sc3b],[sc1, sc2, sc3]], 1, value = bvalu)

sc5 = padimage(bytscl(big_im, min=d_min, max=d_max), 1, value=bvalu)
sc6 = padimage(bytscl(big_res, min=r_min, max=r_max), 1, value=bvalu)
sc7 = padimage([sc5, sc6], 2, value=bvalu)

xoff_1 = (size(sc7))[1] + (size(modim))[1] / 2
xoff_2 = (size(sc7))[1] + (size(modim))[1] + 2 + (size(modres))[1] + 2 + (size(srcim))[1] / 2
yoff_1 = (size(modim2))[2]/2 + 2
yoff_2 = (size(modim2))[2] + 2 + 2 + (size(modim))[2] / 2


sc8 = [sc7, sc4[2:*,*]]
;sc8 = [sc7, sc4]

nxd = (size(sc8))[1]
nyd = (size(sc8))[2]
aspectrat = float(nyd) / float(nxd)

xsiz = 8.
ysiz = xsiz * aspectrat
ofile = 'vppic' + sysname + '.eps'
print, ofile
spawn, 'rm -f ' + ofile
set_plot, 'ps'
device, filename=ofile, /encapsulated, /inches, $
 xsize=xsiz, ysize=ysiz, bits_per_pixel=8 ;, /color
;loadct, 34
tv, 255B-sc8
;tv, sc8
;loadct, 0
plot, [0], [0], /nodata, position=[0,0,1,1], $
 xrange=[-0.5,nxd-0.5], yrange=[-0.5,nyd-0.5], xstyle=5, ystyle=5, /noerase
oplot, [10.,10.+40.], [nyd-10., nyd-10.], thick=4
oplot, [10.,10.], [nyd-14., nyd-6.], thick=4
oplot, [10.+40.,10.+40.], [nyd-14., nyd-6.], thick=4
;oplot, x_crit+xoff_1, y_crit+yoff_2, ps=3, color=255
;oplot, x_cau+xoff_2, y_cau+yoff_2, ps=3, color=255
oplot, x_crit+xoff_1, y_crit+yoff_2, color=255, thick=3
oplot, x_cau+xoff_2, y_cau+yoff_2, color=255, thick=3
wh = where((abs(x_cut) le (imx/2)) and (abs(y_cut) le (imy/2)), nwh)
;stop
if (nwh gt 1) then begin
  idiff = wh[1:*] - wh[0:nwh-2]
  chunk_lo = where(idiff gt 1, nbreak)
  chunk_lo = (nbreak gt 0) ? [0L, chunk_lo+1] : 0L
  chunk_hi = (nbreak gt 0) ? [chunk_lo[1:*]-1, nwh-1] : nwh-1
  for i = 0L, nbreak do oplot, $
   (x_cut+xoff_2)[wh[chunk_lo[i]:chunk_hi[i]]], $
   (y_cut+yoff_2)[wh[chunk_lo[i]:chunk_hi[i]]], thick=4, color=255
endif
xyouts, 0.025, 0.05, sysname, charsize=1.5, charthick=4, /normal ;, color=0
xyouts, 0.032, 0.86, '1"', charsize=1.5, charthick=4, /normal
xyouts, 4./7.+.01, 0.92, 'SIE', charsize=1.2, charthick=4, /normal
xyouts, 4./7.+.01, 0.42, 'LTM', charsize=1.2, charthick=4, /normal
xyouts, 6./7.+.01, 0.92, 'x2', charsize=1.2, charthick=4, /normal
xyouts, 6./7.+.01, 0.42, 'x2', charsize=1.2, charthick=4, /normal
device, /close
set_plot, 'x'
if keyword_set(showplot) then spawn, 'display ' + ofile + ' &'

return
end
