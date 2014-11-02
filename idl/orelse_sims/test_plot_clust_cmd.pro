pro test_plot_clust_cmd, gal, central_halo, rmax, rprojmax
  
; Given a array of halo structures, plots all objects that fall within
; a given distance of the position of a central halo

; Start by getting the size of the input array and creating a
; duplicate array to hold the output

ngal = size(gal,/n_elements)
tmpV = fltarr(ngal)
tmpI = fltarr(ngal)
tmpx = fltarr(ngal)
tmpy = fltarr(ngal)
tmpz = fltarr(ngal)
tmpdx = fltarr(ngal)
tmpdy = fltarr(ngal)
tmpdz = fltarr(ngal)
tmpx1 = fltarr(ngal)
tmpy1 = fltarr(ngal)
tmpz1 = fltarr(ngal)

print, ''
print, 'Number of objects in input catalog: ',ngal
xcent = central_halo.Pos[0]
ycent = central_halo.Pos[1]
zcent = central_halo.Pos[2]
print, 'Central Halo Parameters'
print, '-----------------------'
print, '  ID:    ',central_halo.HaloID
print, '  x:     ',xcent
print, '  y:     ',ycent
print, '  z:     ',zcent

; Loop through the input catalog, selecting the objects that are
; within the requested distance, in comoving h^-1 Mpc

count = 0L
count1 = 0L
for i = 0L,ngal-1L do begin
   dx = gal[i].Pos[0] - xcent
   dy = gal[i].Pos[1] - ycent
   dz = gal[i].Pos[2] - zcent
   dpos = sqrt(dx*dx + dy*dy + dz*dz)
   if dpos lt rmax then begin
      tmpV[count] = gal[i].Mag[5]
      tmpI[count] = gal[i].Mag[6]
      tmpx[count] = gal[i].Pos[0]
      tmpy[count] = gal[i].Pos[1]
      tmpz[count] = gal[i].Pos[2]
      tmpdx[count] = dx
      tmpdy[count] = dy
      tmpdz[count] = dz
      if (gal[i].Type eq 1) or (gal[i].Type eq 0) then begin
         tmpx1[count1] = gal[i].Pos[0]
         tmpy1[count1] = gal[i].Pos[1]
         tmpz1[count1] = gal[i].Pos[2]
         count1 = count1 + 1
      endif
      count = count+1
   end
endfor

; Trim the output catalogs
print, ''
print, 'Volume defined as a sphere of radius',rmax
print, 'Number of objects in this volume:  ',count
x = tmpx[0:count-1]
y = tmpy[0:count-1]
z = tmpz[0:count-1]
hstV = tmpV[0:count-1]
hstI = tmpI[0:count-1]
dxf = tmpdx[0:count-1]
dyf = tmpdy[0:count-1]
dzf = tmpdz[0:count-1]
print, 'Number of Type 0 or Type 1 objects in this volume',count1
x1 = tmpx1[0:count1-1]
y1 = tmpy1[0:count1-1]
y1 = tmpy1[0:count1-1]
delvar, tmpx,tmpy,tmpz,tmpx1,tmpy1,tmpz1,tmpdx,tmpdy,tmpdz

; New loop for plotting
tmppx = fltarr(count)
tmppy = fltarr(count)
tmppV = fltarr(count)
tmppI = fltarr(count)
pcount = 0L
for i=0L,count-1L do begin
   pdp = sqrt(dxf[i]*dxf[i] + dyf[i]*dyf[i])
   if pdp lt rprojmax then begin
      tmppx[pcount] = x[i]
      tmppy[pcount] = y[i]
      tmppV[pcount] = hstV[i]
      tmppI[pcount] = hstI[i]
      pcount += 1
   end
endfor
print, ''
print, 'Cylinder defined as:'
print, '  radius: ',rprojmax
print, '  line-of-sight length: ',rmax
print, 'Found',pcount-1,' objects in projected cylinder'
px = tmppx[0:pcount-1]
py = tmppy[0:pcount-1]
pV = tmppV[0:pcount-1]
pI = tmppI[0:pcount-1]
pVI = pV - pI

!p.multi = [0, 1, 2]
plot, pI, pVI, /ynozero, psym=1
!p.multi = 0
;loadct, 4
;oplot, x1, y1, psym=4, symsize=4., thick=3, color=255

end
