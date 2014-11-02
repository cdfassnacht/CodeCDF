pro get_clust_env, gal, central_halo, rmax, outfits
  
; Given a array of halo structures, returns all halos that are
; within rmax comoving h^-1 Mpc of the central halo.  The halos are
; returned as an array of structures

; Start by getting the size of the input array and creating a
; duplicate array to hold the output

ngal = size(gal,/n_elements)
tmpout = replicate(gal[0],ngal)

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
dist1 = fltarr(ngal)
for i = 0L,ngal-1L do begin
   dx = gal[i].Pos[0] - xcent
   dy = gal[i].Pos[1] - ycent
   dz = gal[i].Pos[2] - zcent
   dpos = sqrt(dx*dx + dy*dy + dz*dz)
   if dpos lt rmax then begin
      tmpout[count] = gal[i]
      dist1[count] = dpos
      count = count+1
      if (gal[i].Type eq 1) or (gal[i].Type eq 0) then count1 = count1 + 1
   end
endfor

; Trim and sort the output catalog
print, ''
print, 'Volume defined as a sphere of radius',rmax
print, 'Number of objects in this volume:  ',count
print, 'Number of Type 0 or Type 1 objects in this volume',count1
tmp2 = tmpout[0:count-1]
dist2 = dist1[0:count-1]
delvar, tmpout
delvar, dist1
sortind = sort(dist2)
outgals=tmp2[sortind]
delvar, tmp2
delvar, dist2
print, 'After sorting index 0 corresponds to ID:',outgals[0].HaloID
print, ''

; Write the output list to a binary table fits file
mwrfits, outgals, outfits

end
