pro plot_clust_cmd, gal, haloid

;

; Start by getting the size of the input array and creating a
; duplicate array to hold the output

ngal = size(gal,/n_elements)
tmp606 = fltarr(ngal)
tmp814 = fltarr(ngal)

; Loop through the input catalog, selecting the objects that match the
; requested cluster ID

print, ''
print, 'Number of objects in input catalog: ',ngal
print, 'Searching for cluster ID: ',haloid
count = 0
for i = 0L,ngal-1L do begin
   if (gal[i].HaloID eq haloid) then begin
      print, i, gal[i].HaloID, gal[i].Mvir, gal[i].Type,gal[i].Mag[6]
      tmp606[count] = gal[i].Mag[5]
      tmp814[count] = gal[i].Mag[6]
      count = count+1
   end
endfor

; Trim the output catalog
print, 'Number of objects associated with this cluster:  ',count+1
f606 = tmp606[0:count-1]
f814 = tmp814[0:count-1]
hstcolor = f606 - f814

; Temporary printouts to check error
print, ''
print, ''
for i = 0,count-1 do begin
   print, i, f814[i], hstcolor[i]
endfor

!p.multi = [0, 1, 2]
plot, f814, hstcolor, /ynozero, psym=1
!p.multi = 0
end
