pro test_plot_xy, gal

; This procedure is mostly just to check what the different HaloIDs
; and other paramters mean.  It will go through an array of structures
; and plot all objects that we believe are associated with one central
; halo.  This assignment is done by finding an object in the array
; with a type==0 and then assigning all halos after that object in the
; array to that central halo until the next occurence of a Type 0 object.

; Start by getting the size of the input array and creating a
; duplicate array to hold the output

ngal = size(gal,/n_elements)
tmpx = fltarr(ngal)
tmpy = fltarr(ngal)
tmpx1 = fltarr(ngal)
tmpy1 = fltarr(ngal)

; Loop through the input catalog, selecting the objects that follow
; the first occurance of a Type 0 object
; requested cluster ID

print, ''
print, 'Number of objects in input catalog: ',ngal
print, ''
i = 0
count  = 0
count1 = 0
samehalo = 1
foundfirst = 0
; Loop over gal, assuming that the first entry has Type 0
while samehalo eq 1 do begin
   if (gal[i].Type eq 0) and (count gt 0) then begin
      print, 'Found next object with Type 0 at i=',i
      samehalo = 0
   endif else begin
      tmpx[count] = gal[i].Pos[0]
      tmpy[count] = gal[i].Pos[1]
      count = count+1
      if (gal[i].Type eq 1) or (gal[i].Type eq 0) then begin
         tmpx1[count1] = gal[i].Pos[0]
         tmpy1[count1] = gal[i].Pos[1]
         count1 = count1 + 1
      endif
      i = i+1
   endelse
endwhile

; Trim the output catalogs
print, 'Number of objects associated with this central halo:  ',count
x = tmpx[0:count-1]
y = tmpy[0:count-1]
print, 'Number of Type 1 objects associated with the central halo',count1
x1 = tmpx1[0:count1-1]
y1 = tmpy1[0:count1-1]

plot, x, y, /ynozero, psym=1
;loadct, 4
oplot, x1, y1, psym=4, symsize=4., thick=3, color=255
end
