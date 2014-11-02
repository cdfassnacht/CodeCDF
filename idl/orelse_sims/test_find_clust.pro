pro test_find_clust, gal, mclust, outcat

;
; test_find_clust.pro
;
; Description: 
;  Takes an array of galaxy structures, read in by Risa's 
;  read_galaxy_subtree_rich9.pro and returns all of the objects that
;  are associated with a central halo that has a mass greater than the
;  passed mclust value.
; 
; Passed parameters
;  gal      - Array of galaxy structures produced by Risa's
;             read_galaxy_subtree_rich9
;  mclust   - Lower limit for defining a central halo as a cluster
;             (e.g., 1.0e14)
;  outcat   - Output catalog variable name.  The catalog gets created
;             by this procedure and gets returned at the end

; Start by getting the size of the input array and creating a
; duplicate array to hold the output

ngal = size(gal,/n_elements)
tmpcat = replicate(gal[0],ngal)

; Loop through the input catalog, selecting the objects with large
; central Mvirs

print, ''
print, 'Number of objects in input catalog: ',ngal
print, 'Minimum CentralMvir for cluster membership: ',mclust
count = 0
for i = 0L,ngal-1L do begin
   if (gal[i].CentralMvir gt mclust) then begin
      if gal[i].Type eq 0 then begin
         print, i, gal[i].HaloID, gal[i].CentralMvir, gal[i].Mvir, gal[i].Type
      end
      tmpcat[count] = gal[i]
      count = count+1
   end
endfor

; Trim the output catalog
print, 'Number of objects associated with clusters:  ',count
outcat = tmpcat[0:count-1]

end
