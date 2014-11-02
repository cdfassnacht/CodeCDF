; clean the four WIRC quadrants for residual background tilt
pro wirc_cleanquadrants,in,out, FIXEDGES=fixedges

; yank the quadrants
;
;   3 4
;   1 2
;
quad1=in[0:1023,0:1023]
quad2=in[1024:2047,0:1023]
quad3=in[0:1023,1024:2047]
quad4=in[1024:2047,1024:2047]



;----- fit the slope model to each quadrant------

model=fltarr(1024)
model[*]=0


; clean quadrant1
for i=0,1023 do begin
    model[i]=median(quad1[i,*])
endfor

for i=0,1023 do begin
    quad1[*,i]=quad1[*,i]-model
endfor


; clean quadrant2
for i=0,1023 do begin
    model[i]=median(quad2[*,i])
endfor

for i=0,1023 do begin
    quad2[i,*]=quad2[i,*]-model
endfor

; clean quadrant3
for i=0,1023 do begin
    model[i]=median(quad3[*,i])
endfor

for i=0,1023 do begin
    quad3[i,*]=quad3[i,*]-model
endfor


; clean quadrant4
for i=0,1023 do begin
    model[i]=median(quad4[i,*])
endfor

for i=0,1023 do begin
    quad4[*,i]=quad4[*,i]-model
endfor

;-----------------------------------------------


; fix the edges to all have the same median value
if keyword_set(FIXEDGES) then begin

print,'Fixing edges.'

; horizontal edges
edge1=(median(quad1[20:1023,1000:1023])+median(quad1[1000:1023,20:1023]))/2
edge2=(median(quad2[0:1020,1000:1023])+median(quad2[0:20,20:1023]))/2
edge3=(median(quad3[20:1023,0:20])+median(quad3[1000:1023,0:1020]))/2
edge4=(median(quad4[0:1020,0:20])+median(quad4[0:20,0:1020]))/2
print,edge1,edge2,edge3,edge4


avgedge=(edge1+edge2+edge3+edge4)/4
quad1=quad1+(avgedge-edge1)
quad2=quad2+(avgedge-edge2)
quad3=quad3+(avgedge-edge3)
quad4=quad4+(avgedge-edge4)

endif





; dump the quadrants into the output array
out[0:1023,0:1023]=quad1
out[1024:2047,0:1023]=quad2
out[0:1023,1024:2047]=quad3
out[1024:2047,1024:2047]=quad4

end