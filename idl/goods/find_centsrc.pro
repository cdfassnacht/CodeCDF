pro find_centsrc, image, secat, best
;
; Procedure find_centsrc
;
; Description: Given an image and a secat array, finds the source closest
;  to the center of the image
;
; Inputs: image      (fltarr)      image
;         secat      (structarr)   catalog info
;         best       (structure)   secat structure containing position
;                                   of best match
;
; Revision history:
;  2003Jul02 Chris Fassnacht -- Moved from goods_er
;

; Check input format

if n_params() lt 3 then begin
    print, ''
    print, 'syntax: find_centsrc, image, secat, best'
    print, ''
    return
endif

; Get size of input image

nx = (size(image,/dimen))[0]
ny = (size(image,/dimen))[1]
print, ''
print, 'find_centsrc: Input image dimensions: ', nx, ny
centx = nx / 2.0
centy = ny / 2.0

; Go through input catalog and find object closest to the image center.

ncat = n_elements(secat)
besti = 0
bestoffset = nx > ny
for i=0,ncat-1 do begin
   dx = centx - secat[i].x
   dy = centy - secat[i].y
   offset = sqrt(dx*dx + dy*dy)
   if offset lt bestoffset then begin
      bestoffset = offset
      besti = i
   endif
endfor

; Set up structure for holding best object info

def_secat, best
best = secat[besti]

print, ''
print, 'find_centsrc: Object closest to center is ',besti
print, 'find_centsrc: 	with x,y = ',best.x,best.y
print, 'find_centsrc: 	and an offset of ',bestoffset

end
