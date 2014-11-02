function stamp_radec , fits, ra, dec, size, pa=pa, outhead=outhead

;
; Function stamp_radec
;
; Description: Takes an input image and cuts out a smaller "postage
;  stamp" image of size nx by ny centered at (ra, dec).  Returns
;  the postage stamp image, with the proper rotation so that north is
;  up and east is to the left.  A complete copy of Lexi's eidol_gim.
;
; Inputs: fits       (string)      name of input fits file
;         ra         (double)      RA to be
;                                   used as the center of the postage stamp
;         dec        (double)      Dec position (in input image) to be
;                                   used as the center of the postage stamp
;         size       (int)         output image size
;         [pa=]      (float)       optional additional rotation to perform.
;         [outhead=] (stringarr)   optional string array to hold header
;                                   card info
;
; Output: [none for right now]
;
; Revision history:
;  2003May04 Chris Fassnacht -- First hack-up of Lexi's eidol_gim
;  2003May06 Chris Fassnacht -- Added outhead optional parameter and 
;                                associated header modification to make
;                                more like poststamp.pro
;

    if n_elements(pa) eq 0 then pa=0.
    size=2l*size
    hd=headfits(fits)
    sz=fxpar(hd,'NAXIS*')
    getrot,hd,rotn
    adxy, hd, ra, dec, x0, y0
    img=make_array(size,size,/float)
    x1 = long(x0 - size/2.) > 0l
    x2 = long(x0 + size/2.) < sz[0]-1l
    y1 = long(y0 - size/2.) > 0l
    y2 = long(y0 + size/2.) < sz[1]-1l
    print,x1,x2,y1,y2
    if (x1 gt sz[0]-1l or x2 lt 0 or y1 gt sz[1]-1l or y2 lt 0) then begin
        print,'stamp_radec: error: pointing not within mosaic boundaries'
        return,-1
    endif
    fxread,fits,imgtmp,hd,x1,x2,y1,y2
    cx1 = long(size-x2-1l) > 0l
    cx2 = long(size-1l) < long(sz[0]-x1-1l) < long(size-1l)
    cy1 = long(size-y2-1l) > 0l
    cy2 = long(size-1l) < long(sz[1]-y1-1l) < long(size-1l)
    img[cx1:cx2,cy1:cy2] = imgtmp;[x1:x2,y1:y2]
    img=rot(img,rotn+pa,missing=0.)
    size=size/2l
    nx1=size/2
    ny1=size/2
    nx2=size/2*3
    ny2=size/2*3
    img=img[nx1:nx2,ny1:ny2]
    outhead = hd
    fxaddpar, outhead, 'naxis1', size
    fxaddpar, outhead, 'naxis2', size
    histstr = '  Postage stamp cutout centered at (' + strtrim(x0,1)
    histstr = histstr + ',' + strtrim(y0,1) + ')'
    fxaddpar, outhead, 'history', histstr
    return,img

end
