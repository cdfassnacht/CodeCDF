pro circ_phot, infile, xphot, yphot, rphot, rbkgd, photstat, $
                zp=zp, dodisplay=dodisplay
;
; Procedure circ_phot
;
; Description: Given a fits file, a central position, aperture radius/radii,
;               and background region radii, compute the background-subtracted
;               counts inside the aperture(s).
;
; Inputs: infile     (string)      base name of input FITS file (root+id number)
;         xphot      (float)       x position for center of aperture(s)
;         yphot      (float)       y position for center of aperture(s)
;         rphot      (floatarr)    array of aperture radii (can contain only
;                                   one value)
;         rbkgd      (floatarr)    inner and outer radii for background region
;         photstat   (structure)   structure containing photometric info
;                                   (created by this procedure)
;         [zp=]      (float)       photometric zero point (Default is 25.0)
;         [/dodisplay]             display images if set
;
;
; Revision history:
;  2003Apr29 Chris Fassnacht -- First working version.
;  2003Jun14 Chris Fassnacht -- Added more parameters to output structure.
;

; Check input format

if n_params() lt 6 then begin
    print, ''
    print, 'syntax: circ_phot, infile, xphot,yphot,rphot,rbkgd,photstat,'
    print, '          [,/dodisplay]'
    print, ''
    return
endif

; Set up defaults

if n_elements(zp) eq 0 then begin
   zp = 25.0
   print,'circ_phot: Setting zeropoint to 25.0'
endif else begin
   print,'circ_phot: Using passed zeropoint of',zp
endelse

; Set up structure for holding best object info

best =  create_struct('x',0.0,$
                      'y',0.0,$
                      'a',0.0,$
                      'b',0.0,$
                      'pa',0.0,$
                      'other',0.0)


; Read in files

print, ''
print, 'circ_phot: Reading input files:'
print, 'circ_phot:   ',infile
im = mrdfits(infile, 0, inhead)

; Get size of input image

sz = (size(im,/dimen))
nx = (size(im,/dimen))[0]
ny = (size(im,/dimen))[1]
print, ''
print, 'circ_phot: Input image dimensions: ', nx, ny

; Display input image

if keyword_set(dodisplay) then dispim, im

; Get distances from center of aperture

dist_circle, r, sz, xphot, yphot, /double

; Get background mean counts

if n_elements(rbkgd) eq 1 then begin
   print, ''
   print, $
    '**** Only one value for rbkgd -- taking background as all r > rbkgd ****'
   bkgd = where(r gt rbkgd[0])
   if keyword_set(dodisplay) then $
      tvcircle, rbkgd[0], xphot, yphot, color=65535, /data
endif else begin
   print, ''
   print, 'Taking background in annulus defined by first 2 elements of rbkgd'
   bkgd = where(r gt rbkgd[0] and r lt rbkgd[1])
   if keyword_set(dodisplay) then begin
      tvcircle, rbkgd[0], xphot, yphot, color=65535, /data
      tvcircle, rbkgd[1], xphot, yphot, color=65535, /data
   endif
endelse
nbkgd = n_elements(bkgd)
tmp = im[bkgd]
bkgdstat = moment(tmp,/double)
bkgdmean = bkgdstat[0]
print, 'Mean background is ',bkgdmean

; Loop through apertures, computing count-rate 

print, ''
print, 'Calculating net counts, m_inst, and m_phot'
print, '--------------------------------------------
for i=0,(n_elements(rphot))-1 do begin
   if keyword_set(dodisplay) then $
      tvcircle, rphot[i], xphot, yphot, color=255, /data
   good = where(r lt rphot[i])
   ngood = n_elements(good)
   rawcts = total(im[good])
   netcts = rawcts - (1.0d * bkgdmean * ngood)
   instmag = -2.5 * alog10(netcts)
   photmag = instmag + zp
   if i eq 0 then begin
      raw = rawcts
      net = netcts
      aparea = ngood
      m_inst = instmag
      m_phot = photmag
      zp_phot = zp
   endif else begin
      raw = [raw, rawcts]
      net = [net, netcts]
      aparea = [aparea, ngood]
      m_inst = [m_inst, instmag]
      m_phot = [m_phot, photmag]
      zp_phot = [zp_phot, zp]
   endelse
   print, netcts, instmag, photmag
endfor

; Put all the info into a structure

photstat = create_struct( $
            'raw', raw, $
            'net', net, $
            'bkgdmean', bkgdmean, $
            'aparea', aparea, $
            'm_inst', m_inst, $
            'm_phot', m_phot, $
            'zp', zp_phot $
           )

end
