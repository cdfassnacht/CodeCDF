pro eval_resid, image, resid, xcent, ycent, rap, rbkgd, outstruct, $
                indir=indir, rsddir=rsddir, dodisplay=dodisplay

;
; Procedure eval_resid
;
; Description: Given two fits files containing an image and a residual created 
;               by subtracting a fit from that image, compares the residual 
;               level to the input level within an aperture.  Information is
;               returned in a structure.
;
; Inputs: image      (string)      fits filename for input image
;         resid      (string)      fits filename for residual image
;         xcent      (float)       x position for center of aperture(s)
;         ycent      (float)       y position for center of aperture(s)
;         rap        (floatarr)    array of aperture radii (can contain only
;                                   one value)
;         rbkgd      (floatarr)    inner and outer radii for background region
;         outstruct  (structure)   output structure containing image info
;         outbase    (string)      base name for output files
;         [indir=]   (string)      name of directory containing residuals
;                                   default is 'z-band'
;         [rsddir=]  (string)      name of directory containing residuals
;                                   default is 'Resid'
;         [shape=]   (floatarray)  a three element vector describing
;                                   the shape of the image, as found
;                                   by SExtractor.  The elements are:
;                                   shape[0] = "size" (fwhm for now)
;                                   shape[1] = ellipticity (1 - b/a)
;                                   shape[2] = PA (PA of major axis,
;                                    astronomical convention)
;         [/dodisplay]             display images if set
;
;
; Revision history:
;  2003Jun14 Chris Fassnacht -- First version.
;

; Check input format

if n_params() lt 7 then begin
    print, ''
    print, 'syntax: eval_resid, image,resid,xcent,ycent,rap,rbkgd,outstruct'
    print, '          [indir=indir,rsddir=rsddir,/dodisplay]'
    print, ''
    print, 'Default values:'
    print, "  indir:  'z-band'"
    print, "  rsddir: 'Resid'"
    print, ''
    return
endif

; Set directory names and defaults

if not keyword_set(dodisplay) then dodisplay=0b
if not keyword_set(indir) then indir = 'z-band'
if not keyword_set(rsddir) then rsddir = 'Resid'

; Set filenames

infits = indir+'/'+image
rsdfits = rsddir+'/'+resid

; Evaluate the images by calling circ_phot

circ_phot, infits,xcent,ycent,rap,rbkgd,instat,dodisplay=dodisplay
circ_phot, rsdfits,xcent,ycent,rap,rbkgd,rsdstat,dodisplay=dodisplay

; Extract the relevant info from the photometric structures

intot = instat.net
rsdtot = rsdstat.raw
rsdrat = rsdtot / intot

; Create the output structure

outstruct = create_struct( $
             'intot', instat.raw, $
             'innet', instat.net, $
             'rsdtot', rsdstat.raw, $
             'rsdrat', rsdtot / instat.net, $
             'r_ap', rap, $
             'r_bkgd', rbkgd, $
             'inbkgd', instat.bkgdmean, $
             'rsdbkgd', rsdstat.bkgdmean, $
             'rsdnet', rsdstat.net $
           )

print, ''
print, 'eval_resid: Aperture size = ',outstruct.r_ap
print, 'eval_resid: Raw counts in input image aperture = ',outstruct.intot
print, 'eval_resid: Net counts in input image aperture = ',outstruct.innet
print, 'eval_resid: Raw counts in resid image aperture = ',outstruct.rsdtot
print, 'eval_resid: Relative residual level = ',outstruct.rsdrat

end
