pro er_phot, infile, xlens, ylens, er, outstat, dodisplay=dodisplay
;
; Procedure er_phot
;
; Description: A driver to call circ_phot for each of the four bands for
;               a GOODS ACS object.
;
; Inputs: infile     (string)      base name of input FITS file (root+id number)
;         xlens      (float)       x position of the lens, in pixels
;         ylens      (float)       y position of the lens, in pixels
;         er         (float)       Einstein ring radius, in pixels
;         outstat    (structure)   Output structure for image info
;         [/dodisplay]             display images if set
;
;
; Revision history:
;  2002Apr28 Chris Fassnacht -- First working version.
;

; Check input format

if n_params() lt 4 then begin
    print, ''
    print, 'syntax: er_phot, infile, xlens, ylens, er ,outstat [,/dodisplay]'
    print, ''
    return
endif

; Initialize some variables and set optional variables to defaults if
;  not set by the function call

if not keyword_set(dodisplay) then dodisplay = 0

; Set zeropoints

bm0=25.65288d0
vm0=26.49341d0
im0=25.64053d0
zm0=24.84315d0

; Set directory names

bdir = 'b-band/'
vdir = 'v-band/'
idir = 'i-band/'
zdir = 'z-band/'
catdir = 'Stampcat/'
residdir = 'Resid/'
fitdir = 'Fit/'
maskdir = 'Mask/'

; Set filenames

bfits = bdir+infile+'_b.fits'
vfits = vdir+infile+'_v.fits'
ifits = idir+infile+'_i.fits'
zfits = zdir+infile+'_z.fits'

; Call circ_phot for each of the bands

rbkgd = [60, 100]

circ_phot, bfits,xlens,ylens,er,rbkgd,bstat,zp=bm0,dodisplay=dodisplay
circ_phot, vfits,xlens,ylens,er,rbkgd,vstat,zp=vm0,dodisplay=dodisplay
circ_phot, ifits,xlens,ylens,er,rbkgd,istat,zp=im0,dodisplay=dodisplay
circ_phot, zfits,xlens,ylens,er,rbkgd,zstat,zp=zm0,dodisplay=dodisplay

; Create the structure to contain important info

outstat = create_struct('object',infile,$
                        'xlens',xlens,$
                        'ylens',ylens,$
                        'theta_e',er,$
                        'b_e',bstat.m_phot,$
                        'v_e',vstat.m_phot,$
                        'i_e',istat.m_phot,$
                        'z_e',zstat.m_phot)

end
