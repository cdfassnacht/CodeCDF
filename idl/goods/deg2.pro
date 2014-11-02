pro deg2, idfile, catfile

;
; Procedure deg2
;
; Description: Converts a list of RA and Dec, each listed in decimal
;               degrees, into the standard h:m:s format.
;
; Inputs: idfile      (string)     name of input file.  The file is
;                                   assumed to have (at least) three
;                                   columns, with the first column being the
;                                   id number, the second the RA, and the
;                                   third the Dec.
;
; Revision history:
;  2003Apr22 Chris Fassnacht -- First working version.
;

; Check input format

if n_params() lt 2 then begin
    print, ''
    print, $
      'syntax: deg2hms, infile, catfile'
    print, ''
    return
endif

; Read in input files with readcol

print, ''
print, 'Reading id list from ',idfile
readcol, idfile, id, format='l'
print,''
print, 'Reading catalog from ',catfile
readcol, catfile, idc,xc,yc,alpha,delta,format='l,f,f,d,d'

; Set up arrays

neid = n_elements(id)
ra = dblarr(neid)
dec = dblarr(neid)
sra = strarr(neid)
sdec = strarr(neid)
idcheck = lonarr(neid)

; Run the loop

print,''
for i=0, neid-1 do begin 
    idcheck[i] = idc[id[i]-1]
    ra[i] = alpha[id[i]-1]
    dec[i] = delta[id[i]-1]
    sra[i]  = dec2hms(ra[i] / 15.0d)
    sdec[i] = dec2hms(dec[i])
    print, id[i], idcheck[i], sra[i], sdec[i]
endfor

print,''
forprint, id, idcheck, sra, sdec, text='hms.out'

end 


