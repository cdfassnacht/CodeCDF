pro deg2hms, infile

;
; Procedure deg2hms
;
; Description: Converts a list of RA and Dec, each listed in decimal
;               degrees, into the standard h:m:s format.
;
; Inputs: infile      (string)     name of input file.  The file is
;                                   assumed to have (at least) three
;                                   columns, with the first column being the
;                                   id number, the second the RA, and the
;                                   third the Dec.
;
; Revision history:
;  2003Apr22 Chris Fassnacht -- First working version.
;

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, $
      'syntax: deg2hms, infile'
    print, ''
    return
endif

; Read in input file with readcol

print, ''
print, 'Reading catalog from ',infile
readcol, infile, id, radeg, decdeg, format='l,d,d'

; Set up string arrays

neid = n_elements(id)
sra = strarr(neid)
sdec = strarr(neid)

; Run the loop

for i=0, neid-1 do begin 
    sra(i)  = dec2hms(radeg(i) / 15.0d)
    sdec(i) = dec2hms(decdeg(i))
    print, id(i), sra(i), sdec(i)
endfor
        
end 


