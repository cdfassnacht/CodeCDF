; cull information from an input list of raw wirc files
;
; some of this information is needed for future workflow
;
pro wirc_logsheet,inlist,$
                  OUTFILE=outfile ; default is logsheet.txt

if (keyword_set(OUTFILE) NE 1) then outfile='logsheet.txt'

readcol,inlist,infile,format="A"
nfiles=n_elements(infile)


openw,1,outfile


for i=0,nfiles-1 do begin

head=headfits(infile[i])

    filter1=strcompress(sxpar(head,'FORE'),/remove_all)
    filter2=strcompress(sxpar(head,'AFT'),/remove_all)
    ra=sxpar(head,'RA')
    dec=sxpar(head,'DEC')
    airmass=sxpar(head,'AIRMASS')
    exptime=sxpar(head,'EXPTIME')
    object=sxpar(head,'OBJECT')
    obstype=sxpar(head,'OBSTYPE')
    time=sxpar(head,'TIME')
    coadds=sxpar(head,'COADDS')


printf,1,format='(I3," ",A-18,A22,A10,A13,A13,"   ",F4.1,A16,A16,"    ",I2,F8.2)',i,infile[i],object,obstype,ra,dec,airmass,filter1,filter2,coadds,exptime

endfor

free_lun,1


end
