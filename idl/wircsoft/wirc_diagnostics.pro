; Some diagnostic output
;  shows the time history and calibration statistics
;  calculated by wirc_fixwcs
;
pro wirc_diagnostics,inlist,$
                     PLOTROOT=PLOTROOT ; rootname for storing plotfiles



readcol,inlist,infile,format="A"
nfiles=n_elements(infile)
magzpt=fltarr(nfiles)
time=fltarr(nfiles)
crval1=fltarr(nfiles)
crval2=fltarr(nfiles)
pnterr=fltarr(nfiles)
raoff=fltarr(nfiles)
decoff=fltarr(nfiles)
pntstars=fltarr(nfiles)
seeing=fltarr(nfiles)
rawmed=fltarr(nfiles)

 ; load the zeropoint array
 print,'Loading magnitude zeropoints.'
 for i=0,nfiles-1 do begin
     hd=headfits(infile[i])
     magzpt[i]=sxpar(hd,'MAGZPT')
     crval1[i]=sxpar(hd,'CRVAL1')
     crval2[i]=sxpar(hd,'CRVAL2')
     pnterr[i]=sxpar(hd,'PNTERR')
     raoff[i]=sxpar(hd,'RAOFF')
     decoff[i]=sxpar(hd,'DECOFF')
     pntstars[i]=sxpar(hd,'PNTSTARS')
     seeing[i]=sxpar(hd,'SEEING')
     lst=sxpar(hd,'LST')
     rawmed[i]=sxpar(hd,'RAWMED')
     get_coords,timetemp,instring=lst+' 0:0:0'
     time[i]=timetemp[0]
     
 endfor
 
 medianzpt=median(magzpt)
 print,'Median magnitude zeropoint is ',medianzpt
 
 k=moment(magzpt,/nan)
 print,'Average is ',k[0],' +/-',sqrt(k[1])
 
 print,'Plotting median magnitude zeropoint vs. LST.'
 window,xsize=640,ysize=480
 plot,time,magzpt,yrange=[min(magzpt)-3.*sqrt(k[1]),max(magzpt)+3.*sqrt(k[1])],psym=2,xtitle='LST',ytitle='Magnitude Zeropoint'
 
 
 if keyword_set (PLOTROOT) then begin
     outname=strcompress(PLOTROOT+'_magzpt.gif',/REMOVE_ALL)
     write_gif,outname,tvrd()
    endif
 
 print,'Press any key to continue...'
 result=get_kbrd()
 
print,'Plotting mean RA offset vs. LST'
plot,time,raoff,psym=2,xtitle='LST',ytitle='RA Offset (arcsec)'
 
if keyword_set (PLOTROOT) then begin
   outname=strcompress(PLOTROOT+'_raoffset.gif',/REMOVE_ALL)
   write_gif,outname,tvrd()
 endif
 
 print,'Press any key to continue...'
 result=get_kbrd()
 
print,'Plotting mean DEC offset vs. LST'
plot,time,decoff,psym=2,xtitle='LST',ytitle='DEC Offset (arcsec)'
 
if keyword_set (PLOTROOT) then begin
   outname=strcompress(PLOTROOT+'_decoffset.gif',/REMOVE_ALL)
   write_gif,outname,tvrd()
  endif
 
 
 print,'Press any key to continue...'
 result=get_kbrd()
 
print,'Plotting mean pointing error vs. LST.'
plot,time,pnterr,psym=2,xtitle='LST',ytitle='Mean Pointing Error (arcsec)'
 
if keyword_set (PLOTROOT) then begin
     outname=strcompress(PLOTROOT+'_pnterr.gif',/REMOVE_ALL)
     write_gif,outname,tvrd()
  endif
 
  print,'Press any key to continue...'
  result=get_kbrd()

  
print,'Plotting median seeing vs. LST.'
plot,time,seeing,psym=2,xtitle='LST',ytitle='Seeing (arcsec)'
 
print,'Median seeing was ',median(seeing),' arcseconds.'
k=moment(seeing)
print,'Average was ',k[0],'+/-',sqrt(k[1]),' arcseconds.'
 
if keyword_set (PLOTROOT) then begin
     outname=strcompress(PLOTROOT+'_seeing.gif',/REMOVE_ALL)
     write_gif,outname,tvrd()
  endif
    
    
  print,'Press any key to continue...'
  result=get_kbrd()    
    
print,'Plotting median background vs. LST'
plot,time,rawmed,psym=2,xtitle='LST',ytitle='Median Raw Counts (DN)'

if keyword_set (PLOTROOT) then begin
     outname=strcompress(PLOTROOT+'_rawcounts.gif',/REMOVE_ALL)
     write_gif,outname,tvrd()
  endif  
    
    
    
 end