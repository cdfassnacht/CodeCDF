; Some diagnostic output
;  shows the time history and calibration statistics
;  calculated by wirc_fixwcs
;
pro wirc_medzpt,inlist



readcol,inlist,infile,format="A"
nfiles=n_elements(infile)
magzpt=fltarr(nfiles)
time=fltarr(nfiles)
crval1=fltarr(nfiles)
crval2=fltarr(nfiles)
pnterr=fltarr(nfiles)
raoff=fltarr(nfiles)
decoff=fltarr(nfiles)

 ; load the zeropoint array
 print,'Loading magnitude zeropoints.'
 for i=0,nfiles-1 do begin
     hd=headfits(infile[i])
     magzpt[i]=sxpar(hd,'MAGZPT')
     crval1[i]=sxpar(hd,'CRVAL1')
     crval2[i]=sxpar(hd,'CRVAL2')
     pnterr[i]=sxpar(hd,'PNTERR')
     raoff[i]=sxpar(hd,'RAOFF')
     decoff[i]=sxpar(hd,'RAOFF')
     lst=sxpar(hd,'LST')
     get_coords,timetemp,instring=lst+' 0:0:0'
     time[i]=timetemp[0]
     
 endfor
 
 medianzpt=median(magzpt)
 print,'Median magnitude zeropoint is ',medianzpt
 
 k=moment(magzpt,/nan)
 print,'Average is ',k[0],' +/-',sqrt(k[1])
 
 print,'Plotting median magnitude zeropoint vs. LST.'
 window,xsize=640,ysize=480
 plot,time,magzpt,yrange=[medianzpt-1.,medianzpt+1.],psym=2
 
 print,'Plotting mean RA offset vs. LST'
 plot,time,raoff,psym=2 
 
  print,'Plotting mean DEC offset vs. LST'
 plot,time,decoff,psym=2 
 
 print,'Plotting mean pointing error vs. LST.'
 plot,time,pnterr,psym=2
 
 
 stop
 
 end