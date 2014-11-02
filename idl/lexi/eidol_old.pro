pro eidol,xx,yy,$
          mosaic=mosaic,datloc=datloc,scale=scale,$
;          color=color,band=band,$
          title=title,$
          minsig=minsig,maxsig=maxsig,$
          size=size,pxsc=pxsc,$
          annotate=annotate,$
          fits=fits,postscript=postscript,$
          segmentation=segmentation,$
          unsharpmask=unsharpmask,$
          xoffset=xoffset,yoffset=yoffset

    if n_elements(datloc) eq 0 then datloc='/data/raid1/acs_archive/mosaics/'
    if n_elements(mosaic) eq 0 then mosaic='s123'
    if n_elements(scale) eq 0 then scale='1' else scale=strcompress(scale,/remove_all)
    if n_elements(band) eq 0 then band ='z'

    if n_elements(minsig) eq 0 then minsig=1.5
    if n_elements(maxsig) eq 0 then maxsig=15.
    if n_elements(pxsc) eq 0 then pxsc=0.05d0
    if n_elements(size) eq 0 then size=100l

    if n_elements(xoffset) eq 0 then xoffset=0.d0
    if n_elements(yoffset) eq 0 then yoffset=0.d0
    
    print,'eidol: output will be a size of '+$
      strcompress(string(size),/remove_all)+' pixels, or '+$
      strcompress(string(size*pxsc,format='(f6.2)'),/remove_all)+' arcsec'

    if n_elements(xx) eq 0 or n_elements(yy) eq 0 then begin
        print,'need to specify coordinates, as ra and dec'
        return
    endif 

    if size(xx,/type) eq 7 then begin
        rahms=xx
        xx=hms2dec(xx)*15.d0-xoffset
    endif else rahms=dec2hms(xx/15.d0-xoffset)
    if size(yy,/type) eq 7 then begin
        dechms=yy
        yy=hms2dec(yy)*1.d0-yoffset
    endif else dechms=dec2hms(yy-yoffset)
  
    if n_elements(title) eq 0 then $
      title = 'GDS_'+$
      strcompress(string(xx),/remove_all)+$
      strcompress(string(yy),/remove_all)

    if n_elements(postscript) ne 0 then begin
        ps_open,title,/color,/portrait;,/ps_fonts
        device,xsize=8,ysize=8,/inches,xoffset=0.25,yoffset=2.75,/portrait,/times
;       XSIZE           FLOAT           8.00000
;       XOFF            FLOAT          0.250000
;       YSIZE           FLOAT           8.00000
;       YOFF            FLOAT           2.75000
;       FILENAME        STRING    '/home/leonidas/nwork/goods/catalogs/testpsfor'...
;       INCHES          INT              1
;       COLOR           INT              1
;       BITS_PER_PIXEL  INT              4
    endif
  
    ra = xx
    dec = yy

    bfits=datloc+'s1b_mos_scale'+scale+'_drz.fits'
    vfits=datloc+mosaic+'v_mos_scale'+scale+'_drz.fits'
    ifits=datloc+mosaic+'i_mos_scale'+scale+'_drz.fits'
    zfits=datloc+mosaic+'z_mos_scale'+scale+'_drz.fits'

    bim=eidol_gim(bfits,ra,dec,size)
    vim=eidol_gim(vfits,ra,dec,size)
    iim=eidol_gim(ifits,ra,dec,size)
    zim=eidol_gim(zfits,ra,dec,size)

;    bimsm=gaussconv(bim,17,50)-bim
;    vimsm=gaussconv(vim,17,50)-vim
;    iimsm=gaussconv(iim,17,50)-iim
;    zimsm=gaussconv(zim,17,50)-zim
        
    sz=size(zim,/dimen)
    literstat,bim,bst,/silent
    literstat,vim,vst,/silent
    literstat,iim,ist,/silent
    literstat,zim,zst,/silent

;    literstat,bimsm,bsmst,/silent
;    literstat,vimsm,vsmst,/silent
;    literstat,iimsm,ismst,/silent
;    literstat,zimsm,zsmst,/silent

    bimsc=bytscl(bim,$
                 min=bst.median-minsig*bst.sigma,$
                 max=bst.median+maxsig*bst.sigma)
    vimsc=bytscl(vim,$
                 min=vst.median-minsig*vst.sigma,$
                 max=vst.median+maxsig*vst.sigma)
    iimsc=bytscl(iim,$
                 min=ist.median-minsig*ist.sigma,$
                 max=ist.median+maxsig*ist.sigma)
    zimsc=bytscl(zim,$
                 min=zst.median-minsig*zst.sigma,$
                 max=zst.median+maxsig*zst.sigma)
    
    bvi=bytarr(3,sz[0],sz[1])
    bvi[0,*,*]=bytscl(iim,$
                      min=ist.median-minsig*ist.sigma,$
                      max=ist.median+maxsig*ist.sigma)
    bvi[1,*,*]=bytscl(vim,$
                      min=vst.median-minsig*vst.sigma,$
                      max=vst.median+maxsig*vst.sigma)
    bvi[2,*,*]=bytscl(bim,$
                      min=bst.median-minsig*bst.sigma,$
                      max=bst.median+maxsig*bst.sigma)
    
    bvz=bytarr(3,sz[0],sz[1])
    bvz[0,*,*]=bytscl(zim,$
                      min=zst.median-minsig*zst.sigma,$
                      max=zst.median+maxsig*zst.sigma)
    bvz[1,*,*]=bytscl(vim,$
                      min=vst.median-minsig*vst.sigma,$
                      max=vst.median+maxsig*vst.sigma)
    bvz[0,*,*]=bytscl(bim,$
                      min=bst.median-minsig*bst.sigma,$
                      max=bst.median+maxsig*bst.sigma)

    biz=bytarr(3,sz[0],sz[1])
    biz[0,*,*]=bytscl(zim,$
                      min=zst.median-minsig*zst.sigma,$
                      max=zst.median+maxsig*zst.sigma)
    biz[1,*,*]=bytscl(iim,$
                      min=ist.median-minsig*ist.sigma,$
                      max=ist.median+maxsig*ist.sigma)
    biz[0,*,*]=bytscl(bim,$
                      min=bst.median-minsig*bst.sigma,$
                      max=bst.median+maxsig*bst.sigma)

    viz=bytarr(3,sz[0],sz[1])
    viz[0,*,*]=bytscl(zim,$
                      min=zst.median-minsig*zst.sigma,$
                      max=zst.median+maxsig*zst.sigma)
    viz[1,*,*]=bytscl(iim,$
                      min=ist.median-minsig*ist.sigma,$
                      max=ist.median+maxsig*ist.sigma)
    viz[2,*,*]=bytscl(vim,$
                      min=vst.median-minsig*vst.sigma,$
                      max=vst.median+maxsig*vst.sigma)

;    vizsm=bytarr(3,sz[0],sz[1])
;    vizsm[0,*,*]=bytscl(zimsm,$
;                        min=zsmst.median-minsig*zsmst.sigma,$
;                        max=zsmst.median+maxsig*zsmst.sigma)
;    vizsm[1,*,*]=bytscl(iimsm,$
;                        min=ismst.median-minsig*ismst.sigma,$
;                        max=ismst.median+maxsig*ismst.sigma)
;    vizsm[2,*,*]=bytscl(vimsm,$
;                        min=vsmst.median-minsig*vsmst.sigma,$
;                        max=vsmst.median+maxsig*vsmst.sigma)
    
;      position=[0.1,0.55,0.45,0.9],title=title,$

    subcellarray, [1,1,1,1],[1,1], newpan, newsubpan
    
    plotimage,bimsc,/preserve,_extra=_extra,panel=newpan[0,1,*],subpan=newsubpan[0,1,*],$
      imgxr=[-1,1]*size*pxsc/2.,imgyr=[-1,1]*size*pxsc/2.,xtit='arcsec',ytit='arcsec',tit=title
    xyouts,-0.9*size*pxsc/2.,0.8*size*pxsc/2.,'B',/data
;    xyouts,-0.9*size*pxsc/2.,-0.8*size*pxsc/2.,title,/data
    plotimage,vimsc,/preserve,_extra=_extra,panel=newpan[1,1,*],subpan=newsubpan[1,1,*],$
      imgxr=[-1,1]*size*pxsc/2.,imgyr=[-1,1]*size*pxsc/2.,xtit='arcsec',/noerase,ytickn=replicate(' ',20),tit=rahms
    xyouts,-0.9*size*pxsc/2.,0.8*size*pxsc/2.,'V',/data
;    xyouts,-0.9*size*pxsc/2.,-0.8*size*pxsc/2.,rahms+dechms,/data
    plotimage,iimsc,/preserve,_extra=_extra,panel=newpan[2,1,*],subpan=newsubpan[2,1,*],$
      imgxr=[-1,1]*size*pxsc/2.,imgyr=[-1,1]*size*pxsc/2.,xtit='arcsec',/noerase,ytickn=replicate(' ',20),tit=dechms
    xyouts,-0.9*size*pxsc/2.,0.8*size*pxsc/2.,'I',/data
    plotimage,zimsc,/preserve,_extra=_extra,panel=newpan[3,1,*],subpan=newsubpan[3,1,*],$
      imgxr=[-1,1]*size*pxsc/2.,imgyr=[-1,1]*size*pxsc/2.,xtit='arcsec',/noerase,ytickn=replicate(' ',20)
    xyouts,-0.9*size*pxsc/2.,0.8*size*pxsc/2.,'Z',/data

    plotimage,bvi,/preserve,_extra=_extra,panel=newpan[0,0,*],subpan=newsubpan[0,0,*],$
      imgxr=[-1,1]*size*pxsc/2.,imgyr=[-1,1]*size*pxsc/2.,xtit='arcsec',ytit='arcsec',/noerase
    xyouts,-0.9*size*pxsc/2.,0.8*size*pxsc/2.,'BVI',/data
    plotimage,bvz,/preserve,_extra=_extra,panel=newpan[1,0,*],subpan=newsubpan[1,0,*],$
      imgxr=[-1,1]*size*pxsc/2.,imgyr=[-1,1]*size*pxsc/2.,xtit='arcsec',/noerase,ytickn=replicate(' ',20)
    xyouts,-0.9*size*pxsc/2.,0.8*size*pxsc/2.,'BVZ',/data
    plotimage,biz,/preserve,_extra=_extra,panel=newpan[2,0,*],subpan=newsubpan[2,0,*],$
      imgxr=[-1,1]*size*pxsc/2.,imgyr=[-1,1]*size*pxsc/2.,xtit='arcsec',/noerase,ytickn=replicate(' ',20)
    xyouts,-0.9*size*pxsc/2.,0.8*size*pxsc/2.,'BIZ',/data
    plotimage,viz,/preserve,_extra=_extra,panel=newpan[3,0,*],subpan=newsubpan[3,0,*],$
      imgxr=[-1,1]*size*pxsc/2.,imgyr=[-1,1]*size*pxsc/2.,xtit='arcsec',/noerase,ytickn=replicate(' ',20)
    xyouts,-0.9*size*pxsc/2.,0.8*size*pxsc/2.,'VIZ',/data

;    plotimage,vizsm,/preserve,_extra=_extra,panel=newpan[3,0,*],subpan=newsubpan[3,0,*],$
;      imgxr=[-1,1]*size*pxsc/2.,imgyr=[-1,1]*size/2.*pxsc/2.,xtit='arcsec',/noerase,ytickn=replicate(' ',20)

;    plotpan,[-1,1]*size*pxsc/2.,[-1,1]*size*pxsc/2.,/nodata,panel=newpan[3,0,*],subpan=newsubpan[3,0,*],$
;      /noerase,xtickn=replicate(' ',20),ytickn=replicate(' ',20)

    if keyword_set(annotate) then begin 
        plot,[-1,1]*size*pxsc/2.,imgyr=[-1,1]*size*pxsc/2.,xtit='arcsec',ytit='arcsec',/noerase,$
          /xst,/yst,/nodata
        srcor,$
          scat.alpha_j2000,scat.delta_j2000,$
          [ra0],[dec0],7.0,pix1,pix2,spherical=2
        if total(pix1) gt 0 then begin
            tra=(scat[pix1].alpha_j2000-ra0[0])*3600.d0*cos(dec0[0]/!radeg)
            tdec=(scat[pix1].delta_j2000-dec0[0])*3600.d0
            for j=0,n_elements(pix1)-1 do begin
                ncrd=[tra[j],tdec[j]]#rotmat
                xyouts,-ncrd[0]+0.8,ncrd[1]-0.2,strcompress(long(scat[pix1[j]].number)),$
                  /data,charsize=1.3,color=redidx
            endfor
        endif
    endif
        
    if n_elements(postscript) ne 0 then ps_close

end 
