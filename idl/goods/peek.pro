pro peek,mag=mag,autoscale=autoscale,sigma=sigma

  if n_elements(mag) eq 0 then mag = 2
  if n_elements(sigma) eq 0 then sigma=1.

  print,'peek: using magnification of '+$
    strcompress(string(mag),/remove_all)+' for the zoom'

  mdir='/data/ubu1/goods/mosaics/'

  bimage = mdir+'s1b_mos4_blk.fits'
  iimage = mdir+'s1i_mos4_blk.fits'
  vimage = mdir+'s1v_mos4_blk.fits'
  zimage = mdir+'s1z_mos4_blk.fits'

  print,'peek: reading in the images'
  print,'      ',bimage
  fxread,bimage,bim,bhd
  print,'      ',vimage
;  fxread,vimage,vim,vhd
;  print,'      ',iimage
  fxread,iimage,iim,ihd
  print,'      ',zimage
  fxread,zimage,zim,zhd

  imwin = 10
  window,imwin,xsize=400,ysize=640,title='data window'

  sz=size(zim,/dimen)
  mltx=nint(float(sz[0])/!d.x_size)
  mlty=nint(float(sz[1])/!d.y_size)
  nxy = max([mltx,mlty])

  newsz=sz/nxy
  print,'peek: rebinning and rescaling the viewing images'
  bimsc=rebin(bytscl(bim,min=-0.00175111,max=0.0088532),newsz[0],newsz[1])
;  vimsc=rebin(bytscl(vim,min=-0.00706771,max=0.0413516),newsz[0],newsz[1])
  iimsc=rebin(bytscl(iim,min=-0.00643398,max=0.0376096),newsz[0],newsz[1])
  zimsc=rebin(bytscl(zim,min=-0.00324947,max=0.0222154),newsz[0],newsz[1])

  loadct,0 & invct
  plotimage,zimsc,imgxr=[0,sz[0]-1],imgyr=[0,sz[1]-1]
;  tvscl,zimsc
;  tvscl,vimsc,channel=2
;  tvscl,bimsc,channel=3

  print,'peek: setting up zoom window'
  zoomsize = 200
  zmwin = 11
  window,zmwin,xsize=zoomsize,ysize=zoomsize,title='zoom window'
  
; ** Structure !MOUSE, 4 tags, length=16:
;    X               LONG               511
;    Y               LONG               252
;    BUTTON          LONG                 4
;    TIME            LONG        1428829775
;   X and Y: Contain the location (in device coordinates) of the 
;   cursor when the mouse button was pressed.
;   BUTTON: Contains 
;   - 1 (one) if the left mouse button was pressed,
;   - 2 (two) if the middle mouse button was pressed
;   - 4 (four) if the right mouse button was pressed. 
;   TIME: Contains the number of milliseconds since a base time.

; CURSOR, X, Y [, Wait | [, /CHANGE | , /DOWN | , /NOWAIT | , /UP | ,
; /WAIT]] [, /DATA | , /DEVICE, | , /NORMAL] 

  print,'peek: starting interactivity'
  !mouse.button = 0
  mousepress = !mouse.button
  while mousepress ne 4 do begin
      mousepress = !mouse.button
      wset,imwin
;      tvrdc,xcoord,ycoord,2,/data
      cursor,xcoord,ycoord,0,/data
      xcoord = (xcoord > 0 < (sz[0]-1))
      ycoord = (ycoord > 0 < (sz[1]-1))

; return alpha/dec for cursor position
      xyad,zhd,xcoord,ycoord,alpha,dec
;      print,form="($,2(f7.2),a)",xcoord,ycoord,carrret()
      print,form="($,2(f7.2),f12.8,2x,f12.8,a)",xcoord,ycoord,alpha,dec,carrret()
;      print,form="($,2(f7.2),g,a)",xcoord,ycoord,zim[xcoord,ycoord],carrret()

; determine the size of the window and set the size of the zoom window.
      wset,zmwin
      xmin = nint(xcoord-!d.x_size/2/mag > 0           )
      xmax = nint(xcoord+!d.x_size/2/mag < (sz[0]-1))
      ymin = nint(ycoord-!d.y_size/2/mag > 0           )
      ymax = nint(ycoord+!d.y_size/2/mag < (sz[1]-1))

;      print,!d.x_size,!d.y_size
;      print,!d.x_size/2/mag,!d.y_size/2/mag
;      print,xcoord,ycoord,xmin,xmax,ymin,ymax

;      print,'... updating zoom window ',xmin,xmax,ymin,ymax

      if n_elements(autoscale) ne 0 then begin
          tvscl,rebin(zim[xmin:xmax,ymin:ymax],$
                      nint(xmax-xmin+1)*mag,nint(ymax-ymin+1)*mag),channel=1
          tvscl,rebin(iim[xmin:xmax,ymin:ymax],$
                      nint(xmax-xmin+1)*mag,nint(ymax-ymin+1)*mag),channel=2
          tvscl,rebin(bim[xmin:xmax,ymin:ymax],$
                      nint(xmax-xmin+1)*mag,nint(ymax-ymin+1)*mag),channel=3
      endif else begin 
          tvscl,rebin(bytscl(zim[xmin:xmax,ymin:ymax],min=-0.00324947,max=sigma*0.0222154),$
                      nint(xmax-xmin+1)*mag,nint(ymax-ymin+1)*mag),channel=1
          tvscl,rebin(bytscl(iim[xmin:xmax,ymin:ymax],min=-0.00643398,max=sigma*0.0376096),$
                      nint(xmax-xmin+1)*mag,nint(ymax-ymin+1)*mag),channel=2
          tvscl,rebin(bytscl(bim[xmin:xmax,ymin:ymax],min=-0.00175111,max=sigma*0.00885321),$
                      nint(xmax-xmin+1)*mag,nint(ymax-ymin+1)*mag),channel=3
      endelse 
          
;      if n_elements(autoscale) ne 0 then $
;        tvscl,rebin(zim[xmin:xmax,ymin:ymax],$
;                    nint(xmax-xmin+1)*mag,nint(ymax-ymin+1)*mag) else $
;        tvscl,rebin(bytscl(zim[xmin:xmax,ymin:ymax],min=-0.00324947,max=sigma*0.0222154),$
;                           nint(xmax-xmin+1)*mag,nint(ymax-ymin+1)*mag)

      wset,imwin
      plotimage,zimsc,imgxr=[0,sz[0]-1],imgyr=[0,sz[1]-1],/nodata,/noerase
  endwhile
  print,' '
  print,'peek: all done with peek'

end


;      tvscl,rebin(zimsc[xmin:xmax,ymin:ymax],$
;                  nint(xmax-xmin+1)*mag,nint(ymax-ymin+1)*mag),channel=1
;      tvscl,rebin(rimsc[xmin:xmax,ymin:ymax],$
;                  nint(xmax-xmin+1)*mag,nint(ymax-ymin+1)*mag),channel=2
;      tvscl,rebin(bimsc[xmin:xmax,ymin:ymax],$
;                  nint(xmax-xmin+1)*mag,nint(ymax-ymin+1)*mag),channel=3

;      tvcircle,10,!d.x_size/2,!d.y_size/2 ;,color=-1


;  bimage = mdir+'s1b_mos1_drz.fits'
;  iimage = mdir+'s1i_mos1_drz.fits'
;  vimage = mdir+'s1v_mos1_drz.fits'
;  zimage = mdir+'s1z_mos1_drz.fits'
;   bimsc=bytscl(bim,min=-0.00175111,max=0.00885321)
;   vimsc=bytscl(vim,min=-0.00706771,max=0.0413516)
;   iimsc=bytscl(iim,min=-0.00643398,max=0.0376096)
;   zimsc=bytscl(zim,min=-0.00324947,max=0.0222154)
;   xxs = float((size(zimsc,/dimen))[0])/!d.x_size
;   yys = float((size(zimsc,/dimen))[1])/!d.y_size
;   tvscl,rebin(zimsc,sz[0]/nxy,sz[1]/nxy),channel=1
;   tvscl,rebin(vimsc,sz[0]/nxy,sz[1]/nxy),channel=2
;   tvscl,rebin(bimsc,sz[0]/nxy,sz[1]/nxy),channel=3

;      x0=10000l & y0=10000l
;      dsz=200l
;      
;      fxread,bimage,bim,bhd,x0-dsz,x0+dsz,y0-dsz,y0+dsz
;      fxread,vimage,vim,vhd,x0-dsz,x0+dsz,y0-dsz,y0+dsz
;      fxread,iimage,iim,ihd,x0-dsz,x0+dsz,y0-dsz,y0+dsz
;      fxread,zimage,zim,zhd,x0-dsz,x0+dsz,y0-dsz,y0+dsz

