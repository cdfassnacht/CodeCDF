function acsproc_datafile, subexposure, cutout=cutout, rectify=rectify, combine=combine, mask=mask, photmod=photmod, biz=biz, old=old, test=test, par=par, hst_root=hst_root, ext=ext, subext=subext, eps=eps 
  
  common acsproc_photmod_rhalf_data, acs_filter, dataroot

  if keyword_set(ext) then begin
    ftype = ''
    extension = '.' + ext 
  endif else begin
    ftype = '_flt'
    extension = '.fits'
  endelse
  
  if keyword_set(subext) then subext = '_' + subext else subext = ''

  if keyword_set(mask) then extension = '_mask' + extension $
  else if keyword_set(old) and keyword_set(photmod) then extension = '_photmod_old' + extension $
  else if keyword_set(test) and keyword_set(photmod) then extension = '_photmod_test' + extension $
  else if keyword_set(par) and keyword_set(photmod) then extension = '_photmod_par' + extension $
  else if keyword_set(eps) and keyword_set(photmod) then extension = '_photmod' $
  else if keyword_set(photmod) then extension = '_photmod' + extension $
  else if keyword_set(biz) then extension = '_biz' + extension
  if keyword_set(old) and keyword_set(photmod) then old=0L

  if keyword_set(mask) or keyword_set(photmod) or keyword_set(biz) then combine = 1L

  if keyword_set(hst_root) then begin
    datafile = dataroot + subexposure + extension
  endif else if keyword_set(combine) then begin
    filter = acs_filter[subexposure[0].filtcode]
    datafile = dataroot + subexposure[0].visit + '/'+ subexposure[0].uniqname + '_' + $
               strtrim(string(subexposure[0].proposid, format='(i)'),2) + '_' + $
               strtrim(string(subexposure[0].expstart, format='(i)'),2) + '_' + $
               filter + '_' +  strtrim(subexposure[0].nexp,2) + subext + extension 
    if keyword_set(old) then datafile = dataroot  + 'old/'+ subexposure[0].uniqname + '_' + $
               strtrim(string(subexposure[0].proposid, format='(i)'),2) + '_' + $
               strtrim(string(subexposure[0].expstart, format='(i)'),2) + '_' + $
               filter + '_' +  strtrim(subexposure[0].nexp,2) + subext + extension 
  endif else begin
    datafile = dataroot + subexposure.visit + '/'+ subexposure.rootname + ftype
    if keyword_set(cutout) then datafile +=  '_cutout' + subext + extension $
    else if keyword_set(rectify) then datafile += '_recti' + subext + extension
  endelse
  
  return, datafile              
  
end

;---------------------------------------------------------------------
function acsproc_load_maskstruc, fits=fits
    
  if keyword_set(fits) then begin
    moffat = {modimID:fits.moffat_modimID, name:fits.moffat_name, modim:fits.moffat_modim, pars:fits.moffat_pars, pars_err:fits.moffat_pars_err, chi2:fits.moffat_chi2}
    bsplinecen = {modimID:fits.bsplinecen_modimID, name:fits.bsplinecen_name, modim:fits.bsplinecen_modim, dmodim:fits.bsplinecen_dmodim, pars:fits.bsplinecen_pars, pars_err:fits.bsplinecen_pars_err, chi2:fits.bsplinecen_chi2, nthetaID:fits.bsplinecen_nthetaID, rbkptID:fits.bsplinecen_rbkptID}
    bsplinefit = [{modimID:fits.bsplinecen_modimID, name:fits.bsplinecen_name, modim:fits.bsplinecen_modim, chi2:fits.bsplinecen_chi2, nthetaID:fits.bsplinecen_nthetaID, rbkptID:fits.bsplinecen_rbkptID}, $
                  {modimID:fits.bsplinefit_02_modimID, name:fits.bsplinefit_02_name, modim:fits.bsplinefit_02_modim, chi2:fits.bsplinefit_02_chi2, nthetaID:fits.bsplinefit_02_nthetaID, rbkptID:fits.bsplinefit_02_rbkptID}, $
                  {modimID:fits.bsplinefit_04_modimID, name:fits.bsplinefit_04_name, modim:fits.bsplinefit_04_modim, chi2:fits.bsplinefit_04_chi2, nthetaID:fits.bsplinefit_04_nthetaID, rbkptID:fits.bsplinefit_04_rbkptID}, $
                  {modimID:fits.bsplinefit_024_modimID, name:fits.bsplinefit_024_name, modim:fits.bsplinefit_024_modim, chi2:fits.bsplinefit_024_chi2, nthetaID:fits.bsplinefit_024_nthetaID, rbkptID:fits.bsplinefit_024_rbkptID}, $
                  {modimID:fits.bsplinefit_01_modimID, name:fits.bsplinefit_01_name, modim:fits.bsplinefit_01_modim, chi2:fits.bsplinefit_01_chi2, nthetaID:fits.bsplinefit_01_nthetaID, rbkptID:fits.bsplinefit_01_rbkptID}, $
                  {modimID:fits.bsplinefit_012_modimID, name:fits.bsplinefit_012_name, modim:fits.bsplinefit_012_modim, chi2:fits.bsplinefit_012_chi2, nthetaID:fits.bsplinefit_012_nthetaID, rbkptID:fits.bsplinefit_012_rbkptID}, $
                  {modimID:fits.bsplinefit_014_modimID, name:fits.bsplinefit_014_name, modim:fits.bsplinefit_014_modim, chi2:fits.bsplinefit_014_chi2, nthetaID:fits.bsplinefit_014_nthetaID, rbkptID:fits.bsplinefit_014_rbkptID}, $
                  {modimID:fits.bsplinefit_0124_modimID, name:fits.bsplinefit_0124_name, modim:fits.bsplinefit_0124_modim, chi2:fits.bsplinefit_0124_chi2, nthetaID:fits.bsplinefit_0124_nthetaID, rbkptID:fits.bsplinefit_0124_rbkptID}]
    maskstruc = {modimID:fits.modimID, xc:0L, yc:0L, has_modim:fits.has_modim, has_bizzle:fits.has_bizzle, has_photmod:fits.has_photmod, nthetaID:fits.nthetaID, lens_status:fits.lens_status, fmask:fits.fmask, jmask:fits.jmask, moffat:moffat, bsplinecen:bsplinecen, bsplinefit:bsplinefit}
  endif else begin
    mask = bytarr(imsize,imsize)
    modim = fltarr(2*hw+1,2*hw+1)
    moffat = {modimID:0L, name:'Moffat Fit', modim:modim, pars:fltarr(7), pars_err:fltarr(7), chi2:0.0}
    bsplinecen = {modimID:1L, name:'Bspline Cen', modim:modim, dmodim:modim, pars:fltarr(4), pars_err:fltarr(4), chi2:0.0, nthetaID:0L, rbkptID:0L}
    bsplinefit = replicate({modimID:2L, name:'Bspline Fit', modim:modim, pars:fltarr(4), pars_err:fltarr(4), chi2:0.0, nthetaID:0L, rbkptID:0L},8)
    maskstruc = {modimID:0L, xc:0L, yc:0L, nthetaID:0L, lens_status:0L, has_modim:0L, has_bizzle:0L, has_photmod:0L, fmask:mask, jmask:mask, moffat:moffat, bsplinecen:bsplinecen, bsplinefit:bsplinefit}
  endelse
  
  return, maskstruc
  
 end
;---------------------------------------------------------------------
function acsproc_get_maskstruc, maskfile, header=mask_header

  if file_test(maskfile) then begin
    maskstruc =  struct_addtags(acsproc_load_maskstruc(fits=mrdfits(maskfile, 1, mask_header, /silent)), {saved:1L})
  endif else begin
    maskstruc = struct_addtags(acsproc_load_maskstruc(), {saved:0L})
  endelse
  
  return, maskstruc
      
 end
;---------------------------------------------------------------------
pro acsproc_maskstruc_write, maskstruc, maskfile, header=mask_header

  ; flatten the maskstruc so it can be saved as a fits file.
    
  flat_maskstruc = {modimID:maskstruc.modimID, xc:maskstruc.xc, yc:maskstruc.yc, nthetaID:maskstruc.nthetaID, lens_status:maskstruc.lens_status, has_modim:maskstruc.has_modim, has_bizzle:maskstruc.has_bizzle, has_photmod:maskstruc.has_photmod, fmask:maskstruc.fmask, jmask:maskstruc.jmask, $
    moffat_modimID:maskstruc.moffat.modimID, moffat_name:maskstruc.moffat.name, moffat_modim:maskstruc.moffat.modim, moffat_pars:maskstruc.moffat.pars, moffat_pars_err:maskstruc.moffat.pars_err, moffat_chi2:maskstruc.moffat.chi2, $
    bsplinecen_modimID:maskstruc.bsplinecen.modimID, bsplinecen_name:maskstruc.bsplinecen.name, bsplinecen_modim:maskstruc.bsplinecen.modim, bsplinecen_dmodim:maskstruc.bsplinecen.dmodim, bsplinecen_pars:maskstruc.bsplinecen.pars, bsplinecen_pars_err:maskstruc.bsplinecen.pars_err, bsplinecen_chi2:maskstruc.bsplinecen.chi2, bsplinecen_nthetaID:maskstruc.bsplinecen.nthetaID, bsplinecen_rbkptID:maskstruc.bsplinecen.rbkptID, $
    bsplinefit_02_modimID:maskstruc.bsplinefit[1].modimID, bsplinefit_02_name:maskstruc.bsplinefit[1].name, bsplinefit_02_modim:maskstruc.bsplinefit[1].modim, bsplinefit_02_chi2:maskstruc.bsplinefit[1].chi2, bsplinefit_02_nthetaID:maskstruc.bsplinefit[1].nthetaID, bsplinefit_02_rbkptID:maskstruc.bsplinefit[1].rbkptID, $
    bsplinefit_04_modimID:maskstruc.bsplinefit[2].modimID, bsplinefit_04_name:maskstruc.bsplinefit[2].name, bsplinefit_04_modim:maskstruc.bsplinefit[2].modim, bsplinefit_04_chi2:maskstruc.bsplinefit[2].chi2, bsplinefit_04_nthetaID:maskstruc.bsplinefit[2].nthetaID, bsplinefit_04_rbkptID:maskstruc.bsplinefit[2].rbkptID, $
    bsplinefit_024_modimID:maskstruc.bsplinefit[3].modimID, bsplinefit_024_name:maskstruc.bsplinefit[3].name, bsplinefit_024_modim:maskstruc.bsplinefit[3].modim, bsplinefit_024_chi2:maskstruc.bsplinefit[3].chi2, bsplinefit_024_nthetaID:maskstruc.bsplinefit[3].nthetaID, bsplinefit_024_rbkptID:maskstruc.bsplinefit[3].rbkptID, $
    bsplinefit_01_modimID:maskstruc.bsplinefit[4].modimID, bsplinefit_01_name:maskstruc.bsplinefit[4].name, bsplinefit_01_modim:maskstruc.bsplinefit[4].modim, bsplinefit_01_chi2:maskstruc.bsplinefit[4].chi2, bsplinefit_01_nthetaID:maskstruc.bsplinefit[4].nthetaID, bsplinefit_01_rbkptID:maskstruc.bsplinefit[4].rbkptID, $
    bsplinefit_012_modimID:maskstruc.bsplinefit[5].modimID, bsplinefit_012_name:maskstruc.bsplinefit[5].name, bsplinefit_012_modim:maskstruc.bsplinefit[5].modim, bsplinefit_012_chi2:maskstruc.bsplinefit[5].chi2, bsplinefit_012_nthetaID:maskstruc.bsplinefit[5].nthetaID, bsplinefit_012_rbkptID:maskstruc.bsplinefit[5].rbkptID, $
    bsplinefit_014_modimID:maskstruc.bsplinefit[6].modimID, bsplinefit_014_name:maskstruc.bsplinefit[6].name, bsplinefit_014_modim:maskstruc.bsplinefit[6].modim, bsplinefit_014_chi2:maskstruc.bsplinefit[6].chi2, bsplinefit_014_nthetaID:maskstruc.bsplinefit[6].nthetaID, bsplinefit_014_rbkptID:maskstruc.bsplinefit[6].rbkptID, $
    bsplinefit_0124_modimID:maskstruc.bsplinefit[7].modimID, bsplinefit_0124_name:maskstruc.bsplinefit[7].name, bsplinefit_0124_modim:maskstruc.bsplinefit[7].modim, bsplinefit_0124_chi2:maskstruc.bsplinefit[7].chi2, bsplinefit_0124_nthetaID:maskstruc.bsplinefit[7].nthetaID, bsplinefit_0124_rbkptID:maskstruc.bsplinefit[7].rbkptID}

  splog, "CREATE: "+maskfile 
  mwrfits, flat_maskstruc, maskfile, mask_header, /create, /silent

 end
;---------------------------------------------------------------------
pro acsproc_database_update_photmod_test, program_number, visit_number, filename, message=message, confirm=confirm
  ;==============================================================================
  ; Database Function: insert exposure redux_name and receive message response from database
  ;==============================================================================

  cmd = 'photmod_test_update'
  arg = " -p "+strtrim(program_number,2) + " -v "+visit_number + " -f "+filename
  spawn, cmd+arg, response
  pos = strlen(cmd+" > ") & message = ""
  start = cmd+" > START*" & started = 0
  done = cmd+" > DONE*" & confirm = 0
  for i=0,n_elements(response)-1 do begin
    if not keyword_set(started) and strmatch(response[i],start) then started=1 $
    else if not keyword_set(confirm) and strmatch(response[i],done) then confirm=1 $
    else if keyword_set(started) and not keyword_set(confirm) then message += strmid(response[i],pos)
  endfor
end
;---------------------------------------------------------------------

;================================================================================
pro plot_dmr, MYFUNCT, p, iter, fnorm, functargs=functargs, I_image=I_image, $
              parinfo=parinfo, quiet=quiet, DOF=dof, $
    	      PFORMAT=pformat, UNIT=unit

  common tmp, I_model
  loadct, 25, /silent
  print, '# of Iteration: ', iter
  print, 'chi2= ', fnorm
  print, p
  window, xsize=960, ysize=320
  tvscl_image, I_image, position=[0.06, 0.1, 0.34, 0.94], max_image=0.5*max(I_model), /erase
  tvscl_image, I_model, position=[0.36, 0.1, 0.64, 0.94], max_image=0.5*max(I_model)
  xyouts, 0.5, 0.95, '# of Iteration: '+string(iter, format='(I4)'), alignment=0.5, /normal, charsize=2.0
  tvscl_image, (I_image-I_model), position=[0.66, 0.1, 0.94, 0.94]
end
;---------------------------------------------------------------------

function SB_approx, x
  common par, pars

  I_0=pars[0] 
  r_e=pars[1]
  nser=pars[6]
  r_b=pars[7]
  gama=pars[9]
  b_ser=(1.9992*nser-0.3271) ; empirical relation for b_ser 
  if (b_ser le 0.0) then b_ser=1.0 
  ;b_ser=1.0D 

  if ((size(x))[0] gt 1) then splog, 'Input X should be a scalar or 1D array.'
  result=0.0*x
  nlen=(size(x))[1]
  wh_PL=where(x le r_b, c_PL)
  wh_Sersic=where(x gt r_b, c_Sersic)
  if (c_PL gt 0) then result[wh_PL]=I_0*exp(-b_ser*(r_b/r_e)^(1.0D/nser))*(r_b/x[wh_PL])^gama
  if (c_Sersic gt 0) then result[wh_Sersic]=I_0*exp(-b_ser*(x[wh_Sersic]/r_e)^(1.0D/nser))
   
  return, result
end

function SB, x
  common par, pars

  I_0=pars[0] 
  r_e=pars[1]
  nser=pars[6]
  r_b=pars[7]
  alpha=pars[8]
  gama=pars[9]
  b_ser=(1.9992*nser-0.3271) ; empirical relation for b_ser
  if (b_ser le 0.0) then b_ser=1.0 
  ;b_ser=1.0D
  
  if (I_0 le 0.0) then begin &$
     splog, 'Negative amplitude' &$
     I_0=10.0^(-12) &$
  endif
  
  if (alpha lt 50) then begin &$
     yy=alog(I_0)+(gama/alpha)*alog(1.0+(r_b/x)^alpha)-b_ser*((x^alpha+r_b^alpha)/r_e^alpha)^(1.0D/(alpha*nser)) &$
;     yy[where(yy lt -88.0)]=-88.0 &$
;     yy[where(yy gt 88.0)]=88.0 &$
     y=exp(yy) &$
  endif else begin &$
     y=SB_approx(x)
  endelse
  return, y
end

function integrand, x
  return, 2.0*!pi*x*SB(x) 
end

function Tlight, x, status=status
  status=0
  catch, err_status
  if err_status ne 0 then begin &$
     splog, 'Error message: ', !ERR_STRING &$
     status=-1 &$
     return, -1 &$
  endif else begin &$
     return, qsimp('integrand', 10.0^(-6), x, eps=10.^(-2))
  endelse
end

function Half_light, x
  qsimp, 'integrand', 10.0^(-6), x, integral, eps=10.^(-2)
  return, integral-0.5*Tlight(1000.0)
end


;================================================================================


pro acsproc_photmod_rhalf, visit, flatweight=flatweight, cycle=cycle,program=program, doplot=doplot

;common acsproc_photmod_rhalf_data

common acsproc_photmod_rhalf_data, acs_filter, dataroot
common par, pars1
  !Except = 0

  acs_filter = ['','F435W','F555W','F814W']
  acs_filtwave = [0,4350,5550,8140] ;Angstroms

  ; pixel size:
  dpix = 0.05D
  ;modim half-width
  hw = 140
  ; PSF half-width:
  phw = 18
  ;size of cutout
  imsize = 1024L
 
  if not keyword_set(cycle) then cycle='18'
  if not keyword_set(program) then program='12210'
  hst_dataroot=getenv('HST_DATAROOT')
  cycle_program = strtrim(cycle,2) + '/' + strtrim(program,2)
  dataroot = djs_filepath('', root=hst_dataroot, subdir=cycle_program)
  metafile = djs_filepath('.metafile.fits', root=dataroot)
  metastruc = mrdfits(metafile, 1)

  nrows = n_elements(metastruc)
  maskflag = 4
  selected = where(metastruc.visit eq visit and metastruc.subexposure eq 1,nvis) 
      
  if (nvis gt 0) then begin
  
    metastruc_visit = metastruc[selected]
    for q = 0,nvis-1 do begin

      combfile = acsproc_datafile(metastruc_visit[q],/combine)
      photmodfile = acsproc_datafile(metastruc_visit[q],/photmod,/test)
      ;photmodfile = acsproc_datafile(metastruc_visit[q],/photmod,old=1-keyword_set(flatweight))
      epsfile = acsproc_datafile(metastruc_visit[q],/photmod,/eps)

      bizfile = acsproc_datafile(metastruc_visit[q],/biz)
          
      maskfile = acsproc_datafile(metastruc_visit[q],/mask)
      maskstruc = acsproc_get_maskstruc(maskfile, header=mask_header)
      pixel_mask = (1B-maskstruc.fmask)*(1B-maskstruc.jmask)

      xc = round(sxpar(mask_header, 'CRPIX1'))-1L
      yc = round(sxpar(mask_header, 'CRPIX2'))-1L
   
      bspline_img = maskstruc.bsplinecen.modim

      x_bspline = dpix * (dindgen(2*hw+1) # replicate(1.d0, 2*hw+1) - hw)
      y_bspline = dpix * (replicate(1.d0, 2*hw+1) # dindgen(2*hw+1) -  hw)

      x = dpix * (findgen(imsize) - xc) # replicate(1.0, imsize)
      y = dpix * replicate(1.0, imsize) # (findgen(imsize) - yc)

      psf = mrdfits(combfile, 4)
      nx_psf = (size(psf))[1]
      junk = max(psf, wmax)
      pxc = wmax mod nx_psf
      pyc = wmax / nx_psf
      psf = psf[pxc-phw:pxc+phw,pyc-phw:pyc+phw]
      psf = psf / total(psf)

      img = mrdfits(combfile, 1, hdr)
;      img=0.0 
        ;photmodstruc = {fitfunc: fitfunc[i], img:photimg, pars: pars, perror: perror,  covar: covar, chisq: bestnorm, ndof: ndof, xc: xc, yc: yc}
        photmodstruc = mrdfits(photmodfile, 1, /silent)
        splog, "READING "+photmodfile
        photimg=photmodstruc.img
        pars = photmodstruc.pars
	perror=photmodstruc.perror

        ; Determine the half-light radius r_half using Bisection method
        pars1=pars
        r_cutoff=[1000.0D, 500.0D, 400.0D, 300.0D, 200.0D, 100.0D, 50.0D, 40.0D, 30.0D, 20.0D, 10.0D]
        ind=0
        total_light=Tlight(r_cutoff[0], status=status)
        ; Decide the cut-off radius iteratively until the total light converges
        while (status eq -1 and ind le 9) do begin &$
           ind=ind+1 &$
           total_light=Tlight(r_cutoff[ind], status=status) &$
        endwhile
           
        eps=10.0^(-4)
        n_step=1
        splog, "total_light="+strtrim(total_light,2)
        if (total_light gt 0.0) then begin &$
    	   x_l=0.1D & x_r=r_cutoff[ind] & x_m=0.5*(x_l+x_r) &$
           if (Tlight(x_l)-0.5*total_light gt 0.0) then splog, "x_l is incorrect!" &$
           if (Tlight(x_r)-0.5*total_light lt 0.0) then splog, "x_r is incorrect!" &$
  
           rtbis=x_l &$
           dx=x_r-x_l &$
           fmid=Tlight(x_r)-0.5*total_light lt 0.0 &$
           for i_step=0L, 39 do begin &$
              dx=dx*0.5 &$
              x_m=rtbis+dx &$
              fmid=Tlight(x_m)-0.5*total_light &$
              if (fmid lt 0.0) then rtbis=x_m &$
              if (abs(dx) gt eps and fmid ne 0.) then continue &$
              ;print, rtbis, fmid &$
              ;print, abs(Tlight(x_m)-0.5*total_light)/total_light &$
              ;print, x_r-x_l &$
           endfor &$
           r_half=rtbis &$
        endif else r_half=-1.0

        splog, 'r_half= ', r_half        

        ; Generate plots
        Radius_bspline=10.^(findgen(201)*0.01-1.0)
        L_R_bspline=0.0*Radius_bspline
        for i=0L, n_elements(L_R_bspline)-1 do begin &$
           wh_bspline=where(sqrt(x_bspline^2.0+y_bspline^2.0) le Radius_bspline[i]) &$
           L_R_bspline[i]=total(bspline_img[wh_bspline])*dpix^2. &$
        endfor

        Radius2=10.^(findgen(401)*0.01-1.0)
        L_R2=0.0*Radius2
        for i=0L, n_elements(L_R2)-1 do begin &$
           wh_2=where(sqrt(x^2.0+y^2.0) le Radius2[i]) &$
           L_R2[i]=total(photimg[wh_2])*dpix^2. &$
        endfor

        R_min=1.0 & R_max=r_cutoff[ind]
        Radius=10.^(findgen(fix((alog10(R_max)+1)/0.01)+1)*0.01-1.0)
        L_R=Tlight(Radius)

       epsfile=metastruc_visit.targname
       set_plot, 'ps'
       !p.font=0
       !p.thick=4
       !x.thick=4
       !y.thick=4
        device, /encapsulated, /color, xsize=20, ysize=15, filename=epsfile+'_light.eps', bits=8, /times
        loadct, 0, /silent
        plot, radius, L_R, xrange=[0.0, 400.0], xtitle='Radius [arcsec]', ytitle='L(R)', $
              charthick=4, symsize=1.5
        oplot, [r_half, r_half], !y.crange, thick=4, linestyle=2
        plot, Radius, L_R, xrange=[0.0, 10.0], xtitle='Radius [arcsec]', ytitle='L(R)', $
              charthick=4, symsize=1.5, position=[0.4, 0.2, 0.9, 0.7], /noerase
        oplot, Radius_bspline, L_R_bspline, color=128, thick=4
        legend, ['core Sersic', 'bspline'], colors=[0, 128], linestyle=[0, 0]
        device, /close

        loadct, 25, /silent
        device, /color, /encapsulated, filename=epsfile+'_model.eps', xsize=27, ysize=9, bits=8, /times
        tvscl_image, img[imsize/2-hw:imsize/2+hw, imsize/2-hw:imsize/2+hw], position=[0.06, 0.1, 0.34, 0.94], min_image=10.^(-8), max_image=0.5*max(photimg), /erase
        tvscl_image, photimg[imsize/2-hw:imsize/2+hw, imsize/2-hw:imsize/2+hw], position=[0.36, 0.1, 0.64, 0.94], min_image=10.^(-8), max_image=0.5*max(photimg)
        tvscl_image, (img-photimg)[imsize/2-hw:imsize/2+hw, imsize/2-hw:imsize/2+hw], position=[0.66, 0.1, 0.94, 0.94]
        device, /close
        set_plot, 'x'
        splog, 'eps files '+epsfile+'_model.eps'+' created'

       
      if (NOT keyword_set(doplot)) then begin 
        ; Calculate the uncertainty in half-light radius r_half
        r_half_arr=fltarr(100)
        flag_arr=intarr(100)+1
        i=0L
        count=0L
        while (i lt 100 and count le 10000) do begin &$
           count=count+1 &$
           ;epsilon=randomu(seed, 10) &$
           ;pars_tmp=pars+(epsilon-0.5)*2.0*perror &$
           epsilon=randomn(seed, 10) &$
           pars_tmp=pars+(epsilon*perror) &$
           if pars_tmp[0] lt 0.0 then continue &$
           if pars_tmp[1] lt 0.0 then continue &$ 
           if pars_tmp[6] lt 0.0 then continue &$ 
           if pars_tmp[7] lt 0.0 then continue &$ 
           if pars_tmp[8] lt 0.0 then continue &$  
           if pars_tmp[9] lt 0.0 then continue &$ 
           
           pars1=pars_tmp &$
            
           n_step=1 &$
           if (total_light gt 0.0) then begin &$
    	      x_l=0.1D & x_r=r_cutoff[ind] & x_m=0.5*(x_l+x_r) &$
              if (Tlight(x_l)-0.5*total_light gt 0.0) then begin &$ 
                 splog, "x_l is incorrect!" &$
                 continue &$
              endif
             
              if (Tlight(x_r)-0.5*total_light lt 0.0) then begin &$
                 splog, "x_r is incorrect!" &$
                 continue &$
              endif
  
              rtbis=x_l &$
              dx=x_r-x_l &$
              fmid=Tlight(x_r)-0.5*total_light lt 0.0 &$
              for i_step=0L, 39 do begin &$
                 dx=dx*0.5 &$
                 x_m=rtbis+dx &$
                 fmid=Tlight(x_m)-0.5*total_light &$
                 if (fmid lt 0.0) then rtbis=x_m &$
                 if (abs(dx) gt eps and fmid ne 0.) then continue &$
                 ;print, rtbis, fmid &$
                 ;print, abs(Tlight(x_m)-0.5*total_light)/total_light &$
                 ;print, x_r-x_l &$
              endfor &$
              r_half_tmp=rtbis &$
           endif else continue &$

           splog, 'Iteration ', i+1, '  ', 'r_half= ', r_half_tmp &$
           r_half_arr[i]=r_half_tmp &$
           i=i+1 &$
        endwhile
        r_halferr=-1.0
        if (count le 10000) then r_halferr=stddev(r_half_arr[where(r_half_arr gt 0.1 and r_half_arr lt r_cutoff[ind])])
        splog, 'r_halferr= ', r_halferr

         ; Get the deVaucouleurs fit parameters and the photometric
         ; header keyword parameters:
         photflam = sxpar(hdr, 'PHOTFLAM')
         photplam = sxpar(hdr, 'PHOTPLAM')

         ; Convert core sersic model amplitude to counts, and to AB magnitude:
         cps_cs= total_light/(dpix^2)
         abmag_zpt=-2.5*alog10(photflam)-21.10-5.*alog10(photplam)+18.6921
         abmag_cs=-2.5*alog10(cps_cs)+abmag_zpt
    
         cps_cserr=cps_cs*(perror[0]/pars[0])

         ;endfor 

         photmodstruc2 = mrdfits(photmodfile, 2)

         photmodstruc2.acs_rhalfcs = r_half
         photmodstruc2.acs_rhalfcserr = r_halferr
         photmodstruc2.abmag_cs=abmag_cs
         photmodstruc2.cps_cs=cps_cs
         photmodstruc2.cps_cserr=cps_cserr

         
         mwrfits, photmodstruc, photmodfile, /silent, /create
         mwrfits, photmodstruc2, photmodfile, /silent
               
         ;has_photmod = maskstruc.has_photmod
         ;maskstruc.has_photmod = file_test(photmodfile)
         ;if has_photmod ne maskstruc.has_photmod then acsproc_maskstruc_write, maskstruc, maskfile, header=mask_header
         
         acsproc_database_update_photmod_test, metastruc_visit[q].proposid, metastruc_visit[q].visit, file_basename(photmodfile), confirm=confirm, message=message
         if confirm then splog, 'Database > '+message else splog, 'update photmod_test failed.'

      endif

    endfor
  
  endif
        
end
