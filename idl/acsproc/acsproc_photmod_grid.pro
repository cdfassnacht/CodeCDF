function acsproc_datafile, subexposure, cutout=cutout, rectify=rectify, combine=combine, mask=mask, photmod=photmod, biz=biz, old=old, test=test, grid=grid, par=par, hst_root=hst_root, ext=ext, subext=subext, eps=eps
  
  common acsproc_photmod_grid_data, acs_filter, dataroot, scratch_dir
  
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
  else if keyword_set(eps) and keyword_set(photmod) then extension = '_photmod.eps' $
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
    if keyword_set(grid) then datafile = scratch_dir + subexposure[0].uniqname + '_' + $
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
pro acsproc_database_update_photmod_test_grid, program_number, visit_number, filename, grid_size=grid_size, grid_index=grid_index, difficulty=difficulty, message=message, confirm=confirm
  ;==============================================================================
  ; Database Function: insert exposure redux_name and receive message response from database
  ;==============================================================================

  cmd = 'photmod_test_grid_update'
  arg = " -p "+strtrim(program_number,2) + " -v "+visit_number + " -f "+filename + " -gs "+strtrim(grid_size,2)+ " -gi "+strtrim(grid_index,2)+ " -d "+strtrim(difficulty,2)
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
  common par, pars1

  I_0=pars1[0] 
  r_e=pars1[1]
  nser=pars1[6]
  r_b=pars1[7]
  gama=pars1[9]
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
  common par, pars1

  I_0=pars1[0] 
  r_e=pars1[1]
  nser=pars1[6]
  r_b=pars1[7]
  alpha=pars1[8]
  gama=pars1[9]
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
  return, x*SB(x) 
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


pro acsproc_photmod_grid, visit, flatweight=flatweight, cycle=cycle,program=program, grid_index=grid_index, grid_size=grid_size, skip_pars=skip_pars

common acsproc_photmod_grid_data, acs_filter, dataroot, scratch_dir
common par, pars1
  !Except = 0


  if not keyword_set(grid_size) then grid_size = 12
  if (grid_index gt grid_size) then message, 'grid_index must be less/equal to grid_size'
  if not keyword_set(cycle) or not keyword_set(program) then message, 'cycle/program keywords required'

  ;modim half-width
  hw = 140
  ;size of cutout
  imsize = 1024L
  
  ; PSF half-width:
  phw = 18
  ; pixel size:
  dpix = 0.05D

  acs_filter = ['','F435W','F555W','F814W']
  acs_filtwave = [0,4350,5550,8140] ;Angstroms
  fitfunc = ['asb_devauc_multi','asb_sersic_multi', 'core_sersic']
  fitfunc_label = ['de Vauc','Sersic', 'core_sersic']

  hst_dataroot=getenv('HST_DATAROOT')
  hst_scratch_dir=getenv('HST_SCRATCH_DIR')
  cycle_program = strtrim(cycle,2) + '/' + strtrim(program,2)
  dataroot = djs_filepath('', root=hst_dataroot, subdir=cycle_program)
  scratch_dir = djs_filepath('', root=hst_scratch_dir, subdir='photmod_test_grid-'+string(grid_size,format='(i2.2)')+'/grid-'+string(grid_index,format='(i2.2)'))
  if not file_test(scratch_dir) then begin
    splog, "CREATE: "+scratch_dir
    file_mkdir, scratch_dir
  endif
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
      photmod_gridfile = acsproc_datafile(metastruc_visit[q],/photmod,/test,/grid)
      bizfile = acsproc_datafile(metastruc_visit[q],/biz)

      maskfile = acsproc_datafile(metastruc_visit[q],/mask)
      maskstruc = acsproc_get_maskstruc(maskfile, header=mask_header)
  
      xc = round(sxpar(mask_header, 'CRPIX1'))-1L
      yc = round(sxpar(mask_header, 'CRPIX2'))-1L
      
      img = mrdfits(combfile, 1, hdr)
      err = mrdfits(combfile, 2)
      nim = mrdfits(combfile, 3)
      
      ivar = (nim gt 0) * (err gt 0.) / (err^2 + (nim le 0) + (err le 0.))
      
      pixel_mask = (1B-maskstruc.fmask)*(1B-maskstruc.jmask)
      ivar = ivar * pixel_mask
      
      psf = mrdfits(combfile, 4)
      nx_psf = (size(psf))[1]
      junk = max(psf, wmax)
      pxc = wmax mod nx_psf
      pyc = wmax / nx_psf
      psf = psf[pxc-phw:pxc+phw,pyc-phw:pyc+phw]
      psf = psf / total(psf)
  
      ymask = 0B * byte(ivar)
      ymask[xc-hw:xc+hw,yc-hw:yc+hw] = 1B
      pmask = replicate(0B, imsize, imsize)
      pmask[phw:imsize-(phw+1),phw:imsize-(phw+1)] = 1B
      ivar = ivar * pmask
      pmask = 0

      ; Get the right pixel statistics:
      djs_iterstat, img * sqrt(ivar), sigma=errsig
      ;errsig=median(img*sqrt(ivar)*(1.0-ymask))/2.0
      ivar = ivar / errsig^2

      ; Make RA and Dec images:
      aimg = -dpix * (findgen(imsize) - xc) # replicate(1.0, imsize)
      dimg = dpix * replicate(1.0, imsize) # (findgen(imsize) - yc)

      photmodstruc = mrdfits(photmodfile, 1, /silent)
      photmodstruc2 = mrdfits(photmodfile, 2, /silent)

      splog, "READING "+photmodfile
      
      pars = photmodstruc.pars
      perror=photmodstruc.perror
      ivar = ivar * (photmodstruc2.outermask eq 0)

      pars_tmp=pars
      if (perror[1] lt pars[1]) then begin
         pars_tmp[1]=pars[1]-perror[1]+(grid_index-1)*(2*perror[1]/(grid_size-1.0))
         difficulty = 1
      endif else if (perror[1] gt pars[1] and perror[1] lt 10.*pars[1]) then begin
         pars_tmp[1]=0.2+(grid_index-1)*(pars[1]+perror[1])/(grid_size-1.0)
         difficulty = 2
      endif else begin
         pars_tmp[1]=0.2+(grid_index-1)*(10.0*pars[1])/(grid_size-1.0)         
         difficulty = 3
      endelse
      
      if not keyword_set(skip_pars) then begin

          parinfo=replicate({value: 0.D, fixed: 0, limited:[0, 0], limits: [0.0, 0.0]}, n_elements(pars_tmp))
          parinfo.value=pars_tmp
          parinfo[1].fixed=1
          if difficulty gt 1 then begin
              parinfo[0].limited[0]=1
              parinfo[0].limits[0]=10.0^(-12)
              parinfo[2].fixed=1
              parinfo[3].fixed=1
              parinfo[4].fixed=1
              parinfo[5].fixed=1
          endif

          i=2
          splog, "MPFIT > SMARTMASK: "+fitfunc_label[i]+" **************************************************************"
          if keyword_set(plot) then $
            pars1 = mpfit2dfun(fitfunc[i], dimg, aimg, $
            img, 0.*img, pars_tmp, weights=ivar, functargs={psf: psf, rmin: 0.1*dpix}, parinfo=parinfo, $
            iterproc='plot_dmr', iterargs={I_image:img},$ ;visualize the fit
            bestnorm=bestnorm, perror=perror, niter=niter, yfit=yfit, covar=covar) $
          else $
            pars1 = mpfit2dfun(fitfunc[i], dimg, aimg, $
            img, 0.*img, pars_tmp, weights=ivar, functargs={psf: psf, rmin: 0.1*dpix}, parinfo=parinfo, $
            bestnorm=bestnorm, perror=perror, niter=niter, yfit=yfit, covar=covar)
            
          pars1[1] = abs(pars1[1])        
          photmodstruc.pars = pars1
          photmodstruc.perror = perror
          photmodstruc.covar = covar
          photmodstruc.chisq = bestnorm

          photmodstruc2.acs_rcs=pars1[1]
          photmodstruc2.acs_rcserr=perror[1]
          photmodstruc2.acs_qcs=pars1[4]
          photmodstruc2.acs_qcserr=perror[4]
          photmodstruc2.acs_pacs=pars1[5]
          photmodstruc2.acs_pacserr=perror[5]
          photmodstruc2.acs_ncs=pars1[6]
          photmodstruc2.acs_ncserr=perror[6]
          photmodstruc2.acs_rcorecs=pars1[7]
          photmodstruc2.acs_rcorecserr=perror[7]
          photmodstruc2.acs_alphacs=pars1[8]
          photmodstruc2.acs_alphacserr=perror[8]
          photmodstruc2.acs_gammacs=pars1[9]
          photmodstruc2.acs_gammacserr=perror[9]
          photmodstruc2.acs_chi2cs=bestnorm
          
          if (photmodstruc2.acs_qcs gt 1.) then begin
                photmodstruc2.acs_qcs=1./photmodstruc2.acs_qcs
                photmodstruc2.acs_qcserr=photmodstruc2.acs_qcserr*photmodstruc2.acs_qcs^2
                photmodstruc2.acs_pacs+=90.
          endif
          pa=photmodstruc2.acs_pacs mod 180.
          ptest=pa ge 0.
          pa=pa*ptest+(180.-abs(pa))*(1B-ptest)
          photmodstruc2.acs_pacs=pa


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
          
          ; Convert core sersic model amplitude to counts, and to AB magnitude:
          photmodstruc2.cps_cs= 2.0*!pi*total_light/(dpix^2)
          photmodstruc2.abmag_zpt=-2.5*alog10(photmodstruc2.photflam)-21.10-5.*alog10(photmodstruc2.photplam)+18.6921
          photmodstruc2.abmag_cs=-2.5*alog10(photmodstruc2.cps_cs)+photmodstruc2.abmag_zpt
          photmodstruc2.cps_cserr=photmodstruc2.cps_cs*(perror[0]/pars1[0])
          photmodstruc2.acs_rhalfcs = r_half
          photmodstruc2.acs_rhalfcserr = -1.0


          splog, "CREATE: "+photmod_gridfile
          mwrfits, photmodstruc, photmod_gridfile, /silent, /create
          mwrfits, photmodstruc2, photmod_gridfile, /silent
          
      endif
            
      acsproc_database_update_photmod_test_grid, metastruc_visit[q].proposid, metastruc_visit[q].visit, file_basename(photmodfile), grid_size=grid_size, grid_index=grid_index, difficulty=difficulty, confirm=confirm, message=message
      if confirm then splog, 'Database > '+message else splog, 'update photmod_test failed.'
      

    endfor
  
  endif
        
end
