pro acsproc_07_bizzle, process, statuscode, count, selected=selected, pmasking=pmasking

  common acsproc_data
  common acsproc_metastruc
  common acsproc_maskstruc

  ; Get the info on each visit (reference exposure) to be processed (correct statuscode) checking for a comb & mask file:
  
  nrows = n_elements(metastruc)
  void = where(metastruc.subexposure eq 1,nvisits)
  maskflag = 3
  
  if not keyword_set(selected) then selected = where(metastruc.statuscode eq statuscode and metastruc.combflag eq 1 and metastruc.maskflag eq maskflag-1 and (metastruc.fmaskflag eq 1 or metastruc.jmaskflag eq 1) and metastruc.subexposure eq 1,nvis) else nvis = n_elements(selected)
  
  if n_elements(metastruc) gt 0 then steps=(metastruc[0].nexp+1)*nvis else steps=nvis
  
  process->Start, steps=steps
  process->Setmessage, 'Processing Proposal ID '+program+' ('+strtrim(nvis,2)+'/'+strtrim(nvisits,2)+' Visits)', /label
  process->Setmessage, ''  
  metastruc_changed = 0L
  
  if (nvis gt 0) then begin
  
    metastruc_visit = metastruc[selected]
    
    for q = 0,nvis-1 do begin
    
      if process->CheckCancel() then begin
        process->Step, steps
        process->Setmessage, 'Cancelling...', /append
        if metastruc_changed then acsproc_metastruc_update
        process->Destroy
        return
      endif else process->Increment, 1.0
    
      combfile = acsproc_datafile(metastruc_visit[q],/combine)
      maskfile = acsproc_datafile(metastruc_visit[q],/mask)
      bizfile = acsproc_datafile(metastruc_visit[q],/biz)
      maskstruc = acsproc_get_maskstruc(maskfile, header=mask_header)
      xc = round(sxpar(mask_header, 'CRPIX1'))-1L
      yc = round(sxpar(mask_header, 'CRPIX2'))-1L

      fmask = 1B-maskstruc.fmask[xc-hw:xc+hw,yc-hw:yc+hw]
      jmask = 1B-maskstruc.jmask[xc-hw:xc+hw,yc-hw:yc+hw]
  
      dmodim = maskstruc.bsplinecen.dmodim
      selected_modim = acsproc_modim(/fit_nthetaID)
      modim = selected_modim.modim

      hdr0 = headfits(combfile, exten=0)
      img = (mrdfits(combfile, 1))[xc-hw:xc+hw,yc-hw:yc+hw]
      err = (mrdfits(combfile, 2))[xc-hw:xc+hw,yc-hw:yc+hw]
      nim = (mrdfits(combfile, 3))[xc-hw:xc+hw,yc-hw:yc+hw]
      ivar = (nim gt 0) / (err^2 + (nim le 0))
  
      psf = mrdfits(combfile, 4)
      nx_psf = (size(psf))[1]
      junk = max(psf, wmax)
      pxc = wmax mod nx_psf
      pyc = wmax / nx_psf
      bpsf = psf[pxc-phw:pxc+phw,pyc-phw:pyc+phw]
      bpsf = bpsf / total(bpsf)
      
      if keyword_set(pmasking) then begin
         np_side = 2*hw+1
         pmask = replicate(0B, np_side, np_side)
         pmask[phw:np_side-(phw+1),phw:np_side-(phw+1)] = 1B
         ivar = ivar * pmask
      endif

      ; x and y pixel ofsets of centroid from central pixel:
      xdiff = -maskstruc.bsplinecen.pars[0]/dpix
      ydiff = maskstruc.bsplinecen.pars[1]/dpix
      xcen = float(hw) + xdiff
      ycen = float(hw) + ydiff
      
      hdrx = ['']
      sxaddpar, hdrx, 'XCENTER', xcen, '(zero-based centered pixels)'
      sxaddpar, hdrx, 'YCENTER', ycen, '(zero-based centered pixels)'
      sxaddpar, hdrx, 'AXISRAT', maskstruc.bsplinecen.pars[2], '(b-spline ellipsoid)'
      sxaddpar, hdrx, 'AXISANGL', maskstruc.bsplinecen.pars[3], '(E from N)'
            
      file_delete, bizfile, /allow_nonexistent, /quiet
      
      mwrfits, img, bizfile, hdr0, /create, /silent
      mwrfits, ivar, bizfile, /silent
      mwrfits, nim, bizfile, /silent
      mwrfits, psf, bizfile, /silent
      mwrfits, fmask, bizfile, /silent
      mwrfits, jmask, bizfile, /silent
      mwrfits, fmask, bizfile, /silent
      mwrfits, modim, bizfile, hdrx, /silent
      mwrfits, dmodim, bizfile, hdrx, /silent
      mwrfits, bpsf, bizfile, /silent
      
      splog, "CREATE: "+bizfile
      
      has_bizzle = maskstruc.has_bizzle
      maskstruc.has_bizzle = file_test(bizfile)
      if has_bizzle ne maskstruc.has_bizzle then acsproc_maskstruc_write, maskfile

      selected_modim_name = selected_modim.name
      if selected_modim.modimID eq 2 then begin
        ntheta_str = ''
        void = acsproc_ntheta(selected_modim,str=ntheta_str)
        selected_modim_name += '['+ntheta_str+']'
      endif
       
      message = 'Selected MODIM is '+selected_modim_name
      splog, message
      acsproc_database_update_table, 'bspline', metastruc_visit[q].proposid, metastruc_visit[q].visit, file_basename(bizfile), confirm=confirm, message=message
      if confirm then splog, 'Database > '+message else splog, 'update failed.'

      process->setMessage, "Target "+metastruc_visit[q].uniqname+" Visit "+metastruc_visit[q].visit+' -> '+message, /APPEND
             
      visit = where(metastruc.visit eq metastruc_visit[q].visit,nexp)
      
      for s=0,nexp-1 do begin          
        process->Increment, 1.0
        r = visit[s]
        if (not metastruc_changed) then if (metastruc[r].maskflag ne maskflag) then metastruc_changed = 1L
        metastruc[r].maskflag = maskflag
      endfor
                          
    endfor
  
  endif
      
  if metastruc_changed then acsproc_metastruc_update
  wait, 2
  process->Destroy  
  
end
