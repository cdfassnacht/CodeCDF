pro acsproc_06_bspline, process, statuscode, count, selected=selected, pmasking=pmasking

  common acsproc_data
  common acsproc_metastruc
  common acsproc_maskstruc

  ; Get the info on each visit (reference exposure) to be processed (correct statuscode) checking for a comb & mask file:
  
  nrows = n_elements(metastruc)
  void = where(metastruc.subexposure eq 1,nvisits)
  maskflag = 2
  if not keyword_set(selected) then selected =  where(metastruc.statuscode eq statuscode and metastruc.combflag eq 1 and metastruc.maskflag eq maskflag-1 and (metastruc.fmaskflag eq 1 or metastruc.jmaskflag eq 1) and metastruc.subexposure eq 1,nvis) else nvis = n_elements(selected)

  if n_elements(metastruc) gt 0 then steps=(metastruc[0].nexp+1+8)*nvis else steps=nvis
  
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
      maskstruc = acsproc_get_maskstruc(maskfile, header=mask_header)
      void = where(maskstruc.fmask ne 0B,nfmask)      
  
      xc = round(sxpar(mask_header, 'CRPIX1'))-1L
      yc = round(sxpar(mask_header, 'CRPIX2'))-1L

      modim = acsproc_modim('bsplinecen')
      ntheta = acsproc_ntheta(modim)
      rbkpt = acsproc_rbkpt(modim)
      
      img = (mrdfits(combfile, 1))[xc-hw:xc+hw,yc-hw:yc+hw]
      err = (mrdfits(combfile, 2))[xc-hw:xc+hw,yc-hw:yc+hw]
      nim = (mrdfits(combfile, 3))[xc-hw:xc+hw,yc-hw:yc+hw]
      ivar = (nim gt 0) / (err^2 + (nim le 0))
      void = where(maskstruc.fmask ne 0B,nfmask)
      void = where(maskstruc.jmask ne 0B,njmask)
      pixel_mask = (1B-maskstruc.fmask[xc-hw:xc+hw,yc-hw:yc+hw])*(1B-maskstruc.jmask[xc-hw:xc+hw,yc-hw:yc+hw])
      void = where(pixel_mask ne 0L,npixels)
      invvar = ivar * pixel_mask

      aimg = -dpix * (dindgen(2*hw+1) # replicate(1.d0, 2*hw+1) - hw)
      dimg = dpix * (replicate(1.d0, 2*hw+1) # dindgen(2*hw+1) -  hw)
  
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

      ft = {x: dimg, y: aimg, psf: bpsf, rbkpt: rbkpt, ntheta: ntheta, $
            data: img, invvar: invvar, deviates: 1B}

      ; The centering fit:
      bspar = [0.01D, 0.01D, .9D, 45.D]
      ; bspar = [0.001, 0.001, ms[j].mofpar_1[5], ms[j].mofpar_1[6]]

      message = 'BSPLINE: '+maskfile
      process->setMessage, "Target "+metastruc_visit[q].uniqname+" Visit "+metastruc_visit[q].visit+' -> '+message, /APPEND
      
      bspars = mpfit('bspline_ellip', bspar, functargs=ft, perror=perror)
      dof = (npixels-n_elements(bspars)-n_elements(ntheta)*n_elements(rbkpt))
      
      maskstruc.bsplinecen.name = 'Bspline Cen'
      maskstruc.bsplinecen.pars = bspars
      maskstruc.bsplinecen.pars_err = perror
      maskstruc.bsplinecen.modim = bspline_ellip(bspars, x=ft.x, y=ft.y, data=ft.data, ntheta=ft.ntheta, invvar=ft.invvar, psf=ft.psf, rbkpt=ft.rbkpt, sset=sset_bspsf)
      maskstruc.bsplinecen.dmodim = bspline_ellip_valu(ft.x, ft.y, bspars, sset_bspsf)
      maskstruc.bsplinecen.chi2 = total((((img - maskstruc.bsplinecen.modim) * sqrt(invvar))[*])^2) / dof
      minchi2 = maskstruc.bsplinecen.chi2
      maskstruc.nthetaID =  0
      ntheta_str = '0'
      maskstruc.bsplinecen.nthetaID = modim.nthetaID
      maskstruc.bsplinecen.rbkptID = modim.rbkptID
      maskstruc.nthetaID =  1
      
      best_ntheta = ntheta_str
      modim = acsproc_modim('bsplinefit')
      for nthetaID = 0,n_elements(modim)-1 do begin  
        process->Increment, 1.0      
        modim[nthetaID].nthetaID = nthetaID
        ntheta = acsproc_ntheta(modim[nthetaID], str=ntheta_str)
        rbkpt = acsproc_rbkpt(modim[nthetaID])
        dof = (npixels-n_elements(bspars)-n_elements(ntheta)*n_elements(rbkpt))

        maskstruc.bsplinefit[nthetaID].name = 'Bspline Fit ['+ntheta_str+']'
        maskstruc.bsplinefit[nthetaID].modim = bspline_ellip(bspars, x=ft.x, y=ft.y, data=ft.data, ntheta=ntheta, invvar=ft.invvar, psf=ft.psf, rbkpt=rbkpt, sset=sset)
        maskstruc.bsplinefit[nthetaID].chi2 = total((((img - maskstruc.bsplinefit[nthetaID].modim) * sqrt(invvar))[*])^2) / dof
        maskstruc.bsplinefit[nthetaID].nthetaID = modim[nthetaID].nthetaID
        maskstruc.bsplinefit[nthetaID].rbkptID = modim[nthetaID].rbkptID
        if maskstruc.bsplinefit[nthetaID].chi2 lt minchi2 then begin
          minchi2 = maskstruc.bsplinefit[nthetaID].chi2
          best_ntheta = ntheta_str
        endif
        
        message = 'BSPLINEFIT: [ntheta='+ntheta_str+'] has  reduced chi2 ='+strtrim(maskstruc.bsplinefit[nthetaID].chi2,2)
        process->setMessage, "Target "+metastruc_visit[q].uniqname+" Visit "+metastruc_visit[q].visit+' -> '+message, /APPEND
        
      endfor
      
      message = 'BSPLINE ['+best_ntheta+'] has lowest reduced chi2 ='+strtrim(minchi2,2)
      process->setMessage, "Target "+metastruc_visit[q].uniqname+" Visit "+metastruc_visit[q].visit+' -> '+message, /APPEND

      acsproc_maskstruc_write, maskfile
             
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
  process->Destroy  
  
end
