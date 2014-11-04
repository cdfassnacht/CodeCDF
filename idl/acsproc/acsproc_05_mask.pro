pro acsproc_05_mask, process, statuscode, count, selected=selected, old=old

  common acsproc_data
  common acsproc_metastruc
  common acsproc_maskstruc
  
  ;old = 1L
  
  ; Get the info on each visit (reference exposure) to be processed (correct statuscode) checking for a comb file:
  
  nrows = n_elements(metastruc)
  void = where(metastruc.subexposure eq 1,nvisits)
  if not keyword_set(selected) then selected = where(metastruc.statuscode eq statuscode and metastruc.combflag eq 1 and metastruc.subexposure eq 1,nvis)  else nvis = n_elements(selected)
  
  if n_elements(metastruc) gt 0 then steps=(metastruc[0].nexp+1)*nvis else steps=nvis
  
  process->Start, steps=steps
  process->Setmessage, 'Processing Proposal ID '+program+' ('+strtrim(nvis,2)+'/'+strtrim(nvisits,2)+' Visits)', /label
  process->Setmessage, ''  
  metastruc_changed = 0L
  
  if (nvis gt 0) then begin

    maskstruc = acsproc_load_maskstruc()
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
      
      if keyword_set(old) then begin
        oldmaskfile = acsproc_datafile(metastruc_visit[q],/mask,/old)
        oldmaskstruc = acsproc_load_maskstruc(fits=mrdfits(oldmaskfile, 1))
      endif
      
      
      mask_header = headfits(combfile, exten=1)
      xc = round(sxpar(mask_header, 'CRPIX1'))-1L
      yc = round(sxpar(mask_header, 'CRPIX2'))-1L
        
      if not keyword_set(old) then begin
        img = (mrdfits(combfile, 1))[xc-hw:xc+hw,yc-hw:yc+hw]
        err = (mrdfits(combfile, 2))[xc-hw:xc+hw,yc-hw:yc+hw]
        nim = (mrdfits(combfile, 3))[xc-hw:xc+hw,yc-hw:yc+hw]
        ivar = (nim gt 0) / (err^2 + (nim le 0))
        xmap = findgen(2*hw+1) # replicate(1., 2*hw+1) - hw
        ymap = transpose(xmap)
        
        ; Moffat profile parameters, in pixel basis:
        moffati = 10.              ; central intensity
        moffatr = 3.               ; scale radius
        moffatg = 1.                ; "gamma" exponent
        moffatq = 0.85              ; axis ratio
        moffata = 15.               ; position angle

        moffatpar = fltarr(7)
        moffatpar[0] = moffati
        moffatpar[1] = moffatr
        moffatpar[2] = moffatg
        ; moffatpar[3,*] = decstack[fhw,fhw,*]
        moffatpar[3] = 0.1
        ; moffatpar[4,*] = rastack[fhw,fhw,*]
        moffatpar[4] = 0.1
        moffatpar[5] = moffatq
        moffatpar[6] = moffata
        
        functargs={data: img, ivar: ivar, x: xmap, y: ymap, deviates: 1B}
        fpars = mpfit('mpmoffat', moffatpar, functargs=functargs, perror=perror)
        maskstruc.moffat.name = 'Moffat Fit'
        maskstruc.moffat.modim = mpmoffat(fpars, x=xmap, y=ymap)
        maskstruc.moffat.pars = fpars
        maskstruc.moffat.pars_err = perror      
        maskstruc.moffat.chi2 = total((((img - maskstruc.moffat.modim) * sqrt(ivar))[*])^2) / (n_elements(img)-n_elements(fpars))
        
      endif
      
      maskstruc.xc = xc
      maskstruc.yc = yc
      
      if keyword_set(old) then begin
        maskstruc.fmask = oldmaskstruc.fmask
        maskstruc.jmask = oldmaskstruc.jmask
        maskstruc.modimID = oldmaskstruc.modimID
        maskstruc.nthetaID = oldmaskstruc.nthetaID
        maskstruc.has_bizzle = oldmaskstruc.has_bizzle
        maskstruc.has_photmod = 0L; oldmaskstruc.has_photmod
        maskstruc.has_modim = oldmaskstruc.has_modim
        maskstruc.lens_status = oldmaskstruc.lens_status
        maskstruc.moffat =  oldmaskstruc.moffat
        maskstruc.bsplinecen =  oldmaskstruc.bsplinecen
        maskstruc.bsplinefit[0].modimID =  oldmaskstruc.bsplinefit[0].modimID
        maskstruc.bsplinefit[0].name =  oldmaskstruc.bsplinefit[0].name
        maskstruc.bsplinefit[0].modim =  oldmaskstruc.bsplinefit[0].modim
        maskstruc.bsplinefit[0].chi2 =  oldmaskstruc.bsplinefit[0].chi2
        maskstruc.bsplinefit[0].nthetaID =  oldmaskstruc.bsplinefit[0].nthetaID
        maskstruc.bsplinefit[0].rbkptID =  oldmaskstruc.bsplinefit[0].rbkptID
        maskstruc.bsplinefit[1].modimID =  oldmaskstruc.bsplinefit[1].modimID
        maskstruc.bsplinefit[1].name =  oldmaskstruc.bsplinefit[1].name
        maskstruc.bsplinefit[1].modim =  oldmaskstruc.bsplinefit[1].modim
        maskstruc.bsplinefit[1].chi2 =  oldmaskstruc.bsplinefit[1].chi2
        maskstruc.bsplinefit[1].nthetaID =  oldmaskstruc.bsplinefit[1].nthetaID
        maskstruc.bsplinefit[1].rbkptID =  oldmaskstruc.bsplinefit[1].rbkptID
        maskstruc.bsplinefit[2].modimID =  oldmaskstruc.bsplinefit[2].modimID
        maskstruc.bsplinefit[2].name =  oldmaskstruc.bsplinefit[2].name
        maskstruc.bsplinefit[2].modim =  oldmaskstruc.bsplinefit[2].modim
        maskstruc.bsplinefit[2].chi2 =  oldmaskstruc.bsplinefit[2].chi2
        maskstruc.bsplinefit[2].nthetaID =  oldmaskstruc.bsplinefit[2].nthetaID
        maskstruc.bsplinefit[2].rbkptID =  oldmaskstruc.bsplinefit[2].rbkptID
        maskstruc.bsplinefit[3].modimID =  oldmaskstruc.bsplinefit[3].modimID
        maskstruc.bsplinefit[3].name =  oldmaskstruc.bsplinefit[3].name
        maskstruc.bsplinefit[3].modim =  oldmaskstruc.bsplinefit[3].modim
        maskstruc.bsplinefit[3].chi2 =  oldmaskstruc.bsplinefit[3].chi2
        maskstruc.bsplinefit[3].nthetaID =  oldmaskstruc.bsplinefit[3].nthetaID
        maskstruc.bsplinefit[3].rbkptID =  oldmaskstruc.bsplinefit[3].rbkptID
        maskstruc.bsplinefit[4].modimID =  oldmaskstruc.bsplinefit[4].modimID
        maskstruc.bsplinefit[4].name =  oldmaskstruc.bsplinefit[4].name
        maskstruc.bsplinefit[4].modim =  oldmaskstruc.bsplinefit[4].modim
        maskstruc.bsplinefit[4].chi2 =  oldmaskstruc.bsplinefit[4].chi2
        maskstruc.bsplinefit[4].nthetaID =  oldmaskstruc.bsplinefit[4].nthetaID
        maskstruc.bsplinefit[4].rbkptID =  oldmaskstruc.bsplinefit[4].rbkptID
        maskstruc.bsplinefit[5].modimID =  oldmaskstruc.bsplinefit[5].modimID
        maskstruc.bsplinefit[5].name =  oldmaskstruc.bsplinefit[5].name
        maskstruc.bsplinefit[5].modim =  oldmaskstruc.bsplinefit[5].modim
        maskstruc.bsplinefit[5].chi2 =  oldmaskstruc.bsplinefit[5].chi2
        maskstruc.bsplinefit[5].nthetaID =  oldmaskstruc.bsplinefit[5].nthetaID
        maskstruc.bsplinefit[5].rbkptID =  oldmaskstruc.bsplinefit[5].rbkptID
        maskstruc.bsplinefit[6].modimID =  oldmaskstruc.bsplinefit[6].modimID
        maskstruc.bsplinefit[6].name =  oldmaskstruc.bsplinefit[6].name
        maskstruc.bsplinefit[6].modim =  oldmaskstruc.bsplinefit[6].modim
        maskstruc.bsplinefit[6].chi2 =  oldmaskstruc.bsplinefit[6].chi2
        maskstruc.bsplinefit[6].nthetaID =  oldmaskstruc.bsplinefit[6].nthetaID
        maskstruc.bsplinefit[6].rbkptID =  oldmaskstruc.bsplinefit[6].rbkptID
        maskstruc.bsplinefit[7].modimID =  oldmaskstruc.bsplinefit[7].modimID
        maskstruc.bsplinefit[7].name =  oldmaskstruc.bsplinefit[7].name
        maskstruc.bsplinefit[7].modim =  oldmaskstruc.bsplinefit[7].modim
        maskstruc.bsplinefit[7].chi2 =  oldmaskstruc.bsplinefit[7].chi2
        maskstruc.bsplinefit[7].nthetaID =  oldmaskstruc.bsplinefit[7].nthetaID
        maskstruc.bsplinefit[7].rbkptID =  oldmaskstruc.bsplinefit[7].rbkptID
      endif

      acsproc_maskstruc_write, maskfile
      
      maskflag = file_test(maskfile) 
      if maskflag then message = "MOFFAT has reduced chi2 ="+strtrim(maskstruc.moffat.chi2,2) else message = "MASK MISSING"
      process->setMessage, "Target "+metastruc_visit[q].uniqname+" Visit "+metastruc_visit[q].visit+' -> '+message, /APPEND
      
      splog, "Target "+metastruc_visit[q].uniqname+" Visit "+metastruc_visit[q].visit+' -> '+message
      
             
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
