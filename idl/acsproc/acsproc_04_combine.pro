;BEGIN ORIGINAL DOCS
; combine_them.pro
;
; combine the rectified and registered images into
; single co-added frames!
;
; abolton@cfa 2007apr
;END ORIGINAL DOCS

pro acsproc_04_combine, process, statuscode, count, selected=selected

  common acsproc_data
  common acsproc_metastruc
  
  void = where(metastruc.statuscode eq statuscode and metastruc.rectiflag eq 1,nexp)
  
  process->Start, steps=nexp
  metastruc_changed = 0L
    
  ; Get the info on each visit (reference exposure) to be processed (correct statuscode) checking for a recti file:
  
  nrows = n_elements(metastruc)
  void = where(metastruc.subexposure eq 1,nvisits)
  if not keyword_set(selected) then selected = where(metastruc.statuscode eq statuscode and metastruc.rectiflag eq 1 and metastruc.subexposure eq 1,nvis) else nvis = n_elements(selected)
  
  process->Setmessage, 'Processing Proposal ID '+program+' ('+strtrim(nvis,2)+'/'+strtrim(nvisits,2)+' Visits)', /label
  process->Setmessage, ''  
  
  if (nvis gt 0) then begin

    metastruc_visit = metastruc[selected]
    
    for q = 0,nvis-1 do begin
          
      if process->CheckCancel() then begin
        process->Step, nexp
        process->Setmessage, 'Cancelling...', /append
        if metastruc_changed then acsproc_metastruc_update
        process->Destroy
        return
      endif

      visit = metastruc_visit[q].visit
      good = (strtrim(metastruc.expflag, 2) eq 'NORMAL') * (metastruc.exptime gt 1.)
      maybe = ((strtrim(metastruc.expflag, 2) eq 'EXCESSIVE DOWNTIME') or (strtrim(metastruc.expflag, 2) eq 'INTERRUPTED')) * (metastruc.exptime gt 1.)
      this_visit = metastruc.visit eq visit
            
      wh = where((good or maybe) and this_visit,nexp)
      
      if (nexp gt 0) then begin
      
        subexposure = metastruc[wh]
              
        rectifile = acsproc_datafile(subexposure,/rectify)
        combfile = acsproc_datafile(subexposure,/combine) 
        maskfile = acsproc_datafile(subexposure,/mask)

        ; Get the individual frame info:
        hdr0 = headfits(rectifile[0])
        hdr1 = headfits(rectifile[0], exten=1)
        nx = sxpar(hdr1, 'NAXIS1')
        ny = sxpar(hdr1, 'NAXIS2')
        phdr = headfits(rectifile[0], exten=4)
        nxp = sxpar(phdr, 'NAXIS1')
        nyp = sxpar(phdr, 'NAXIS2')
        
        ; Get exposure times:
        exptimes = subexposure.exptime
        
        ; Different handling for multiple vs. single exposures:
        
        if (nexp gt 1) then begin
    
            data = fltarr(nx,ny,nexp)
            mdim = fltarr(nx,ny,nexp)
            errs = fltarr(nx,ny,nexp)
            mask = bytarr(nx,ny,nexp)
            psfs = fltarr(nxp,nyp,nexp)
    
    
            medbox = 11             ; to help w/ more CR flagging
            
            for j  = 0L, nexp-1 do begin
                process->setMessage, "Target "+subexposure[j].uniqname+" Visit "+subexposure[j].visit+' COMBINE: Subexposure '+strtrim(subexposure[j].subexposure,2)+' of '+strtrim(nexp,2), /APPEND
                process->Increment, 1.0
            
                data[*,*,j] = mrdfits(rectifile[j],1,/silent) / exptimes[j]
                mdim[*,*,j] = median(data[*,*,j], medbox)
                errs[*,*,j] = mrdfits(rectifile[j],2,/silent) / exptimes[j]
                mask[*,*,j] = mrdfits(rectifile[j],3,/silent)
                psfs[*,*,j] = mrdfits(rectifile[j],4,/silent)
            endfor
            
           ; individual "pretty" images
            pimg = data * (mask ne 0) + mdim * (mask eq 0)
    
            ; Overall median image:
            tmed = median(pimg, dimen=3, /even)
    
            ; Sigma-thresholding for missed CRs:
            sigdiff = 0. * data
            for j = 0L, nexp-1 do sigdiff[*,*,j] = $
              ((data[*,*,j] - tmed) * mask[*,*,j]) / errs[*,*,j]
            siglim = 4.
            mask = mask * (sigdiff lt siglim)
    
            ; Total CR weight:
            mask_tot = total(mask, 3)
    
            ; Combined images and errors:
            data_tot = total(data * mask, 3) / (mask_tot > 1.)
            errs_tot = sqrt(total(errs^2*mask, 3)) / (mask_tot > 1.)
    
            ; Combined PSF:
            opsf = total(psfs, 3)
            opsf = opsf / total(opsf)
    
        endif else begin            ; what we do for single exposures:
            process->setMessage, "Target "+subexposure[0].uniqname+" Visit "+subexposure[0].visit+' COMBINE: Subexposure '+strtrim(subexposure[0].subexposure,2)+' of '+strtrim(nexp,2), /APPEND
            process->Increment, 1.0
            
            mask_tot = byte(mrdfits(rectifile[0],3,/silent))
            data_tot = float(mask_tot * mrdfits(rectifile[0],1,/silent) / exptimes[0])
            errs_tot = float(mask_tot * mrdfits(rectifile[0],2,/silent) / exptimes[0])
            opsf = float(mrdfits(rectifile[0],4,/silent))
        endelse
    
        ; Update header zero:
        ttime = total(exptimes)
        sxaddpar, hdr0, 'TOT_TIME', ttime, 'total co-added exposure time'
        for j = 0, nexp-1 do sxaddpar, hdr0, 'SUBIM_'+string(j,format='(i2.2)'), $
          subexposure[0].rootname, 'component images'
        sxdelpar, hdr0, 'SUBR_XLO'
        sxdelpar, hdr0, 'SUBR_XHI'
        sxdelpar, hdr0, 'SUBR_YLO'
        sxdelpar, hdr0, 'SUBR_YHI'          
        
        ; Write out:
        mwrfits, 0, combfile, hdr0, /create, /silent
        mwrfits, data_tot, combfile, hdr1, /silent
        mwrfits, errs_tot, combfile, /silent
        mwrfits, fix(mask_tot), combfile, /silent
        mwrfits, opsf, combfile, /silent
        
        combflag = file_test(combfile)
        maskflag = file_test(maskfile)
        
        for s=0,nexp-1 do begin          
          r = wh[s]
          if (not metastruc_changed) then if (metastruc[r].combflag ne combflag) then metastruc_changed = 1L
          metastruc[r].combflag = combflag
          if combflag then begin
            spos = strpos(combfile, '/', /reverse_search)
            acsproc_database_update, 'exposure', 'redux_name', strmid(combfile, spos+1), metastruc[r].visit, strtrim(metastruc[r].exposure,2), confirm=confirm, message=message
            if confirm then message += 'Database > update exposure redux_name successful.' else message += 'update exposure redux_name failed.'
            process->setMessage, message, /append
           acsproc_database_update_table, 'exposure', metastruc[r].proposid, metastruc[r].visit, file_basename(combfile), confirm=confirm, message=message
           if confirm then splog, 'Database > '+message else splog, 'update failed.'
          endif
          if (not metastruc_changed) then if (metastruc[r].maskflag ne maskflag) then metastruc_changed = 1L
          metastruc[r].maskflag = maskflag 
        endfor

      endif
      
    endfor
    
  endif
  
  if metastruc_changed then acsproc_metastruc_update
  process->Destroy  

end