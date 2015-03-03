;BEGIN ORIGINAL DOC:
; singlerect.pro
;
; make rectified versions of individual frames, along with
; their error and mask images, AND along with their PSFs.
;
; abolton@cfa 2007apr
;END ORIGINAL DOC

pro acsproc_03_rectify, process, statuscode, count, selected=selected

  common acsproc_data
  common acsproc_metastruc

  ; The size of images to make:
  xsiz = imsize
  ysiz = imsize
  dpix = 0.05d0
  ra_grid = dpix * (reverse(dindgen(xsiz)) # replicate(1.d0, ysiz) - xsiz / 2)
  dec_grid = dpix * (replicate(1.d0, xsiz) # dindgen(ysiz) - ysiz / 2 + 1)
  
  whref = where((ra_grid eq 0) and (dec_grid eq 0))
  xref = whref[0] mod xsiz
  yref = whref[0] / xsiz
  
    good = (strtrim(metastruc.expflag, 2) eq 'NORMAL') * (metastruc.exptime gt 1.)
    maybe = ((strtrim(metastruc.expflag, 2) eq 'EXCESSIVE DOWNTIME') or (strtrim(metastruc.expflag, 2) eq 'INTERRUPTED')) * (metastruc.exptime gt 1.)
    ready = metastruc.statuscode eq statuscode and metastruc.cutflag eq 1
  void = where((good or maybe) and ready,nexp)
  
  process->Start, steps=nexp
  metastruc_changed = 0L
                        
  ; Get the info on each visit (reference exposure) to be processed (correct statuscode) with a cutout file:
  
  nrows = n_elements(metastruc)
  if not keyword_set(selected) then selected = where(metastruc.statuscode eq statuscode and metastruc.cutflag eq 1 and metastruc.subexposure eq 1,nvis)  else nvis = n_elements(selected)
  
  process->Setmessage, 'Processing Proposal ID '+program+' ('+strtrim(nexp,2)+'/'+strtrim(nrows,2)+' exposures).', /label
  process->Setmessage, ''  
  
  if (nvis gt 0) then begin

    metastruc_visit = metastruc[selected]
    
    for q = 0,n_elements(metastruc_visit)-1 do begin
          
      if process->CheckCancel() then begin
        process->Step, nexp
        process->Setmessage, 'Cancelling...', /append
        if metastruc_changed then acsproc_metastruc_update
        process->Destroy
        return
      endif

      visit = metastruc_visit[q].visit
      this_visit = metastruc.visit eq visit
            
      wh = where((good or maybe) and this_visit,nexp)
      
      if (nexp gt 0) then begin
      
        subexposure = metastruc[wh]
              
        cutoutfile = acsproc_datafile(subexposure, /cutout)
                  
        rectfiles = strarr(nexp)
        
        ; Loop over all exposures per visit:
        
        for i = 0L, nexp - 1 do begin
        
        process->setMessage, "Target "+subexposure[i].uniqname+" Visit "+subexposure[i].visit+' RECTIFY: Subexposure '+strtrim(subexposure[i].subexposure,2)+' of '+strtrim(nexp,2), /APPEND
        process->Increment, 1.0
        
        ; Get the data:
            hdr0 = headfits(cutoutfile[i],exten=0)
            data = mrdfits(cutoutfile[i],1,hdr1)
            errs = mrdfits(cutoutfile[i],2)
            mask = mrdfits(cutoutfile[i],4)
            nx = (size(data))[1]
            ny = (size(data))[2]
                
        ; Get the pixels for the rectified map:
            acssip_ad2xy, ra_grid + subexposure[i].ra_tan_targ, $
              dec_grid + subexposure[i].dec_tan_targ, hdr1, xpix, ypix
        
        
        ; Make the interpolated images
            rect_data = bilinear(data, xpix, ypix, missing=0)
            rect_errs = bilinear(errs, xpix, ypix, missing=0)
            rect_mask = bilinear(float(mask gt 0), xpix, ypix, missing=4.)
        
        
        ; Update the WCS:
            sxaddpar, hdr1, 'CRPIX1', float(xref+1)
            sxaddpar, hdr1, 'CRPIX2', float(yref+1)
            sxaddpar, hdr1, 'CRVAL1', sxpar(hdr0, 'RA_TARG')
            sxaddpar, hdr1, 'CRVAL2', sxpar(hdr0, 'DEC_TARG')
            sxaddpar, hdr1, 'CD1_2', 0.d0
            sxaddpar, hdr1, 'CD2_1', 0.d0
            sxaddpar, hdr1, 'CD1_1', -dpix / 3600.d0
            sxaddpar, hdr1, 'CD2_2', dpix / 3600.d0
            sxaddpar, hdr1, 'ORIENTAT', 0.d0
            sxaddpar, hdr1, 'VAFACTOR', 1.d0
            sxdelpar, hdr1, 'RA_APER'
            sxdelpar, hdr1, 'DEC_APER'
            sxdelpar, hdr1, 'PA_APER'
            sxdelpar, hdr1, 'A_ORDER'
            sxdelpar, hdr1, 'B_ORDER'
            for ii = 0, 4 do for jj = 0, 4 do sxdelpar, hdr1, $
              'A_' + string(ii, format='(i1)') + '_' + string(jj, format='(i1)')
            for ii = 0, 4 do for jj = 0, 4 do sxdelpar, hdr1, $
              'B_' + string(ii, format='(i1)') + '_' + string(jj, format='(i1)')
        
        
            ; A good test:
            ;acssip_ad2xy, ra_grid, dec_grid, hdr1, xtest, ytest
        
            psf_file = getenv('HST_PSFDIR') + '/acs_' + acs_filter[subexposure[0].filtcode] + '_lens_1x.fits'
            psf_full = mrdfits(psf_file,/silent)
        
            nxp = (size(psf_full))[1]
            nyp = (size(psf_full))[2]
            junk = max(psf_full, pmax)
            cxp = pmax mod nxp
            cyp = pmax / nxp
        
            xmap_psf = xpix[xref-nxp/2:xref+nxp/2,yref-nyp/2:yref+nyp/2,*] - xpix[xref,yref] + cxp
            ymap_psf = ypix[xref-nxp/2:xref+nxp/2,yref-nyp/2:yref+nyp/2,*] - ypix[xref,yref] + cyp
        
            psf_out = bilinear(psf_full, xmap_psf, ymap_psf, missing=0)
            psf_out = psf_out / total(psf_out)
        
        ; Write it out:
            outroot = (strsplit(cutoutfile[i], '_cutout.fits', /extract, /regex))[0]
            ofile = outroot + '_recti.fits'
            rectfiles[i] = ofile
            mwrfits, 0, ofile, hdr0, /create, /silent
            mwrfits, rect_data, ofile, hdr1, /silent
            mwrfits, rect_errs, ofile, /silent
            mwrfits, (rect_mask eq 0), ofile, /silent
            mwrfits, psf_out, ofile, /silent
            
            r = wh[i]
            rectiflag = file_test(acsproc_datafile(metastruc[r], /rectify))
            if (not metastruc_changed) then if (metastruc[r].rectiflag ne rectiflag) then metastruc_changed = 1L
            metastruc[r].rectiflag = rectiflag
        
        endfor
        
      endif
      
    endfor
    
  endif
  
  if metastruc_changed then  acsproc_metastruc_update
  process->Destroy  
  
end
