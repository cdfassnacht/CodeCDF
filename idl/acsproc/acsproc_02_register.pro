;
; acsproc_02_register.pro
;
; registr images through constrained simultaneous MPFITting
; of a suitably general image profile to the lens galaxy.
;
; abolton@cfa 2007apr
; 
;       Refactored into a subprogram of acsproc for GUI management
;       By: Joel Brownstein, October, 2010
; 
;

pro acsproc_02_register, process, statuscode, count, selected=selected

  common acsproc_data
  common acsproc_metastruc  

  chw = 100  ; cross-correlation half-width
  fhw = 50   ; fitting half-width
  shwin = 30 ; window of allowable pixel shifts
     
  nrows = n_elements(metastruc)
            
  wh = where(metastruc.statuscode eq statuscode and metastruc.cutflag eq 1,nexp)
  
  process->Start, steps=nexp
  metastruc_changed = 0L

  process->Setmessage, 'Processing Proposal ID '+program+' ('+strtrim(nexp,2)+'/'+strtrim(nrows,2)+' exposures).', /label
  process->Setmessage, ''
  
  if (nexp gt 0) then begin
  
    n=1

    ; Get the data stacks for the cross-correlation:
    nx = 1500
    ny = 1500
    xfull = dindgen(nx) # replicate(1.d0, ny)
    yfull = replicate(1.d0, nx) # dindgen(ny)
    
    expj = 0L
    for j = 0L, nrows-1 do begin
    
      if process->CheckCancel() then begin
        process->Step, nexp
        process->Setmessage, 'Cancelling...', /append
        if metastruc_changed then acsproc_metastruc_update
        process->Destroy
        return
      endif
    
      err = 0
      if (metastruc[j].statuscode gt statuscode and metastruc[j].cutflag eq 1 and keyword_set(metastruc0)) then begin ;check to see if this row has been already been processed
        wh = where(metastruc.visit eq metastruc0.visit and $
                   metastruc.exposure eq metastruc0.exposure and $
                   metastruc.subexposure eq metastruc0.subexposure, c)
        if (c eq 1) then metastruc[j] = metastruc0[(wh)[0]] else err = 1
        
      endif else if (metastruc[j].statuscode eq statuscode and metastruc[j].cutflag eq 1) then begin ;process (set) normally
        if (metastruc[j].subexposure eq 1) then begin ;reference exposure
          process->Increment, 1.0
          process->setMessage, "Target "+metastruc[j].uniqname+" Visit "+metastruc[j].visit+' MPREGISTER: '+strtrim(n,2)+' of '+strtrim(nexp,2)+' [reference]', /APPEND
          cname1 = dataroot +  metastruc[j].visit + '/' + metastruc[j].rootname + '_flt_cutout.fits'
          head = headfits(cname1, exten=1)
          this_xc = fix(sxpar(head, 'LCOORDX'))
          this_yc = fix(sxpar(head, 'LCOORDY'))
          if (this_xc ne 0 and this_yc ne 0) then begin
            good = (strtrim(metastruc.expflag, 2) eq 'NORMAL') * (metastruc.exptime gt 1.)
            maybe = ((strtrim(metastruc.expflag, 2) eq 'EXCESSIVE DOWNTIME') or (strtrim(metastruc.expflag, 2) eq 'INTERRUPTED')) * (metastruc.exptime gt 1.)
            this_visit =  metastruc.visit eq metastruc[j].visit and metastruc.exposure eq metastruc[j].exposure
            expset = where(this_visit and (good or maybe), nf) ;get full expset for this reference exposure
            if (nf gt 0) then begin
              this_ref = (where(metastruc[expset].subexposure eq 1))[0] ;identify the refence/set
              files = dataroot +  metastruc[expset].visit + '/' + metastruc[expset].rootname + '_flt_cutout.fits'
              
              ; Do the cross-correlation registration to the nearest pixel (if we have more than one image!):
              if (nf gt 1) then begin
                  datastack = fltarr(2*chw+1, 2*chw+1, nf)
                  crmstack = fltarr(2*chw+1, 2*chw+1, nf)
                  medstack = fltarr(2*chw+1, 2*chw+1, nf)
                  for i = 0, nf-1 do begin
                      splog, "READING: "+strtrim(files[i],2)
                      datastack[*,*,i] = (mrdfits(files[i], 1, /silent))[this_xc-chw:this_xc+chw,this_yc-chw:this_yc+chw]
                      crmstack[*,*,i] = (mrdfits(files[i], 4, /silent))[this_xc-chw:this_xc+chw,this_yc-chw:this_yc+chw]
                      medstack[*,*,i] = (mrdfits(files[i], 5, /silent))[this_xc-chw:this_xc+chw,this_yc-chw:this_yc+chw]
                  endfor
                  imgstack = datastack * (crmstack eq 0) + medstack * (crmstack ne 0)
                  fft_register, imgstack, shwin, xsh, ysh, /nofrac
                  xc_all = this_xc + xsh[*,0]
                  yc_all = this_yc + ysh[*,0]
              endif else begin
                  xc_all = this_xc
                  yc_all = this_yc
              endelse
          
          
              ; Now get the subimages for centroiding:
              datastack = fltarr(2*fhw+1, 2*fhw+1, nf)
              errstack = fltarr(2*fhw+1, 2*fhw+1, nf)
              crmstack = fltarr(2*fhw+1, 2*fhw+1, nf)
              medstack = fltarr(2*fhw+1, 2*fhw+1, nf)
              rastack = fltarr(2*fhw+1, 2*fhw+1, nf)
              decstack = fltarr(2*fhw+1, 2*fhw+1, nf)
              for i = 0, nf-1 do begin
                  datastack[*,*,i] = (mrdfits(files[i], 1, hdr, /silent))[xc_all[i]-fhw:xc_all[i]+fhw,yc_all[i]-fhw:yc_all[i]+fhw]
                  errstack[*,*,i] = (mrdfits(files[i], 2, /silent))[xc_all[i]-fhw:xc_all[i]+fhw,yc_all[i]-fhw:yc_all[i]+fhw]
                  crmstack[*,*,i] = (mrdfits(files[i], 4, /silent))[xc_all[i]-fhw:xc_all[i]+fhw,yc_all[i]-fhw:yc_all[i]+fhw]
                  medstack[*,*,i] = (mrdfits(files[i], 5, /silent))[xc_all[i]-fhw:xc_all[i]+fhw,yc_all[i]-fhw:yc_all[i]+fhw]
                  xsub = xfull[xc_all[i]-fhw:xc_all[i]+fhw,yc_all[i]-fhw:yc_all[i]+fhw]
                  ysub = yfull[xc_all[i]-fhw:xc_all[i]+fhw,yc_all[i]-fhw:yc_all[i]+fhw]
                  acssip_xy2ad, xsub, ysub, hdr, asub, dsub
                  rastack[*,*,i] = asub
                  decstack[*,*,i] = dsub
              endfor
          
              ivarstack = (crmstack eq 0) / errstack^2
          
          
              ; Initialize some moffat parameters:
          
              moffati = 500.              ; central intensity
              moffatr = .25               ; scale radius
              moffatg = 1.                ; "gamma" exponent
              moffatq = 0.85              ; axis ratio
              moffata = 15.               ; position angle
          
              moffatpar = fltarr(7, nf)
              moffatpar[0,*] = moffati
              moffatpar[1,*] = moffatr
              moffatpar[2,*] = moffatg
              moffatpar[3,*] = decstack[fhw,fhw,*]
              moffatpar[4,*] = rastack[fhw,fhw,*]
              moffatpar[5,*] = moffatq
              moffatpar[6,*] = moffata
          
              if (nf gt 1) then begin
                  tie_vec = 'P[' + strtrim(sindgen(7), 2) + ']'
                  tied = strarr(7, nf)
                  for i = 0, nf-1 do tied[*,i] = tie_vec
                  tied[*,0] = ''
                  tied[3:4,*] = ''
                  pinfo = replicate({tied: ''}, 7L*nf)
                  pinfo.tied = tied[*]
              endif else pinfo = replicate({tied: ''}, 7L)
          
          
              par = moffatpar[*]
              fpar = mpfit('mpmoffatstack', par, functargs = $
                           {x: decstack, y: rastack, data: datastack, ivar: ivarstack, deviates: 1B}, $
                           parinfo = pinfo, perror = perror, /quiet)
              modim = mpmoffatstack(fpar, x=decstack, y=rastack)
          
              ; Select out the params of interest:
              mo_idx = 7L * lindgen(nf)
              if (not metastruc_changed) then metastruc_changed = 1L
              metastruc[expset].ra_tan_targ = fpar[mo_idx + 4]
              metastruc[expset].dec_tan_targ = fpar[mo_idx + 3]
              metastruc[expset].moffati_targ = fpar[mo_idx]
              metastruc[expset].moffatr_targ = fpar[mo_idx + 1]
              metastruc[expset].moffatg_targ = fpar[mo_idx + 2]
              metastruc[expset].moffatq_targ = fpar[mo_idx + 5]
              metastruc[expset].moffata_targ = fpar[mo_idx + 6]                                
              
            endif else err = 2
                        
          endif else err = 3
          
        endif else begin
          ; skip over exposures other than the reference (since they are co-treated in the stack)
          process->setMessage, "Target "+metastruc[j].uniqname+" Visit "+metastruc[j].visit+' MPREGISTER: '+strtrim(n,2)+' of '+strtrim(nexp,2)+' [stacked on reference]', /APPEND
        endelse

        n+=1

      endif else err = 4
      
      if err then print, "***** error ("+strtrim(err,2)+") processing proposid "+program+" target "+metastruc[j].targname+" visit "+metastruc[j].visit+" subexposure "+strtrim(metastruc[j].subexposure,2)
      
    endfor
    
  endif
       
  if metastruc_changed then acsproc_metastruc_update    
  process->Destroy  

end
