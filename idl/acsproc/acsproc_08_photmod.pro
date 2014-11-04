pro acsproc_08_photmod, process, statuscode, count, flatweight=flatweight, selected=selected

  common acsproc_data
  common acsproc_metastruc
  common acsproc_maskstruc

  smartguess = keyword_set(flatweight)   
  fitfunc = ['asb_devauc_multi','asb_sersic_multi']
  fitfunc_label = ['de Vauc','Sersic']

  ; Get the info on each visit (reference exposure) to be processed (correct statuscode) checking for a comb & mask file:
  
  nrows = n_elements(metastruc)
  void = where(metastruc.subexposure eq 1,nvisits)
  maskflag = 4
  if not keyword_set(selected) then begin
    selected = where(metastruc.statuscode eq statuscode and metastruc.combflag eq 1 and metastruc.maskflag eq maskflag-1 and (metastruc.fmaskflag eq 1 or metastruc.jmaskflag eq 1) and metastruc.subexposure eq 1,nvis) 
  endif else begin
    ;selected -= 1
    nvis = n_elements(selected)
  endelse
  
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
      photmodfile = acsproc_datafile(metastruc_visit[q],/photmod,old=1-keyword_set(flatweight))
      bizfile = acsproc_datafile(metastruc_visit[q],/biz)
            
      ;bz2file = strmid(combfile, 0,strpos(combfile,'.fits')) + '.bz2'

      maskfile = acsproc_datafile(metastruc_visit[q],/mask)
      maskstruc = acsproc_get_maskstruc(maskfile, header=mask_header)
  
      xc = round(sxpar(mask_header, 'CRPIX1'))-1L
      yc = round(sxpar(mask_header, 'CRPIX2'))-1L
      
      bspars = maskstruc.bsplinecen.pars

      img = mrdfits(combfile, 1, hdr)
      err = mrdfits(combfile, 2)
      nim = mrdfits(combfile, 3)
      
      ivar = (nim gt 0) * (err gt 0.) / (err^2 + (nim le 0) + (err le 0.))
      
      pixel_mask = (1B-maskstruc.fmask)*(1B-maskstruc.jmask)
      ivar = ivar * pixel_mask
      
      ; Make RA and Dec images:
      aimg = -dpix * (dindgen(2*hw+1) # replicate(1.d0, 2*hw+1) - hw)
      dimg = dpix * (replicate(1.d0, 2*hw+1) # dindgen(2*hw+1) -  hw)
  
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
      ivar = ivar / errsig^2
      
      ; A fit with flat weighting:
      if keyword_set(flatweight) then begin
        med_ivar = median(ivar[where(ivar gt 0.)])
        ivar_weight = med_ivar * (ivar gt 0.)
      endif else ivar_weight = ivar
        
      ; Do the initial deVaucouleurs fit:
      dspars = [10., 2., bspars]
      if keyword_set(smartguess) then begin
          ivar_weight0 = alog(maskstruc.bsplinecen.modim/min(maskstruc.bsplinecen.modim))
          img0 = maskstruc.bsplinecen.modim
          parinfo = replicate({value:0.D, fixed:1}, n_elements(dspars))
          parinfo[*].value = dspars
          parinfo[0:1].fixed = 0
          splog, "MPFIT > SMARTGUESS **************************************************************"
          devpars = mpfit2dfun(fitfunc[0], dimg, aimg, img0, 0.*img0, dspars, weights=ivar_weight0, functargs={psf: psf, rmin: 0.1*dpix}, bestnorm=bestnorm, perror=perror, niter=niter, yfit=yfit00, covar=covar, parinfo=parinfo, dof=dof, /silent)
          dspars = devpars
          dspars[1] = abs(dspars[1])
          splog, "MPFIT > SMARTGUESS completed with REDUCED CHI^2 ="+strtrim(bestnorm/dof,2)
      endif else dspars = [10., 2., bspars]
            
      ivar0 = ivar[xc-hw:xc+hw,yc-hw:yc+hw]
      ivar_weight0 = ivar_weight[xc-hw:xc+hw,yc-hw:yc+hw]
      img0 = img[xc-hw:xc+hw,yc-hw:yc+hw]

      splog, "MPFIT > Initial deVaucouleurs **************************************************************"
      devpars = mpfit2dfun(fitfunc[0], dimg, aimg, $
        img0, 0.*img0, dspars, weights=ivar_weight0, functargs={psf: psf, rmin: 0.1*dpix}, $
        bestnorm=bestnorm, perror=perror, niter=niter, yfit=yfit0, covar=covar)
                
      ndof = long(total(ivar0 gt 0.)) - n_elements(devpars)
      splog, "MPFIT > Initial deVaucouleurs fit completed with REDUCED CHI^2 ="+strtrim(bestnorm/ndof,2)
      yfit = img*0.
      yfit[xc-hw:xc+hw,yc-hw:yc+hw] = yfit0

      ; Make RA and Dec images:
      aimg = -dpix * (findgen(imsize) - xc) # replicate(1.0, imsize)
      dimg = dpix * replicate(1.0, imsize) # (findgen(imsize) - yc)

      if keyword_set(flatweight) then process->setMessage, "PHOTMOD: FLATWEIGHT", /APPEND else process->setMessage, "PHOTMOD: OLD WEIGHT", /APPEND
      message = 'PHOTMOD: '+photmodfile
      process->setMessage, "Target "+metastruc_visit[q].uniqname+" Visit "+metastruc_visit[q].visit+' -> '+message, /APPEND

      ndof = long(total(ivar gt 0.)) - n_elements(devpars)
 
      message = 'PHOTMOD: Initial deVaucouleurs fit chi2/dof='+strtrim(bestnorm/ndof,2)
      process->setMessage, "Target "+metastruc_visit[q].uniqname+" Visit "+metastruc_visit[q].visit+' -> '+message, /APPEND

      ; Make a mask for any outlying junk:
      sigthresh = 4.0
            
      outermask = ((img-yfit) * sqrt(ivar)) ge sigthresh
      outermask = outermask * (ymask eq 0)
      ivar = ivar * (outermask eq 0)

      for i=0,n_elements(fitfunc)-1 do begin
      
        if (i eq 0) then pars = devpars $
        else if (i eq 1) then pars = [devpars, 4.]

        ; Do the smart mask fit for each fitfunc
        splog, "MPFIT > SMARTMASK **************************************************************"
        pars = mpfit2dfun(fitfunc[i], dimg, aimg, $
          img, 0.*img, pars, weights=ivar_weight, functargs={psf: psf, rmin: 0.1*dpix}, $
          bestnorm=bestnorm, perror=perror, niter=niter, yfit=yfit, covar=covar)
            
        message = 'PHOTMOD: Smart mask '+fitfunc[i]+' chi2/dof='+strtrim(bestnorm/ndof,2)
        process->setMessage, "Target "+metastruc_visit[q].uniqname+" Visit "+metastruc_visit[q].visit+' -> '+message, /APPEND
  
        ;;; If you dare...
        ; Subsample the grid and do one last MPFIT:
        subs = 5
        xbase = 0.5 * (2.*findgen(imsize*subs) + 1. - float(subs))/float(subs)
        ybase = 0.5 * (2.*findgen(imsize*subs) + 1. - float(subs))/float(subs)
        abig = bilinear(aimg, xbase, ybase)
        dbig = bilinear(dimg, xbase, ybase)
        xbase = 0
        ybase = 0
        
        splog, "MPFIT > "+fitfunc_label[i]+" **************************************************************"
        pars = mpfit2dfun(fitfunc[i], dbig, abig, $
          img, 0.*img, pars, weights=ivar_weight, $
          functargs={psf: psf, rmin: 0.1*dpix, subs: subs}, $
          bestnorm=bestnorm, perror=perror, niter=niter, yfit=yfit, covar=covar)
        ; Keep it positive!
        pars[1] = abs(pars[1])
          
        message = 'PHOTMOD: Subsampled (x'+strtrim(subs,2)+') grid '+fitfunc[i]+' chi2/dof='+strtrim(bestnorm/ndof,2)
        process->setMessage, "Target "+metastruc_visit[q].uniqname+" Visit "+metastruc_visit[q].visit+' -> '+message, /APPEND
        
        ;-------------------------------------------------------------------------
              
        splog, "REDUCED CHI^2 of "+fitfunc_label[i]+" fit="+strtrim(bestnorm/ndof,2)

        if fitfunc[i] eq 'asb_devauc_multi' then photimg = asb_devauc_multi(dimg,aimg,pars,rmin=0.1*dpix,psf=psf) $
        else if fitfunc[i] eq 'asb_sersic_multi' then photimg = asb_sersic_multi(dimg,aimg,pars,rmin=0.1*dpix,psf=psf) $
        else photimg = img*0.0
   
        photmodstruc = {fitfunc: fitfunc[i], img:photimg, pars: pars, perror: perror, $
          covar: covar, chisq: bestnorm, ndof: ndof, xc: xc, yc: yc}
  
        abig = 0
        dbig = 0
        
        if (i eq 0) then begin
        
          ; Turn this all into proper band photometry:
          
          ; Get the deVaucouleurs fit parameters and the photometric
          ; header keyword parameters:
          photflam = sxpar(hdr, 'PHOTFLAM')
          photplam = sxpar(hdr, 'PHOTPLAM')
          
          ; Convert deVaucouleurs model amplitude to counts, and to AB magnitude:
          k_dev = 7.66925001
          cps_dev = ((2. * !pi * pars[0] * pars[1]^2) * 4. * gamma(8.) / k_dev^8.) / (dpix^2)
          abmag_zpt = -2.5 * alog10(photflam) - 21.10 - 5. * alog10(photplam) + 18.6921
          abmag_dev = -2.5 * alog10(cps_dev) + abmag_zpt
          
          cps_deverr = cps_dev * (perror[0] / pars[0])
          acs_rdev = pars[1]
          acs_rdeverr = perror[1]
          acs_qdev = pars[4]
          acs_qdeverr = perror[4]
          acs_padev = pars[5]
          acs_padeverr = perror[5]
          
          if (acs_qdev gt 1.) then begin
            acs_qdev = 1./acs_qdev
            acs_qdeverr = acs_qdeverr * acs_qdev^2
            acs_padev = acs_padev + 90.
          endif
          pa = acs_padev mod 180.
          ptest = pa ge 0.
          pa = pa * ptest + (180. - abs(pa)) * (1B - ptest)
          acs_padev = pa
                    
          mwrfits, photmodstruc, photmodfile, /create, /silent
          
        endif else if (i eq 1) then begin
        
          acs_nser = pars[6]
          acs_nsererr = perror[6]

          mwrfits, photmodstruc, photmodfile, /silent
          
        endif
       
      endfor 

      ;-------------------------------------------------------------------------
      ; Breakpoints for b-spline photometry:
      rbkpt = dpix*[1., 2., 4., 6., 8., 12., 16., 22., 30., 50., 60., 100., 150., 200., 300., 400., 500., 600., 800.]
      
      message = 'PHOTMOD: B-Spline Photometry'
      process->setMessage, "Target "+metastruc_visit[q].uniqname+" Visit "+metastruc_visit[q].visit+' -> '+message, /APPEND
      ; Do the b-spline photometry:
      bfit = bspline_ellip(dspars[2:*], x=dimg, y=aimg, invvar=ivar, data=img, psf=psf, ntheta=0, rbkpt=rbkpt, sset=sset)
      
      ; Store the bspline photometry of interest.
      ; Really we're just evaluating the ellipsoidal
      ; b-spline surface brightness on a baseline.
      ; Anything else we want to do can be done later.

      bsr_arcs = (findgen(800) + 1.) / 100.
      bsflux = bsr_arcs
      bsflux = bspline_valu(bsr_arcs, sset)

      ;-------------------------------------------------------------------------

      photmodstruc = {photflam: photflam, photplam: photplam, abmag_zpt:abmag_zpt, abmag_dev: abmag_dev, $
         cps_dev: cps_dev, cps_deverr: cps_deverr, acs_rdev: acs_rdev, acs_rdeverr: acs_rdeverr, acs_qdev: acs_qdev, $
         acs_qdeverr: acs_qdeverr, acs_padev: acs_padev, acs_padeverr: acs_padeverr, acs_nser: acs_nser, acs_nsererr: acs_nsererr, $
         bsr_arcs: bsr_arcs, bsflux: bsflux, outermask:outermask}
      
      mwrfits, photmodstruc, photmodfile, /silent
      
      acsproc_database_update_table, 'photmod', metastruc_visit[q].proposid, metastruc_visit[q].visit, file_basename(photmodfile), confirm=confirm, message=message
       if confirm then splog, 'Database > '+message else splog, 'update failed.'
      
      ;pos = strpos(combfile,'/',/reverse_search) - 2
      ;dir = strmid(combfile,0,pos)
      ;cd, dir
      ;spawn, 'tar -cjf '+bz2file+' '+strmid(combfile,pos)+' '+strmid(maskfile,pos)+' '+strmid(photmodfile,pos)+' '+strmid(bizfile,pos)

      has_photmod = maskstruc.has_photmod
      maskstruc.has_photmod = file_test(photmodfile)
      if has_photmod ne maskstruc.has_photmod then acsproc_maskstruc_write, maskfile

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
