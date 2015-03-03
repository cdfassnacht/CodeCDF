;+
;
; NAME: acsproc_01_cutout
;
; PURPOSE: make SLACS ACS-image cutouts, with computed CR masks
;   and sky subtraction.
;
; WRITTEN: abolton@cfa 2006apr.
;   Took out sky subtraction, took off "slacs_" prefix 2006july
;   Put sky subtraction back in, 2007mar
;   
; Refactored into a subprogram of acsproc for GUI management
; By: Joel Brownstein, October, 2010
;-

PRO acsproc_01_cutout, process, statuscode, count, selected=selected
  
  common acsproc_data
  common acsproc_metastruc  

  ; Keep with the snapshot convention:
  rad = 750
  xrange = 2048 + [0, 2*rad-1] - rad
  yrange = 1024 + [0, 2*rad-1] - rad
  
  medrad=10
  xr2 = [xrange[0]-2*medrad,xrange[1]+2*medrad]
  yr2 = [yrange[0]-2*medrad,yrange[1]+2*medrad]
  nx2 = xr2[1]-xr2[0]+1
  ny2 = yr2[1]-yr2[0]+1
  
  filename = dataroot + metastruc.visit + '/' + metastruc.filename
  
  froot = acsproc_datafile(metastruc, /froot)

  cutoutfile = acsproc_datafile(metastruc, /cutout)  
     
  wh = where(metastruc.statuscode eq statuscode,nexp)
  
  process->Start, steps=nexp
  metastruc_changed = 0L
    
  nrows = n_elements(metastruc)
   
  process->Setmessage, 'Processing Proposal ID '+program+' ('+strtrim(nexp,2)+'/'+strtrim(nrows,2)+' exposures).', /label
  process->Setmessage, ''
  
  if (nexp gt 0) then begin
  
    n = 1
    for j = 0L, nrows-1 do begin
        
      if (metastruc[j].statuscode eq statuscode) then begin
    
        if process->CheckCancel() then begin
          process->Step, nexp
          process->Setmessage, 'Cancelling...', /append
          if metastruc_changed then acsproc_metastruc_update    
          process->Destroy
          return
        endif else begin
          good = (strtrim(metastruc[j].expflag, 2) eq 'NORMAL') * (metastruc[j].exptime gt 1.)
          maybe = ((strtrim(metastruc[j].expflag, 2) eq 'EXCESSIVE DOWNTIME') or (strtrim(metastruc[j].expflag, 2) eq 'INTERRUPTED')) * (metastruc[j].exptime gt 1.)
          if good then checkflag = 'PROCESSING' else if maybe then checkflag = 'PROCESSING (?)' else CHECKFLAG = 'SKIPPING'
          process->Setmessage, checkflag + ' Visit:'+metastruc[j].visit + $
                               ' ['+strtrim(metastruc[j].expflag,2) + $
                               ' '+strtrim(metastruc[j].exptime,2) + $
                               's] '+filename[j]+' ('+strtrim(n,2)+' of '+strtrim(nexp,2)+')', /append
          process->Increment, 0.25
          if not maybe and (strtrim(metastruc[j].expflag, 2) eq 'EXCESSIVE DOWNTIME') then splog, "WARNING: "+strtrim(metastruc[j].expflag,2) + $
                               ' '+strtrim(metastruc[j].exptime,2) + $
                               's for '+filename[j]
        endelse
    
        wf2 = metastruc[j].aperture eq "WFC2"
        extoff = Keyword_set(wf2) ? 1 : 4
        splog, "Looking for aperture="+metastruc[j].aperture+" in ext="+strtrim(extoff,2)
        
        if (good or maybe) then begin
        
          
          splog, "READ: "+filename[j]
          hdr0 = headfits(filename[j])
          img = (mrdfits(filename[j], extoff, hdr1))[xr2[0]:xr2[1],yr2[0]:yr2[1]]
          err = (mrdfits(filename[j], extoff+1))[xr2[0]:xr2[1],yr2[0]:yr2[1]]
          dqm = (mrdfits(filename[j], extoff+2))[xr2[0]:xr2[1],yr2[0]:yr2[1]]
          
          ; Change the headers:
          crpix1 = sxpar(hdr1, 'CRPIX1') - xrange[0]
          crpix2 = sxpar(hdr1, 'CRPIX2') - yrange[0]
          ocrpix1 = sxpar(hdr1, 'OCRPIX1') - xrange[0]
          ocrpix2 = sxpar(hdr1, 'OCRPIX2') - yrange[0]

          sxaddpar, hdr1, 'CRPIX1', crpix1
          sxaddpar, hdr1, 'CRPIX2', crpix2
          sxaddpar, hdr1, 'OCRPIX1', ocrpix1
          sxaddpar, hdr1, 'OCRPIX2', ocrpix2
          
          ; Mask the negative pixels:
          negsig = 4.5
          djs_iterstat, img, mean=imgmean, median=imgmed, sigma=imgsig
          errmed = Median(err)
          mfimg = Median(img, 2*medrad+1)
          nmask = (img - imgmed) LT (-Abs(negsig)*errmed)
          
          img2 = img * (1B-nmask) + mfimg * nmask
           
          ; Write out a file for la_cosmic to work on:
          oimg = img2[medrad:nx2-medrad-1,medrad:ny2-medrad-1]
          tempname = froot[j]+'_temp.fits'
          crname = froot[j]+'_temp-mask.fits'
          mwrfits, oimg, tempname, /create
          
          ; Run la_cosmic:
          gain = Float(sxpar(hdr0, 'CCDGAIN'))
          rn_a = sxpar(hdr0, 'READNSEA')
          rn_b = sxpar(hdr0, 'READNSEB')
          rn_c = sxpar(hdr0, 'READNSEC')
          rn_d = sxpar(hdr0, 'READNSED')
          readn = Float(Round(mean([rn_a, rn_b, rn_c, rn_d])))
          
          ;sigfrac = 0.6  for 10587
          ;sigfrac = 0.2  for 10174
          sigfrac = Round(gain) EQ 2 ? 0.6 : 0.2
          objlim = 2.5
          niter = 5
          sigclip = 4.5
          
          process->Increment, 0.25
          process->Setmessage, 'running la_cosmic...', /append
          la_cosmic, tempname, gain=gain, readn=readn, niter=niter, sigfrac=sigfrac, objlim=objlim, sigclip=sigclip, verbose=0
          process->Increment, 0.25
          process->Setmessage, 'Exit la_cosmic', /append
            
          crmask = mrdfits(crname, /silent)
          
          ;atv, oimg
          ;atvplot, where(crmask) mod 321, where(crmask) / 321, ps=1

          nx3 = (Size(crmask))[1]
          ny3 = (Size(crmask))[2]
          
          ;for cycle > 17, get mdsky fron extention 0 (previously was extension 1)
          
          mdsky = sxpar(hdr0, 'MDRIZSKY', count=has_mdsky)
          if not keyword_set(has_mdsky) then mdsky = sxpar(hdr1, 'MDRIZSKY', count=has_mdsky)
          if not keyword_set(has_mdsky) then message, 'MDRIZSKY not found in ext=0 or ext='+strtrim(extoff)
          if (not metastruc_changed) then if (metastruc[j].mdrizsky ne mdsky) then metastruc_changed = 1L
          metastruc[j].mdrizsky = mdsky
          
          sxaddpar, hdr0, 'SUBR_XLO', xrange[0], 'ASB 0-based subraster xlo'
          sxaddpar, hdr0, 'SUBR_XHI', xrange[1], 'ASB 0-based subraster xhi'
          sxaddpar, hdr0, 'SUBR_YLO', yrange[0], 'ASB 0-based subraster ylo'
          sxaddpar, hdr0, 'SUBR_YHI', yrange[1], 'ASB 0-based subraster yhi'
          sxaddpar, hdr0, 'LA_SIGF', sigfrac, 'la_cosmic sigfrac'
          sxaddpar, hdr0, 'LA_OBJL', objlim, 'la_cosmic objlim'
          sxaddpar, hdr0, 'LA_SIGC', sigclip, 'la_cosmic sigclip'
          sxaddpar, hdr0, 'LA_GAIN', gain, 'la_cosmic gain'
          sxaddpar, hdr0, 'LA_READN', readn, 'la_cosmic readn'
          sxaddpar, hdr0, 'LA_NITER', niter, 'la_cosmic niter'
          sxaddpar, hdr0, 'MEDRAD', medrad, 'ASB median-filter radius'
          sxaddpar, hdr0, 'NEGSIG', negsig, 'ASB negative-pixel thresh'
          sxaddpar, hdr0, 'ASB_SSUB', 'YES', 'Sky subtracted using MDRIZSKY'
          
          img = img[2*medrad:nx2-2*medrad-1,2*medrad:ny2-2*medrad-1]
          err = err[2*medrad:nx2-2*medrad-1,2*medrad:ny2-2*medrad-1]
          dqm = dqm[2*medrad:nx2-2*medrad-1,2*medrad:ny2-2*medrad-1]
          mfimg = mfimg[2*medrad:nx2-2*medrad-1,2*medrad:ny2-2*medrad-1]
          nmask = nmask[2*medrad:nx2-2*medrad-1,2*medrad:ny2-2*medrad-1]
          crmask = crmask[medrad:nx3-medrad-1,medrad:ny3-medrad-1]
          omask = nmask + 2B*crmask
          
          img = img - mdsky
          mfimg = mfimg - mdsky
          
          process->Increment, 0.25
          process->Setmessage, 'writing out/cleaning up...', /append
          mwrfits, 0, cutoutfile[j], hdr0, /create
          mwrfits, Float(img), cutoutfile[j], hdr1
          mwrfits, Float(err), cutoutfile[j]
          mwrfits, dqm, cutoutfile[j]
          mwrfits, omask, cutoutfile[j]
          mwrfits, Float(mfimg), cutoutfile[j]
          
          Spawn, 'rm -f ' + tempname
          Spawn, 'rm -f ' + crname
          Spawn, 'rm -f ' + froot[j] + '_temp-out.fits'
          
          cutflag = file_test(cutoutfile[j])
          if (not metastruc_changed) then if (metastruc[j].cutflag ne cutflag) then metastruc_changed = 1L
          metastruc[j].cutflag = cutflag
          
          
          metastruc[j].rectiflag = file_test(acsproc_datafile(metastruc[j], /rectify)) 
          metastruc[j].combflag = file_test(acsproc_datafile(metastruc[j], /combine))
          maskfile = acsproc_datafile(metastruc[j], /mask)
          if (file_test(maskfile)) then begin
            maskstruc = acsproc_get_maskstruc(maskfile)
            if maskstruc.has_photmod and file_test(acsproc_datafile(metastruc[j], /photmod)) then metastruc[j].maskflag = 4 $
            else if maskstruc.has_bizzle and file_test(acsproc_datafile(metastruc[j], /biz)) then metastruc[j].maskflag = 3 $
            else if maskstruc.has_modim then metastruc[j].maskflag = 2 $
            else metastruc[j].maskflag = 1
            print , metastruc[j].visit,metastruc[j].maskflag
            void = where(maskstruc.fmask ne 0,fmaskflag)
            void = where(maskstruc.jmask ne 0,jmaskflag)
            metastruc[j].fmaskflag = (fmaskflag gt 0)
            metastruc[j].jmaskflag = (jmaskflag gt 0)
            metastruc[j].lens = acsproc_get_lens_status(maskstruc.lens_status)
          endif else begin
            metastruc[j].maskflag = file_test(maskfile)
            metastruc[j].fmaskflag = 0L
            metastruc[j].jmaskflag = 0L
            metastruc[j].lens = acsproc_get_lens_status(0L)
            print , metastruc[j].visit," --> found no mask" 
          endelse


        endif ; else acsproc_database_update, entity, column, value, visit_number, exposure_number, message=message, confirm=confirm
  
        n+=1
  
      endif

    endfor
    
  endif
    
  if metastruc_changed then acsproc_metastruc_update    
  process->Destroy  
      
end
