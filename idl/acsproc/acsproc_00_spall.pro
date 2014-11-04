function acsproc_00_spall_specprimary, ra, dec, spec=spec

  common acsproc_spall, spall, search_radius
  
  wh = where(spall.specprimary and abs(spall.plug_ra - ra) le search_radius and abs(spall.plug_dec - dec) le search_radius, specprimary)
  
  if (specprimary eq 1) then spec = spall[wh]
        
  return, specprimary
  
end

;---------------------------------------------------------------------

function acsproc_00_spall_galfitstruc, empty=empty, create=create
  

  galfitstrucfile = acsproc_datafile('galfit_targets',/hst_root)
  
  if not keyword_set(empty) then begin
    galfitstruc = [{uniqname:'BELLSJ0151+0049',galfit_sersic_xcen:511.92,galfit_sersic_ycen:511.93,galfit_sersic_r_eff:18.64,galfit_sersic_n:3.93,galfit_sersic_q:0.60,galfit_sersic_pa:-70.45,galfit_sersic_chi2:1.079}, $
    {uniqname:'BELLSJ0747+4448',galfit_sersic_xcen:511.97,galfit_sersic_ycen:512.05,galfit_sersic_r_eff:36.29,galfit_sersic_n:4.60,galfit_sersic_q:0.62,galfit_sersic_pa:-39.29,galfit_sersic_chi2:1.165}, $
    {uniqname:'BELLSJ0801+4727',galfit_sersic_xcen:512.00,galfit_sersic_ycen:512.03,galfit_sersic_r_eff:12.34,galfit_sersic_n:4.27,galfit_sersic_q:0.98,galfit_sersic_pa:40.26,galfit_sersic_chi2:1.108}, $
    {uniqname:'BELLSJ0830+5116',galfit_sersic_xcen:512.06,galfit_sersic_ycen:512.13,galfit_sersic_r_eff:22.87,galfit_sersic_n:4.11,galfit_sersic_q:0.76,galfit_sersic_pa:-55.34,galfit_sersic_chi2:1.080}, $
    {uniqname:'BELLSJ1221+3806',galfit_sersic_xcen:512.02,galfit_sersic_ycen:512.08,galfit_sersic_r_eff:55.93,galfit_sersic_n:8.54,galfit_sersic_q:0.84,galfit_sersic_pa:70.65,galfit_sersic_chi2:0.546}, $
    {uniqname:'BELLSJ1352+3216',galfit_sersic_xcen:511.94,galfit_sersic_ycen:512.03,galfit_sersic_r_eff:31.36,galfit_sersic_n:6.83,galfit_sersic_q:0.94,galfit_sersic_pa:-74.41,galfit_sersic_chi2:1.083}, $
    {uniqname:'BELLSJ1452+3323',galfit_sersic_xcen:511.97,galfit_sersic_ycen:511.92,galfit_sersic_r_eff:39.37,galfit_sersic_n:6.48,galfit_sersic_q:0.79,galfit_sersic_pa:82.53,galfit_sersic_chi2:1.072}, $
    {uniqname:'BELLSJ1631+1854',galfit_sersic_xcen:511.94,galfit_sersic_ycen:511.93,galfit_sersic_r_eff:49.86,galfit_sersic_n:5.26,galfit_sersic_q:0.94,galfit_sersic_pa:-1.87,galfit_sersic_chi2:1.079}, $
    {uniqname:'BELLSJ1637+1439',galfit_sersic_xcen:511.90,galfit_sersic_ycen:512.03,galfit_sersic_r_eff:14.62,galfit_sersic_n:2.26,galfit_sersic_q:0.53,galfit_sersic_pa:-58.04,galfit_sersic_chi2:1.115}, $
    {uniqname:'BELLSJ2122+0409',galfit_sersic_xcen:511.96,galfit_sersic_ycen:511.92,galfit_sersic_r_eff:56.79,galfit_sersic_n:6.65,galfit_sersic_q:0.82,galfit_sersic_pa:-48.85,galfit_sersic_chi2:1.092}, $
    {uniqname:'BELLSJ2125+0411',galfit_sersic_xcen:511.96,galfit_sersic_ycen:511.89,galfit_sersic_r_eff:42.10,galfit_sersic_n:5.62,galfit_sersic_q:0.75,galfit_sersic_pa:-69.29,galfit_sersic_chi2:1.078}, $
    {uniqname:'BELLSJ1541+1812',galfit_sersic_xcen:511.90,galfit_sersic_ycen:511.96,galfit_sersic_r_eff:14.21,galfit_sersic_n:3.04,galfit_sersic_q:0.78,galfit_sersic_pa:-32.45,galfit_sersic_chi2:1.108}, $
    {uniqname:'BELLSJ2303+0037',galfit_sersic_xcen:511.96,galfit_sersic_ycen:511.90,galfit_sersic_r_eff:44.40,galfit_sersic_n:5.10,galfit_sersic_q:0.78,galfit_sersic_pa:86.24,galfit_sersic_chi2:1.093}]
  endif else begin
    galfitstruc = {galfit_devauc_xcen:0., galfit_devauc_ycen:0., galfit_devauc_r_eff:0., galfit_devauc_q:0., galfit_devauc_pa:0., galfit_devauc_chi2:0., galfit_sersic_xcen:0., galfit_sersic_ycen:0., galfit_sersic_r_eff:0., galfit_sersic_n:0., galfit_sersic_q:0., galfit_sersic_pa:0., galfit_sersic_chi2:0.}
  endelse
  
  if keyword_set(create) then mwrfits, galfitstruc, galfitstrucfile, /create, /silent
  
  return, galfitstruc
          
end

;---------------------------------------------------------------------

pro acsproc_00_spall, process, selected=selected, photmod=photmod, galfit=galfit

  common acsproc_data
  common acsproc_metastruc
  common acsproc_spall, spall, search_radius

  search_radius = 0.001
  ebv2ext = [5.155, 3.793, 2.751, 2.086, 1.479]
  ifilt = 3
  rv = 3.1
  e2e = (rv * ext_ccm(acs_filtwave[ifilt], rv))[0]
  
  void = where(metastruc.subexposure eq 1,nvisits)
  if not keyword_set(selected) then begin
    if keyword_set(photmod) then selected = where(metastruc.subexposure eq 1,nvis) $
    else selected = where(metastruc.subexposure eq 1,nvis)
  endif else begin
    selected -= 1
    nvis = n_elements(selected)
  endelse
  
  steps=nvis + 1
    
  process->Start, steps=steps
  process->Setmessage, 'Processing Proposal ID '+program+' ('+strtrim(nvis,2)+'/'+strtrim(nvisits,2)+' Visits)', /label
  process->Setmessage, ''  
  metastruc_changed = 0L
  
  
  boss_spectro_redux = getenv('BOSS_SPECTRO_REDUX')
  print , boss_spectro_redux
  if strpos(boss_spectro_redux,'/',strlen(boss_spectro_redux)-1) lt 0 then boss_spectro_redux+='/'
  
  ver = getenv('RUN2D')
  spallver = 'spAll-' + ver
  spallfile = boss_spectro_redux + spallver+'.fits'
    
  if (nvis gt 0) and file_test(spallfile) then begin
  
    metastruc_visit = metastruc[selected]
    
    process->Increment, 1.0
    
    spallstrucfile = acsproc_datafile(spallver+'_targets',/hst_root)
    
    message = 'SPALL: Reading '+spallfile
    process->setMessage, message, /APPEND
    
    spall = mrdfits(spallfile,1,/silent)
        
    spallstruc = []

    galfitstrucfile = acsproc_datafile('galfit_targets',/hst_root)
    if keyword_set(galfit) and file_test(galfitstrucfile) then begin
      message = 'SPALL: Reading '+galfitstrucfile
      process->setMessage,+message, /APPEND
      ;galfitstruc = acsproc_00_spall_galfitstruc(/create)
      galfitstruc = mrdfits(galfitstrucfile,1,/silent)
      galfitstruc_empty = acsproc_00_spall_galfitstruc(/empty)
    endif

    for q = 0,nvis-1 do begin
    
      if process->CheckCancel() then begin
        process->Step, steps
        process->Setmessage, 'Cancelling...', /append
        if metastruc_changed then acsproc_metastruc_update
        process->Destroy
        return
      endif else process->Increment, 1.0
                
      specprimary = acsproc_00_spall_specprimary(metastruc_visit[q].ra_targ, metastruc_visit[q].dec_targ, spec=sp)
      
      process->Increment, 1.0
      
      if specprimary eq 1 then begin
      
        message = 'SPALL: spectrum found at '+strtrim(sp.plate,2)+'-'+strtrim(sp.mjd,2)+'-'+strtrim(sp.fiberid,2)
        process->setMessage, "Target "+metastruc_visit[q].uniqname+" Visit "+metastruc_visit[q].visit+' -> '+message, /APPEND
        
        targ_data = acsproc_get_target_data(metastruc_visit[q].targname)
        
        this_spallstruc = { $
          plate:sp.plate, $
          mjd:sp.mjd, $
          fiberid:sp.fiberid, $
          proposid: metastruc_visit[q].proposid, $
          targname: metastruc_visit[q].targname, $
          uniqname: metastruc_visit[q].uniqname, $
          visit: metastruc_visit[q].visit, $
          ra_targ: metastruc_visit[q].ra_targ, $
          dec_targ: metastruc_visit[q].dec_targ, $
          equinox: metastruc_visit[q].equinox, $
          fileroot: metastruc_visit[q].uniqname + '_' + $
                         strtrim(string(metastruc_visit[q].proposid, format='(i)'),2) + '_' + $
                         strtrim(string(metastruc_visit[q].expstart, format='(i)'),2) + '_' + $
                         acs_filter[metastruc_visit[0].filtcode] + '_' +  strtrim(metastruc_visit[q].nexp,2), $
           date_obs: metastruc_visit[q].date_obs, $
           time_obs: metastruc_visit[q].time_obs, $
           expstart: metastruc_visit[q].expstart, $
           expend: metastruc_visit[q].expend, $
           exptime: metastruc_visit[q].exptime, $
           lens: strtrim(metastruc_visit[q].lens,2), $
           z_lens:targ_data.z_lens, $
           z_source:targ_data.z_source, $
           d_lens:targ_data.d_lens, $
           d_source:targ_data.d_source, $
           devmag:sp.devmag, $
           devmagerr:sp.devmagerr, $
           devflux:sp.devflux, $
           devflux_ivar:sp.devflux_ivar, $
           theta_dev:sp.theta_dev, $
           theta_deverr:sp.theta_deverr, $
           ab_dev:sp.ab_dev, $
           ab_deverr:sp.ab_deverr, $
           acs_filter:acs_filter[ifilt], $
           extin_corr:e2e * sp.extinction[ifilt] / ebv2ext[ifilt]}
           
         photmodfile = acsproc_datafile(metastruc_visit[q],/photmod)
         if keyword_set(photmod) and file_test(photmodfile) then begin
  
           message = 'SPALL: Reading '+photmodfile
           process->setMessage, "Target "+metastruc_visit[q].uniqname+" Visit "+metastruc_visit[q].visit+' -> '+message, /APPEND
           
           photmod = mrdfits(photmodfile,3,/silent)
           
           photmodstruc = {photflam:photmod.photflam, $
             photplam:photmod.photplam, $
             abmag_zpt:photmod.abmag_zpt, $
             abmag_dev:photmod.abmag_dev, $
             cps_dev:photmod.cps_dev, $
             cps_deverr:photmod.cps_deverr, $
             acs_rdev:photmod.acs_rdev, $
             acs_rdeverr:photmod.acs_rdeverr, $
             acs_qdev:photmod.acs_qdev, $
             acs_qdeverr:photmod.acs_qdeverr, $
             acs_padev:photmod.acs_padev, $
             acs_padeverr:photmod.acs_padeverr, $
             acs_nser:photmod.acs_nser, $
             acs_nsererr:photmod.acs_nsererr }
             
           this_spallstruc = struct_addtags(this_spallstruc, photmodstruc)
                      
         endif
                                              
           
         if not keyword_set(photmod) or file_test(photmodfile) then begin
         
           if keyword_set(galfit) and (n_elements(galfitstruc) gt 0) then begin
           
             wh = where(galfitstruc.uniqname eq metastruc_visit[q].uniqname, found)
             
             if found then begin
               message = 'SPALL: Adding galfit'
               process->setMessage, "Target "+metastruc_visit[q].uniqname+" Visit "+metastruc_visit[q].visit+' -> '+message, /APPEND
               this_galfitstruc = struct_trimtags(galfitstruc[wh],except_tags=['uniqname'])
               this_spallstruc = struct_addtags(this_spallstruc, this_galfitstruc)
             endif else begin
               message = 'SPALL: ERROR Cannot find galfit'
               process->setMessage, "Target "+metastruc_visit[q].uniqname+" Visit "+metastruc_visit[q].visit+' -> '+message, /APPEND
               this_spallstruc = struct_addtags(this_spallstruc, galfitstruc_empty)
             endelse
    
           endif
           
           spallstruc = [spallstruc,this_spallstruc]

         endif
         
                        
       endif else begin
          
         message = 'SPALL: specprimary='+strtrim(specprimary,2)+'. Cannot find spectrum found at ra='+strtrim(metastruc_visit[q].ra_targ,2)+'-'+strtrim(metastruc_visit[q].dec_targ,2)
         process->setMessage, "Target "+metastruc_visit[q].uniqname+" Visit "+metastruc_visit[q].visit+' -> '+message, /APPEND
              
       endelse
                          
     endfor
    
     if n_elements(spallstruc) gt 0 then mwrfits, spallstruc, spallstrucfile, /create, /silent
    
  endif else print, 'missing spAll file:',spallfile
      
  process->Destroy  
  
end

;---------------------------------------------------------------------

