; Copyright 2005-2014, Joel R. Brownstein, Adam S. Bolton, and Scott Burles
; All rights reserved.  Permission required for use.

common acsproc_state, state, data, offline, cycles, programs
common acsproc_data, hst_dataroot, dataroot, metafile, cycle, program, acs_filter, acs_filtwave, run_process, subexposure_count, hw, imsize, phw, dpix
common acsproc_metastruc, metastruc, header
common acsproc_maskstruc, maskstruc, mask_header, cornerpix
;---------------------------------------------------------------------

function acsproc_data_table, newdata=newdata

; function to return the contents of the acsproc table

  common acsproc_state
  common acsproc_data  
  common acsproc_metastruc
  
  ;data_table = {uniqname:'',visit:'',subexposure:0L,timestamp:'',status:'',lens:''}
  data_table = {uniqname:'',visit:'',timestamp:'',status:'',lens:''}
  
  if keyword_set(newdata) then metastruc=0 
  
  if size(metastruc, /tname) eq 'STRUCT' and state.limit_exposures then begin
    wh = where(metastruc.subexposure eq 1,nexp)
    if nexp gt 0 then exposure = metastruc[wh]
  endif else exposure = metastruc
  
  
  if size(exposure, /tname) eq 'STRUCT' then begin
    data_table = replicate(data_table,n_elements(exposure))
    for q=0,n_elements(exposure)-1 do begin
      data_table[q].uniqname = exposure[q].uniqname 
      data_table[q].visit = exposure[q].visit 
     ; data_table[q].subexposure = exposure[q].subexposure
      data_table[q].status = exposure[q].status
      data_table[q].timestamp = exposure[q].timestamp
      data_table[q].lens = exposure[q].lens
    endfor
  endif else data_table = replicate(data_table,256)
  
  return, data_table
   
 end

;---------------------------------------------------------------------

function acsproc_data_header, newdata=newdata

; function to return the contents of the acsproc table

  common acsproc_state
  common acsproc_data  
  common acsproc_metastruc
  
  blank_header = '                                                                                           '
  if keyword_set(newdata) then header=0 
    
    if size(header, /tname) eq 'STRING' then begin
      timestamp = sxpar(header, 'TIMESTMP')
      proposid = strtrim(sxpar(header, 'PROPOSID'),2)
    
      data_header = 'Phase II ID '+proposid
      if timestamp ne '' then begin
        data_header_right = 'UPDATED: '+timestamp 
        data_header = strmid(data_header + blank_header,0,strlen(blank_header)-strlen(data_header_right))+data_header_right
      endif else begin
        data_header = strmid(data_header + blank_header,0,strlen(blank_header))
      endelse
      
   endif else data_header = blank_header
      
   return, data_header
   
 end
;---------------------------------------------------------------------

function acsproc_status_table, status_table=status_table

; function to return a blank status table

  common acsproc_state
  common acsproc_data  
  common acsproc_metastruc

; function to return a blank status table
  if not keyword_set(status_table) then begin
    status_data = {statuscode:'',program_count:'',run_process:''} 
    status_table = replicate(status_data,n_elements(run_process))
  
    for r=0,n_elements(run_process)-1 do begin
      status_table[r].statuscode = run_process[r].statuscode
      status_table[r].run_process= run_process[r].title
    endfor
  
  endif
  
  ; function to return a status table
  if size(metastruc, /tname) eq 'STRUCT' and keyword_set(status_table) then begin
    nexp = metastruc[0].nexp
    nmetastruc = strtrim(n_elements(metastruc)/nexp,2)
    nvis = strtrim(subexposure_count,2)
    for r=0,n_elements(run_process)-1 do begin
      if r gt 0 then begin
        if r eq 7 then void = where(metastruc.statuscode eq r and metastruc.subexposure eq 1 and not (metastruc.fmaskflag ne 0 and metastruc.jmaskflag ne 0),c) $
        else void = where(metastruc.statuscode eq r and metastruc.subexposure eq 1,c)
        c = strtrim(c,2)
        if r lt n_elements(run_process)-1 then begin
          if c eq '0' then c=''
          status_table[r].program_count = c
        endif else status_table[r].program_count = c+"/"+nmetastruc
       endif else status_table[r].program_count = nmetastruc+"/"+nvis
    endfor
  endif
   
  return, status_table
   
 end
;---------------------------------------------------------------------

function acsproc_status_header

; function to return the contents of the acsproc table

  common acsproc_state
  common acsproc_data  
  
    blank_header = '                                                                                           '
    
      timestamp = systime()
    
      status_header = 'Status Report'
      if timestamp ne '' then begin
        status_header_right = 'UPDATED: '+timestamp 
        status_header = strmid(status_header + blank_header,0,strlen(blank_header)-strlen(status_header_right))+status_header_right
      endif else begin
        status_header = strmid(status_header + blank_header,0,strlen(blank_header))
      endelse
            
   return, status_header
   
 end
 
;--------------------------------------------------------------------

function acsproc_datafile, subexposure, froot=froot, cutout=cutout, rectify=rectify, combine=combine, mask=mask, photmod=photmod, biz=biz, old=old, test=test, hst_root=hst_root, ext=ext, subext=subext

  common acsproc_data
  
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
  endif else begin
    datafile = dataroot + subexposure.visit + '/'+ subexposure.rootname + ftype
    if keyword_set(cutout) then datafile +=  '_cutout' + subext + extension $
    else if keyword_set(rectify) then datafile += '_recti' + subext + extension $
    else if not keyword_set(froot) then datafile +=  subext + extension
  endelse
    
  return, datafile              
  
end

;--------------------------------------------------------------------

function acsproc_filtcode, subexposure

  common acsproc_data

  filtcode = 0L
  for i=0,n_elements(acs_filter)-1 do begin
    if subexposure.filter1 eq acs_filter[i] or subexposure.filter2 eq acs_filter[i] then begin
      filtcode = i
      break
    endif
  endfor
  
  return, i
  
end

;--------------------------------------------------------------------

function get_num_subdir, root, ascending=ascending, descending=descending

  subdir=0
  if strpos(root,'/',strlen(root)-1) lt 0 then root+='/'
  subdir_search = file_search(root+'*')
  for i=0,n_elements(subdir_search)-1 do begin
    subdir_test_pos0 = strpos(subdir_search[i],'/',/reverse_search)+1
    subdir_test_pos1 = strlen(subdir_search[i])
    subdir_test = strmid(subdir_search[i],subdir_test_pos0,subdir_test_pos1-subdir_test_pos0)
    if valid_num(subdir_test) then begin
        if keyword_set(subdir) then subdir=[subdir,subdir_test] $
        else subdir=[subdir_test]
    endif
  endfor
  
  if n_elements(subdir) gt 1 then begin
    if keyword_set(ascending) then subdir = subdir[sort(subdir)]  else if keyword_set(descending) then subdir = subdir[reverse(sort(subdir))]
  endif


  return, subdir
  
end
  
;---------------------------------------------------------------------
function acsproc_database_webget, url, POST=post
  ;==============================================================================
  ; wrapper for webget that catches errors in case user is offline/site is down.
  ;==============================================================================
  servers = n_elements(url)
  server = 0
  catch, webget_error
  if (webget_error ne 0) then begin
    print, "unable to connect to host"
    server = server+1
    if (server eq servers) then begin
      catch,/cancel
      response =  {Text:'NULL'}
    endif else response = webget(url[server],POST=post,/silent)
  endif
  if (server lt servers) then begin
    response = webget(url[server],POST=post)
    if (response.Text[0] eq '' and n_elements(response.Text) eq 2) then response.Text[0] = response.Text[1]
    if (response.Text[0] eq '') then begin
      if (server lt servers-1) then begin
        print, "unable to get response from host (1)"
        response = webget(url[server+1],POST=post,/silent)
      endif
    endif
  endif else return, ''
  return, response.Text[0]
end
;---------------------------------------------------------------------
pro acsproc_subexposure_count
  ;==============================================================================
  ; Database Function: communicate with database URL by HTTP POST
  ;==============================================================================
  common acsproc_state
  common acsproc_data
  
  if not offline then begin
    ;url = ['http://slacs.astro.utah.edu/','http://slacs.astro.utah.edu/']
    url = ['http://slacs.astro.utah.edu/']
    cmdurl =  url + 'mysql.php'
    post = {func:'subexposure_count',mysql:'zwicky',ver:'101120101139',program_number:program}
    subexposure_count = acsproc_database_webget(cmdurl,POST=post)
  endif else subexposure_count=''
end
;---------------------------------------------------------------------
function acsproc_get_target_data, target_name
  ;==============================================================================
  ; Database Function: communicate with database URL by HTTP POST
  ;==============================================================================
  common acsproc_data
  post = {func:'target',mysql:'zwicky',ver:'101120101139',action:'get_data',program_number:program,target_name:repstr(target_name, '+', '<PLUS/>')}
  item = {z_lens:0.0,z_source:0.0,d_lens:0.0,d_source:0.0,confirmed:0L}
  acsproc_database_select, post, '', item, select
  if (n_elements(select) eq 1) then begin
    if select[0].confirmed eq 1L then return, select[0] else return, item
  endif else return, item
end
;---------------------------------------------------------------------
pro acsproc_program_timestamp
  ;==============================================================================
  ; Database Function: communicate with database URL by HTTP POST
  ;==============================================================================
  common acsproc_state
  common acsproc_data
  
  if not offline then begin
    post = {func:'program',mysql:'zwicky',ver:'101120101139',action:'update_modified',program_number:program}
    item = {modified:'',confirmed:0L}
    acsproc_database_select, post, '', item, select
    if (n_elements(select) eq 1) then begin
      if select[0].confirmed eq 1L then print, 'PROGRAM TIMESTAMP: '+strtrim(select[0].modified,2) else print, 'ERROR - WILL NOT UPDATE PROGRAM'
    endif else print, 'ERROR - CANNOT UPDATE PROGRAM'
  endif else print, 'OFFLINE - CANNOT UPDATE PROGRAM'
end
;---------------------------------------------------------------------
pro acsproc_set_target_lens, target_name, lens
  ;==============================================================================
  ; Database Function: communicate with database URL by HTTP POST
  ;==============================================================================
  common acsproc_state
  common acsproc_data
  
  if not offline then begin
    post = {func:'target_lens',mysql:'zwicky',ver:'101120101139',program_number:program,target_name:repstr(target_name, '+', '<PLUS/>'),lens:lens}
    item = {confirmed:0L}
    acsproc_database_select, post, '', item, select
    if (n_elements(select) eq 1) then begin
      update_lens='???'
      if (lens eq 0) then update_lens = 'NEITHER YES/NO/MAYBE?' $
      else if (lens eq 3) then update_lens = 'NO' $
      else if (lens eq 2) then update_lens = 'YES = DETECTED WITH COUNTER-IMAGE (MODELABLE)' $
      else if (lens eq 1) then update_lens = 'MAYBE = DETECTED WITHOUT COUNTER-IMAGE (UNMODELABLE)' 
      if select[0].confirmed eq 1L then print, 'UPDATE LENS '+update_lens else print, 'ERROR - WILL NOT UPDATE LENS'
    endif else begin
	print, 'ERROR - CANNOT UPDATE LENS'
    endelse
  endif else print, 'OFFLINE - CANNOT UPDATE LENS'
end
;---------------------------------------------------------------------
pro acsproc_database_post, cmd, post, sid, response=response
  ;==============================================================================
  ; Database Function: communicate with database URL by HTTP POST
  ;==============================================================================
  url = ['http://slacs.astro.utah.edu/','http://slacs.astro.utah.edu/']
  if (sid ne '') then cmdurl = url + cmd + '.php?PHPSESSID='+sid else cmdurl = url + cmd + '.php'
  response = acsproc_database_webget(cmdurl,POST=post)
end

;---------------------------------------------------------------------

pro acsproc_database_select, post, sid, item, select
  ;==============================================================================
  ; Database Function: parse response from database on selected item
  ;==============================================================================
  cmd = 'mysql'
  acsproc_database_post, cmd, post, sid, response=response
  responseText=''
  if (n_elements(response) eq 1) then responseText = response
  if (responseText ne 'NULL') and (responseText ne '') then begin
    escape = '`'
    responses = STRSPLIT(responseText, ESCAPE=escape, /EXTRACT,';')
    nresponses = N_ELEMENTS(responses)
    select = replicate(item,nresponses)
    if (nresponses gt 0) then begin
      select = replicate(item,nresponses)
      for i = 0,nresponses-1 do begin
        pairs = STRSPLIT(responses[i], ESCAPE=escape, /EXTRACT,',')
        npairs = N_ELEMENTS(pairs)
        itemtag = TAG_NAMES(item)
        if (npairs eq N_ELEMENTS(itemtag)) then begin
          for j=0,npairs-1 do begin
            pair = STRSPLIT(pairs[j], ESCAPE=escape, /EXTRACT,'=')
            if (N_ELEMENTS(pair) eq 2) then begin
              key = pair[0]
              value = pair[1]
              if (key = itemtag[j]) then select[i].(j)=value
            endif
          endfor
        endif else begin
          uumessage = "bad response from database, please contact administrator." & print, uumessage
        endelse
      endfor
    endif else begin
      uumessage = "null response from database, please contact administrator." & print, uumessage
    endelse
  endif else begin
  endelse
end

;---------------------------------------------------------------------

pro acsproc_database_member, action
  ;==============================================================================
  ; Database Function: authenticate member (login) and get recent comments
  ;==============================================================================
  common uuEmlineBase_state, uuState, uuData, uuComment, recentcommentlist
  post = {func:'member',mysql:'emline',ver:'053120101324',siteID:'1',action:action,username:uuState.username,password:uuState.password,remote:uuState.remote}
  item = {loggedin:'',sid:'',memberid:0L}
  acsproc_database_select, post, uuState.sid, item, select
  if (n_elements(select) eq 1) then begin
    uuState.loggedin=fix(select[0].loggedin)
    uuState.sid=select[0].sid
    uuState.memberid=select[0].memberid
  endif else begin
    uuState.loggedin=0
    uuState.sid=''
    uuState.memberid=0
  endelse
  if (uuState.loggedin) then acsproc_database_select_recentcommentlist
end

;---------------------------------------------------------------------

pro acsproc_database_update_status, program_number, visit_number, exposure_number, status, message=message, confirm=confirm
  ;==============================================================================
  ; Database Function: insert exposure redux_name and receive message response from database
  ;==============================================================================
  common acsproc_data
    
  redux_name = StrJoin( StrSplit(redux_name, '+', /Extract), '%2b')
  post = {func:'status',mysql:'zwicky',ver:'101120101139',action:'update',program_number:program_number,visit_number:visit_number,exposure_number:exposure_number,status:status}
  item = {message:'',confirm:''}
  acsproc_database_select, post, '', item, select
  confirm = 0
  if (n_elements(select) eq 1) then begin
    message=select[0].message
    confirm=select[0].confirm
  endif else begin
    message = 'Database Error > '
    confirm = 0L
  endelse
end

;---------------------------------------------------------------------

pro acsproc_database_update, entity, column, value, visit_number, exposure_number, message=message, confirm=confirm
  ;==============================================================================
  ; Database Function: insert exposure redux_name and receive message response from database
  ;==============================================================================
  common acsproc_data
    
  func = entity+'_'+column
  redux_name = StrJoin( StrSplit(value, '+', /Extract), '%2b')
  post = {func:func,mysql:'zwicky',ver:'101120101139',action:'update',program_number:program,visit_number:visit_number,exposure_number:exposure_number,redux_name:redux_name}
  item = {message:'',confirm:''}
  acsproc_database_select, post, '', item, select
  confirm = 0
  if (n_elements(select) eq 1) then begin
    message=select[0].message
    confirm=select[0].confirm
  endif else begin
    message = 'Database Error > '
    confirm = 0L
  endelse
end
;---------------------------------------------------------------------

pro acsproc_database_update_table, table, program_number, visit_number, filename, message=message, confirm=confirm
  ;==============================================================================
  ; Database Function: insert exposure redux_name and receive message response from database
  ;==============================================================================

  cmd = table+'_update'
  arg = " -p "+strtrim(program_number,2) + " -v "+visit_number + " -f "+filename
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

pro acsproc_set_cycle_program, refresh=refresh

  common acsproc_state
  common acsproc_data
  common acsproc_metastruc

   dataroot = hst_dataroot+cycle+'/'+program+'/'
   metafile = dataroot+'.metafile.fits'
   acsproc_subexposure_count
      
   status_table = acsproc_status_table()
   status_ncols = n_elements(status_table[0])

   newdata = 1 - file_test(metafile) 
   
   if newdata then begin
     data_table = acsproc_data_table(/newdata)
     data_header = acsproc_data_header(/newdata)
   endif else begin
     metastruc = mrdfits(metafile, 1, header, /silent)
     data_table = acsproc_data_table()
     data_header = acsproc_data_header()
     status_table = acsproc_status_table(status_table=status_table)
   endelse   
   
   status_header = acsproc_status_header()
   
   if keyword_set(refresh) then begin
     widget_control, data.header_id, set_value=data_header
     widget_control, data.table_id, set_value=data_table
     widget_control, state.header_id, set_value=status_header
     widget_control, state.table_id, set_value=status_table
   endif else begin
     table_nrows = 255 ; the maximum number of subexposures in any program (for format consistency)
     table_format = make_array(5,table_nrows,/string,value='')
     ;for r=1,table_nrows-1 do table_format[2,r] = '(i)'
     blue_textcolor = [0,0,222]
     blue_cellcolor = [210,210,222]
     red_textcolor = [255,25,25]
     white_cellcolor = [222,222,255]
     yellow_cellgrad  = [-4,-4,-3]
     data_table_col = data_table[0]
     table_ncols = n_tags(data_table_col)
     new_data_table = replicate(data_table_col,table_nrows) ; give lots of room 
     table_textcolor = make_array(3,table_ncols,table_nrows,/integer,value=0)
     table_cellcolor = make_array(3,table_ncols,table_nrows,/integer,value=222)
     for j=0,table_nrows-1 do begin
       table_textcolor[*,0,j] = blue_textcolor[*]
       table_textcolor[*,table_ncols-1,j] = blue_textcolor[*]
     endfor
     for j=0,table_nrows-1,2 do begin
      for i = 0,2 do table_cellcolor[i,0:table_ncols-1,j] = blue_cellcolor[i]
     endfor
     status_ncols = n_tags(status_table[0])
     status_nrows = n_elements(status_table)
     status_format = make_array(3,status_nrows,/string,value='')
     status_textcolor = make_array(3,status_ncols,status_nrows,/integer,value=0)
     status_cellcolor = make_array(3,status_ncols,status_nrows,/integer,value=222)
     for j=0,status_nrows-1 do begin
       status_textcolor[*,0,j] = blue_textcolor[*]
       status_textcolor[*,status_ncols-1,j] = red_textcolor[*]
       status_cellcolor[*,status_ncols-1,j] = white_cellcolor[*]
     endfor
     for j=1,status_nrows-1 do begin
      for i = 0,2 do begin
        status_cellcolor[i,0:status_ncols-1,j] = white_cellcolor[i]+j*yellow_cellgrad[i]
      endfor
     endfor
     status_textcolor[*,1,0] = blue_textcolor[*]
     status_textcolor[*,1,status_nrows-1] = blue_textcolor[*]
     
     table_widths = make_array(5,/integer,value=150)
     table_widths[1] = 50
     table_widths[4] = 50
     status_table_widths = make_array(3,/integer,value=186)

     data.header_id=widget_label(state.data_base_id,value=data_header,frame=1)
     data.table_id=widget_table(state.data_base_id,value=new_data_table,/row_major,/no_row_headers,y_scroll_size=26,column_labels=['Name','Visit','Reduction Date','Process Waiting','Lens'],alignment=1,format=table_format,column_widths=table_widths,uvalue='data_table',background_color=table_cellcolor,foreground_color=table_textcolor,/all_events)
     state.header_id=widget_label(state.data_base_id,value=status_header,frame=1)
     state.table_id=widget_table(state.data_base_id,value=status_table,/row_major,/no_row_headers,y_scroll_size=n_elements(status_table),column_labels=['Process', 'Number of Visits', 'Run Process NOW'],alignment=1,format=status_format,column_widths=status_table_widths,uvalue='run_process',background_color=status_cellcolor,foreground_color=status_textcolor,/all_events)
   endelse
   
end

;---------------------------------------------------------------------

function acsproc_load_maskstruc, fits=fits

  common acsproc_data
    
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

pro acsproc_maskstruc_write, maskfile

  ; flatten the maskstruc so it can be saved as a fits file.
  
  common acsproc_maskstruc
  
  flat_maskstruc = {modimID:maskstruc.modimID, xc:maskstruc.xc, yc:maskstruc.yc, nthetaID:maskstruc.nthetaID, lens_status:maskstruc.lens_status, has_modim:maskstruc.has_modim, has_bizzle:maskstruc.has_bizzle, has_photmod:maskstruc.has_photmod, fmask:maskstruc.fmask, jmask:maskstruc.jmask, $
    moffat_modimID:maskstruc.moffat.modimID, moffat_name:maskstruc.moffat.name, moffat_modim:maskstruc.moffat.modim, moffat_pars:maskstruc.moffat.pars, moffat_pars_err:maskstruc.moffat.pars_err, moffat_chi2:maskstruc.moffat.chi2, $
    bsplinecen_modimID:maskstruc.bsplinecen.modimID, bsplinecen_name:maskstruc.bsplinecen.name, bsplinecen_modim:maskstruc.bsplinecen.modim, bsplinecen_dmodim:maskstruc.bsplinecen.dmodim, bsplinecen_pars:maskstruc.bsplinecen.pars, bsplinecen_pars_err:maskstruc.bsplinecen.pars_err, bsplinecen_chi2:maskstruc.bsplinecen.chi2, bsplinecen_nthetaID:maskstruc.bsplinecen.nthetaID, bsplinecen_rbkptID:maskstruc.bsplinecen.rbkptID, $
    bsplinefit_02_modimID:maskstruc.bsplinefit[1].modimID, bsplinefit_02_name:maskstruc.bsplinefit[1].name, bsplinefit_02_modim:maskstruc.bsplinefit[1].modim, bsplinefit_02_chi2:maskstruc.bsplinefit[1].chi2, bsplinefit_02_nthetaID:maskstruc.bsplinefit[1].nthetaID, bsplinefit_02_rbkptID:maskstruc.bsplinefit[1].rbkptID, $
    bsplinefit_04_modimID:maskstruc.bsplinefit[2].modimID, bsplinefit_04_name:maskstruc.bsplinefit[2].name, bsplinefit_04_modim:maskstruc.bsplinefit[2].modim, bsplinefit_04_chi2:maskstruc.bsplinefit[2].chi2, bsplinefit_04_nthetaID:maskstruc.bsplinefit[2].nthetaID, bsplinefit_04_rbkptID:maskstruc.bsplinefit[2].rbkptID, $
    bsplinefit_024_modimID:maskstruc.bsplinefit[3].modimID, bsplinefit_024_name:maskstruc.bsplinefit[3].name, bsplinefit_024_modim:maskstruc.bsplinefit[3].modim, bsplinefit_024_chi2:maskstruc.bsplinefit[3].chi2, bsplinefit_024_nthetaID:maskstruc.bsplinefit[3].nthetaID, bsplinefit_024_rbkptID:maskstruc.bsplinefit[3].rbkptID, $
    bsplinefit_01_modimID:maskstruc.bsplinefit[4].modimID, bsplinefit_01_name:maskstruc.bsplinefit[4].name, bsplinefit_01_modim:maskstruc.bsplinefit[4].modim, bsplinefit_01_chi2:maskstruc.bsplinefit[4].chi2, bsplinefit_01_nthetaID:maskstruc.bsplinefit[4].nthetaID, bsplinefit_01_rbkptID:maskstruc.bsplinefit[4].rbkptID, $
    bsplinefit_012_modimID:maskstruc.bsplinefit[5].modimID, bsplinefit_012_name:maskstruc.bsplinefit[5].name, bsplinefit_012_modim:maskstruc.bsplinefit[5].modim, bsplinefit_012_chi2:maskstruc.bsplinefit[5].chi2, bsplinefit_012_nthetaID:maskstruc.bsplinefit[5].nthetaID, bsplinefit_012_rbkptID:maskstruc.bsplinefit[5].rbkptID, $
    bsplinefit_014_modimID:maskstruc.bsplinefit[6].modimID, bsplinefit_014_name:maskstruc.bsplinefit[6].name, bsplinefit_014_modim:maskstruc.bsplinefit[6].modim, bsplinefit_014_chi2:maskstruc.bsplinefit[6].chi2, bsplinefit_014_nthetaID:maskstruc.bsplinefit[6].nthetaID, bsplinefit_014_rbkptID:maskstruc.bsplinefit[6].rbkptID, $
    bsplinefit_0124_modimID:maskstruc.bsplinefit[7].modimID, bsplinefit_0124_name:maskstruc.bsplinefit[7].name, bsplinefit_0124_modim:maskstruc.bsplinefit[7].modim, bsplinefit_0124_chi2:maskstruc.bsplinefit[7].chi2, bsplinefit_0124_nthetaID:maskstruc.bsplinefit[7].nthetaID, bsplinefit_0124_rbkptID:maskstruc.bsplinefit[7].rbkptID}
  
  splog, "CREATE: "+maskfile 
  mwrfits, flat_maskstruc, maskfile, mask_header, /create, /silent

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

function acsproc_rbkpt, modim

  common acsproc_data
  
  case modim.rbkptID of
  
    0:  rbkpt= dpix * [1.D, 2.D, 4.D, 6.D, 8.D, 12.D, 16.D, 22.D, 30.D, 50.D, 60.D, 100.D, 150.D, 200.D]
    else: begin
      rbkpt = dpix * [1.D, 2.D, 4.D, 6.D, 8.D, 12.D, 16.D, 22.D, 30.D, 50.D, 60.D, 100.D, 150.D, 200.D]
      print, 'unknown rbkpt'
    end
    
  endcase
  
  return, rbkpt
  
 end
;---------------------------------------------------------------------

function acsproc_ntheta, modim, str=str
    
  if keyword_set(str) then begin
  
    case modim.nthetaID of
    
      0: str = '0'
      1: str = '2'
      2: str = '04'
      3: str = '024'
      4: str = '01
      5: str = '012'
      6: str = '014'
      7: str = '0124'
      else: begin
        str = '0'
        print, 'unknown nthetaID'
      end
      
    endcase

  endif
  
  case modim.nthetaID of
  
    0: ntheta = [0]
    1: ntheta = [0,-2,2]
    2: ntheta = [0,-4,4]
    3: ntheta = [0,-2,2,-4,4]
    4: ntheta = [0,-1,1]
    5: ntheta = [0,-1,1,-2,2]
    6: ntheta = [0,-1,1,-4,4]
    7: ntheta = [0,-1,1,-2,2,-4,4]
    else: begin
      ntheta = [0]
      print, 'unknown nthetaID'
    end
    
  endcase
 
  
  return, ntheta
  
 end
 
;---------------------------------------------------------------------

function acsproc_modim, modim_name, fit_nthetaID=fit_nthetaID

  common acsproc_maskstruc
  
  if n_elements(modim_name) gt 0 then begin
    case modim_name of
    
      'moffat': maskstruc.modimID = 0
      'bsplinecen': maskstruc.modimID = 1
      'bsplinefit': maskstruc.modimID = 2
      else: begin
        maskstruc.modimID = 0
        print, 'unknown modim name: ',modim_name
      end
      
    endcase
  endif
    
  case maskstruc.modimID of
  
    0: modim = maskstruc.moffat
    1: modim = maskstruc.bsplinecen
    2: begin
      if keyword_set(fit_nthetaID) then modim=maskstruc.bsplinefit[maskstruc.nthetaID] else modim=maskstruc.bsplinefit
    end
    else: begin
      modim = maskstruc.moffat
      print, 'unknown modimID: ',maskstruc.modimID
    end
    
  endcase
  
  return, modim
  
 end

;---------------------------------------------------------------------

pro acsproc_startup, cycle=this_cycle,program=this_program,offline=this_offline

; This routine initializes the acsproc internal variables, and creates and
; realizes the window widgets.  It is only called by the acsproc main
; program level, when there is no previously existing acsproc window.

; Initialize the common blocks

  common acsproc_state
  common acsproc_data

 
  hst_dataroot=getenv('HST_DATAROOT')
  if hst_dataroot eq '' then begin
    ok = DIALOG_MESSAGE('Please restart IDL with your HST_DATAROOT environment variable set',/INFORMATION)
    return
  endif else if strpos(hst_dataroot,'/',strlen(hst_dataroot)-1) lt 0 then hst_dataroot+='/'
  
  hst_psfdir=getenv('HST_PSFDIR')
  if hst_psfdir eq '' then begin
    ok = DIALOG_MESSAGE('Please restart IDL with your HST_PSFDIR environment variable set',/INFORMATION)
    return
  endif else if strpos(hst_psfdir,'/',strlen(hst_psfdir)-1) lt 0 then hst_psfdir+='/'  
  
  offline = keyword_set(this_offline)

  state = {base_id: 0L, $                 ; id of top-level base
           header_id: 0L, $        ; widget id for Status Table Header
           table_id: 0L, $                ; widget id for Status Table
           data_base_id: 0L, $            ; base id for Data Table
           limit_data: 0L, $              ; value of limit data widget (all = 0, new = 1)
           limit_exposures: 1L, $         ; value of exposures limit widget (all = 0, primary = 1)
           cycle_droplist_id:0L, $        ; id of cycle droplist widget
           cycle_droplist_index:0L, $     ; current index of cycle droplist widget
           program_droplist_id:0L, $      ; id of program droplist widget
           program_droplist_index:0L, $   ; current index of program droplist widget
           version: 1.0}                 ; acsproc version number
           
           
  data = {header_id: 0L, $
           table_id: 0L}
           
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
  
  run_process = [{title:'REFRESH', statuscode:'Downloaded/Archived'}, $
                {title:'CUTOUT', statuscode:'CUTOUT'}, $
                {title:'LCOORDS', statuscode:'LCOORDS'}, $
                {title:'REGISTER', statuscode:'REGISTER'}, $
                {title:'RECTIFY', statuscode:'RECTIFY'}, $
                {title:'COMBINE', statuscode:'COMBINE'}, $
                {title:'INITIAL MASK', statuscode:'INITIAL MASK'}, $
                {title:'FEATURE/JUNK MASK', statuscode:'FEATURE/JUNK MASK'}, $
                {title:'BSPLINE', statuscode:'BSPLINE'}, $
                {title:'BIZZLE', statuscode:'BIZZLE'}, $
                {title:'PHOTMOD', statuscode:'PHOTMOD'}, $
                {title:'RESET ALL', statuscode:'REDUCED'}]

  cycles = get_num_subdir(hst_dataroot,/descending)
  
  if (keyword_set(this_cycle)) then begin
    print, 'new cycle=',this_cycle
    wh = where(strtrim(this_cycle[0],2) eq cycles,found)
    if (found gt 0) then state.cycle_droplist_index = wh[0] else state.cycle_droplist_index=0
  endif else state.cycle_droplist_index=0
  cycle = cycles[state.cycle_droplist_index]
  
  programs = get_num_subdir(hst_dataroot+cycle,/ascending)

  if (keyword_set(this_program)) then begin
    wh = where(strtrim(this_program[0],2) eq programs,found)
    if (found gt 0) then state.program_droplist_index = wh[0] else state.program_droplist_index=0
  endif else state.program_droplist_index=0
  program = programs[state.program_droplist_index]
       
  base = widget_base(title = 'acsproc > Cycle:'+cycle+' Program: '+program, $
                   /column, /base_align_right, $
                   uvalue = 'acsproc_base', $
                   /tlb_size_events)
                   
   control_base = Widget_Base(base,/row,/base_align_right) 
   state.data_base_id = Widget_Base(base,/column,/base_align_left) 
   option_base = Widget_Base(base,/row,/nonexclusive,/base_align_right) 
   button_base = Widget_Base(base,/row,/base_align_right) 
         
   state.cycle_droplist_id = widget_droplist(control_base, title = 'Cycle:', uvalue = 'cycle', value = cycles)
   widget_control, state.cycle_droplist_id, set_droplist_select=state.cycle_droplist_index
   state.program_droplist_id = widget_droplist(control_base, title = 'Program:', uvalue = 'program', value = programs)
   widget_control, state.program_droplist_id, set_droplist_select=state.program_droplist_index
   void = widget_label(control_base, value = ' ')     

   table_widths = make_array(5,/integer,value=150)
   table_widths[1:2] = 50
   status_table_widths = make_array(3,/integer,value=186)
      
   status_table = acsproc_status_table()
   status_ncols = n_elements(status_table[0])

    acsproc_set_cycle_program
      
   export_button = widget_button(button_base, value = 'Export FITS',  uvalue = 'export')     
   ;spall_button = widget_button(button_base, value = 'Gen SpAll',  uvalue = 'spall')     
   timestamp_button = widget_button(button_base, value = 'Timestamp',  uvalue = 'program_timestamp')     
   done_button = widget_button(button_base, value = 'Done',  uvalue = 'done')     
   void = widget_label(button_base, value = '  ')     
   
   widget_control, base, /realize
    
   state.base_id = base
   
   acsproc_get_metadata, /newdata
end

;---------------------------------------------------------------------

pro acsproc_shutdown, windowid

; routine to kill the acsproc window(s) and clear variables to conserve
; memory when quitting acsproc.  The windowid parameter is used when
; acsproc_shutdown is called automatically by the xmanager, if acsproc is
; killed by the window manager.

  common acsproc_state

  ; Kill top-level base if it still exists
  if (xregistered ('acsproc')) then widget_control, state.base_id, /destroy
  
  ; Kill ATV 
  acsproc_atv_shutdown
  
  delvarx, state
  
  return    
end

;---------------------------------------------------------------------

pro acsproc_run_process, run_process_id, uniqname=uniqname

; Main event loop for acsproc top-level base, and for all the buttons.

  common acsproc_state
  common acsproc_data
  common acsproc_metastruc
  
  if keyword_set(uniqname) then begin
    wh = where(metastruc.uniqname eq uniqname and metastruc.subexposure eq 1, nselected)
    if nselected gt 0 then selected = wh else selected = 0L
  endif else selected = 0L
  
  case run_process_id of 
    0: acsproc_get_metadata, /newdata
    1: acsproc_cutout, process=Obj_New("acsproc_process", state.base_id, Title=title, Message_xsize='120', Message_ysize='10', /CancelButton), selected=selected
    2: acsproc_atv_imagemode, /cutout, /lcoords, /newmode
    3: acsproc_register, process=Obj_New("acsproc_process", state.base_id, Title=title, Message_xsize='120', Message_ysize='10', /CancelButton), selected=selected
    4: acsproc_rectify, process=Obj_New("acsproc_process", state.base_id, Title=title, Message_xsize='120', Message_ysize='10', /CancelButton), selected=selected
    5: acsproc_combine, process=Obj_New("acsproc_process", state.base_id, Title=title, Message_xsize='120', Message_ysize='10', /CancelButton), selected=selected
    6: acsproc_mask, process=Obj_New("acsproc_process", state.base_id, Title=title, Message_xsize='120', Message_ysize='10', /CancelButton), selected=selected
    7: begin
       acsproc_atv_imagemode, /combine, /newmode
       acsproc_atv_maskmode, 3, /newmode
       acsproc_atv_mode, 0, /newmode
    end
    8: acsproc_bspline, process=Obj_New("acsproc_process", state.base_id, Title=title, Message_xsize='120', Message_ysize='10', /CancelButton), selected=selected
    9: acsproc_bizzle, process=Obj_New("acsproc_process", state.base_id, Title=title, Message_xsize='120', Message_ysize='10', /CancelButton), selected=selected
    10: acsproc_photmod, process=Obj_New("acsproc_process", state.base_id, Title=title, Message_xsize='120', Message_ysize='10', /CancelButton), selected=selected
    11: acsproc_get_metadata, /newdata, /reset
    else:  print, 'No match for uvalue....' 
  endcase
  
end

;---------------------------------------------------------------------

pro acsproc_event, event

; Main event loop for acsproc top-level base, and for all the buttons.

  common acsproc_state
  common acsproc_data
  
  
  widget_control, event.id, get_uvalue = uvalue
    
  case uvalue of
  
    'acsproc_base': begin
      c = where(tag_names(event) EQ 'ENTER', count)
      if (count EQ 0) then widget_control, state.base_id, /clear_events
    end
      
    'run_process':  begin
      run_process_id = -1L
      if event.type eq 4 and (event.SEL_TOP ge 0) then begin
        if (event.SEL_LEFT eq 2 and event.SEL_RIGHT eq 2 ) and (event.SEL_TOP eq event.SEL_BOTTOM) then begin
          run_process_id = event.SEL_TOP
          if run_process_id ge 0 and run_process_id lt n_elements(run_process) then title = run_process[run_process_id].title else title=''
          acsproc_run_process, run_process_id
        endif
      endif
    end
    
    'data_table':  begin
      
      if event.type eq 4 and (event.SEL_TOP ge 0) then begin
        if (event.SEL_TOP eq event.SEL_BOTTOM) then begin
          selected = event.SEL_TOP
          if (event.SEL_LEFT eq 3 and event.SEL_RIGHT eq 3) then begin
            widget_control, data.table_id, get_value=data_table
            status = strtrim(data_table[selected].status,2)
            run_process_id = where(run_process.title eq status,found)
            if found then begin
              run_process_id = run_process_id[0]
              acsproc_run_process, run_process_id, uniqname=data_table[selected].uniqname
            endif else acsproc_atv_metastruc_nav, selected
          endif else acsproc_atv_metastruc_nav, selected
        endif
      endif
    end
    
    'limit_data': begin
      if (event.value) then state.limit_data = event.select
      acsproc_atv_limit_data, state.limit_data
    end
    
    'limit_exposures': begin
      limit_exposures = event.select
      if (event.value) then limit_exposures = event.select
      acsproc_atv_limit_exposures, limit_exposures
      acsproc_limit_exposures, limit_exposures
    end
    
    'cycle':  begin
      state.cycle_droplist_index = event.index
      cycle = cycles[state.cycle_droplist_index]
      programs = get_num_subdir(hst_dataroot+cycle,/ascending)
      state.program_droplist_index = 0   
      program = programs[state.program_droplist_index]
      widget_control, state.program_droplist_id, set_value=programs     
      widget_control, state.program_droplist_id, set_droplist_select=state.program_droplist_index      
      acsproc_set_cycle_program, /refresh
      acsproc_get_metadata
    end
    
    'program':  begin
      state.program_droplist_index = event.index
      program = programs[state.program_droplist_index]
      acsproc_set_cycle_program, /refresh
      acsproc_get_metadata
    end
    
    'program_timestamp':  begin
      acsproc_program_timestamp
    end
    
    'spall':  begin
      acsproc_00_spall, Obj_New("acsproc_process", state.base_id, Title='SPALL', Message_xsize='120', Message_ysize='10', /CancelButton)
      ;acsproc_00_spall, Obj_New("acsproc_process", state.base_id, Title='SPALL', Message_xsize='120', Message_ysize='10', /CancelButton), /photmod, /galfit
    end
    
    'export':  begin
      acsproc_00_export, Obj_New("acsproc_process", state.base_id, Title='SPALL', Message_xsize='120', Message_ysize='10', /CancelButton)
    end
    
    'done':  acsproc_shutdown

    else:  print, 'No match for uvalue: ',uvalue  ; bad news if this happens
      
  endcase
  
end
 

;---------------------------------------------------------------------

function acsproc_goodones, candidate
  
  return, strtrim(candidate.lens,2) eq 'YES' or strtrim(candidate.lens,2) eq 'MAYBE'
  
end
;---------------------------------------------------------------------

pro acsproc_set_lens_status, status

  common acsproc_maskstruc

  case strtrim(status,2) of

    '': maskstruc.lens_status = 0
    'YES': maskstruc.lens_status = 1
    'MAYBE': maskstruc.lens_status = 2
    'NO': maskstruc.lens_status = 3
  
  endcase
      
end
;---------------------------------------------------------------------

function acsproc_get_lens_status, status

  common acsproc_maskstruc
  
  case status of

    0: lens_status = ''
    1: lens_status = 'YES'
    2: lens_status = 'MAYBE'
    3: lens_status = 'NO'
  
  endcase
  
  return, lens_status  
      
end
;---------------------------------------------------------------------

pro acsproc_update_maskflag, visit

; routine to update the acsproc metastruc on disk
;
  common acsproc_metastruc
  
  wh = where(metastruc.proposid eq visit.proposid and metastruc.visit eq visit.visit and metastruc.exposure eq visit.exposure, c)

  if (c gt 0) then begin
    metastruc[wh].maskflag = visit.maskflag
    metastruc[wh].fmaskflag = visit.fmaskflag
    metastruc[wh].jmaskflag = visit.jmaskflag
    acsproc_metastruc_update
  endif
      
end
 
;---------------------------------------------------------------------

pro acsproc_metastruc_update

; routine to update the acsproc metastruc on disk
;
  common acsproc_state
  common acsproc_data 
  common acsproc_metastruc 
  
   
  sxaddpar, header, 'TIMESTMP', systime(), 'metafile timestamp'
    
  mwrfits,  metastruc, metafile, header, /create
   
end
 
;---------------------------------------------------------------------

pro acsproc_limit_exposures, limit_exposures

; routine to limit the exposures shown in the data tables

  common acsproc_state
  common acsproc_data 
  common acsproc_metastruc 
   
  state.limit_exposures = limit_exposures
  
  metastruc = mrdfits(metafile, 1, /silent)
  widget_control, data.table_id, set_value = acsproc_data_table()
  
end

;---------------------------------------------------------------------

pro acsproc_cutout, process=process, selected=selected

; routine to generate the flt_cutouts

  common acsproc_state
  common acsproc_data  
  
  widget_control, state.table_id, get_value = status_table
  statuscode = (where(run_process.title eq 'CUTOUT'))[0]
  if keyword_set(selected) then count = n_elements(selected) else count = fix(status_table[statuscode].program_count)
  acsproc_01_cutout, process, statuscode, count, selected=selected
  acsproc_get_metadata

end

;---------------------------------------------------------------------

pro acsproc_register, process=process, selected=selected

; routine to generate the flt_cutouts

  common acsproc_state
  common acsproc_data  

  widget_control, state.table_id, get_value = status_table
  statuscode = (where(run_process.title eq 'REGISTER'))[0]
  if keyword_set(selected) then count = n_elements(selected) else count = fix(status_table[statuscode].program_count)
  acsproc_02_register, process, statuscode, count, selected=selected
  acsproc_get_metadata

end

;---------------------------------------------------------------------
pro acsproc_rectify, process=process, selected=selected

; routine to generate the flt_recti

  common acsproc_state
  common acsproc_data  

  widget_control, state.table_id, get_value = status_table
  statuscode = (where(run_process.title eq 'RECTIFY'))[0]
  if keyword_set(selected) then count = n_elements(selected) else count = fix(status_table[statuscode].program_count)
  acsproc_03_rectify, process, statuscode, count, selected=selected
  acsproc_get_metadata

end

;---------------------------------------------------------------------
pro acsproc_combine, process=process, selected=selected

; routine to generate the flt_comb

  common acsproc_state
  common acsproc_data  

  widget_control, state.table_id, get_value = status_table
  statuscode = (where(run_process.title eq 'COMBINE'))[0]
  if keyword_set(selected) then count = n_elements(selected) else count = fix(status_table[statuscode].program_count)
  acsproc_04_combine, process, statuscode, count, selected=selected
  acsproc_get_metadata

end

;---------------------------------------------------------------------
pro acsproc_mask, process=process, selected=selected

; routine to generate the flt_mask

  common acsproc_state
  common acsproc_data  

  widget_control, state.table_id, get_value = status_table
  statuscode = (where(run_process.title eq 'INITIAL MASK'))[0]
  if keyword_set(selected) then count = n_elements(selected) else count = fix(status_table[statuscode].program_count)
  acsproc_05_mask, process, statuscode, count, selected=selected
  acsproc_get_metadata

end

;---------------------------------------------------------------------

pro acsproc_bspline, process=process, selected=selected

; routine to generate bsplinecen

  common acsproc_state
  common acsproc_data  

  widget_control, state.table_id, get_value = status_table
  statuscode = (where(run_process.title eq 'BSPLINE'))[0]
  if keyword_set(selected) then count = n_elements(selected) else count = fix(status_table[statuscode].program_count)
  acsproc_06_bspline, process, statuscode, count, selected=selected
  acsproc_get_metadata

end

;---------------------------------------------------------------------
pro acsproc_bizzle, process=process, selected=selected

; routine to generate bsplinecen

  common acsproc_state
  common acsproc_data  
  
  widget_control, state.table_id, get_value = status_table
  statuscode = (where(run_process.title eq 'BIZZLE'))[0]
  if keyword_set(selected) then count = n_elements(selected) else count = fix(status_table[statuscode].program_count)
  acsproc_07_bizzle, process, statuscode, count, selected=selected
  acsproc_get_metadata

end

;---------------------------------------------------------------------

pro acsproc_photmod, process=process, selected=selected

; routine to generate bsplinecen

  common acsproc_state
  common acsproc_data  
  
  widget_control, state.table_id, get_value = status_table
  statuscode = (where(run_process.title eq 'PHOTMOD',c))[0]
  if keyword_set(selected) then count = n_elements(selected) else count = fix(status_table[statuscode].program_count)
  acsproc_08_photmod, process, statuscode, count, selected=selected, /flatweight
  ;acsproc_08_photmod, process, statuscode, count, selected=selected
  acsproc_get_metadata

end

;---------------------------------------------------------------------

function acsproc_get_metastruc, allfiles, metastruc0=metastruc0,redo=redo

nf = n_elements(allfiles)

metastruc = { $
 proposid: 0L, $
 targname: ' ', $
 uniqname: '', $
 visit: ' ', $
 exposure: 0L, $
 subexposure: 0L, $
 nexp: 0L, $
 ra_targ: 0.d0, $
 dec_targ: 0.d0, $
 equinox: 0., $
 filename: ' ', $
 instrume: ' ', $
 rootname: ' ', $
 gyromode: ' ', $
 date_obs: ' ', $
 time_obs: ' ', $
 expstart: 0.d0, $
 expend: 0.d0, $
 exptime: 0.d0, $
 expflag: ' ', $
 pa_v3: 0.d0, $
 filter1: ' ', $
 filter2: ' ', $
 filtcode: 0B, $
 aperture: ' ', $
 opus_ver: ' ', $
 cal_ver: ' ', $
 cutflag: 0L, $
 rectiflag: 0L, $
 combflag: 0L, $
 maskflag: 0L, $
 fmaskflag: 0L, $
 jmaskflag: 0L, $
 ra_tan_targ: 0., $
 dec_tan_targ: 0., $
 moffati_targ: 0., $
 moffatr_targ: 0., $
 moffatg_targ: 0., $
 moffatq_targ: 0., $
 moffata_targ: 0., $
 statuscode: 0L, $
 status: 'New Data', $
 lens: '', $
 timestamp: systime(), $
 mdrizsky: 0.}


metastruc = replicate(metastruc, nf)

for i = 0L, nf-1 do begin
    if file_test(allfiles[i]) then begin
        hdr = headfits(allfiles[i])
        hdr_filename = strtrim(sxpar(hdr, 'FILENAME'),2)
        if (keyword_set(metastruc0)) then wh = where(metastruc0.filename eq hdr_filename,processed) else processed = 0
        if (not processed) then begin
          metastruc[i].proposid = sxpar(hdr, 'PROPOSID')
          if metastruc[i].proposid eq 12209 then survey='BELLS' else survey='SLACS' 
          linenum = strtrim(sxpar(hdr, 'LINENUM'),2)
          spos = strpos(linenum, '.')
          metastruc[i].visit = strmid(linenum, 0, spos)
          if (metastruc[i].proposid eq 10494) then begin
             sdss_name = strarr(54) & sdss_name[1] = 'SDSSJ1627-0053' & sdss_name[2] = 'SDSSJ0912+0029' & sdss_name[3] = 'SDSSJ0737+3216' & sdss_name[4] = 'SDSSJ0959+0410' & sdss_name[5] = 'SDSSJ1402+6321' & sdss_name[6] = 'SDSSJ1630+4520' & sdss_name[7] = 'SDSSJ2321-0939' & sdss_name[8] = 'SDSSJ2300+0022' & sdss_name[9] = 'SDSSJ0216-0813' & sdss_name[10] = 'SDSSJ2238-0754' & sdss_name[11] = 'SDSSJ2303+1422' & sdss_name[12] = 'SDSSJ1420+6019' & sdss_name[13] = 'SDSSJ1250+0523' & sdss_name[14] = 'SDSSJ0956+5100' & sdss_name[15] = 'SDSSJ1205+4910' & sdss_name[52] = 'SDSSJ0912+0029' & sdss_name[53] = 'SDSSJ0956+5100'
             metastruc[i].targname =  sdss_name(fix(metastruc[i].visit))
          endif else if (metastruc[i].proposid eq 10174) then begin
             sdss_name = strarr(49) & sdss_name[0] = 'SDSSJ0819+4534' & sdss_name[1] = 'SDSSJ1205+4910' & sdss_name[2] = 'SDSSJ1259+6134' & sdss_name[3] = 'SDSSJ0803+4536' & sdss_name[4] = 'SDSSJ1432+6317' & sdss_name[5] = 'SDSSJ2202-0846' & sdss_name[6] = 'SDSSJ1330-0148' & sdss_name[7] = 'SDSSJ1246+0440' & sdss_name[8] = 'SDSSJ0953+5205' & sdss_name[9] = 'SDSSJ1636+4707' & sdss_name[10] = 'SDSSJ1451-0239' & sdss_name[11] = 'SDSSJ1251-0208' & sdss_name[12] = 'SDSSJ0948+0416' & sdss_name[13] = 'SDSSJ1547+5719' & sdss_name[14] = 'SDSSJ1333-0055' & sdss_name[15] = 'SDSSJ1443+0304' & sdss_name[16] = 'SDSSJ1420+6019' & sdss_name[17] = 'SDSSJ0959+0410' & sdss_name[18] = 'SDSSJ1416+5136' & sdss_name[19] = 'SDSSJ1306+0600' & sdss_name[20] = 'SDSSJ0802+4504' & sdss_name[21] = 'SDSSJ1155+6237' & sdss_name[22] = 'SDSSJ2251-0926' & sdss_name[23] = 'SDSSJ0935-0003' & sdss_name[24] = 'SDSSJ2343-0030' & sdss_name[25] = 'SDSSJ1702+3320' & sdss_name[26] = 'SDSSJ0216-0813' & sdss_name[27] = 'SDSSJ2300+0022' & sdss_name[28] = 'SDSSJ0813+4518' & sdss_name[29] = 'SDSSJ1025-0035' & sdss_name[30] = 'SDSSJ1402+6321' & sdss_name[31] = 'SDSSJ2238-0754' & sdss_name[32] = 'SDSSJ1250-0135' & sdss_name[33] = 'SDSSJ1250+0523' & sdss_name[34] = 'SDSSJ1204+0358' & sdss_name[35] = 'SDSSJ2303+1422' & sdss_name[36] = 'SDSSJ1618+4353' & sdss_name[37] = 'SDSSJ1535-0038' & sdss_name[38] = 'SDSSJ2321-0939' & sdss_name[39] = 'SDSSJ0737+3216' & sdss_name[40] = 'SDSSJ1627-0053' & sdss_name[41] = 'SDSSJ1630+4520' & sdss_name[42] = 'SDSSJ0912+0029' & sdss_name[43] = 'SDSSJ0037-0942' & sdss_name[44] = 'SDSSJ0956+5100' & sdss_name[45] = 'SDSSJ1117+0534' & sdss_name[46] = 'SDSSJ2052+0001' & sdss_name[47] = 'SDSSJ2347-0005' & sdss_name[48] = 'SDSSJ1718+6424'
              metastruc[i].targname =  sdss_name(fix(metastruc[i].visit)-1)
          endif else begin
            targname = strtrim(sxpar(hdr, 'TARGNAME'),2)
            pos = strpos(targname,'-COPY')
            if (pos gt 0) then targname = strmid(targname,0,pos)
            metastruc[i].targname = targname
          endelse
          pos = strpos(metastruc[i].targname,'SDSS')
          if (pos ge 0) then metastruc[i].uniqname= survey + strmid(metastruc[i].targname, pos+4) 
          metastruc[i].exposure = fix(strmid(linenum, spos+1))
          pattstep  = sxpar(hdr, 'PATTSTEP') 
          nexp  = sxpar(hdr, 'P1_NPTS') 
          if (pattstep gt 0 and nexp gt 0) then begin
            metastruc[i].subexposure = pattstep
            metastruc[i].nexp = nexp
          endif else begin
            metastruc[i].subexposure = 1
            metastruc[i].nexp = 1
          endelse
          metastruc[i].ra_targ = sxpar(hdr, 'RA_TARG')
          metastruc[i].dec_targ = sxpar(hdr, 'DEC_TARG')
          metastruc[i].equinox = sxpar(hdr, 'EQUINOX')
          metastruc[i].filename = hdr_filename
          metastruc[i].instrume = strtrim(sxpar(hdr, 'INSTRUME'),2)
          metastruc[i].rootname = strtrim(sxpar(hdr, 'ROOTNAME'),2)
          metastruc[i].gyromode = strtrim(sxpar(hdr, 'GYROMODE'),2)
          metastruc[i].date_obs = strtrim(sxpar(hdr, 'DATE-OBS'),2)
          metastruc[i].time_obs = strtrim(sxpar(hdr, 'TIME-OBS'),2)
          metastruc[i].expstart = sxpar(hdr, 'EXPSTART')
          metastruc[i].expend = sxpar(hdr, 'EXPEND')
          metastruc[i].exptime = sxpar(hdr, 'EXPTIME')
          metastruc[i].expflag = strtrim(sxpar(hdr, 'EXPFLAG'),2)
          metastruc[i].pa_v3 = sxpar(hdr, 'PA_V3')
          metastruc[i].filter1 = strtrim(sxpar(hdr, 'FILTER1'),2)
          metastruc[i].filter2 = strtrim(sxpar(hdr, 'FILTER2'),2)
          metastruc[i].filtcode = acsproc_filtcode(metastruc[i])
          metastruc[i].aperture = strtrim(sxpar(hdr, 'APERTURE'),2)
          metastruc[i].opus_ver = strtrim(sxpar(hdr, 'OPUS_VER'),2)
          metastruc[i].cal_ver = strtrim(sxpar(hdr, 'CAL_VER'),2)
          metastruc[i].mdrizsky = strtrim(sxpar(hdr, 'MDRIZSKY', count=has_mdsky),2)
          if not keyword_set(has_mdsky) then metastruc[i].mdrizsky = strtrim(sxpar(headfits(allfiles[i], exten=1), 'MDRIZSKY', count=has_mdsky),2)
          metastruc[i].cutflag = file_test(acsproc_datafile(metastruc[i], /cutout)) 
          metastruc[i].rectiflag = file_test(acsproc_datafile(metastruc[i], /rectify)) 
          metastruc[i].combflag = file_test(acsproc_datafile(metastruc[i], /combine))
          maskfile = acsproc_datafile(metastruc[i], /mask)
          if (file_test(maskfile)) then begin
            maskstruc = acsproc_get_maskstruc(maskfile)
            if maskstruc.has_photmod and file_test(acsproc_datafile(metastruc[i], /photmod)) then metastruc[i].maskflag = 4 $
            else if maskstruc.has_bizzle and file_test(acsproc_datafile(metastruc[i], /biz)) then metastruc[i].maskflag = 3 $
            else if maskstruc.has_modim then metastruc[i].maskflag = 2 $
            else metastruc[i].maskflag = 1
            print , metastruc[i].visit,metastruc[i].maskflag
            void = where(maskstruc.fmask ne 0,fmaskflag)
            void = where(maskstruc.jmask ne 0,jmaskflag)
            metastruc[i].fmaskflag = (fmaskflag gt 0)
            metastruc[i].jmaskflag = (jmaskflag gt 0)
            metastruc[i].lens = acsproc_get_lens_status(maskstruc.lens_status)
          endif else begin
            metastruc[i].maskflag = file_test(maskfile)
            metastruc[i].fmaskflag = 0L
            metastruc[i].jmaskflag = 0L
            metastruc[i].lens = acsproc_get_lens_status(0L)
            print , metastruc[i].visit," --> found no mask" 
          endelse
        endif else if (processed) then metastruc[i] = metastruc0[wh]
    endif else splog, "Cannot find FLT files at: "+allfiles[i]
    
endfor


return, metastruc
end

;---------------------------------------------------------------------

pro acsproc_set_metastruc, reset=reset

  common acsproc_data
        
  if not keyword_set(reset) and not file_test(metafile) then reset=1L
      
  if not keyword_set(reset) then metastruc0 = mrdfits(metafile, 1, /silent) else metastruc0=0L
            
  allfiles = file_search(dataroot+'*/*flt.fits')

  metastruc = acsproc_get_metastruc(allfiles,metastruc0=metastruc0)
              
  sxaddpar, header, 'TIMESTMP', systime(), 'metafile timestamp'
  sxaddpar, header, 'PROPOSID', program, 'Phase II ID'

  ofile = metafile
  
  mwrfits,  metastruc, ofile, header, /create
        

end

;---------------------------------------------------------------------

pro acsproc_set_metadata, newdata=newdata, refresh=refresh

; routine to compute the acsproc processes according to the state of the gui

  common acsproc_state
  common acsproc_data   
  common acsproc_metastruc  
  
  status_changed = 0L
    
  cutoutfile = acsproc_datafile(metastruc, /cutout)
  for q=0,n_elements(metastruc)-1 do begin
    statuscode = metastruc[q].statuscode
    if not file_test(cutoutfile[q]) then metastruc[q].cutflag = 0L
    if (metastruc[q].cutflag) then begin
      if (metastruc[q].subexposure eq 1) then cutout_header = headfits(cutoutfile[q],exten=1)
      cutout_needs_lcoord = (fix(sxpar(cutout_header, 'LCOORDX')) eq 0 or fix(sxpar(cutout_header, 'LCOORDY')) eq 0)
      has_register = (metastruc[q].moffati_targ ne 0.)
      initial_complete = (metastruc[q].cutflag and not cutout_needs_lcoord and metastruc[q].rectiflag and metastruc[q].combflag)
      ;splog,metastruc[q].visit, initial_complete, metastruc[q].cutflag, cutout_needs_lcoord, metastruc[q].rectiflag, metastruc[q].combflag, metastruc[q].maskflag, metastruc[q].fmaskflag, metastruc[q].jmaskflag
      if (initial_complete and metastruc[q].maskflag eq 4) then metastruc[q].statuscode = 11 $
      else if (initial_complete and metastruc[q].maskflag eq 3) then metastruc[q].statuscode = 10 $
      else if (initial_complete and metastruc[q].maskflag eq 2) then metastruc[q].statuscode = 9 $
      else if (initial_complete and metastruc[q].fmaskflag) or (initial_complete and metastruc[q].jmaskflag) then metastruc[q].statuscode = 8 $
      else if (initial_complete and metastruc[q].maskflag eq 1) then metastruc[q].statuscode = 7 $
      else if (metastruc[q].combflag) then metastruc[q].statuscode = 6 $
      else if (metastruc[q].rectiflag) then metastruc[q].statuscode = 5 $
      else if (has_register) then metastruc[q].statuscode = 4 $
      else if (not cutout_needs_lcoord) then metastruc[q].statuscode = 3 $
      else metastruc[q].statuscode = 2
    endif else metastruc[q].statuscode = 1
    
    metastruc[q].status = run_process[metastruc[q].statuscode].statuscode
    if keyword_set(refresh) and metastruc[q].maskflag then begin
      if (metastruc[q].subexposure eq 1) then begin
          maskfile = acsproc_datafile(metastruc[q], /mask)
          maskstruc = acsproc_get_maskstruc(maskfile)
          lens = acsproc_get_lens_status(maskstruc.lens_status)
      endif
      if (strtrim(metastruc[q].lens,2) ne strtrim(lens,2)) then begin
        metastruc[q].lens = lens
        if not status_changed then status_changed = 1L
      endif
    endif
    if not status_changed and (metastruc[q].statuscode ne statuscode) then status_changed = 1L
  endfor
  if status_changed then acsproc_metastruc_update

  data_table = acsproc_data_table()
  widget_control, data.table_id, set_value = data_table
  widget_control, data.table_id, table_ysize = n_elements(data_table)
  widget_control, data.header_id, set_value = acsproc_data_header()
  status_table = acsproc_status_table(status_table=status_table)
  
  if (n_elements(metastruc) gt 0) then begin
    widget_control, state.table_id, set_value = status_table
    if not keyword_set(refresh) then begin
      acsproc_atv, metastruc, dataroot=dataroot, cycle=cycle, program=program, reset=newdata
      acsproc_limit_exposures, state.limit_exposures
      acsproc_atv_limit_exposures, state.limit_exposures
    endif
  endif
  
end

;---------------------------------------------------------------------
pro acsproc_get_metadata, newdata=newdata, reset=reset

; routine to compute the acsproc processes according to the state of the gui

  common acsproc_state
  common acsproc_data   
  common acsproc_metastruc  

  if not file_test(metafile) then begin
    newdata=1L
    reset=1L
  endif
  
  if keyword_set(newdata) then acsproc_set_metastruc, reset=reset
  
  status_table = 0L
            
  if file_test(metafile) then begin

    metastruc = mrdfits(metafile, 1, header, /silent)
    acsproc_set_metadata, newdata=newdata
       
  endif
           
end

;--------------------------------------------------------------------
;    acsproc main program.  needs to be last in order to compile.
;---------------------------------------------------------------------

; Main program routine for acsproc.  If there is no current acsproc session,
; then run acsproc_startup to create the widgets.  

pro acsproc, cycle=this_cycle,program=this_program,offline=this_offline
         
  common acsproc_state
  common acsproc_data
  
  if (!d.name NE 'X' AND !d.name NE 'WIN' AND !d.name NE 'MAC') then begin
      print, 'Graphics device must be set to X, WIN, or MAC to work.'
      retall
  endif
  
  if (not (xregistered('acsproc', /noshow))) then begin
      acsproc_startup, cycle=this_cycle,program=this_program,offline=this_offline
      xmanager, 'acsproc', state.base_id, no_block = 1, cleanup = 'acsproc_shutdown'
      wset, !d.window
      acsproc_atv_resize, setting=1
      ;acsproc_atv_imagemode, /newmode
      acsproc_atv_imagemode, /combine, /newmode
      acsproc_atv_maskmode, 3, /newmode
      acsproc_atv_mode, 0, /newmode
  endif

end
