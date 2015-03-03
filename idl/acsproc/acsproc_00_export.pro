;---------------------------------------------------------------------

pro acsproc_00_export, process, selected=selected

  common acsproc_data
  common acsproc_metastruc
  
  cmd = "hst_export -p "+program+" -m merge"
    
  process->Start, steps=steps
  process->Setmessage, 'Exporting Proposal ID '+program, /label
  process->Setmessage, '', /APPEND
  
  spawn, cmd, result
  
  for i=1,n_elements(result)-1 do begin
    process->Setmessage, result[i], /APPEND
    print, result[i]
  endfor
 
  wait, 2
  process->Destroy  
  
end

;---------------------------------------------------------------------

