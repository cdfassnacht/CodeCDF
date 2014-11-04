pro photmodupdate, cycle=cycle, program=program
  common acsproc_maskstruc
  if not keyword_set(cycle) then cycle='18'
  if not keyword_set(program) then program='12210'
  hst_dataroot=getenv('HST_DATAROOT')
  cycle_program = strtrim(cycle,2) + '/' + strtrim(program,2)
  fileroot = ['57/SLACSJ1040+3626_12210_56040_F814W_1','65/SLACSJ1103+3625_12210_56054_F814W_1','72/SLACSJ1127+2312_12210_56066_F814W_1','82/SLACSJ1213+2930_12210_56086_F814W_1']
  photmodfile = djs_filepath(fileroot+'_photmod.fits', root=hst_dataroot, subdir=cycle_program)
  maskfile = djs_filepath(fileroot+'_mask.fits', root=hst_dataroot, subdir=cycle_program)
  for i=0,n_elements(maskfile)-1 do begin
	maskstruc = acsproc_get_maskstruc(maskfile[i], header=mask_header)
	help,maskstruc,/struc
	has_photmod = maskstruc.has_photmod
	maskstruc.has_photmod = file_test(photmodfile[i])
	if has_photmod ne maskstruc.has_photmod then begin
		print, i,maskfile[i]
		print , 'photmod manual update mask',maskstruc.has_photmod
		acsproc_maskstruc_write, maskfile[i]
	endif
endfor
end
