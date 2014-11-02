pro def_pltsetup, setup
;
; Procedure def_pltsetup
;
; Description: Creates a structure that is used to pass plotting info
;               to several of the plotting routines.
;
; Inputs: setup (struct)          output structure containing setup info
;
; Revision history:
;  2003Apr11 Chris Fassnacht -- First version.
;

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, 'syntax: def_pltsetup, setup'
    print, ''
    return
endif

; Set up output structure

setup =  create_struct('row',0,$
                       'column',0,$
                       'fullxsz',0.0,$
                       'fullysz',0.0,$
                       'imgxsz',0.0,$
                       'imgysz',0.0,$
                       'fullxsz2',0.0,$
                       'fullysz2',0.0,$
                       'imgxsz2',0.0,$
                       'imgysz2',0.0,$
                       'labx',0.0,$
                       'laby',0.0,$
                       'labx2',0.0,$
                       'laby2',0.0,$
                       'axiscolor',16777215L,$
                       'namecolor',16777215L,$
                       'tcolor1',16777215L,$
                       'tcolor2',16777215L,$
                       'tcolor3',16777215L,$
                       'tcolor4',16777215L,$
                       'namex',0.0,$
                       'namey',0.0,$
                       'namex2',0.0,$
                       'namey2',0.0,$
                       'name','',$
                       'name2','',$
                       'name3','',$
                       'name4','',$
                       'lab1','',$
                       'lab2','',$
                       'lab3','',$
                       'lab4','',$
                       'lab5','',$
                       'lab6','',$
                       'lab7','',$
                       'doticks',1B,$
                       'doxlab',1B,$
                       'doylab',1B)

end

