pro eidolwrap,eidolfile,istart,_extra=_extra

; left mouse button to move forward in your list
; right mouse button to move back in your list
; middle mouse button to stop

    readcol,eidolfile,rahms,dechms,title,format='a,a,a'
    neid=n_elements(title)

    if n_elements(istart) eq 0 then istart=0l
    if istart lt 0 or istart gt neid then begin
        print,'give a i-value in range'
        return
    endif

    i=istart
    mousebutt = 0
    while mousebutt ne 2 do begin
        print,strcompress(i)+': '+title[i]
        eidol,rahms[i],dechms[i],title=title[i],_extra=_extra
        cursor,x,y,/wait
        mousebutt = !mouse.button
        if mousebutt eq 1 and i lt neid then i=i+1
        if mousebutt eq 4 and i gt -1 then i=i-1
    endwhile
        
end 


