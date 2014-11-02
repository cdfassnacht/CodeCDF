pro get_photoz, id, photcat=photcat, matchcat=matchcat, catdir=catdir

;
; Procedure get_photoz
;
; Description: Given an ACS id and the location of the photo-z catalog
;               and its associated "match" catalog (which matches photo-z id's
;               with ACS id's), extracts the photo-z information.
;               a galaxy and residual images created by subtracting off
;               some form of fit to the galaxy light distribution.
;
; Inputs: id          (intarr)     ACS id(s)
;         [photcat=]  (string)     name of photo-z catalog -- default is
;                                   'photz_only_s_v2.cat'
;         [matchcat=] (string)     name of mathc catalog -- default is
;                                   'CDFS_assoc_v1.7.cat'
;         [catdir=]   (string)     name of directory containing the catalogs
;                                   -- default is 'Photoz'
;
; Revision history:
;  2003Apr05 Chris Fassnacht -- First working version.
;

; Check input format

if n_params() lt 1 then begin
    print, ''
    print, $
      'syntax: get_photoz, id [,photcat=, matchcat=, catdir=]
    print, ''
    print, " Default value for photcat is 'photz_only_s_v2.cat'"
    print, " Default value for matchcat is 'CDFS_assoc_v1.7.cat'"
    print, " Default value for catdir is 'Photoz'"
    print, ''
    return
endif

; More input checking

nid = n_elements(id)
if nid eq 0 then begin
   print, ''
   print, '**** ERROR: The passed id parameter has no elements'
   return
endif

; Set directory names and defaults

if not keyword_set(photcat) then photcat = 'photz_only_s_v2.cat'
if not keyword_set(matchcat) then matchcat = 'CDFS_assoc_v1.7.cat'
if not keyword_set(catdir) then catdir = 'Photoz'
pzcat = catdir+'/'+photcat
mcat = catdir+'/'+matchcat

; Read in catalogs with readcol

; readcol, pzcat, pid,zphot,minz,maxz,gtype,odds, format='l,f,f,f,f,f'
;npz = n_elements(pid)
readcol, mcat, aid,zid,dist,format='l,x,x,l,x,x,f'
nm = n_elements(aid)

; Find ids in the "match" catalog


for i=0,nid-1 do begin
   print, i
   for j=0,nm-1 do begin
   endfor
endfor

end 
