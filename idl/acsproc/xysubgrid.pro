;+
;
; NAME: xysubgrid
;
; PURPOSE: return subsampled x and y grid images.
;
; USAGE:
;   xysubgrid, nx, ny, xsub, ysub, ximg, yimg
;
; INPUTS:
;   nx, ny: extent of grid images in WHOLE pixels.
;   xsub, ysub: subsampling factors in x and y.
;
; RETURNS:
;   ximg, yimg: nx*xsub X ny*ysub grid coordinate images.
;
; WRITTEN: 2006may abolton@cfa
;
;-

pro xysubgrid, nx, ny, xsub, ysub, ximg, yimg

nx_loc = fix(abs(nx)) > 1
ny_loc = fix(abs(ny)) > 1
xsub_loc = fix(abs(xsub)) > 1
ysub_loc = fix(abs(ysub)) > 1

xbase = 0.5 * (2.*findgen(nx_loc*xsub_loc) + 1. - float(xsub_loc))/float(xsub_loc)
ybase = 0.5 * (2.*findgen(ny_loc*ysub_loc) + 1. - float(ysub_loc))/float(ysub_loc)

ximg = xbase # replicate(1., ny_loc*ysub_loc)
yimg = replicate(1., nx_loc*xsub_loc) # ybase

return
end
