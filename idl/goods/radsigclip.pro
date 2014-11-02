pro radsigclip, xin, yin, nsig, mask, xbin=xbin, radrms=radrms, ngood=ngood 
;
; Procedure radsigclip
;
; Description: Takes input x and y vectors and computes a the RMS
;               scatter in regularly-spaced bins in x.  It is these
;               RMS values that are used to create the bad pixel mask,
;               after doing an iterative clipping in each bin.
;
; Inputs: xin       (floatarray)  input x vector
;         xin       (floatarray)  input x vector
;         nsig      (float)       level of sigma clipping.  e.g., for
;                                  3-sigma clipping, set nsig = 3.0
;         [xbin=]   (float)       binsize for bins in x
;
; Outputs accessed through passed parameters:
;         mask      (floatarray)  mask of bad pixels (bad pixels set
;                                  to 0)
;         [radrms=] (double)      radially-binned rms of clipped image
;         [ngood=]  (int)         number of good pixels after clipping
;                                  converges
;
; Revision history:
;  2002Apr23 Chris Fassnacht -- First working version.
;

; Check input format

if n_params() lt 4 then begin
    print, ''
    print, $
      'syntax: radsigclip, xin, yin, nsig, mask [, xbin=, radrms=, ngood=]'
    print, ''
    return
endif

; Find the limits for the x binning

if not keyword_set(xbin) then xbin=3.0
xmax = max(xin,min=xmin)
xstepmin = floor(xmin)
xstepmax = ceil(xmax)
nbin = ceil((xmax - xmin)/xbin)

; Set up bin boundaries

xstep = findgen(nbin+1)*xbin + xstepmin

; Set up vectors to be filled

rmsbin = fltarr(nbin)
rmsall = 0.0 * yin
mask = 0.0 * yin + 1.0


; For each bin, get rms value

if not keyword_set(nsig) then nsig=5.0
for i=0L,nbin-1 do begin
    wsel = where (xin ge xstep[i] and xin lt xstep[i+1],nsel)
    case nsel of
        0: rmsbin[i] = 0.0
        1: rmsbin[i] = yin[wsel]
        else: begin 
            rmsbin[i] = stddev(yin[wsel],/double) 
            rmsall[wsel] = rmsbin[i]
        end
    endcase
endfor

; Set mask equal to 0.0 for points greater than nsig times their local
; RMS

wbad = where(abs(yin) gt (nsig * rmsall), nbad)
if nbad gt 0 then mask[wbad] = 0.0
print, 'radsigclip: Masking out ',nbad,' bad pixels.'

end

