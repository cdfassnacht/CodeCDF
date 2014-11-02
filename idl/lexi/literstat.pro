FUNCTION lmode, data
   
; Naive script to determine the mode, only take the first 
; in case of double-valued result...
   
   hist = histogram(data) 
   ifin = n_elements(hist)
   xvals = findgen(ifin)+min(data)+0.5
   xvals = xvals(where(hist eq max(hist)))
   return, xvals(0)
   
END 

PRO literstat, data, out, nsigrej=nsigrej,$
               maxiter=maxiter, lower=lower, upper=upper, silent=silent
   
; Program based on dimsum's iterstat:
; Created lam98jan18ucb
; 
; Converted to idl, and made output into a structure.
; 
; procedure iterstat(image)
; 
; # Script to find image statistics excluding deviant pixels
; # 4 August 1992 by John Ward
; # Minor modifications 4 August 1992 MD
; # Various subsequent variations.
; # Latest revision:  18 Aug 1993 MD
;
   
   IF n_params() LE 0 THEN BEGIN 
      print,' '
      print,'literstat, data, out, [nsigrej=], [maxiter=],'
      print,'           [lower=], [upper=], [silent]'
      print,' '
      print,'  N.b. -- The output array ''out'' is a structure:'
      print,'    out.NPIX .MEAN .MEDIAN .SIGMA .MODE .MIN .MAX'
      print,' '
      return 
   ENDIF 
   
   out =  create_struct('Npix',0.,$
                        'Mean',0.,$
                        'Median',0.,$
                        'Sigma',0.,$
                        'Mode',0.,$
                        'min',0.,$
                        'max',0.)
   
   IF keyword_set(nsigrej) EQ 0. THEN nsigrej=5. ; N * sigma for limits
   IF keyword_set(maxiter) EQ 0  THEN maxiter=10 ; Max # of iterations
   IF NOT keyword_set(lower) THEN lower=min(data)
   IF NOT keyword_set(upper) THEN upper=max(data)
   
   data = float(data)           ; Verify that data are floated
   sz = size(data)
   npx = sz(n_elements(sz)-1)   ; Number of relevant pixels
;   lind = lindgen(npx)
   lind = where(data GE lower AND data LE upper)
;;;   data = data(lind)
;    IF n_elements(where(finite(data) EQ 1)) NE n_elements(data) THEN BEGIN 
;       print,'There are NaNs in the array'
;       lind = where(finite(data) EQ 1)
;    ENDIF 
   mn  =  mean  (data(lind))    ; Mean
   sig =  stdev (data(lind))    ; Standard Deviation
   med =  median(data(lind))    ; Median
   mde =  lmode (data(lind))    ; Mode, cf 'lmode' below
   IF NOT keyword_set(silent) THEN print, ' '
   IF NOT keyword_set(silent) THEN print,'    Iter       Npix        Mean         Sigma        Median       Mode'

   Tpx =  0
   m = 1
   WHILE (m LE maxiter AND Tpx NE npx) DO BEGIN 
      npx = Tpx
      ll =  mn - (nsigrej*sig)
      ul =  mn + (nsigrej*sig)
      IF (ll NE min(data(lind)) AND ll LE lower) THEN ll=lower
      IF (ul NE max(data(lind)) AND ul GE upper) THEN ul=upper
      lind = where(data GE ll AND data LE ul)
      IF total(lind) EQ -1 THEN BEGIN 
         print,'Error -- all the pixels masked...'
         return
      ENDIF 
      Tpx =  n_elements(lind)   ; Number of relevant pixels
      mn  =  mean  (data(lind)) ; Mean
      sig =  stdev (data(lind)) ; Standard Deviation
      med =  median(data(lind)) ; Median
      mde =  lmode (data(lind)) ; Mode, cf 'lmode' below
      IF NOT keyword_set(silent) THEN print,m,tpx,mn,sig,med,mde
      m = m+1
   ENDWHILE 
   
   IF NOT keyword_set(silent) THEN print, ' '
   
   out.npix = npx
   out.mean = mn
   out.median = med
   out.sigma = sig
   out.mode = mde
   out.min = ll
   out.max = ul
   
   IF NOT keyword_set(silent) THEN help,/structure,out
   
END 
