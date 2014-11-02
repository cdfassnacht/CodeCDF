pro read_galaxy_subtree_rich9, filebasename, FirstFile, LastFile, TotNTrees, TotNGals, FirstGalIDs, Gal
;; the magnitudes are
;; sdss-u, sdss-g, sdss-r, sdss-i, sdss-z, WFC2-F606, WFC2-F814, UKIRT-J, UKIRT-Ks
;; mags are in the observed frame.
;; no dust model
Hubble_h=0.7

scale = 1.0
scalemass = scale^3.
scaleenergy = scale^5.
MassUnitinMsun=1.d10
duston=0


  Gstruct = {$
              Type                  : 0L, $
              HaloIndex             : 0L, $
              SnapNum               : 0L, $
		  HaloID			    : 0UL, $
              CentralMvir           : 0.0, $
              Mvir                  : 0.0, $
              Rvir                  : 0.0, $
              Vvir                  : 0.0, $
              Vmax                  : 0.0, $
		  Pos					: fltarr(3), $
		  Vel					: fltarr(3), $
              ColdGas               : 0.0, $
              StellarMass           : 0.0, $
              BulgeMass             : 0.0, $
              HotGas                : 0.0, $
              EjectedMass           : 0.0, $
              BlackHoleMass         : 0.0, $
              MetalsColdGas         : 0.0, $
              MetalsStellarMass     : 0.0, $
              MetalsBulgeMass       : 0.0, $
              MetalsHotGas          : 0.0, $
              MetalsEjectedMass     : 0.0, $
              Sfr                   : 0.0, $
              SfrBulge              : 0.0, $
              XrayLum               : 0.0, $
              DiskRadius            : 0.0, $
              CoolingRadius         : 0.0, $
              Mag                   : fltarr(9), $
              MagBulge              : fltarr(9), $
              MagDust               : fltarr(9) $
            }

    TotNTrees = 0L & TotNGals = 0L
    close, 1
    for fnr = FirstFile, LastFile do begin
      fname = strcompress(filebasename+'_'+string(fnr), /remove_all)
      openr, 1, fname
      Ntrees = 0L     & readu, 1, Ntrees
      NtotGals = 0L   & readu, 1, NtotGals
      print,fname,Ntrees,NtotGals
      TotNTrees = TotNTrees + Ntrees
      TotNGals =  TotNGals  + NtotGals
      close, 1
    endfor
    print, TotNGals
    G = replicate(Gstruct, TotNGals)
    ;G = replicate(Gstruct, 5)
        print, 'Made it this far'
	GGalsPerTree = lonarr(TotNTrees)
	FirstGalIDs = lonarr(TotNTrees)

    Offset = 0L
	Treeoffset = 0L
    for fnr = FirstFile, LastFile do begin
      fname = strcompress(filebasename+'_'+string(fnr), /remove_all)
      openr, 1, fname
      Ntrees = 0L     & readu, 1, Ntrees
      NtotGals = 0L   & readu, 1, NtotGals
      GalsPerTree = lonarr(Ntrees)
      readu, 1, GalsPerTree
      print, fnr, ":   Reading N=", NtotGals, "   galaxies."
      GG = replicate(Gstruct, NtotGals)
      readu, 1, GG
      G(offset:offset+NtotGals-1) = GG(*)
	  GGalsPerTree[Treeoffset:Treeoffset+Ntrees-1] = GalsPerTree[*]
      offset = offset + NtotGals
	  Treeoffset = Treeoffset + Ntrees
      close, 1
    endfor

	for i=1L, TotNTrees-1 do begin
		FirstGalIDs[i]=total(GGalsPerTree[0:i-1])
	endfor
	FirstGalIDs[0] = 0

;;;;;;;; some setup

  G.CentralMvir = G.CentralMvir*scalemass/Hubble_h*MassUnitinMsun

  G.Mvir = G.Mvir*scalemass/Hubble_h*MassUnitinMsun
  G.Rvir = G.Rvir*scale
  G.Vvir = G.Vvir*scale
  G.Vmax = G.Vmax*scale

  G.ColdGas       = G.ColdGas*scalemass/Hubble_h*MassUnitinMsun
  G.StellarMass   = G.StellarMass*scalemass/Hubble_h*MassUnitinMsun
  G.BulgeMass     = G.BulgeMass*scalemass/Hubble_h*MassUnitinMsun
  G.HotGas        = G.HotGas*scalemass/Hubble_h*MassUnitinMsun
  G.EjectedMass   = G.EjectedMass*scalemass/Hubble_h*MassUnitinMsun
  G.BlackHoleMass = G.BlackHoleMass*scalemass/Hubble_h*MassUnitinMsun

  G.MetalsColdGas       = G.MetalsColdGas*scalemass/Hubble_h*MassUnitinMsun
  G.MetalsStellarMass   = G.MetalsStellarMass*scalemass/Hubble_h*MassUnitinMsun
  G.MetalsBulgeMass     = G.MetalsBulgeMass*scalemass/Hubble_h*MassUnitinMsun
  G.MetalsHotGas        = G.MetalsHotGas*scalemass/Hubble_h*MassUnitinMsun
  G.MetalsEjectedMass   = G.MetalsEjectedMass*scalemass/Hubble_h*MassUnitinMsun

  G.Sfr                = G.Sfr*scalemass
  G.SfrBulge           = G.SfrBulge*scalemass
  G.XrayLum            = G.XrayLum + alog10(scaleenergy)
  G.DiskRadius         = G.DiskRadius*scale
  G.CoolingRadius      = G.CoolingRadius*scale
	G.Pos				= G.Pos/1.e3

  if (duston eq 0) then begin
    G.Mag  = G.Mag - 2.5*alog10(scaleenergy)
  endif else begin
    G.Mag  = G.MagDust - 2.5*alog10(scaleenergy)
  endelse
  G.MagBulge = G.MagBulge - 2.5*alog10(scaleenergy)

	Gal = G
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
