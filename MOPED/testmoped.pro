Function modelfn, P
  age = P[0]
  metal = P[1]
  
  spec = bc03spec(age, metal)
  wave = spec->wavelength()
  flux = spec->flux()
  obj_destroy, spec
  w=where(wave gt 3350 and wave lt 4500)
  wave=wave[w]
  flux=flux[w]
  return, flux
end

Function dmodelfn, P
  flux0 = call_function('modelfn', P)
  P1 = P
  P1[0] += 0.01*P[0]
  flux1 = call_function('modelfn', P1)
  dflux1 = (flux1-flux0)/(0.01*P[0])
  P2 = P
  P2[1] += 0.01*P[1]
  flux2 = call_function('modelfn', P2)
  dflux2 = (flux2-flux0)/(0.01*P[1])
  return, [[dflux1],[dflux2]]
end

Pro testmoped
  ;;first we need to create a fake spectrum to try and capture...
  s2n = 10.
  testparams = [1.e9, 0.04]
  testflux = call_function('modelfn', testparams)
  err = double(randomu(seed,n_elements(testflux),/normal)*testflux/s2n)
  ;;got spectrum now

  ivar = (s2n / testflux)^2
  ;; now to fit it...
  fiducialparams = [5.e8, 0.01]
  min = fiducialparams
;  for i=0, 5 do begin
;     moped = obj_new('moped', 'modelfn','dmodelfn', /norm, $
;                     2, ivar, min, testflux+err, /plot, /verbose)
;     ;min = moped->amoeba(0.001, min, [1.e9, 0.01])
;     min = moped->grid(
;     print, min
;     obj_destroy, moped
;
;  endfor

  agerange = [1.e8, 1.e10]
  metalrange = [1.e-3, 5.e-2]
  
  lages = findgen(25)/24.*(alog10(agerange[1])-alog10(agerange[0]))+alog10(agerange[0])
  ages = 10^lages
  lmetals = findgen(25)/24.*(alog10(metalrange[1])-alog10(metalrange[0]))+alog10(metalrange[0])
  metals = 10^lmetals
  moped = obj_new('hybridmoped', 'modelfn','dmodelfn', /norm, $
                  2, ivar, min, testflux+err, /plot, /verbose)
  void = moped->grid([Ptr_New(ages), Ptr_New(metals)], like=like)
  
atv, like

end
