Pro testdrizzle
  deimos_getspec, $
     '/home/jmeyers314/scp1/CL-spec/Keck/CL-R/RA/spec1d.R.000.R-012.fits.gz', $
     wave, flux, ivar
  inbins = wave
  influx=flux[1:*]
  inivar=ivar[1:*]
  wfinite = where(finite(inivar), complement=wnotfinite)
  if wnotfinite[0] ne -1 then inivar[wnotfinite]=0.d
  outbins = dindgen(501)/500.d * 5000.d + 5000.d
  preflux = dblarr(500)*0.d
  preivar = dblarr(500)*0.d
  
  blah=dblarr(10)
  
  blah = drizzle1d( inbins, influx, inivar, outbins, preflux, preivar )
  outflux = blah[0:499]
  outivar = blah[500:999]
  splot, outbins[0:499], outflux
end
