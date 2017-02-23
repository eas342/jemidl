Function delayedexp_SFH, galage, param=param
  return, galage/param.tau^2*exp(-galage^2/(2*param.tau^2))
end

Function exp_SFH, galage, param=param
  return, 1./param.tau*exp(-galage/param.tau)
end

Function cspspec, SFH_func, T_form, metal, param=param, _extra=extra
  ages = sspages(_extra=extra)
  w = where(ages lt T_form, nw)
  ages = [ages[w], T_form]
  ageaves = [0.5*ages[0],0.5*(ages+shift(ages,1))[1:*]]
  galages = T_form - ageaves
  psis = call_function(SFH_func, galages, param=param)
  dts = ages-shift(ages,1)
  dts[0] = ages[0]

  ssps = sspspec(ages, metal, /full, _extra=extra)
  csps = ssps.flux
  ms = total(psis*dts*ssps.ms)

  csphist = csps*transpose(rebin(psis*dts,nw+1,n_elements(ssps.wave)))
  csp = total(csphist,2)
  return, {wave:ssps.wave, flux:csp, ms:ms}
end
