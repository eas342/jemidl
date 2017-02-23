Function bc03mags, $
   age, $
   metallicity, $
   z, $
   filters, $
   mass=mass, $
   zpspec=zpspec, $
   zcluster=zcluster, $
   host_ebv=host_ebv, $
   host_Rv=host_Rv, $
   MW_ebv=MW_ebv, $
   MW_Rv=MW_Rv
  
  resolve_routine, 'spectrum__define', /no_recompile
  resolve_routine, 'dopplershift__define', /no_recompile

  if n_elements(z) eq 0 then z=0.
  if n_elements(zcluster) eq 0 then zcluster=z
  if n_elements(host_ebv) eq 0 then host_ebv=0.
  if n_elements(host_Rv) eq 0 then host_Rv = 3.1
  if n_elements(MW_ebv) eq 0 then MW_ebv = 0.
  if n_elements(MW_Rv) eq 0 then MW_Rv = 3.1
  if n_elements(zpspec) eq 0 then begin
     restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/zeropointspecs.sav'
     resolve_obj, ABspec
     zpspec=ABspec
  endif

  spec = bc03spec( age, metallicity, mass=mass, /erg )
  lumdens = spec->flux()
  wave = spec->wavelength()
  
  wave_obs = wave * (1.0+z)

  if host_ebv ne 0. then lumdens *= ccm( wave, host_ebv, host_Rv )
  if MW_ebv ne 0. then lumdens *= ccm( wave_obs, MW_ebv, MW_Rv )

  if zcluster eq 0 then d_L = 1.e-5 $   ;; 10 kpc in units of Mpc
  else d_L = lumdist(zcluster, /silent) ;;output in units of Mpc
  d_L *= 3.0857d24 ;; Mpc -> cm conversion  (since comparison spectrum is in erg/cm^2/s/Ang)

  out = dblarr(n_elements(filters))
  for ifilter=0, n_elements(filters)-1 do begin
     fw = filters[ifilter]->wavelength()
     ft = filters[ifilter]->throughput()
     
     w = where( wave_obs ge min(fw) and wave_obs le max(fw) )
     wave_obs = wave_obs[w]
     flux = lumdens[w] / (4.*!dpi*d_L^2*(1.0+z))
     
     ft_interp = interpol( ft, fw, wave_obs )
     
     zpflux = zpspec->flux()
     zpwave = zpspec->wavelength()
     zpflux = interpol( zpflux, zpwave, wave_obs )
     
     dwave = wave_obs - shift( wave_obs, 1)
     dwave[0] = dwave[1]
     
     same = wave_obs*ft_interp*dwave
     num = total(flux*same)
     den = total(zpflux*same)
     out[ifilter] = -2.5*alog10(num/den)
  endfor
  
  return, out
end
