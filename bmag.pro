Function bmag, wave, lumdens
  resolve_routine, 'spectrum__define', /no_recompile
  resolve_routine, 'dopplershift__define', /no_recompile
  restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/filters.sav'
  restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/zeropointspecs.sav'
  resolve_obj, filters.acs_f775w
  resolve_obj, ABspec
  
  d_L = 1.e-5
  d_L *= 3.0857d24 ;; Mpc -> cm

  fw = (filters.bessellb)->wavelength()
  ft = (filters.bessellb)->throughput()

  w = where( wave ge min(fw) and wave le max(fw) )
  wave_int = wave[w]
  flux = lumdens[w] / (4.*!dpi*d_L^2)
  ft_interp = interpol( ft, fw, wave_int )
  
  ABflux = ABspec->flux( frame='obs', unit='flambda' )
  ABwave = ABspec->wavelength( frame='obs' )
  ABflux_int = interpol( ABflux, ABwave, wave_int )
  
  dwave = wave_int - shift( wave_int, 1 )
  dwave[0] = dwave[1]

  same = wave_int*ft_interp*dwave
  num = total(flux*same)
  den = total(ABflux_int*same)

  return, -2.5*alog10(num/den)
end
