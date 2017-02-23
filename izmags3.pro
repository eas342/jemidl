Function izmags3, mass, z, age, metal
  common jem$_izmags3, ABspec, filters
  resolve_routine, 'spectrum__define', /no_recompile
  resolve_routine, 'dopplershift__define', /no_recompile
  if n_elements(_izmag_struct) eq 0 then begin
     restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/filters.sav'
     restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/zeropointspecs.sav'
     resolve_obj, filters.acs_f775w
     resolve_obj, ABspec
  endif

  spec = bc03spec( age, metal )
  wave = spec->wavelength(frame='obs')
  flux = spec->flux( frame='obs', unit='flambda' )
  wave *= 1.+z ;; redshift
  flux *= mass
  flux *= 3.826d33 ;; L_sun -> erg conversion
  flux /= 1.+z ;; more redshift...
  d_L = lumdist(z, /silent)
  d_L *= 3.0857d24 ;; Mpc -> cm conversion
  flux /= 4*!dpi*d_L^2

  ifw = (filters.acs_f775w)->wavelength()
  ift = (filters.acs_f775w)->throughput()
  zfw = (filters.acs_f850lp)->wavelength()
  zft = (filters.acs_f850lp)->throughput()

  iw = where( wave ge min(ifw) and wave le max(ifw) )
  iwave = wave[iw]
  iflux = flux[iw]
  zw = where( wave ge min(zfw) and wave le max(zfw) )
  zwave = wave[zw]
  zflux = flux[zw]

  ift_interp = interpol( ift, ifw, iwave )
  zft_interp = interpol( zft, zfw, zwave )

  ABflux = ABspec->flux( frame='obs', unit='flambda' )
  ABwave = ABspec->wavelength( frame='obs' )
  iABflux = interpol(ABflux, ABwave, iwave)
  zABflux = interpol(ABflux, ABwave, zwave)

  idwave = iwave-shift(iwave,1)
  idwave[0] = idwave[1]
  zdwave = zwave-shift(zwave,1)
  zdwave[0] = zdwave[1]
  
  zsame = zwave*zft_interp*zdwave
  znum = total(zflux*zsame)
  zden = total(zABflux*zsame)

  isame = iwave*ift_interp*idwave
  inum = total(iflux*isame)
  iden = total(iABflux*isame)

  return, -2.5*alog10([inum/iden,znum/zden])
  
end
