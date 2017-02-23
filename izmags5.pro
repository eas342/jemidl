;returns the i_775 and z_850 magnitudes for an observer frame spectrum.

Function izmags5, $
   wave, $
   flux

  common jem$_izmags5, _izmag_struct
  resolve_routine, 'spectrum__define', /no_recompile
  resolve_routine, 'dopplershift__define', /no_recompile
  if n_elements(_izmag_struct) eq 0 then begin
     restore, '/home/scpdata01/FILTERS/filters.sav'
     restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/zeropointspecs.sav'
     resolve_obj, filters.f775w
     resolve_obj, ABspec
     resolve_routine, 'dopplershift__define', /no_recompile
     _izmag_struct = { ABspec:ABspec, $
                       F775W: filters.f775w, $
                       F850LP: filters.f850lp }
     tags = tag_names(filters)
     wtags = where(~member(tags, ['F775W','F850LP']))
     for i=0, n_elements(wtags)-1 do obj_destroy, filters.(wtags[i])
     obj_destroy, [vegaspec, STspec]
  endif

  ifw = (_izmag_struct.F775W)->wavelength()
  ift = (_izmag_struct.F775W)->throughput()
  zfw = (_izmag_struct.F850LP)->wavelength()
  zft = (_izmag_struct.F850LP)->throughput()

  ABflux = _izmag_struct.ABspec->flux( frame='obs', unit='flambda' )
  ABwave = _izmag_struct.ABspec->wavelength( frame='obs' )
  if min(wave) lt min(ifw) and max(wave) gt max(ifw) then begin
     iw = where( wave ge min(ifw) and wave le max(ifw) )
     iwave = wave[iw]
     iflux = flux[iw]
     ift_interp = interpol( ift, ifw, iwave )
     iABflux = interpol( ABflux, ABwave, iwave )
     idwave = iwave - shift( iwave, 1 )
     idwave[0] = idwave[1]
     isame = iwave*ift_interp*idwave
     inum = total(iflux*isame)
     iden = total(iABflux*isame)
     imag = -2.5*alog10(inum/iden)
  endif else imag = !values.f_nan

  if min(wave) lt min(zfw) and max(wave) gt max(zfw) then begin
     zw = where( wave ge min(zfw) and wave le max(zfw) )
     zwave = wave[zw]
     zflux = flux[zw]
     zft_interp = interpol( zft, zfw, zwave )
     zABflux = interpol( ABflux, ABwave, zwave )
     zdwave = zwave - shift( zwave, 1 )
     zdwave[0] = zdwave[1]
     zsame = zwave*zft_interp*zdwave
     znum = total(zflux*zsame)
     zden = total(zABflux*zsame)
     zmag = -2.5*alog10(znum/zden)
  endif else zmag = !values.f_nan

  return, [imag, zmag]
end
