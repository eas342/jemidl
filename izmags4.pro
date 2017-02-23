;;wavelength in units of Angstroms
;;luminosity density in units of ergs/s/Angstrom
;;assume standard cosmology (lumdist)

Function izmags4, $
   wave_emit, $           ;;rest(emitted) wavelengths in Angstroms
   lumdens, $             ;;luminosity density of source (erg/s/Ang)
   z, $                   ;;redshift for velocity correction
   zcluster=zcluster, $   ;;redshift for distance correction
   host_ebv=host_ebv, $   ;;E(B-V) for host in rest-frame
   host_Rv=host_Rv, $     ;;total-to-selective extinction coefficient for host
   MW_ebv=MW_ebv, $       ;;E(B-V) for MW in obs-frame
   MW_Rv=MW_Rv            ;;total-to-selective extinction coefficient for MW

  if n_elements(z) eq 0 then z=0.
  if n_elements(zcluster) eq 0 then zcluster=z
  if n_elements(host_ebv) eq 0 then host_ebv=0.
  if n_elements(host_Rv) eq 0 then host_Rv = 3.1
  if n_elements(MW_ebv) eq 0 then MW_ebv=0.
  if n_elements(MW_Rv) eq 0 then MW_Rv=3.1
  if ~finite(z) or ~finite(zcluster) then return, [!values.f_nan, !values.f_nan]

  common jem$_izmags4, _izmag_struct
  resolve_routine, 'spectrum__define', /no_recompile
  resolve_routine, 'dopplershift__define', /no_recompile
  if n_elements(_izmag_struct) eq 0 then begin
     restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/filters.sav'
     restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/zeropointspecs.sav'
     resolve_obj, filters.acs_f775w
     resolve_obj, ABspec
     resolve_routine, 'dopplershift__define', /no_recompile
     _izmag_struct = { ABspec:ABspec, $
                       ACS_F775W: filters.acs_f775w, $
                       ACS_F850LP: filters.acs_f850lp }
     tags = tag_names(filters)
     wtags = where(~member(tags, ['ACS_F775W','ACS_F850LP']))
     for i=0, n_elements(wtags)-1 do obj_destroy, filters.(wtags[i])
     obj_destroy, [vegaspec, STspec]
  endif


  wave_obs = wave_emit * (1.0+z)

  if zcluster eq 0 then d_L = 1.0e-5 $ ;; 10 kpc in units of Mpc
  else d_L = lumdist(zcluster, /silent) ;;output in units of Mpc
  d_L *= ( (1.0+z)/(1.0+zcluster) )^2 ;; peculiar velocity correction
  d_L *= 3.0857d24 ;; Mpc -> cm conversion  (since zeropoint spectrum is in erg/cm^2/s/Ang)

  flux = lumdens / (4.*!dpi*d_L^2*(1.0+z)) 
  ;; (1+z) above is d_lambda(obs) = (1+z)*d_lambda(emit)

  if host_ebv ne 0. then flux *= ccm( wave_emit, host_ebv, host_Rv )
  if MW_ebv ne 0. then flux *= ccm( wave_obs, MW_ebv, MW_Rv )

  ifw = (_izmag_struct.ACS_F775W)->wavelength()
  ift = (_izmag_struct.ACS_F775W)->throughput()
  zfw = (_izmag_struct.ACS_F850LP)->wavelength()
  zft = (_izmag_struct.ACS_F850LP)->throughput()

  iw = where( wave_obs ge min(ifw) and wave_obs le max(ifw) )
  iwave_obs = wave_obs[iw]
  iflux = flux[iw]
  zw = where( wave_obs ge min(zfw) and wave_obs le max(zfw) )
  zwave_obs = wave_obs[zw]
  zflux = flux[zw]

  ift_interp = interpol( ift, ifw, iwave_obs )
  zft_interp = interpol( zft, zfw, zwave_obs )

  ABflux = _izmag_struct.ABspec->flux( frame='obs', unit='flambda' )
  ABwave = _izmag_struct.ABspec->wavelength( frame='obs' )
  iABflux = interpol( ABflux, ABwave, iwave_obs )
  zABflux = interpol( ABflux, ABwave, zwave_obs )

  idwave = iwave_obs - shift( iwave_obs, 1 )
  idwave[0] = idwave[1]
  zdwave = zwave_obs - shift( zwave_obs, 1 )
  zdwave[0] = zdwave[1]
  
  zsame = zwave_obs*zft_interp*zdwave
  znum = total(zflux*zsame)
  zden = total(zABflux*zsame)

  isame = iwave_obs*ift_interp*idwave
  inum = total(iflux*isame)
  iden = total(iABflux*isame)

  return, -2.5*alog10([inum/iden,znum/zden])
end
