Function synphot, $
   wave_emit, $                 ;; rest(emmited) wavelengths in Angstroms
   lum_dens, $                  ;; luminosity density of source (erg/s/Ang)
   z_velocity, $                ;; redshift for velocity correction
   z_distance=z_distance, $     ;; use this redshift for luminosity distance
   l_distance=l_distance, $     ;; or just use this luminosity distance (Mpc)
   host_ebv=host_ebv, $         ;; E(B-V) for host in rest-frame
   host_Rv=host_Rv, $           ;; total-to-selective extinction coefficient for host
   MW_ebv=MW_ebv, $             ;; E(B-V) for MW in obs-frame
   MW_Rv=MW_Rv, $               ;; total-to-selective extinction coefficient for MW
   f_name=f_name, $             ;; name of bandpass filter
   f_wave=f_wave, $             ;; alternatively, wavelengths of filter
   f_throughput=f_throughput, $ ;; and throughput of filter
   zp_name=zp_name, $           ;; zeropoint spectrum : ['Vega', 'AB', 'ST']
   zp_wave=zp_wave, $           ;; alternatively, zeropoint spectrum wavelength (Ang)
   zp_flux=zp_flux              ;; alternatively, zeropoint spectrum flux (erg/s/cm^2/Ang)

  if size(lum_dens, /n_dim) eq 1 then nspec=1 $
  else nspec=(size(lum_dens, /dim))[1]
  nwave=(size(lum_dens, /dim))[0]

  if n_elements(z_velocity) eq 0 then z_velocity=0.
  if n_elements(z_distance) eq 0 then z_distance=z_velocity
  if n_elements(host_ebv) eq 0 then host_ebv=0.
  if n_elements(host_Rv) eq 0 then host_Rv = 3.1
  if n_elements(MW_ebv) eq 0 then MW_ebv=0.
  if n_elements(MW_Rv) eq 0 then MW_Rv=3.1
  if n_elements(l_distance) eq 0 then begin
     if z_distance eq 0.0 then l_distance = 1.0e-5 $ ;; 10pc in units of Mpc
     else l_distance = lumdist( z_distance, H0=70.2, omega_m=0.272, lambda0=0.728, /silent )
  endif
  if ~finite(z_velocity) then return, !values.f_nan

  if n_elements(f_wave) eq 0 then begin
     f = filter(f_name)
     f_wave = f.wave
     f_throughput = f.throughput
  endif

  if n_elements(zp_wave) eq 0 then begin
     zp = zpspec(zp_name)
     zp_wave = zp.wave
     zp_flux = zp.flux
  endif

  wave_obs = wave_emit * (1.0 + z_velocity)

  ;; compute effective luminosity distance, correcting for peculiar velocity
  l_distance_cm = l_distance * 3.0857d24 ;; Mpc -> cm conversion (since zeropoint spectra are in erg/cm^2/s/Ang)
  l_distance_cm_eff = l_distance_cm * ( (1.0 + z_velocity)/(1.0 + z_distance) )^2

  flux = lum_dens / ( 4.0d*!dpi*l_distance_cm_eff^2 * (1.0 + z_velocity) )

  ;; apply dust
  if host_ebv ne 0.0 then begin
     if nspec eq 1 then flux *= ccm( wave_emit, host_ebv, host_Rv ) $
     else $
        flux *= rebin( ccm( wave_emit, host_ebv, host_Rv ), nwave, nspec )
  endif
  if MW_ebv ne 0.0 then begin
     if nspec eq 1 then flux *= ccm( wave_obs, MW_ebv, MW_Rv ) $
     else $
        flux *= rebin( ccm( wave_obs, host_ebv, host_Rv ), nwave, nspec )
  endif

  ;; truncate object spectrum to range of filter
  w = where( wave_obs ge min(f_wave) and wave_obs le max(f_wave) )
  wave_obs_trunc = wave_obs[w]
  if nspec eq 1 then flux_trunc = flux[w] $
  else flux_trunc = flux[w,*]

  ;; interpolate filter throughput and zeropoint spectrum
  f_throughput_interp = interpol( f_throughput, f_wave, wave_obs_trunc )
  zp_spec_interp = interpol( zp_flux, zp_wave, wave_obs_trunc )

  ;; need d_lambda
  dwave_obs_trunc = wave_obs_trunc - shift( wave_obs_trunc, 1 )
  dwave_obs_trunc[0] = dwave_obs_trunc[1]

  ;; integrate and return result
  same = dwave_obs_trunc*wave_obs_trunc*f_throughput_interp
  denominator = total( zp_spec_interp*same )
  if nspec eq 1 then begin
     numerator = total( flux_trunc*same )
  endif else begin
     same = rebin( same, n_elements(same), nspec )
     numerator = total( flux_trunc*same, 1 )
  endelse
  return, -2.5*alog10( numerator / denominator )
end
