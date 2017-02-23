Function onelickflux2, obswave, flux, fluxfit, err, redshift, $
                       band1, band2
  speedOfLight = 299792.458d
  ;; load Lick index definitions

  restwave = obswave / (1d + redshift)
  if min(restwave) gt band1 or max(restwave) lt band2 then $
     return, {flux:!values.d_nan, $
              flux_err:!values.d_nan, $
              continuum:!values.d_nan, $
              EW:!values.d_nan, $
              EW_err:!values.d_nan}

  wwave = where(restwave ge band1-50 and restwave le band2+50)

  j_wave = restwave[wwave]
  j_spec = flux[wwave]
  j_fit = fluxfit[wwave]
  j_err = err[wwave]

  ;; normalize here, remember to multiply result by medspec at the end
  medspec = median(j_spec)
  j_spec /= medspec
  j_fit /= medspec
  j_err /= medspec

  dlambda = (max(j_wave)-min(j_wave))/(n_elements(j_wave)-1)
  lambda = mean(j_wave)

  smoothspec = j_spec
  smooth_err = j_err

  band = where((j_wave gt band1) and (j_wave lt band2))

  theta_band = j_wave[band] - j_wave[band-1]
  bandspec = smoothspec[band]
  bandfit = j_fit[band]

  ;;-----------------------------------------------
  ;; Include fractional pixels at end of bandpass
  ;;-----------------------------------------------

  ;; f_b1 = (j_wave[Min(blue)] - blue1) /               $
  ;;        (j_wave[Min(blue)] - j_wave[Min(blue)-1])
  ;; f_b2 = (blue2 - j_wave[Max(blue)]) /             $
  ;;        (j_wave[Max(blue)+1] - j_wave[Max(blue)])
  ;; last_b = (Size(blue))[1]-1

  ;; f_r1 = (j_wave[Min(red)] - red1) /               $
  ;;        (j_wave[Min(red)] - j_wave[Min(red)-1])
  ;; f_r2 = (red2 - j_wave[Max(red)]) /             $
  ;;        (j_wave[Max(red)+1] - j_wave[Max(red)])
  ;; last_r = (Size(red))[1]-1

  ;;--------------------------------
  ;; Measure Lick Index
  ;;-------------------------------

  f_band1 = (j_wave[Min(band)] - band1) /               $
            (j_wave[Min(band)] - j_wave[Min(band)-1])
  f_band2 = (band2 - j_wave[Max(band)]) /             $
            (j_wave[Max(band)+1] - j_wave[Max(band)])
  last_band = (Size(band))[1]-1

  ;; Main part of the index
  ;; ew = Int_Tabulated(j_wave[band],(1-bandspec/pseudocont))

  ;; Add in the fractional pixels at the edge
  ;; blue_fracpix = 0.5 * f_band1 * $
  ;;                ((1-smoothspec[Min(band)-1] / wholecont[Min(band)-1]) * $
  ;;                 theta_band[0] + $
  ;;                 (1-smoothspec[Min(band)] / wholecont[Min(band)]) * $
  ;;                 theta_band[0])
  ;; red_fracpix = 0.5 * f_band2 * $
  ;;               ((1-smoothspec[Max(band)+1] / wholecont[Max(band)+1]) * $
  ;;                theta_band[last_band] + $
  ;;                (1-smoothspec[Max(band)] / wholecont[Max(band)]) * $
  ;;                theta_band[last_band])
  ;; ew = ew + blue_fracpix + red_fracpix

  intflux = Int_Tabulated(j_wave[band], bandspec-bandfit)
  blue_fracpix = 0.5 * f_band1 * theta_band[0] $
                 * (smoothspec[Min(band)-1]-j_fit[Min(band)-1] $
                    + smoothspec[Min(band)]-j_fit[Min(band)])
  red_fracpix = 0.5 * f_band2 * theta_band[last_band] $
                * (smoothspec[Max(band)+1]-j_fit[Max(band)+1] $
                   + smoothspec[Max(band)]-j_fit[Max(band)])
  intflux += blue_fracpix + red_fracpix

  var_flux = total(smooth_err[band]^2 * theta_band^2)
  blue_fracvar = 0.5 * f_band1 * theta_band[0]^2 $
                 * (smooth_err[Min(band)-1]^2 $
                    + smooth_err[Min(band)]^2)
  red_fracvar = 0.5 * f_band2 * theta_band[last_band]^2 $
                * (smooth_err[Max(band)+1]^2 $
                    + smooth_err[Max(band)]^2)

  var_flux += blue_fracvar + red_fracvar

  continuum = mean(bandfit)
  return, { flux:intflux*medspec, $
            flux_err:sqrt(var_flux)*medspec, $
            continuum:continuum*medspec, $
            EW:-intflux/continuum, $
            EW_err:sqrt(var_flux)/continuum }
end
