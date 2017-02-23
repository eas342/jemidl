Function onelickflux, obswave, flux, err, redshift, vdisp, res, index
  speedOfLight = 299792.458d
  ;; load Lick index definitions
  indexlist = '$EZ_AGES_DIR/lickindexlist.txt'
  readcol, indexlist, band1, band2, blue1, blue2, red1, red2, indres, $
           indmag, indexname, c1, c2, $
           FORMAT='F,F,F,F,F,F,F,I,A,F,F', /silent
  nindex = (size(indexname))[1]

  ;; convert to vacuum wavelengths
  airtovac, band1
  airtovac, band2
  airtovac, blue1
  airtovac, blue2
  airtovac, red1
  airtovac, red2

  ;; which index are we measuring
  j = (where(index eq indexname))[0]
  if j eq -1 then return, 0

  restwave = obswave / (1d + redshift)

  wwave = where(restwave ge blue1[j]-50 and restwave le red2[j]+50)

  j_wave = restwave[wwave]
  j_spec = flux[wwave]
  j_err = err[wwave]
  j_res = res[wwave]

  ;; normalize here, remember to multiply result by medspec at the end
  medspec = median(j_spec)
  j_spec /= medspec
  j_err /= medspec

  dlambda = (max(j_wave)-min(j_wave))/(n_elements(j_wave)-1)
  lambda = mean(j_wave)

  sig_ind = indres[j] / sqrt(8.0*alog(2.0)) / dlambda
  sig_gal = (double(vdisp)/speedOfLight)*(lambda/dlambda)

  sig_res = (mean(j_res) / speedOfLight) * (lambda / dlambda) / sqrt(8.0*alog(2.0))
  sig_tot = sqrt(sig_res^2 + sig_gal^2)

  if (sig_ind le sig_tot) then begin
     smoothspec = j_spec
     smooth_err = j_err
  endif else begin
     sigma = sqrt(sig_ind^2 - sig_tot^2)
     filtersize=(fix(sigma*8)+(fix(sigma*8+1) mod 2))
     gaussfilter = psf_gaussian(npixel=filtersize, $
                                st_dev=sigma, ndimen=1, /normalize)
     smoothspec = convol(j_spec, gaussfilter, /center, $
                         /edge_truncate)

     smooth_err = convol((j_err)^2, gaussfilter, /center, $
                         /edge_truncate)
     smooth_err = sqrt(smooth_err)
  endelse

  blue = where((j_wave gt blue1[j]) and (j_wave lt blue2[j]))
  red = where((j_wave gt red1[j]) and (j_wave lt red2[j]))
  band = where((j_wave gt band1[j]) and (j_wave lt band2[j]))

  lambda_b = (blue1[j] + blue2[j]) / 2.0
  lambda_r = ( red1[j] +  red2[j]) / 2.0

  theta_blue = j_wave[blue] - j_wave[blue-1]
  theta_red  = j_wave[ red] - j_wave[ red-1]
  theta_band = j_wave[band] - j_wave[band-1]

  bandspec = smoothspec[band]

  ;;-----------------------------------------------
  ;; Include fractional pixels at end of bandpass
  ;;-----------------------------------------------

  f_b1 = (j_wave[Min(blue)] - blue1[j]) /               $
         (j_wave[Min(blue)] - j_wave[Min(blue)-1])
  f_b2 = (blue2[j] - j_wave[Max(blue)]) /             $
         (j_wave[Max(blue)+1] - j_wave[Max(blue)])
  last_b = (Size(blue))[1]-1

  f_r1 = (j_wave[Min(red)] - red1[j]) /               $
         (j_wave[Min(red)] - j_wave[Min(red)-1])
  f_r2 = (red2[j] - j_wave[Max(red)]) /             $
         (j_wave[Max(red)+1] - j_wave[Max(red)])
  last_r = (Size(red))[1]-1

  ;;-----------------------------------------------------
  ;; Find "average" point in blue and red pseudocontinua
  ;;-----------------------------------------------------

  s_blue = (Int_Tabulated(j_wave[blue], smoothspec[blue])      +  $
            f_b1 * smoothspec[Min(blue)-1] * theta_blue[0]     +  $
            f_b2 * smoothspec[Max(blue)] * theta_blue[last_b]) /  $
           (blue2[j] - blue1[j])
  s_red  = (Int_Tabulated(j_wave[red] , smoothspec[red])       +  $
            f_r1 * smoothspec[Min(red)-1] * theta_red[0]      +  $
            f_r2 * smoothspec[Max(red)] * theta_red[last_r] ) /  $
           (red2[j]  -  red1[j])

  ;;--------------------------------------------------
  ;; Make pseudocontinuum.
  ;;--------------------------------------------------

  pseudoslope = (s_red-s_blue) / $
                (lambda_r-lambda_b)
  pseudointercept = s_blue-pseudoslope*lambda_b

  pseudocont = pseudoslope*j_wave[band] + pseudointercept
  wholecont  = pseudoslope*j_wave + pseudointercept

  ;;--------------------------------
  ;; Measure Lick Index
  ;;-------------------------------

  f_band1 = (j_wave[Min(band)] - band1[j]) /               $
            (j_wave[Min(band)] - j_wave[Min(band)-1])
  f_band2 = (band2[j] - j_wave[Max(band)]) /             $
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

  flux = Int_Tabulated(j_wave[band], bandspec-pseudocont)
  blue_fracpix = 0.5 * f_band1 * theta_band[0] $
                 * (smoothspec[Min(band)-1]-wholecont[Min(band)-1] $
                    + smoothspec[Min(band)]-wholecont[Min(band)])
  red_fracpix = 0.5 * f_band2 * theta_band[last_band] $
                * (smoothspec[Max(band)+1]-wholecont[Max(band)+1] $
                   + smoothspec[Max(band)]-wholecont[Max(band)])
  flux += blue_fracpix + red_fracpix

  var_flux = total(smooth_err[band]^2 * theta_band^2)
  blue_fracvar = 0.5 * f_band1 * theta_band[0]^2 $
                 * (smooth_err[Min(band)-1]^2 $
                    + smooth_err[Min(band)]^2)
  red_fracvar = 0.5 * f_band2 * theta_band[last_band]^2 $
                * (smooth_err[Max(band)+1]^2 $
                    + smooth_err[Max(band)]^2)

  var_flux += blue_fracvar + red_fracvar

  ;; plot, j_wave, j_spec, title=index, $
  ;;       xrange=[min(j_wave[blue])-50, max(j_wave[red])+50], $
  ;;       yrange=[min(j_spec)-0.15*(max(j_spec)-min(j_spec)), $
  ;;               max(j_spec)], ps=10
  ;; oplot, [red1[j], red2[j]], [s_red, s_red], $
  ;;        thick=2, color=ctable.green
  ;; oplot, [blue1[j], blue2[j]], [s_blue, s_blue], $
  ;;        thick=2, color=ctable.green
  ;; oplot, j_wave, smoothspec, thick=2, color=ctable.red
  ;; oplot, [lambda_b, lambda_r], [s_blue, s_red], $
  ;;        linestyle=2, color=ctable.blue
  ;; oplot, j_wave[band], pseudocont, thick=4, color=ctable.blue
  ;; polyfill, [j_wave[band[0]], j_wave[band], j_wave[reverse(band)]], $
  ;;           [pseudocont[0], smoothspec[band], reverse(pseudocont)], $
  ;;           /line_fill, color=ctable.gold

  return, {flux:flux*medspec, var_flux:var_flux*medspec^2}
end
