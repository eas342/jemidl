; This function returns an IDL structure containing cleaned sextractor ACS i/z catalogs.

Function import_sexcat, ifile, zfile, clusterid
  fmt = 'I,F,F,D,D,F,F,F,F,F,F,F,F,F,F,D,I,F,F,F,F,F,F,F,F,F' ;through magerr_petro
  for i=0,31 do fmt += ',F,F'   ;take care of aperture mags....
  print, 'reading F775W catalog'
  myreadcol, ifile, format=fmt, /silent, $
             galnum, xwin_image, ywin_image, alphawin_j2000, deltawin_j2000, $
             iflux_radius, iclass_star, a_image, b_image, $
             theta_image, ellipticity, kron_radius, petro_radius, $
             imu_max, imu_threshold, ibackground, $
             iisoarea_image, ifwhm_image, $
             imag_best, imagerr_best, $
             imag_auto, imagerr_auto, $
             imag_iso, imagerr_iso, $
             imag_isocor, imagerr_isocor, $
             imag_petro, imagerr_petro, $
             iap1,  iap2,  iap3,  iap4, $
             iap5,  iap6,  iap7,  iap8, $
             iap9,  iap10, iap11, iap12, $
             iap13, iap14, iap15, iap16, $
             iap17, iap18, iap19, iap20, $
             iap21, iap22, iap23, iap24, $
             iap25, iap26, iap27, iap28, $
             iap29, iap30, iap31, iap32, $
             iaperr1,  iaperr2,  iaperr3,  iaperr4, $
             iaperr5,  iaperr6,  iaperr7,  iaperr8, $
             iaperr9,  iaperr10, iaperr11, iaperr12, $
             iaperr13, iaperr14, iaperr15, iaperr16, $
             iaperr17, iaperr18, iaperr19, iaperr20, $
             iaperr21, iaperr22, iaperr23, iaperr24, $
             iaperr25, iaperr26, iaperr27, iaperr28, $
             iaperr29, iaperr30, iaperr31, iaperr32
  iap=[[iap1],  [iap2],  [iap3],  [iap4], $
       [iap5],  [iap6],  [iap7],  [iap8], $
       [iap9],  [iap10], [iap11], [iap12], $
       [iap13], [iap14], [iap15], [iap16], $
       [iap17], [iap18], [iap19], [iap20], $
       [iap21], [iap22], [iap23], [iap24], $
       [iap25], [iap26], [iap27], [iap28], $
       [iap29], [iap30], [iap31], [iap32] ]
  iaperr=[[iaperr1],  [iaperr2],  [iaperr3],  [iaperr4], $
          [iaperr5],  [iaperr6],  [iaperr7],  [iaperr8], $
          [iaperr9],  [iaperr10], [iaperr11], [iaperr12], $
          [iaperr13], [iaperr14], [iaperr15], [iaperr16], $
          [iaperr17], [iaperr18], [iaperr19], [iaperr20], $
          [iaperr21], [iaperr22], [iaperr23], [iaperr24], $
          [iaperr25], [iaperr26], [iaperr27], [iaperr28], $
          [iaperr29], [iaperr30], [iaperr31], [iaperr32] ]
  print, 'reading F850LP catalog'
  myreadcol, zfile, format=fmt, /silent, $
             galnum, xwin_image, ywin_image, alphawin_j2000z, deltawin_j2000z, $
             zflux_radius, zclass_star, a_image, b_image, $
             theta_image, ellipticity, kron_radius, petro_radius, $
             zmu_max, zmu_threshold, zbackground, $
             zisoarea_image, zfwhm_image, $
             zmag_best, zmagerr_best, $
             zmag_auto, zmagerr_auto, $
             zmag_iso, zmagerr_iso, $
             zmag_isocor, zmagerr_isocor, $
             zmag_petro, zmagerr_petro, $
             zap1,  zap2,  zap3,  zap4, $
             zap5,  zap6,  zap7,  zap8, $
             zap9,  zap10, zap11, zap12, $
             zap13, zap14, zap15, zap16, $
             zap17, zap18, zap19, zap20, $
             zap21, zap22, zap23, zap24, $
             zap25, zap26, zap27, zap28, $
             zap29, zap30, zap31, zap32, $
             zaperr1,  zaperr2,  zaperr3,  zaperr4, $
             zaperr5,  zaperr6,  zaperr7,  zaperr8, $
             zaperr9,  zaperr10, zaperr11, zaperr12, $
             zaperr13, zaperr14, zaperr15, zaperr16, $
             zaperr17, zaperr18, zaperr19, zaperr20, $
             zaperr21, zaperr22, zaperr23, zaperr24, $
             zaperr25, zaperr26, zaperr27, zaperr28, $
             zaperr29, zaperr30, zaperr31, zaperr32
  zap=[[zap1],  [zap2],  [zap3],  [zap4], $
       [zap5],  [zap6],  [zap7],  [zap8], $
       [zap9],  [zap10], [zap11], [zap12], $
       [zap13], [zap14], [zap15], [zap16], $
       [zap17], [zap18], [zap19], [zap20], $
       [zap21], [zap22], [zap23], [zap24], $
       [zap25], [zap26], [zap27], [zap28], $
       [zap29], [zap30], [zap31], [zap32] ]
  zaperr=[[zaperr1],  [zaperr2],  [zaperr3],  [zaperr4], $
          [zaperr5],  [zaperr6],  [zaperr7],  [zaperr8], $
          [zaperr9],  [zaperr10], [zaperr11], [zaperr12], $
          [zaperr13], [zaperr14], [zaperr15], [zaperr16], $
          [zaperr17], [zaperr18], [zaperr19], [zaperr20], $
          [zaperr21], [zaperr22], [zaperr23], [zaperr24], $
          [zaperr25], [zaperr26], [zaperr27], [zaperr28], $
          [zaperr29], [zaperr30], [zaperr31], [zaperr32] ]

  for i=0, n_elements(galnum)-1 do begin
     acscat1 = {sexcat, $
                clusterid:clusterid, $
                galid:galnum[i], $
                xwin_image:xwin_image[i], $
                ywin_image:ywin_image[i], $
                alphawin_j2000:alphawin_j2000[i], $
                deltawin_j2000:deltawin_j2000[i], $
                iflux_radius:iflux_radius[i], $
                iclass_star:iclass_star[i], $
                zflux_radius:zflux_radius[i], $
                zclass_star:zclass_star[i], $
                a_image:a_image[i], $
                b_image:b_image[i], $
                theta_image:theta_image[i], $
                ellipticity:ellipticity[i], $
                kron_radius:kron_radius[i], $
                petro_radius:petro_radius[i], $
                imu_max:imu_max[i], $
                imu_threshold:imu_threshold[i], $
                ibackground:ibackground[i], $
                zmu_max:zmu_max[i], $
                zmu_threshold:zmu_threshold[i], $
                zbackground:zbackground[i], $

                iisoarea_image:iisoarea_image[i], $
                zisoarea_image:zisoarea_image[i], $
                ifwhm_image:ifwhm_image[i], $
                zfwhm_image:zfwhm_image[i], $

                imag_best:imag_best[i], $
                imagerr_best:imagerr_best[i], $
                zmag_best:zmag_best[i], $
                zmagerr_best:zmagerr_best[i], $
                
                imag_auto:imag_auto[i], $
                imagerr_auto:imagerr_auto[i], $
                zmag_auto:zmag_auto[i], $
                zmagerr_auto:zmagerr_auto[i], $
                
                imag_iso:imag_iso[i], $
                imagerr_iso:imagerr_iso[i], $
                zmag_iso:zmag_iso[i], $
                zmagerr_iso:zmagerr_iso[i], $
                
                imag_isocor:imag_isocor[i], $
                imagerr_isocor:imagerr_isocor[i], $
                zmag_isocor:zmag_isocor[i], $
                zmagerr_isocor:zmagerr_isocor[i], $
                
                imag_petro:imag_petro[i], $
                imagerr_petro:imagerr_petro[i], $
                zmag_petro:zmag_petro[i], $
                zmagerr_petro:zmagerr_petro[i], $
                
                imag_aper:reform(iap[i,*], 32), $
                imagerr_aper:reform(iaperr[i,*], 32), $
                zmag_aper:reform(zap[i,*], 32), $
                zmagerr_aper:reform(zaperr[i,*], 32) }
     if n_elements(acscat) eq 0 then acscat=acscat1 else $
        acscat=[acscat, acscat1]
  endfor
  return, acscat
end
