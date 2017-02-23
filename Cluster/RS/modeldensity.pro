Function denslikenobkg, data, CMRparams
  sigpart = CMRparams.sigma^2 + data.csigma^2 + CMRparams.slope*data.msigma^2
  logpart = -0.5*alog(2.*!pi)-0.5*alog(sigpart)
  constpart = -(data.color - CMRparams.slope*(data.mag - 23) - CMRparams.intercept23)^2
  constpart /= 2.*sigpart
  return, logpart + constpart
end

Function denslikebkg, data, bkgparams
  return, bkgparams.d + bkgparams.e*data.mag + bkgparams.f*data.color + bkgparams.g*data.mag^2
end

Function densfieldintegrand, mag, color
  common JEM$_modeldensity, bkgparams, CMRparams, Schechterparams, iskysigma, zskysigma, colors, icolor
  return, denslikebkg({mag:mag, color:color}, bkgparams)
end

Function densclusterintegrand, mag, color
  common JEM$_modeldensity
  return, exp(denslikenobkg( {mag:mag, msigma:0., color:color, csigma:0.03}, $
                             CMRparams ))*SchechterM( mag, $
                                                      Schechterparams.phi_star, $
                                                      Schechterparams.M_star, $
                                                      Schechterparams.alpha )
end

Function col_limits, x
  common JEM$_modeldensity
  return, [colors[icolor], colors[icolor+1]]
end

Function modeldensity, CMRparams1, bkgparams1, Schechterparams1
  common JEM$_modeldensity
  CMRparams=CMRparams1
  bkgparams=bkgparams1
  Schechterparams=Schechterparams1
  mags = findgen(18)/17*6+19.
  colors = findgen(20)/19.*2-0.5
  dens = fltarr( n_elements(mags)-1, n_elements(colors)-1 )
;  colors = findgen(50)/49*2.-0.5
;  mags = findgen(100)/99.*6.+19.
;  dens = fltarr( n_elements(mags)-1, n_elements(colors)-1 )
  for icolor=0, n_elements(colors)-2 do begin
     for imag=0, n_elements(mags)-2 do begin
        dens[imag,icolor] = int_2d( 'densfieldintegrand', [mags[imag], mags[imag+1]], 'col_limits', 6 )
        dens[imag,icolor] += int_2d( 'densclusterintegrand', [mags[imag], mags[imag+1]], 'col_limits', 6 )
     endfor
  endfor
  return, dens
end
