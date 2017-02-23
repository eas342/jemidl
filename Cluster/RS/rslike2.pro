Function likenobkg, data, CMRparams
  sigpart = CMRparams.sigma^2 + data.csigma^2 + CMRparams.slope*data.msigma^2
  constpart = 1./sqrt(2.*!pi*sigpart)
  exppart = -(data.color - CMRparams.slope*(data.mag - 23) - CMRparams.intercept23)^2
  exppart /= 2.*sigpart
  return, constpart * exp(exppart)
end

Function likebkg, data, bkgparams
  return,  bkgparams.d $
           + bkgparams.e*data.mag $
           + bkgparams.f*data.color $
           + bkgparams.g*data.mag^2
end

Function pq_limits, x
  common JEM$_rslike, magrange, colorrange, $
     bkgparams, CMRparams, Schechterparams, $
     clusterdata, bkgdata, $
     iskysigma, zskysigma, $
     clusterarea, bkgarea ;;colorrange
  return, [colorrange[0], colorrange[1]]
end

;;scalar mag, vector color
Function fieldintegrand, mag, color
  common JEM$_rslike ;; bkgparams
  return, likebkg( {mag:replicate(mag, n_elements(color)), $
                    color:color}, $
                   bkgparams )
end

Function colorerror, zmag, color
  common JEM$_rslike
  common JEM$_colorerror, ifits, zfits, radii
  if n_elements(ifits) eq 0 then $
     restore, '/home/scpdata02/joshimages/skystats/sky2apsig.sav'
  ifluxerr = ifits[8,0]+ifits[8,1]*iskysigma
  zfluxerr = zfits[8,0]+zfits[8,1]*zskysigma
  imag = zmag+color
  iflux = 10^(-0.4*(imag-25.67849))
  zflux = 10^(-0.4*(zmag-24.86663))
  iexptime=3000
  zexptime=8000
  ifluxerr = sqrt(ifluxerr^2+iflux/iexptime)
  zfluxerr = sqrt(zfluxerr^2+zflux/zexptime)
  idfof = ifluxerr/iflux
  zdfof = zfluxerr/zflux
  ierr = 2.5/2*alog10((1+idfof)/(1-idfof))
  zerr = 2.5/2*alog10((1+zdfof)/(1-zdfof))
  return, sqrt(ierr^2+zerr^2)
end

Function clusterintegrand, mag, color
  common JEM$_rslike
  csigma = colorerror( mag, color )
  return, likenobkg( {mag:mag, msigma:0., color:color, csigma:csigma}, $
                     CMRparams )*SchechterM( mag, $
                                             Schechterparams.phi_star, $
                                             Schechterparams.M_star, $
                                             Schechterparams.alpha )
end

Function loglikebkg, abscissa
  common JEM$_rslike ;;data, magrange
  d=abscissa[0]
  e=abscissa[1]
  f=abscissa[2]
  g=abscissa[3]
  bkgparams = { d:d, e:e, f:f, g:g }
  integral = int_2d( 'fieldintegrand', magrange, 'pq_limits', 6 )
  out = total(alog(likebkg( bkgdata, bkgparams )))-integral
  print, integral
 return, -out
end

Function loglike, abscissa
  common JEM$_rslike ;;clusterdata, bkgdata
  CMRparams = { intercept23:abscissa[0], $
                slope:abscissa[1], $
                sigma:abscissa[2] }
  phi_star=abscissa[3]
  m_star=abscissa[4]
  alpha=abscissa[5]
  Schechterparams = { phi_star:phi_star, $
                      m_star:m_star, $
                      alpha:alpha }
  bkgparams = { d:abscissa[6], $
                e:abscissa[7], $
                f:abscissa[8], $
                g:abscissa[9] }
  
  bkgintegral = int_2d( 'fieldintegrand', magrange, 'pq_limits', 6 )
  clusterintegral = int_2d( 'clusterintegrand', magrange, 'pq_limits', 6 )
  clusterout = total(alog(likebkg( clusterdata, bkgparams ) $
                          +likenobkg( clusterdata, CMRparams ) $
                          *SchechterM( clusterdata.mag, phi_star, m_star, alpha ) )) $
               - bkgintegral - clusterintegral + alog(clusterarea)
  bkgout = total(alog(likebkg( bkgdata, bkgparams ))) - bkgintegral + alog(bkgarea)
  print, -(clusterout+bkgout)
  return, -(clusterout+bkgout)
end

Pro rslike2, clustercmd, start=start, scale=scale
  common JEM$_rslike
  if n_elements(start) eq 0 then start = [0.946, -0.028, 0.028, 10., 22., -1., -35.5, 3.4, -15.3, -0.036]
  if n_elements(scale) eq 0 then scale = [0.01, 0.001, 0.001, 1., 0.3, 0.005, 0.2, 0.01, 0.2, 0.002]
  if n_elements(clustercmd) eq 0 then begin
     restore, '/home/scpdata02/cluster3/Y/Y_2.sav'
     resolve_obj, obj
     clustercmd = cmd( obj, /xpsf, /nofit )
     s=obj->summary()
     iskysigma = mean( s.isky_sigma, /nan )
     zskysigma = mean( s.zsky_sigma, /nan )
     clusterarea = 2.56
     bkgarea = 1.69*30
  endif

  magrange = [20., 24.]
  colorrange = [0.6, 1.2]

  cgals = *clustercmd.gals
  clusterdata = {mag:0., color:0., msigma:0., csigma:0.}
  clusterdata = replicate( clusterdata, n_elements(cgals) )
  clusterdata.mag = cgals.zmag
  clusterdata.color = cgals.iz
  clusterdata.msigma = cgals.zmagerr
  clusterdata.csigma = cgals.iz_err
  wcluster = where( clusterdata.mag gt magrange[0] $
                    and clusterdata.mag le magrange[1] $
                    and clusterdata.color gt colorrange[0] $
                    and clusterdata.color le colorrange[1] )
  clusterdata = clusterdata[wcluster]

  restore, '/home/scpdata03/goods/cmd/cmds.sav'
  bkggals = *(*cmds[0]).gals
  for i=1, n_elements(cmds)-1 do begin
     bkggals=[bkggals, *(*cmds[i]).gals]
  endfor
  bkgdata = {mag:0., color:0.}
  bkgdata = replicate( bkgdata, n_elements(bkggals) )
  bkgdata.mag = bkggals.zmag
  bkgdata.color = bkggals.iz
  wbkg = where( bkgdata.mag gt magrange[0] $
                and bkgdata.mag le magrange[1] $
                and bkgdata.color gt colorrange[0] $
                and bkgdata.color le colorrange[1] )
  bkgdata = bkgdata[wbkg]

  result = downhillsimplex( 'loglike', start, scale, $
                            tol=0.00001, fval=fval, miniter=100 )
  
  print, "intercept23:", result[0]
  print, "slope:      ", result[1]
  print, "sigma:      ", result[2]
  print, "phi_star:   ", result[3]
  print, "m_star:     ", result[4]
  print, "alpha:      ", result[5]
  print, "d:          ", result[6]
  print, "e:          ", result[7]
  print, "f:          ", result[8]
  print, "g:          ", result[9]
end
