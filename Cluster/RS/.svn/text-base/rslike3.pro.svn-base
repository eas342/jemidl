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

Function fieldintegrand, mag, color
  common JEM$_rslike, bkgparams, CMRparams, Schechterparams, $
     iskysigma, zskysigma
  return, likebkg({mag:mag, color:color}, bkgparams)
end

Function colorerror, zmag, color
  common JEM$_rslike
  common JEM$_rslike2, ifits, zfits, radii
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

Function pq_limits, x
  return, [0.7, 1.1]
end

Function rslike, clusterdata, bkgdata, $
                 clusterarea, bkgarea, $
                 magrange, colorrange, $
                 iskysigma1, zskysigma1, $
                 CMRparams1, bkgparams1, Schechterparams1
  common JEM$_rslike
  CMRparams=CMRparams1
  bkgparams=bkgparams1
  Schechterparams=Schechterparams1
  iskysigma=iskysigma1
  zskysigma=zskysigma1
  wcluster = where( clusterdata.mag gt magrange[0] $
                    and clusterdata.mag lt magrange[1] $
                    and clusterdata.color gt colorrange[0] $
                    and clusterdata.color lt colorrange[1] )
  wbkg = where( bkgdata.mag gt magrange[0] $
                    and bkgdata.mag lt magrange[1] $
                    and bkgdata.color gt colorrange[0] $
                    and bkgdata.color lt colorrange[1] )

  likenobkgpart = likenobkg( clusterdata[wcluster], CMRparams )
  Schechterpart = SchechterM( clusterdata[wcluster].mag, $
                              Schechterparams.phi_star, $
                              Schechterparams.M_star, $
                              Schechterparams.alpha )
  clusterlikebkgpart = likebkg( clusterdata[wcluster], bkgparams )
  bkglikebkgpart = likebkg( bkgdata[wbkg], bkgparams )
  bkgintegral = int_2d( 'fieldintegrand', magrange, 'pq_limits', 6 )
  clusterintegral = int_2d( 'clusterintegrand', magrange, 'pq_limits', 6 )
  return, total(alog(likenobkgpart*Schechterpart+clusterlikebkgpart))+alog(clusterarea) $
          +total(alog(bkglikebkgpart))+alog(bkgarea)-bkgintegral*(clusterarea+bkgarea) $
          -clusterintegral*clusterarea
end

