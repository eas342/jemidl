Function RSlike::Init, $
   _extra=extra

  self->SetProperty, _extra=extra
  return, 1
end

Function RSlike::likenobkg, data, P
  CMRint23 = P[0]
  CMRslope = P[1]
  CMRsigma = P[2]

  sigpart = CMRsigma^2 + data.csigma^2 + CMRslope^2*data.msigma^2
  constpart = 1./sqrt(2.*!pi*sigpart)
  exppart = -(data.color - CMRslope*(data.mag-23.) - CMRint23)^2
  exppart /= 2.*sigpart
  
  return, constpart * exp(exppart)
end

Function RSlike::likebkg, data, P
  d = P[6]
  e = P[7]
  f = P[8]
  g = P[9]
  return, d + e*data.mag $
          + f*data.color $
          + g*data.mag^2
end

Function RSlike::BkgIntegrand, mag, color
  return, self.bkgarea*self->likebkg( { mag:replicate(mag, n_elements(color)), $
                                        color:color }, self.param )
end

Function RSlike::colorerror, zmag, color
  common JEM$_RSlike_colorerror, ifits, zfits
  if n_elements(ifits) eq 0 then $
     restore, '/home/scpdata02/joshimages/skystats/sky2apsig.sav'
  ifluxerr = ifits[8,0]+ifits[8,1]*self.iskysigma
  zfluxerr = zfits[8,0]+zfits[8,1]*self.zskysigma
  imag = zmag+color
  iflux = 10^(-0.4*(imag-25.67849))
  zflux = 10^(-0.4*(zmag-24.86663))
  iexptime=self.iexptime
  zexptime=self.zexptime
  ifluxerr = sqrt(ifluxerr^2+iflux/iexptime)
  zfluxerr = sqrt(zfluxerr^2+zflux/zexptime)
  idfof = ifluxerr/iflux
  zdfof = zfluxerr/zflux
  ierr = 2.5/2*alog10((1+idfof)/(1-idfof))
  zerr = 2.5/2*alog10((1+zdfof)/(1-zdfof))
  return, sqrt(ierr^2+zerr^2)
end

Function RSlike::ClusterIntegrand, mag, color
  csigma = self->colorerror( mag, color )
  phi_star = (self.param)[3]
  m_star = (self.param)[4]
  alpha = (self.param)[5]
  abscissa = { mag:replicate(mag, n_elements(color)), $
               color:color, $
               msigma:replicate(0., n_elements(color)), $
               csigma:csigma }
  return, (self->likenobkg( abscissa, self.param ) $
           * SchechterM( mag, phi_star, m_star, alpha ) $
           + self->likebkg( abscissa, self.param )) $
          * self.clusterarea
end

Function RSlike::PQ_limits, mag
  return, self.colorrange
end

Function RSlike::LogLikelihood, P
  self.param = P
  phi_star = P[3]
  m_star = P[4]
  alpha = P[5]
  bkglike = total( $
            alog( $
            self.bkgarea*self->likebkg( *self.bkgdata, P ) ) )
  clusterlike = total( $
                alog( $
                self.clusterarea*self->likenobkg( *self.clusterdata, P ) $
                * SchechterM( (*self.clusterdata).mag, $
                              phi_star, m_star, alpha ) $
                + self.clusterarea*self->likebkg( *self.clusterdata, P ) ) )
  bkgintegral = jem_int_2d( 'bkgintegrand', self.magrange, $
                            'pq_limits', 20, object=self )
  clusterintegral = jem_int_2d( 'clusterintegrand', self.magrange, $
                                'pq_limits', 20, object=self )
  loglikelihood = bkglike+clusterlike-bkgintegral-clusterintegral
;  print, format='(%"bkglike:          %11.4f")', bkglike
;  print, format='(%"clusterlike:      %11.4f")', clusterlike
;  print, format='(%"bkgintegral:      %11.4f")', bkgintegral
;  print, format='(%"clusterintegral:  %11.4f")', clusterintegral
;  print, format='(%"loglikelihood:  %13.4f")', loglikelihood
  return, -loglikelihood
end

Pro RSlike::CleanUp
  if ptr_valid(self.clusterdata) then ptr_free, self.clusterdata
  if ptr_valid(self.bkgdata) then ptr_free, self.bkgdata
end

Pro RSlike__Define, struct
  struct = { RSlike, $
             Inherits JEMobject, $
             clusterdata:ptr_new(), $
             clusterarea:0., $
             bkgdata:ptr_new(), $
             bkgarea:0., $
             magrange:fltarr(2), $
             colorrange:fltarr(2), $
             iexptime:0., $
             zexptime:0., $
             iskysigma:0., $
             zskysigma:0., $
             param:fltarr(10) $
           }
end
