Function warpSky::Init, $
   vacwave=vacwave, $
   skyflux=skyflux, $
   skyMask=skyMask, $
   _extra=extra
  self->SetProperty, _extra=extra
  self.vacwave = ptr_new(vacwave)
  self.skyflux = ptr_new(skyflux)
  self.skyMask = ptr_new(skyMask)
  a = mrdfits('/Users/josh/skylines/sky.fwhm.3.0.fits',1,/silent)
  self.modelVacWave = ptr_new(a.vacwave)
  self.modelSkyFlux = ptr_new(a.flux)
  return, 1
end

Function warpSky::warpedSky, param
  waveRange = param.waveRange
  waveCoeffs = param.waveCoeffs
  fluxCoeffs = param.fluxCoeffs
  waveOrder = n_elements(waveCoeffs)-1
  fluxOrder = n_elements(fluxCoeffs)-1
  x = (*self.vacWave-waveRange[0])/(waveRange[1]-waveRange[0])
  waveOut = *self.vacWave - (waveCoeffs##fchebyshev(x, waveOrder+1))
  warpFactor = abs(fluxCoeffs##fchebyshev(x, fluxOrder+1))
  fluxOut = *self.skyFlux*warpFactor
  return, {vacWave:waveOut, flux:fluxOut}
end

Function warpSky::modelDeviates, param
  s = self->warpedSky( param )
  modelComparison = interpol(*self.modelSkyFlux, *self.modelVacWave, s.vacWave)
  return, (s.flux - modelComparison)*(*self.skyMask)
end

Function warpSky::model
  return, {vacWave:*self.modelVacWave, flux:*self.modelSkyFlux}
end

Pro warpSky__define, struct
  struct = { warpSky, $
             Inherits JEMobject, $
             vacWave:ptr_new(), $
             skyFlux:ptr_new(), $
             skyMask:ptr_new(), $
             modelVacWave:ptr_new(), $
             modelSkyFlux:ptr_new() $
           }
end
