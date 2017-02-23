Function Filter::Init, $
   wavelength, $
   throughput

  if n_elements(wavelength) ne n_elements(throughput) then $
     message, 'Wavelength and throughput arrays must have the same length!'
  
  sort = sort(wavelength)
  self.wavelength = Ptr_New(wavelength[sort])
  self.throughput = Ptr_New(throughput[sort])
  return, 1
end

Function Filter::Interpolate, $
   wavelength

  throughput = dblarr( n_elements(wavelength) )
  lower = where( wavelength le min(*self.wavelength) )
  middle = where( wavelength gt min(*self.wavelength) and wavelength lt max(*self.wavelength) )
  upper = where( wavelength ge max(*self.wavelength) )

  if lower[0] ne -1 then throughput[lower] = 0.d
  if middle[0] ne -1 then throughput[middle] = interpol(*self.throughput, *self.wavelength, wavelength[middle], /spline)
  if upper[0] ne -1 then throughput[upper] = 0.d
  return, throughput

end

Function Filter::Wavelength
  return, *self.wavelength
end

Function Filter::Throughput, $
   wavelength=wavelength

  if n_elements(wavelength) ne 0 then $
     throughput = self->interpolate(wavelength) $
  else $
     throughput = *self.throughput
  return, throughput
end

Pro Filter::CleanUp
  ptr_free, self.wavelength
  ptr_free, self.throughput
end

Pro Filter__Define
  struct = { Filter, $
             wavelength: Ptr_New(), $
             throughput: Ptr_New() }
end
