pro testprospect
  resolve_routine, 'spectrum__define', /com
  resolve_routine, 'filter__define', /com
  resolve_routine, 'dopplershift__define', /com
  restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/filters.sav'
  restore, '/home/jmeyers314/olivetos/jemidl/ProSpect/data/zeropointspecs.sav'

  ivegazp = obj_new( 'zeropoint', vegaspec, filters.ACS_F775W )
  zvegazp = obj_new( 'zeropoint', vegaspec, filters.ACS_F850LP )
  iABzp = obj_new( 'zeropoint', ABspec, filters.ACS_F775W )
  zABzp = obj_new( 'zeropoint', ABspec, filters.ACS_F850LP )
  iSTzp = obj_new( 'zeropoint', STspec, filters.ACS_F775W )
  zSTzp = obj_new( 'zeropoint', STspec, filters.ACS_F850LP )
  print, 'These should be the same'
  print, -(ivegazp->base_mag() - zvegazp->base_mag())
  print, 25.291 - 24.347
  print

  print, -(iABzp->base_mag() - zABzp->base_mag())
  print, 25.678 - 24.867
  print

  print, -(iSTzp->base_mag() - zSTzp->base_mag())
  print, 26.417 - 25.955
end
