Pro ZeroPointSpecs_sav
  sf=obj_new('spectrumfactory')
  vegafile = '/home/scpdata01/CALSPEC/ascii_2009-01/alpha_lyr_stis_004.ascii'
  vegaspec = sf->createFromFile( vegafile )
  wave = vegaspec->wavelength()
  STflux = wave*0+3.63d-9
  STspec = obj_new( 'spectrum',wave, STflux, unit='flambda' )
  c_ang_s=2.99792d18
  freq = c_ang_s/wave
  ABflux_nu = freq*0+3.63d-20
  ABflux_lambda = ABflux_nu*c_ang_s/wave^2
  ABspec = obj_new( 'spectrum', wave, ABflux_lambda, unit='flambda' )

  
  save, vegaspec, STspec, ABspec, filename='zeropointspecs.sav'
end
