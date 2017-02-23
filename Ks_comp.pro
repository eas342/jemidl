Function Ks_comp, age, metal, logmass
  sun = mrdfits('/Users/josh/reference/CALSPEC/fits_2010-02/sun_reference_stis_001.fits', 1, /silent)
  sunwave=sun.wavelength
  sunlumdens=sun.flux*(1.5e13*1.5e13*4*!pi)
  Ks_sun = synphot(sunwave, sunlumdens, f_name='KS_2mass', zp_name='AB')
  spec = sspspec(age, metal, /full)
  Ks_spec = synphot(spec.wave, spec.flux*10^(logmass)/spec.ms[0], f_name='KS_2mass', zp_name='AB')
  return, -0.4*(Ks_spec - Ks_sun)
end
