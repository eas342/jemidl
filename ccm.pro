Function ccm_ir, x, rv
  a0 = 0.574
  b0 = -0.527
  c0 = rv * a0 + b0
  return, c0 * ( x^1.61 )
end

Function ccm_optical, x, rv
  y = x - 1.82
  ;; a coeffs
  a = [ 1.0, 0.17699, -0.50447, -0.02427, $
        0.72085, 0.01979, -0.77530, 0.32999 ]
  ;; b coeffs
  b = [ 0.0, 1.41338, 2.28305, 1.07233, $
        -5.38434, -0.62251, 5.30260, -2.09002 ]
  ;; total coeffs
  c = rv * a + b

  ;; finally result
  return, poly(y, c)
end

Function ccm_uv1, x, rv
  y  = x - 5.90
  fa = (-0.04473*y^2 - 0.009779*y^3) * (x ge 5.90)
  fb = ( 0.21300*y^2 + 0.120700*y^3) * (x ge 5.90)
  aa = 1.752 - 0.316 * x - 0.104 / ( ( x - 4.67 )^2 + 0.341 ) + fa
  bb = -3.09 + 1.825 * x + 1.206 / ( ( x - 4.62 )^2 + 0.263 ) + fb
  return, rv * aa + bb
end

Function ccm_uv2, x, rv
  y = x - 8.0
  ; a coeffs
  a = [ -1.073, -0.628, 0.137, -0.070]
  ; b coeffs
  b = [ 13.670, 4.257, -0.420, 0.374]
  ; total coeffs
  c = rv * a + b
  return, poly( y, c )
end

Function ccm, wavelength, ebv, rv
  if n_elements(rv) eq 0 then rv=3.1
  x = 1.0 / (wavelength * 1.e-4) ;;Angstroms to inverse microns
  ex1 = ccm_ir( x, rv )      * ( x lt 1.1 )
  ex2 = ccm_optical( x, rv ) * ( x lt 3.3 and x ge 1.1 )
  ex3 = ccm_uv1( x, rv )     * ( x lt 8.0 and x ge 3.3 )
  ex4 = ccm_uv2( x, rv )     * ( x ge 8.0 )
  extinction = ex1 + ex2 + ex3 + ex4
  return, 10.0^(-extinction * ebv / 2.5 )
end
