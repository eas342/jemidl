Function ldist_f, z, w0=w0, wa=wa
  return, (1d + z)^(3d * (1d + w0 + wa))*exp(-3d * wa * z / (1d + z))
end

Function ldist2, z, omega_m=omega_m, omega_l=omega_l, w0=w0, wa=wa
  return, 1d/sqrt((omega_m*(1d + z)^3 + omega_l*ldist_f(z, w0=w0, wa=wa)))
end

Function lumdist2, z, h=h, omega_m=omega_m, omega_l=omega_l, w0=w0, wa=wa
  if n_elements(h) eq 0 then h=0.7
  if n_elements(omega_m) eq 0 then omega_m=0.3
  if n_elements(omega_l) eq 0 then omega_l=0.7
  if n_elements(w0) eq 0 then w0=-1.0
  if n_elements(wa) eq 0 then wa=0.0

  cLight = 299792.458d ;; km/s
  h0 = 100*h;; km/s/Mpc
  if omega_l EQ 0 then begin
     denom = sqrt(1+omega_m*z) + 1 + omega_m*0.5*z
     dlum = (cLight*z/h0)*(1 + z*(1-omega_m*0.5)/denom)
     return,dlum
  endif

  if n_elements(z) eq 1 then begin
     qsimp, 'ldist2', 0, z, lz, omega_m=omega_m, omega_l=omega_l, w0=w0, wa=wa
     return, cLight*(1d + z)/h0*lz
  endif else begin
     lz=dblarr(n_elements(z))
     for i=0, n_elements(z)-1 do begin
        qsimp, 'ldist2', 0d, z[i], lz1, omega_m=omega_m, omega_l=omega_l, w0=w0, wa=wa
        lz[i] = cLight*(1d + z[i])/h0*lz1
     endfor
     return, lz
  endelse
end
