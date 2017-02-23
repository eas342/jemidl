Function Mpc2arcmin, rMpc, z
  rarcmin = rMpc / (lumdist( z, /silent )/(1.+z)^2) * 180./!pi * 60.
  return, rarcmin
end
