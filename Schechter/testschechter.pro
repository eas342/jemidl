Pro test, seed
  window, 0, xsize=500, ysize=200
  window, 1, xsize=500, ysize=200
  ntrials = 1000
  ngal=55
  m_star=21.
  alpha=-0.75
  binbounds = findgen(9)*0.5+20
  ngals = intarr(n_elements(binbounds)-1)
  ngals_fiducial = intarr(n_elements(binbounds)-1)

  ngal_star = schechtermbin( 1., m_star, alpha, $
                             min(binbounds), max(binbounds) )
  phi_star = ngal/ngal_star

  fitalphas = dblarr(ntrials)
  fitm_stars = dblarr(ntrials)
  for ibin=0, n_elements(ngals)-1 do begin
     ngals_fiducial[ibin] = schechtermbin( phi_star, m_star, alpha, $
                                           binbounds[ibin], binbounds[ibin+1] )
  endfor

  for itrial=0, ntrials-1 do begin
     for ibin=0, n_elements(ngals)-1 do begin
        ngals[ibin] = randomn( seed, poisson=ngals_fiducial[ibin] )
     endfor
     fitschechter, binbounds, ngals, cs=cs, alphas=alphas, m_stars=m_stars
     junk = min(cs, m)
     fitalphas[itrial] = alphas[m mod n_elements(alphas)]
     fitm_stars[itrial] = m_stars[m / n_elements(alphas)]
     if itrial gt 3 then begin
        wset, 0
        plothist, fitalphas[0:itrial], bin=0.1, xrange=[-1.3, 0.], xstyle=1
        wset, 1
        plothist, fitm_stars[0:itrial], bin=0.2, xrange=[18., 22.5], xstyle=1
     endif
  endfor
end
