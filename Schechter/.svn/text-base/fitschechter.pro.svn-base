Pro FitSchechter, binbounds, ngals, cs=cs, alphas=alphas, m_stars=m_stars
  if n_elements(binbounds)-1 ne n_elements(ngals) then begin
     message, 'binbounds-1 != ngals'
     return
  endif
  alphas = findgen(25)*0.05-1.25
  m_stars = findgen(25)*0.1+19.5

  cs=dblarr( n_elements(alphas), n_elements(m_stars) )
  for ialpha=0, n_elements(alphas)-1 do begin
     for im_star=0, n_elements(m_stars)-1 do begin
        alpha=alphas[ialpha]
        m_star=m_stars[im_star]
        ngal_star = schechtermbin( 1., m_star, alpha, $
                                   min(binbounds), max(binbounds) )
        phi_star = total(ngals)/ngal_star
        c=0
        for ibin=0, n_elements(ngals)-1 do begin
           ni=ngals[ibin]
           mi = schechtermbin( phi_star, m_star, alpha, $
                               binbounds[ibin], binbounds[ibin+1] )
           c += 2*(ni*alog(mi) - mi - alog(factorial(ni)))
        endfor
        cs[ialpha, im_star] = -c
     endfor
  endfor
  
end
