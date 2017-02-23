Function ACSsegmask, image, thresh
  acskern = [[0.0322, 0.0718, 0.0322], $
             [0.0718, 0.1600, 0.0718], $
             [0.0322, 0.0718, 0.0322]]
  acskern /= total(acskern)
  detect = convol( image, acskern, /edge_truncate )
  return, segment( detect gt thresh )
end
