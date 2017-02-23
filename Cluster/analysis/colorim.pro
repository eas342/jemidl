Pro colorim, obj, filehead
  resolve_obj, obj
;  obj->LoadImage
  R=*(obj->extract('image'))[0]
  mwrfits, R, filehead+'R.fits', /create
  delvarx, R
  G=0.5*(*(obj->extract('image'))[0]+*(obj->extract('image'))[1])
  mwrfits, G, filehead+'G.fits', /create
  delvarx, G
  B=*(obj->extract('image'))[1]
  mwrfits, B, filehead+'B.fits', /create
  delvarx, B
  spawn, 'stiff '+filehead+'R.fits '+filehead+'G.fits '+filehead+'B.fits -c stiff.conf'
  spawn, 'mv stiff.tif '+filehead+'.tif'
  spawn, 'rm -rf '+filehead+'R.fits '+filehead+'G.fits '+filehead+'B.fits'
end
