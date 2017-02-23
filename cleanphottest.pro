Function CLEANPhotTest, image, psf, radius, loopgain=loopgain
  if n_elements(loopgain) eq 0 then loopgain=0.003
  coord=size(image,/dim)/2
  clean=obj_new('CLEAN',image,psf,1.5,loopgain=loopgain)
  iter=-1
  mag=0
  mag0=apphot(image, coord[0], coord[1], radius)
  for i=0,99 do begin
     counter, i+1, 100
     void=clean->Clean(10000L)
     mag = clean->phot(coord, radius)
     stats = clean->stats()
     if n_elements(out) eq 0 then out = struct_addtags(stats, {mag:mag}) $
     else out = [out, struct_addtags(stats, {mag:mag})]
  endfor
  obj_destroy, clean
  return, out
end
