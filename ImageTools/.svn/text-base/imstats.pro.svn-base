Pro ImStats, image, p=p
  if n_elements(p) ne 0 then begin
     im=image[p[0]:p[1],p[2]:p[3]]
     mom = moment(im)
  endif else mom = moment(image)
  print, string(format='(%"av:%f  stdev:%f")', mom[0], sqrt(mom[1]))
end
