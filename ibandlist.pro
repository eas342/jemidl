Function ibandlist, image_list, zlist=zlist
  for i=0, n_elements(image_list)-1 do begin
     if strmid(sxpar(headfits(image_list[i]),'FILTER1'),0,5) eq 'F775W' then begin
        if n_elements(outlist) eq 0 then outlist = image_list[i] $
        else outlist=[outlist,image_list[i]]
     endif else begin
        if n_elements(zlist) eq 0 then zlist = image_list[i] $
        else zlist=[zlist,image_list[i]]
     endelse
  endfor
  return, outlist
end
