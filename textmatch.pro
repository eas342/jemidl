Pro textmatch, list1, list2, match1, match2
  n1 = n_elements(list1)
  n2 = n_elements(list2)
  maxmatch = min([n1, n2])
  match1 = intarr(maxmatch)
  match2 = intarr(maxmatch)
  imatch = 0
  for i=0, n_elements(list1)-1 do begin
     w=where(list2 eq list1[i])
     if w[0] ne -1 then begin
        match1[imatch] = i
        match2[imatch] = w
        imatch += 1
     endif
  endfor
  if imatch eq 0 then begin
     match1 = -1
     match2 = -1
  endif else begin
     match1 = match1[0:imatch-1]
     match2 = match2[0:imatch-1]
  endelse
end
