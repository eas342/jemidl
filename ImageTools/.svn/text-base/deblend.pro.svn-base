Function deblend, image1, analysisthresh
  acskern = [[0.0322, 0.0718, 0.0322], $
             [0.0718, 0.1600, 0.0718], $
             [0.0322, 0.0718, 0.0322]]
  acskern /= total(acskern)
  image = convol( image1, acskern, /edge_truncate )
  background = image le analysisthresh
  
  npeaks = 30
  mm=max(image)
  ncontours = 64
  lmax = alog10(mm)
  lmin = lmax - 2.
  lcontours = findgen(64)/63*(lmax-lmin)+lmin
  contours = 10^lcontours

  sz = size( image, /dim )
  levels=bytarr( sz[0], sz[1] )
  for i=0, ncontours-1 do begin
     levels += image gt contours[i]
  endfor

  seg=image*0
  mask = seg eq 0 and ~background
  for i=0, npeaks-1 do begin
     max = max(image*mask, m)
     level = levels[m]
     seg1 = segment( levels eq level)
     obj = seg1 eq seg1[m]
     level -= 1
     while level ge 0 do begin
        seg2 = segment( (levels eq level) or obj )
        if total(seg2 eq seg2[m] and seg) gt 0 then break
        if level lt analysisthresh then break
        obj = seg2 eq seg2[m]
        level -= 1
     endwhile
     seg += (i+1)*obj
     atv, seg mod 10, min=0, max=9
     mask = seg eq 0 and ~background
  endfor
  return, seg-background
end
