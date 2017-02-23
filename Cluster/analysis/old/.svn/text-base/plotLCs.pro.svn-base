Pro plotLCs, LCcat, band=band
  if n_elements(band) eq 0 then band='z'
  for igal=0, n_elements(LCcat)-1 do begin
     soplot, LCcat[igal].zmag, color=(igal mod 7)+1, linestyle=0 ;solid
     soplot, LCcat[igal].zmag_PSFmatch, color=(igal mod 7)+1, linestyle=1;dotted
     soplot, LCcat[igal].zmag_CLEAN, color=(igal mod 7)+1, linestyle=2;dashed
  endfor
end
