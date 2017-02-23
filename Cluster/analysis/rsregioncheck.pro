Function RSRegionCheck, zmag_auto, iz, magrange, colorrange, bluelimit, slope=slope
  if n_elements(slope) eq 0 then slope=0.
  rsregioncheck = zmag_auto ge magrange[0] $
                  and zmag_auto le magrange[1] $
                  and iz ge (colorrange[0] + (zmag_auto-23)*slope) $
                  and iz le (colorrange[1] + (zmag_auto-23)*slope) $
                  and iz ge bluelimit
  return, rsregioncheck
end
