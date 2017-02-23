function isnot_telluric, wave,aggressive=aggressive
  ; if want to use aggressive definition of telluric features
  if (keyword_set(aggressive)) then begin
      return, wave lt 6272 or $
        (wave gt 6301 and wave lt 6863) or $
        (wave gt 6915 and wave lt 7162) or $
        (wave gt 7340 and wave lt 7590) or $
        (wave gt 7700 and wave lt 8123) or $
        (wave gt 8362 and wave lt 8926) or $
        (wave gt 9218 and wave lt 9267) or $
        (wave gt 9674)
  endif else begin
      ;else don't use features that make sensfunc ratty
      return, wave lt 6272 or $
        (wave gt 6301 and wave lt 6863) or $
        (wave gt 6915 and wave lt 7590) or $
        (wave gt 7700 and wave lt 9267) or $
        (wave gt 9674)
  endelse
end
