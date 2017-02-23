Function mode, data, res=res, sig=sig, minmax=minmax
  if n_elements(res) eq 0 then res=0.1
  if n_elements(sig) eq 0 then sig=1.0
  if n_elements(minmax) eq 0 then minmax=minmax(data)

  center = minmax[0]+res
  statmax = -1.
  while center lt minmax[1] do begin
     w = where(abs(data-center) lt 4.*sig)
     stat = total(exp(-((data[w]-center)/sig)^2))
     if stat gt statmax then begin
        statmax = stat
        mode = center
     endif
     center += res
  endwhile
  return, mode
end
