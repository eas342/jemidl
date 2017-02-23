Function jem_cosmic1, im
  sz = size(im, /dim)
  mask = byte(im*0)
  ;; loop over all pixels
  for i=0, sz[0]-1 do begin
     for j=0, sz[1]-1 do begin
        v1 = im[i,j,*] ;;time evolution of pixel
        ngood = total(finite(v1))
        case ngood of 
           0: begin
              mask[i,j,*] = 1
           end
           1: begin
              w = where(~finite(v1) or v1 lt -1. or v1 gt 10.)
              if w[0] ne -1 then mask[i,j,w] = 1
           end
           2: begin
              w = where(finite(v1), complement=wnot)
              dif = abs(v1[w[0]] - v1[w[1]])
              if dif lt 0.01 then begin
                 if wnot[0] ne -1 then mask[i,j,wnot] = 1
              endif else begin
                 if wnot[0] ne -1 then mask[i,j,wnot] = 1
                 mx = max(v1, imx, /nan)
                 mask[i,j,imx] = 1
              endelse
           end
           else: begin
              mx = max(v1, imx, /nan)
              v2 = rm_elt(v1, imx)
              if ngood gt 4 then begin
                 mx2 = max(v2, imx2, /nan)
                 v2 = rm_elt(v2, imx2)
              endif
              mmt = moment( v2, /nan )
              s2 = sqrt(mmt[1])
              mn2 = mmt[0]
              w = where( ~finite(v1) or v1 lt -1. or v1 gt 10. or v1 gt 6.*s2+mn2 )
              if w[0] ne -1 then mask[i,j,w]=1
           end
        end
     endfor
  endfor
  return, mask
end
