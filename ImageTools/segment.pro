;; Segments the 1-bit image 'image' into 8-contiguous regions.

Function segment, image
  s=size(image,/dim)
  nx=s[0]
  ny=s[1]
  isegment = 1
  segmap = intarr(s)
  identical = [[0,0]]
  ;;do first row separately
  ;;do first pixel separately
  if image[0,0] then begin
     segmap[0,0] = isegment
     isegment += 1
  endif
  for ix=1, nx-1 do begin
     if ~image[ix,0] then continue
     if ~image[ix-1,0] then begin
        segmap[ix,0] = isegment
        isegment += 1
     endif else segmap[ix,0]=segmap[ix-1,0]
  endfor
  for iy=1, ny-1 do begin
     ;;do first column separately
     if image[0,iy] then begin
        if segmap[0,iy-1] ne 0 then begin
           segmap[0,iy]=segmap[0,iy-1]
        endif else if segmap[1,iy-1] ne 0 then begin
           segmap[0,iy]=segmap[1,iy-1]
        endif else begin
           segmap[0,iy] = isegment
           isegment += 1
        endelse
     endif
     ;;do middle columns
     for ix=1, nx-2 do begin
        if ~image[ix,iy] then continue
        if segmap[ix,iy-1] ne 0 then begin ;;check lower
           segmap[ix,iy]=segmap[ix,iy-1]
           continue
        endif
        if segmap[ix-1,iy-1] ne 0 then begin ;;check lower-left
           if segmap[ix+1,iy-1] ne 0 $
              and segmap[ix+1,iy-1] ne segmap[ix-1,iy-1] then begin
              identical=[[identical],[[segmap[ix+1,iy-1],segmap[ix-1,iy-1]]]]
           endif
           segmap[ix,iy]=segmap[ix-1,iy-1]
           continue
        endif
        if segmap[ix-1,iy] ne 0 then begin ;;check left
           if segmap[ix+1,iy-1] ne 0 $
              and segmap[ix+1,iy-1] ne segmap[ix-1,iy] then begin
              identical=[[identical],[[segmap[ix+1,iy-1],segmap[ix-1,iy]]]]
           endif
           segmap[ix,iy]=segmap[ix-1,iy]
           continue
        endif
        if segmap[ix+1,iy-1] ne 0 then begin ;;check lower-right
           segmap[ix,iy]=segmap[ix+1,iy-1]
           continue
        endif
        segmap[ix,iy] = isegment
        isegment += 1
     endfor
     ;;do final column
     if image[nx-1, iy] then begin
        if segmap[nx-1,iy-1] ne 0 then begin ;;check lower
           segmap[nx-1,iy]=segmap[nx-1,iy-1]
           continue
        endif
        if segmap[nx-2,iy-1] ne 0 then begin ;;check lower-left
           segmap[nx-1,iy]=segmap[nx-2,iy-1]
           continue
        endif
        if segmap[nx-2,iy] ne 0 then begin ;;check left
           segmap[nx-1,iy]=segmap[nx-2,iy]
           continue
        endif
        segmap[nx-1,iy]=isegment
        isegment += 1
     endif
  endfor

  if size(identical, /n_dim) ne 1 then begin
     nident = (size(identical, /dim))[1]
     for iident=1, nident-1 do begin
        w=where(segmap eq identical[0,iident])
        if w[0] ne -1 then begin
           segmap[w] = identical[1,iident]
           w1 = where(identical eq identical[0,iident])
           if w1[0] ne -1 then identical[w1] = identical[1,iident]
        endif
     endfor
     newident=1
     for iident=1, max(segmap) do begin
        w=where(segmap eq iident)
        if w[0] ne -1 then begin
           segmap[w] = newident
           newident += 1
        endif
     endfor
  endif
  return, segmap
end
