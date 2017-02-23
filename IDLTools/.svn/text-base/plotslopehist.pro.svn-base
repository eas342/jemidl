Pro plotslopehist, $
   xhist, $
   yhist, $
   slope, $
   yhistogram=yhistogram, $
   fill=fill, $
   xhist2=xhist2, $
   yhist2=yhist2, $
   noplot=noplot, $
   _extra=extra

  xhist1 = xhist+slope*yhist
  xrange = minmax([xhist1,xhist])

  xhist2 = dblarr(n_elements(xhist)*2-1)
  yhist2 = dblarr(n_elements(yhist)*2-1)
  if n_elements(xhist) eq 1 then begin
     xhist2[0]=xhist[0]
     yhist2[0]=yhist[0]
  endif else begin
     for i=0, n_elements(xhist)-1 do begin
        if i eq 0 then begin
           xhist2[0]=xhist[0]
           yhist2[0]=yhist[0]
           yhist2[1]=yhist[0]
        endif else if i eq n_elements(xhist)-1 then begin
           xhist2[2*i-1]=xhist[i]
           xhist2[2*i]=xhist[i]
           yhist2[2*i]=yhist[i]
        endif else begin
           xhist2[2*i-1]=xhist[i]
           xhist2[2*i]=xhist[i]
           yhist2[2*i]=yhist[i]
           yhist2[2*i+1]=yhist[i]
        endelse
     endfor
  endelse
  ;;now add final 2 points and first point
  xhist2=[xhist2[0],xhist2]
  yhist2=[0,yhist2]

  xhist2=[xhist2,xhist[1]-xhist[0]+xhist[n_elements(xhist)-1]]
  yhist2=[yhist2,yhist2[n_elements(yhist2)-1]]

  xhist2=[xhist2,xhist2[n_elements(xhist2)-1]]
  yhist2=[yhist2,0]

  xhist2 += yhist2*slope

  if keyword_set(yhistogram) then begin
     tmp = xhist2
     xhist2 = yhist2
     yhist2 = tmp
  endif
  if ~keyword_set(noplot) then begin
     if keyword_set(fill) then begin
        polyfill, xhist2, yhist2, _extra=extra
     endif else begin
        oplot, xhist2, yhist2, ps=0, _extra=extra
     endelse
  endif
end
