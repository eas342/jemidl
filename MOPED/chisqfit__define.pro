Function ChiSqFit::Init, $
   modelfn, $
   nparams, $
   ivar, $
   data, $
   normalize=normalize, $
   verbose=verbose, $
   plot=plot, $
   win=win

  self.normalize=keyword_set(normalize)
  self.verbose=keyword_set(verbose)
  self.plot=keyword_set(plot)
  if n_elements(win) eq 0 then self.win=-1 else self.win=win
  self.pixwin=-1

  self.modelfn=modelfn
  self.nparams=nparams
  self.ivar=ptr_new(ivar)
  self.data=ptr_new(data)
  
  if self.plot then self->wininit

  return, 1
end

Pro ChiSqFit::WinInit
  if self.win eq -1 then begin
     window, /free, xsize=500, ysize=300
     self.win = !d.window
  endif
  if self.pixwin eq -1 then begin
     window, /free, /pixmap, xsize=500, ysize=300
     self.pixwin = !d.window
  endif
  yrange = percentile(*self.data,[0.02,0.98])
  yrange[0] = yrange[0] < 0.
  yrange += [-0.1,0.3]*(yrange[1]-yrange[0])
  wset, self.pixwin
  plot, *self.data, xstyle=1, ystyle=1, yrange=yrange
end

Function ChiSqFit::Likelihood, P
  mu = call_function(self.modelfn, P)
  if self.normalize then begin
     junk = computechi2(*self.data, sqrt(*self.ivar), mu, acoeff=acoeff)
     mu *= acoeff[0]
  endif
  like = 0.5*total((mu-*self.data)^2*(*self.ivar))
  if self.verbose then print, P, like
  if self.plot then begin
     wset, self.win
     device, copy = [0,0,500,300,0,0,self.pixwin]
     oplot, mu, color=fsc_color('red', !d.table_size-2)
     for i=0, self.nparams-1 do $
        xyouts, 0.1, 0.9-0.02*i, string(P[i]), /norm
  endif
  return, like
end

Function ChiSqFit::Amoeba, tol, P0, scale, _extra=extra
  return, amoebaobj(tol, P0=P0, scale=scale, function_name='Likelihood', object=self, _extra=extra)
end

Function ChiSqFit::Grid, $
   gridpts, $ ; an nparam array of ptrs to arrays of gridpoints to search
   like=like ; the output likelihood map

  dims = intarr(self.nparams)
  for i=0, self.nparams-1 do dims[i] = n_elements(*gridpts[i])
  like = make_array(dims, /double)
  for i=0, n_elements(like)-1 do begin
     P = dblarr(self.nparams)
     indices = array_indices(like, i)
     for j=0, self.nparams-1 do P[j] = (*gridpts[j])[indices[j]]
     like[i] = self->Likelihood(P)
  endfor
  junk = min(like, m)
  index = array_indices(like, m)
  for j=0, self.nparams-1 do P[j] = (*gridpts[j])[index[j]]
  return, P
end

Pro ChiSqFit::Cleanup
  ptr_free, self.ivar, self.data
  if self.pixwin ne -1 then wdelete, self.pixwin
end

Pro ChiSqFit__Define, struct
  struct = { ChiSqFit, $
             Inherits JEMobject, $
             modelfn: "", $
             nparams: 0L, $
             ivar: ptr_new(), $
             data: ptr_new(), $
             normalize: 0, $
             verbose: 0, $
             plot: 0, $
             win: -1, $
             pixwin: -1 }
end
