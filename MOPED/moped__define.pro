;; MOPED class
;; The MOPED algorithm is describe in Heavens et al. MNRAS 2000
;; IDL code written by JEMeyers LBL 06/02/08

;; Moped initializer.  
Function MOPED::Init, $
   modelfn, $  ;;IDL function of one variable modeling data
   dmodelfn, $  ;;derivatives of above
   nparams, $ ;;number of model parameters
   ivar, $ ;;inverse variance vector of data
   fiducialparams, $ ;;fiducial parameters where derivatives are computed
   data, $ ;;the data to be fitted
   limits=limits, $ ;;valid ranges for variables
   normalize=normalize, $ ;;optional normalization
   verbose=verbose, $ ;;term output
   plot=plot, $
   win=win;;graphical output

  self.normalize=Keyword_Set(normalize)
  self.verbose=Keyword_Set(verbose)
  self.plot=Keyword_Set(plot)
  if n_elements(win) eq 0 then self.win=-1 else self.win=win
  self.pixwin=-1
  if n_elements(limits) ne 0 then self.limits = ptr_new(limits)

  self.modelfn = modelfn
  self.dmodelfn = dmodelfn
  self.nparams = nparams
  self.ivar = Ptr_New(ivar)
  self.data = Ptr_New(data)
  self.fiducialparams = Ptr_New(fiducialparams)

  self->FiducialInit
  return, 1
end

;; You must initialize MOPED at a fiducial point to compute derivatives
;; for weight vectors.
Pro MOPED::FiducialInit
  derivs = call_function(self.dmodelfn, *self.fiducialparams)
  
  n = dblarr(self.nparams) ;; normalizations
  b = dblarr(n_elements(*self.ivar), self.nparams) ; weight vectors
  ybar = dblarr(self.nparams) ;; compressed data; to be compared to compressed model

  n[0] = sqrt(total((derivs[*,0])*(*self.ivar)*(derivs[*,0])))
;  n[0] = sqrt(derivs[*,0]#*self.cinv#derivs[*,0])
  b[*,0] = (*self.ivar*derivs[*,0])/n[0]
;  b[*,0] = (*self.cinv#derivs[*,0])/n[0]
  ybar[0] = (b[*,0]##transpose(*self.data))[0]
  for i=1, self.nparams-1 do begin
     sum0=0
     sum1=0
     for j=0, i do begin
        term = ((transpose(derivs[*,i])#b[*,j])[0])
        sum0 += term^2
        sum1 += term*b[*,j]
     endfor
     n[i] = sqrt(total((derivs[*,i])*(*self.ivar)*(derivs[*,i])) - sum0)
;     n[i] = sqrt(derivs[*,i]#*self.cinv#derivs[*,i] - sum0)
     b[*,i] = (*self.ivar*derivs[*,i] - sum1)/n[i]
;     b[*,i] = (*self.cinv#derivs[*,i] - sum1)/n[i]
     ybar[i] = (b[*,i]##transpose(*self.data))[0]
  endfor
  self.b = Ptr_New(b)
  self.ybar = Ptr_New(ybar)

  if self.plot then begin
     if self.win eq -1 then begin
        window, /free, xsize=500, ysize=300
        self.win = !d.window
     endif
     if self.pixwin eq -1 then begin
        window, /free, /pixmap, xsize=500, ysize=300
        self.pixwin = !d.window
     endif
     xrange=[0,max(where(b[*,0] ne 0))]
     yrange=percentile(*self.data, [0.02, 0.98])
     yrange[0] = yrange[0] < 0.
     yrange += [-0.1,0.3]*(yrange[1]-yrange[0])
     wset, self.pixwin
     plot, *self.data, xstyle=1, ystyle=1, xrange=xrange, yrange=yrange
  endif
end

;;  For debugging
Pro MOPED::SpillGuts
  print, self.modelfn
  print, self.dmodelfn
  print, self.nparams
  help, *self.ivar
  help, *self.fiducialparams
  help, *self.b
end

;; Returns the log likelihood of the model given parameters P
Function MOPED::Likelihood, P
  if ptr_valid(self.limits) then begin
     for i=0, n_elements(p)-1 do begin
        if p[i] lt (*self.limits)[0,i] or p[i] gt (*self.limits)[1,i] then return, !values.f_infinity
     endfor
  endif
  mu = call_function(self.modelfn, P)
  if self.normalize then begin
     junk = computechi2(*self.data, sqrt(*self.ivar), mu, acoeff=acoeff)
     mu *= acoeff[0]
  endif
  
  y = dblarr(self.nparams)
  for i=0, self.nparams-1 do $
     y[i] = ((*self.b)[*,i]##transpose(mu))[0]
  like = 0.5*total((y - *self.ybar)^2)
  if self.verbose then print, format='(%"%18.7e  %18.7e  %18.7e")',P, like
  if self.plot then begin
     wset, self.win
     device, copy = [0,0,500,300,0,0,self.pixwin]
     oplot, mu, color=fsc_color('red', !d.table_size-2)
     for i=0, self.nparams-1 do $
        xyouts, 0.1, 0.9-0.02*i, string(P[i]), /norm
  endif
  return, like
end

Function MOPED::b
  return, *self.b
end

;; Downhill simplex minimization of the log-likelihood function
Function MOPED::Amoeba, tol, P0, scale, _extra=extra
  return, amoebaobj(tol, P0=P0, scale=scale, function_name='Likelihood', object=self, _extra=extra)
end

;; Josh's downhill simplex algorithm...
Function MOPED::Downhillsimplex, tol, P0, scale, _extra=extra
  return, downhillsimplex( 'Likelihood', P0, scale, object=self, _extra=extra )
end

;; Grid search minimization of the log-likelihood function
Function MOPED::Grid, $
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

;; Moped destructor
Pro MOPED::Cleanup
  Ptr_Free, self.ivar
  Ptr_Free, self.fiducialparams
  Ptr_Free, self.data
  Ptr_Free, self.b
  Ptr_Free, self.ybar
  ptr_free, self.limits
  if self.pixwin ne -1 then wdelete, self.pixwin
end

;; Moped class definition
Pro MOPED__Define, struct
  struct = { MOPED, $
             Inherits JEMobject, $
             modelfn: "", $  ; name of function that returns model
             dmodelfn: "", $ ; name of function that returns dmodel/dparam for 
             nparams: 0L, $ ; n different parameters
             ivar: Ptr_New(), $ ; inverse covariance matrix to use
             fiducialparams: Ptr_New(), $ ; fiducial guess at which derivatives are evaluated
             data: Ptr_New(), $ ; data to be fit
             b: Ptr_New(), $ ; weight vectors
             ybar: Ptr_New(), $ ; compressed data to be fit
             limits: Ptr_New(), $
             normalize: 0, $ ; normalization flag
             verbose: 0, $ ; term output
             plot: 0, $
             win: -1, $
             pixwin: -1}         ; window output
end
