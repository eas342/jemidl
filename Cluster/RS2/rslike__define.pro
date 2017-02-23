Function RSlike::Init, $
   _extra=extra

  self->SetProperty, _extra=extra
  return, 1
end

Function RSlike::Like_RS, data, P
  int23 = P[0]
  sigma = P[1]

  sigpart = sigma^2 + data.izerr^2
  constpart = 1./sqrt(2*!pi*sigpart)
  exppart = -(data.rsresid - int23)^2
  exppart /= 2.*sigpart
  return, constpart * exp(exppart)
end

Function RSlike::Like_Field, data, P
  return, poly(data.rsresid-1., P[3:*]) > 1.e-5
end

Function RSlike::BkgIntegrand, rsresid
  return, self->like_field( {rsresid:rsresid}, *self.param )*self.bkgarea
end

Function RSlike::ClusterIntegrand, rsresid
  N = (*self.param)[2]/self.clusterarea
  return, self->like_RS( {rsresid:rsresid, izerr:self.izerr}, *self.param )*N*self.clusterarea $
          + self->like_field( {rsresid:rsresid}, *self.param )*self.clusterarea
end

Function RSlike::SpecialClusterIntegrand, rsresid
  return, self->like_field( {rsresid:rsresid}, *self.param )*self.clusterarea
end

Function RSlike::RSIntegrand, rsresid
  N = (*self.param)[2]/self.clusterarea
  return, self->like_RS( {rsresid:rsresid, izerr:self.izerr}, *self.param )*N*self.clusterarea
end

Function RSlike::LogLikelihood, P, verbose=verbose
  ptr_free, self.param
  self.param = ptr_new(P)
  N = P[2]/self.clusterarea

  bkglike = total( $
            alog( $
            self.bkgarea*self->like_field(*self.bkgdata, P ) ), /nan )

  clusterlike = total( $
                alog( $
                self.clusterarea * ( $
                self->like_rs( *self.clusterdata, P )*N $
                + self->like_field( *self.clusterdata, P ) ) ) )
  bkgexpect = qpint1d( 'BkgIntegrand', object=self, $
                       self.rsresidrange[0], self.rsresidrange[1] )
  if P[1] lt 0.01 then begin
     clusterexpect = qpint1d( 'SpecialClusterIntegrand', object=self, $
                              self.rsresidrange[0], self.rsresidrange[1] ) $
                     + N*self.clusterarea
  endif else begin
     clusterexpect = qpint1d( 'ClusterIntegrand', object=self, $
                              self.rsresidrange[0], self.rsresidrange[1] )
  endelse
  loglikelihood = bkglike+clusterlike-bkgexpect-clusterexpect

  if keyword_set(verbose) then begin
     print, format='(%"bkglike:          %11.4f")', bkglike
     print, format='(%"clusterlike:      %11.4f")', clusterlike
     print, format='(%"bkgexpect:        %11.4f")', bkgexpect
     print, format='(%"clusterexpect:    %11.4f")', clusterexpect
     print, format='(%"loglikelihood:  %13.4f")', loglikelihood
  endif

  return, -loglikelihood
end

Function RSlike::LogLikelihoodNoBkg, P, verbose=verbose
  ptr_free, self.param
  self.param = ptr_new(P)
  N = P[2]/self.clusterarea
  
  clusterlike = total( $
                alog( $
                self.clusterarea * ( $
                self->like_rs( *self.clusterdata, P )*N ) ) )
  if P[1] lt 0.01 then begin
     clusterexpect = N*self.clusterarea
  endif else begin
     clusterexpect = qpint1d( 'RSIntegrand', object=self, $
                              self.rsresidrange[0], self.rsresidrange[1] )
  endelse
  loglikelihood = clusterlike-clusterexpect

  if keyword_set(verbose) then begin
     print, format='(%"clusterlike:      %11.4f")', clusterlike
     print, format='(%"clusterexpect:    %11.4f")', clusterexpect
     print, format='(%"loglikelihood:  %13.4f")', loglikelihood
  endif

  return, -loglikelihood
end

Pro RSlike::CleanUp
  if ptr_valid(self.clusterdata) then ptr_free, self.clusterdata
  if ptr_valid(self.bkgdata) then ptr_free, self.bkgdata
  if ptr_valid(self.param) then ptr_free, self.param
end

Pro RSlike__Define, struct
  struct = { RSlike, $
             Inherits JEMobject, $
             clusterdata:ptr_new(), $
             clusterarea:0., $
             izerr:0., $
             bkgdata:ptr_new(), $
             bkgarea:0., $
             rsresidrange:fltarr(2), $
             param:ptr_new() $
           }
end
