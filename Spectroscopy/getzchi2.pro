;+
; NAME:
;   getzchi2
;
; PURPOSE:
;   The template fit statistics for object of known redshift.
;
; CALLING SEQUENCE:
;    result = getzchi2( objwave, objflux, objivar, $
;       starflux=starflux, starwave=starwave, [ z=z, npoly=npoly,
;       synflux=synflux ] )
;
; INPUTS:
;   objwave    - Object wavelengths [NOBJPIX].
;   objflux    - Object flux [NOBJPIX]
;   objivar    - Object inverse variance [NOBJPIX]
;
; REQUIRED KEYWORDS:
;   starflux   - Eigenspectra [NPIXSTAR,NSTAR]
;   starwave   - Wavelengths for starflux [NSTARPIX]
;
; OPTIONAL KEYWORDS:
;   z          - Obect redshift. Default z=0.
;   npoly      - Number of polynomial terms to append to eigenspectra;
;                default to none.
;
; OUTPUTS:
;   result     - Structure with fit information.
;
; OPTIONAL OUTPUTS:
;   synflux    - Best fit synthetic flux. [NOBJPIX] 
;   
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
; computechi2()
;
; INTERNAL SUPPORT ROUTINES:
; create_zans()
; cheb_array()
;
; REVISION HISTORY:
;   16-Feb-2007  Written by J. Meyers, LBNL
;------------------------------------------------------------------------------

function create_zans, nstar, nfind
  zans1 = create_struct( $
          name = 'ZANS'+strtrim(string(nstar),1), $
          'z', 0.0, $
          'z_err', 0.0, $
          'chi2',  0.0, $
          'dof',   0L, $
          'theta', fltarr(nstar), $
          'theta_covar', fltarr(nstar,nstar))
  return, replicate(zans1,nfind)
end

function getzchi2, objwave, objflux, objivar, z=z, npoly=npoly, $
                   starflux=starflux, starwave=starwave, synflux=synflux, $
                   acoeff=acoeff, covar=covar
  
  if (n_elements(z) eq 0) then z=0.
  if (n_elements(npoly) eq 0) then npoly=0
  
  if (npoly ne 0) then fixed_template=cheb_array(npoly-1,objwave)

  nobjpix=(size(objflux,/dimens))[0]
  ndim=(size(starflux,/n_dim))
  if ndim eq 1 then nstar=1 else $
     nstar=(size(starflux,/dimens))[1]
  zans1=create_zans(nstar,1)
  
  starwavez=starwave*(1.0+z)
  sqobjivar=sqrt(objivar)
  starfluxz=dblarr(nobjpix,nstar)
  for istar=0,nstar-1 do $
     starfluxz[*,istar]=interpol(starflux[*,istar],starwavez,objwave, $
                                 spline=spline)
  
  idx=where(objwave le max(starwavez) and objwave ge min(starwavez))
  if (idx[0] eq -1) then return,zans1
  if (keyword_set(fixed_template)) then begin
     chi2=computechi2(objflux[idx],sqobjivar[idx], $
                      ([[starfluxz],[fixed_template]])[idx,*], $
                      acoeff=acoeff, covar=covar, dof=dof, yfit=yfit)
     zans1.theta=acoeff[0:nstar-1]
     zans1.theta_covar=covar[0:nstar-1,0:nstar-1]
  endif else begin
     chi2=computechi2(objflux[idx],sqobjivar[idx],starfluxz[idx,*], $
                      acoeff=acoeff, covar=covar,dof=dof, yfit=yfit)
     zans1.theta=acoeff
     zans1.theta_covar=covar
  endelse
  synflux=objwave*0.
  synflux[idx]=yfit
  zans1.z=z
  zans1.z_err=-1.
  zans1.chi2=chi2
  zans1.dof=dof
  return,zans1
end
