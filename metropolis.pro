Function metropolis_logp, param, j=j
  gamma = 0.57721
  logP0 = param[0]
  jstar = param[1]
  alpha = param[2]
  term = (jstar/j)^alpha
  model=logP0 + alog(term/(1+term)) - gamma
  return, model
end


Function metropolis_devlogp, param, j=j, pj=pj
  model = metropolis_logp(param, j=j)
  return, alog(pj)-model
end


Pro metropolis_plotchains, paramchain, valchain, trial=trial
  common metropolisblock, fcn, obj, iter, nparam, covar, $
     scale, silent, nupdate, prefix, parinfo

  nsum=floor(iter/2000L) > 1

  filename=prefix+'.chains.ps'
  if keyword_set(trial) then filename=prefix+'.trialchains.ps'
  thisDevice=!D.Name
  set_plot, 'ps'
  device, xsize=8.5, ysize=11.0, /inches, xoffset=0, yoffset=0
  device, encapsulated=0
  device, /times, /isolatin1
  device, filename=filename
  thisFont = !p.font
  !p.font=0
  !p.multi=[0,2,3]
  for i=0, nparam-1 do begin
     if parinfo[i].fixed then begin
        val=paramchain[i,0]
        if val gt 0 then yrange=[0, val*2] $
        else if val lt 0 then yrange=[-2*val,0] $
                                     else yrange=[-1,1]
        plot, paramchain[i,*], $
              title=parinfo[i].name, /ynozero, nsum=nsum, charsize=1.5, $
              yrange=yrange
     endif else begin
        plot, paramchain[i,*], $
              title=parinfo[i].name, /ynozero, nsum=nsum, charsize=1.5
     endelse
  endfor
  plot, valchain-max(valchain), $
        title='Relative Log Likelihood', yrange=[-15,0], charsize=1.5
  !p.multi=0
  !p.font=thisFont
  device, /close
  set_plot, thisDevice
end


Pro metropolis_plotcorrelations, paramchain, trial=trial
  common metropolisblock

  skip = floor((size(paramchain,/dim))[1]/1000) > 1
  filename=prefix+'.correlations.ps'
  if keyword_set(trial) then filename=prefix+'.trialcorrelations.ps'

  loadct, 0
  color95 = fsc_color('red', !D.Table_Size-2)
  color68 = fsc_color('blue', !D.Table_Size-3)
  thisDevice=!D.Name
  set_plot, 'ps', /copy
  device, xsize=8.5, ysize=11.0, /inches, xoffset=0, yoffset=0
  device, encapsulated=0, color=1
  device, /times, /isolatin1
  device, filename=filename
  thisFont = !p.font
  !p.font = 0
  !p.multi=[0,2,3]
  for i=0, nparam-1 do begin
     if parinfo[i].fixed then continue
     for j=i+1, nparam-1 do begin
        if parinfo[j].fixed then continue
        xrange = percentile(paramchain[i,*], [0.01, 0.99])
        xrange += [-0.1, 0.1]*(xrange[1]-xrange[0])
        yrange = percentile(paramchain[j,*], [0.01, 0.99])
        yrange += [-0.1, 0.1]*(yrange[1]-yrange[0])

        ndim = size(paramchain, /n_dim)
        if ndim eq 1 then nchain = n_elements(paramchain) $
        else nchain = (size(paramchain, /dim))[1]
        value=fltarr(nchain)+1.0
        nx = 5 > floor(sqrt(nchain)) < 30
        ny = 5 > floor(sqrt(nchain)) < 30
        xpos = (paramchain[i,*]-xrange[0])*nx/(xrange[1]-xrange[0])
        ypos = (paramchain[j,*]-yrange[0])*ny/(yrange[1]-yrange[0])
        w=where(xpos gt 0 and xpos le nx and ypos gt 0 and ypos le ny)
        field = cic(value[w], xpos[w], nx, ypos[w], ny, /isolated, /no_message)
        xbinsize=(xrange[1]-xrange[0])/nx
        ybinsize=(yrange[1]-yrange[0])/ny
        xbins=findgen(nx)*xbinsize+xbinsize/2.+xrange[0]
        ybins=findgen(ny)*ybinsize+ybinsize/2.+yrange[0]
        dens=field/total(field)
        denst=dens[sort(dens)]
        cs=total(denst, /cum)
        test=min(abs(cs-(1.0-0.68)),ind68)
        test=min(abs(cs-(1.0-0.95)),ind95)
        levels=[denst[ind95], denst[ind68]]
        delvarx, xpout, ypout
        surf=min_curve_surf(dens, nx=nx, ny=ny, xpout=xpout, ypout=ypout)
        xpout=xbins
        ypout=ybins
        contour, surf, xpout, ypout, levels=levels, /nodata, $
                 xrange=minmax(xpout), xstyle=1, $
                 yrange=minmax(ypout), ystyle=1, $
                 c_colors=[color68, color95], $
                 xtitle=parinfo[i].name, ytitle=parinfo[j].name, $
                 thick=3, xthick=3, ythick=3, $
                 charsize=1.5
        image = bytscl(-surf, top=!D.Table_Size-4)
        tv, image, $
            !X.WINDOW[0], !Y.WINDOW[0], $
            XSIZE = !X.WINDOW[1] - !X.WINDOW[0], $
            YSIZE = !Y.WINDOW[1] - !Y.WINDOW[0], /NORM
        contour, surf, xpout, ypout, levels=levels, /over, $
                 xrange=minmax(xpout), xstyle=1, $
                 yrange=minmax(ypout), ystyle=1, $
                 c_colors=[color68, color95], $
                 xtitle=parinfo[i].name, ytitle=parinfo[j].name, $
                 thick=3, xthick=3, ythick=3, $
                 charsize=1.5
        axis, xaxis=0, xrange=minmax(xpout), xstyle=1, xthick=3, charsize=1.5
        axis, xaxis=1, xrange=minmax(xpout), xstyle=1, xthick=3, xtickname=replicate(' ',30)
        axis, yaxis=0, yrange=minmax(ypout), ystyle=1, ythick=3, charsize=1.5
        axis, yaxis=1, yrange=minmax(ypout), ystyle=1, ythick=3, ytickname=replicate(' ',30)
     endfor
  endfor
  !p.multi=0
  !p.font=thisFont
  device, /close
  set_plot, thisDevice
end

Pro metropolis_plothist, paramchain
  common metropolisblock

  thisDevice=!D.Name
  set_plot, 'ps'
  device, xsize=8.5, ysize=11.0, /inches, xoffset=0, yoffset=0
  device, encapsulated=0
  device, /times, /isolatin1
  device, filename=prefix+'.hist.ps'
  thisFont = !p.font
  !p.font=0
  !p.multi=[0,2,3]
  for i=0, nparam-1 do begin
     if parinfo[i].fixed then continue
     plothist, paramchain[i,*], title=parinfo[i].name, $
               bin=(sqrt(diag_matrix(covar)))[i]/5., $
               thick=3, charsize=1.5
  endfor
  !p.multi=0
  !p.font=thisFont
  device, /close
  set_plot, thisDevice
end


Pro metropolis_plotpower, powerstruct
  common metropolisblock

  thisDevice=!D.Name
  set_plot, 'ps'
  device, xsize=8.5, ysize=11.0, /inches, xoffset=0, yoffset=0
  device, encapsulated=0
  device, /times, /isolatin1
  device, filename=prefix+'.power.ps'
  thisFont = !p.font
  !p.font=0
  !p.multi=[0,2,3]
  for i=0, nparam-1 do begin
     if parinfo[i].fixed then continue
     k =        powerstruct[i].k
     j =        powerstruct[i].j
     chainlen = powerstruct[i].chainlen
     power =    powerstruct[i].power
     parms =    powerstruct[i].parms
     jstar =    powerstruct[i].jstar
     r =        powerstruct[i].r
     plotmax = floor( 50 > jstar*20 < (2000 < chainlen/2) )
     plot, alog10(k[1:plotmax]), power[1:plotmax], /ylog, $
           title=parinfo[i].name, $
           thick=3, xthick=3, ythick=3, charsize=1.5
     oplot, alog10(k[1:*]), exp(metropolis_logp(parms, j=j[1:*]))
     xyouts, !X.Window[0] + 0.1 * (!X.Window[1] - !X.Window[0]), $
             !Y.Window[0] + 0.2 * (!Y.Window[1] - !Y.Window[0]), $
             textoidl('j^* = '+str(jstar)), /normal
     xyouts, !X.Window[0] + 0.1 * (!X.Window[1] - !X.Window[0]), $
             !Y.Window[0] + 0.1 * (!Y.Window[1] - !Y.Window[0]), $
             textoidl('r = '+str(r)), /normal
  endfor
  !p.multi=0
  !p.font=thisFont
  device, /close
  set_plot, thisDevice
end


Function metropolis_power, chains
  common metropolisblock

  chainlen = (size(chains, /dim))[1]
  j = findgen(floor(chainlen/2+1))
  k = 2*!pi/chainlen*j
  if chainlen mod 2 eq 1 then begin ;;ensure even
     chains = chains[*,0:chainlen-2]
     chainlen -= 1
  endif
  powerstruct = { k:k, $
                  j:j, $
                  chainlen:chainlen, $
                  power:fltarr(chainlen)+!values.f_nan, $
                  parms:fltarr(3)+!values.f_nan, $
                  jstar:!values.f_nan, $
                  r:!values.f_nan $
                }
  powerstruct = replicate(powerstruct, nparam)
  for i=0, nparam-1 do begin
     if parinfo[i].fixed then continue
     chain = chains[i,*]
     chain -= mean(chain)
     chain /= stdev(chain)
     power = fft( chain, 1 )
     power /= sqrt(chainlen)
     power = abs(power)^2

     start = [alog(10), 20, 2]
     jmax = floor(1000 < chainlen/2)
     functargs = {j:j[1:jmax], pj:power[1:jmax]}
     parms = mpfit( 'metropolis_devlogp', start, $
                    functargs=functargs, /quiet )
     jstar = parms[1]

     ;;iterate once
     jmax = floor(chainlen/2 < jstar*10. > 50)
     functargs = {j:j[1:jmax], pj:power[1:jmax]}
     parms = mpfit( 'metropolis_devlogp', start, $
                    functargs=functargs, /quiet )
     jstar = parms[1]
     r = exp(parms[0])/chainlen
     powerstruct[i].power=power
     powerstruct[i].parms=parms
     powerstruct[i].jstar=jstar
     powerstruct[i].r=r
  endfor
  return, powerstruct
end


Function metropolis_getval, indepvar
  common metropolisblock

  if ~obj_valid(obj) then $
     return, call_function( fcn, indepvar ) $
  else $
     return, call_method( fcn, obj, indepvar )
end


Pro metropolis_iterate, start, $
                         niter=niter, $
                         paramchain=paramchain, valchain=valchain
  common metropolisblock

  currentparam = start
  currentval = metropolis_getval( currentparam )

  paramchain = fltarr(nparam, niter)
  valchain = fltarr(niter)

  for i=0L, niter-1 do begin
     if ~keyword_set(silent) then $
        counter, i+1+iter, niter+iter, 'Metropolis iterations '
     proposedparam = currentparam + mrandomn( seed, covar*scale )
     w=where(parinfo.fixed)
     if w[0] ne -1 then begin
        proposedparam[w] = start[w]
     endif
     proposedval = metropolis_getval( proposedparam )
     if (proposedval gt currentval) $
        || (exp(proposedval-currentval) gt randomu(seed)) then begin
        currentval = proposedval
        currentparam = proposedparam
        nupdate += 1
     endif
     paramchain[*,i] = currentparam
     valchain[i] = currentval
  endfor
  if ~keyword_set(silent) then print
  iter += niter
end


Function metropolis, fcn1, start, covar1, $
                      object=object, miniter=miniter, $
                      maxiter=maxiter, fval=fval, $
                      paramchain=paramchain, valchain=valchain, $
                      trialparamchain=trialparamchain, $
                      trialvalchain=trialvalchain, $
                      trialniter=trialniter, ntrial=ntrial, $
                      updatefrac=updatefrac, time=time, $
                      silent=silent1, parinfo=parinfo1, $
                      prefix=prefix1, nstop=nstop
  common metropolisblock

  time0 = systime(1)

  if n_elements(object) eq 0 then obj='' else obj=object
  if n_elements(prefix1) ne 0 then prefix=prefix1
  fcn = fcn1
  covar = covar1
  if n_elements(ntrial) eq 0 then ntrial=0
  silent = keyword_set(silent1)
  if n_elements(nstop) eq 0 then nstop=20
  if nstop lt 2 then nstop=2


  iter = 0L
  nupdate=0L
  ;;need currentparam, currentval
  nparam = n_elements(start)

  parinfo = { name:'', $
              jstargoal:20., $
              rgoal:0.01, $
              fixed:0 }
  parinfo = replicate( parinfo, nparam )
  parinfo.name = 'Var'+strtrim(str(indgen(nparam)+1), 2)
  if n_elements(parinfo1) ne 0 then begin
     for i=0, nparam-1 do begin
        tmp = parinfo[i]
        struct_assign, parinfo1[i], tmp, /nozero
        parinfo[i] = tmp
     endfor
  endif

  nfixed = total(parinfo.fixed)
  scale = 2.4^2/(nparam-nfixed) ;; proposal distribution rescaling of the covariance matrix

  ;; optionally run trial chain(s) to determine an efficient covariance matrix
  currentparam = start

  for j=0L, ntrial-1 do begin
     if ~keyword_set(silent) then $
        print, 'Estimating Covariance '+str(j+1)+' of '+str(ntrial)
     metropolis_iterate, currentparam, niter=trialniter, $
                          paramchain=trialparamchain, valchain=trialvalchain
     bestval = max(trialvalchain, wbest)
     wlikely=where(trialvalchain gt bestval-alog(100.))
     if n_elements(wlikely) lt 30 then begin
        s=reverse(sort(trialvalchain))
        wlikely=[s[0:29]]
     endif
     newcovar = covariance(trialparamchain[*,wlikely])
     wfixed = where(parinfo.fixed)
     if wfixed[0] ne -1 then begin
        newcovar[wfixed,*]=0.
        newcovar[*,wfixed]=0.
        newcovar[wfixed, wfixed]=1.
     endif
     if ~keyword_set(silent) then begin
        print
        print, 'Update efficiency: ', nupdate*1.0/trialniter
        print, 'Updated covariance matrix from ', n_elements(wlikely), ' samples.'
        print, '    Variable      Old sigma      New sigma'
        for i=0, nparam-1 do begin
           print, parinfo[i].name, $
                  (sqrt(diag_matrix(covar)))[i], $
                  (sqrt(diag_matrix(newcovar)))[i], $
                  f='(A12,2F15.4)'
        endfor
     endif
     covar = newcovar
     ;; set currentparam to best found param in trial
     currentparam = trialparamchain[*,wbest]
     if n_elements(prefix) ne 0 then begin
        metropolis_plotchains, trialparamchain, $
                                trialvalchain, /trial
        metropolis_plotcorrelations, trialparamchain[*,wlikely], $
                                      /trial
     endif
  endfor

  ;; breakup into 20 possible stopping points
  stoppoint = floor([0,exp(findgen(nstop)/(nstop-1) $
                           *(alog(1.0*maxiter/miniter))+alog(miniter))])
  ;; find stoppoints whose prime factors are 2,3,5
  candstop = powers([2,3,5], maxiter)
  for j=1, nstop-1 do begin
     stoppoint[j] = min(candstop[where(candstop gt stoppoint[j])])
  endfor
  ;; eliminate duplications
  stoppoint = stoppoint[uniq(stoppoint)]
  nstop = n_elements(stoppoint)-1

  ;; Start actual iterations
  iter = 0L
  nupdate = 0L
  paramchain = fltarr(nparam, maxiter)
  valchain = fltarr(maxiter)-(!values.f_nan)
  jstar = fltarr(nparam)
  r = fltarr(nparam)

  for j=0,nstop-1 do begin
     niter = stoppoint[j+1]-stoppoint[j]
     metropolis_iterate, currentparam, niter=niter, $
                          paramchain=paramchain1, valchain=valchain1
     paramchain[*,stoppoint[j]:stoppoint[j+1]-1]=paramchain1
     valchain[stoppoint[j]:stoppoint[j+1]-1]=valchain1
     wlikely = where(valchain[0:stoppoint[j+1]-1] gt max(valchain)-alog(100.))
     burn = min(wlikely)
     bestval = max(valchain,m)
     bestparam = transpose(paramchain[*,m])

     powerstruct = metropolis_power(paramchain[*,burn:stoppoint[j+1]-1])
     if n_elements(prefix) ne 0 then begin
        metropolis_plotchains, paramchain[*,0:stoppoint[j+1]-1], $
                                valchain[0:stoppoint[j+1]-1]
        metropolis_plotcorrelations, paramchain[*,wlikely]
        metropolis_plotpower, powerstruct
        metropolis_plothist, paramchain[*,burn:stoppoint[j+1]-1]
     endif

     if ~keyword_set(silent) then begin
        print
        print
        print, '    variable          value          error          jstar              r'
        for i=0, nparam-1 do begin
           print, parinfo[i].name, $
                  bestparam[i], $
                  biweight_scale(paramchain[i,burn:stoppoint[j+1]-1]-bestparam[i], /zero), $
                  powerstruct[i].jstar, $
                  powerstruct[i].r, $
                  f='(A12,2F15.5,F15.2,F15.5)'
        endfor
     endif

     if (total(powerstruct.jstar lt parinfo.jstargoal) eq 0) $
        and (total(powerstruct.r gt parinfo.rgoal) eq 0) then break
  endfor
  paramchain = paramchain[*,0:iter-1]
  valchain = valchain[0:iter-1,*]

  updatefrac = nupdate*1.0/iter
  fval = max(valchain, m)
  bestparam = transpose(paramchain[*,m])
  results = replicate({parammean:0., parambest:0., paramerr:0.}, nparam)
  for i=0, nparam-1 do begin
     results[i].parammean = biweight_location(paramchain[i,wlikely])
     results[i].paramerr = biweight_scale(paramchain[i,wlikely])
     results[i].parambest = bestparam[i]
  endfor
  time = systime(1)-time0

  parinfo1 = struct_addtags(parinfo, struct_trimtags(powerstruct, select=['R','JSTAR']))
  parinfo1 = struct_addtags(parinfo1, results)

  return, bestparam
end
