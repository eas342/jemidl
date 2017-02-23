Function stex_icolor, color

   if (n_elements(color) EQ 0) then color='white'

   ncolor = N_elements(color)

   ; If COLOR is a string or array of strings, then convert color names
   ; to integer values
   if (size(color,/tname) EQ 'STRING') then begin ; Test if COLOR is a string

      ; Detemine the default color for the current device
      if (!d.name EQ 'X') then defcolor = 7 $ ; white for X-windows
       else defcolor = 0 ; black otherwise

      icolor = 0 * (color EQ 'black') $
             + 1 * (color EQ 'red') $
             + 2 * (color EQ 'green') $
             + 3 * (color EQ 'blue') $
             + 4 * (color EQ 'cyan') $
             + 5 * (color EQ 'magenta') $
             + 6 * (color EQ 'yellow') $
             + 7 * (color EQ 'white') $
             + defcolor * (color EQ 'white')

;     if (!d.N_colors EQ 16777216) then begin
;        red   = [0, 1, 0, 0, 0, 1, 1, 1]
;        green = [0, 0, 1, 0, 1, 0, 1, 1]
;        blue  = [0, 0, 0, 1, 1, 1, 0, 1]
;        colors = 255L*red + ishft(255L*green,8) + ishft(255L*blue,16)
;        icolor = colors[icolor]
;     endif

   endif else begin
      icolor = long(color)
   endelse

   return, icolor
end

;------------------------------------------------------------------------------

Pro inserttreenode, base, levels
  widget_control, base, get_uvalue=base_uvalue
  if n_elements(levels) gt 1 then begin
     if ptr_valid(base_uvalue.children) then begin ;; check children
        children = *(base_uvalue.children)
        for i=0, n_elements(children)-1 do begin
           child_ID = children[i]
           widget_control, child_ID, get_uvalue=child_uvalue
           if child_uvalue.subpath eq levels[0] then begin ;; child matches
              break
           endif
        endfor
        if i eq n_elements(children) then begin ;; no matching children
           wtNode = widget_tree( base, value=levels[0], $
                                 uvalue={ filterbranch, $
                                          subpath:levels[0], $
                                          parent:base, $
                                          children:ptr_new() }, $
                                 /folder )
           children = [children, wtNode]
           child_ID = wtNode
           ptr_free, base_uvalue.children
           base_uvalue.children = ptr_new(children)
           widget_control, base, set_uvalue=base_uvalue
        endif
     endif else begin ;; add first child
        wtNode = widget_tree( base, value=levels[0], $
                              uvalue={ filterbranch, $
                                       subpath:levels[0], $
                                       parent:base, $
                                       children:ptr_new() }, $
                              /folder )
        child_ID = wtNode
        base_uvalue.children=ptr_new([wtNode])
        widget_control, base, set_uvalue=base_uvalue
     endelse
     inserttreenode, child_ID, levels[1:*]  ;; go recursive
  endif else begin ;; add leaf, no recursion
     if ptr_valid(base_uvalue.children) then begin
        children = *(base_uvalue.children)
        for i=0, n_elements(children)-1 do begin
           child_ID = children[i]
           widget_control, child_ID, get_uvalue=child_uvalue
           if child_uvalue.subpath eq levels[0] then begin ;; child matches
              break
           endif
        endfor
        if i eq n_elements(children) then begin ;; no matching children
           wtNode = widget_tree( base, value=levels[0], $
                                 uvalue={ filterleaf, $
                                          subpath:levels[0], $
                                          parent:base, $
                                          children:ptr_new() } )
           children = [children, wtNode]
           ptr_free, base_uvalue.children
           base_uvalue.children = ptr_new(children)
           widget_control, base, set_uvalue=base_uvalue
        endif else begin
           message, 'duplicate leaf'
        endelse
     endif else begin
        wtNode = widget_tree( base, value=levels[0], $
                              uvalue={ filterleaf, $
                                       subpath:levels[0], $
                                       parent:base, $
                                       children:ptr_new() } )
        children = [wtNode]
        base_uvalue.children = ptr_new(children)
        widget_control, base, set_uvalue=base_uvalue
     endelse
  endelse
end

;------------------------------------------------------------------------------

Function filterchoice, group_leader=groupleader, cancel=cancel
  filterdir = getenv('JEM_REF')+'FILTERS/'
  files=file_search(filterdir, '*.dat', count=nfiles)

  catch, theerror
  if theerror ne 0 then begin
     catch, /cancel
     ok = dialog_message(!error_state.msg)
     if destroy_groupleader then widget_control, groupleader, /destroy
     cancel = 1
     return, ""
  endif

  if n_elements(groupleader) eq 0 then begin
     groupleader = widget_base(map=0)
     widget_control, groupleader, /realize
     destroy_groupleader = 1
  endif else destroy_groupleader = 0

  wTLB = widget_base(/column, title='Filters', uvalue='', group_leader=groupleader)
  wTree = widget_tree(wTLB)
  wtRoot = widget_tree( wTree, $
                        value='Root', $
                        uvalue={ filterbranch, $
                                 subpath:'', $
                                 parent:0, $
                                 children:ptr_new() }, $
                        /folder, $
                        /expanded )
  for i=0, nfiles-1 do begin
     path=strmid(files[i],strlen(filterdir))
     inserttreenode, wtRoot, strsplit(path, '/', /extract)
  endfor
  wDone = widget_button(wTLB, value='Done', uvalue='DONE')
  wCancel = widget_button(wTLB, value='Cancel', uvalue='CANCEL')

  widget_control, wTLB, /realize
  ptr = ptr_new({filter:'', cancel:1})
  info = {filter:'', cancel:0, wTree:wTree, wtRoot:wtRoot, wDone:wDone, wCancel:wCancel}
  infoptr = ptr_new(info, /no_copy)
  widget_control, wTLB, set_uvalue=infoptr
  xmanager, 'filterchoice', wTLB
  filter = (*infoptr).filter
  cancel = (*infoptr).cancel
  ptr_free, infoptr
  if destroy_groupleader then widget_control, groupleader, /destroy
  return, filter
end

;------------------------------------------------------------------------------

Pro filterchoice_event, ev
  compile_opt hidden
  widget_control, ev.ID, get_uvalue=uvalue
  widget_control, ev.top, get_uvalue=infoptr
  t = size(uvalue, /tname)
  case t of
     'STRING': begin
        if uvalue eq 'CANCEL' then (*infoptr).cancel=1
        widget_control, ev.top, /destroy
     end
     'STRUCT': begin
        if tag_names(uvalue, /structure_name) eq 'FILTERLEAF' then begin
           uvalue1=uvalue
           path=''
           repeat begin
              path = [path, uvalue1.subpath]
              parent=uvalue1.parent
              widget_control, parent, get_uvalue=uvalue1
           endrep until parent eq (*infoptr).wtRoot
           (*infoptr).filter=strjoin(reverse(path[1:*]),'/')
        endif
     end
  endcase
end

;------------------------------------------------------------------------------

Function colorbox, group_leader=groupleader, cancel=cancel

  catch, theerror
  if theerror NE 0 then begin
     catch, /cancel
     ok = dialog_message(!error_state.msg)
     if destroy_groupleader then widget_control, groupleader, /destroy
     cancel = 1
     return, ""
  endif

  if n_elements(groupleader) eq 0 then begin
     groupleader = widget_base(map=0)
     widget_control, groupleader, /realize
     destroy_groupleader = 1
  endif else destroy_groupleader = 0

  filterdir = getenv('JEM_REF')+'FILTERS/'
  filters=file_search(filterdir, '*.dat', count=nfiles)
  filters=strmid(filters, strlen(filterdir))
  specdir = getenv('JEM_REF')+'SPECTRA/'
  spectra=file_search(specdir, '*.dat', count=nfiles)
  spectra=strmid(spectra, strlen(specdir))

  tlb = widget_base(title='Color filters', column=1, /modal, $
                    /base_align_center, group_leader=groupleader)

  label1_base = widget_base(tlb, row=1)
  label1ID = widget_label(label1_base, value='filter 1')
  combobox1_base = widget_base(tlb, row=1)
  combobox1ID = widget_combobox(combobox1_base, value=filters)

  label2_base = widget_base(tlb, row=1)
  label2ID = widget_label(label2_base, value='zeropoint 1')
  combobox2_base = widget_base(tlb, row=1)
  combobox2ID = widget_combobox(combobox2_base, value=spectra)

  label3_base = widget_base(tlb, row=1)
  label3ID = widget_label(label3_base, value='filter 2')
  combobox3_base = widget_base(tlb, row=1)
  combobox3ID = widget_combobox(combobox3_base, value=filters)

  label4_base = widget_base(tlb, row=1)
  label4ID = widget_label(label4_base, value='zeropoint 1')
  combobox4_base = widget_base(tlb, row=1)
  combobox4ID = widget_combobox(combobox4_base, value=spectra)

  redshift_base = widget_base(tlb, row=1)
  redshiftID = cw_field(redshift_base, $
                        /float, $
                        value = 0.0, $
                        title = 'Redshift:', $
                        xsize = 12)

  buttonbase = widget_base(tlb, row=1)
  cancelID = widget_button(buttonbase, value='Cancel')
  acceptID = widget_button(buttonbase, value='Accept')

  widget_control, tlb, /realize
  ptr = ptr_new({filter1:"", filter2:"", zp1:"", zp2:"", redshift:0.d, cancel:1})

  info = {ptr:ptr, combobox1ID:combobox1ID, combobox2ID:combobox2ID, $
          combobox3ID:combobox3ID, combobox4ID:combobox4ID, $
          redshiftID:redshiftID, cancelID:cancelID, acceptID:acceptID }
  widget_control, tlb, set_uvalue = info, /no_copy

  xmanager, 'colorbox', tlb

  color = {filter1:(*ptr).filter1, filter2:(*ptr).filter2, $
           zp1:(*ptr).zp1, zp2:(*ptr).zp2, $
           redshift:(*ptr).redshift}
  cancel = (*ptr).cancel
  ptr_free, ptr
  if destroy_groupleader then widget_control, groupleader, /destroy

  return, color
end

;------------------------------------------------------------------------------

Pro colorbox_event, event
  widget_control, event.top, get_uvalue=info
  case event.id of
     info.cancelID: widget_control, event.top, /destroy
     info.acceptID: begin
        filter1 = widget_info(info.combobox1ID, /combobox_gettext)
        zp1 = widget_info(info.combobox2ID, /combobox_gettext)
        filter2 = widget_info(info.combobox3ID, /combobox_gettext)
        zp2 = widget_info(info.combobox4ID, /combobox_gettext)
        widget_control, info.redshiftID, get_value=redshift
        (*info.ptr).filter1=filter1[0]
        (*info.ptr).zp1=zp1[0]
        (*info.ptr).filter2=filter2[0]
        (*info.ptr).zp2=zp2[0]
        (*info.ptr).redshift=redshift
        (*info.ptr).cancel = 0
        widget_control, event.top, /destroy
     end
     else:
  endcase
  return
end

;------------------------------------------------------------------------------

Pro stex_makecolor

  common stex_common, nstex, istex, stexbases
  widget_control, /hourglass
  widget_control, stexbases[istex], get_uvalue=state

  color = colorbox( group_leader=stexbases[istex], cancel=cancelled )
  if cancelled then return

  curves = ptrarr(state.ndata)

  filterdir=getenv('JEM_REF')+'FILTERS/'
  readcol, filterdir+color.filter1, f1wave, f1tput, /silent
  readcol, filterdir+color.filter2, f2wave, f2tput, /silent
  specdir=getenv('JEM_REF')+'SPECTRA/'
  readcol, specdir+color.zp1, zp1wave, zp1flux, /silent
  readcol, specdir+color.zp2, zp2wave, zp2flux, /silent

  ;; for idata=0, state.ndata-1 do begin
  ;;    w=(*state.data[idata]).wave*(1.+color.redshift)
  ;;    day=(*state.data[idata]).day
  ;;    flux=(*state.data[idata]).flux
  ;;    ww=rebin(w, n_elements(w), n_elements(day))

  ;;    f1tputinterp = interpol(f1tput, ftwave, w)
  ;;    filter1 = rebin(filter1, n_elements(w), n_elements(day))
  ;;    filter2 = make_filter(w, color.filter2)
  ;;    filter2 = rebin(filter2, n_elements(w), n_elements(day))

  ;;    dlambda = w-shift(w,1)
  ;;    dlambda[0] = w[1]-w[0]
  ;;    dlambda = rebin(dlambda, n_elements(w), n_elements(day))

  ;;    product1 = flux * ww * dlambda
  ;;    product2 = product1

  ;;    product1 *= filter1
  ;;    product2 *= filter2

  ;;    LCflux1 = total(product1,1)
  ;;    LCflux2 = total(product2,1)

  ;;    c = -2.5*alog10(LCflux1) - vega_zeropoint(color.filter1) + $
  ;;        2.5*alog10(LCflux2) + vega_zeropoint(color.filter2)
  ;;    curve = {day:day, curve:c}
  ;;    curves[idata]=ptr_new(curve)
  ;; endfor

  for idata=0, state.ndata-1 do begin
     w=(*state.data[idata]).wave
     day=(*state.data[idata]).day
     flux=(*state.data[idata]).flux
     c= synphot(w,flux,color.redshift,f_wave=f1wave,f_throughput=f1tput,zp_wave=zp1wave,zp_flux=zp1flux) $
        - synphot(w,flux,color.redshift,f_wave=f2wave,f_throughput=f2tput,zp_wave=zp2wave,zp_flux=zp2flux)
     curve = {day:day, curve:c}
     curves[idata] = ptr_new(curve)
  endfor

  stex_make_curve, color.filter1+' - '+color.filter2+' at z='+string(color.redshift), curves
  stex_plot_curves

  return
end

;------------------------------------------------------------------------------

Pro stex_makeKcorrection
  common stex_common
  widget_control, /hourglass
  widget_control, stexbases[istex], get_uvalue=state

  cancelled=0
  color = colorbox( group_leader=stexbases[istex], cancel=cancelled)
  if cancelled then return

  curves = ptrarr(state.ndata)

  for idata=0, state.ndata-1 do begin
     wave=(*state.data[idata]).wave
     w=wave*(1.+color.redshift)
     day=(*state.data[idata]).day
     flux=(*state.data[idata]).flux

     filter1 = make_filter(w, color.filter1)
     filter1 = rebin(filter1, n_elements(w), n_elements(day))
     filter2 = make_filter(wave, color.filter2)
     filter2 = rebin(filter2, n_elements(wave), n_elements(day))

     dlambda = wave-shift(wave,1)
     dlambda[0] = wave[1]-wave[0]
     dlambda = rebin(dlambda, n_elements(wave), n_elements(day))

     product1 = filter1*flux*dlambda*(1.+color.redshift)*rebin(w, n_elements(w), n_elements(day))
     product2 = filter2*flux*dlambda*rebin(wave, n_elements(wave), n_elements(day))

     LCflux1 = total(product1, 1)
     LCflux2 = total(product2, 1)

     c = -2.5*alog10(LCflux1) - vega_zeropoint(color.filter1) $
         + 2.5*alog10(LCflux2) + vega_zeropoint(color.filter2) $
         + 2.5*alog10(1.+color.redshift)
     curves[idata] = ptr_new({day:day, curve:c})
  endfor

  stex_make_curve, 'Kcorrection : '+color.filter1+' to '+color.filter2+' at z = ' $
                   + string(color.redshift, format='(f9.5)'), curves
  stex_plot_curves
  return
end

;------------------------------------------------------------------------------

Function magbox, group_leader=groupleader, cancel=cancel

  catch, theerror
  if theerror NE 0 then begin
     catch, /cancel
     ok = dialog_message(!error_state.msg)
     if destroy_groupleader then widget_control, groupleader, /destroy
     cancel = 1
     return, ""
  endif

  if n_elements(groupleader) eq 0 then begin
     groupleader = widget_base(map=0)
     widget_control, groupleader, /realize
     destroy_groupleader = 1
  endif else destroy_groupleader = 0

  filterdir = getenv('JEM_REF')+'FILTERS/'
  filters=file_search(filterdir, '*.dat', count=nfiles)
  filters=strmid(filters, strlen(filterdir))
  specdir = getenv('JEM_REF')+'SPECTRA/'
  spectra=file_search(specdir, '*.dat', count=nfiles)
  spectra=strmid(spectra, strlen(specdir))

  tlb = widget_base(title='Magnitude filters', column=1, /modal, $
                    /base_align_center, group_leader=groupleader)

  label1_base = widget_base(tlb, row=1)
  label1ID = widget_label(label1_base, value='filter 1')
  combobox1_base = widget_base(tlb, row=1)
  combobox1ID = widget_combobox(combobox1_base, value=filters)

  label2_base = widget_base(tlb, row=1)
  label2ID = widget_label(label2_base, value='zeropoint 1')
  combobox2_base = widget_base(tlb, row=1)
  combobox2ID = widget_combobox(combobox2_base, value=spectra)

  redshift_base = widget_base(tlb, row=1)
  redshiftID = cw_field(redshift_base, $
                        /float, $
                        value = 0.0, $
                        title = 'Redshift:', $
                        xsize = 12)

  buttonbase = widget_base(tlb, row=1)
  cancelID = widget_button(buttonbase, value='Cancel')
  acceptID = widget_button(buttonbase, value='Accept')

  widget_control, tlb, /realize
  ptr = ptr_new({filter1:"", zp1:"", redshift:0.d, cancel:1})

  info = {ptr:ptr, combobox1ID:combobox1ID, combobox2ID:combobox2ID, $
          redshiftID:redshiftID, cancelID:cancelID, acceptID:acceptID }
  widget_control, tlb, set_uvalue = info, /no_copy

  xmanager, 'magbox', tlb

  mag = {filter1:(*ptr).filter1, $
           zp1:(*ptr).zp1, $
           redshift:(*ptr).redshift}
  cancel = (*ptr).cancel
  ptr_free, ptr
  if destroy_groupleader then widget_control, groupleader, /destroy

  return, mag
end

;------------------------------------------------------------------------------

Pro magbox_event, event
  widget_control, event.top, get_uvalue=info
  case event.id of
     info.cancelID: widget_control, event.top, /destroy
     info.acceptID: begin
        filter1 = widget_info(info.combobox1ID, /combobox_gettext)
        zp1 = widget_info(info.combobox2ID, /combobox_gettext)
        widget_control, info.redshiftID, get_value=redshift
        (*info.ptr).filter1=filter1[0]
        (*info.ptr).zp1=zp1[0]
        (*info.ptr).redshift=redshift
        (*info.ptr).cancel = 0
        widget_control, event.top, /destroy
     end
     else:
  endcase
  return
end

;------------------------------------------------------------------------------

Pro stex_makefluxLC, filterfile
  common stex_common
  widget_control, /hourglass
  widget_control, stexbases[istex], get_uvalue=state

  cancelled = 0
  if n_elements(filterfile) le 0 then $
     filterfile = filterchoice( group_leader=stexbases[istex], cancel=cancelled )
  if cancelled then return

  curves = ptrarr(state.ndata)

  filterdir=getenv('JEM_REF')+'FILTERS/'
  readcol, filterdir+filterfile, fwave, ftput0, /silent

  for idata=0, state.ndata-1 do begin
     wave=(*state.data[idata]).wave
     day=(*state.data[idata]).day
     flux=(*state.data[idata]).flux

     ftput = interpol(ftput0, fwave, wave)
     w=where(wave gt max(fwave) or wave lt min(fwave))
     ftput[w]=0
     ftput = rebin(ftput, n_elements(wave), n_elements(day))

     dlambda = wave-shift(wave,1)
     dlambda[0] = wave[1]-wave[0]
     dlambda = rebin(dlambda, n_elements(wave), n_elements(day))

     product = ftput*flux*dlambda*rebin(wave, n_elements(wave), n_elements(day))
     LCflux = total(product,1)

     curves[idata] = ptr_new({day:day, curve:LCflux})
  endfor

  stex_make_curve, filterfile+' flux lightcurve', curves
  stex_plot_curves
  return
end

;------------------------------------------------------------------------------

Pro stex_makemagLC, filterfile
  common stex_common
  widget_control, /hourglass
  widget_control, stexbases[istex], get_uvalue=state

  cancelled = 0
  if n_elements(filterfile) le 0 then $
     mag = magbox( group_leader=stexbases[istex], cancel=cancelled )
  if cancelled then return

  curves = ptrarr(state.ndata)

  filterdir=getenv('JEM_REF')+'FILTERS/'
  readcol, filterdir+mag.filter1, f1wave, f1tput, /silent
  specdir=getenv('JEM_REF')+'SPECTRA/'
  readcol, specdir+mag.zp1, zp1wave, zp1flux, /silent

  ;; for idata=0, state.ndata-1 do begin
  ;;    wave=(*state.data[idata]).wave
  ;;    day=(*state.data[idata]).day
  ;;    flux=(*state.data[idata]).flux

  ;;    ftput = interpol(ftput0, fwave, wave)
  ;;    w=where(wave gt max(fwave) or wave lt min(fwave))
  ;;    ftput[w]=0
  ;;    ftput = rebin(ftput, n_elements(wave), n_elements(day))

  ;;    dlambda = wave-shift(wave,1)
  ;;    dlambda[0] = wave[1]-wave[0]
  ;;    dlambda = rebin(dlambda, n_elements(wave), n_elements(day))

  ;;    product = ftput*flux*dlambda*rebin(wave, n_elements(wave), n_elements(day))
  ;;    LCflux = total(product,1)
  ;;    LCmag = -2.5*alog10(LCflux); - vega_zeropoint(filterfile)

  ;;    curves[idata] = ptr_new({day:day, curve:LCmag})
  ;; endfor

  for idata=0, state.ndata-1 do begin
     w=(*state.data[idata]).wave
     day=(*state.data[idata]).day
     flux=(*state.data[idata]).flux
     m = synphot(w,flux,mag.redshift,f_wave=f1wave,f_throughput=f1tput,zp_wave=zp1wave,zp_flux=zp1flux)
     curve = {day:day, curve:m}
     curves[idata] = ptr_new(curve)
  endfor

  stex_make_curve, mag.filter1+' magnitude lightcurve', curves, /mag
  stex_plot_curves
  return
end

;------------------------------------------------------------------------------

Pro stex_make_curve, name, curves, mag=mag
  common stex_common
  widget_control, stexbases[istex], get_uvalue=state

  if n_elements(mag) eq 0 then mag=0

  draw_window_size = [600L, 200L]
  base = widget_base(title = 'STEX '+name, $
                     /column, /base_align_right, /tlb_size_events, $
                     /tlb_kill_request_events, group_leader = stexbases[istex])

  base1 = widget_base(base, /row, /base_align_right)
  base2 = widget_base(base, /row, /base_align_right)

  zoomone_button = widget_button(base1, value='Zoom1')

  tmp_string = string('',format='(a30)')
  location_bar_id = widget_label(base1, value=tmp_string, $
                                 uvalue = 'location', frame=1)

  keyboard_text_id = widget_text(base2, $
                                 /all_events, $
                                 scr_xsize = 1, $
                                 scr_ysize = 1, $
                                 units = 0, $
                                 uvalue = 'keyboard_text', $
                                 value='')
  draw_base_id = widget_base(base2, $
                             /column, /base_align_left, $
                             /tracking_events, $
                             uvalue='draw_base', $
                             frame=2)
  draw_widget_id = widget_draw(draw_base_id, $
                               uvalue='draw_window', $
                               /motion_events, $
                               /button_events, $
                               scr_xsize = draw_window_size[0], $
                               scr_ysize = draw_window_size[1])
  widget_control, base, /realize
  widget_control, draw_widget_id, get_value = tmp_value
  draw_window_id = tmp_value

  str = { window_id: draw_window_id, $
          curves: curves, $
          mag: mag, $
          location_bar_id: location_bar_id, $
          draw_widget_id: draw_widget_id, $
          draw_base_id: draw_base_id, $
          draw_window_size: draw_window_size, $
          keyboard_text_id: keyboard_text_id, $
          zoomone_button: zoomone_button, $
          mouse:[0L, 0L], $
          xrange:[0.d0, 1.d0], $
          yrange:[0.d0, 1.d0], $
          position:[0.15d0, 0.10d0, 0.95d0, 0.95d0], $
          mphys: [0.0, 0.0], $
          base_pad: [0L, 0L], $
          pad: [0L, 0L], $
          base_min_size: [600L, 225L] }

  xrange = minmax(*state.dayptr,/nan)
  extra = (max(xrange)-min(xrange))*0.1
  xrange += [-extra, extra]

  if str.mag then begin
     yrange = [-!values.f_infinity, !values.f_infinity]
     for icurve=0,n_elements(str.curves)-1 do begin
        yrange[0] = yrange[0] > max((*str.curves[icurve]).curve,/nan)
        yrange[1] = yrange[1] < min((*str.curves[icurve]).curve,/nan)
     endfor
     extra = (max(yrange)-min(yrange))*0.1
     yrange += [extra, -extra]
  endif else begin
     yrange = [!values.f_infinity, -!values.f_infinity]
     for icurve=0,n_elements(str.curves)-1 do begin
        yrange[0] = yrange[0] < min((*str.curves[icurve]).curve,/nan)
        yrange[1] = yrange[1] > max((*str.curves[icurve]).curve,/nan)
     endfor
     extra = (max(yrange)-min(yrange))*0.1
     yrange += [-extra, extra]
  endelse
  str.xrange=xrange
  str.yrange=yrange

  basegeom = widget_info(base, /geometry)
  drawbasegeom = widget_info(draw_base_id, /geometry)
  str.pad[0] = basegeom.xsize - str.draw_window_size[0]
  str.pad[1] = basegeom.ysize - str.draw_window_size[1]
  str.base_pad[0] = basegeom.xsize - drawbasegeom.xsize $
                    + (2 * basegeom.margin)
  str.base_pad[1] = basegeom.ysize - drawbasegeom.ysize $
                    + (2 * basegeom.margin)

  widget_control, base, set_uvalue = str

  xmanager, 'stex_curve', base, /no_block

  state.curveids[state.ncurves] = base
  state.ncurves += 1

  widget_control, stexbases[istex], set_uvalue=state

  stex_plot_curves

  return
end

;------------------------------------------------------------------------------

Pro stex_curve_event, event
  common stex_common
  widget_control, stexbases[istex], get_uvalue=state
  widget_control, event.top, get_uvalue=str

  if tag_names(event, /structure_name) eq 'WIDGET_KILL_REQUEST' then begin
     i = where(state.curveids eq event.id)
     if i eq state.ncurves-1 then begin
        state.curveids[i]=0
        widget_control, event.id, /destroy
        state.ncurves -= 1
     endif else begin
        state.curveids[i:state.ncurves-2]=state.curveids[i+1:state.ncurves-1]
        state.curveids[state.ncurves-1] = 0
        state.ncurves -= 1
        widget_control, event.id, /destroy
     endelse
     widget_control, stexbases[istex], set_uvalue=state
     return
  endif

  state.icurve=(where(event.top eq state.curveids))[0]
  widget_control, stexbases[istex], set_uvalue=state

  case event.id of
     str.draw_widget_id: begin
        if (event.type eq 2) then begin ;motion event
           temp_event = [event.x, event.y]
           str.mouse = temp_event
           widget_control, event.top, set_uvalue=str
           stex_curve_gettrack
        endif
        if (event.type eq 0) then begin
           case event.modifiers of
              1: xycase = 'yonly'
              4: xycase = 'xandy'
              else: xycase = 'xonly'
           endcase
           xycase = keyword_set(event.modifiers) ? 'xandy' : 'xonly'
           case event.press of
              1: stex_curve_zoom, xycase+'-in', /recenter
              2: stex_curve_zoom, xycase+'-recen', /recenter
              4: stex_curve_zoom, xycase+'-out', /recenter
              else: print, 'trouble in stex_curve_event, mouse zoom'
           endcase
           stex_curve_gettrack
        endif
        widget_control, str.keyboard_text_id, /input_focus
     end
     str.zoomone_button: begin
        xrange = minmax(*state.dayptr,/nan)
        extra = (max(xrange)-min(xrange))*0.1
        xrange += [-extra, extra]
        if str.mag then begin
           yrange = [-!values.f_infinity, !values.f_infinity]
           for icurve=0,n_elements(str.curves)-1 do begin
              yrange[0] = yrange[0] > max((*str.curves[icurve]).curve,/nan)
              yrange[1] = yrange[1] < min((*str.curves[icurve]).curve,/nan)
           endfor
           extra = (max(yrange)-min(yrange))*0.1
           yrange += [extra, -extra]
        endif else begin
           yrange = [!values.f_infinity, -!values.f_infinity]
           for icurve=0,n_elements(str.curves)-1 do begin
              yrange[0] = yrange[0] < min((*str.curves[icurve]).curve,/nan)
              yrange[1] = yrange[1] > max((*str.curves[icurve]).curve,/nan)
           endfor
           extra = (max(yrange)-min(yrange))*0.1
           yrange += [-extra, extra]
        endelse
        str.xrange=xrange
        str.yrange=yrange
        widget_control, event.top, set_uvalue=str
        stex_plot_curves
     end
     event.top: begin
        stex_curve_resize, event
        stex_plot_curves
     end
     else:
  endcase
  return
end

;------------------------------------------------------------------------------

Pro stex_curve_zoom, zchange, recenter=recenter
  common stex_common
  widget_control, stexbases[istex], get_uvalue=state
  widget_control, state.curveids[state.icurve], get_uvalue=str

   case zchange of
      'xonly-in': begin
         str.xrange = str.mphys[0] $
          + [-0.25, 0.25] * (str.xrange[1] - str.xrange[0])
         widget_control, state.curveids[state.icurve], set_uvalue=str
      end
      'yonly-in': begin
         str.yrange = str.mphys[1] $
          + [-0.25, 0.25] * (str.yrange[1] - str.yrange[0])
         widget_control, state.curveids[state.icurve], set_uvalue=str
      end
      'xandy-in': begin
         str.xrange = str.mphys[0] $
          + [-0.25, 0.25] * (str.xrange[1] - str.xrange[0])
         str.yrange = str.mphys[1] $
          + [-0.25, 0.25] * (str.yrange[1] - str.yrange[0])
         widget_control, state.curveids[state.icurve], set_uvalue=str
      end
      'xonly-out': begin
         str.xrange = str.mphys[0] $
          + [-1.0, 1.0] * (str.xrange[1] - str.xrange[0])
         widget_control, state.curveids[state.icurve], set_uvalue=str
      end
      'yonly-out': begin
         str.yrange = str.mphys[1] $
          + [-1.0, 1.0] * (str.yrange[1] - str.yrange[0])
         widget_control, state.curveids[state.icurve], set_uvalue=str
      end
      'xandy-out': begin
         str.xrange = str.mphys[0] $
          + [-1.0, 1.0] * (str.xrange[1] - str.xrange[0])
         str.yrange = str.mphys[1] $
          + [-1.0, 1.0] * (str.yrange[1] - str.yrange[0])
         widget_control, state.curveids[state.icurve], set_uvalue=str
      end
      'one': begin
         stex_auto_zoom
      end
      'xonly-recen': begin ; no change to zoom level: X recenter on mouse pos'n
         str.xrange = str.mphys[0] $
          + [-0.5, 0.5] * (str.xrange[1] - str.xrange[0])
         widget_control, state.curveids[state.icurve], set_uvalue=str
      end
      'yonly-recen': begin ; no change to zoom level: X recenter on mouse pos'n
         str.yrange = str.mphys[1] $
          + [-0.5, 0.5] * (str.yrange[1] - str.yrange[0])
         widget_control, state.curveids[state.icurve], set_uvalue=str
      end
      'xandy-recen': begin ; no change to zoom level: Y recenter on mouse pos'n
         str.xrange = str.mphys[0] $
          + [-0.5, 0.5] * (str.xrange[1] - str.xrange[0])
         str.yrange = str.mphys[1] $
          + [-0.5, 0.5] * (str.yrange[1] - str.yrange[0])
         widget_control, state.curveids[state.icurve], set_uvalue=str
      end
      else: print, 'Problem in splot_zoom!'
   endcase

   stex_plot_curves
   return
end

;------------------------------------------------------------------------------

Pro stex_curve_resize, event
  widget_control, event.top, get_uvalue=str

  tmp_event = [event.x, event.y]
  window = (str.base_min_size > tmp_event)

  newbase = window - str.base_pad
  newsize = window - str.pad

  widget_control, str.draw_base_id, $
                  xsize = newbase[0], ysize = newbase[1]
  widget_control, str.draw_widget_id, $
                  xsize = newsize[0], ysize = newsize[1]

  str.draw_window_size = newsize

  widget_control, event.top, set_uvalue=str
  return
end
;------------------------------------------------------------------------------

Pro stex_curve_gettrack
  common stex_common
  widget_control, stexbases[istex], get_uvalue=state
  widget_control, state.curveids[state.icurve], get_uvalue=str

  xphysize = str.xrange[1] - str.xrange[0]
  xdevsize = str.draw_window_size[0] $
             * (str.position[2] - str.position[0])
  xdev0 = str.draw_window_size[0] * str.position[0]
  str.mphys[0] = $
     (str.mouse[0] - xdev0) * xphysize / xdevsize + str.xrange[0]

  yphysize = str.yrange[1] - str.yrange[0]
  ydevsize = str.draw_window_size[1] $
             * (str.position[3] - str.position[1])
  ydev0 = str.draw_window_size[1] * str.position[1]
  str.mphys[1] = $
     (str.mouse[1] - ydev0) * yphysize / ydevsize + str.yrange[0]

  loc_string = strcompress( string(str.mphys[0], str.mphys[1]) )
  widget_control, str.location_bar_id, set_value=loc_String
  widget_control, state.curveids[state.icurve], set_uvalue=str
  return
end

;------------------------------------------------------------------------------

Pro stex_plot_curves
  common stex_common
  widget_control, stexbases[istex], get_uvalue=state

  pposition = !p.position
  !p.position = state.position

  for i=0, state.ncurves-1 do begin
     widget_control, state.curveids[i], get_uvalue=str

     wset, str.window_id
     widget_control, state.slider, get_value=iflux

     options = (*state.rawdata[0]).options
     c = where(tag_names(options) eq 'COLOR', ct)
     if (ct eq 1) then begin
        curvecolor=fsc_color(options.color)
        options.color = fsc_color('white')
     endif

     plot, [0], [0], /nodata, xrange=str.xrange, yrange=str.yrange, xstyle=1, $
           ystyle=1, _EXTRA=options

     for icurve=0, n_elements(str.curves)-1 do begin
        options = (*state.rawdata[icurve]).options
        c = where(tag_names(options) eq 'COLOR', ct)
        if (ct eq 1) then options.color = fsc_color(options.color)

        day = (*str.curves[icurve]).day
        curve = (*str.curves[icurve]).curve
        plot, day, curve, xrange=str.xrange, yrange=str.yrange, /noerase, $
              xstyle=5, ystyle=5, _EXTRA=options
        iday = (*state.dayptr)[iflux]
        if iday le max(day,/nan) and iday ge min(day,/nan) then begin
           w=(where(iday eq day))[0]
           plot, [day[w]], [curve[w]], /noerase, psym=2, symsize=2, $
                 xrange=str.xrange, yrange=str.yrange, xstyle=5, ystyle=5, $
                 color=options.color
           ;loc_string = strcompress( string(day[w], curve[w]))
           ;widget_control, str.location_id, set_value=loc_string
        endif
     endfor
  endfor
  !p.position = pposition
  return
end
;------------------------------------------------------------------------------

Pro stex_startup
  common stex_common

  ;load a simple color table
  ;loadct, 0, /silent
  ;if (!d.n_colors LE 16777216 AND !d.table_size LT 12) then $
  ;   message, 'Too few colors available for color table'
  ;if (!d.n_colors GT 256) then device, decomposed=0
  ;loadct, 0, /silent, bottom=8
  ;red=[0, 1, 0, 0, 0, 1, 1, 1]
  ;green = [0, 0, 1, 0, 1, 0, 1, 1]
  ;blue  = [0, 0, 0, 1, 1, 1, 0, 1]
  ;tvlct, 255*red, 255*green, 255*blue

  keylist = {key:' ', x:0.0, y:0.0}

  ;compare up to maxdata spectral templates
  maxdata = 10
  data = ptrarr(maxdata)
  rawdata = ptrarr(maxdata)
  dayptr = ptr_new()

  ;for lightcurve, color, k-correction plots...
  maxcurves = 50
  curveids = lonarr(maxcurves)


  state = { $
          version: '0.1a', $
          base_min_size: [600L, 400L], $
          draw_base_id: 0L, $
          draw_window_id: 0L, $
          draw_widget_id: 0L, $
          location_bar_id: 0L, $
          location_bar2_id: 0L, $
          xmin_text_id: 0L, $
          xmax_text_id: 0L, $
          mode_droplist_id: 0L, $
          ymin_text_id: 0L, $
          ymax_text_id: 0L, $
          keyboard_text_id: 0L, $
          norm_button:0L, $
          zoomone_button: 0L, $
          nkey: 0, $
          keylist: replicate(keylist,5), $
          xrange: [0.d0,1.d0], $
          yrange: [0.d0,1.d0], $
          yfix: 'fix', $
          position: [0.15d0, 0.10d0, 0.95d0, 0.95d0], $
          draw_window_size: [600L,512L], $
          mouse: [0L, 0L], $
          mphys: [0.0, 0.0], $
          base_pad: [0L, 0L], $
          pad: [0L, 0L],$
          slider: 0L, $
          ncurves: 0L, $
          curveIDs:curveIDs, $
          icurve:0L, $
          dayptr:dayptr, $
          data:data, $
          rawdata:rawdata, $
          ndata:1L}



  base = widget_base(title=strcompress('Spectral Template Explorer : ' $
                                       +string(nstex)), $
                     /column, /base_align_right, /tlb_size_events, $
                    /tlb_kill_request_events)

  ;store base id in common block
  stexbases[nstex] = base
  nstex += 1
  istex = nstex-1

  button_base1 =  widget_base(base, /row, /base_align_right)
  button_base2 =  widget_base(base, /row, /base_align_right)
  button_base3 =  widget_base(base, /row, /base_align_right)
  button_base4 =  widget_base(base, /row, /base_align_right)

  state.xmin_text_id = cw_field(button_base1, $
                                /float,  $
                                title = 'XMIN=', $
                                value = string(state.xrange[0]),  $
                                /return_events, $
                                xsize = 12)

  state.xmax_text_id = cw_field(button_base1, $
                                /float,  $
                                title = 'XMAX=', $
                                value = string(state.xrange[1]),  $
                                /return_events, $
                                xsize = 12)

  fixlist = ['Fix', 'Float']
  state.mode_droplist_id = widget_droplist(button_base2, $
                                           title= ' ', $
                                           uvalue = 'yfix', $
                                           value = fixlist)


  state.ymin_text_id = cw_field(button_base2, $
                                uvalue = 'ymin_text', /float,  $
                                title = 'YMIN=', $
                                value = string(state.yrange[0]),  $
                                /return_events, $
                                xsize = 12)

  state.ymax_text_id = cw_field(button_base2, $
                                uvalue = 'ymax_text', /float,  $
                                title = 'YMAX=', $
                                value = string(state.yrange[1]),  $
                                /return_events, $
                                xsize = 12)

  tmp_string = string('',format='(a30)')
  state.location_bar_id = widget_label(button_base1, $
                                       value = tmp_string,  $
                                       uvalue = 'location_bar', frame=1)

  state.slider = widget_slider(button_base2, /drag, $
                               uvalue = 'phase_slider', $
                               xsize = 123, $
                               ysize = 20, $
                               /suppress_value)

  tmp_string = string('',format='(a9)')
  state.location_bar2_id = widget_label(button_base2, $
                                       value = tmp_string,  $
                                       uvalue = 'location_bar2', frame=1)

  state.norm_button = widget_button(button_base3, $
                                    value = 'Normalize' )

  state.zoomone_button = widget_button(button_base3, $
                                       value='Zoom1' )

  state.keyboard_text_id =  widget_text(button_base4, $
                                        /all_events, $
                                        scr_xsize = 1, $
                                        scr_ysize = 1, $
                                        units = 0, $
                                        uvalue = 'keyboard_text', $
                                        value = '')

  state.draw_base_id = widget_base(base, $
                                   /column, /base_align_left, $
                                   /tracking_events, $
                                   uvalue = 'draw_base', $
                                   frame = 2)

  state.draw_widget_id = widget_draw(state.draw_base_id, $
                                     uvalue = 'draw_window', $
                                     /motion_events,  /button_events, $
                                     scr_xsize = state.draw_window_size[0], $
                                     scr_ysize = state.draw_window_size[1])

  widget_control, base, /realize

  widget_control, state.draw_widget_id, get_value = tmp_value
  state.draw_window_id = tmp_value

  basegeom = widget_info(base, /geometry)
  drawbasegeom = widget_info(state.draw_base_id, /geometry)
  state.pad[0] = basegeom.xsize - state.draw_window_size[0]
  state.pad[1] = basegeom.ysize - state.draw_window_size[1]
  state.base_pad[0] = basegeom.xsize - drawbasegeom.xsize $
                      + (2 * basegeom.margin)
  state.base_pad[1] = basegeom.ysize - drawbasegeom.ysize $
                      + (2 * basegeom.margin)

  widget_control, base, set_uvalue = state

  xmanager, strcompress('stex'), base, /no_block
  return
end

;------------------------------------------------------------------------------


Pro stex_event, event
  common stex_common
  ;if a stex window is closed, delete from common block
  if tag_names(event, /structure_name) eq 'WIDGET_KILL_REQUEST' then begin
     i = where(stexbases eq event.id)
     if i eq nstex-1 then begin
        stexbases[i]=0
        widget_control, event.id, /destroy
        nstex -= 1
     endif else begin
        stexbases[i:nstex-2]=stexbases[i+1:nstex-1]
        stexbases[nstex-1] = 0
        nstex -= 1
        widget_control, event.id, /destroy
     endelse
     if istex ge nstex then istex=nstex-1
     return
  endif

  ;now determine which instance of stex is active...
  istex=(where(event.top eq stexbases))[0]

  ;now we have to determine which sub-widget spawned the event...
  widget_control, event.top, get_uvalue=state
  case event.id of
     state.slider: begin
        stex_plot
        stex_plot_curves
        stex_slider_label
     end
     state.mode_droplist_id: begin
        case event.index of
           0: state.yfix = 'fix'
           1: state.yfix = 'float'
        endcase
        widget_control, event.top, set_uvalue=state
     end
     state.xmin_text_id: begin
        reads, event.value, tmp1
        state.xrange[0] = tmp1
        widget_control, event.top, set_uvalue=state
        stex_set_minmax
        stex_plot
     end
     state.xmax_text_id: begin
        reads, event.value, tmp1
        state.xrange[1] = tmp1
        widget_control, event.top, set_uvalue=state
        stex_set_minmax
        stex_plot
     end
     state.ymin_text_id: begin
        reads, event.value, tmp1
        state.yrange[0] = tmp1
        widget_control, event.top, set_uvalue=state
        stex_set_minmax
        stex_plot
     end
     state.ymax_text_id: begin
        reads, event.value, tmp1
        state.yrange[1] = tmp1
        widget_control, event.top, set_uvalue=state
        stex_set_minmax
        stex_plot
     end
     state.draw_widget_id: begin
        if (event.type eq 2) then begin ;motion event
           temp_event = [event.x, event.y]
           state.mouse = temp_event
           widget_control, event.top, set_uvalue=state
           stex_gettrack
        endif

        if (event.type eq 0) then begin
           case event.modifiers of
              1: xycase = 'yonly'
              4: xycase = 'xandy'
              else: xycase = 'xonly'
           endcase
           xycase = keyword_set(event.modifiers) ? 'xandy' : 'xonly'
           case event.press of
              1: stex_zoom, xycase+'-in', /recenter
              2: stex_zoom, xycase+'-recen', /recenter
              4: stex_zoom, xycase+'-out', /recenter
              else: print, 'trouble in stex_event, mouse zoom'
           endcase
        endif
        widget_control, state.keyboard_text_id, /input_focus
     end
     event.top: begin
        stex_resize, event
        stex_plot
     end
     state.draw_base_id: begin
        case event.enter of
           0: begin
              widget_control, state.keyboard_text_id, set_value=''
           end
           1: begin
              widget_control, state.keyboard_text_id, /input_focus
           end
        endcase
     end
     state.keyboard_text_id: begin
        eventchar = string(event.ch)

        if (state.nkey LT N_elements(state.keylist)) then begin
           state.keylist[state.nkey].key = eventchar
           state.keylist[state.nkey].x = state.mphys[0]
           state.keylist[state.nkey].y = state.mphys[1]
           state.nkey += 1
           widget_control, event.top, set_uvalue=state
        endif else begin
           stex_clearkeylist
        endelse

        case eventchar of
           '4' : stex_retract_phase
           '6' : stex_advance_phase
           'f' : stex_makefluxLC
           'm' : stex_makemagLC
           'c' : stex_makecolor
           'k' : stex_makeKcorrection
           'u' : stex_makemagLC, 'BessellUx'
           'b' : stex_makemagLC, 'BessellBx'
           'v' : stex_makemagLC, 'BessellV'
           'r' : stex_makemagLC, 'BessellR'
           'i' : stex_makemagLC, 'BessellI'
           else:
        endcase
        widget_control, state.keyboard_text_id, /clear_events
     end
     state.norm_button: begin
        stex_normalize_Bmax
     end
     state.zoomone_button: begin
        stex_auto_zoom
        stex_plot
     end
     else:
  endcase
  return
end

;------------------------------------------------------------------------------

Pro stex_zoom, zchange, recenter=recenter
  common stex_common
  widget_control, stexbases[istex], get_uvalue=state

   case zchange of
      'xonly-in': begin
         state.xrange = state.mphys[0] $
          + [-0.25, 0.25] * (state.xrange[1] - state.xrange[0])
         widget_control, stexbases[istex], set_uvalue=state
      end
      'yonly-in': begin
         state.yrange = state.mphys[1] $
          + [-0.25, 0.25] * (state.yrange[1] - state.yrange[0])
         widget_control, stexbases[istex], set_uvalue=state
      end
      'xandy-in': begin
         state.xrange = state.mphys[0] $
          + [-0.25, 0.25] * (state.xrange[1] - state.xrange[0])
         state.yrange = state.mphys[1] $
          + [-0.25, 0.25] * (state.yrange[1] - state.yrange[0])
         widget_control, stexbases[istex], set_uvalue=state
      end
      'xonly-out': begin
         state.xrange = state.mphys[0] $
          + [-1.0, 1.0] * (state.xrange[1] - state.xrange[0])
         widget_control, stexbases[istex], set_uvalue=state
      end
      'yonly-out': begin
         state.yrange = state.mphys[1] $
          + [-1.0, 1.0] * (state.yrange[1] - state.yrange[0])
         widget_control, stexbases[istex], set_uvalue=state
      end
      'xandy-out': begin
         state.xrange = state.mphys[0] $
          + [-1.0, 1.0] * (state.xrange[1] - state.xrange[0])
         state.yrange = state.mphys[1] $
          + [-1.0, 1.0] * (state.yrange[1] - state.yrange[0])
         widget_control, stexbases[istex], set_uvalue=state
      end
      'one': begin
         stex_auto_zoom
      end
      'xonly-recen': begin ; no change to zoom level: X recenter on mouse pos'n
         state.xrange = state.mphys[0] $
          + [-0.5, 0.5] * (state.xrange[1] - state.xrange[0])
         widget_control, stexbases[istex], set_uvalue=state
      end
      'yonly-recen': begin ; no change to zoom level: X recenter on mouse pos'n
         state.yrange = state.mphys[1] $
          + [-0.5, 0.5] * (state.yrange[1] - state.yrange[0])
         widget_control, stexbases[istex], set_uvalue=state
      end
      'xandy-recen': begin ; no change to zoom level: Y recenter on mouse pos'n
         state.xrange = state.mphys[0] $
          + [-0.5, 0.5] * (state.xrange[1] - state.xrange[0])
         state.yrange = state.mphys[1] $
          + [-0.5, 0.5] * (state.yrange[1] - state.yrange[0])
         widget_control, stexbases[istex], set_uvalue=state
      end
      else: print, 'Problem in splot_zoom!'
   endcase

   stex_set_minmax
   stex_plot
   return
end

;------------------------------------------------------------------------------

Pro stex_auto_zoom
  common stex_common
  widget_control, stexbases[istex], get_uvalue=state
  for idata=0, state.ndata-1 do begin
     ;; why [0] and not [idata]?
     xr = [min((*state.data[0]).wave,/nan), max((*state.data[0]).wave,/nan)]
     yr = [min((*state.data[0]).flux,/nan) < 0., max((*state.data[0]).flux, /nan)]
     if idata eq 0 then begin
        state.xrange = xr
        state.yrange = yr
     endif
     state.xrange[0] = state.xrange[0] < xr[0]
     state.xrange[1] = state.xrange[1] > xr[1]
     state.yrange[0] = state.yrange[0] < yr[0]
     state.yrange[1] = state.yrange[1] > yr[1]
  endfor
  widget_control, stexbases[istex], set_uvalue=state
  stex_set_minmax
  return
end

;------------------------------------------------------------------------------

Pro stex_resize, event
; routine to resize the draw window when a top-level resize event occurs.

  common stex_common
  widget_control, stexbases[istex], get_uvalue=state

  tmp_event = [event.x, event.y]
  window = (state.base_min_size > tmp_event)

  newbase = window - state.base_pad
  newsize = window - state.pad

  widget_control, state.draw_base_id, $
                  xsize = newbase[0], ysize = newbase[1]
  widget_control, state.draw_widget_id, $
                  xsize = newsize[0], ysize = newsize[1]

  state.draw_window_size = newsize
  widget_control, stexbases[istex], set_uvalue=state
  return
end

;------------------------------------------------------------------------------

Pro stex_gettrack
  common stex_common
  widget_control, stexbases[istex], get_uvalue=state

  xphysize = state.xrange[1] - state.xrange[0]
   xdevsize = state.draw_window_size[0] $
    * (state.position[2] - state.position[0])
   xdev0 = state.draw_window_size[0] * state.position[0]
   state.mphys[0] = $
    (state.mouse[0] - xdev0) * xphysize / xdevsize + state.xrange[0]

   yphysize = state.yrange[1] - state.yrange[0]
   ydevsize = state.draw_window_size[1] $
    * (state.position[3] - state.position[1])
   ydev0 = state.draw_window_size[1] * state.position[1]
   state.mphys[1] = $
    (state.mouse[1] - ydev0) * yphysize / ydevsize + state.yrange[0]

   loc_string = strcompress( string(state.mphys[0], state.mphys[1]) )
   widget_control, state.location_bar_id, set_value=loc_String
   widget_control, stexbases[istex], set_uvalue=state
   return
end


;------------------------------------------------------------------------------

Pro stex_slider_label
  common stex_common
  widget_control, stexbases[istex], get_uvalue=state
  widget_control, state.slider, get_value=iflux
  label = (*state.dayptr)[iflux]
  if label lt 1000 then $
     loc_string = strcompress( string(label, format='(f7.2)' ) ) $
  else $
     loc_string = strcompress( string(label, format='(e9.3)' ) )
  widget_control, state.location_bar2_id, set_value=loc_string
  return
end

;------------------------------------------------------------------------------

Pro stex_advance_phase
  common stex_common
  widget_control, stexbases[istex], get_uvalue=state

  widget_control, state.slider, get_value=phase
  day=(*state.data[0]).day
  if phase+1 le n_elements(day)-1 then $
     widget_control, state.slider, set_value=phase+1
  widget_control, stexbases[istex], set_uvalue=state
  stex_plot
  stex_plot_curves
  stex_slider_label
  return
end

;------------------------------------------------------------------------------

Pro stex_retract_phase
  common stex_common
  widget_control, stexbases[istex], get_uvalue=state

  widget_control, state.slider, get_value=phase
  if phase-1 ge 0 then $
     widget_control, state.slider, set_value=phase-1
  widget_control, stexbases[istex], set_uvalue=state
  stex_plot
  stex_plot_curves
  stex_slider_label
  return
end

;------------------------------------------------------------------------------

Pro stex_clearkeylist
  common stex_common
  widget_control, stexbases[istex], get_uvalue=state
  state.nkey=0
  state.keylist.key = ' '
  state.keylist.x = 0.0
  state.keylist.y = 0.0
  widget_control, stexbases[istex], set_uvalue=state

  return
end

;------------------------------------------------------------------------------

Pro stex_plot
  common stex_common

  widget_control, stexbases[istex], get_uvalue=state
  wset, state.draw_window_id

  pposition = !p.position
  !p.position = state.position
  widget_control, state.slider, get_value=iflux

  if state.ndata lt 1 then begin
     plot, [0], [0], /nodata, xrange=state.xrange, yrange=state.yrange, $
           ystyle=1, xstyle=1
     return
  endif

  if state.yfix eq 'float' then begin
     state.yrange = [!values.f_infinity,-!values.f_infinity]
     for idata=0, state.ndata-1 do begin
        wave=(*state.data[idata]).wave
        day=(*state.data[idata]).day
        flux=(*state.data[idata]).flux

        if (*state.dayptr)[iflux] gt max(day) $
           or (*state.dayptr)[iflux] lt min(day,/nan) then continue
        jflux = (where((*state.dayptr)[iflux] eq day))[0]

        w=where(wave ge state.xrange[0] and wave le state.xrange[1])
        extra = 0.1*(max(flux[w,jflux[0]], /nan)-min(flux[w,jflux], /nan))
        yr=dblarr(2)
        yr[0] = min(flux[w,jflux], /nan)-extra < 0.
        yr[1] = max(flux[w,jflux], /nan)+extra
        state.yrange[0] = state.yrange[0] < yr[0]
        state.yrange[1] = state.yrange[1] > yr[1]
     endfor
     widget_control, stexbases[istex], set_uvalue=state
     stex_set_minmax
  endif

  options = (*state.rawdata[0]).options
  c = where(tag_names(options) eq 'COLOR', ct)
  if (ct eq 1) then options.color = fsc_color('white')

  plot, [0], [0], /nodata, xrange=state.xrange, yrange=state.yrange, xstyle=1, $
        ystyle=1, _EXTRA=options

  for idata=0, state.ndata-1 do begin
     wave=(*state.data[idata]).wave
     day=(*state.data[idata]).day
     flux=(*state.data[idata]).flux
     options=(*state.rawdata[idata]).options

     if (*state.dayptr)[iflux] gt max(day,/nan) $
        or (*state.dayptr)[iflux] lt min(day,/nan) then continue
     jflux = (where((*state.dayptr)[iflux] eq day))[0]

     c = where(tag_names(options) eq 'COLOR', ct)
     if (ct eq 1) then options.color = fsc_color(options.color)

     plot, wave, flux[*,jflux], xrange=state.xrange, yrange=state.yrange, $
           /noerase, xstyle=5, ystyle=5, _EXTRA=options
  endfor
  !p.position = pposition
  return
end

;------------------------------------------------------------------------------

Pro stex_set_minmax
  common stex_common
  widget_control, stexbases[istex], get_uvalue=state

  widget_control, state.xmin_text_id, set_value = state.xrange[0]
  widget_control, state.xmax_text_id, set_value = state.xrange[1]
  widget_control, state.ymin_text_id, set_value = state.yrange[0]
  widget_control, state.ymax_text_id, set_value = state.yrange[1]

  return
end

;------------------------------------------------------------------------------
; stex_interpolate_data
;
; Because we can only display one age at a time, we interpolate all of
; the spectral templates to a common age array.  The age array
; consists of the union of all of the input age arrays.

Pro stex_interpolate_data
  common stex_common
  widget_control, stexbases[istex], get_uvalue=state

  ;determine the union of the day datasets
  for idata=0, state.ndata-1 do begin
     if idata eq 0 then dayunion=(*state.rawdata[idata]).day $
     else dayunion=[dayunion,(*state.rawdata[idata]).day]
  endfor
  dayunion=(dayunion[sort(dayunion)])[uniq(dayunion[sort(dayunion)])]

  ;now interpolate onto these dates...
  for idata=0, state.ndata-1 do begin

     wdayout = where(dayunion ge min((*state.rawdata[idata]).day,/nan) and $
                     dayunion le max((*state.rawdata[idata]).day,/nan), ndayout)
     dayout = dayunion[wdayout]
     nwave = n_elements((*state.rawdata[idata]).wave)

     flux = dblarr(nwave, ndayout)

     progressBar = Obj_New("PROGRESSBAR", Color='red', $
                           Text='data '+strtrim(idata+1,2)+' of '+strtrim(state.ndata,2), $
                           /fast_loop, /nocancel)
     progressBar -> Start
     for iwave=0, nwave-1 do begin
        progressBar -> Update, 100.*iwave/(nwave-1)
        flux[iwave,*] = interpol((*state.rawdata[idata]).flux[iwave,*], $
                                 (*state.rawdata[idata]).day, $
                                 dayout, /spline)

     endfor
     progressBar -> Destroy

     state.dayptr=ptr_new(dayunion)
     wave = (*state.rawdata[idata]).wave
     data={day:dayout, wave:wave, flux:flux}
     state.data[idata] = ptr_new(data)
  endfor
  widget_control, state.slider, $
                  set_slider_min=0, $
                  set_slider_max=(n_elements(dayunion)-1), $
                  set_value=(n_elements(dayunion)/2.)
  widget_control, stexbases[istex], set_uvalue=state
  stex_slider_label
  return
end

;------------------------------------------------------------------------------

Pro stex_interpolate_data2
  common stex_common
  widget_control, stexbases[istex], get_uvalue=state

  ;;take the union of all 'day' arrays
  for idata=0, state.ndata-1 do begin
     if idata eq 0 then dayunion=(*state.rawdata[idata]).day $
     else dayunion = [dayunion,(*state.rawdata[idata]).day]
  endfor
  dayunion=(dayunion[sort(dayunion)])[uniq(dayunion[sort(dayunion)])]

  ;now start filling in the flux...
  for idata=0, state.ndata-1 do begin
     nwave = n_elements((*state.rawdata[idata]).wave)
     wdayout = where(dayunion ge min((*state.rawdata[idata]).day) and $
                     dayunion le max((*state.rawdata[idata]).day), ndayout)
     dayout = dayunion[wdayout]
     flux = dblarr(nwave, ndayout)

     donedays = where(member(dayout, (*state.rawdata[idata]).day), complement=needdays)
     flux[*,donedays] = (*state.rawdata[idata]).flux
     ; now interpolate the rest...
     if needdays[0] ne -1 then begin
        progressBar = Obj_New("PROGRESSBAR", Color='red', $
                              Text='data '+strtrim(idata+1,2)+' of '+strtrim(state.ndata,2), $
                              /fast_loop, /nocancel)
        progressBar -> Start
        for iwave=0, nwave-1 do begin
           flux[iwave,needdays] = interpol((*state.rawdata[idata]).flux[iwave,*], $
                                           (*state.rawdata[idata]).day, $
                                           dayout[needdays])
           progressBar -> Update, 100.*iwave/(nwave-1)
        endfor
        progressBar -> Destroy
     endif

     ptr_free, state.dayptr
     state.dayptr=ptr_new(dayunion)
     wave = (*state.rawdata[idata]).wave
     data = {day:dayout, wave:wave, flux:flux}
     state.data[idata] = ptr_new(data)
  endfor
  widget_control, state.slider, $
                  set_slider_min=0, $
                  set_slider_max=(n_elements(dayunion)-1), $
                  set_value=(n_elements(dayunion)/2.)
  widget_control, stexbases[istex], set_uvalue=state
  stex_slider_label
  return
end

;------------------------------------------------------------------------------

Pro stex_normalize_Bmax
  common stex_common
  widget_control, /hourglass
  widget_control, stexbases[istex], get_uvalue=state

  Bmax = dblarr(state.ndata)
  for idata=0, state.ndata-1 do begin
     wave=(*state.rawdata[idata]).wave
     day=(*state.rawdata[idata]).day
     flux=(*state.rawdata[idata]).flux
     filter = make_filter(wave, 'BessellBx')
     filter = rebin(filter, n_elements(wave), n_elements(day))

     dlambda = wave-shift(wave,1)
     dlambda[0] = wave[1]-wave[0]
     dlambda = rebin(dlambda, n_elements(wave), n_elements(day))

     product = filter*flux*dlambda*rebin(wave, n_elements(wave), n_elements(day))

     LCflux = total(product,1)

     ;Bmax[idata] = interpol(LCflux, day, 0, /spline)
     Bmax[idata] = max(LCflux)
  endfor

  for idata=0, state.ndata-1 do begin
     rawdata = *state.rawdata[idata]
     rawdata.flux *= Bmax[0]/Bmax[idata]
     ptr_free, state.rawdata[idata]
     state.rawdata[idata] = ptr_new(rawdata)

     data = *state.data[idata]
     data.flux *= Bmax[0]/Bmax[idata]
     ptr_free, state.data[idata]
     state.data[idata] = ptr_new(data)
  endfor
  widget_control, stexbases[istex], set_uvalue=state
  stex_plot
  stex_plot_curves
  return
end

;------------------------------------------------------------------------------

Pro sterase, nerase
  common stex_common
  widget_control, stexbases[istex], get_uvalue=state

  ndata = (state.ndata)

  if (N_params() lt 1) then nerase = ndata $
  else if (nerase gt ndata) then nerase = ndata

  for idata=ndata-nerase, ndata-1 do begin
     ptr_free, state.data[idata]
     ptr_free, state.rawdata[idata]

  endfor
  state.ndata -= nerase

  widget_control, stexbases[istex], set_uvalue=state

  stex_plot
  stex_plot_curves
  return
end

;------------------------------------------------------------------------------

Pro ostex, d, w, f, istex=is, _EXTRA=options
  common stex_common
  widget_control, /hourglass

  if xregistered('stex') eq 0 then $
     message, 'You need to start stex first!'

  if n_elements(is) ne 0 then begin
     if is gt nstex-1 then message, 'Not that many stexes!!' $
     else istex=is
  endif

  widget_control, stexbases[istex], get_uvalue=state

  if (n_elements(options) eq 0) then $
     options = {color:'white'}

  state.ndata+=1
  rawdata={day:double(d), wave:double(w), flux:double(f), options:options}
  state.rawdata[state.ndata-1]=ptr_new(rawdata)

  widget_control, stexbases[istex], set_uvalue=state

  stex_interpolate_data2

  stex_plot

  return
end

;------------------------------------------------------------------------------

Pro stex, d, w, f, _EXTRA=options
  common stex_common
  if n_elements(stexbases) eq 0 then begin
     nstex=0L
     istex=-1L
     stexbases=lonarr(10)
  endif

  if (N_params() LT 3) then begin
     print, 'Too few parameters for STEX'
     return
  endif

  stex_startup
  widget_control, stexbases[istex], get_uvalue=state

  if (N_elements(options) eq 0) then $
     options = {color:'white'}

  data = {day:double(d), wave:double(w), flux:double(f),options:options}
  rawdata = data
  state.data[0]=ptr_new(data)
  state.rawdata[0]=ptr_new(rawdata)
  state.dayptr = ptr_new(d)

  widget_control, stexbases[istex], set_uvalue=state
  widget_control, state.slider, $
                  set_slider_min=0, $
                  set_slider_max=(n_elements(d)-1), $
                  set_value=(n_elements(d)/2.)
  stex_slider_label
  stex_auto_zoom

  stex_plot

  return
end
