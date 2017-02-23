Function readWaveCal, fileName
  MAXOUT = 20L
  openr, lun, fileName, /get_lun
  outProto = {isBox:0, $
              isOpt:0, $
              iObj:0, $
              waveCoeffs:ptr_new(), $
              fluxCoeffs:ptr_new(), $
              waveRange:dblarr(2) }
  outs = replicate(outProto, MAXOUT)
  i = 0
  while ~eof(lun) do begin
     out1 = outProto
     a=''
     readf, lun, a
     a=strsplit(a, ' ', /extract)
     if a[1] eq 'BOX' then out1.isBox = 1
     if a[1] eq 'OPT' then out1.isOpt = 1
     out1.iObj = fix(a[2])
     a=''
     readf, lun, a
     a=strsplit(a, ' ', /extract)
     waveCoeffs = dblarr(fix(a[0]))
     for j=0, fix(a[0])-1 do begin
        waveCoeffs[j] = double(a[j+1])
     endfor
     waveRange = dblarr(2)
     waveRange[0] = double(a[j+1])
     waveRange[1] = double(a[j+2])
     a=''
     readf, lun, a
     a=strsplit(a, ' ', /extract)
     if fix(a[0]) gt 0 then begin
        fluxCoeffs = dblarr(fix(a[0]))
        for j=0, fix(a[0])-1 do begin
           fluxCoeffs[j] = double(a[j+1])
        endfor
        out1.fluxCoeffs = ptr_new(fluxCoeffs)
     endif
     out1.waveCoeffs = ptr_new(waveCoeffs)
     out1.waveRange = waveRange
     outs[i] = out1
     i += 1
  endwhile
  close, lun
  free_lun, lun
  return, outs[0:i-1]
end
