Function zcheck
  out1={ SNname:'', $
         reduction:'', $
         setup:'', $
         MJD:0.d, $
         vdisp:0.d, $
         vdisp_err:0.d, $
         z:0.d, $
         z_err:0.d $
       }
  search='~/scp4/LowZHosts/GGreduction/*/*/*/Final/*galstruct.fits'
  f=file_search(search, count=nfiles)
  out = replicate(out1, nfiles)
  out.reduction='GG'
  for i=0, nfiles-1 do begin
     split=strsplit(f[i],'/',/extract)
     out[i].setup=split[7]
     out[i].SNname=repstr(split[9], '_galstruct.fits', '')
     h=headfits(repstr(f[i], '_galstruct.fits', '_F.fits'))
     out[i].mjd=double(sxpar(h,'MJD-OBS'))
     a=mrdfits(f[i], 1, /silent)
     out[i].vdisp=a[0].vdisp
     out[i].vdisp_err=a[0].vdisp_err
     out[i].z=a[0].z
     out[i].z_err=a[0].z_err
  endfor

  search='~/scp4/LowZHosts/JMreduction/*/*/Final/*galstruct.fits'
  f=file_search(search, count=nfiles2)
  out2 = replicate(out1, nfiles2)
  out2.reduction='JM'
  for i=0, nfiles2-1 do begin
     split=strsplit(f[i],'/',/extract)
     out2[i].setup=split[6]
     out2[i].SNname=repstr(split[8], '_galstruct.fits', '')
     h=headfits(repstr(f[i], '_galstruct.fits', '_F.fits'))
     out2[i].mjd=double(sxpar(h,'MJD-OBS'))
     a=mrdfits(f[i], 1, /silent)
     out2[i].vdisp=a[0].vdisp
     out2[i].vdisp_err=a[0].vdisp_err
     out2[i].z=a[0].z
     out2[i].z_err=a[0].z_err
  endfor
  out=[out,out2]
  out.setup = repstr(out.setup, 'red600_5000', 'R5000')
  out.setup = repstr(out.setup, 'red600_7500', 'R7500')
  out.setup = repstr(out.setup, 'blue600_4000_d460', 'B460')
  out.setup = repstr(out.setup, 'blue600_4000_d560', 'B560')

  out.snname = repstr(out.snname, 'Ellip2', 'SNF20080614-010')
  out.snname = repstr(out.snname, 'Ellip3', 'SNF20060521-008')
  out.snname = repstr(out.snname, 'Ellip4', 'SNF20070417-002')
  out.snname = repstr(out.snname, 'Ellip5', 'SNF20070727-016')
  out.snname = repstr(out.snname, 'Ellip6', 'SNF20080731-000')
  out.snname = repstr(out.snname, 'SNF20080512', 'SNF20080512-010')
  out.snname = repstr(out.snname, 'SNF20080512-010-010', 'SNF20080512-010')
  out.snname = repstr(out.snname, 'sn', 'SN')
  w=where(strmid(out.snname,0,2) ne 'SN')
  out[w].snname='SN'+out[w].snname
  s=sort(out.snname)
  u=out[uniq(out.snname, s)].snname
  for i=0, n_elements(u)-1 do begin
     print, u[i]
     for j=0, n_elements(s)-1 do begin
        if out[j].snname ne u[i] then continue
        print, out[j].reduction, out[j].setup, out[j].z, out[j].vdisp, $
               format='(A4,A7,F9.5,F12.5)'
     endfor
  endfor

  return, out
end
