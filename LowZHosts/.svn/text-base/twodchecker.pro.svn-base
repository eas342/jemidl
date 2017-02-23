Pro twodchecker, snname, setup, i_indiv
  dir='~/scp4/LowZHosts/JMreduction/'
  setups = ['B460','B560','R5000','R7500','R1000B','R2000B','R2500R']
  longsetups = ['blue600_4000_d460', $
                'blue600_4000_d560', $
                'red600_5000', $
                'red600_7500', $
                'R1000B', $
                'R2000B', $
                'R2500R']
  w=(where(setup eq setups))[0]
  checkfile=file_search(dir+'*/'+longsetups[w]+'/FinalCR/'+snname+'_check.fits')
  fluxfile=file_search(dir+'*/'+longsetups[w]+'/FinalCR/'+snname+'_F.fits')
  ck=mrdfits(checkfile, 1, /silent)
  n_indiv = (size(ck.indiv, /dim))[0]
  infiles=strarr(n_indiv)
  for i=0, n_indiv-1 do begin
     infiles[i] = sxpar(headfits(fluxfile), 'INFILE'+strn(i,F='(I02)'))
  endfor
  split = strsplit(checkfile, '/', /extract)
  dir=strjoin(split[0:n_elements(split)-3], '/')

  filename='/'+dir+'/Science_box/'+infiles[i_indiv]
  filename = repstr(filename, '.newwave', '')
  sci=mrdfits(filename, 0, /silent)
  sky=mrdfits(filename, 2, /silent)
  print, infiles[i_indiv]
  atv, sci-sky
end
