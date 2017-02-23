Pro preparedir, dir, version, epoch=epoch, outdir=outdir
  if n_elements(epoch) ne 0 then epochstr='*'+str(epoch,f='(I02)') $
  else epochstr=''
  files = file_search(dir+'CL-*_v'+str(version)+epochstr+'*_drz.fits')
  for ifile=0, n_elements(files)-1 do begin
     preparefile, files[ifile], outdir=outdir
  endfor
end
