Pro preparefile, filename, outdir=outdir
  if n_elements(outdir) eq 0 then cd, current=outdir ;;set an output directory
  outdir = directoryify(outdir)
  file_mkdir, outdir
  filesplit = strsplit(filename, '/', /extract)  ;;split fullfilename into path+filename
  if n_elements(filesplit) ge 2 then begin
     dir = '/'+strjoin(filesplit[0:n_elements(filesplit)-2],'/')+'/'
     file = filesplit[n_elements(filesplit)-1]
  endif else begin
     dir = './'
     file = filesplit[0]
  endelse
  fileroot = strmid(file, 0, strlen(file)-5)

  ;; this is a hack!
  ostr=strpos(fileroot, 'CL-O') ;;change 'O' to 'D'
  if ostr ne -1 then fileroot = strmid(fileroot, 0, ostr)+'CL-D'+strmid(fileroot,ostr+4)
  ;; end hack!

  print, dir
  print, file
  print, fileroot

  file_delete, outdir+fileroot+'_SCI.fits', /allow_nonexistent  ;;imcopy sci and wht extensions
  file_delete, outdir+fileroot+'_WHT.fits', /allow_nonexistent
  print, 'copying multidrizzle science extension'
  spawn, 'imcopy "'+filename+'[1]" '+outdir+fileroot+'_SCI.fits'
  print, 'copying multidrizzle weight extension'
  spawn, 'imcopy "'+filename+'[2]" '+outdir+fileroot+'_WHT.fits'

  head = headfits( filename )  ;;construct sky image
  print, 'reading in science image'
  image = mrdfits( filename, 1, /silent )
  print, 'reading in effective exposure time image'
  time  = mrdfits( filename, 2, /silent )
  print, 'reading in context exposures image'
  ctx   = mrdfits( filename, 3, /silent )

  print, 'computing effective number of exposures'
  nexp  = nexp( ctx )
  mwrfits, nexp, outdir+fileroot+'_NXP.fits', head, /create

  print, 'computing estimated sky image'
  fltdir = '/home/scpdata01/acs/acs_wcsfixed3/'
  sky = float( skyim( ctx, head, fltdir, looksubdirs='**/' ) )
  mwrfits, sky, outdir+fileroot+'_SKY.fits', head, /create

  gain  = 1.0
  rdnoise = 5.0
  dark  = 0.0

  print, 'computing estimated error image'
  error = float( errim( image, time, gain, rdnoise, sky, dark, nexp ) )
  mwrfits, error, outdir+fileroot+'_ERR.fits', head, /create

  print, 'gzipping results'
  spawn, 'gzip -f '+outdir+fileroot+'_SCI.fits'
  spawn, 'gzip -f '+outdir+fileroot+'_WHT.fits'
  spawn, 'gzip -f '+outdir+fileroot+'_SKY.fits'
  spawn, 'gzip -f '+outdir+fileroot+'_ERR.fits'
  spawn, 'gzip -f '+outdir+fileroot+'_NXP.fits'
end
