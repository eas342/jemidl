; This procedure executes sextractor on the i/z files given using the sextractor
; configuration file SEXFILE.  Output catalogs are prefixed with PREFIX.
; zeropoints are for 77K AB magnitudes and are taken from the web.

Pro sextract_iz, ifile, zfile, sexfile, prefix, outdir=outdir, segment=segment
  cd, current=pwd
  if n_elements(outdir) eq 0 then outdir=pwd
  outdir=directoryify(outdir)
  ;first get the correct exposure times...
  z_h=headfits(ifile)
  i_h=headfits(zfile)
  z_exptime = sxpar( z_h, 'EXPTIME' )
  i_exptime = sxpar( i_h, 'EXPTIME' )

  file_mkdir, outdir+prefix+'tmp/'

  ;now copy appropriate extensions to temp files...
  spawn, 'imcopy "'+ifile+'[1]" '+outdir+prefix+'tmp/ifile.fits'
  spawn, 'imcopy "'+zfile+'[1]" '+outdir+prefix+'tmp/zfile.fits'
  spawn, 'imcopy "'+ifile+'[2]" '+outdir+prefix+'tmp/iweight.fits'
  spawn, 'imcopy "'+zfile+'[2]" '+outdir+prefix+'tmp/zweight.fits'

  ;now do the z-band sextraction
  cmd = 'sex '+outdir+prefix+'tmp/zfile.fits -c ' + sexfile $
        + ' -MAG_ZEROPOINT 24.867' $
        + ' -GAIN '+string(z_exptime) $
        + ' -WEIGHT_IMAGE '+outdir+prefix+'tmp/zweight.fits' $
        + ' -WEIGHT_TYPE MAP_WEIGHT' $
        + ' -CATALOG_NAME ' + outdir+prefix+'_z.cat'
  if keyword_set(segment) then cmd += ' -CHECKIMAGE_TYPE SEGMENTATION' $
                                      + ' -CHECKIMAGE_NAME ' + outdir+prefix+'segment.fits'
  print, cmd
  spawn, cmd

  ;now do the i-band sextraction
  cmd = 'sex '+outdir+prefix+'tmp/zfile.fits,' $
        + outdir+prefix+'tmp/ifile.fits -c ' + sexfile $
        + ' -MAG_ZEROPOINT 25.678' $
        + ' -GAIN '+string(i_exptime) $
        + ' -WEIGHT_IMAGE '+outdir+prefix+'tmp/zweight.fits,' $
        + outdir+prefix+'tmp/iweight.fits' $
        + ' -WEIGHT_TYPE MAP_WEIGHT' $
        + ' -CATALOG_NAME ' + outdir+prefix+'_i.cat'
  print, cmd
  spawn, cmd
  spawn, 'rm -rf '+outdir+prefix+'tmp/'

end
