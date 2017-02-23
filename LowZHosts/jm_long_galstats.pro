;+
; NAME:
;   long_galstats
;
; PURPOSE:
;   Compute galaxy stats for LRIS galaxy spectra.
;
; CALLING SEQUENCE:
;   res   =  long_galstats(ffile,resfile,inter=inter)
;
; INPUTS:
;   ffile      - Fluxed galaxy spectrum.
;   resfile    - Resolution file generated by long_resvec.pro
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   res      - Output structure with redshift and velocity dispersion
;              data
;
;   We currently mask within +/- 280 km/sec of the following wavelengths
;   that could have emission lines:
;     linelist = [3725.94, 3727.24, 3970.072, 4101.73, 4340.46, $
;      4861.3632, 4958.911, 5006.843, 6300.32, 6548.05, 6562.801, $
;      6583.45, 6716.44, 6730.82]
;
;   The constructed over-sampled and smoothed eigenspectra are stored
;   in a common block between calls if the eigen-vector file is the same.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/templates/spEigenVdisp*.fits
;
; PROCEDURES CALLED:
;   long_vdispfit
;   zfind
;   airtovac
;   combine1fiber
;   computechi2()
;   djs_filepath()
;   find_nminima
;   mrdfits()
;   poly_array()
;   splog
;   sxpar()
;
; REVISION HISTORY:
;   11-Nov-2009  Written by J. Hennawi based on SDSS spreduce1d.pro
;------------------------------------------------------------------------------
; Create output structure
PRO MASK_TELLURIC, loglam, ivar
  telpix = WHERE(10.0^loglam GE 7585.8D AND 10.0^loglam LE 7703.0D, ntel)
  IF ntel NE 0 THEN ivar[telpix] = 0.0D
  telpix = WHERE(10.0^loglam GE 6862.1D AND 10.0^loglam LE 6964.6D, ntel)
  IF ntel NE 0 THEN ivar[telpix] = 0.0D
  ;;telpix = WHERE(10.0^loglam GE 7143.3D AND 10.0^loglam LE 7398.2D, ntel)
  ;;IF ntel NE 0 THEN ivar[telpix] = 0.0D
  ;; TESTING masking
  ;;telpix = WHERE(10.0^loglam GE 5885.0D AND 10.0^loglam LE 5900.0D, ntel)
  ;;IF ntel NE 0 THEN ivar[telpix] = 0.0D
  ;;telpix = WHERE(10.0^loglam GE 5573.0D AND 10.0^loglam LE 5585.0D, ntel)
  ;;IF ntel NE 0 THEN ivar[telpix] = 0.0D
  ;;telpix = WHERE(10.0^loglam LE 4800, ntel)
  ;;IF ntel NE 0 THEN ivar[telpix] = 0.0D

RETURN
END

FUNCTION JM_LONG_GALSTATS, influx, inivar, resvec1, inloglam $
                         , TELLURIC = TELLURIC, INTER = INTER $
                         , SPLINEMASK = SPLINEMASK, NRES = NRES $
                         , VDISPMODEL = VDISPMODEL, RESET=reset $
                         , HARDCOPY = HARDCOPY


   IF KEYWORD_SET(INTER) THEN debug = 1
   IF (keyword_set(debug)) THEN doplot = 1
   nobj = 1

   clight = 2.99792458e5
   min_log = min(inloglam)
   max_log = max(inloglam)
   npixobj = n_elements(inloglam)
   COEFF0 = min_log
   COEFF1 = (max_log-min_log)/double(npixobj-1L)
   med_res = djs_median(resvec1)
   ;; TESTING!!!
   ;;COEFF1 = dvpix_in/alog(10.0d)/clight
   ;;COEFF1 = 1d-4 ;; SDSS grid for TESTING only
   ;;npixobj = ceil((max_log-min_log)/COEFF1)

   ;; Mask Telluric absorption
   IF KEYWORD_SET(TELLURIC) THEN mask_telluric, inloglam, inivar
   IF KEYWORD_SET(SPLINEMASK) THEN BEGIN
      IF NOT KEYWORD_SET(NRES) THEN NRES = 1.0
      bkspace = NRES*(med_res/alog(10.0d)/clight)
      spec_set = bspline_iterfit(inloglam, influx $
                                 , invvar = inivar*(inivar GT 0) $
                                 , bkspace = bkspace $
                                 , yfit = flux_spline, maxiter = 10L $
                                 , upper = 7.0d, lower = 7.0d, nord = 3 $
                                 , maxrej = 10, outmask = outmask1 $
                                 , /silent, /sticky)
      ;; Grow the mask by two pixels
      outmask = long_grow_mask(outmask1, 2)
      ind_bad  = WHERE(outmask EQ 0 OR inivar EQ 0, nrej)
      IF nrej GT 0 THEN inivar[ind_bad] = 0
      IF KEYWORD_SET(INTER) THEN BEGIN
         var = (inivar GT 0.0)/(inivar + (inivar EQ 0.0))
         sig = (inivar GT 0.0)/sqrt((var + (var EQ 0.0)))
         x_splot, 10.0D^inloglam, influx, psym1 = 10 $
                  , ytwo = sig, psym2 = 10 $
                  , xthr = 10.0D^inloglam[ind_bad] $
                  , ythr = influx[ind_bad], psym3 = 6 $
                  , /block, title = 'Masked pixels' $
                  , yfou = flux_spline $
                  , XMNX = 10.0D^[min_log, max_log] $
                  , YMNX = [-0.1, max(flux_spline)*1.2]
      ENDIF
   ENDIF

   ;; Make a fake hdr for zfind
   mkhdr, hdr, 0, /image
   sxaddpar, hdr, 'COEFF0', COEFF0
   sxaddpar, hdr, 'COEFF1', COEFF1
   objloglam0 = coeff0
   objdloglam = coeff1
   objloglam = objloglam0 + lindgen(npixobj) * objdloglam

   objflux = interpol(influx, inloglam, objloglam)
   objivar = interpol(inivar, inloglam, objloglam)
   resvec  = interpol(resvec1, inloglam, objloglam)
   mask    = interpol(double(inivar GT 0.0), inloglam, objloglam)
   ;; This masks pixels in the interpolated ivar which were previously masked
   objivar = objivar*(mask GT 0.5D)
   ;; Now compute the disp Sigma of the intrumental response
   ;; (dispersion) in units of the km/s
   sigres = resvec/2.35482D
   ;; Bin template onto same grid
   if (n_elements(eigendir) EQ 0) then $
      eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')
   eigenfile = 'spEigenGal-*.fits'

   allfiles = findfile(djs_filepath(eigenfile, root_dir = eigendir), count = ct)
   if (ct EQ 0) then $
      message, 'Unable to find EIGENFILE matching '+eigenfile
   thisfile = allfiles[ (reverse(sort(allfiles)))[0] ]
   splog, 'Selecting EIGENFILE=' + thisfile
   if (keyword_set(columns)) then $
      splog, 'Selecting columns=', columns
   ;;----------
   ;; Read the template file, and optionally trim to only those columns
   ;; specified by COLUMNS.
   ;; Assume that the wavelength binning is the same as for the objects
   ;; in log-wavelength.
   starflux1 = readfits(thisfile, shdr, /silent)
   if n_elements(starflux1) LE 1 then begin
      splog, 'Looking for a bspline structure ', thisfile
      bspline_set = mrdfits(thisfile, 1, shdr, /silent)
   endif else begin
      starloglam0 = sxpar(shdr, 'COEFF0')
      stardloglam0 = sxpar(shdr, 'COEFF1')
   endelse
   dims = size(starflux1, /dim)
   npixstar = dims[0]
   starloglam = starloglam0 + lindgen(npixstar)*stardloglam0

   min_star = min(starloglam)
   max_star = max(starloglam)
   npixnew = ceil((max_star-min_star)/objdloglam)

   starflux = fltarr(npixnew, dims[1])
   starnewloglam = starloglam0 + lindgen(npixnew)*objdloglam
   FOR ii = 0L, dims[1]-1L DO $
      starflux[*, ii] = interpol(starflux1[*, ii], starloglam, starnewloglam)

   npoly = 3
   zmin = 0.00 ; 0 km/sec
   zmax = 0.25 ; Max z
   pspace = 2
   nfind = 5
   plottitle = 'Galaxy Redshift'

   splog, 'Compute GALAXY redshifts:', $
    ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
   t0 = systime(1)
   res_gal = zfind(objflux, objivar, hdr = hdr $
                   , npoly = npoly $
                   , zmin = zmin, zmax = zmax, pspace = pspace $
                   , nfind = nfind, width = 5*pspace $
                   , plottitle = plottitle, doplot = doplot $
                   , debug = debug, /verbose $
                   , starflux = starflux, starloglam0 = starloglam0)
   splog, 'CPU time to compute GALAXY redshifts = ', systime(1)-t0

   splog, 'Locally re-fitting GALAXY redshifts'
   t0 = systime(1)
   res_gal = zrefind(objflux, objivar, hdr=hdr $
                     , pwidth = 5, pspace = 1, width = 5, zold = res_gal $
                     , plottitle = plottitle, doplot = doplot, debug = debug $
                     , starflux = starflux, starloglam0 = starloglam0)

   nstruct = n_elements(res_gal)
   res_gal.tfile = strmid(thisfile,strpos(thisfile,'/spE')+1)
   proto = create_struct('VDISP_TOT', 0.0)
   res_gal = struct_addtags(res_gal,  replicate(proto, nstruct))
   splog, 'CPU time to re-fit GALAXY redshifts = ', systime(1)-t0
   ; Only solve for velocity dispersions for the best-fit
   splog, 'Find velocity dispersions for galaxies'
   t0 = systime(1)
   ifind = 0
   if keyword_set(HARDCOPY) then begin
      set_plot, 'ps'
      device, filename=HARDCOPY+'.vdisp.eps'
      device, /encapsulated, /color, bits_per_pixel=8
      device, xsize=10, ysize=8, /inches, xoffset=0, yoffset=0
      device, /times, /isolatin1
      !p.font=0
      bangp=!p
      debug=1
   endif

;   vdans = long_vdispfit(objflux, objivar, hdr = hdr $
;                         , zobj = res_gal[ifind, *].z $
;                         , eigenfile = 'spEigenElodie.fits' $
;                         , columns = lindgen(24), yfit = dispflux $
;                         , debug = debug, reset=reset)
   vdans = long_vdispfit(objflux, objivar, hdr = hdr $
                         , zobj = res_gal[ifind, *].z $
                         , eigenfile = 'spEigenElodie.fits' $
                         , columns = lindgen(24), yfit = dispflux $
                         , debug = 0, reset=reset)
   if n_elements(HARDCOPY) ne 0 then begin
      !p=bangp
      device, /close
   endif

   VDISPMODEL = interpol(dispflux, objloglam, inloglam)
   ;; No longer pass sigres to the long_vdispfit, so it will not
   ;; smooth to LRIS resolution. If you require true velocity
   ;; dispersion, you'll need to fix whatever is wrong with
   ;; our resolution model.
;;                         , sigres = sigres

   res_gal[ifind, *].vdisp = reform([vdans.vdisp], 1, nobj)
   res_gal[ifind, *].vdisp_tot = reform([vdans.vdisp_tot], 1, nobj)
   res_gal[ifind, *].vdisp_err = reform([vdans.vdisp_err], 1, nobj)
   res_gal[ifind,*].vdispchi2 = reform([vdans.vdispchi2],1,nobj)
   res_gal[ifind,*].vdispnpix = reform([vdans.vdispnpix],1,nobj)
   res_gal[ifind,*].vdispdof = reform([vdans.vdispdof],1,nobj)
   splog, 'CPU time to fit GALAXY velocity dispersions = ', systime(1)-t0

   res_gal.class = 'GALAXY'
   res_gal.subclass = ' '

   res_all = res_gal ; Append results
   RETURN, res_all
   END
