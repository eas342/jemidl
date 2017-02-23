pro testgalfit
  a=mrdfits('/home/jmeyers314/olivetos/jemidl/Galfit/A0044z.fits',1,/silent)
  restore, '/home/jmeyers314/olivetos/jemidl/Galfit/Apsf.sav'
  psf = zpsf/max(zpsf)
  overpsf = zoverpsf/max(zoverpsf)
  overpsf = undersampleimage(overpsf,3)

  params = { psffactor:1, $
             fit_xmin:1, $
             fit_ymin:1, $
             fit_xmax:31, $
             fit_ymax:30, $
             conv_x:31, $
             conv_y:31, $
             display: 'regular', $
             interactive:0, $
             options:0, $
             platescalex:0.05, $
             platescaley:0.05, $
             zeropoint:24.867 }
  
  constraint1 = obj_new('galfitabsconstraint', 1, 're', 0.18, 300)
  constraint2 = obj_new('galfitrelconstraint', 1, 'mag', -5, 5)
  constraint3 = obj_new('galfitabsconstraint', 1, 'n', 0.5, 8)
  constraint4 = obj_new('galfitabsconstraint', 2, 're', 0.18, 300)
  constraint5 = obj_new('galfitrelconstraint', 2, 'mag', -5, 5)
  constraint6 = obj_new('galfitabsconstraint', 2, 'n', 0.5, 8)
  
  constraints = [constraint1, constraint2, constraint3, constraint4, constraint5, constraint6]
  
  component1 = obj_new('galfitsersic', xpos=15.2681, ypos=14.41, $
                    mag=26.2236, r_e=1.1476, n=1.5, $
                    boa=0.8891, pa=-129.1, out=0)
  component2 = obj_new('galfitsersic', xpos=17.0569, ypos=22.5490, $
                    mag=25.2214, r_e=3.1107, n=1.5, $
                    boa=0.4402, pa=-71.3, out=0)
  component3 = obj_new('galfitsky', sky=-0.0009, dskydx=0, dskydy=0, $
                    out=0, fitsky=0, fitdskydx=0, fitdskydy=0)
  components=[component1, component2, component3]

  outdir='/home/hstdata/users/jmeyers314/svn/scp/members/jmeyers314/IDL/Galfit'

  obj = obj_new('galfit', image=a.image, errim=a.errim, $
                mskim=a.mskim, params=params, constraints=constraints, $
                components=components, psf=psf, prefix='test', outdir = outdir)

  obj->Execute
  results = obj->Results()
  help, results, /st
  obj_destroy, obj

  stop
  ;;now try with an oversampled PSF

  params = { psffactor:3, $
             fit_xmin:1, $
             fit_ymin:1, $
             fit_xmax:31, $
             fit_ymax:30, $
             conv_x:31, $
             conv_y:31, $
             display: 'regular', $
             interactive:0, $
             options:0, $
             platescalex:0.05, $
             platescaley:0.05, $
             zeropoint:24.867 }

  constraint1 = obj_new('galfitabsconstraint', 1, 're', 0.18, 300)
  constraint2 = obj_new('galfitrelconstraint', 1, 'mag', -5, 5)
  constraint3 = obj_new('galfitabsconstraint', 1, 'n', 0.5, 8)
  constraint4 = obj_new('galfitabsconstraint', 2, 're', 0.18, 300)
  constraint5 = obj_new('galfitrelconstraint', 2, 'mag', -5, 5)
  constraint6 = obj_new('galfitabsconstraint', 2, 'n', 0.5, 8)
  
  constraints = [constraint1, constraint2, constraint3, constraint4, constraint5, constraint6]
  
  component1 = obj_new('galfitsersic', xpos=15.2681, ypos=14.41, $
                    mag=26.2236, r_e=1.1476, n=1.5, $
                    boa=0.8891, pa=-129.1, out=0)
  component2 = obj_new('galfitsersic', xpos=17.0569, ypos=22.5490, $
                    mag=25.2214, r_e=3.1107, n=1.5, $
                    boa=0.4402, pa=-71.3, out=0)
  component3 = obj_new('galfitsky', sky=-0.0009, dskydx=0, dskydy=0, $
                    out=0, fitsky=0, fitdskydx=0, fitdskydy=0)
  components=[component1, component2, component3]

  outdir='/home/hstdata/users/jmeyers314/svn/scp/members/jmeyers314/IDL/Galfit'

  obj = obj_new('galfit', image=a.image, errim=a.errim, $
                mskim=a.mskim, params=params, constraints=constraints, $
                components=components, psf=overpsf, prefix='overtest', outdir=outdir)
  obj->Execute
  results = obj->Results()
  help, results, /st


  stop
  ;;now try just modeling a galaxy

  params = { psffactor:3, $
             fit_xmin:1, $
             fit_ymin:1, $
             fit_xmax:31, $
             fit_ymax:30, $
             conv_x:31, $
             conv_y:31, $
             display: 'regular', $
             interactive:0, $
             options:0, $
             platescalex:0.05, $
             platescaley:0.05, $
             zeropoint:24.867 }

  constraint1 = obj_new('galfitabsconstraint', 1, 're', 0.18, 300)
  constraint2 = obj_new('galfitrelconstraint', 1, 'mag', -5, 5)
  constraint3 = obj_new('galfitabsconstraint', 1, 'n', 0.5, 8)
  constraint4 = obj_new('galfitabsconstraint', 2, 're', 0.18, 300)
  constraint5 = obj_new('galfitrelconstraint', 2, 'mag', -5, 5)
  constraint6 = obj_new('galfitabsconstraint', 2, 'n', 0.5, 8)
  
  constraints = [constraint1, constraint2, constraint3, constraint4, constraint5, constraint6]
  
  component1 = obj_new('galfitsersic', xpos=15.2681, ypos=14.41, $
                    mag=26.2236, r_e=1.1476, n=1.5, $
                    boa=0.8891, pa=-129.1, out=0)
  component2 = obj_new('galfitsersic', xpos=17.0569, ypos=22.5490, $
                    mag=25.2214, r_e=3.1107, n=1.5, $
                    boa=0.4402, pa=-71.3, out=0)
  component3 = obj_new('galfitsky', sky=-0.0009, dskydx=0, dskydy=0, $
                    out=0, fitsky=0, fitdskydx=0, fitdskydy=0)
  components=[component1, component2, component3]

  outdir='/home/hstdata/users/jmeyers314/svn/scp/members/jmeyers314/IDL/Galfit'

  obj = obj_new('galfit', image=a.image, errim=a.errim, $
                mskim=a.mskim, params=params, constraints=constraints, $
                components=components, psf=overpsf, prefix='model', outdir=outdir)
  atv, obj->model()

end
