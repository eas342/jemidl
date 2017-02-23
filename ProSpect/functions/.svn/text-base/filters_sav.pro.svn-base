Pro Filters_sav
  ff=obj_new('filterfactory')
  filterdir='/home/scpdata01/FILTERS/'

  BessellUx=ff->createFromFile( filterdir+'STANDARDS/Bessell/BessellUx.txt', /johnson )
  BessellBx=ff->createFromFile( filterdir+'STANDARDS/Bessell/BessellBx.txt', /johnson )
  BessellB=ff->createFromFile(  filterdir+'STANDARDS/Bessell/BessellB.txt', /johnson )
  BessellV=ff->createFromFile(  filterdir+'STANDARDS/Bessell/BessellV.txt', /johnson )
  BessellR=ff->createFromFile(  filterdir+'STANDARDS/Bessell/BessellR.txt', /johnson )
  BessellI=ff->createFromFile(  filterdir+'STANDARDS/Bessell/BessellI.txt', /johnson )


  NIC_F110W=ff->createFromFile( filterdir+'HST/NICMOS/F110W.dat' )
  

  ACS_F775W=ff->createFromFile(  filterdir+'HST/ACS/HST_ACS_F775W_T77.dat' )
  ACS_F850LP=ff->createFromFile( filterdir+'HST/ACS/HST_ACS_F850LP_T77.dat' )

  SDSSu=ff->createFromFile( newfilterdir+'SDSS7/APO_SDSS_U.dat' )
  SDSSg=ff->createFromFile( newfilterdir+'SDSS7/APO_SDSS_G.dat' )
  SDSSr=ff->createFromFile( newfilterdir+'SDSS7/APO_SDSS_R.dat' )
  SDSSi=ff->createFromFile( newfilterdir+'SDSS7/APO_SDSS_I.dat' )
  SDSSz=ff->createFromFile( newfilterdir+'SDSS7/APO_SDSS_Z.dat' )

  filters = { $
            BessellUx  : BessellUx, $
            BessellBx  : BessellBx, $
            BessellB   : BessellB, $
            BessellV   : BessellV, $
            BessellR   : BessellR, $
            BessellI   : BessellI, $
            ACS_F775W  : ACS_F775W, $
            ACS_F850LP : ACS_F850LP, $
            NIC_F110W  : NIC_F110W, $
            SDSSu      : SDSSu, $
            SDSSg      : SDSSg, $
            SDSSr      : SDSSr, $
            SDSSi      : SDSSi, $
            SDSSz      : SDSSz }
  save, filters, filename='filters.sav'
end
