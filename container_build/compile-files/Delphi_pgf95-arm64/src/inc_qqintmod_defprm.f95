      !  ATTENTION!  This file is part of qqinttmod module.
      !==================================================================
      subroutine defprm
      !  logical parameters used in qdiff, mainly from parameter file
      !-------------------------------------------------------
      ! 2011-06-14 Declarations added due to IMPLICIT NONE
      integer :: i
      !!c       data conversion to kT/charge at 25 celsius
      fpi=12.566370614359 !  data fpi /12.566/ data statement discouraged in F95
      
      ideveloper=.false. 
      !   changed for better readability of the code on 2011-04-14
      toplbl="qdiffxas: qdiffxs4 with an improved surfacing routine"
      isch=.false. ; isen=.false. ; iacs=.false. 
      irea=.false. ; iautocon=.true.
      iper=.false.     ! changed on 2011-04-14 to array operation    
      iconc=.false. ; ibios=.false. ;	isite=.false. 
      iatout=.false.;  diff=.false. ;	isph=.false. 
      phiwrt=.false. ; logs=.false. ;     logc=.false. 
      loga=.false. ;	logg=.false.
      !!c b++++++++++++++++++++++++++++++w Oct 2000
      logions=.false.
      !!c flag for energy calculation of contribution from the solvent
      lognl=.false.
      !!c flag for non linear energy calculation
      icreapdb=.false.
      !!c flag for automatic insertion of objects
      imanual=.false.
      !!c flag for manual assignment of relaxation parameterr
      verbose=.true.
      !!c flag for removing part of the standard output
      !!c e++++++++++++++++++++++++++++++++++++++++
      logas=.false.; ipdbrd=.false. ; epswrt=.false. ; iacent=.false.
      ipdbwrt=.false. ; ifrcwrt=.false. ; ifrcrd=.false.; ipoten=.false.
      igraph=.false. ; imem=.false.; ihome=.false. ; ibufz=.false.
      inrgwrt=.false. ; iwgcrg=.false.; iuspec=.false. ; icheb=.false.
      idbwrt=.false. ; ixphird=.false.; iyphird=.false.; izphird=.false.
      ixphiwt=.false. ; iyphiwt=.false. ; izphiwt=.false.; isita=.false.
      isitq=.false. ;	isitp=.false.; 
      isitap=.false.; isitmd=.false. ; isittf=.false. ; isitsf=.false.
      isitf=.false. ;	isitr=.false.;	isitc=.false.;	isitx=.false.
      isiti=.false. ;	iself=.false.;	isitrf=.false. ; isitcf=.false.
      isitt=.false. ;	isolv=.true. ;	isrf=.false.
      
      icon1=10 ;	icon2=1
      
      epsnam="fort.17" ;   phinam="fort.14" ;  srfnam="grasp.srf"
      frcnam="fort.16" ;  mpdbnam="fort.19" ; updbnam="fort.20"
      ufrcnam="fort.21" ; centnam="fort.15" ;  pdbnam="fort.13"
      crgnam="fort.12"  ;  siznam="fort.11" ; phiinam="fort.18"
      frcinam="fort.15" ;  prmnam="fort.10";  scrgnam="scrg.dat"
      nrgnam="energy.dat";gcrgnam="crg.dat";    dbnam="db.dat"
      xphinam="fort.31" ; yphinam="fort.32" ;  zphinam="fort.33"
      xphonam="fort.34" ; yphonam="fort.35" ;  zphonam="fort.36"
      
      epslen=7 ;  philen=7 ;  srflen=9 ;  frclen=7 ; mpdblen=7  
      updblen=7 ; ufrclen=7 ; centlen=7 ;  pdblen=7 ;  crglen=7 
      sizlen=7 ; phiilen=7 ; frcilen=7 ;	scrglen=8 ;  nrglen=10 
      prmlen=7 ; gcrglen=7 ;   dblen=6 ;	xphilen=7 ; yphilen=7
      zphilen=7 ;	xpholen=7 ; ypholen=7 ;	zpholen=7 
      
      phifrm=0 ;	epsfrm=0 ; frcfrm=0 ; mpdbfrm=0 ; updbfrm=0
      ufrcfrm=0;	pdbfrm=0 ; crgfrm=0 ;  sizfrm=0 ;  prmfrm=0
      phiifrm=0; frcifrm=0 ; scrgfrm=0;  nrgfrm=0 ; gcrgfrm=0
      
      radprb(1)=1.4 ; scale=10000. ; exrad=2.0 ;	perfil=10000.
      !!c b++++++++++++++++++
      radpolext=1.0 ; radprb(2)=-1.0
      conc=0.0        ! changed 2011-04-14 to array syntax
      rionst=0.0 ; relpar=1.0 
      !!c ival are the valencies of salts
      ival=1   ! changed 2011-04-14 to array syntax
      ival2=0  ! changed 2011-04-14 to array syntax
      
      res1=0.0 ;  res2=0.0 ;  res5=0.0
      vdrop=coord(0.,0.,0.)
      atompotdist=0.5;   temperature=297.3342119
      
      realsiz=4
      realsiz=realsiz+4
      resnummax=0
      !! e++++++++++++++++++
      repsout=80 ; epsout=80 ; repsin=2 ;	epsin=2
      gten=0.0 ;	uspec=0.9975
      ! 2011-05-09  Changed to coord type variables            
      offset=coord(0.,0.,0.)
      acent=coord(0.,0.,0.)
      bufz%min=int_coord(0,0,0); bufz%max=int_coord(0,0,0)   
      igrid=0 ;	nlit=0 ;  nnit=0 ; ibctyp=2

      !######Lin Li:Gaussian:
          cutoff=1.0
          sigma=1.0
          inhomo=0
          srfcut=20.0
          gaussian=0 
      end subroutine defprm  ! changed syntax 2011-04-14
