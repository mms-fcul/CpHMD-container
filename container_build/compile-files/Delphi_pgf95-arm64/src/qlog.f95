!#####################################################################
! change to F95 by PK from file qlog.h     
! changed 2011-04-14 to module style
! changed 2011-04-15 to F95 syntax
!#####################################################################

      !logical parameters used in qdiff, mainly from parameter file
      module qlog

      logical :: iautocon,iper(6),iconc,ibios,isite,iatout,diff,isph
      logical :: ipdbwrt,ifrcwrt,ifrcrd,ipoten,igraph,imem,ihome,icheb
      logical :: phiwrt,logs,logc,loga,logg,ipdbrd,logas,epswrt,iacent
      logical :: isen,isch,iacs,irea,ibufz,inrgwrt,iexun,iwgcrg,iuspec
      logical :: idbwrt,ixphiwt,iyphiwt,izphiwt,ixphird,iyphird,izphird
      logical :: isita,isitq,isitp,isitf,isitr,isitc,isitx,isiti,iself
      logical :: isitrf,isitcf,isitt,isolv,isrf,ibem
      logical :: lognl,logions,icreapdb,imanual,isitap,isitmd,isittf
      logical ::isitsf,ionlymol,debug,isitpot,verbose,isitdeb
      logical ::ivertnorm,idebwrt,ideveloper !ideveloper is for output format

      !2011-05-10 Transferred from local declaration in main program
      logical :: uniformdiel

      character(80) :: epsnam,phinam,frcnam,mpdbnam,updbnam,ufrcnam
      character(80) :: srfnam,debnam
      character(80) :: centnam,siznam,crgnam,pdbnam,frcinam,phiinam
      character(80) :: prmnam,usernam,uhomenam,scrgnam,nrgnam,gcrgnam
      character(80) :: dbnam
      character(80) :: xphonam,yphonam,zphonam,xphinam,yphinam,zphinam
      character(60) :: toplbl
      real :: scale,prbrad(2),exrad,perfil,repsout,repsin,cutoff,sigma,srfcut,ergsgaussian
      real, parameter :: Pi=3.14159265359

      real :: rionst,conc(2),chi1,chi2,chi3,chi4,chi5,gten,uspec
      real :: radpolext,relpar,epkt,fpi,deblen
      real :: start,timetot

      real :: res1,res2,res5,tol,ergestout
      !real :: res1,res2,res5,tol,ergestout,vdropx,vdropy,vdropz
      real :: atompotdist, temperature
      integer :: ival(2),ival2(2),resnummax,realsiz,npotenziali,nvincfit

      !2011-05-09 acent and offset are now coord type variables
      !           declared in pointers module
      real :: epsout,epsin,radprb(2)
      !real :: offset(3),epsout,epsin,radprb(2),acent(3)
      !integer :: igrid,nlit,nnit,ibctyp,icon1,icon2,bufz(2,3)
      integer :: igrid,nlit,nnit,ibctyp,icon1,icon2
      integer :: epslen,philen,srflen,frclen,mpdblen,updblen,ufrclen
      integer :: phifrm,epsfrm,frcfrm,mpdbfrm,updbfrm,ufrcfrm
      integer :: centlen,sizlen,crglen,pdblen,frcilen,phiilen
      integer :: centfrm,sizfrm,crgfrm,pdbfrm,frcifrm,phiifrm
      integer :: prmfrm,scrgfrm,scrglen,prmlen,userlen,uhomelen
      integer :: nrgfrm,nrglen,gcrglen,gcrgfrm,dblen,xphilen,yphilen
      integer :: zphilen,xphipos,yphipos,zphipos,xphopos,yphopos,zphopos
      integer :: xpholen,ypholen,zpholen,ngp,nhgp,nbgp,debnamlen
      integer :: nobject,nmedia,natom,numbmol,ndistr
      integer :: ibnum,ibmx, ibnumsurf

      !----2011-05-17 Added for module epsmakmod (vwtms subroutine)---
      real :: radpmax,zeta,axdist
      integer :: extot, iall

      !--------------Transferred from acc2.h file---------------------
      real :: grdi,rdmx, cbai, sideinter, sidemin 
      real*4 :: tary(2)
      integer :: lcb, mcb, ncb, lcb1, mcb1, ncb1
      integer, parameter :: idmax=50

      !----2011-05-30 Added for module crgarrmod (crgarr subroutine)--
      integer :: nqass, nqgrd,extracrg
      real :: qmin,qnet,qplus

      !-----2011-05-30 Added for module setcrgmod --------------------
      integer :: nsp,icount2a,icount2b,idirectalg,ngrid,inhomo,gaussian
      integer :: icount1a,icount1b,ibc
      real :: dbval(0:1,0:6,0:1),sfd(5,0:1)

      !-----2011-05-30 Added for module setbcmod ---------------------
      real :: spec

      !-----2011-05-30 Added for module setrc and ass-----------------
      integer :: norecord

      !-----2011-07-18 Added for timing-------------------------------
      character(10) :: day,time
      integer :: values(8)
      
!*********************************************************************
!Common blocks become unnecessary in module architecture 
!removed on 2011-04-14
!	common
!     &  /log1/iautocon,iper,iconc,ibios,isite,iatout,diff,isph,
!     &  ipdbwrt,ifrcwrt,ifrcrd,ipoten,igraph,imem,ihome,icheb,
!     &  phiwrt,logs,logc,loga,logg,ipdbrd,logas,epswrt,iacent,
!     &  isen,isch,iacs,irea,ibufz,inrgwrt,iexun,iwgcrg,iuspec,
!     &  idbwrt,ixphiwt,iyphiwt,izphiwt,ixphird,iyphird,izphird,
!     &  isita,isitq,isitp,isitf,isitr,isitc,isitx,isiti,iself,
!     &  isitrf,isitcf,isitt,isolv,isrf,ibem,ionlymol,debug,
!     &  lognl,logions,icreapdb,imanual,isitap,isittf,isitmd,isitsf,
!     &  isitpot,verbose,isitdeb
!
!	common
!     &  /val1/scale,prbrad,exrad,perfil,rionst,repsout,repsin,
!     &  gten,offset,epsout,epsin,radprb,acent,epkt,fpi,deblen,uspec,
!     &  conc,chi1,chi2,chi3,chi4,chi5,radpolext,relpar,res1,res2,tol,
!     &  ergestout,res5,vdropx,vdropy,vdropz,atompotdist,temperature
!
!	common
!     &  /ival1/igrid,nnit,nlit,ibctyp,epslen,philen,srflen,frclen,
!     &  mpdblen,updblen,ufrclen,centlen,sizlen,crglen,pdblen,
!     &  frcilen,phiilen,sizfrm,crgfrm,pdbfrm,frcifrm,phiifrm,
!     &  phifrm,epsfrm,frcfrm,mpdbfrm,updbfrm,ufrcfrm,prmfrm,
!     &  prmlen,icon1,icon2,userlen,uhomelen,scrglen,scrgfrm,bufz,
!     &  nrgfrm,nrglen,gcrglen,gcrgfrm,dblen,xphilen,yphilen,zphilen,
!     &  xphipos,yphipos,zphipos,xphopos,yphopos,zphopos,
!     &  xpholen,ypholen,zpholen,ngp,nhgp,nbgp,ival,ival2,resnummax,
!     &  realsiz,npotenziali,nvincfit
!
!     &  /icar1/epsnam,phinam,frcnam,mpdbnam,updbnam,ufrcnam,centnam,
!     &  uhomenam,usernam,scrgnam,nrgnam,gcrgnam,dbnam,srfnam,
!     &  xphonam,yphonam,zphonam,xphinam,yphinam,zphinam,
!     &  prmnam,pdbnam,siznam,crgnam,phiinam,frcinam,toplbl
!
      end module qlog
