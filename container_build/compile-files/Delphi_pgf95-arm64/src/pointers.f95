!#####################################################################
!Created by PK 2011-04-21 from files qdiffpar5.h and pointer.h
!allocatable arrays for qdiff4.f and subroutines
!
!2011-04-21  some parameters are now obsolete due to allocatable arrays
!parameter (nclist = 1500)
!ncmax= maximum number of entries in charge file 
!parameter (ncmax = 1500)
!ncrgmx = maximum number of charged atoms
!parameter (ncrgmx = 50000)
!ngrid= maximum grid size 
!parameter (ngrid = 1025)   ***********IRIX64*****
!parameter (ngrid = 513)
!ngp = ngrid**3 = number of grid points 
!parameter (ngp = 274626)
!nhgp= half ngp
!parameter (nhgp = 137313)
!natmax = maximum number of atoms 
!parameter (natmax =100000)
!nmediamax = maximum number of different media
!parameter (nmediamax = 1000)
!nobjectmax = maximum number of objects
!parameter (nobjectmax = 1000)
!ndistrmax = maximum number of charge distributions    
!parameter (ndistrmax = 1000)
!ibcmax = maximum no. grid points charged and at boundary
!parameter (ibcmax = 30000)
!ngcrg = maximum number of grid points that may be assigned charge 
!parameter (ngcrg =110000)
!nbgp = number of points at one box boundary = ngrid**2
!parameter (nbgp = 4226)
!nsp = maximum number of dielectric boundary points, appox.= 5*nbgp 
!parameter (nsp = 40000)
!---------------------------------------------------
!2011-04-24  Derived variable type for records in radii/charge table
!            and corresponding allocatable array
!nrmax=  number of entries in radius file (found in rdh module) 
!parameter (nrmax = 1500)
!parameter (nrlist = 1500)   !  not necessary
!-------- atom names for radii table
!dimension atnam(nrmax)
!-------- residue names for radii table
!dimension rnam(nrmax)		
!-------- radius table 
!dimension radt(nrmax)		
!-------- chain names for radii table
!dimension rchn(nrmax)		
!-------- residue number for radii table
!dimension rrnum(ncmax)		
!-------- links for radius hash table
!dimension irlink(nrlist)	
!--------- radii id numbers in hash table
!dimension irnumb(nrlist)	
!-------- links for radius hash table
!dimension irlink(nrlist)	
!-------- radii id numbers in hash table
!dimension irnumb(nrlist)	
!-------- atom names   for charge table
!dimension catnam(ncmax)		
!-------- chain names for charge table
!dimension cchn(ncmax)		
!-------- residue names for charge table
!dimension crnam(ncmax)		
!-------- residue number for charge table
!dimension crnum(ncmax)		
!-------- charge table 
!dimension chrgvt(ncmax)		
!-------- links for charge hash table
!dimension iclink(nclist)	
!-------- charge entry id numbers in hash table
!dimension icnumb(nclist)	
!nrmax=nrdrec=  number of entries in radius file (found in rdh module) 
!ncmax=nchrec=  number of entries in charge file (found in rdh module) 
!irtot=  number of entries in hash radius file (found in rdh module) 
!ictot=  number of entries in hash charge file (found in rdh module) 
!character*1 chn,cchn,rchn
!character*3 crnam,rnam,res
!character*4 rnum,crnum,rrnum
!character*6 atm,catnam,atm
!real, allocatable :: chgpos(:,:) !dimension chgpos(3,ncrgmx),schrg(1)
!pointer (i_chgpos,chgpos)
!pointer (i_phimap,phimap),(i_phimap1,phimap1)
!pointer (i_phimap2,phimap2),(i_phimap3,phimap3)
!pointer (i_db,db),(i_idpos,idpos),(i_sf1,sf1),(i_sf2,sf2)
!pointer (i_qmap1,qmap1),(i_qmap2,qmap2)
!pointer (i_debmap1,debmap1),(i_debmap2,debmap2)
!pointer (i_bndx1,bndx1),(i_bndx2,bndx2),(i_bndx3,bndx3)
!pointer (i_bndx4,bndx4)
!pointer (i_ibndx,ibndx),(i_ibndy,ibndy),(i_ibndz,ibndz)
!pointer (i_neps,neps),(i_keps,keps)
!pointer (i_iepsmp,iepsmp),(i_idebmap,idebmap)
!integer limeps(2,3)
!pointer (i_limgunit,limgunit)
!real oldmid(3),oldmid1(3),pmid(3)
!pointer (i_limobject,limobject)
!pointer (i_cbal,cbal),(i_icbn,icbn) 
!pointer (i_iatmobj,iatmobj) 
!pointer (i_expos,expos),(i_pls,pls),(i_ast,ast),(i_ast2,ast2)
!pointer (i_coi,coi)
!integer, allocatable :: cbn1(:,:,:), cbn2(:,:,:) ! in vwtms2.f 
!pointer (i_cbn1,cbn1),(i_cbn2,cbn2)
!pointer (i_r0,r0),(i_r02,r02),(i_rs2,rs2)
!pointer (i_ibnd,ibnd),(i_ibgrd,ibgrd),(i_bndeps,bndeps) 
!pointer (i_atndx,atndx),(i_scspos,scspos),(i_atsurf,atsurf)
!pointer (i_scsnor,scsnor),(i_vnorm,vnorm)
!integer, allocatable :: iab1(:,:,:),iab2(:,:,:),icume(:) !in vwtms2.f
!pointer (i_iab1,iab1),(i_iab2,iab2),(i_icume,icume)
!pointer (i_ioff,ioff)
!pointer (i_atmeps,atmeps)
!pointer (i_nqgrdtonqass,nqgrdtonqass)
!pointer (i_cgrid,cgrid),(i_spdiv,spdiv),(i_sen,sen)
!pointer (i_schrg,schrg),(i_crgatn,crgatn)
!pointer (i_spot,spot),(i_sqs,sqs),(i_iepsv,iepsv)
!pointer (i_chrgv2,chrgv2),(i_cqs,cqs)
!pointer (i_sitephi,sitephi) !in react2.f
!pointer (i_rfield,rfield)
!pointer (i_atmforce,atmforce)
!pointer (i_iexpos,iexpos)
!real, allocatable :: atpos(:),rad3(:), chrgv4(:) !in getatm2.f
!character(15), allocatable :: atinf(:) !in getatm2.f
!pointer (i_iatmmed,iatmmed) !in getatm2.f 
!pointer (i_atpos,atpos),(i_xn2,xn2),(i_rad3,rad3) ! xn2 ???
!pointer (i_chrgv4,chrgv4),(i_atinf,atinf)                       
!pointer (i_medeps,medeps)
!pointer (i_dataobject,dataobject)
!pointer (i_datadistr,datadistr)
!pointer (i_internal,internal)
!pointer (i_tmpmap,tmpmap)
!pointer (i_gchrgtmp,gchrgtmp)
!pointer (i_cgbp,cgbp)
!integer, allocatable :: gchrgp(:,:) !in setcrg.f and setfcrg.f
!pointer (i_qval,qval),(i_gchrg,gchrg),(i_gchrgp,gchrgp)
!pointer (i_iqpos,iqpos)
!pointer (i_gchrgd,gchrgd),(i_gchrg2,gchrg2),(i_gval,gval)
!real, allocatable :: atmcrg(:,:)  Now deruved type ! in crdarr.f
!pointer (i_atmcrg,atmcrg) 
!real, allocatable :: vert(:,:), vnorm(:,:), vnorm2(:,:) !in msrf.f
!pointer (i_vert,vert),(i_vindx,vindx),(i_vtemp,vtemp)
!pointer (i_nsel,nsel),(i_vtlen,vtlen),(i_vtlst,vtlst)
!pointer (i_tmlst,tmlst),(i_vtpnt,vtpnt)
!pointer (i_sout,sout)   
!pointer (i_phimap4,phimap4)
!pointer (i_polariz,polariz) !used exclusively in rforcenew subroutine
!pointer (i_epsmap,epsmap)
!pointer (i_debmap,debmap),(i_qmap,qmap)

!other arrays
!atmcrg contians the postions of all chrages in grid units, and the 
!charge itself in the fourth field. chgpos contians the same but in 
!angstroms schrg contains the induced surface charges in electrons
!dimension atmcrg(4,ncrgmx),chgpos(3,ncrgmx),schrg(nsp)
!dimension phimap(ngrid,ngrid,ngrid),phimap1(nhgp),phimap2(nhgp)
!dimension phimap3(ngp)
!integer iepsmp(ngrid,ngrid,ngrid,3),idebmap(ngrid,ngrid,ngrid)
!integer atsurf(nsp)
!integer atsurf(1)
!dimension schrg(1)
!dimension cgbp(5,20000),atmcrg(4,ncrgmx)
!logical logtab(100)
!#####################################################################

      module pointers

      integer, parameter :: ncrgmx = 50000
      integer :: nrmax, nrdrec ,ncmax, nchrec !makes those variables 
                                              !public
      integer :: irtot, ictot !makes those variables public 

      type :: coord
         real :: x,y,z
      end type coord

      type :: object_min_max
         type(coord) :: min, max
      end type object_min_max

      type :: int_coord
         integer :: i,j,k
      end type int_coord

      type :: int_extrema
         type(int_coord) :: min, max
      end type int_extrema

      !---------------------------------------------------
      !2011-04-24 Instead for 10=5x2 separate arrays 
      !           atnam, rnam, radt, rchn, rrnumb
      !           catnam, crnam, crnum, cchn, chrgvt                 
      !           for records in radii/charge files
      !---------------------------------------------------
      type parameter_file_record
         character(6) :: atnam !atom name in radii/charge table
         character(3) :: rnam  !residue name in radii/charge table   
         character    :: chn   !chain name in radii/chatge table
         character(4) :: rnum  !residue number in radii/charge table
         real         :: value !value of atom charge/radius 
      end type parameter_file_record

      type(parameter_file_record), allocatable :: radii(:), charge(:)

      !---------------------------------------------------
      !2011-04-24 Instead for 4=2x2 separate hash arrays 
      !           irlink,irnumb, iclink, icnumb
      !           for records in radii/charge files
      !---------------------------------------------------
      type hash_table       
         integer :: ilink !links for hash table
         integer :: inumb !id numbers in hash table
      end type hash_table

      type(hash_table), allocatable :: rhash(:), chash(:)

      !---------------------------------------------------
      type delphi_pdb_file_record
         real         :: rad3   !atom radius instead of rad3 array  
         real         :: chrgv4 !atom charge instead of chrgv4
         type(coord)  :: xyz    !atom coordinates
         !real        :: x      !  
         !real        :: y      !instead of atpos array  
         !real        :: z      !  
         character(15) :: atinf !instead of atinf array
      end type delphi_pdb_file_record

      type(delphi_pdb_file_record), allocatable :: delphipdb(:)

      !2011-05-30 New generic variable type for grid positions and value
      !           at these positions
      type grid_value
           type(coord) :: xyz
           real :: value
      end type grid_value

      type int_grid_value
         type(int_coord):: ijk
         real :: value
      end type int_grid_value

      type double_grid_value
         type(coord) xyz
         real :: value1
         real :: value2
      end type double_grid_value

      !temporarly character varibales used in the assignment of radii
      character(6) :: atm      !  atom name 
      character(3) :: res      !  residue name 
      character ::  chn        !  chain name 
      character(4) ::  rnum    !  residue number 

      !******** ALLOCATABLE ARRAYS USED IN MAIN PROGRAM (qdiff4v)
      !------------------------------------------------------------------
      !2011-04-21 Allocatable arrays from pointer.h file used in main 
      !           program
      !------------------------------------------------------------------
      real, allocatable :: phimap(:,:,:), phimap1(:)
      real, allocatable :: phimap2(:), phimap3(:) !also in mkdbsf.f and 
                                                  !wrtphi.f
      real, allocatable :: db(:,:), sf1(:), sf2(:)
      integer, allocatable :: idpos(:)
      real, allocatable :: qmap1(:), qmap2(:), debmap1(:), debmap2(:)
      real, allocatable :: bndx1(:), bndx2(:), bndx3(:), bndx4(:)

      type(int_coord), allocatable :: neps(:,:,:)
      integer, allocatable :: ibndx(:),ibndy(:),ibndz(:),keps(:,:,:)

      type(int_extrema) :: limeps !in epsmak.f
      type(object_min_max), allocatable :: limobject(:) !in main,extrm 
                                                        !extrmobject, 
      type(object_min_max), allocatable :: limgunit(:)  !in epsmak.f
      type(int_coord), allocatable :: iepsmp(:,:,:)     !in epsmak.f
      type(coord), allocatable :: gepsmp(:,:,:),gepsmp2(:,:,:) ! in epsmak.f for gaussian eps

      logical, allocatable :: idebmap(:,:,:)            !in epsmak.f

      !2011-05-09 Using coord type variables
      type(coord) :: oldmid,oldmid1,pmid,cmin,cmax,cran,xo,offset
      type(coord) :: acent,vdrop,xl,xr

      !------------------------------------------------------------------
      !2011-04-23  Allocatable arrays from pointer.h file used in cube.f
      !------------------------------------------------------------------
      integer, allocatable :: cbal(:), icbn(:,:) !in vwtms2.f, cube.f 
      integer, allocatable :: iatmobj(:)         !in cube.f

      !------------------------------------------------------------------
      !2011-04-23  Allocatable arrays from pointer.h file used in sas.f
      !------------------------------------------------------------------
      type(int_coord), allocatable :: pls(:) !in sas.f
      type(coord), allocatable :: expos(:) !in sas.f and vwtms2.f

      real, allocatable :: coi(:,:) ! in sas.f

      !------------------------------------------------------------------
      !2011-04-23 Allocatable arrays from pointer.h file used in 
      !           vwtms2.f
      !------------------------------------------------------------------
      integer, allocatable :: cbn1(:), cbn2(:) !in vwtms2.f 
      integer, allocatable :: ast(:) !in vwtms2.f and sas.f
      real, allocatable :: r0(:),r02(:),rs2(:) !in vwtms2.f
      type(int_coord), allocatable :: ibnd(:),ibgrd(:) !in vwtms2.f
      integer, allocatable :: bndeps(:,:,:,:) !in vwtms2.f
      type(coord), allocatable :: scspos(:) !in vwtms2.f
      integer, allocatable :: atsurf(:),atndx(:) !in vwtms2.f
      type(coord), allocatable :: scsnor(:) !in vwtms2.f
      integer, allocatable :: iab1(:),iab2(:),icume(:) !in vw
      type(int_coord), allocatable :: egrid(:,:,:)

      !2011-05-18 Transferred from acc2.h file
      type(coord) :: mnxyz, xyzo, mxxyz
      type(int_coord) :: lmncb1, lmncb

      !2011-04-23 Allocatable arrays from pointer.h file used in 
      !           setout.f
      type(int_coord), allocatable :: ioff(:) !in setout.f

      !2011-05-30  Allocatable arrays allocated and variables used in 
      !crgarr.f
      type(grid_value), allocatable :: atmcrg(:),chrgv2(:)
      type(coord), allocatable :: chgpos(:)
      integer, allocatable :: crgatn(:), nqgrdtonqass(:)
      real, allocatable :: atmeps(:)
      type(coord) :: cqplus,cqmin,cmid

      !2011-04-23 Allocatable arrays from pointer.h file used in 
      !           react2.f
      real, allocatable :: cgrid(:,:,:), spdiv(:), sen(:) !in react2.f 
      real, allocatable :: schrg(:) !in react2.f
      real, allocatable :: spot(:),sqs(:) !in react2.f
      real, allocatable :: cqs(:) !in react2.f
      real, allocatable :: sitephi(:,:)

      !2011-04-23 Allocatable arrays from pointer.h file used in 
      !           mkdbsf.f
      integer, allocatable :: iepsv(:) !in mkdbsf.f

      !2011-04-23 Non-standard variables and allocatable arrays
      !           used in encalc subroutine
      type(int_extrema) :: bufz

      !2011-04-23 Allocatable arrays from pointer.h file used in 
      !           wrtsit4.f
      type(coord), allocatable :: rfield(:), atmforce(:) !in wrtsit4.f

      !2011-04-23 Allocatable arrays from pointer.h file used in 
      !           epsmak.f

      !2011-04-23 Allocatable arrays from pointer.h file used in 
      !           indver.f
      type(int_coord), allocatable :: iexpos(:) ! in indver.f

      !2011-04-23 Allocatable arrays from pointer.h file used in 
      !           getatm2.f

      !2011-05-01 Derived type to replace 4 arrays used in getatm2.f
      type(coord), allocatable :: xn1(:), xn2(:) !in main and vwtms2.f
      real, allocatable :: medeps(:) !in getatm2.f
      character(96), allocatable :: dataobject(:,:) !in getatm2.f
      character(80), allocatable :: datadistr(:) !in getatm2.f
      integer, allocatable :: iatmmed(:), tmpiatmmed(:) !in getatm2.f 

      !2011-04-23 Allocatable arrays from pointer.h file used in scale.f
      logical, allocatable :: internal(:) !in scale.f

      !2011-04-23 Allocatable arrays from pointer.h file used in 
      !           setcrg.f and setfcrg.f
      real, allocatable :: tmpmap(:,:,:) !in setcrg.f and setfcrg.f
      real, allocatable :: gchrgtmp(:) !in setcrg.f and setfcrg.f 
      type(double_grid_value), allocatable :: cgbp(:) !in setcrg.f and 
                                                      !setfcrg.f
      real, allocatable :: qval(:),gchrg(:) !in setcrg.f and setfcrg.f
      type(int_coord), allocatable :: gchrgp(:) !in setcrg.f and 
                                                !setfcrg.f
      integer, allocatable :: iqpos(:) !in setcrg.f and setfcrg.f

      real, allocatable :: gchrgd(:), gval(:) !in setcrg.f and setfcrg.f
      integer, allocatable :: gchrg2(:) !in setcrg.f and setfcrg.f

      !2011-04-23 Allocatable arrays from pointer.h file used in 
      !         crdarr.f

      !2011-04-23 Allocatable arrays from pointer.h file used in msrf.f
      type(int_coord), allocatable :: vindx(:)
      integer, allocatable :: vtlen(:), vtlst(:),tmlst(:,:) !in msrf.f
      integer, allocatable :: vtpnt(:) !in msrf.f
      type(coord), allocatable :: vert(:), vnorm(:), vnorm2(:) !in 
                                                               !msrf.f

      !2011-04-23 Allocatable arrays from pointer.h file used in clb.f 
      !and nlener.f
      type(grid_value), allocatable :: sout(:) !in clb.f and nlener.f

      !2011-04-23 Allocatable arrays from pointer.h file used in clb.f 
      !           and wrtphi.f
      real, allocatable :: phimap4(:) !in wrtphi.f

      !2011-04-23 Allocatable arrays from pointer.h file with 
      !           unidentified location

!---------------------------------------------------
!   2011-04-21 common block removed due to module architecture
!	common
!     &	/link/  irlink,irnumb,iclink,icnumb,irtot,ictot
!     &	/name/  atnam,rnam,catnam,cchn,rchn,crnam,crnum,rrnum
!     &	/value/ radt,chrgvt
!c    &	/maps/  phimap,phimap1,phimap2,phimap3
!c    &  /array/ atmcrg,chgpos,schrg
!c    &  /imaps/ iepsmp,iepsmp2,idebmap,atsurf
!c    &	/scale/ oldmid,scale1,oldmid1,ibc,cgbp,gval,rmmin,rmmax
!     &	/scale/ oldmid,scale1,oldmid1,rmmin,rmmax,pmid
!     &	/scaleint/ ibc
!c    &  /array/atmcrg,chgpos
!c    &	/log/   logtab
 
      end module pointers
