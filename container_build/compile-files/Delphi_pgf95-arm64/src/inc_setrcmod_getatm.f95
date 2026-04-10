      !  ATTENTION!  This file is part of setrcmod module.
      !==================================================================
      ! 2011-04-29  Other arguments to subroutine now passed via pointers, qlog 
      !*****************************************************************************
      subroutine getatm(fname,nam1)
      !*****************************************************************************
      ! read a pdb file and extract out the coordinates, the atom info
      ! and optionally the extra two fields
      
      ! 2011-04-29  Parameters removed due to allocatable arrays
      !2011-04-29     Those variables decleared in the beginning of the module
      ! b+++++++++++++++++++
      ! 2011-05-01 allocatable arrays just for this sub (obsolete?)
      ! iatmmed :vector containing internal media-number per atom and object
      ! nobject= number of objects, any molecule is counted as an object
      ! medeps :vector containing correspondence media<->epsilon/epkt
      ! dataobject: vector containing string with object data, and pre-elab data
      ! ionlymol: flag : .false. => there are objects other than molecules
      ! objecttype : 1-sphere, 2-cylinder --
      ! datadistr : vector containing string with charge distribution data
      ! e+++++++++++++++++++
      
      character(80) ::  fname,fname2,asci
      integer :: atot,objecttype
      logical :: flag
      character(15) :: atminf
      character(10) ::  cnum='0123456789'
      character(6) :: headfirst=" "
      !-------------------------------------------------------------------
      ! 2011-06-14 Declarations added due to IMPLICIT NONE
      integer :: nmediamax,nobjectmax,ndistrmax,imedianumb
      integer :: ios,ios1,ios2,i,j,io,ii,nam1
      real :: rad,chrgv,radx,crgx,repsintmp
      !-------------------------------------------------------------------
      fname(nam1+1:len(fname))=' '
      atot=0
      
      ! assume a formatted file..
      ! open and determine format of file, i.e. formatted or unformatted..
      
      call form(fname,nam1,ifrm)
      
      ! ifrm=.true. then this is a formatted file
      ! 2011-05-01  New piece of code since the old one is too messy
      !----------------------------------------------------------------------------------------
      !   2011-05-01   Setting initial/default values of parameters
      !----------------------------------------------------------------------------------------
      nmediamax=1; nmedia=1
      nobjectmax=1; nobject=1
      objecttype=0
      imedianumb=1
      ndistrmax=0
      Natom=0
      atot=0
      repsintmp=repsin
      iatinf=.false. ; iatrad=.false.;  iatcrg=.false.
      ifrm=.true.; ionlymol=.true.
      i=0
      idfrm=0      ! 2011-05-01; in case FORMAT = idfm is not present
      !===========================================================================
      !   2011-05-01   Formatted (true) or unformatted file (false)
      !===========================================================================
      if(ifrm)  then          ! 2011-05-01; Reading formatted file
         write(6,*) 'opening formatted file:',trim(fname)
         open(9,file=fname(:nam1),form='formatted')
         !------------------------------------------------------------------------------
         ! 2011-05-01 If DELPHI keyword is present in first line, read FORMAT label idfm
         !------------------------------------------------------------------------------
         read(9,'(a80)',iostat=ios) line
         if (ios.ne.0) then
            write(6,*) 'An error occurred in reading first &
            &line of this formatted file'
            return
         end if
         if(index(line,"DELPHI").ne.0) then
            if(index(line,"PDB").eq.0) then
               write(6,*) "this is not a delphi pdb file! &
               &Check the file type!"
               return                 ! 2011-04-30 instead of  goto 999
            end if
            write(6,*) 'Reading a Delphi-style pdb file'
            read(9,'(a80)',iostat=ios)line
            if (ios.ne.0) then
               write(6,*) &
               &'An error occurred in reading this formatted file'
               return
            end if
            !  look for FORMAT NUMBERa=
            j=index(line,"=")
            read(line(j+1:80),'(i5)',iostat=io) idfrm
            if(io.ne.0) then
               write(6,*)'Error in reading Delphi Format number,'
               idfrm=0
               write(6,*)'assuming Delphi Format number =',idfrm
            end if
            write(6,*) 'Read Delphi Format Number = ',idfrm
         else                     !  2011-05-01 No DELPHI in first line
            rewind 9              !          start reading from beginning
            write(6,*)&
            &'No DELPHI keyword, assuming Delphi Format number =',idfrm
         end if
         
         !_________________Lin Li: for modpdb4
         if(pdbfrm.ge.40) idfrm=pdbfrm
         if(pdbfrm.eq.3) idfrm=13 !just for pqr file,can not use #3 becuase it has been used.
         
         !------------------------------------------------------------------------------
         ! 2011-05-01 Setting logical parameters according to idfrm
         !------------------------------------------------------------------------------
         if (idfrm.eq.0) then
            iatinf=.true.; iatrad=.false.; iatcrg=.false.
         else
            iatinf=.true.; iatrad=.true.; iatcrg=.true.
         end if
         !------------------------------------------------------------------------------
         ! 2011-05-01 First reading of file, looking for different keywords in pdb file 
         !------------------------------------------------------------------------------
         do
            read(9,'(a80)',iostat=ios) line
            if(ios.ne.0) exit     ! 2011-05-01 Instead of err= end= labels
            head=line(1:6)
            call up(head,6)       ! 2011-05-01 Converts to capital letters
            select case(head)
            case('MEDIA ')        ! 2011-05-01 Reads media numbers
               read(line(8:10),'(i3)') nmedia
               if(nmedia.gt.nmediamax) nmediamax=nmedia
            case('OBJECT')        ! 2011-05-01 Reads object numbers
               read(line(8:10),'(i3)') nobject
               if(nobject.gt.nobjectmax) nobjectmax=nobject
            case('CRGDST')        ! 2011-05-01 Reads charge 
               read(line(8:10),'(i3)') ndistr
               if(ndistr.gt.ndistrmax) ndistrmax=ndistr
            case('ATOM  ','HETATM')          ! Counts number of atoms
               Natom=Natom+1
            end select
         end do                     ! 2011-05-01 End of first file reading
         rewind 9
         !-------------------------------------------------------------------------
         ! 2011-05-01 allocation of necessary arrays of above determined sizes
         !-------------------------------------------------------------------------
         allocate(delphipdb(Natom))
         allocate(medeps(0:nmediamax))
         allocate(dataobject(nobjectmax,2))
         allocate(datadistr(Ndistrmax))
         allocate(tmpiatmmed(nobjectmax))
         allocate(iatmmed(Natom+Nobjectmax))
         !-------------------------------------------------------------------------
         ! 2011-05-01 Second file reading , putting values in proper arrays,
         !-------------------------------------------------------------------------
         line=' '
         atot=0; i=0
         read(9,'(a80)',iostat=ios) line
         headfirst=line(1:6)
         DREAD: do
            i=i+1
            head=line(1:6)
            call up(head,6)         ! 2011-05-01 Converts to capital letters
            select case(head)
            case('MEDIA ')
               read(line(8:10),'(i3)') nmedia
            case('OBJECT')
               read(line(8:10),'(i3)',iostat=ios)nobject   
               read(line(12:14),'(i3)',iostat=ios)objecttype
               read(line(16:18),'(i3)',iostat=ios)imedianumb
               read(line(20:27),'(f8.3)',iostat=ios)repsintmp
               if(ios.ne.0) then
                  write(6,*)' Error reading OBJECT line ',i ; cycle DREAD
               end if
               medeps(imedianumb)=repsintmp/epkt
               read(9,'(a80)',iostat=ios) line
               if (objecttype.ne.0) then
                  ionlymol= .false. 
                  dataobject(nobject,1)=line
               else
                  dataobject(nobject,1)="is a molecule    0"
                  tmpiatmmed(nobject)=imedianumb
               end if
               cycle DREAD
            case('CRGDST')
               ionlymol= .false.
               read(line(8:10),'(i3)',iostat=ios)ndistr    
               if(ios.ne.0) then
                  write(6,*)' Error reading CRGDST line ',i ; cycle DREAD
               end if
               datadistr(ndistr)=line
            case('ATOM  ','HETATM')          ! Atoms
               atot=atot+1
               delphipdb(atot)%atinf=line(12:26)
               read(line(31:54),'(3f8.3)',iostat=ios) delphipdb(atot)%xyz
               iatmmed(atot)=imedianumb
               if(ios.ne.0) then
                  write(6,*)' Error reading coordinates in line ',i ; cycle DREAD
               end if
               !------------------------------------------------------------------------------
               ! 2011-05-01 Different values of DELPHI FORMAT label idfrm
               !------------------------------------------------------------------------------
               select case(idfrm)                      
                  !------------------------------------------------------------------------------
                  ! 2011-05-01 DELPHI FORMAT label idfm=0, standard pdb file, no charges and radii
                  !------------------------------------------------------------------------------
               case(0)
                  radx=0. ; crgx=0. 
                  !----------------------------------------------------------------------------------------
                  ! 2011-05-01 DELPHI FORMAT label idfm=1,2,3,4, charges and radii presented in
                  !----------------------------------------------------------------------------------------
               case(1)                        ! 2011-05-01 idfm=1
                  read(line(55:60),'(f6.2)',iostat=ios) radx
                  read(line(61:67),'(f7.3)',iostat=ios) crgx
               case(2,3)                        ! 2011-05-01 idfm=2
                  read(line(55:80),'(2f12.7)',iostat=ios) radx,crgx
               case(4)                        ! 2011-05-01 idfm=4
                  read(line(63:70),'(f8.3)',iostat=ios) radx
                  read(line(55:62),'(f8.4)',iostat=ios) crgx
               case(13)                        ! 2011-05-01 idfm=2
                  read(line(55:68),'(2f7.3)',iostat=ios) crgx,radx
               case(42)                        ! 2011-05-01 idfm=2
                  read(line(55:69),'(f7.4,f8.4)',iostat=ios) radx,crgx
               case(43)                        ! 2011-05-01 idfm=2
                  read(line(55:69),'(f8.4,f7.4)',iostat=ios) crgx,radx
                  
                  
               end select
               delphipdb(atot)%rad3=radx
               delphipdb(atot)%chrgv4=crgx
            end select
            read(9,'(a80)',iostat=ios) line
            if(ios.ne.0) exit DREAD
         end do DREAD          ! 2011-05-01 End of second file reading
         close(9)
         write(6,*) "number of atoms read in = ",atot," formatted"
         
         !==============================================================================
         !   2011-05-01  Reading unformatted file (ifrm=.false.)
         !==============================================================================
      else                       ! 2011-05-01 Reading unformatted file
         !-----------------------------------------------------------------
         !   211-05-01   Reading first two lines
         !-----------------------------------------------------------------
         idfrm=0; atot=0; flag=.true.
         open(9,file=trim(fname),form='unformatted')
         read(9,iostat=ios1) line
         if(ios1.ne.0) then
            write(6,*)'Error in reading first line of this unformatted file'
            flag=.false.
         end if
         if(index(line,"DELPHI").ne.0)then
            if(index(line,"PDB").eq.0) then
               write(6,*) "this is not a delphi pdb file!"
               write(6,*) "check the file type"
               return
            end if
         end if
         write(6,*) "Reading a Delphi-style pdb file"
         read(9,iostat=ios2)idfrm
         if(ios2.ne.0) then
            write(6,*)'Error in reading second line of this unformatted file'
            flag=.false.
         end if
         !-----------------------------------------------------------------
         !   211-05-01   Reading rest of the file
         !-----------------------------------------------------------------
         select case(idfrm)
         case(0)
            rewind(9)
            iatinf=.false.; iatrad=.true.; iatcrg=.true.; Natom=0
            !-----------------------------------------------------------------
            ! 2011-05-01   First reading of file to determine number of atoms
            !-----------------------------------------------------------------
            if (flag) then   !  skipping first two lines in count
               read(9,iostat=ios1) line
               read(9,iostat=ios2)idfrm
            end if
            do               ! 2011-05-01  Reading rest of file
               read(9,iostat=io) xo,radx,crgx
               if(io.ne.0) exit
               Natom=Natom+1
            end do
            rewind(9)
            allocate(delphipdb(Natom)) 
            if (flag) then   !  skipping first two lines in count
               read(9,iostat=ios1) line
               read(9,iostat=ios2)idfrm
            end if
            !-----------------------------------------------------------------
            ! 2011-05-01   Second reading of file to read actual data
            !-----------------------------------------------------------------
            do i=1,Natom      
               read(9,iostat=io) delphipdb(i)%xyz,radx,crgx
               if(io.ne.0) exit
               delphipdb(i)%atinf=' '
               delphipdb(i)%chrgv4=crgx
               delphipdb(i)% rad3=radx
            end do
            !-----------------------------------------------------------------
            ! 2011-05-01   Another form of unformatted file with 
            !-----------------------------------------------------------------
         case(1)          
            read(9,iostat=ios) Natom
            iatinf=.true.
            iatrad=.true.
            iatcrg=.true.
            allocate(delphipdb(Natom))
            do i=1,Natom
               read(9,iostat=ios) atminf,delphipdb(i)%xyz,radx,crgx
               delphipdb(i)%atinf=atminf
               delphipdb(i)%chrgv4=crgx
               delphipdb(i)% rad3=radx
            end do
         end select
         write(6,*) "number of atoms read in = ",Natom,' unformatted'
      end if
      !-----------------------------------------------------------------
      ! 2011-05-01   Some post-processing
      !-----------------------------------------------------------------
      if((ifrm.and.headfirst.ne.'MEDIA').or.(.not.ifrm.and.idfrm.eq.0)) then
         if(.not.allocated(medeps))allocate(medeps(0:nmediamax))
         if(.not.allocated(dataobject))allocate(dataobject(nobjectmax,2))
         if(.not.allocated(datadistr))allocate(datadistr(Ndistrmax))
         if(.not.allocated(tmpiatmmed))allocate(tmpiatmmed(nobjectmax))
         if(.not.allocated(iatmmed))allocate(iatmmed(Natom+Nobjectmax))
         write(6,*)'You are not reading from an objectfile!'
         write(6,*)'Assuming having only molecules, and one medium'
         objecttype=0
         medeps(1)=repsin/epkt
         dataobject(1,1)="is a molecule    0"
         tmpiatmmed(1)=imedianumb
      end if
      do ii=1,nobject
         ! regardless of being a molecule or not, iatmmed has a field to say
         ! which is its medium
         iatmmed(Natom+ii)=tmpiatmmed(ii)        
      end do
      
      !-----------------------------------------------------------------
      ! 2011-05-01 End new piece of code
      !-----------------------------------------------------------------
      end subroutine getatm
