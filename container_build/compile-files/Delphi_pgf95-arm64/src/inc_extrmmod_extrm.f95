      !  ATTENTION!  This file is part of extrmmod module.
      !==================================================================
      !-------------------------------------------------------------
      ! 2011-05-08  Arguments are not necessary, transfered via modules
      subroutine extrm
      
      ! 2011-05-08 Arrays and variables are declared in pointers module  
      !!c b++++++++++++++++++++++++++++++++++++++++
      !!c e++++++++++++++++++++++++++++++++++++++++
      
      ! 2011-05-08 atpos replaced by delphipdb%xyz
      !-------------------------------------------------------------------
      ! 2011-06-14 Declarations added due to IMPLICIT NONE
      integer :: ix,ii,ios,i
      real :: frcnum
      !-------------------------------------------------------------------
      
      !!c find extrema and calculate scale according to them and
      !!c to the percent box fill
      ! 2011-05-08  cmin and cmax are now variables of coord type                
      cmin=coord(6000.,6000.,6000)
      cmax=coord(-6000.,-6000.,-6000)
      
      do ix=1,natom
         ! 2011-05-07  Using new operations on coord type variables  defined
         cmin=min(cmin,delphipdb(ix)%xyz-delphipdb(ix)%rad3)
         cmax=max(cmax,delphipdb(ix)%xyz+delphipdb(ix)%rad3)
      end do
      !!c b++++++++++++++++++++
      !!c find global extrema, both, molecule and objects, are considered
      do ii=1,nobject
         strtmp=dataobject(ii,1)
         if (strtmp(1:4).ne.'is a') then
            ! 2011-05-07  Using new operations on coord type variables  defined
            cmin=min(cmin,limobject(ii)%min)
            cmax=max(cmax,limobject(ii)%max)
         end if
      end do
      return
      end subroutine extrm
	  
      !----------------------------------------------------------------------
      ! 2011-05-09 Hash table is not necessary in module architecture
      subroutine off(oldmid,pmid)
      ! 2011-05-09 Replaced by coord type variables
      type(coord) :: oldmid,pmid,summid,tempmid
      character(80) ::  line,fn
      character(6) :: head
      character(24) :: crdstr
      integer :: un
      !------------------------------------------------------------------------------ 
      ! 2011-06-14 Declarations added due to IMPLICIT NONE
      integer :: ios
      real :: frcnum
      !-------------------------------------------------------------------

      if(iacent) then
         oldmid=acent
         return
      end if
      
      !  2011-05-09  Using operations on coord type variables defined
      oldmid=pmid-offset/scale
      
      !  2011-05-09 Changed to coord type variable
      if((offset%x.eq.999).or.(offset%x.eq.777)) then
         ! 2011-09-05 changed a logical flow of the program
         if(offset%x.eq.999) then
            write(6,*) 'modifying midpoints using frc input file'
            un=15 ; fn=trim(centnam)
         else
            write(6,*) 'modifying midpoints using fort.27'
            un=27 ; fn='fort.27'
         end if
         !  2011-05-09 Changed to coord type variable
         summid=coord(0.,0.,0.)
         frcnum=0.0
         open(un,file=trim(fn),status='old',iostat=ios)
         if(ios.ne.0) then
            write(6,*)'Nonexistence of atom file, for calculating midpoints'
            return 
         end if
         
         do
            read(un,'(a)',iostat=ios) line
            if(ios.ne.0) exit
            head = line(1:6); call up(head,6)
            ! skip header lines
            if((head.ne.'ATOM  ').and.(head.ne.'HETATM')) cycle
            frcnum=frcnum + 1.0
            crdstr = line(31:54)
            read(crdstr,'(3f8.3)')tempmid
            !  2011-05-09  Using operations on coord type variables defined
            summid=summid+tempmid
            if((offset%y.eq.999).or.(offset%y.eq.777)) exit
         end do
         
         if(frcnum.gt.0.0) then
            !  2011-05-09  Using operations on coord type variables defined
            oldmid=summid/frcnum
         else
            if(offset%x.eq.999) &
            &write(6,*) 'frc file empty of atoms for'
            if(offset%x.eq.777) &
            &write(6,*) 'unit 27 empty of atoms for'
            write(6,*) 'midpoint determination therefore &
            &assuming zero offsets'
            !  2011-05-09  Using operations on coord type variables defined
            oldmid=oldmid + offset/scale
         end if
         close(un)
      end if     
      end subroutine off
      
      subroutine warning
      write(6,*)&
      & "WARNING! geometric parameter too small &
      &to be 'resolved' by the program!"   
      end subroutine warning
