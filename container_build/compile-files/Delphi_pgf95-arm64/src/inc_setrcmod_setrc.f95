      !  ATTENTION!  This file is part of setrcmod module.
      !==================================================================
      !   Argument to subroutine are now passed via module qlog
      subroutine setrc
      
      !  2011-04-28   Removed due to the module architecture
      !	include "qdiffpar4.h"
      !	include "qlog.h"
      
      !!c b++++++++++++++++++++
      !!c e++++++++++++++++++++
      !-------------------------------------------------------------------
      ! 2011-06-14 Declarations added due to IMPLICIT NONE
      integer :: i
      real :: rad,chrgv
      !-------------------------------------------------------------------
      write(6,*) 'assigning charges and radii...'
      write(6,*) 
      natom=0
      !!c begin read
      !!c b++++++++++++++++++++
      !    2011-04-28   Other arguments to subroutine are now passed 
      call getatm(pdbnam,pdblen)
      ! 2011-04-28   Obsolete due to allocatable arrays
      if(.not.iatrad) then
         do i=1,Natom
            ! 2011-05-02  Changed to derived-type array delphipdb
            atm = delphipdb(i)%atinf(1:5)
            res = delphipdb(i)%atinf(7:9)
            rnum = delphipdb(i)%atinf(12:15)
            chn = delphipdb(i)%atinf(11:11)
            if(atm.ne.' ') then
               call up(atm,6) ;  call elb(atm,6)
            end if
            if(res.ne.' ') then
               call up(res,3) ;  call elb(res,3)
            end if
            if(rnum.ne.' ') then
               call up(rnum,4);  call elb(rnum,4)
            end if
            if(chn.ne.' ') then
               call up(chn,1);   call elb(chn,1)
            end if
            
            ! assign radius, searching for decreasingly specific specification
            ! ending with generic atom type
            ! note all atoms must have an assignment
            
            !  2011-04-28   atm, res, rnum and chn variables declared in pointers module
            !  2011-05-03  radass and crgass subroutinse are now combined
            call ass(rad,radii,rhash,Nrmax)
            
            !  2011-05-03 Assignment to rad=0 if record is not found is made in
            if(norecord.eq.1) then
               !2011-04-28  changed a15 format specification 
               write(6,'(''!!! WARNING: no radius record for '',a)')&
               &trim(delphipdb(i)%atinf)
               ! need stop here
               !
            elseif(rad.lt.1.e-6.and.(atm(1:1).ne.'H'.and.atm(2:2).ne.'H')) then 
               ! 2011-04-28  changed a15 format specification 
               
               write(6,'(''!!! WARNING: radius of heavy atom'',a,''&
               & is set to zero'')') trim(delphipdb(i)%atinf)
            end if
            
            ! store rad,xn in rad3,  for later use
            !  2011-04-28, Array declared in module pointers and allocated in getatm subroutine.
            delphipdb(i)%rad3=rad       
            ! scale and assign charge to grid
            
            ! 2011-04-28 Other arguments to subroutine are now passed via module pointer 
            call ass(chrgv,charge,chash,Ncmax)
            !  2011-04-28, Array declared in module pointers and allocated in getatm subroutine.
            delphipdb(i)%chrgv4=chrgv
            !debug                write(6,*)'--->',i,Natom,rad,chrgv
         end do
         
         ! write record to new coordinate file if required, with
         ! occupancy and temperature factor fields replaced by radius and charge
         
      end if
      
      ! check charge assignments (Sri April 2, 1996)
      !    2011-04-28   Arguments to subroutine are now passed via modules pointer and qlog
      call chkcrg
      ! unformatted write...
      if(ipdbwrt) then
         !  2011-05-02   updblen obsolete, used trim instead
         open(20,file=trim(updbnam),form='unformatted')
         do i=1,natom
            !  2011-04-28    Changed to derived-type array delphipdb declared 
            write(20)  delphipdb(i)%xyz,delphipdb(i)%rad3,&
            &delphipdb(i)%chrgv4
         end do
         close(20)
      end if
      ! 2011-04-28 iatout declared in qlog module 
      if(iatout) then
         ! 2011-04-28 used trim intrinsic function to get rid of trailing spaces  
         open(19,file=trim(mpdbnam))
         filnam = ' '
         inquire(19,name = filnam)
         ! 2011-04-28 used trim intrinsic function to get rid of trailing spaces  
         write(6,*)&
         & 'atomic coordinates, charges and radii written to file ',&
         & trim(filnam)
         write(6,*)
         if (mpdbfrm.eq.0) then
            write(19,*)'DELPHI PDB FILE'
            write(19,*)'FORMAT = 1', mpdbfrm
            write(19,*)'HEADER output from qdiff'
            write(19,*)'HEADER atom radii in columns 55-60'
            write(19,*)'HEADER atom charges in columns 61-67'
            ! 2011-04-28  Joined two pieces under the same if condition
            line=' '
            line(1:11)="ATOM       "
            do i=1,natom
               !  2011-04-28    Changed to derived-type array delphipdb declared 
               ! 2011-04-28    Removed format labels 205 and 206
               write(line(7:11),'(i5)')i
               line(12:30)=delphipdb(i)%atinf
               
               write(line(55:67),'(F6.2,F7.3)')&
               &delphipdb(i)%rad3,delphipdb(i)%chrgv4
               write(line(31:54),'(3f8.3)') delphipdb(i)%xyz
               line(27:30)='    '
               write(19,'(a)')trim(line)
            end do
         elseif(mpdbfrm.eq.1) then
            write(19,*)'DELPHI PDB FILE'
            write(19,*)'FORMAT = PQR'
            write(19,*)'HEADER output from qdiff'
            write(19,*)'HEADER atom charges in columns 56-61'
            write(19,*)'HEADER atom radii   in columns 63-68'
            ! 2011-04-28  Joined two pieces under the same if condition
            line=' '
            line(1:6)="ATOM  "
            do i=1,natom
               write(line(7:11),'(i5)')i
               line(12:30)=delphipdb(i)%atinf
               write(line(31:54),'(3f8.3)') delphipdb(i)%xyz
               line(55:55)=' '
               write(line(56:61),'(f6.3)')delphipdb(i)%chrgv4
               line(62:62)=' '
               write(line(63:68),'(f6.3)')delphipdb(i)%rad3
               
               !	  write(crgstr,207)chrgv,rad
               !	  line(55:70) = crgstr
               write(19,'(a)')trim(line)
            end do
         end if
         
         
         if (mpdbfrm.eq.40) then
            write(19,*)'DELPHI PDB FILE'
            write(19,*)'FORMAT = 1', mpdbfrm
            write(19,*)'4 digits precison'
            write(19,*)'HEADER output from qdiff'
            write(19,*)'HEADER atom radii in columns 55-61'
            write(19,*)'HEADER atom charges in columns 62-69'
            ! 2011-04-28  Joined two pieces under the same if condition
            line=' '
            line(1:11)="ATOM       "
            do i=1,natom
               write(line(7:11),'(i5)')i
               line(12:30)=delphipdb(i)%atinf
               
               write(line(55:69),'(F7.4,F8.4)')&
               &delphipdb(i)%rad3,delphipdb(i)%chrgv4
               write(line(31:54),'(3f8.3)') delphipdb(i)%xyz
               line(27:30)='    '
               write(19,'(a)')trim(line)
            end do
         elseif(mpdbfrm.eq.41) then
            write(19,*)'DELPHI PDB FILE'
            write(19,*)'FORMAT = PQR'
            write(19,*)'4 digits precison'
            write(19,*)'HEADER output from qdiff'
            write(19,*)'HEADER atom charges in columns 56-61'
            write(19,*)'HEADER atom radii   in columns 63-68'
            ! 2011-04-28  Joined two pieces under the same if condition
            line=' '
            line(1:6)="ATOM  "
            do i=1,natom
               write(line(7:11),'(i5)')i
               line(12:30)=delphipdb(i)%atinf
               write(line(31:54),'(3f8.3)') delphipdb(i)%xyz
               write(line(55:62),'(f8.4)')delphipdb(i)%chrgv4
               write(line(63:69),'(f7.4)')delphipdb(i)%rad3
               
               write(19,'(a)')trim(line)
            end do
         end if
         !!c ++++++++++++++WWW++++++++++
         close(19)
      end if
      !  2011-04-28, label 903 removed
      if((natom.eq.0).and.(nobject.eq.0)) then
         write(6,*) &
         &"exiting due to non-existence of atom file nor object data"
         return 
      end if
      end subroutine setrc