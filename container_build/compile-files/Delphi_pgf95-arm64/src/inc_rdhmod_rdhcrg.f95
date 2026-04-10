      !  ATTENTION!  This file is part of rdhmod module.
      !==================================================================
      ! read charge file and stores in hash table
      subroutine rdhcrg
      !------------------------------------------------------------------------------ 
      ! 2011-06-14 Declarations added due to IMPLICIT NONE
      integer :: io,i
      real :: chrgv
      !------------------------------------------------------------------------------ 
      character(80) :: filnam
      character(80) :: comment
      !    crgnam declared in qlog and assigned in defprm (module qinttot)  
      !    crglen is not necessary due to trim intrinsic function
      
      open(12,file=trim(crgnam),status='old',iostat=io)
      if(io.ne.0) then     ! 2011-04-24 instead of err=901
         call mess1('charge') ; return  
      end if
      
      write(6,*)' '
      write(6,*)'atomic charges read from file ', trim(crgnam)
      write(6,*)' '
      ! skip/print comments (start with !) and one column header line 
      ! 2011-04-24   changed to F95 style to get rid of GOTO statement
      
      do
         read(12,'(a)',iostat=io )  comment  ! 2011-04-24 removed
         if(io.ne.0) then                ! 2011-04-24 instead of end=901
            call mess1('charge') ; return  
         end if
         if(comment(1:1).ne.'!') exit
         write(6,*)comment
      end do
      !    2011-04-24   defines number of records in file and allocates 
      ncmax = 0   
      do
         read(12,'(a)',iostat=io )  comment  
         if (io.ne.0) exit
         ncmax=ncmax+1
      end do
      allocate(charge(ncmax),chash(ncmax))
      if(.not.allocated(charge)) allocate(charge(ncmax))
      if(.not.allocated(chash)) allocate(chash(ncmax))
      
      !   2011-04-23   Rewinds file, skips comments and header line
      rewind(12)
      do
         read(12,'(a)')  comment
         if(comment(1:1).ne.'!') exit
      end do
      !    2011-04-24   Reads actual data from file 
      
      nchrec = 0
      101        continue
      do
         nchrec = nchrec + 1
         ! 2011-04-24   unnecessary due to allocatable arrays
         read(12,'(A6,A3,A4,A1,F8.5)',iostat=io)&
         &atm,res,rnum,chn,chrgv
         if(io.ne.0) exit   !  instead of end=300 label
         if(atm.ne.' ') then
            call up(atm,6);  call elb(atm,6)
         end if
         if(res.ne.' ') then
            call up(res,3);  call elb(res,3)
         end if
         if(rnum.ne.' ') then
            call up(rnum,4); call elb(rnum,4)
         end if
         if(chn.ne.' ') then
            call up(chn,1);  call elb(chn,1)
         end if
         charge(nchrec)%atnam = atm
         charge(nchrec)%rnam  = res
         charge(nchrec)%rnum  = rnum
         charge(nchrec)%chn   = chn
         charge(nchrec)%value = chrgv
         !	  write(6,*) atm,res,rnum,chn,nchrec
      end do
      close(12)
      nchrec = nchrec - 1 
      write(6,*)'# of charge parameter records:',nchrec
      end subroutine rdhcrg
