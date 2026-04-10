      !  ATTENTION!  This file is part of rdhmod module.
      !==================================================================
      !---------------------------------------------------------------------------
      ! read radius file and store in hash table
      subroutine rdhrad     
      
      character(80) :: filnam
      character(80) :: comment
      logical pK
      !------------------------------------------------------------------------------ 
      integer :: io,inc
      real :: rad=0.
      !------------------------------------------------------------------------------ 
      !    siznam declared in qlog and assigned in defprm (module qinttot)  
      !    sizlen is not necessary due to trim intrinsic function
      open(11,file=trim(siznam),status='old',iostat=io)
      if(io.ne.0) then     ! 2011-04-23 instead of err=901
         call mess1('radius') ; return  
      end if
      
      write(6,*)'atom radii read from file ', trim(siznam)
      write(6,*)' '
      
      ! skip/print comments (start with !) and one column header line 
      ! 2011-04-23   changed to F95 style to get rid of GOTO statement
      do
         read(11,'(a)',iostat=io )  comment  ! 2011-04-23 removed
         if(io.ne.0) then    ! 2011-04-23 instead of end=901
            call mess1('radius') ; return  
         end if
         if(comment(1:1).ne.'!') exit
         write(6,*)comment
      end do
      
      inc=index(comment,' ')-1
      if(comment(1:16).eq.'atom__res_radius')then
         pK=.false.
      elseif(comment(1:inc).eq.'atom__resnumbc_radius_')then
         write(6,*)'reading pK style radius file'
         pK=.true.
      else
         write(6,*)'unknown header format in radius file:'
         write(6,*)comment
         stop
      endif
      !    2011-04-23   defines number of records in file and allocates 
      nrmax = 0   
      do
         read(11,'(a)',iostat=io )  comment  
         if (io.ne.0) exit
         nrmax=nrmax+1
      end do
      allocate(radii(nrmax),rhash(nrmax))
      
      !   2011-04-23   Rewinds file, skips comments and header line
      rewind(11); nrdrec=0
      do
         read(11,'(a)')  comment
         if(comment(1:1).ne.'!') exit
      end do
      !    2011-04-23   Reads actual data from file 
      DREAD: do
         nrdrec = nrdrec + 1
         ! 2011-04-23   unnecessary due to allocatable arrays
         if(pK)then
            ! 2011-04-23   removed 202 format label and err=904 end=300 labels 
            read(11,'(A6,A3,A4,A1,F8.3)',iostat=io)&
            &atm,res,rnum,chn,rad
            if(io.ne.0) exit DREAD   !  instead of end=300 label
            if(atm.ne.' ') then
               call up(atm,6); call elb(atm,6)
            end if
            if(res.ne.' ') then
               call up(res,3); call elb(res,3)
            end if
            if(rnum.ne.' ') then
               call up(rnum,4); call elb(rnum,4)
            end if
            if(chn.ne.' ') then
               call up(chn,1);  call elb(chn,1)
            end if
            radii(nrdrec)%atnam = atm
            radii(nrdrec)%rnam  = res
            radii(nrdrec)%rnum  = rnum
            radii(nrdrec)%chn   = chn
            radii(nrdrec)%value = rad
            
         else
            read(11,'(A6,A3,F8.3)',iostat=io)atm,res,rad
            if(io.ne.0) exit DREAD   !  instead of end=300 label
            if(atm.ne.' ') then
               call up(atm,6);  call elb(atm,6)
            end if
            if(res.ne.' ') then
               call up(res,3);  call elb(res,3)
            end if
            radii(nrdrec)%atnam = atm
            radii(nrdrec)%rnam  = res
            radii(nrdrec)%rnum  = '    '
            radii(nrdrec)%chn   = ' '
            radii(nrdrec)%value  = rad
         endif
         ! 2011-04-24   Added three last parameters to make rent
         ! 2011-07-18 Has tables are not necessary anymore        
      end do DREAD
      904          close(11)
      nrdrec = nrdrec - 1
      write(6,*)'# of radius parameter records:',nrdrec
      end subroutine rdhrad