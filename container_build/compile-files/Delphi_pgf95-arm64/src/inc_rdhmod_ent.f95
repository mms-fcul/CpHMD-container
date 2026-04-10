      !  ATTENTION!  This file is part of rdhmod module.
      !==================================================================
      !---------------------------------------------------------------------------
      subroutine ent(atm,res,rnum,chn,nent,hash,Nhash,itot)
      character(6) :: atm      !  atom name 
      character(3) :: res      !  residue name 
      character ::  chn        !  chain name 
      character(4) ::  rnum    !  residue number 
      integer :: nent,Nhash,itot,n,new
      !------------------------------------------------------------------------------ 
      
      
      ! enter character strings res,atm into hash table for radii
      ! by assigning it a number nent
      
      ! check to see if there is room
      ! 2011-04-24  Not necessary due to allocatable array
      
      ! 2011-04-24  either radii or charge hash table 
      type(hash_table) :: hash(Nhash)
      
      n = ichash(atm,res,rnum,chn,Nhash)
      if(hash(n)%inumb.ne.0) then
         
         ! slot filled
         ! run down linked list
         ! 2011-04-24  new F95 code, removed labels 9000 and 9001
         do
            if(hash(n)%ilink.eq.0) exit
            n = hash(n)%ilink
         end do
         !  search for empty slot
         
         new = 1
         !  2011-04-24 new F95 code, removed labels 9002 and 9003 
         do
            if(hash(new)%inumb.eq.0) exit
            new = new + 1
         end do
         ! found one- add link
         
         hash(n)%ilink = new
         n = new
      end if
      
      ! slot empty
      ! fill slot
      
      hash(n)%inumb = nent
      hash(n)%ilink = 0
      itot = itot + 1
      end subroutine ent
      !---------------------------------------------------------------------------
      subroutine mess1(s)   !  2011-04-23 created instead err=901
      character(*) :: s
      write(6,*) "nonexistence or unexpected end of ",trim(s)," file"
      write(6,*)"This is correct in case there are ONLY objects,"
      write(6,*)"or in case some specific delphi pdb format is used!"
      end subroutine mess1