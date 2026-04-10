      !  ATTENTION!  This file is part of qqinttmod module.
      !==================================================================
      subroutine ppi(line,mlen)
      ! changed to F95 syntax on 2011-04-16
      character(400) ::  mline,line             
      character(4) :: alph1="[{]}",alph2="(())" 
      character  ::  chr                        
      integer :: comnum,comb(50),come(50)       
      logical :: iuse,istat,icom(50)            
      !-------------------------------------------------------------------
      ! 2011-06-14 Declarations added due to IMPLICIT NONE
      integer :: mlen,ml,bnum1,bnum2,branum,i,j
      !-------------------------------------------------------------------
      !!c prepare line, removing blanks, changing brackets
      !!c if we run into an explanation mark then turn off the character
      !!c grabber until we reach another one. also ignore blanks
      
      if(mlen.eq.0) return ! 2011-04-16 instead of goto 99
      line(mlen+1:400)=' '
      ml=0 ; iuse=.true.
      do i=1,mlen
         chr=line(i:i)
         if(chr.eq."!") then
            iuse=.not.iuse ; cycle !  2011-04-16 instead of goto 50
         end if
         if(iuse) then
            if(chr.ne." ") then
               ml=ml+1
               j=index(alph1,chr)
               if(j.ne.0) chr=alph2(j:j)
               mline(ml:ml)=chr
            end if
         end if
      end do
      !!c b++++++++++++++++++++++++++++++there was a bug
      if (ml.gt.0.and.mline(ml:ml).eq.",")  ml=ml-1   
      !!c e++++++++++++++++++++++++++++++++
      !!c try and interpret the line. look for statement (= signs), or
      !!c functions ( a closed pair of brackets.)
      !!c if there is trouble, e.g. an unpaired bracket, then print
      
      branum=0
      
      !!c branum = bracket number
      !!c comb(i) = beginning of command i
      !!c come(i) = end of command i
      
      j=1 ;  comb(1)=1
      do i=1,ml
         chr=mline(i:i)
         if(chr.eq."(") branum=branum+1
         if(chr.eq.")") branum=branum-1
         if(((chr.eq.",").or.(chr.eq.":").or.(chr.eq."|")).and.(branum&
         & .eq.0)) then
            come(j)=i-1 ;  j=j+1 ;  comb(j)=i+1
         end if
      end do
      come(j)=i-1 ; comnum=j
      
      ! comnum = number of commands
      ! icom is a logical array, true if the corresponding command is
      ! a statement. once determined, pass to the statement or function
      ! interpreter
      
      do i=1,comnum
         bnum1=0 ; bnum2=0
         istat=.false. ; icom(i)=.false.
         do j=comb(i),come(i)
            if(mline(j:j).eq."(") bnum1=bnum1+1
            if(mline(j:j).eq.")") bnum2=bnum2+1
            if(mline(j:j).eq."=".and.(bnum1.eq.bnum2) ) istat=.true.
         end do
         if(istat) icom(i)=.true.
      end do
      
      do i=1,comnum
         if((come(i)-comb(i)).lt.1) cycle    
         if(icom(i))then
            call statint(mline(comb(i):come(i)),come(i)-comb(i)+1)
         else
            call funcint(mline(comb(i):come(i)),come(i)-comb(i)+1)
         end if
      end do
      
      end subroutine ppi