      !  ATTENTION!  This file is part of qqinttmod module.
      !==================================================================
      subroutine yesno(type,ieq1,quant,ieq2,itf)
      
      character(20) :: quant      ! changed to F95 syntax on 2011-04-17
      character(30) :: type       ! changed to F95 syntax on 2011-04-17
      logical :: itf              ! changed to F95 syntax on 2011-04-17
      logical :: flag=.false.     ! added to avoid GOTO statements on 2011-04-17
      integer :: ieq1,ieq2
      integer :: j,k
      !-------------------------------------------------------------------
      
      itf=.false.
      j=index(quant,"TRUE")
      if(j.ne.0) then
         itf=.true.
         if(quant(1:ieq2).ne."TRUE") call spuchar(quant,type,itf) 
         ! 2011-04-17 removed goto 80 and goto 90
         return                                     ! replacing by if call to subroutine with 
      end if                                        ! message about spurious character
      
      
      j=index(quant,"YES")
      if(j.ne.0) then
         itf=.true.
         if(quant(1:ieq2).ne."YES") call spuchar(quant,type,itf)
         return
      end if
      
      j=index(quant,"ON")
      if(j.ne.0) then
         itf=.true.
         if(quant(1:ieq2).ne."ON") call spuchar(quant,type,itf)
         return
      end if
      
      j=index(quant,"T")
      if(j.ne.0) then
         itf=.true.
         if(quant(1:ieq2).ne."T") call spuchar(quant,type,itf)
         return
      end if
      
      k=index(quant,"FALSE")
      if(k.ne.0) then
         itf=.false.
         if(quant(1:ieq2).ne."FALSE") call spuchar(quant,type,itf)
         return
      end if
      
      k=index(quant,"OFF")
      if(k.ne.0) then
         itf=.false.
         if(quant(1:ieq2).ne."OFF") call spuchar(quant,type,itf)
         return
      end if
      
      k=index(quant,"NO")
      if(k.ne.0) then
         itf=.false.
         if(quant(1:ieq2).ne."NO") call spuchar(quant,type,itf)
         return
      end if
      
      k=index(quant,"F")
      if(k.ne.0) then
         itf=.false.
         if(quant(1:ieq2).ne."F") call spuchar(quant,type,itf)
         return
      end if
      
      write(6,*) "!!!!!!!!!!!!!!"
      write(6,*) "Could not assign a value to the field ",quant(1:ieq2)
      write(6,*) "for the statement type ",type(1:ieq1)
      write(6,*) "!!!!!!!!!!!!!!"
      return !  goto 90
      
      end subroutine yesno
      !----------------------------------------------------------------------------------
      !----------------------------------------------------------------------------------              
      subroutine spuchar(quant,type,itf)
      character(20)  :: quant                   
      character(30)  :: type       
      logical :: itf
      
      write(6,*) "!!!!!!!!!!!!!!"
      write(6,*) "Spurious characters found in the field ",trim(quant)
      write(6,*) "for the statement type ",trim(type)
      if(itf) write(6,*) "which was never the less set to be true"
      if(.not.itf) write(6,*) &
      &"which was never the less set to be false"
      write(6,*) "!!!!!!!!!!!!!!"
      end subroutine spuchar
      
