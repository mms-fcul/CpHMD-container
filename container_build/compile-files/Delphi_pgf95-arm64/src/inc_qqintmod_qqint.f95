      !  ATTENTION!  This file is part of qqinttmod module.
      !==================================================================
      subroutine qqint(i,argnam)
      !   2011-04-14 Argument arglen is obsolete
      !  removed include "qlog.h" 2011-04-14
      integer :: arglen
      character(len=80) :: argnam
      !-------------------------------------------------------
      ! 2011-06-14 Declarations added due to IMPLICIT NONE
      integer :: i
      !-------------------------------------------------------
      call defprm
      if(i.gt.0)then
         prmnam=trim(argnam); prmlen=len_trim(argnam)
      end if
      
      if(i.gt.1)then
         write(6,*) " "
         write(6,*) "WARNING!"
         write(6,*) " too many file names.."
         write(6,*) "only using the first as the parameter file"
         write(6,*) "file name= ",prmnam(:prmlen)
         write(6,*) " "
      end if
      
      call rdprm
      !!c  parameter assessment now because too complex to do elsewhere
      if (radprb(2).eq.-1.) radprb(2)=radprb(1)
      end subroutine qqint  ! changed syntax 2011-04-14