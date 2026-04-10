      !  ATTENTION!  This file is part of setbcmod module.
      !==================================================================
      subroutine conplt(array,title,iclr,iscl,imk,iplt,&
      & symb,ixmin,ixmax,iplot,ymin,ymax)
      
      !!c produces line plot of array
      integer, parameter :: nxran = 60, nyran=20
      real :: array(nxran)
      character :: line(nxran),symb
      character(70) :: iplot(nyran),title
      !------------- Declarations added due t o IMPLICIT NONE
      integer :: iclr,iscl,imk,iplt,ixmin,ixmax,iymin,iymax
      integer :: iyup,iylw,iyran,ibin,i
      real :: ymin,ymax,yminl,ymaxl,temp,temp1
      !------------------------------------------------------
      line='-'
      !!c------------------------------------------------------
      if(iclr.eq.1) then
         !!c clear plot
         ymin = 1.e10;  ymax = 0.0
         iplot=' '
      end if
      if(iscl.eq.1) then
         
         !!c scale plot
         !!c find max, min of data
         
         do i = 1,nxran
            if(array(i).gt.0.0) then
               ymin = min(ymin,array(i))
               ymax = max(ymax,array(i))
            end if
         end do
         !!c find y plot range in log scale
         yminl = log10(ymin);  ymaxl = log10(ymax)
         iyup = (1. + ymaxl);  iylw = (yminl - 1)
         iyran = iyup - iylw
      end if
      if(imk.eq.1) then
         
         !!c make plot
         !!c stick x values in the appropriate bins after clipping
         yminl = log10(ymin);  ymaxl = log10(ymax)
         iyup = (1. + ymaxl);  iylw = (yminl - 1)
         iyran = iyup - iylw
         do i = 1,nxran
            if((array(i).ge.ymin).and.(array(i).le.ymax)) then
               temp = log10(array(i))
               temp1 = (temp - iylw)/iyran
               ibin = temp1*(nyran - 1) + 1
               iplot(ibin)(i:i)= symb
            end if
         end do
      end if
      if(iplt.eq.1) then
         
         !!c draw out plot
         write(6,*)'  '
         write(6,'(5X,A70)')title
         write(6,'(1x,g9.2,'' |-'',60A,''-|'')')ymax,line
         do i = nyran,1,-1
            write(6,'(11X,''| '',A60,'' |'')')iplot(i)
         end do
         write(6,'(1x,g9.2,'' |-'',60A,''-|'')')ymin,line
         write(6,'(11X,''|'',62X,''|'')')
         write(6,'(6X,I5,58X,I5)')ixmin,ixmax
         202           format(6X,I5,58X,I5)
      end if
      end subroutine conplt
