!#####################################################################
!  ATTENTION!  This file is part of crgarrmod module.
!              Do not compile separately!
!#####################################################################

      !==============================================================>
      !find a couple of versors orthogonal to (xa-xb),say xu and xv
      !modul is the intensity of (xa-xb)
      subroutine basisortho(xa,xb,xu,xv,xw,modul)

      type(coord) xb,xa,xu,xv,xw,tmpv
      real a,b,c,tmp,modul

      tmpv=xa-xb
      
      !call diffvect(xa,xb,tmpv)

      if ((abs(tmpv%x).lt.1.e-6).and.(abs(tmpv%y).lt.1.e-6)) then
         xu=coord(0.,0.,0.) ; xv=xu
         modul=abs(xb%z-xa%z)
      else
         modul=tmpv.dot.tmpv
         tmp=modul-tmpv%z*tmpv%z 
         xv%z=-tmp
         tmp=1./sqrt(tmp)
         xu=coord(tmpv%y*tmp, -tmpv%x*tmp, 0.)
         modul=sqrt(modul)
         tmp=tmp/modul
         xv%x=tmpv%x*tmpv%z ; xv%y=tmpv%z*tmpv%y
         xv=xv*tmp
      end if
      
      xw=tmpv/modul
      !call mul(1./modul,tmpvect,xw)
      
      end subroutine basisortho

      !==============================================================>
      !swapping two vectors and their square modulus
      subroutine swap(xa,xb,moda,modb)

      type(coord) :: xb,xa,tmpv
      real :: moda,modb,tmp
        
      tmpv=xa; xa=xb; xb=tmpv

      tmp=moda; moda=modb; modb=tmp

      end subroutine swap
