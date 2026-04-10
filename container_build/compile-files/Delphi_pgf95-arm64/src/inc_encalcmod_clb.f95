!#####################################################################
!ATTENTION!  This file is part of encalcmod module.
!            Do not compile separately!
!#####################################################################

!---------------------------------------------------------------------
      subroutine clb(ergc)

      !2011-06-11 Declarations added due to IMPLICIT NONE
      double precision :: ergc,en,en1
      integer :: izero,i,j
      type(coord) :: txyz
      real :: dist

      izero=0
      en=0.0

      do i=1,nqass-1
         en1=0.0
         do j=i+1,nqass
            txyz=chgpos(i)-chgpos(j); dist=txyz.dot.txyz
            en1 = en1 + atmcrg(j)%value/sqrt(dist)
         end do

         en=en+atmcrg(i)%value*en1
      end do

      ergc=en 

      end subroutine clb

!---------------------------------------------------------------------
      subroutine clbmedia(ergc)

      !2011-06-11 Declarations added due to IMPLICIT NONE
      double precision :: ergc,en,en1
      integer :: izero,i,j
      type(coord) :: txyz
      real :: dist

      en=0.0
      do i=1,nqass
         en1=0.0
         do j=1,nqass
            if (i.ne.j) then
               txyz=chgpos(i)-chgpos(j); dist=txyz.dot.txyz
               en1 = en1 + atmcrg(j)%value/sqrt(dist)
            end if
         end do

         !2011-06-11 Second loop is no longer needed due to IF 
         !           condition in the first loop
         en=en + atmcrg(i)%value*en1/atmeps(i)
      end do
      
      en=en/2.0; ergc=en
      
      end subroutine clbmedia
      
!---------------------------------------------------------------------
      subroutine clbnonl(ergc,ergest,igridout)

      !2011-06-11 Declarations added due to IMPLICIT NONE
      double precision :: ergc,ergest,en,en1,en2
      integer :: izero,i,j,n,igridout
      type(coord) :: txyz
      real :: sqrt8,c,dist

      en=0.0; sqrt8=sqrt(8.0); c=.0006023

      do i=1,nqass
         en1=0.0; en2=0.0
         do j=1,nqass
            if (i.ne.j) then
               txyz=chgpos(i)-chgpos(j); dist=txyz.dot.txyz
               en1 = en1 + atmcrg(j)%value/sqrt(dist)
            end if
         end do

         en=en + atmcrg(i)%value*en1/atmeps(i)

         !calculation for the solvent contribution
         !2011-06-11 Array sout is declared in pointers module and 
         !allocated in nlener subroutine
         if (rionst.gt.1.e-6) then
            do n=1,igridout
               txyz=chgpos(i)-sout(n)%xyz; dist=txyz.dot.txyz
               en2=en2+sout(n)%value/sqrt(dist)
            end do

            ergest=ergest+en2*atmcrg(i)%value
         end if
      end do

      en=en/2.0; ergest=ergest*c/(2.0*epsout); ergc=en

      if(allocated(sout))deallocate(sout)

      end subroutine clbnonl

!---------------------------------------------------------------------
      subroutine clbtot(ergest,ergc)

      real :: cutedgesi,cutedgesj,cutedgesk

      !2011-06-11 Declarations added due to IMPLICIT NONE
      double precision :: ergc,ergest,en,en1,en2
      integer :: izero,i,j,k,n,igridout
      type(coord) :: txyz,gxyz
      real :: goff,dist,dist1,c,carica,phi,tmp

      allocate(sout(igrid*igrid*igrid))

      en=0.0; n=0
      goff=(igrid+1.)/2.
      gxyz=(-goff/scale)+oldmid

      c=scale*scale*scale; tmp=-2.*rionst/c

      do k=1,igrid
         cutedgesk=1.
         if (k.eq.1.or.k.eq.igrid) cutedgesk=.5
         
         do j=1,igrid
            cutedgesj=cutedgesk
            if (j.eq.1.or.j.eq.igrid) cutedgesj=cutedgesk*.5
  
            do i=1,igrid
               cutedgesi=cutedgesj
               if (i.eq.1.or.i.eq.igrid) cutedgesi=cutedgesj*.5 
  
               if (idebmap(i,j,k)) then
                  phi=phimap(i,j,k); carica=phi*tmp

                  !if the gp is in solution and the contribution is 
                  !higher than a threshold, then put this 
                  !information in a list
                  n=n+1
                  sout(n)%xyz=(float(int_coord(i,j,k))/scale)+gxyz
                  sout(n)%value=carica*cutedgesi
               end if
            end do
         end do
      end do
      
      igridout=n
      write(6,*)'number of g.p. in solution contributing to the energy',igridout

      !2011-06-11 Re-sizing of array seems to be unnecessary because
      !           of short life of that array (deallocation at the end
      !           of this subroutine)
      do i=1,nqass
         en1=0.0; en2=0.0
         do j=1,nqass
            if (i.ne.j) then
               txyz=chgpos(i)-chgpos(j); dist=txyz.dot.txyz
               en1 = en1 + atmcrg(j)%value/sqrt(dist)
            end if
         end do

         en=en+atmcrg(i)%value*en1/atmeps(i)

         !calculation for the solvent contribution
         do n=1,igridout
            txyz=chgpos(i)-sout(n)%xyz; dist=txyz.dot.txyz !test
            en2=en2+sout(n)%value/sqrt(dist)
         end do

         ergest=ergest+en2*atmcrg(i)%value
      end do

      en=en/2.0; ergest=ergest*.0006023/(2.0*epsout); ergc=en
      deallocate(sout)

      end subroutine clbtot
