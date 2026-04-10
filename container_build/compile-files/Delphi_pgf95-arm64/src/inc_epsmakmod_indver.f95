!#####################################################################
!ATTENTION! This file is part of epsmakmod module.
!           Do not compile separately!
!
!2011-05-24 Changed to coord type variable, mnxyz and mxxyz, which
!           are accessible through declaration in pointers module
!#####################################################################
 
!---------------------------------------------------------------------
      subroutine indverdata(prbrad,scale)

      implicit none

      real :: jig
      !2011-05-27 Declaration added due to IMPLICIT NONE
      real :: prbrad,scale,cba
      integer :: idtemp

      !cba defines scale, cba and jig the box limits..
      !note that cba is the minumum box size. it can be bigger if need 
      !be. i.e. can rescale to a larger box if need be.
      cba=prbrad

      jig=2.0*(1.0/scale)

      !2011-05-24  Using operations on coord and int_coord type 
      !variables defined in module operators_on_coordinates 
      mnxyz=cmin-(rdmx+prbrad)
      mxxyz=cmax+(rdmx+prbrad)

      !accessible surface points are prbrad away from the actual surface
      !mid points are at most (1/scale) away from the actual surface.
      !midpoints must never be in box zero. accessible point can be in 
      !box zero but not in box -1.
      !
      !e.g. mnx+prbrad==van der Waals surface
      !vanderwaals surface - (1/scale)= minumum midpoint position.
      !therefore mnx+prbrad-(1/scale) gt 1
      !and mnx gt 0. the latter is always true since cba is ge. prbrad.
      !therefore we have...

      mnxyz=mnxyz-jig; mxxyz=mxxyz+jig

      !2011-05-24 Cycle is introduced to remove GOTO 100 statement
      do
         mxxyz=mxxyz+(cba-prbrad); mnxyz=mnxyz-(cba-prbrad)
         lmncb1=int(((mxxyz-mnxyz)/cba)+1.)

         !if points are too widely seperated for the scale and idmax 
         !then rescale..

         !2011-05-24 idmax is declared as parameter in qlog module
         if (lmncb1.vorgt.idmax) then
            idtemp=max(lmncb1)
            cba=cba*real(idtemp+1)/idmax
            write(6,*) "initial cube size too small, "
            write(6,*) "in assigning accessible points to a grid"
            write(6,*) "therefore rescaling..."

            !2011-05-24 CYCLE and EXIT statements instead of GOTO 
            cycle
         else
            exit
         end if
      end do
      
      lcb1=lmncb1%i; mcb1=lmncb1%j; ncb1=lmncb1%k 

      !grdi is just the inverse of cba..,used more...
      grdi=1.0/cba

      end subroutine indverdata

!---------------------------------------------------------------------
      subroutine indver(extot1,iab1,iab2)

      implicit none

      !program to compile the lists iab1,iab2 and icume for use in
      !nearest vertex work. iexpos are box coordinates of vertices
      !and dont need to be passed if comaprisons are done with real  
      !angstroms..
      !but its often more convenient to do so since grid points have to 
      !be converted anyway...

      !2011-05-24 Arrays are asseccible via declaration in pointers 
      !module and allocation in calling vwtms subroutine
      integer :: extot1
      integer iab1(0:lcb1,0:mcb1,0:ncb1),iab2(0:lcb1,0:mcb1,0:ncb1)
      !2011-05-27 Declarations added due  to IMPLICIT NONE
      integer :: i,n,ix,iy,iz,k,j
      
      !2011-05-24 Array iexpos now is of int_coord type and thus 1D
      allocate(iexpos(extot1))

      !initialize grid..
      !2011-05-24 Changed to array operation           
      iab1=1; iab2=0

      !make linear arrays for expos
      !find the number of points in each box, put in iab2, make iexpos
      do i=1,extot1
         !2011-05-24  Using operations on coord and int_coord type 
         !variables defined in module operators_on_coordinates 
         iexpos(i)=int( (expos(i)-mnxyz)*grdi )
         iab2(iexpos(i)%i,iexpos(i)%j,iexpos(i)%k)= &
                         & iab2(iexpos(i)%i,iexpos(i)%j,iexpos(i)%k)+1
      end do

      !check each box for occupancy, using fill number to mark out
      !space in icume, using
      n=0
      do i=0,lcb1
         do j=0,mcb1
            do k=0,ncb1
               !if the box is not empty put start position to n+1 in 
               !iab1 end to n+box occupancy in iab2, overwriting 
               !occupancy..
               if (iab2(i,j,k).ne.0) then
                  iab1(i,j,k)=n+1
                  n=n+iab2(i,j,k)
                  iab2(i,j,k)=n
               end if
            end do
         end do
      end do

      !fill icume using iab1 and iab2, note that iab1 is used to hold 
      !the position in icume, therefore needs to be recalculated..
      do i=1,extot1
         ix=iexpos(i)%i; iy=iexpos(i)%j; iz=iexpos(i)%k
         j=iab1(ix,iy,iz)
         icume(j)=i
         iab1(ix,iy,iz)=iab1(ix,iy,iz)+1
      end do

      !reset iab1 for use in inner loop
      do i=1,extot1
         ix=iexpos(i)%i; iy=iexpos(i)%j; iz=iexpos(i)%k 
         iab1(ix,iy,iz)=iab1(ix,iy,iz)-1
      end do

      !icume now contains pointers to each dot inside a particular box
      !and each box has 2 pointers into icume.,a start pointer and a 
      !finish pointer
      !use however you want...

      !2011-05-24 Array deallocation is done with DEALLOCATE statement
      if(allocated(iexpos)) deallocate(iexpos)

      end subroutine indver
