!#####################################################################
!ATTENTION! This file is part of epsmakmod module.
!           Do not compile separately!
!
!2011-05-26 Mismatch in dimension for egrid array (was 4D in the 
!           calling subroutine)
!           In order to transfer arrays via pointers module
!           names of the following arrays changed:
!           ivert -> vindx , vc -> vert
!2011-05-26 Seems parameters are not used
!parameter (mxtri=200000)
!parameter (mxvtx=100000)
!2011-05-26 Arrays below are transferred via pointers module
!           integer egrid(igrid,igrid,igrid), e2(257,257,2)
!           integer ivert(150000),ipv(5,129,129)
!           dimension vc(3,100000)
!           integer ivert(1),ipv(5,257,257)
!           real vc(3,1) 
!           real oldmid(3)
!#####################################################################

      subroutine ex(ivn,iv,vdat1,lvdat1)
      
      implicit none

      integer ::  e2(257,257,2),ipv(5,257,257)
      integer :: itind(0:256),itrn(2208)
      type(int_coord) :: ivslb(32000)
      character(20) :: uplbl
      character(10) :: nxtlbl
      character(60) :: toplbl
      character(16) :: botlbl
      character(80) :: vdat1,vdat2
      !2011-05-26 Declaration of variables due to implicit none
      integer :: ix,iy,iz,ivn,iv,lvdat1,i1,i2,i,j,k,n,ibox,indx
      real :: fourth,x2,y2,z2

      !era un 4 bytes
      !vdat2(:lvdat1+1)=vdat1(:lvdat1+1)
      !vdat2(lvdat1+2:lvdat1+7)="v3.dat"

      fourth=1.0/4.0

      open(11,file="/usr/local/bin/v3.dat",form="unformatted")
      read(11)itind,itrn
      close(11)

      !loop for times
      ipv=0

      do iy=1,igrid
         do ix=1,igrid
            !2011-05-26 Just added missing dimension to egrid, not sure 
            !if it absolutely correct from algorithmic standpoint
            if (egrid(ix,iy,1)%i.le.0) then
               e2(ix,iy,2)=0
            else
               e2(ix,iy,2)=1
            end if
         end do
      end do

      ibox=0; iv=0; ivn=0

      do iz=1,igrid-1
         do iy=1,igrid
            do ix=1,igrid
               e2(ix,iy,1)=e2(ix,iy,2)
     
               !2011-05-26 Just added missing dimension to egrid, not 
               !sure if it absolutely correct from algorithmic 
               !standpoint
               if (egrid(ix,iy,iz+1)%i.le.0) then
                  e2(ix,iy,2)=0
               else
                  e2(ix,iy,2)=1
               end if
            end do
         end do

         k=0
         do iy=1,igrid-1
            j=e2(1,iy,1)+e2(1,iy+1,1)+e2(1,iy+1,2)+e2(1,iy,2)
            do ix=1,igrid-1
               i2= e2(ix+1,iy  ,1) + e2(ix+1,iy+1,1)+&
                                   & e2(ix+1,iy+1,2) + e2(ix+1,iy,  2)
               i=i2+j; j=i2
               
               if ((i.ne.0).and.(i.ne.8)) then
                  k=k+1
 
                  !determine index of box, 1,254
                  !it was found to be better to calculate the index 
                  !NOW, not later
                  indx=     e2(ix,  iy,  1)+  2*e2(ix,  iy+1,1)+&
                       &  4*e2(ix+1,iy+1,1)+  8*e2(ix+1,iy,  1)+&
                       & 16*e2(ix,  iy,  2)+ 32*e2(ix,  iy+1,2)+&
                       & 64*e2(ix+1,iy+1,2)+128*e2(ix+1,iy,  2)

                  !2011-05-26 Changed to int_coord type
                  ivslb(k)=int_coord(ix,iy,indx)
               end if

               !loop over those verticies in the triangle list for this 
               !index

            !loop to next box
            end do
         end do

         z2=real(2*iz)
         do i=1,k
            ix=ivslb(i)%i; iy=ivslb(i)%j; indx=ivslb(i)%k
            x2=real(2*ix); y2=real(2*iy)

            do i1=itind(indx),itind(indx+1)-1
               iv=iv+1 ;  n=itrn(i1)

               !for each vertex, if it has a number in the vertex slab 
               !use it
               !if not, increment the new vertex counter anf fill ivc 
               !with the right coordinaates
               select case(n)
               case(1)
                  j=ipv(3,ix,iy)

                  if (j.ne.0) then
                     !2011-05-26 Again added missing dimension to array 
                     !vindx (?)
                     vindx(iv)%i=j
                  else
                     ivn=ivn+1
                     ipv(3,ix,iy)=ivn
                     vindx(iv)%i=ivn
                     vert(ivn)=coord(x2,y2,z2+1)
                  end if
               case(2)
                  j=ipv(5,ix,iy)

                  if (j.ne.0) then
                     vindx(iv)%i=j
                  else
                     ivn=ivn+1
                     ipv(5,ix,iy)=ivn
                     vindx(iv)%i=ivn
                     vert(ivn)=coord(x2,y2+1,z2+1)
                  end if
               case(3)
                  j=ipv(3,ix,iy+1)

                  if (j.ne.0) then
                     vindx(iv)%i=j
                  else
                     ivn=ivn+1
                     ipv(3,ix,iy+1)=ivn
                     vindx(iv)%i=ivn
                     vert(ivn)=coord(x2,y2+2,z2+1)
                  end if
               case(4)
                  j=ipv(2,ix,iy)

                  if (j.ne.0) then
                     vindx(iv)%i=j
                  else
                     ivn=ivn+1
                     ipv(2,ix,iy)=ivn
                     vindx(iv)%i=ivn
                     vert(ivn)=coord(x2,y2+1,z2)
                  end if
               case(5)
                  j=ipv(1,ix,iy)

                  if (j.ne.0) then
                     vindx(iv)%i=j
                  else
                     ivn=ivn+1
                     ipv(1,ix,iy)=ivn
                     vindx(iv)%i=ivn
                     vert(ivn)=coord(x2+1,y2,z2)
                  end if
               case(6)
                  j=ipv(4,ix,iy)

                  if (j.ne.0) then
                     vindx(iv)%i=j
                  else
                     ivn=ivn+1
                     ipv(4,ix,iy)=ivn
                     vindx(iv)%i=ivn
                     vert(ivn)=coord(x2+1,y2,z2+2)
                  end if
               case(7)
                  j=ipv(4,ix,iy+1)

                  if (j.ne.0) then
                     vindx(iv)%i=j
                  else
                     ivn=ivn+1
                     ipv(4,ix,iy+1)=ivn
                     vindx(iv)%i=ivn
                     vert(ivn)=coord(x2+1,y2+2,z2+2)
                  end if
               case(8)
                  j=ipv(1,ix,iy+1)

                  if (j.ne.0) then
                     vindx(iv)%i=j
                  else
                     ivn=ivn+1
                     ipv(1,ix,iy+1)=ivn
                     vindx(iv)%i=ivn
                     vert(ivn)=coord(x2+1,y2+2,z2)
                  end if
               case(9)
                  j=ipv(3,ix+1,iy)

                  if (j.ne.0) then
                     vindx(iv)%i=j
                  else
                     ivn=ivn+1
                     ipv(3,ix+1,iy)=ivn
                     vindx(iv)%i=ivn
                     vert(ivn)=coord(x2+2,y2,z2+1)
                  end if
               case(10)
                  j=ipv(5,ix+1,iy)

                  if (j.ne.0) then
                     vindx(iv)%i=j
                  else
                     ivn=ivn+1
                     ipv(5,ix+1,iy)=ivn
                     vindx(iv)%i=ivn
                     vert(ivn)=coord(x2+2,y2+1,z2+2)
                  end if
               case(11)
                  j=ipv(3,ix+1,iy+1)
 
                  if (j.ne.0) then
                     vindx(iv)%i=j
                  else
                     ivn=ivn+1
                     ipv(3,ix+1,iy+1)=ivn
                     vindx(iv)%i=ivn
                     vert(ivn)=coord(x2+2,y2+2,z2+1)
                  end if
               case(12)
                  j=ipv(2,ix+1,iy)

                  if (j.ne.0) then
                     vindx(iv)%i=j
                  else
                     ivn=ivn+1
                     ipv(2,ix+1,iy)=ivn
                     vindx(iv)%i=ivn
                     vert(ivn)=coord(x2+2,y2+1,z2)
                  end if
               end select
            end do
         end do

         do i=1,k
            ix=ivslb(i)%i; iy=ivslb(i)%j
            ipv(3,ix,iy)=0
            ipv(1,ix,iy)=ipv(4,ix,iy)
            ipv(2,ix,iy)=ipv(5,ix,iy)
            ipv(4,ix,iy)=0
            ipv(5,ix,iy)=0
         end do
      end do

      !it loop for timing
      !snow= cputime(tnow)
      !write(6,*) "number of vertices= ",ivn
      !write(6,*) "number of triangles = ",iv/3

      end subroutine ex
