!#####################################################################
!ATTENTION!  This file is part of misc2 module.
!            Do not compile separately!
!#####################################################################

      !==================== From file phintp4.f ======================
      !interpolates the potential at any point inside a cubical volume 
      !using the potential values at the 8 vertices by means of a 
      !trilinear function:
	   !
	   !W = A1.XYZ + A2.XY + A3.XZ + A4.YZ + A5.X +A6.Y + A7.Z + A8
	   ! 
      !where Ai coefficients are linear combinations of Wi at the cube 
      !corners
      subroutine phintp(gp,phi)

      !include 'qdiffpar4.h'
      !include 'qlog.h'

      !dimension gp(3),phimap(igrid,igrid,igrid)
      type(coord) :: gp

      !2011-06-06 Declarations added due to IMPLICIT NONE
      integer :: nx,ny,nz,nx1,ny1,nz1
      real :: phi,rgrid,a1,a2,a3,a4,a5,a6,a7,a8
      real:: xg,yg,zg,xgr,ygr,zgr

      !return 0.0 if outside grid

      rgrid = real(igrid)
      if ((gp.vorlt.1.).or.(gp.vorgt.rgrid) ) then
         !do 9000 k = 1,3
         !   if((gp(k).lt.1).or.(gp(k).gt.rgrid)) then
         phi = 0.0
         write(6,*)"Pay attention, point out of the cube!!"
         write(6,*)"Values:",gp," Igrid:",igrid
         !write(6,*)"i=",k,"Value:",gp(k),"Igrid:",igrid
         return
      end if
      !end if
      !9000        continue
      !end do

      !find lower left bottom grid point
      !nx = int(gp(1))
      nx = int(gp%x); nx1 = nx+1; if(nx1.gt.igrid) nx1=nx
      !ny = int(gp(2))
      ny = int(gp%y); ny1 = ny+1; if(ny1.gt.igrid)ny1=ny
      !nz = int(gp(3))
      nz = int(gp%z); nz1 = nz+1; if(nz1.gt.igrid)nz1=nz

      !calculate cube coordinates of point
      xg = nx; yg = ny; zg = nz
      xgr=gp%x-xg; ygr=gp%y-yg; zgr=gp%z-zg
      !xgr=gp(1)-xg; ygr=gp(2)-yg; zgr=gp(3)-zg
         
      !calculate coefficients of trilinear function
      a8 = phimap(nx,ny,nz)
      a7 = phimap(nx,ny,nz1) - a8
      a6 = phimap(nx,ny1,nz) - a8
      a5 = phimap(nx1,ny,nz) - a8
      a4 = phimap(nx,ny1,nz1) - a8-a7-a6
      a3 = phimap(nx1,ny,nz1) - a8-a7-a5
      a2 = phimap(nx1,ny1,nz) - a8-a6-a5
      a1 = phimap(nx1,ny1,nz1) -a8-a7-a6-a5-a4-a3-a2
            
      !determine value of phi
      phi = a1*xgr*ygr*zgr + a2*xgr*ygr + a3*xgr*zgr &
            &+ a4*ygr*zgr + a5*xgr + a6*ygr + a7*zgr + a8

      end subroutine phintp
