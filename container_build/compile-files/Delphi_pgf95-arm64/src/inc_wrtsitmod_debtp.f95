      !  ATTENTION!  This file is part of wrtsitmod module.
      !==================================================================
      subroutine debtp(gp,deb)
      
      !!c interpolates the debye map at any point inside
      !!c a cubical volume using idebmap values at the
      !!c 8 vertices by means of a trilinear function:
      	
      !!C W = A1.XYZ + A2.XY + A3.XZ + A4.YZ + A5.X +
      !!C		A6.Y + A7.Z + A8
      	
      !!C where Ai coefficients are linear combinations of Wi at
      !!C the cube corners
      
      !!c---------------------------------------------------
      type(coord) :: gp,xyzgr
      type(int_coord) :: nxyz, nxyz1
      real :: debm(8)
      !---------------------------------------------------------------------------
      ! 2011-06-12 Declarations added due to IMPLICIT NONE
      real :: rgrid,xgr,ygr,zgr,a1,a2,a3,a4,a5,a6,a7,a8,deb
      integer :: k,nx,ny,nz,nx1,ny1,nz1
      !---------------------------------------------------------------------------
      !!c return 0.0 if outside grid
      rgrid = real(igrid)
      if((gp.vorlt.1.).or.(gp.vorgt.rgrid)) then
         deb = 0.0
         write(6,*)"Pay attention, point out of the cube!!"
         write(6,*)"Values:",gp,"Igrid:",igrid
         return
      end if
      
      !!c find lower left bottom grid point
      nxyz=int(gp);nxyz1=nxyz+1
      if(nxyz1%i.gt.igrid)nxyz1%i=nxyz%i
      if(nxyz1%j.gt.igrid)nxyz1%j=nxyz%j
      if(nxyz1%k.gt.igrid)nxyz1%k=nxyz%k
      
      !!c calculate cube coordinates of point
      xyzgr=gp-float(nxyz)
      
      !!c convert logical idebmap to real
      debm=0.
      nx=nxyz%i; ny=nxyz%j; nz=nxyz%k
      nx1=nxyz1%i; ny1=nxyz1%j; nz1=nxyz1%k
      if (idebmap(nx,ny,nz)) debm(1)=1.0
      if (idebmap(nx,ny,nz1)) debm(2)=1.0
      if (idebmap(nx,ny1,nz)) debm(3)=1.0
      if (idebmap(nx1,ny,nz)) debm(4)=1.0
      if (idebmap(nx,ny1,nz1)) debm(5)=1.0
      if (idebmap(nx1,ny,nz1)) debm(6)=1.0
      if (idebmap(nx1,ny1,nz)) debm(7)=1.0
      if (idebmap(nx1,ny1,nz1)) debm(8)=1.0
      
      !!c calculate coefficients of trilinear function
      a8 = debm(1)
      a7 = debm(2) - a8
      a6 = debm(3) - a8
      a5 = debm(4) - a8
      a4 = debm(5) - a8 - a7 - a6
      a3 = debm(6) - a8 - a7 - a5
      a2 = debm(7) - a8 - a6 - a5
      a1 = debm(8) - a8 - a7 - a6 - a5 - a4 - a3 - a2
      
      !!c determine value of phi
      xgr=xyzgr%x; ygr=xyzgr%y; zgr=xyzgr%z
      deb = a1*xgr*ygr*zgr + a2*xgr*ygr + a3*xgr*zgr + &
      &a4*ygr*zgr + a5*xgr + a6*ygr + a7*zgr + a8
      
      end subroutine debtp