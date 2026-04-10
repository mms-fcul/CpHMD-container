!#####################################################################
!  ATTENTION!  This file is part of epsmakmod module.
!              Do not compile separately!
!#####################################################################

!---------------------------------------------------------------------
      subroutine cubedata(fac,cbln)

      implicit none

      !2011-05-27 Declarations added due to IMPLICIT NONE
      real :: fac,cbln,off
      type(coord) :: xyzp
      off=0.1

      !2011-05-27 Using operations on coord and int_coord type 
      !           variables defined in module operators_on_coordinates 
      xyzo=cmin-((fac*cbln)+off)
      xyzp=cmax+((fac*cbln)+off)
      !xo=-fac*cbln-off
      !yo=-fac*cbln-off
      !zo=-fac*cbln-off
      !xp=+fac*cbln+off
      !yp=+fac*cbln+off
      !zp=1.+fac*cbln+off
      
      lmncb=int((xyzp-xyzo)/cbln)
      lcb=lmncb%i; mcb=lmncb%j; ncb=lmncb%k
      !bl1=xp-xo
      !bl2=yp-yo
      !bl3=zp-zo
      !lcb=bl1/cbln
      !mcb=bl2/cbln
      !ncb=bl3/cbln
      
      cbai=1./cbln
      
      end subroutine cubedata

!---------------------------------------------------------------------
! 2011-05-27 Other parameters transfered via modules qlog and pointers
!            rda changed to delphipdb()%rad3          
      subroutine cube(crd,prbrd,cbn1,cbn2)
      
      implicit none

      type(coord) :: crd(natom)

      !here rda equals rad3
      !2011-05-27 Arrays declared in pointers and allocated in calling 
      !subroutine
      integer cbn1(0:lcb,0:mcb,0:ncb),cbn2(0:lcb,0:mcb,0:ncb)
      integer newatm,objecttype,itmp,kind

      !2011-05-27 Variables are accessible via qlog and pointers modules
      !           integer nobject,numbmol
      !           integer cbal(1) 

      !creating a set of fictious atoms occupying little cubes
      type(int_coord) :: icbn(natom+(nobject-numbmol)*(lcb+1)*(mcb+1)*(ncb+1))

      !iatmobj connects the fictious atoms to the objects
      integer iatmobj((nobject-numbmol)*(lcb+1)*(mcb+1)*(ncb+1))
      
      character(96) :: strtmp
      type(coord) :: vectx,vecty,vectz,dxyz,tmpvect,xa,xb,xyz2,tmpvect1
      type(coord) :: tmpvect2,xc,xd,xq,xloc,xyz
      type(int_coord) :: ixyz
      real :: modx,mody
      real :: tmp,tmp1,dist,shift
      real modul,mod2,tmp2
      real alpha,tan2,dot
      real cost

      !2011-05-27 Declarations added due to IMPLICIT NONE
      real :: prbrd,cbln
      integer :: i,ii,ix,iy,iz,jx,jy,jz,newatom,icum

      !2011-05-27 Array allocation is not necessary, arrays are used 
      !           only in this subroutine 
      cbln=1./cbai

      !2011-05-27 Changed to array operations
      cbn1=1; cbn2=0

      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            !2011-05-27 Using operations on coord and int_coord type 
            !           variables defined in module 
            !           operators_on_coordinates 
            !           xyzo variable is declared in pointers module, 
            !           cbai in qlog
            xyz=(crd(i)-xyzo)*cbai ; ixyz=int(xyz)
            if(ixyz.vorlt.1) write(6,*)'ix,iy,iz: ',ixyz
            if(ixyz.vorge.lmncb) write(6,*)'ix,iy,iz: ',ixyz

            ix=ixyz%i ; iy=ixyz%j ; iz=ixyz%k

            do jz=iz-1,iz+1
               do jy=iy-1,iy+1
                  do jx=ix-1,ix+1
                     cbn2(jx,jy,jz)=cbn2(jx,jy,jz)+1
                  end do
               end do
            end do

            icbn(i)=ixyz
         end if
      end do

      newatm=natom;  cost=cbln*.87+1./scale;  shift=cost+prbrd

      !icbn will contain also coord center of fictious atoms
      do ii=1,nobject
         strtmp=dataobject(ii,1)
         read(strtmp(16:18),*)kind

         if (strtmp(1:4).ne.'is a'.and.kind.ne.2) then
            itmp=ii+natom

            do iz=0,ncb
               do iy=0,mcb
                  do ix=0,lcb
                     ixyz=int_coord(ix,iy,iz)
                     xq=((float(ixyz)+0.5)*cbln)+xyzo

                     !2011-05-27 Using operations on coord and 
                     !           int_coord type variables defined
                     !           in module operators_on_coordinates 
                     if ((limobject(ii)%min.vandle.(xq+shift)).and.&
                         & (limobject(ii)%max.vandge.(xq-shift))) then
                        call distobj(xq,dist,dxyz,ii,prbrd,.true.)

                        if (dist.le.cost) then
                           newatm=newatm+1
                           cbn2(ix,iy,iz)=cbn2(ix,iy,iz)+1
                           icbn(newatm)=ixyz
                           iatmobj(newatm-natom)=itmp
                        end if
                     end if
                  end do
               end do
            end do
         end if
      end do

      icum=1
      do iz=0,ncb
         do iy=0,mcb
            do ix=0,lcb
               if (cbn2(ix,iy,iz).gt.0) then
                  cbn1(ix,iy,iz)=icum; icum=icum+cbn2(ix,iy,iz)
               end if
            end do
         end do
      end do

      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
         end if
      end do

      do i=natom+1,newatm
         ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
         cbal(cbn1(ix,iy,iz))=iatmobj(i-natom)
         cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
      end do

      !-1,0,0
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            ix=ix-1
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%i=ix
         end if
      end do

      !1,0,0
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            ix=ix+2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%i=ix
         end if
      end do

      !0,-1,0
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            ix=ix-1
            iy=iy-1
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%i=ix
            icbn(i)%j=iy
         end if
      end do

      !0,1,0
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            iy=iy+2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%j=iy
         end if
      end do

      !0,0,-1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            iy=iy-1
            iz=iz-1
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%j=iy
            icbn(i)%k=iz
         end if
      end do

      !0,0,1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            iz=iz+2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%k=iz
         end if
      end do

      !nn=2
      !1,0,1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            ix=ix+1
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%i=ix
         end if
      end do

      !-1,0,1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            ix=ix-2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%i=ix
         end if
      end do

      !0,1,1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            ix=ix+1
            iy=iy+1
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%i=ix
            icbn(i)%j=iy
         end if
      end do

      !0,-1,1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            iy=iy-2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%j=iy
         end if
      end do

      !-1,-1,0
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            ix=ix-1
            iz=iz-1
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%i=ix
            icbn(i)%k=iz
         end if
      end do

      !1,-1,0
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            ix=ix+2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%i=ix
         end if
      end do

      !1,1,0
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            iy=iy+2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%j=iy
         end if
      end do

      !-1,1,0
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            ix=ix-2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%i=ix
         end if
      end do

      !-1,0,-1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            iz=iz-1
            iy=iy-1
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%j=iy
            icbn(i)%k=iz
         end if
      end do

      !1,0,-1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            ix=ix+2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%i=ix
         end if
      end do

      !0,1,-1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            ix=ix-1
            iy=iy+1
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%i=ix
            icbn(i)%j=iy
         end if
      end do

      !0,-1,-1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            iy=iy-2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%j=iy
         end if
      end do

      !nn=3
      !-1,-1,-1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            ix=ix-1
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%i=ix
         end if
      end do

      !1,-1,-1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            ix=ix+2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%i=ix
         end if
      end do

      !1,1,-1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            iy=iy+2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%j=iy
         end if
      end do

      !-1,1,-1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            ix=ix-2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%i=ix
         end if
      end do

      !-1,1,1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            iz=iz+2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%k=iz
         end if
      end do

      !1,1,1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            ix=ix+2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%i=ix
         end if
      end do

      !1,-1,1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            iy=iy-2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%j=iy
         end if
      end do

      !-1,-1,1
      do i=1,natom
         if (delphipdb(i)%rad3.gt.0.0) then
            ix=icbn(i)%i; iy=icbn(i)%j; iz=icbn(i)%k
            ix=ix-2
            cbal(cbn1(ix,iy,iz))=i
            cbn1(ix,iy,iz)=cbn1(ix,iy,iz)+1
            icbn(i)%i=ix
         end if
      end do

      !reset cbn1
      icum=1
      do iz=0,ncb
         do iy=0,mcb
            do ix=0,lcb
               if (cbn2(ix,iy,iz).gt.0) then
                  cbn1(ix,iy,iz)=icum
                  icum=icum+cbn2(ix,iy,iz)
                  cbn2(ix,iy,iz)=icum-1
               end if
            end do
         end do
      end do
      icum=icum-1

      !2011-05-27 Deallocation of arrays is not necessary, they are
      !           used oonly in this subroutine and are deallocated
      !           authomatically

      end subroutine cube
