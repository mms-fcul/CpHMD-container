!#####################################################################
!ATTENTION! This file is part of epsmakmod module.
!           Do not compile separately!
!
!2011-05-20 Parameters to subroutine are transfered via qlog and 
!           pointers modules
!#####################################################################
      subroutine msrf
      
      implicit none

      !2011-05-20 Leaving variables not declared in qlog and pointers 
      !modules
      integer :: itot, vtot
      character(5) :: fixyn
      character(80) :: vdat1, line(5), fname
      type(coord) :: xo,xo2
      logical :: fix
      integer :: epsdim
      type(coord) :: v1,v2,v3,vxyz
      logical :: out,nbe(0:6),exists
      type(int_coord) :: iv123
      !2011-05-27 Declarations added due to IMPLICIT NONE
      integer :: mxvtx,i,ia1,ia2,ia3,iate,imxtri,it,imxvtx
      integer :: ib,iv1,iv2,iv3,j,k,mxtri,ntot,ntot2
      real :: aa,arear,areac,areas,area,cc,bb,rad,rad2,ss
      real :: tar,tne4,vmg

      mxvtx=ibnum*2.2; mxtri=2*mxvtx; epsdim=natom+nobject+2

      !2011-05-20 Using operations on coord and int_coord type 
      !variables defined in module operators_on_coordinates 
      do k=1,igrid
         do j=1,igrid
            do i=1,igrid
               !cambiato da mod a div, mi dovrebbe servire solo il mezzo
               egrid(i,j,k)=iepsmp(i,j,k)/epsdim
            end do
         end do
      end do

      do i=1,igrid
         do j=1,igrid
            egrid(i,1,j)%i = egrid(i,1,j)%j
            egrid(i,j,1)%i = egrid(i,j,1)%k
         end do
      end do

      do k=2,igrid
         do j=2,igrid
            do i=2,igrid
               iate = 0
               if (egrid(i,j,k)%i .gt. 0) iate = iate + 1
               if (egrid(i,j,k)%j .gt. 0) iate = iate + 1
               if (egrid(i,j,k)%k .gt. 0) iate = iate + 1
               if (egrid(i-1,j,k)%i .gt. 0) iate = iate + 1
               if (egrid(i,j-1,k)%j .gt. 0) iate = iate + 1
               if (egrid(i,j,k-1)%k .gt. 0) iate = iate + 1

               if (iate.le.3) then
                  egrid(i,j,k)%i = 0
               else
                  egrid(i,j,k)%i = 1
               end if
            end do
         end do
      end do
           
      !2011-05-20 Arrays are allocated with ALLOCATE F95 statement
      !           Initial values are added (just in case)
      allocate(vindx(mxtri),vert(mxvtx))
      vindx=int_coord(0,0,0); vert=coord(0,0,0)
      print*,'mxtri,mxvtx=',mxtri,mxvtx

      !***************************************************************
      !2011-05-20  In subroutine ex egrid array is three-dimensional 
      !            while in this subroutine this array was four-
      !            dimensional.
      !            Also array vindx was 2D here, but 1D in the EX 
      !            subroutine.
      !            Also it should be clarified why file in /usr/local/
      !            bin is opened and read in ex subroutine.
      vdat1='./'
      call ex(vtot, itot, vdat1, 2 )

      if (vtot.gt.mxvtx) then
         write(6,*)'vtot = ',vtot,' > mxvtx = ',mxvtx
         write(6,*)'increase mxvtx in msrf.f'
         stop
      end if

      do ib=1,vtot
         vert(ib)=vert(ib)/2.
      end do

      itot = itot/3

      !scale boundary grid point positions relative to acc data
      write(6,*)'scaling vertices'
      allocate(vnorm(vtot),vnorm2(vtot))

      call sclbp(vert,vnorm,vtot,iab1,iab2)

      !fix holes and make vertex to triangle arrays
      !allocate hole-fixing arrays next
      !hole-fixing variables
      if (vtot .lt. mxvtx/2) then
         imxvtx = vtot*2
      else
         imxvtx = mxvtx
      end if

      if (itot .lt. mxtri/2) then
         imxtri = itot*2
      else
         imxtri = mxtri
      end if
      
      allocate(vtlen(imxvtx),vtlst(6*imxvtx))
      allocate(vtpnt(imxvtx),tmlst(9,imxvtx))

      !2011-05-25 Arrays to subroutine are transfered via pointers 
      !           module
      call mkvtl(1,vtot, 1,itot, imxtri, imxvtx)

      call fxhl(1,vtot,1,itot,ntot, itot, imxvtx, imxtri)

      call fxhl(1,vtot,1,itot,ntot2, itot, imxvtx, imxtri)

      !write(*,*) "number of triangles added to fix holes= ",ntot+ntot2
 
      if (ntot2 .gt. 0)  then
         fix = .true.
      else
         fix = .false.
      end if

      do while (fix)
         call fxhl(1,vtot,1,itot,ntot2, itot, imxvtx, imxtri)

         !write(*,*) "number of triangles added to fix holes= ",ntot2
         
         if (ntot2 .gt. 0) then
            fix = .true.
         else
            fix = .false.
         end if
      end do
      
      if (itot.gt.mxtri) then
         write(6,*)'itot = ',itot,' > mxtri = ',mxtri
         write(6,*)'increase mxtri in msrf.f'
         stop
      end if
      
      if(allocated(vtlen)) deallocate(vtlen)
      if(allocated(vtlst)) deallocate(vtlst)
      if(allocated(tmlst)) deallocate(tmlst)
      if(allocated(vtpnt)) deallocate(vtpnt)

      write(6,*) 'number of vertices = ',vtot
      write(6,*) 'number of triangles = ',itot

      vnorm2=coord(0.,0.,0)

      !calculate area
      area=0.0; areas=0.0; areac=0.0; arear=0.0

      !2011-05-20  Using operations on coord and int_coord type 
      !variables defined in module operators_on_coordinates 
      do it=1,itot
         iv123=vindx(it)
         v1=vert(iv2)-vert(iv1); v2=vert(iv3)-vert(iv1)

         !2011-05-20 Vector product defined in operators_on_coordinates 
         !           module
         vxyz=v1.x.v2; vmg=sqrt(vxyz.dot.vxyz)
         tar=vmg/2.
         vxyz=vnorm(iv1)+vnorm(iv2)+vnorm(iv3); vmg=sqrt(vxyz.dot.vxyz)

         vnorm2(iv1)=vnorm2(iv1)+(vxyz/vmg)
         vnorm2(iv2)=vnorm2(iv2)+(vxyz/vmg)
         vnorm2(iv3)=vnorm2(iv3)+(vxyz/vmg)

         !calculate spherical triangle area if appropriate
         ia1=atndx(iv1); ia2=atndx(iv2); ia3=atndx(iv3)

         if (ia1.gt.0) then
            if (ia1.eq.ia2.and.ia1.eq.ia3) then
               rad=delphipdb(ia1)%rad3;  rad2=rad*rad
               aa=sum((vert(iv2)-vert(iv1))*(vert(iv2)-vert(iv1)))
               bb=sum((vert(iv3)-vert(iv2))*(vert(iv3)-vert(iv2)))
               cc=sum((vert(iv1)-vert(iv3))*(vert(iv1)-vert(iv3)))

               aa=acos(1.-aa/(2.*rad2))
               bb=acos(1.-bb/(2.*rad2))
               cc=acos(1.-cc/(2.*rad2))
               ss=(aa+bb+cc)*.5
               tne4=sqrt(tan(ss*.5)*tan((ss-aa)*.5)*&
                                     &tan((ss-bb)*.5)*tan((ss-cc)*.5))
               tar=4.*atan(tne4)*rad2
            end if
         end if
         
         area=area+tar
      
      end do
	
      do i=1,vtot
         vmg=sqrt(vnorm2(i).dot.vnorm2(i))
         vnorm2(i)=vnorm2(i)/vmg
      end do

      write(6,*)'MS area                = ',area

      !2011-05-26 Subroutine wrtsurf is short and called only once,
      !           therefore put the body of the subroutine here
      !--------  Begin of  wrtsurf body -----------------------------
      if (.not.ibem) then
         fname='grasp.srf'
         write(6,*) "writing GRASP file to ",trim(fname)
         line=' '
         line(1)='format=2'
         line(2)='vertices,normals,triangles'
         write(line(4),'(3i6,f12.6)') vtot,itot,igrid,scale
         write(line(5),'(3f12.6)') oldmid
         
         !NB, from 6/14/92 output surface coordinates in angstroms..
         open(20,file=trim(fname),form="unformatted")
         do i=1,5
            write(20) line(i)
         end do
         write(6,*) 'writing data for',vtot, ' vertices and', &
                     &	itot, ' triangles'
         write(20) vert ; write(20) vnorm ; write(20) vindx
         close(20)
         write(6,*) 'finished writing ',trim(fname)
      else
         open(7,file="bem.srf")
         write(7,*)vtot,itot
      
         do i=1,vtot
            write(7,*)vert(i)
         end do
      
         do i=1,itot
            write(7,*)vindx(i)
         end do
      
         do i=1,vtot
            write(7,*)vnorm(i)
         end do
      
         close (7)
      end if
      !-------- End of  wrtsurf body --------------------------------

      end subroutine msrf
