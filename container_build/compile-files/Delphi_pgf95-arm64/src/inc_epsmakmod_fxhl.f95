!#####################################################################
!ATTENTION!  This file is part of epsmakmod module.
!            Do not compile separately!
!
! 2011-05-25 Arrays to subroutine are transfered via pointers module
!#####################################################################

      subroutine fxhl(ivbeg,ivend,itbeg,itend,ntot,itot,imxvtx,imxtri)

      implicit none

      integer,parameter :: mxvtx=150000, mxtri=300000
      integer :: intsiz,vtot,itot
      integer itbeg,itend,ivbeg,ivend,ntot,imxvtx,imxtri
      integer is,ir,i,j1,j2,j3,i1,i2,n,it,id,nsel(mxvtx),tv(20),vl(100)
      real a1,b1,c1,a2,b2,c2,disnor,xnor,ynor,znor,p1,p2,p3
      integer j,k,n1
      type(int_coord) :: j123,k123
      type(coord) :: abc1, abc2, xyznor

      ntot=itot; nsel=0

      is=0; ir=0
      do i=ivbeg,ivend
         i1=vtpnt(i); i2=i1+vtlen(i)-1

         n=0
         do j=i1,i2
            j123=vindx(vtlst(j))

            if (j123%i.gt.i) then
               nsel(j123%i)=nsel(j123%i)+1; n=n+1; vl(n)=j123%i
            end if

            if (j123%j.gt.i) then
               nsel(j123%j)=nsel(j123%j)+1; n=n+1; vl(n)=j123%j
            end if

            if (j123%k.gt.i) then
               nsel(j123%k)=nsel(j123%k)+1; n=n+1; vl(n)=j123%k
            end if

         !k= one of the triangles
         end do

         it=0
         do j=1,n
            k=vl(j)
          
            if ((nsel(k).ne.0).and.(nsel(k).ne.2)) then
               it=it+1; tv(it)=k
            end if

            nsel(k)=0
         end do

         if(it.gt.0) is=is+1

         !mend
         if (it.eq.2) then
            n1=tv(1);  id=1
            do j=i1,i2
               j123=vindx(vtlst(j))
               if(j123%i.eq.n1) id=1
               if(j123%k.eq.n1) id=2

               if (j123%j.eq.n1) then
                  if(j123%i.eq.i) id=2
               else
                  if(j123%k.eq.i) id=1
               end if
            end do
                 
            j1=i
            j2=tv(1)
            j3=tv(2)

            !2011-05-25  Using operations on coord and int_coord type 
            !variables defined in module operators_on_coordinates 

            abc1=vert(j2)-vert(j1)
            abc2= vert(j3)-vert(j1)
            xyznor=abc1.x.abc2
            
            disnor=sqrt(xyznor.dot.xyznor)
            xyznor=xyznor/disnor
            
            p1=xyznor.dot.vnorm(j1)
            p2=xyznor.dot.vnorm(j2)
            p3=xyznor.dot.vnorm(j3)

            !2011-05-25 Just compacted the code
            itot=itot+1

            if ((p1.lt.0).and.(p2.lt.0).and.(p3.lt.0)) then
               vindx(itot)=int_coord(tv(2),tv(1),i)
            else
               vindx(itot)=int_coord(i,tv(1),tv(2))
            end if

            ir=ir+1
         end if
      end do

      ntot=itot-ntot

      !remake vertex to triangles listing..
      call mkvtl(ivbeg,ivend,itbeg,itot,imxtri,imxvtx)
           
      end subroutine fxhl
