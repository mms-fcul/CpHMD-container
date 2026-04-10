!#####################################################################
!ATTENTION! This file is part of epsmakmod module.
!           Do not compile separately!
!
!2011-05-25 Arrays to subroutine are transfered via pointers module
!#####################################################################
      subroutine mkvtl(ivbeg,ivend,itbeg,itend,imxtri,imxvtx)

      implicit none

      integer mxvtx, mxtri, intsiz, itbeg, itend, ivbeg, ivend
      integer imxvtx, imxtri
      parameter(mxtri=300000)
      integer vtemp(mxtri), i, j1, j2, j3, k1, k2, k3
      integer im, it1, it2, j
      real tarray1(2),tarray2(2), t1
      type(int_coord) :: j123, k123

      !2011-05-25 Changed to array operations
      vtlen=0

      do i=itbeg,itend
         !2011-05-25 Changed to int_coord type 
         j123=vindx(i)
         vtlen(j123%i)=vtlen(j123%i)+1
         vtlen(j123%j)=vtlen(j123%j)+1
         vtlen(j123%k)=vtlen(j123%k)+1
      end do

      if(ivbeg.ne.1) vtpnt(ivbeg)=vtpnt(ivbeg-1)+vtlen(ivbeg-1)
      if(ivbeg.eq.1) vtpnt(ivbeg)=1

      do i=ivbeg+1,ivend
         vtpnt(i)=vtpnt(i-1)+vtlen(i-1)
      end do
 
      vtpnt(ivend+1)=vtpnt(ivend)+vtlen(ivend)
           
      vtlen=0

      do i=itbeg,itend
         j123=vindx(i)
         k1=vtpnt(j123%i)+vtlen(j123%i)
         k2=vtpnt(j123%j)+vtlen(j123%j)
         k3=vtpnt(j123%k)+vtlen(j123%k)
         vtlen(j123%i)=vtlen(j123%i)+1
         vtlen(j123%j)=vtlen(j123%j)+1
         vtlen(j123%k)=vtlen(j123%k)+1
         vtlst(k1)=i; vtlst(k2)=i; vtlst(k3)=i
      end do

      vtemp=0

      im=0
      do i=itbeg,itend
         j123=vindx(i)

         !which triangle borders edge j1-j2
         do j=vtpnt(j123%i),vtpnt(j123%i)+vtlen(j123%i)-1
            k1=vtlst(j)
            vtemp(k1)=j123%i
         end do

         it1=0
         do j=vtpnt(j123%j),vtpnt(j123%j)+vtlen(j123%j)-1
            k1=vtlst(j)
            if((vtemp(k1).eq.j1).and.(k1.ne.i)) it1=k1
            vtemp(k1)=j123%j
         end do

         !which point completes the triangle
         it2=0
         if (it1.ne.0) then
            k123=vindx(it1)
            if((k123%i.ne.j123%i).and.(k123%i.ne.j123%j)) it2=k123%i
            if((k123%j.ne.j123%i).and.(k123%j.ne.j123%j)) it2=k123%j
            if((k123%k.ne.j123%i).and.(k123%k.ne.j123%j)) it2=k123%k
         end if

         !fill tmlst
         tmlst(1,i)=it1; tmlst(4,i)=it2; tmlst(7,i)=j123%k
         im=im+it1
             
         !which triangle borders edge j2-j3
         it1=0
         do j=vtpnt(j123%k),vtpnt(j123%k)+vtlen(j123%k)-1
            k1=vtlst(j)
            if((vtemp(k1).eq.j123%j).and.(k1.ne.i)) it1=k1
            vtemp(k1)=j123%k
         end do

         !which point completes the triangle
         it2=0
         if (it1.ne.0) then
            k123=vindx(it1)
            if((k123%i.ne.j123%j).and.(k123%i.ne.j123%k)) it2=k123%i
            if((k123%j.ne.j123%j).and.(k123%j.ne.j123%k)) it2=k123%j
            if((k123%k.ne.j123%j).and.(k123%k.ne.j123%k)) it2=k123%k
         end if
  
         !fill tmlst
         tmlst(2,i)=it1; tmlst(5,i)=it2; tmlst(8,i)=j123%i
         im=im+it1

         !which triangle borders edge j3-j1
         it1=0
         do j=vtpnt(j123%i),vtpnt(j123%i)+vtlen(j123%i)-1
            k1=vtlst(j)
            if((vtemp(k1).eq.j3).and.(k1.ne.i)) it1=k1
         end do

         !which point completes the triangle
         it2=0
         if (it1.ne.0) then
            k123=vindx(it1)
            if((k123%i.ne.j123%k).and.(k123%i.ne.j123%i)) it2=k123%i
            if((k123%j.ne.j123%k).and.(k123%j.ne.j123%i)) it2=k123%j
            if((k123%k.ne.j123%k).and.(k123%k.ne.j123%i)) it2=k123%k
         end if

         !fill tmlst
         tmlst(3,i)=it1; tmlst(6,i)=it2;  tmlst(9,i)=j123%j
         im=im+it1

      end do

      end subroutine mkvtl
