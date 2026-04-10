      !  ATTENTION!  This file is part of setcrgmod module.
      !==================================================================
      ! 2011-06-03 All parameters transfered via qlog and pointers modules
      subroutine setcrg
      
      !!c gchrg is the fractional charge in electron units assigned to
      !!c each grid point, gchrgp is the position of each such charge on the
      !!c grid.
      
      ! 2011-06-03 Arrays are declared in pointers module and allocated
      !!c b+++++++++++++++++++++++
      ! 2011-06-03 Variables declared in qlog module
      real :: re,deb
      integer :: ic1,epsdim,inearest(6),iext,ibgp,cont
      logical :: isum(3*igrid)
      !----------------------------------------------------------
      ! 2011-06-03 Declarations added due to IMPLICIT NONE
      integer :: i,j,k,frmnum,i1,i2,ico,ig,isgrid
      integer :: itmp,iv,itemp,ib,ix,iy,iz,kb1,kb2,kb3
      integer :: n,iw,ios
      real :: difeps,cg1,cg2,cg3,fpoh,debfct,sixeps,temp
      ! 2011-06-03 Allocation of cgbp array will be done later after
      !!c e++++++++++++++++++++++
      
      !!c assign fractional charges to grid points, using phimap
      !!c as a dummy array (set to zero again later)
      difeps=epsin-epsout;  epsdim=natom+nobject+2
      sixeps=epsout*6.0;   fpoh=fpi*scale
      debfct = epsout/(deblen*scale)**2
      ! debug tested         write(6,*)'---SETRC:',difeps,epsdim,sixeps,fpoh,debfct
      
      do ig=1,nqgrd 
         kb1=chrgv2(ig)%xyz%x; kb2=chrgv2(ig)%xyz%y; kb3=chrgv2(ig)%xyz%z
         do ix=0,1
            i=kb1+ix
            cg1=kb1-chrgv2(ig)%xyz%x+1-ix
            do iy=0,1
               j=kb2+iy
               cg2=kb2-chrgv2(ig)%xyz%y+1-iy
               do iz=0,1
                  k=kb3+iz
                  cg3=kb3-chrgv2(ig)%xyz%z+1-iz
                  re=abs(cg1*cg2*cg3)*chrgv2(ig)%value
                  phimap(i,j,k)=phimap(i,j,k)+re
               end do
            end do
         end do
      end do
      !debug tested         write(6,*)'---SETRC:',re, cg1,cg2,cg3,ig, size(chrgv2)
      
      n=0
      do k=2,igrid-1
         do j=2,igrid-1
            do i=2,igrid-1
               if(phimap(i,j,k).ne.0.) n=n+1
            end do
         end do
      end do
      
      allocate(gchrgp(n),gchrg(n),gchrgd(n),gchrg2(n))
      allocate(qval(n),iqpos(n),gval(n))
      !debug tested         write(6,*)'---SETRC:',n
      
      !!c set up odd/even logical array 
      ! 2011-06-03 Changed to array operation
      isum=.false. ; isum(1:3*igrid:2)=.true.
      !!c find which grid points have charge assigned to them
      !!c (will use this array later to calculate grid energy)
      
      n=0
      do k=2,igrid-1
         do j=2,igrid-1
            do i=2,igrid-1
               if(phimap(i,j,k).ne.0.) then
                  n=n+1
                  gchrgp(n)=int_coord(i,j,k)
                  gchrg(n)=phimap(i,j,k)
                  phimap(i,j,k)=0.0
               end if
            end do
         end do
      end do
      icount1b=n
      ! 2011-06-03 Parameter ngcrg is no longer needed
      
      !!c b++++++++++++++++++++June 2001++to save memory waste time++++++
      do ig=1,nqgrd
         kb1=chrgv2(ig)%xyz%x; kb2=chrgv2(ig)%xyz%y; kb3=chrgv2(ig)%xyz%z
         ic1=nqgrdtonqass(ig)
         do ix=0,1
            i=kb1+ix
            cg1=kb1-chrgv2(ig)%xyz%x+1-ix
            do iy=0,1
               j=kb2+iy
               cg2=kb2-chrgv2(ig)%xyz%y+1-iy
               do iz=0,1
                  k=kb3+iz
                  cg3=kb3-chrgv2(ig)%xyz%z+1-iz
                  re=abs(cg1*cg2*cg3)*chrgv2(ig)%value
                  deb=0.
                  if (idebmap(i,j,k)) deb=1.
                  phimap(i,j,k)=phimap(i,j,k)+re/(6.*atmeps(ic1)+&
                  &   debfct*deb)
               end do
            end do
         end do
      end do
      allocate(gchrgtmp(icount1b))
      do n=1,icount1b
         i=gchrgp(n)%i; j=gchrgp(n)%j;  k=gchrgp(n)%k
         gchrgtmp(n)=phimap(i,j,k)
      end do
      
      ! 2011-06-03 Changed to array operation
      phimap(2:igrid-1,2:igrid-1,2:igrid-1)=0.
      
      !!c e+++++++++++++++++++++++++++++++++++++++++++++
      if(iwgcrg) then
         ! 2011-06-03 Subroutine wrtcrg is short and called only once, thus
         !---------------------------- Begin of wrtcrg.f -----------------------------
         open(20,file=trim(gcrgnam))
         frmnum=1
         write(20,*) "DELPHI OUTPUT FILE: GRID CHARGE"
         write(20,*) "FORMAT NUMBER=",frmnum
         write(20,*) "NUMBER OF CHARGES=",icount1b
         !!c b+++++++++++++++
         do i = 0,nmedia
            write(20,*)'DIELECTRIC IN MEDIUM NUMBER ',i,' :',medeps(i)*epkt
         end do
         !!c e+++++++++++++++                       
         write(20,*) "GRID SCALE=",scale
         do i=1,icount1b
            write(20,*)gchrg(i),gchrgp(i)
         end do
         close(20)
         !---------------------------- End of wrtcrg.f -----------------------------
      end if
      !!c determine how many charged grid points are odd
      
      icount1a=0
      do i=1,icount1b
         !  2011-06-03  Using operations on coord and int_coord type variables defined
         itemp=sum(gchrgp(i))
         if(isum(itemp)) icount1a=icount1a+1
      end do
      
      !!c set up odd/even pointer array, to be used in making qval
      !!c and iqpos
      
      i1=0 ; i2=icount1a
      do  i=1,icount1b
         itemp=sum(gchrgp(i))
         if(isum(itemp)) then
            i1=i1+1;  gchrg2(i)=i1
         else
            i2=i2+1; gchrg2(i)=i2
         end if
      end do
      
      !!c determine denominator at all charged grid points
      !---------------------------------------------------------
      ! 2011-06-03 Finding actual size for the array cgbp
      ib=0
      do i=1,icount1b
         ix=gchrgp(i)%i; iy=gchrgp(i)%j; iz=gchrgp(i)%k
         iext=0;  ibgp=0
         inearest(1)=iepsmp(ix,iy,iz)%i/epsdim
         inearest(2)=iepsmp(ix,iy,iz)%j/epsdim
         inearest(3)=iepsmp(ix,iy,iz)%k/epsdim
         inearest(4)=iepsmp(ix-1,iy,iz)%i/epsdim
         inearest(5)=iepsmp(ix,iy-1,iz)%j/epsdim
         inearest(6)=iepsmp(ix,iy,iz-1)%k/epsdim
         if(inearest(1).eq.0) iext=1
         if(inearest(1).ne.inearest(6)) ibgp=1
         do cont=2,6
            if(inearest(cont).eq.0) iext=1
            if(inearest(cont).ne.inearest(cont-1)) ibgp=1
         end do
         if ((ibgp.eq.1).or.(iext.eq.1)) ib=ib+1
      end do
      !---------------------------------------------------------
      ibc=ib
      allocate(cgbp(ibc))
      ib=0;  ico=0
      do i=1,icount1b
         !!c     let us scan all the charged grid points
         ix=gchrgp(i)%i; iy=gchrgp(i)%j; iz=gchrgp(i)%k
         
         !!c b++++++++++++++++++++
         iext=0;  ibgp=0
         inearest(1)=iepsmp(ix,iy,iz)%i/epsdim
         inearest(2)=iepsmp(ix,iy,iz)%j/epsdim
         inearest(3)=iepsmp(ix,iy,iz)%k/epsdim
         inearest(4)=iepsmp(ix-1,iy,iz)%i/epsdim
         inearest(5)=iepsmp(ix,iy-1,iz)%j/epsdim
         inearest(6)=iepsmp(ix,iy,iz-1)%k/epsdim
         if(inearest(1).eq.0) iext=1
         if(inearest(1).ne.inearest(6)) ibgp=1
         do cont=2,6
            if(inearest(cont).eq.0) iext=1
            if(inearest(cont).ne.inearest(cont-1)) ibgp=1
         end do
         deb=0.
         if (idebmap(ix,iy,iz)) deb=1.
         if (idirectalg.eq.0)  then
            !!c           to be updated
            gchrgd(i)=itemp*difeps + debfct*deb + sixeps 
            if(nmedia.gt.1) then
               write(6,*)'in setcrg this part needs to be changed'
               stop
            end if
         else
            if (gaussian.eq.0) then
            temp=0.0
            itmp=iepsmp(ix,iy,iz)%i/epsdim
            temp=temp+medeps(itmp)
            itmp=iepsmp(ix,iy,iz)%j/epsdim
            temp=temp+medeps(itmp)
            itmp=iepsmp(ix,iy,iz)%k/epsdim
            temp=temp+medeps(itmp)
            itmp=iepsmp(ix-1,iy,iz)%i/epsdim
            temp=temp+medeps(itmp)
            itmp=iepsmp(ix,iy-1,iz)%j/epsdim
            temp=temp+medeps(itmp)
            itmp=iepsmp(ix,iy,iz-1)%k/epsdim
            temp=temp+medeps(itmp)
            !!c here temp=sum(eps of 6 closest midpoints)
            gchrgd(i)=temp+ debfct*deb
!            print *,'LinLi in setcrg:temp,debfct,deb',temp,debfct,deb

            else if(gaussian.eq.1) then
               if(uniformdiel)then
                  temp=epsin*6
                  gchrgd(i)=temp+debfct*deb
               else
!                 temp=0.0
                  temp=gepsmp2(ix,iy,iz)%x+gepsmp2(ix,iy,iz)%y+gepsmp2(ix,iy,iz)%z &
                  & +gepsmp2(ix-1,iy,iz)%x+gepsmp2(ix,iy-1,iz)%y+gepsmp2(ix,iy,iz-1)%z
!                print *,'Lin Li,temp,gepsmp2(ix,iy,iz)',temp,gepsmp2(ix,iy,iz)
                  temp=temp/epkt
                  gchrgd(i)=temp+ debfct*deb
!                print *,'Lin Li,i,gchrgd(i),temp,epkt',i,gchrgd(i),temp,epkt
!                print *,'Lin Li,gepsmp2(ix,iy,iz),',gepsmp2(ix,iy,iz)
               endif
            else
            endif


         end if
         !!c e++++++++++++++++++++
         if ((ibgp.eq.1).or.(iext.eq.1)) then
            ib=ib+1
            if(ibgp.eq.0) ico=ico+1
            cgbp(ib)%xyz=coord(real(ix),real(iy),real(iz))
            !!c b++++++++++++++++++++
            cgbp(ib)%value1=gchrgtmp(i)*fpoh
            !!c e++++++++++++++++++++
            cgbp(ib)%value2=gchrg2(i)
         end if
      end do
      if(allocated(gchrgtmp)) deallocate(gchrgtmp)
      if(verbose)&
      &  write(6,*) 'no. charged boundary grid points =',ibc
      ! 2011-06-03 Parameter ibcmax is no longer needed
      if(ico.ne.0) write(6,*) "## out of them,",ico," charged points &
      &are in solution ## Increased resolution is needed"
      
      !!c make qval, fpoh term so potentials will be in kt/e
      
      do i=1,icount1b
         j=gchrg2(i)
         qval(j)=gchrg(i)*(fpoh/gchrgd(i))
         gval(j)=gchrg(i)
      end do
      if(allocated(gchrgd)) deallocate(gchrgd) 
      
      !!c make iqpos
      isgrid=igrid*igrid
      do i=1,icount1b
         j=gchrg2(i)
         ix=gchrgp(i)%i; iy=gchrgp(i)%j; iz=gchrgp(i)%k
         iw=1+ix+igrid*(iy-1)+isgrid*(iz-1)
         iv=iw/2; iqpos(j)=iv
      end do
      if(allocated(gchrg2)) deallocate(gchrg2,stat=ios)
      
      !!c end of chrgup, return with qval,iqpos and gchrgp and gchrg
      !!c also icount1a, icount1b 
      ! 2011-06-03 Did not understand meaning of this statement
      end subroutine setcrg
