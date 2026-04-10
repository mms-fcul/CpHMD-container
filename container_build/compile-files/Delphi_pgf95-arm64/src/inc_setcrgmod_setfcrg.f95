      !  ATTENTION!  This file is part of setcrgmod module.
      !==================================================================
      ! 2011-06-02 All parameters transfered via qlog and pointers modules
      subroutine setfcrg
      
      !!c      Mike Gilson's fancy charge distribution embedded into qdiffx.
      !!c      A combination of Kim Sharp's chrgit1 and Anthony Nicholls' chrgup.
      !!c      Grafting by Richard Friedman.
      !!c	 Latest version:6/2/89. First Version:4/25/89.
      real re,deb
      integer :: ic1,epsdim 
      logical isum(3*igrid),flag
      
      real :: ddist2(500), rsave(500)
      integer ::  nlist(500)
      type(coord) :: dg,xyzdist,xyzsum,ddxyz
      type(int_coord) :: kb,kxyz
      type(int_grid_value) :: clist(50,500)
      real, parameter :: crit = .0001
      
      real :: radmin=1.2, radmax=2.0,radstep=0.005
      !--------------------------------------------------------
      ! 2011-06-03 Declarations added due to IMPLICIT NONE
      integer :: i,j,k,frmnum,i1,i2,ico,ig,isgrid
      integer :: itmp,iv,itemp,ib,ix,iy,iz,kb1,kb2,kb3
      integer :: n,iw,inearest(6),kx,ky,kz,irad,cont
      integer :: ibgp,ichoice,iext,iii,il,ilist,ir,irmax
      real :: difeps,cg1,cg2,cg3,fpoh,debfct,sixeps,temp
      real :: dist,dmin,epsins6,rmax2,rmax,ssum
      real :: xdist,ydist,zdist,dist2
      
      irad=(radmax-radmin)/radstep
      !!c b+++++++++++++++++++++++
      ! 2011-06-03 Allocation of cgbp array will be done later after
      !!c e++++++++++++++++++++++
      !!c assign fractional charges to grid points, using phimap
      !!c as a dummy array (set to zero again later)
      
      epsdim=natom+nobject+2; difeps=epsin-epsout; sixeps=epsout*6.0
      fpoh=fpi*scale;  debfct = epsout/(deblen*scale)**2
      do ig=1,nqgrd
         !!c truncate to nearest grid point
         !  2011-05-16  Using operations on coord and int_coord type variables defined
         kb = int(chrgv2(ig)%xyz)
         
         !!c find position of charge in box in terms
         !!c of fractional distance along edges
         dg=chrgv2(ig)%xyz-float(kb)
         !!c loop over increasing radius to find suitable one
         
         irmax = 0; flag=.true.
         do iii=0,irad
            rmax=radmin+iii*radstep
            rmax2 =rmax*rmax;  irmax=irmax+1; rsave(irmax)=rmax; ilist=0
            ir=rmax+dg%x+1
            
            !!c loop over all grid points within range of rmax:
            !!c place base point of grid box containing the charge, at the origin, for now.
            
            do i = -ir, ir
               xdist = i - dg%x
               do j = -ir, ir
                  ydist = j - dg%y
                  do  k = -ir, ir
                     zdist = k - dg%z
                     
                     !!c calculate distance from current grid point to charge:
                     xyzdist=coord(xdist,ydist,zdist)
                     dist2=xyzdist.dot.xyzdist
                     dist = sqrt(dist2)
                     
                     !!c if grid point is closer than R, store the point and its distance:
                     
                     if (dist2.le. Rmax2) then
                        ilist = ilist + 1
                        clist(ilist,irmax)%ijk=int_coord(i,j,k)
                        clist(ilist,irmax)%value=dist
                     endif
                  end do
               end do
            end do
            nlist(irmax) = ilist
            
            !!c generate weighting factors for this rmax:
            !!c sum normalizes the weighting to one:
            
            ssum = 0
            do il = 1, nlist(irmax)
               ssum = ssum + Rmax - clist(il, irmax)%value
            end do
            !!c normalize the weighting factors to sum to one:
            do il = 1, nlist(irmax)
               clist(il, irmax)%value=clist(il,irmax)%value/ssum
               !9007              continu
            end do
            
            !!c calculate center of charge for this rmax:
            
            xyzsum=coord(0.,0.,0.)
            do il = 1, nlist(irmax)
               xyzsum=xyzsum+(clist(il,irmax)%value*float(clist(il,irmax)%ijk))
            end do
            !!c check whether criterion is satisfied, and if so, exit rmax loop.
            
            ddxyz=dg-xyzsum
            ddist2(irmax)=ddxyz.dot.ddxyz
            ! 2011-06-03 Removing GOTO and statements labels
            if (ddist2(irmax) .le.crit) then
               ichoice = irmax
               dmin = sqrt(ddist2(ichoice))
               flag=.false.; exit
            end if
            !!c otherwise, try another cutoff radius:		
         end do
         
         !!c if loop gets finished without a radius yielding good enough
         !!c results, print warning,and use the best cutoff radius found:
         ! 2011-06-03 If exit from 9002 cycle was normal
         if (flag) then
            dmin = 100000
            do i = 1, irmax
               if (ddist2(i).lt.dmin) then
                  dmin = ddist2(i); ichoice = i
               endif
            end do
            !!c---------------------------------------------------
            dmin = sqrt (dmin)
            rmax = rsave(ichoice)
         end if
         
         
         !!c now we know what set of grids we're distributing the charge over
         !!c (clist(1-3, 1-ilist(ichoice), ichoice) for ichoice), and the weighting
         !!c (clist(4,1-ilist, ichoice)...
         !!c now, distribute the charge:
         
         do ilist = 1, nlist(ichoice)
            
            !!c get grid point by adding offset to it
            
            kxyz=kb+clist(ilist,ichoice)%ijk
            
            !!c make sure grid point is within the big box (should be a problem only
            !!c in cases where box edge cuts through
            !!c or very near the protein):
            kx=kxyz%i ; ky=kxyz%j ; kz=kxyz%k
            flag = .false.
            if (kx.lt.1 ) then
               kx = 1; flag = .true.
            endif
            if (ky.lt.1) then
               ky = 1; flag = .true.
            endif
            if (kz .lt.1) then
               kz = 1; flag = .true.
            endif
            if (kx.gt.igrid) then
               kx = igrid;  flag = .true.
            endif
            if (ky.gt.igrid)then
               ky = igrid;  flag = .true.
            endif
            if (kz.gt.igrid) then
               kz = igrid; flag = .true.
            endif
            
            if(flag)then
               write(6,*)' problem for charge at', chrgv2(ig)
            end if
            
            re=clist(ilist,ichoice)%value*chrgv2(ig)%value
            phimap(kx,ky,kz)= phimap(kx,ky,kz) + re
            
         end do
         
      end do
      !!c set up odd/even logical array
      ! 2011-06-03 Changed to array operation
      isum=.false.; isum(1:3*igrid:2)=.true. ; 
      !10          isum(i)=.true.
      !20             isum(i)=.false.
      
      !!c find which grid points have charge assigned to them
      !!c (will use this array later to calculate grid energy)
      
      n=0
      do k=2,igrid-1
         do j=2,igrid-1
            do i=2,igrid-1
               if(phimap(i,j,k).ne.0) n=n+1
            end do
         end do
      end do
      allocate(gchrgp(n),gchrg(n),gchrgd(n),gchrg2(n))
      allocate(qval(n),iqpos(n),gval(n))
      
      n=0
      do  k=2,igrid-1
         do  j=2,igrid-1
            do  i=2,igrid-1
               if(phimap(i,j,k).ne.0) then
                  n=n+1
                  gchrgp(n)=int_coord(i,j,k)
                  gchrg(n)=phimap(i,j,k); phimap(i,j,k)=0.0
               end if
            end do
         end do
      end do
      icount1b=n
      ! 2011-06-03 Parameter ngcrg is no longer needed
      !!c b++++++++++++++++++++June 2001++to save memory waste time++++++
      do  ig=1,nqgrd
         ic1=nqgrdtonqass(ig)
         !!c truncate to nearest grid point
         kb=int(chrgv2(ig)%xyz)
         !!c find position of charge in box in terms
         !!c of fractional distance along edges
         dg=chrgv2(ig)%xyz-float(kb)
         
         !!c loop over increasing radius to find suitable one
         
         irmax = 0 ; irad=(radmax-radmin)/radstep
         do iii=0,irad
            rmax=radmin+iii*radstep; rmax2=rmax*rmax
            irmax=irmax+1; rsave(irmax)=rmax
            ilist=0; ir=rmax+dg%x+1
            
            !!c loop over all grid points within range of rmax:
            !!c place base point of grid box containing the charge, at the origin, for now.
            
            do i = -ir, ir
               xdist = i - dg%x
               do j = -ir, ir
                  ydist = j - dg%y
                  do  k = -ir, ir
                     zdist = k - dg%z
                     
                     !!c calculate distance from current grid point to charge:
                     
                     xyzdist=coord(xdist,ydist,zdist)
                     dist2=xyzdist.dot.xyzdist
                     dist = sqrt(dist2)
                     
                     !!c if grid point is closer than R, store the point and its distance:
                     
                     if (dist2.le. Rmax2) then
                        ilist = ilist + 1
                        clist(ilist,irmax)%ijk=int_coord(i,j,k)
                        clist(ilist,irmax)%value=dist
                     endif
                  end do
               end do
            end do
            nlist(irmax) = ilist
            
            !!c generate weighting factors for this rmax:
            !!c sum normalizes the weighting to one:
            
            ssum = 0
            do il = 1, nlist(irmax)
               ssum = ssum + Rmax - clist(il, irmax)%value
            end do
            
            !!c normalize the weighting factors to sum to one:
            
            do il = 1, nlist(irmax)
               clist(il, irmax)%value=clist(il,irmax)%value/ssum
            end do
            
            !!c calculate center of charge for this rmax:
            
            xyzsum=coord(0.,0.,0.)
            do il = 1, nlist(irmax)
               xyzsum=xyzsum+(clist(il,irmax)%value*float(clist(il,irmax)%ijk))
            end do
            !!c check whether criterion is satisfied, and if so, exit rmax loop.
            ddxyz=dg-xyzsum
            ddist2(irmax)=ddxyz.dot.ddxyz
            ! 2011-06-03 Removed GOTO and statement labels
            if (ddist2(irmax) .le.crit) then
               ichoice = irmax
               dmin = sqrt(ddist2(ichoice))
               flag=.false.; exit
            end if
            !!c otherwise, try another cutoff radius:
         end do
         
         !!c if loop gets finished without a radius yielding good enough
         !!c results, print warning,and use the best cutoff radius found:
         
         if(flag) then
            dmin = 100000
            do i = 1, irmax
               if (ddist2(i).lt.dmin) then
                  dmin = ddist2(i)
                  ichoice = i
               endif
            end do
            !!c---------------------------------------------------
            dmin = sqrt (dmin)
            rmax = rsave(ichoice)
         end if
         !!c now we know what set of grids we're distributing the charge over
         !!c (clist(1-3, 1-ilist(ichoice), ichoice) for ichoice), and the weighting
         !!c (clist(4,1-ilist, ichoice)...
         !!c now, distribute the charge:
         
         do ilist = 1, nlist(ichoice)
            !!c get grid point by adding offset to it
            kxyz=kb+clist(ilist,ichoice)%ijk
            !!c make sure grid point is within the big box (should be a problem only
            !!c in cases where box edge cuts through
            !!c or very near the protein):
            kx=kxyz%i ; ky=kxyz%j ; kz=kxyz%k
            flag = .false.
            if (kx.lt.1 ) then
               kx = 1; flag = .true.
            endif
            if (ky .lt.1) then
               ky = 1; flag = .true.
            endif
            if (kz .lt.1) then
               kz = 1; flag = .true.
            endif
            if (kx.gt.igrid) then
               kx = igrid; flag = .true.
            endif
            if (ky.gt.igrid)then
               ky = igrid; flag = .true.
            endif
            if (kz.gt.igrid) then
               kz = igrid;  flag = .true.
            endif
            
            if(flag)then
               write(6,*)' problem for charge at', &
               & chrgv2(ig)%xyz
            end if
            
            re=clist(ilist,ichoice)%value*chrgv2(ig)%value
            deb=0.
            if (idebmap(kx,ky,kz)) deb=1.
            phimap(kx,ky,kz)=phimap(kx,ky,kz)+&
            & re/(6.*atmeps(ic1)+debfct*deb)
         end do
         
      end do
      !!c e++++++++++++++++++++++++++++++++++++++++++++
      allocate(gchrgtmp(icount1b))
      !!c determine how many charged grid points are odd
      icount1a=0
      do i=1,n
         itemp=sum(gchrgp(i))
         if(isum(itemp)) icount1a=icount1a+1
      end do
      
      !!c set up odd/even pointer array, to be used in making qval
      !!c and iqpos
      i1=0; i2=icount1a
      do i=1,n
         itemp=sum(gchrgp(i))
         if(isum(itemp)) then
            i1=i1+1; gchrg2(i)=i1
         else
            i2=i2+1; gchrg2(i)=i2
         end if
      end do
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
      
      !!c determine denominator at all charged grid points
      ib=0;  epsins6= 6.0*epsin 
      do i=1,n
         ix=gchrgp(i)%i; iy=gchrgp(i)%j; iz=gchrgp(i)%k
         itemp=0
         !!c cambiato da mod a div, mi dovrebbe servire solo il mezzo
         if((iepsmp(ix,iy,iz)%i/epsdim).ne.0) itemp=itemp+1
         if((iepsmp(ix,iy,iz)%j/epsdim).ne.0) itemp=itemp+1
         if((iepsmp(ix,iy,iz)%k/epsdim).ne.0) itemp=itemp+1
         if((iepsmp(ix-1,iy,iz)%i/epsdim).ne.0) itemp=itemp+1
         if((iepsmp(ix,iy-1,iz)%j/epsdim).ne.0) itemp=itemp+1
         if((iepsmp(ix,iy,iz-1)%k/epsdim).ne.0) itemp=itemp+1
         !!c itemp=number of internal closest midpoints
         !!c b++++++++++++++++++++
         deb=0.
         if (idebmap(ix,iy,iz)) deb=1.
         if (idirectalg.eq.0)  then
            gchrgd(i)=itemp*difeps + debfct*deb + sixeps
         else
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
         end if
         !!c e++++++++++++++++++++
         if(itemp.ne.6) then
            ib=ib+1
            cgbp(ib)%xyz=coord(real(ix),real(iy),real(iz))
            !!c b++++++++++++++++++++
            cgbp(ib)%value1=gchrgtmp(i)*fpoh
            !!c e++++++++++++++++++++
            cgbp(ib)%value2=gchrg2(i)
         end if
      end do
      if(allocated(gchrgtmp)) deallocate(gchrgtmp)
      write(6,*) '# grid points charged and at boundary=',ib
	  
      !!c make qval, fpoh term so potentials will be in kt/e
      do  i=1,n
         j=gchrg2(i)
         qval(j)=gchrg(i)*fpoh/gchrgd(i)
      end do
      if(allocated(gchrgd)) deallocate(gchrgd) 
      
      !!c make iqpos
      isgrid=igrid**2
      do i=1,n
         j=gchrg2(i)
         ix=gchrgp(i)%i; iy=gchrgp(i)%j; iz=gchrgp(i)%k
         iw=1+ix+igrid*(iy-1)+isgrid*(iz-1)
         iv=iw/2
         iqpos(j)=iv
      end do
      if(allocated(gchrg2)) deallocate(gchrg2)
      
      !!c end of chrgup, return with qval,iqpos and gchrgp and gchrg
      !!c also icount1a, icount1b
      !!c b+++++++++++++++++++++++++++++++
      !!c       i_tmpmap=memalloc(i_tmpmap,0)
      !!c       i_tmpmap1=memalloc(i_tmpmap1,0)
      !!c e+++++++++++++++++++++++++++++++
      end subroutine setfcrg