      !  ATTENTION!  This file is part of setbcmod module.
      !==================================================================
      ! 2011-06-05 All parameters are accessible via qlog and pointers
      subroutine setbc
      
      !!c assigns either zero boundary conditions (ibctyp=1)
      !!c quasi-coulombic based on debye dipole qplus at cqplus
      !!c and qmin at cqmin (ibctyp=2)
      !!c focussing (ibctyp=3)(continuance if scales are the same)
      !!c full quasi-coulombic (ibctyp=4)
      !!c or constant external field (ibctyp=5)
      !!c option 2 will be appropriately modified by periodic
      !!c boundary condition flags (iper)
      !!c--------------------------------------------------------------
      ! 2011-06-05 Arrays are accessible via pointers module
      type(coord) :: g,go,c,cqnet,cqnetA,xyz,temp
      character(20) :: toblbl
      character(16) :: botlbl
      character(60) :: title
      character(10) :: label
      character(80) :: filnam
      !!c--------------------------------------------------------------
      ! 2011-06-05 Declarations added due to IMPLICIT NONE
      real :: h,ergestout,dist,fact,cutedgesx,cutedgesy,cutedgesz
      real :: cutedges,dist1,dist2,dist3,ergestemp,goff,gmid
      real :: radius,scale1,phiv,tempp,tempn,tmp,tempd,subt
      integer :: ix,iy,iz,ios,i,midg,ic,igrid1,iout,isgrid
      !-----------------------------------------------------------------
      h=1./scale; ergestout=0.0
      !!c zero option, clear boundary values
      
      if(verbose) then
         write(6,*) " "
         write(6,*) " setting boundary conditions"
         write(6,*) " "
      end if
      
      ! 2011-06-05 Changed to array operations
      phimap(1,1:igrid,1:igrid)=0.; phimap(igrid,1:igrid,1:igrid)=0.
      phimap(1:igrid,1,1:igrid)=0.; phimap(1:igrid,igrid,1:igrid)=0.
      phimap(1:igrid,1:igrid,1)=0.; phimap(1:igrid,1:igrid,igrid)=0.
      
      !!c end of zero option
      ! 2011-06-05 Variables ibctyplogions,nnit and qnet are 
      if(ibctyp.eq.6.and.(.not.(logions.and.nnit.eq.0.and.qnet.ne.0.)))&
      &   ibctyp=2
      if(ibctyp.eq.7.and.(.not.(logions.and.nnit.eq.0.and.qnet.ne.0.)))&
      &   ibctyp=4
      ! 2011-06-05 Multiple IFs changed to SELECT CASE construct
      select case(ibctyp)
      case(2)    !if(ibctyp.eq.2) then
         
         !!c quasi coulombic dipole option
         
         qnet=qmin+qplus
         do iz=1,igrid
            do iy=1,igrid
               do ix=1,igrid,igrid-1
                  !  2011-05-26  Using operations on coord and int_coord type variables defined
                  xyz=coord(real(ix),real(iy),real(iz))
                  temp=cqplus-xyz ; dist=sqrt(temp.dot.temp)/scale
                  tempp = qplus*exp(-dist/deblen )/(dist*epsout) 
                  temp=cqmin-xyz ; dist=sqrt(temp.dot.temp)/scale
                  tempn = qmin*exp(-dist/deblen )/(dist*epsout) 
                  phimap(ix,iy,iz) =  tempp + tempn
                  
               end do
            end do
         end do
         !--------------------------------------------------------------------------
         do iz=1,igrid
            do  iy=1,igrid,igrid-1
               do ix=1,igrid
                  xyz=coord(real(ix),real(iy),real(iz))
                  temp=cqplus-xyz ; dist=sqrt(temp.dot.temp)/scale
                  tempp = qplus*exp(-dist/deblen )/(dist*epsout) 
                  temp=cqmin-xyz ; dist=sqrt(temp.dot.temp)/scale
                  tempn = qmin*exp(-dist/deblen )/(dist*epsout) 
                  phimap(ix,iy,iz) =  tempp + tempn
               end do
            end do
         end do
         !--------------------------------------------------------------------------
         do iz=1,igrid,igrid-1
            do iy=1,igrid
               do ix=1,igrid
                  xyz=coord(real(ix),real(iy),real(iz))
                  temp=cqplus-xyz ; dist=sqrt(temp.dot.temp)/scale
                  tempp = qplus*exp(-dist/deblen )/(dist*epsout) 
                  temp=cqmin-xyz ; dist=sqrt(temp.dot.temp)/scale
                  tempn = qmin*exp(-dist/deblen )/(dist*epsout) 
                  phimap(ix,iy,iz) =  tempp + tempn
               end do
            end do
         end do
         
         !!c end of quasi coulombic dipole option
         
         !--------------------------------------------------------------------------
      case(4)   !        if(ibctyp.eq.4) then
         
         !!c a summation of the potential resulted from each point of charge 
         
         qnet=qmin+qplus
         do iz=1,igrid
            do iy=1,igrid
               do ix=1,igrid,igrid-1
                  xyz=coord(real(ix),real(iy),real(iz))
                  do ic = 1,nqass
                     temp=atmcrg(ic)%xyz-xyz
                     dist=sqrt(temp.dot.temp)/scale
                     tempd = atmcrg(ic)%value*exp(-dist/deblen)/(dist*epsout)
                     phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempd 
                  end do
               end do
            end do
         end do
         !--------------------------------------------------------------------------
         do iz=1,igrid
            do iy=1,igrid,igrid-1
               do ix=1,igrid
                  xyz=coord(real(ix),real(iy),real(iz))
                  do ic = 1,nqass
                     temp=atmcrg(ic)%xyz-xyz
                     dist=sqrt(temp.dot.temp)/scale
                     tempd = atmcrg(ic)%value*exp(-dist/deblen)/(dist*epsout)
                     phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempd 
                  end do
               end do
            end do
         end do
         !--------------------------------------------------------------------------
         do  iz=1,igrid,igrid-1
            do  iy=1,igrid
               do  ix=1,igrid
                  xyz=coord(real(ix),real(iy),real(iz))
                  do  ic = 1,nqass
                     temp=atmcrg(ic)%xyz-xyz
                     dist=sqrt(temp.dot.temp)/scale
                     tempd = atmcrg(ic)%value*exp(-dist/deblen)/(dist*epsout)
                     phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempd 
                  end do
               end do
            end do
         end do
         
         !!c end of the option for the complete summation of potential
         
         !--------------------------------------------------------------------------
      case(3) !         if(ibctyp.eq.3) then
         
         !!c focussing option-bc's come from a previous phimap
         
         open(18,file=phiinam(:phiilen),status='old',&
         & iostat=ios,form='unformatted')
         if(ios.ne.0) then
            write(6,*)' no potl map for focussing boundary conditions'
            stop
         end if
         write(6,*)' '
         write(6,*)'focussing boundary condition '
         write(6,*)'read from file'
         write(6,*) trim(phiinam)
         write(6,*)' '
         
         !!c read in old potential map
         
         read(18)toblbl
         read(18)label,title
         read(18)phimap
         read(18)botlbl
         read(18)scale1,oldmid1,igrid1
         close(18)
         
         if(scale1.eq.scale) then
            write(6,*) 'scales are the same.' 
            write(6,*) 'therefore assuming this to be a continuence'
         else
            write(6,*)' '
            write(6,*)' focussing potential map:'
            write(6,'(a60)') trim(title)
            write(6,*)'original scale (grids/A)      : ',scale1
            write(6,*)'object centre at (A) : ',oldmid1
            write(6,*)' '
            !!c check to see that new grid lies within old one that is going to
            !!c provide bc's
            
            iout = 0
            !!c             goff = (ngrid+1.)/2.
            goff = (igrid1+1.)/2.
            do  iz=1,igrid,igrid-1
               do  iy=1,igrid,igrid-1
                  do  ix=1,igrid,igrid-1
                     g=coord(real(ix),real(iy),real(iz))
                     !!c for each new grid corner, calculate old grid coords
                     c=((g-goff)/scale)+oldmid
                     temp=((c-oldmid1)*scale1)+goff
                     if((temp.vorle.1.).or.(temp.vorge.real(igrid1))) iout=1
                     !!c line below replace 9/16/91 cos i think it should be between
                     !!c 1 and ngrid, not 2 and ngrid-1
                     !!c             if((gold.le.2.).or.(gold.ge.ngrid-1.))  iout = 1
                     !!c             if((gold.le.1).or.(gold.ge.ngrid))  iout = 1
                  end do
               end do
            end do
            if(iout.ne.0) then
               write(6,*)'part of new grid lies outside old grid'
               write(6,*)'check scaling of both grids'
               write(6,*)'old scale:'
               write(6,*)'scale (grids/A)      : ',scale1
               write(6,*)'object centre at (A) : ',oldmid1
               write(6,*)'new scale:'
               write(6,*)'scale (grids/A)      : ',scale
               write(6,*)'object centre at (A) : ',oldmid
               stop
            end if
            
            !!c for each boundary point
            !!c convert to real coordinates
            !!c convert to old grid coordinates
            !!c interpolate potential
            !!c note that can use same potential array for boundaries
            !!c since old potentials at boundary are not used for new ones
            
            !!c save new grid size, and set temporarily to 65
            isgrid = igrid
            igrid = igrid1
            gmid = (isgrid + 1.)/2.
            write(6,*)'pulling boundary values out of old potential map...'
            do iz=2,isgrid-1
               do iy=2,isgrid-1
                  do  ix=1,isgrid,isgrid-1
                     g=coord(real(ix),real(iy),real(iz))
                     c=((g-gmid)/scale)+oldmid
                     go=((c-oldmid1)*scale1)+goff
                     !!c for each new grid side, calculate old grid coords
                     !!c find potential
                     call phintp(go,phiv)
                     phimap(ix,iy,iz) = phiv
                  end do
               end do
            end do
            !--------------------------------------------------------------------------
            do iz=2,isgrid-1
               do  iy=1,isgrid,isgrid-1
                  do  ix=2,isgrid-1
                     
                     !!c for each new grid side, calculate old grid coords
                     g=coord(real(ix),real(iy),real(iz))
                     c=((g-gmid)/scale)+oldmid
                     go=((c-oldmid1)*scale1)+goff
                     !!c find potential
                     call phintp(go,phiv)
                     phimap(ix,iy,iz) = phiv
                  end do
               end do
            end do
            !--------------------------------------------------------------------------
            do  iz=1,isgrid,isgrid-1
               do iy=2,isgrid-1
                  do  ix=2,isgrid-1
                     !!c for each new grid side, calculate old grid coords
                     g=coord(real(ix),real(iy),real(iz))
                     c=((g-gmid)/scale)+oldmid
                     go=((c-oldmid1)*scale1)+goff
                     !!c find potential
                     call phintp(go,phiv)
                     phimap(ix,iy,iz) = phiv
                  end do
               end do
            end do
            !!c restore new grid size
            igrid = isgrid
            !511         end if
         end if
         
         !!c end of focussing option
         
         !--------------------------------------------------------------------------
      case(5)     !  if(ibctyp.eq.5) then
         fact=0.
         subt=vdrop%z*(fact-1.)/2.
         !!c     fixed potential over z direction
         do  iz=2,igrid-1
            do  ix=1,igrid
               do  iy=1,igrid
                  phimap(ix,iy,iz)=vdrop%z*fact/(1.+&
                  & exp((-2.*iz+igrid+1.)/35.5))-subt
               end do
            end do
         end do
         
         do ix=1,igrid
            do iy=1,igrid
               phimap(ix,iy,1)    = 0.0
               phimap(ix,iy,igrid)= vdrop%z
            end do
         end do
         
         
         !!c end external field option
         !--------------------------------------------------------------------------
      case(6)    !  if(ibctyp.eq.6) then
         !!c quasi coulombic modified dipole option
         qnet=qmin+qplus
         !!c b+++++++++++++++++++++++++++++++++++++++++++++++++++++
         i=2
         if (qmin.eq.0.or.qplus.eq.0) i=1
         goff=(igrid+1.)/2.
         cqnet=(cqmin+cqplus)/real(i)
         cqnetA=((cqnet-goff)/scale)+oldmid
         !!c     calcolo  raggio medio molecola
         tmp=0.
         do i=1,ibnum
            temp=scspos(i)-cqnetA
            tmp=tmp+(temp.dot.temp)
         end do
         radius=sqrt(tmp/ibnum)+exrad
         !------------------------------------------------------------
         ergestemp=0.0
         !!c setting bc and calculating direct ionic contribution
         do  iz=1,igrid
            cutedgesz=-2.5
            if (iz.eq.1.or.iz.eq.igrid) cutedgesz=2.
            do iy=1,igrid
               cutedges=cutedgesz
               if (iy.eq.1.or.iy.eq.igrid) cutedges=cutedgesz+1
               ! 2011-06-05 Added additional loop over ix= 1 or igrid
               do ix=1,igrid,igrid-1
                  xyz=coord(real(ix),real(iy),real(iz))
                  temp=cqplus-xyz; temp%x=-temp%x
                  dist=sqrt(temp.dot.temp)/scale
                  tempp = qplus*exp(-dist/deblen)/(dist*epsout)
                  
                  temp=cqmin-xyz; temp%x=-temp%x
                  dist=sqrt(temp.dot.temp)/scale
                  tempn = qmin*exp(-dist/deblen )/(dist*epsout)
                  phimap(ix,iy,iz) =  tempp + tempn
                  
                  
                  temp=cqnet-xyz; temp%x=-temp%x
                  dist=sqrt(temp.dot.temp)/scale
                  fact=exp(-(dist-radius)/deblen)/(1+radius/deblen)
                  dist1=temp%x*(3.5-cutedges)/6.
                  tempd=qnet*fact/(dist*epsout)
                  if(ix.eq.1) then
                     ergestemp=ergestemp+qnet*tempd*dist1/dist**2
                  else
                     ergestemp=ergestemp-qnet*tempd*dist1/dist**2
                  end if
               end do
            end do
         end do
         if (.not.iper(1)) ergestout=ergestout+ergestemp
         !-----------------------------------------------------------
         ergestemp=0.0
         do iz=1,igrid
            cutedgesz=-2.5
            if (iz.eq.1.or.iz.eq.igrid) cutedgesz=2.
            do ix=1,igrid
               cutedges=cutedgesz
               if (ix.eq.1.or.ix.eq.igrid) cutedges=cutedgesz+1
               ! 2011-06-05 Added additional loop over iy = 1 or igrid
               do iy=1,igrid,igrid-1
                  xyz=coord(real(ix),real(iy),real(iz))
                  temp=cqplus-xyz; temp%y=-temp%y
                  dist=sqrt(temp.dot.temp)/scale
                  tempp = qplus*exp(-dist/deblen)/(dist*epsout)
                  
                  temp=cqmin-xyz; temp%y=-temp%y
                  dist=sqrt(temp.dot.temp)/scale
                  tempn = qmin*exp(-dist/deblen )/(dist*epsout)
                  phimap(ix,iy,iz) =  tempp + tempn
                  
                  fact=exp(-(dist-radius)/deblen)/(1+radius/deblen)
                  dist2=temp%y*(3.5-cutedges)/6.
                  tempd=qnet*fact/(dist*epsout)
                  if(iy.eq.1) then
                     ergestemp=ergestemp+qnet*tempd*dist2/dist**2
                  else
                     ergestemp=ergestemp-qnet*tempd*dist2/dist**2
                  end if
                  
                  
               end do
            end do
         end do
         if (.not.iper(2)) ergestout=ergestout+ergestemp
         !--------------------------------------------------------------------------
         ergestemp=0.0
         do iy=1,igrid
            cutedgesy=-2.5
            if (iy.eq.1.or.iy.eq.igrid) cutedgesy=2.
            do ix=1,igrid
               cutedges=cutedgesy
               if (ix.eq.1.or.ix.eq.igrid) cutedges=cutedgesy+1
               do iz=1,igrid,igrid-1
                  xyz=coord(real(ix),real(iy),real(iz))
                  temp=cqplus-xyz; temp%z=-temp%z
                  dist=sqrt(temp.dot.temp)/scale
                  tempp = qplus*exp(-dist/deblen )/(dist*epsout)
                  
                  temp=cqmin-xyz; temp%z=-temp%z
                  dist=sqrt(temp.dot.temp)/scale
                  tempn = qmin*exp(-dist/deblen )/(dist*epsout)
                  phimap(ix,iy,iz) = tempp + tempn
                  
                  fact=exp(-(dist-radius)/deblen)/(1+radius/deblen)
                  dist3=temp%z*(3.5-cutedges)/6.
                  tempd=qnet*fact/(dist*epsout)
                  if(iz.eq.1) then
                     ergestemp=ergestemp+qnet*tempd*dist3/dist**2
                  else
                     ergestemp=ergestemp-qnet*tempd*dist3/dist**2
                  end if
               end do
            end do
         end do
         if (.not.iper(3)) ergestout=ergestout+ergestemp
         
         ergestout=ergestout*rionst*deblen*.0006023/(epsout*scale**3)
         !!c end of quasi coulombic modified dipole option
         !--------------------------------------------------------------------------
      case(7)   !        if(ibctyp.eq.7) then
         !!c a summation of the potential resulted from each point of charge 
         qnet=qmin+qplus
         !!c here sets bc and calculates ionic direct contribution at the same time
         i=2
         if (qmin.eq.0.or.qplus.eq.0) i=1
         goff=(igrid+1.)/2.
         cqnet=(cqmin+cqplus)/real(i)
         cqnetA=((cqnet-goff))/scale+oldmid
         !!c     calcolo  raggio medio molecola
         tmp=0.
         do i=1,ibnum
            temp=scspos(i)-cqnetA
            tmp=tmp+(temp.dot.temp)
         end do
         radius=sqrt(tmp/ibnum)+exrad
         
         !--------------------------------------------------------------------------
         ergestemp=0.0
         do  iz=1,igrid
            cutedgesz=-2.5
            if (iz.eq.1.or.iz.eq.igrid) cutedgesz=2.
            do  iy=1,igrid
               cutedges=cutedgesz
               if (iy.eq.1.or.iy.eq.igrid) cutedges=cutedgesz+1
               do ix=1,igrid,igrid-1
                  xyz=coord(real(ix),real(iy),real(iz))
                  do  ic = 1,nqass
                     temp=atmcrg(ic)%xyz-xyz; temp%x=-temp%x
                     dist=sqrt(temp.dot.temp)/scale
                     tempd=atmcrg(ic)%value*exp(-dist/deblen)/(dist*epsout)
                     phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempd
                  end do
                  
                  temp=cqnet-xyz ; temp%x=-temp%x
                  dist=sqrt(temp.dot.temp)/scale
                  fact=exp(-(dist-radius)/deblen)/(1+radius/deblen)
                  dist1=temp%x*(3.5-cutedges)/6.
                  if(ix.eq.1) then
                     ergestemp=ergestemp+qnet*tempd*dist1/dist**2
                  else
                     ergestemp=ergestemp-qnet*tempd*dist1/dist**2
                  end if
                  
               end do
            end do
         end do
         if (.not.iper(1)) ergestout=ergestout+ergestemp
         !--------------------------------------------------------------------------
         ergestemp=0.0
         do  iz=1,igrid
            cutedgesz=-2.5
            if (iz.eq.1.or.iz.eq.igrid) cutedgesz=2.
            do  ix=1,igrid
               cutedges=cutedgesz
               if (ix.eq.1.or.ix.eq.igrid) cutedges=cutedgesz+1
               do iy=1,igrid,igrid-1
                  xyz=coord(real(ix),real(iy),real(iz))
                  do ic = 1,nqass
                     temp=atmcrg(ic)%xyz-xyz; temp%y=-temp%y
                     dist=sqrt(temp.dot.temp)/scale
                     tempd=atmcrg(ic)%value*exp(-dist/deblen)/(dist*epsout)
                     phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempd
                     
                  end do
                  temp=cqnet-xyz ; temp%y=-temp%y
                  dist=sqrt(temp.dot.temp)/scale
                  fact=exp(-(dist-radius)/deblen)/(1+radius/deblen)
                  dist2=temp%y*(3.5-cutedges)/6.
                  if(iy.eq.1) then
                     ergestemp=ergestemp+qnet*tempd*dist2/dist**2
                  else
                     ergestemp=ergestemp-qnet*tempd*dist2/dist**2
                  end if
                  1                     dist2=dist2*(3.5-cutedges)/6.
               end do
            end do
         end do
         if (.not.iper(2)) ergestout=ergestout+ergestemp
         !--------------------------------------------------------------------------
         ergestemp=0.0
         do iy=1,igrid
            cutedgesy=-2.5
            if (iy.eq.1.or.iy.eq.igrid) cutedgesy=2.
            do ix=1,igrid
               cutedges=cutedgesy
               if (ix.eq.1.or.ix.eq.igrid) cutedges=cutedgesy+1
               do iz=1,igrid,igrid-1
                  xyz=coord(real(ix),real(iy),real(iz))
                  do ic = 1,nqass
                     temp=atmcrg(ic)%xyz-xyz; temp%z=-temp%z
                     dist=sqrt(temp.dot.temp)/scale
                     tempd=atmcrg(ic)%value*exp(-dist/deblen)/(dist*epsout)
                     phimap(ix,iy,iz) = phimap(ix,iy,iz) + tempd
                     
                  end do
                  temp=cqnet-xyz ; temp%z=-temp%z
                  dist=sqrt(temp.dot.temp)/scale
                  fact=exp(-(dist-radius)/deblen)/(1+radius/deblen)
                  dist3=temp%z*(3.5-cutedges)/6.
                  if(iz.eq.1) then
                     ergestemp=ergestemp+qnet*tempd*dist3/dist**2
                  else
                     ergestemp=ergestemp-qnet*tempd*dist3/dist**2
                  end if
               end do
            end do
         end do
         if (.not.iper(3)) ergestout=ergestout+ergestemp
         ergestout=ergestout*rionst*deblen*.0006023/(epsout*scale**3)
      end select
      !!c end of the option for the complete summation of potential
      !!c e+++++++++++++++++++++++++++++++++++++++++++++++++++++++
      midg = (igrid+1)/2
      if(verbose) then
         
         if(ideveloper) then
            write(6,*)' some initial phi values: '
            write(6,*)' midg,midg,1; midg,midg,igrid '
            write(6,"(2e23.15)")phimap(midg,midg,1), phimap(midg,midg,igrid)
            write(6,*)' midg,1,midg; midg,igrid,midg '
            write(6,"(2e23.15)")phimap(midg,1,midg),phimap(midg,igrid,midg)
            write(6,*)' 1,midg,midg; igrid,midg,midg '
            write(6,"(2e23.15)")phimap(1,midg,midg),phimap(igrid,midg,midg)
         else
            write(6,*)' some initial phi values: '
            write(6,*)' midg,midg,1; midg,midg,igrid '
            write(6,"(2e13.5)")phimap(midg,midg,1), phimap(midg,midg,igrid)
            write(6,*)' midg,1,midg; midg,igrid,midg '
            write(6,"(2e13.5)")phimap(midg,1,midg),phimap(midg,igrid,midg)
            write(6,*)' 1,midg,midg; igrid,midg,midg '
            write(6,"(2e13.5)")phimap(1,midg,midg),phimap(igrid,midg,midg)
         end if
      end if
      !!C---------------------------------------------------------------
      !900      write(6,*)' no potl map for focussing boundary conditions'
      end subroutine setbc
