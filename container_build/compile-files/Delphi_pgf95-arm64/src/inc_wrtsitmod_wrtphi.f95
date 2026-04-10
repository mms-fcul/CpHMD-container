      !  ATTENTION!  This file is part of wrtsitmod module.
      !==================================================================
      subroutine wrtphi
      
      !!cNB uses phimap3 as a temp storage for phimap in case
      !!c want to write not potentials bu salt concentraions
      integer, parameter :: mgrid=65
      !!c b+++++++++++for back compatibility with single precision phi readers
      real :: scalesingle,scalesingle1
      type(coord) :: oldmidsingle,oldmidsingle1,xyzstart,xyzend
      type(coord) :: origin
      real :: minim,maxim,somma,average
      integer :: igrid1      
      character(80) :: filnam
      character(10) :: nxtlbl
      character(20) :: uplbl
      character(16) :: botlbl
      !---------------------------------------------------------------------------
      ! 2011-06-15 Declarations added due to IMPLICIT NONE
      integer :: i,j,k,ivary,ix,iy,iz,ii,intx,inty,intz,intdat,nbyte
      real :: extent,range,xmax,xang,yang,zang,stepsize,coeff
      !---------------------------------------------------------------------------
      scalesingle=scale
      oldmidsingle=oldmid
      
      if (realsiz.ne.4.and.phifrm.ne.2) then
         allocate(phimap4(igrid*igrid*igrid))
         i=1
         do iz=1,igrid
            do iy=1,igrid
               do ix=1,igrid
                  phimap4(i)=phimap(ix,iy,iz)
                  i=i+1
               end do
            end do
         end do
      end if
      write(6,*) phimap(5,5,5)
      write(6,*) phimap4(5+(5-1)*igrid+(5-1)*igrid*igrid)
      !!c e++++++++++++++++++++++++++++
      if(iconc) then
         allocate(phimap3(ngp))
         i=1
         do iz=1,igrid
            do iy=1,igrid
               do ix=1,igrid
                  phimap3(i)=phimap(ix,iy,iz)
                  i=i+1
               end do
            end do
         end do
         call phicon
      end if
      
      if(ibios) then
         
         !!c write phimap in insight format
         open(14,file=trim(phinam),form="unformatted")
         filnam = ' '
         inquire(14,name = filnam)
         write(6,*)'potential map written in INSIGHT format to file'
         write(6,*)trim(filnam)
         write(6,*)'  '
         ivary = 0; nbyte = 4; intdat = 0
         xang = 90.; yang = 90.; zang = 90.
         intx = igrid - 1
         inty = igrid - 1
         intz = igrid - 1
         !  2011-06-15  Using operations on coord type variables defined
         xmax=max(oldmid)
         range = (igrid-1.)/(2.*scale)
         extent = range + xmax
         xyzstart=(oldmid-range)/extent
         xyzend=  (oldmid+range)/extent
         write(14)toplbl
         write(14)ivary,nbyte,intdat,extent,extent,extent,&
         & xang,yang,zang,xyzstart%x,xyzend%x,&
         &xyzstart%y,xyzend%y,xyzstart%z,xyzend%z,intx,inty,intz
         if (realsiz.ne.4) then
            call wrtphimap(igrid,phimap4,1)
         else
            do k = 1,igrid
               do j = 1,igrid
                  write(14)(phimap(i,j,k),i=1,igrid)
               end do
            end do
         end if
      elseif(phifrm.eq.2)then
         !!c GRASP phimap - output a 65^3 grid and leave out ngrid spec.
         
         write(6,*)'  '
         write(6,*)'writing potential map in GRASP format'
         write(6,*)'  '
         
         allocate(phimap4(mgrid*mgrid*mgrid))
         call expand(mgrid,phimap4)
         open(14,file=trim(phinam),form="unformatted")
         filnam = ' '
         inquire(14,name = filnam)
         write(6,*)'potential map written to file'
         write(6,*)trim(filnam)
         write(6,*)'  '
         if(iconc.and.(rionst.ne.0)) then
            nxtlbl="concentrat"
         else
            nxtlbl="potential "
         end if
         write(14)'now starting phimap '
         write(14)nxtlbl,toplbl
         !!c        write(14)phimap4
         call wrtphimap(mgrid,phimap4,0)
         write(14)' end of phimap  '
         write(14)scalesingle,oldmidsingle
         close(14)
         deallocate(phimap4)
         !!c b++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      elseif(phifrm.eq.3)then
         !!c CCP4 phimap - output 
         
         write(6,*)'  '
         write(6,*)'writing potential map in CCP4 format'
         write(6,*)'  ' 
         !!c       '(A20)'
         open(14,file=trim(phinam),form="unformatted")
         filnam = ' '
         inquire(14,name = filnam)
         write(6,*)'potential map written to file'
         write(6,*)trim(filnam)
         write(6,*)'  '
         !!c NC,NR,NS,MODE,NCSTART,NRSTART,NSSTART,NX,NY,NZ
         write(14)igrid,igrid,igrid,2,1,1,1,igrid-1,igrid-1,igrid-1
         !!c X Y Z lengths, alpha
         write(14)(igrid-1)/scale,(igrid-1)/scale,(igrid-1)/scale,90.0
         !!c beta, gamma,MAPC,MAPR,MAPS 
         write(14)90.0,90.0,1,2,3
         
         minim=10000.0;  maxim=-10000.0
         if (realsiz.ne.4) then
            i=1
            do iz=1,igrid
               do iy=1,igrid
                  do ix=1,igrid
                     if (phimap4(i).lt.minim) minim=phimap4(i)
                     if (phimap4(i).gt.maxim) maxim=phimap4(i)
                     somma=somma+phimap4(i)
                     i=i+1
                  end do
               end do
            end do
            average=somma/(igrid*igrid*igrid)
         else
            do iz=1,igrid
               do iy=1,igrid
                  do ix=1,igrid
                     if (phimap(ix,iy,iz).lt.minim) minim=phimap(ix,iy,iz)
                     if (phimap(ix,iy,iz).gt.maxim) maxim=phimap(ix,iy,iz)
                     somma=somma+phimap(ix,iy,iz)
                  end do
               end do
            end do
            average=somma/(igrid*igrid*igrid)
         end if
         !!c amin, amax, amean,ispg,nsymbt,LSKFLG
         write(14)minim,maxim,average,1,1,0
         !!c skwmat(3,3)
         write(14)0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
         !!c skwtrn(3),future use
         write(14)0.0,0.0,0.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
         write(14)'MAP '
         !!c MACHST,ARMS,NLABL
         write(14)1,0.0,0
         do ii=1,200
            write(14)0
         end do
         if (realsiz.ne.4) then
            call wrtphimap(igrid,phimap4,0)
         else
            write(14)phimap
         end if
         close(14)
         
      elseif(phifrm.eq.4)then
         !!c diffrential phimap - 
         
         write(6,*)'  '
         write(6,*)'reading previous.phi file'
         write(6,*) realsiz, phimap(5,5,5)
         write(6,*) phimap4(5+(5-1)*igrid+(5-1)*igrid*igrid)
         write(6,*)'  '
         
         open(14,file="previous.phi",form="unformatted")
         read(14)uplbl
         read(14)nxtlbl,toplbl
         write(6,*)uplbl,nxtlbl,toplbl
         call rdphimap(igrid,phimap4)
         read(14)botlbl
         read(14)scalesingle1,oldmidsingle1,igrid1
         write(6,*)botlbl,scalesingle1,oldmidsingle1
         close(14)
         if ((scalesingle1.ne.scalesingle).or.&
         &(oldmidsingle1.vorne.oldmidsingle).or.&
         &(igrid1.ne. igrid)) then
            write(6,*) 'Error: the two potential maps do not mach'
            return
         endif
         write(6,*)phimap4(5+(5-1)*igrid+(5-1)*igrid*igrid)
         if (realsiz.ne.4) then
            call submaps(igrid,phimap,phimap4,2)
         else
            call submaps(igrid,phimap,phimap4,1)
         end if
         
         open(14,file=trim(phinam),form="unformatted")
         filnam = ' '
         inquire(14,name = filnam)
         write(6,*)'potential map written to file'
         write(6,*)trim(filnam)
         write(6,*)'  '
         if(iconc.and.(rionst.ne.0)) then
            nxtlbl="concentrat"
         else
            nxtlbl="potential "
         end if
         !!c         the following is uplabel, 20 char
         write(14)'now starting phimap '
         write(14)nxtlbl,toplbl
         if (realsiz.ne.4) then
            call wrtphimap(igrid,phimap4,0)
         else
            write(14)phimap
         end if
         !!c         this is botlbl
         write(14)' end of phimap  '
         write(14)scalesingle,oldmidsingle,igrid
         close(14)
         !!c e++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! 2011-07-23 New piece of code from Maxim
      elseif (phifrm.eq.5) then
         open(14,file=phinam(:philen), form='formatted')
         filnam = ' '
         inquire(14,name = filnam)
         write(6,*)' Potential map in cube format '
         write(6,*)'written to file',filnam
         
         write(14,*)'qdiffxs4 with an improved surfacing routine'
         write(14,*) 'Gaussian cube format phimap'
         coeff=0.5291772108
         stepsize=1.0/scale
         origin=(oldmid-stepsize*(igrid-1)/2)/coeff
         write(14,'(i5,3f14.6)') 1, origin !06/12/12 Lin: 12.6 to 14.6, for large grids.
         write(14,'(i5,3f12.6)') igrid, stepsize/coeff,0.0,0.0
         write(14,'(i5,3f12.6)') igrid, 0.0,stepsize/coeff,0.0
         write(14,'(i5,3f12.6)') igrid, 0.0,0.0,stepsize/coeff
         write(14,'(i5,4f12.6)') 1,0.0,0.0,0.0,0.0
         do i = 1,igrid
            do j = 1,igrid
               write(14,'(6E13.5)')(phimap(i,j,k),k=1,igrid)
            end do
         end do
         close(14)
         ! 2011-07-23 End new piece of code from Maxim
         !--------------------------------------------------------------------
      else
         write(6,*)'  '
         write(6,*)'writing potential map in DELPHI format'
         write(6,*) realsiz, phimap(5,5,5)
         write(6,*) phimap4(5+(5-1)*igrid+(5-1)*igrid*igrid)
         write(6,*)'  '
         
         open(14,file=trim(phinam),form="unformatted")
         filnam = ' '
         inquire(14,name = filnam)
         write(6,*)'potential map written to file'
         write(6,*)trim(filnam)
         write(6,*)'  '
         if(iconc.and.(rionst.ne.0)) then
            nxtlbl="concentrat"
         else
            nxtlbl="potential "
         end if
         !!c         the following is uplabel, 20 char
         write(14)'now starting phimap '
         write(14)nxtlbl,toplbl
         !debug tested              write(6,*)nxtlbl,len(nxtlbl)
         ! 2011-07-21 Just causing trouble with F95
         write(14)phimap
         !!c         this is botlbl
         write(14)' end of phimap  '
         write(14)scalesingle,oldmidsingle,igrid
         close(14)
      end if
      	
      if(realsiz.ne.4.and.phifrm.ne.2) deallocate(phimap4)
      
      if(iconc) then
         i=1
         do iz=1,igrid
            do iy=1,igrid
               do ix=1,igrid
                  phimap(ix,iy,iz)=phimap3(i)
                  i=i+1
               end do
            end do
         end do
         deallocate(phimap3)
      end if
      
      end subroutine wrtphi
      !=======================================================================
      subroutine wrtphimap(mgrid,phimap,opt)
      real :: phimap(mgrid,mgrid,mgrid)
      integer :: opt, mgrid,i,j,k
      
      if (opt.eq.1) then
         do  k = 1,mgrid
            do  j = 1,mgrid
               write(14)(phimap(i,j,k),i=1,mgrid)
            end do
         end do
      else
         write(14)phimap
      endif
      end subroutine wrtphimap
      !=======================================================================
      subroutine rdphimap(mgrid,phimap)
      integer :: mgrid
      real :: phimap(mgrid,mgrid,mgrid)
      read(14)phimap
      write(6,*) 'inside rdphimap',phimap(5,5,5)
      return
      end subroutine rdphimap
      !=======================================================================
      subroutine submaps(mgrid,phimap,phimap4,opt)
      integer :: mgrid,opt,ix,iy,iz
      real :: phimap(mgrid,mgrid,mgrid)
      real :: phimap4(mgrid,mgrid,mgrid)
      
      if(opt.eq.1) then
         do iz=1,mgrid
            do iy=1,mgrid
               do ix=1,mgrid
                  phimap(ix,iy,iz)=phimap(ix,iy,iz)-phimap4(ix,iy,iz)
               end do
            end do
         end do
      else
         do iz=1,mgrid
            do iy=1,mgrid
               do ix=1,mgrid
                  phimap4(ix,iy,iz)=phimap(ix,iy,iz)-phimap4(ix,iy,iz)
               end do
            end do
         end do
      endif
      
      end subroutine submaps
      !=======================================================================
      
      subroutine phicon
      !---------------------------------------------------------------------------
      ! 2011-06-15 Declarations added due to IMPLICIT NONE
      integer :: ix,iy,iz
      real :: sixth, tmp,temp,phi,phisq
      !---------------------------------------------------------------------------
      !!c convert potentials to concentrations
      
      sixth = 1./6.; tmp=abs(chi2*chi4)
      if(rionst.gt.0.0) then
         write(6,*)'  '
         write(6,*)'converting potentials to '
         write(6,*)'net charge concentrations...'
         write(6,*)'  '
         !!c IF nonlinear equation is used then use the exponential form
         !!c otherwise use the linear form
         !!c NB use same number of terms in expansion
         !!c of sinh as for iteration in itit.f
         write(6,*)'PHICON: this option has not been tested yet'
         
         if(nnit.ne.0) then
            if (tmp.lt.1.e-6) then
               do iz = 1,igrid
                  do  iy = 1,igrid
                     do ix = 1,igrid
                        !!c use first three terms of sinh
                        if (idebmap(ix,iy,iz)) then
                           phi = phimap(ix,iy,iz)
                           phisq = phi*phi
                           !!c Horner scheme for charge and osmotic term
                           temp = phisq*chi5 + chi3
                           temp = temp*phisq + chi1
                           phimap(ix,iy,iz)=temp*phi     
                        else
                           phimap(ix,iy,iz)=0.0
                        end if
                     end do
                  end do
               end do
            else
               !!c b+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               !!c asymmetric salt
               do iz = 1,igrid
                  do iy = 1,igrid
                     do ix = 1,igrid
                        if (idebmap(ix,iy,iz)) then
                           phi = phimap(ix,iy,iz)
                           !!c Horner scheme for charge and osmotic term
                           temp = phi*chi5 + chi4
                           temp = phi*temp + chi3
                           temp = phi*temp + chi2
                           temp = temp*phi + chi1
                           phimap(ix,iy,iz)=temp*phi   
                        else
                           phimap(ix,iy,iz)=0.0
                        end if
                     end do
                  end do
               end do
            end if
            !!c e++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         else
            do  iz = 1,igrid
               do  iy = 1,igrid
                  do ix = 1,igrid
                     if (idebmap(ix,iy,iz)) then
                        phi = phimap(ix,iy,iz)
                        phimap(ix,iy,iz)=chi1*phi
                     else
                        phimap(ix,iy,iz)=0.0
                     end if
                  end do
               end do
            end do
         end if
      end if
      
      if(rionst.eq.0.0) then
         write(6,*) "cannot convert from potentials to concentrations"
         write(6,*) "if the ionic strenth is zero!"
      end if
      end subroutine phicon
	  
      !===========================================================================
      subroutine expand(mgrid,phimap4)
      !!c expands igrid**3 grid to ngrid**3
      !!c grid (65**3 default) to give compatibilty with previous
      !!c phimap and epsmap formats using trilinear interpolation
      real :: phimap4(mgrid,mgrid,mgrid)
      type(coord) :: gc
      ! 2011-06-15 Array temp is not used
      !---------------------------------------------------------------------------
      ! 2011-06-15 Declarations added due to IMPLICIT NONE
      integer :: ix,iy,iz,mgrid
      real :: phiv,rscale
	  
      rscale = (igrid-1.)/(mgrid-1.)
      !!c do high end first to prevent overwriting
      !!c find small grid values and interpolate into big grid
      
      if(igrid.eq.mgrid)then
         do iz=1,igrid
            do iy=1,igrid
               do ix=1,igrid
                  phimap4(ix,iy,iz) = phimap(ix,iy,iz)
               end do
            end do
         end do
      else
         if(.not.ibios) then
            do  iz = mgrid,1,-1
               gc%z = (iz-1)*rscale + 1.
               do  iy = mgrid,1,-1
                  gc%y = (iy-1)*rscale + 1.
                  do  ix = mgrid,1,-1
                     gc%x = (ix-1)*rscale + 1.
                     call phintp(gc,phiv)
                     phimap4(ix,iy,iz) = phiv
                  end do
               end do
            end do
         end if
      end if
      
      scale = scale/rscale
      write(6,*)'new scale is ',scale,' grids/ang'
      end subroutine expand