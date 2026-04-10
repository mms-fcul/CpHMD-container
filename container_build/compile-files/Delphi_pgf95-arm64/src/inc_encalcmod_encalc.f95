!#####################################################################
!ATTENTION!  This file is part of encalcmod module.
!            Do not compile separately!
!
! 2011-06-11 All parameters are accessible via qlog and pointers
!            modules
! 2011-06-11 Arrays and variables below are declared in qlog and
!            pointers modules
!            real phimap(igrid,igrid,igrid)
!            dimension gchrg(icount1b)
!            dimension rad3(natom),xn2(3,natom),scspos(3,ibnum)
!            integer gchrgp(3,icount1b),ibgrd(3,ibnum),nqgrd
!            real medeps(0:nmedia)
!
! comments on energies
!    logs= solvation energy
!    logc= coulombic energy
!    loga= analytic energy
!    logg= grid energy
!    logions= ionic contribution
!    lognl=ionic contribution in non linear case
!#####################################################################

      subroutine encalc

      logical :: ido
      integer :: igridout
      !real*8 :: ergs,ergas,ergnl,ergc,ergg
      type(int_extrema) :: lim
      type(int_coord) :: ixyz
      !06/12/12Lin: changed ergs,ergas from double to real*8, new 
      !compiler doesn't allow real*16 to real*8.
      real*8 :: ergs,ergas 
      
      !Declarations added due to IMPLICIT NONE
      double precision :: ergnl,ergc,ergg,erga,ergest
      double precision :: ergtot
      real :: deltaphi,finish,phidown,zfieldown,zfieldup
      real :: phiup,enface
      integer :: i,ix,iy,iz,iisitpot
      real*4  timediff, timearray(2), dtime

      if (inrgwrt) open(42,file=trim(nrgnam))

      if (loga) then
         erga=0.0
         write(6,*) 'analytic grid energy is no longer available'
         stop
         if(inrgwrt) write(42,*) 'analytic grid energy is',erga,' kt'
      end if

      if (logg) then
         ergg=0.0

         !2011-06-11  Using operations on coord and int_coord type 
         !variables defined in module operators_on_coordinates 
         lim%min=2+bufz%min ; lim%max=igrid-1-bufz%max 

         do i=1,icount1b
            ixyz=gchrgp(i)
            ido=.true.
            if((ixyz.vorlt.lim%min).or.(ixyz.vorgt.lim%max)) ido=.false.
            if(ido) ergg=ergg + phimap(ixyz%i,ixyz%j,ixyz%k)*gchrg(i)
         end do

         ergg=ergg/2.0

         if(gaussian.eq.1.and.inhomo.eq.1.and.logs) then
            ergsgaussian=ergg
         else
            if(ideveloper) then
               write(6,*) 'total grid energy          :     ',ergg,' kt'
            else
               write(6,"(a,f20.4,a)") ' total grid energy               :',ergg,' kt'
            endif
         endif

         if(inrgwrt) write(42,*)'total grid energy: ',ergg,' kt'

      end if
        
      if(gaussian.eq.1.and.inhomo.eq.0.and.logs)then
        write(6,"(a,f20.4,a)"),        ' corrected reaction field energy :',ergg-ergsgaussian,' kt'
      endif
        
      if(gaussian.eq.1) then !LinLi skip other energy terms
         goto 1212
      endif

      if (logg.and.loga) then
         write(6,*) 'difference energy, in kt, is',(ergg-erga)
         write(6,*) 'difference energy, in kcals, is',(ergg-erga)*0.6
      end if

      if (ibctyp.eq.5) then
         write(6,*) "WARNING!!!Not completely tested routine &
                                            &for polarized membrane!!"

         if (lognl.or.logas) then
            write(6,*) "This option is not yet working with&
                                        & fixed potential difference!"
         ergas=0.0; ergnl=0.0
         end if

         !calcolo contributo facce "condensatore" assunte facce zeta!!
         !uso formula 0.5*(|q| media sulle due facce)*deltaV
         deltaphi=0.0; zfieldown=0.;  zfieldup=0.
         open(52,file='fields.txt',form='formatted')
         do ix=2,igrid-1
            do iy=2,igrid-1
               phiup    = phimap(ix,iy,igrid)
               phidown  = phimap(ix,iy,1)
               deltaphi = deltaphi+ (phiup-phimap(ix,iy,igrid-1))
               zfieldup=zfieldup-(phiup-phimap(ix,iy,igrid-1))*scale

               write(52,*)ix,iy,phimap(ix,iy,igrid-1),&
                                                &phimap(ix,iy,igrid-2)

               !averaging over upper and lower cube faces in order to 
               !cancel numerical errors
               deltaphi = deltaphi-(phidown-phimap(ix,iy,2))
               zfieldown=zfieldown+(phidown-phimap(ix,iy,2))*scale
            end do
         end do

         enface=deltaphi*vdrop%z*epsout*((igrid-1.)/(igrid-2.))**2/&
                                                       &(4.*fpi*scale)

         !ci sarebbe *((igrid-1)/(igrid-4))**2 ma ininfluente
         write(6,*)"Energy contribution from voltage drop=",enface,"kt"
         close (52)

         open(52,file='potcen.txt',form='formatted')
         write(6,*)"fieldup medio: ",zfieldup/(igrid-2.)**2
         write(6,*)"fieldown medio: ",zfieldown/(igrid-2.)**2

         do iz=1,igrid
            write(52,*)iz,phimap((igrid+1)/2,(igrid+1)/2,iz)
         end do
         close (52)
      end if

      if (irea.or.logs.or.lognl.or.logas.or.isen.or.isch) then
         !ergest=interaction energy of the solvent and the fixed charges
	      ergs=0.0; ergas=0.0; ergnl=0.0; ergest=0.0

         if (diff) ibc=0
         iisitpot=0
         if (isitpot) iisitpot=1

         !2011-06-11 All other parameters are accessible via qlog 
         !and pointers modules
         call react(ergs,ergas,iisitpot)
      end if

1212  continue

      if (logc.and.(.not.logions.or..not.lognl)) then
         ergc=0.0
            
         if (logions) then
            !linear case
            ergest=0.0

            call clbtot(ergest,ergc)

            write(6,*)'solvent contribution to fixed charges'
            write(6,*)'respectively inside and outside the cube :',&
                      &ergest,'kt',ergestout,'kt'
            write(6,*)'total ionic direct contribution :',&
                                                &ergest+ergestout,'kt'
         else
            if (nmedia.eq.1) then
               call clb(ergc)
               ergc=ergc/epsin
            else
               call clbmedia(ergc)
            end if
         end if

         if(ideveloper) then
            write(6,*) 'coulombic energy                :',ergc,' kt'
         else
            write(6,"(a,f20.4,a)") ' coulombic energy                :',ergc,' kt'
         end if
         if(inrgwrt) write(42,*) 'total coulombic energy:',ergc,' kt'
      end if

      if (lognl) then
         call nlener(ergnl,igridout)
         ergc=0.0; ergest=0.0

         if (logions) then
            call clbnonl(ergc,ergest,igridout)

            write(6,*) 'direct ionic contrib. inside the box:',&
                       &ergest,' kt'
            write(6,*) 'coulombic energy:                     ',&
                       &ergc,' kt'
            if(inrgwrt) write(42,*) 'total coulombic energy:',ergc,&
                                    &' kt'
         end if
      end if
      if(gaussian.eq.1)then
         goto 1213
         print *,"gaussian==1"
      endif
      if (logs.and.logions) then
         write(6,*)'Energy arising from solvent and boundary pol.',&
                   &ergnl+ergs+ergest+ergestout,' kt'
      end if

      if (lognl.and.logg) then
         if(ideveloper) then
            write(6,*)'total non linear grid energy    :',&
               &ergg+ergnl,' kt'
         else
            write(6,"(a,f20.4,a)")&
               &' total non linear grid energy    :',ergg+ergnl,' kt'
         end if
      end if

      ergtot=ergnl+ergc+ergs+ergest+ergestout

      if (logs.or.logc) then
         if(ideveloper) then
            write(6,*)'all required energy terms but grid and self_react.:',ergtot,'kt'
         else
            write(6,"(a,f20.4,a)")&
               &' all required energy terms but grid and self_react.:',ergtot,'kt'
         end if
         if(inrgwrt) write(42,*)'total required energy (everything calculated but grid and self_reaction energies: ',ergtot,'kt'
      end if

      if (logas.and.loga.and.logg) then
         write(6,*) "excess grid energy= ",ergg-ergas-erga
      end if

      call date_and_time(VALUES=values)
      finish=values(5)*3600.+values(6)*60.+values(7)+values(8)/1000.
      timediff = dtime( timearray)
      timetot=timetot+timediff

      if(ideveloper) then
         if(verbose)write(6,*)"energy calculations done at",timetot," s"
      else
         if(verbose)write(6,"(a,f10.2,a)")" energy calculations done at",timetot," s"
      end if

1213  continue
      if(inrgwrt) close(42)
      end subroutine encalc
