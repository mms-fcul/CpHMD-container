      !  ATTENTION!  This file is part of extrmmod module.
      !==================================================================
      ! 2011-05-10 Created by PK from file wrtprm.f
      subroutine wrtprm(cmin,cmax,cran)
      ! 2011-05-10  Include statement removed due to module architecture
      character(6) :: head
      character(12), dimension(5) :: bclab=&
      &(/ 'zero        ','dipolar     ',&
      &'focussing   ','coulombic   ','Voltage Drop' /)
      character(8) :: testchar
      !! b+++++++++++++++++++++++
      ! 2011-05-10 medeps declared in pointers module and allocated
      ! 2011-05-10 Changed to coord type variable
      logical :: is_lt_one
      type(coord) :: cmin,cmax,cran
      !------------------------------------------------------------------------------ 
      ! 2011-06-14 Declarations added due to IMPLICIT NONE
      integer :: i
      !------------------------------------------------------------------------------ 
      !2011-05-10  Data assignment made in declaration statement
      write(6,*)'  '      
      write(6,*)'grid size                  :',igrid
      if(ideveloper) then
         write(6,*)'percent of box to be filled:',perfil

      write(6,*)'scale,in grids/A :',scale
      ! 2011-05-10 Changed to coord type variable
      write(6,*)'xmin,xmax     (A):',cmin%x,cmax%x
      write(6,*)'ymin,ymax     (A):',cmin%y,cmax%y
      write(6,*)'zmin,zmax     (A):',cmin%z,cmax%z
      write(6,*)'x,y,z range   (A): ',cran
      write(6,*)'system geometric center in (A):',pmid
      write(6,*)' grid  box is centered in (A) :',oldmid
      write(6,*)'object centre offset (gu)  :',offset

      if(isolv)then
         write(6,*)'outer dielectric           :',repsout
         !! b+++++++++++++++
         do i = 0,nmedia
            write(6,*)'dielectric in medium number',i,' :',medeps(i)*epkt
            is_lt_one = (medeps(i)*epkt.lt.1.0)
         end do
         if (is_lt_one) then
            write(6,*) 'DelPhi is not able to deal with eps < 1.0'
            write(6,*) 'Therefore Exiting' 
            stop
         end if
         
         write(6,*)'first kind salt concentration (M)   :',conc(1)
         write(6,*)'valences salt 1 are        ',ival(1),' and',ival(2)
         write(6,*)'second kind salt concentration (M)  :',conc(2)
         write(6,*)'valences salt 2 are        ',ival2(1),' and',ival2(2)
         !! e+++++++++++++++
         write(6,*)'ionic strength (M)         :',rionst
         write(6,*)'debye length (A)           :',deblen
         write(6,*)'absolute temperature (K)   :',temperature
         write(6,*)'ion exclusion radius (A)   :',exrad
         write(6,*)'probe radius facing water(A:',radprb(1)
         write(6,*)'probe radius, internal (A) :',radprb(2)
         write(6,*)'boundary conditions        : ',bclab(ibctyp)
         write(6,*)'x,y,z periodic bc. and volt. drop flags:',iper
         if(iper(4).or.iper(5).or.iper(6)) then
            ! 2011-05-10   Changed to coord type variable
            write(6,*)"Voltage drops along x,y,z :",vdrop
         end if
         if(iautocon) then
            if(gten.gt.0.)then
               write(6,*)'convergence by grid energy :',gten,' kt'
            else
               write(6,*)'# of linear iterations     : automatic'
            endif
         else
            write(6,*)'# of linear iterations     :',nlit
         end if
         !! b+++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if (res1.gt.0..or.res2.gt.0..or.res5.gt.0.) then
            write(6,*)'convergence by rms  change :',res1,' kt'
            write(6,*)'convergence by max  change :',res2,' kt'
            write(6,*)'convergence by norm change :',res5,' kt'
         end if




         if (rionst.lt.1.e-6.and.nnit.gt.0) then
            write(6,*)'ionic strength=0 ==> only linear iterations'
         else
            write(6,*)'# of non-linear iterations :',nnit
            write(6,*)'non-linear energy calculat.:',lognl
            write(6,*)'manual relaxation parameter :',imanual
         end if
         !! e+++++++++++++++++++++++++++++++++++++++++++++++++++++++
         write(6,*)'ionic direct energy contribution:',logions
         write(6,*)'concentration map output   :',iconc
         write(6,*)'spherical charge distbn.   :',isph
         write(6,*)'INSIGHT format output      :',ibios
         write(6,*)'site potential output      :',isite
      endif



      else
         write(6,"(a,f12.2)")' percent of box to be filled:',perfil

      write(6,"(a,f7.2)")' scale,in grids/A :',scale
      ! 2011-05-10 Changed to coord type variable
      write(6,"(a,2f10.3)")' xmin,xmax     (A):',cmin%x,cmax%x
      write(6,"(a,2f10.3)")' ymin,ymax     (A):',cmin%y,cmax%y
      write(6,"(a,2f10.3)")' zmin,zmax     (A):',cmin%z,cmax%z
      write(6,"(a,3f10.3)")' x,y,z range   (A):',cran
      write(6,"(a,3f10.3)")' system geometric center in (A):',pmid
      write(6,"(a,3f10.3)")' grid  box is centered in (A)  :',oldmid
      write(6,"(a,3f10.3)")' object centre offset (gu)     :',offset

      if(isolv)then
         write(6,"(a,f7.2)")' outer dielectric              :',repsout
         !! b+++++++++++++++
         do i = 0,nmedia
            write(6,"(a,i3,a,f7.2)")' dielectric in medium number',i,':',medeps(i)*epkt
            is_lt_one = (medeps(i)*epkt.lt.1.0)
         end do
         if (is_lt_one) then
            write(6,*) 'DelPhi is not able to deal with eps < 1.0'
            write(6,*) 'Therefore Exiting' 
            stop
         end if
         
         write(6,"(a,f7.4)")' first kind salt concentration (M)       :',conc(1)
         write(6,*)'valences salt 1 are        ',ival(1),' and',ival(2)
         write(6,"(a,f7.4)")' second kind salt concentration (M)      :',conc(2)
         write(6,*)'valences salt 2 are        ',ival2(1),' and',ival2(2)
         !! e+++++++++++++++
         write(6,"(a,f7.4)")' ionic strength (M)            :',rionst
         write(6,"(a,f7.4)")' debye length (A)              :',deblen
         write(6,"(a,f7.2)")' absolute temperature (K)      :',temperature
         write(6,"(a,f7.4)")' ion exclusion radius (A)      :',exrad
         write(6,"(a,f7.4)")' probe radius facing water(A)  :',radprb(1)
         write(6,"(a,f7.4)")' probe radius, internal (A)    :',radprb(2)
         write(6,*)'boundary conditions           : ',bclab(ibctyp)
         write(6,*)'x,y,z periodic bc. and volt. drop flags:',iper
         if(iper(4).or.iper(5).or.iper(6)) then
            ! 2011-05-10   Changed to coord type variable
            write(6,*)"Voltage drops along x,y,z :",vdrop
         end if
         if(iautocon) then
            if(gten.gt.0.)then
               write(6,*)'convergence by grid energy :',gten,' kt'
            else
               write(6,*)'# of linear iterations        : automatic'
            endif
         else
            write(6,*)'# of linear iterations     :',nlit
         end if
         !! b+++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if (res1.gt.0..or.res2.gt.0..or.res5.gt.0.) then
            write(6,"(a,f10.6,a)")' convergence by rms  change    :',res1,' kt'
            write(6,"(a,f10.6,a)")' convergence by max  change    :',res2,' kt'
            write(6,"(a,f10.6,a)")' convergence by norm change    :',res5,' kt'
         end if

         if (rionst.lt.1.e-6.and.nnit.gt.0) then
            write(6,*)'ionic strength=0 ==> only linear iterations'
         else
            write(6,"(a,i10)")' # of non-linear iterations    :',nnit
            write(6,*)'non-linear energy calculat.   :',lognl
            write(6,*)'manual relaxation parameter   :',imanual
         end if
         !! e+++++++++++++++++++++++++++++++++++++++++++++++++++++++
         write(6,*)'ionic direct energy contribution:',logions
         write(6,*)'concentration map output      :',iconc
         write(6,*)'spherical charge distbn.      :',isph
         write(6,*)'INSIGHT format output         :',ibios
         write(6,*)'site potential output         :',isite
      endif


      end if


      write(6,*)'modified atom file output     :',iatout
      write(6,*)'map file label                :'
      write(6,'(a)')toplbl
      
      if(ipdbrd.or.ipdbwrt) then
         if(ipdbrd) write(6,*) 'set to  read unformatted pdb file'
         if(ipdbrd.and.ipdbwrt) then
            write(6,*) ' '
            write(6,*) ' WARNING: can not write an unformatted pdb'
            write(6,*) ' file, while reading in an unformatted pdb file'
            write(6,*) ' Therefore the write option is disabled'
         else
            if(ipdbwrt) write(6,*) ' set to write unformatted pdb file'
         end if
      end if



      
      if(ifrcrd.or.ifrcwrt) then
         if(ifrcrd.and.isite) write(6,*) 'set to read unformatted frc.pdb &
         &   file'
         if(ifrcrd.and.ifrcwrt) then
            write(6,*) ' '
            write(6,*) ' WARNING: can not write an unformatted frc'
            write(6,*) ' file, while reading in an unformatted frc file'
            write(6,*) ' Therefore the write option is disabled'
         else
            if(ifrcwrt) write(6,*) ' set to write unformatted frc.pdb file'
         end if
      end if
      
      if(.not.igraph) write(6,*) ' convergence graph turned off'
      if(.not.ipoten) write(6,*) ' potential listings turned off'
      if((icon1.ne.10).or.(icon2.ne.1)) then
         write(6,*) 'convergence test interval is every',icon1,'loops'
         write(6,*) 'testing',100/icon2,'%'
      end if
      write(6,*)' '
      end subroutine wrtprm
	  
      !---------------------------------------------------------------------------
      subroutine wrt(ipgn)
      integer :: ipgn
      !!c subroutine to handle log file writes
      if(ipgn.eq.1) then
         write(6,*)'   '
         write(6,*)' ______________DelPhi V. 6.1 (f95)_______________  '
         write(6,*)'|                                                | '
         write(6,*)'| A program to solve the PB equation             | '
         write(6,*)'| in 3D, using non-linear form, incorporating    | '
         write(6,*)'| many dielectric regions, multisalt ionic       | '
         write(6,*)'| strength, different probe radii, periodic      | '
         write(6,*)'| and focussing boundary conditions, utilizing   | '
         write(6,*)'| stripped optimum successive over-relaxation    | '
         write(6,*)'| and an improved algorithm for mapping the      | '
         write(6,*)'| Mol. Surface to the finite-Difference grid     | '
         write(6,*)'|                                                | '
         write(6,*)'|    If there is any question, please go to:     | '
         write(6,*)'|       http://compbio.clemson.edu/forum/        | '
         write(6,*)'|     June 2012,by DelPhi Development Team       | '
         write(6,*)'|_________________             __________________| '
         write(6,*)'                  DelPhi V. 6.1                    '
         write(6,*)'  '
         
         call date_and_time(DATE=day,TIME=time,VALUES=values)
         write(6,*)'Program started on ',day(1:4),'-',day(5:6),'-',&
         &day(7:8),' at ',time(1:2),':',time(3:4),':',time(5:6)
      end if
      return
      end subroutine wrt
      
