!#####################################################################
!  ATTENTION!  This file is part of epsmakmod module.
!              Do not compile separately!
!
! 2011-05-10 Parameters transferred via modules qlog and pointers
!#####################################################################

      subroutine epsmak
      
      implicit none

      type(coord) :: amin,amax
      !2011-05-27 Declarations added due to IMPLICIT NONE
      real :: rmid,start,finish,rmax,cputime
      integer :: ii
      !2011-05-10 Seems bndeps is not used here, moved to vwtms 
      !           subroutine
      
      !2011-05-10 Array is allocated using conventional F95 statement
      !           limgunit is array of object_min_max derived type
      allocate(limgunit(Nobject))

      !Lin Li: ibmx value used to be 1000000, sometimes not enough
      ibmx=50000000 

      !here limobject is expressed in grid units
      rmid=real((igrid+1)/2)

      do ii=1,nobject
         !2011-05-10  Using operations on coord type variables defined
         !            in module operators_on_coordinates 
         limgunit(ii)%min=(limobject(ii)%min-oldmid)*scale+rmid
         limgunit(ii)%max=(limobject(ii)%max-oldmid)*scale+rmid
      end do

      !save file datadelphi, containing some variables, to be used by 
      !GUI
      !open(52,file='datadelphi',form='unformatted')
      !write(52)igrid
      !write(6,*)limgunit
      !write(6,*)oldmid
      !write(6,*)scale
      !write(52)nobject
      !write(6,*)limobject
      !write(52)dataobject
      !close(52)
      !
      !write(6,*)'limobject'
      !write(6,*)limobject
      !write(6,*)'limgunit'
      !write(6,*)limgunit
      
      if (uniformdiel) then
         write(6,*) "not going to calculate boundary elements since"
         write(6,*) "uniform dielectric"
         ibnum=0
         return
      end if

      !lepsx.y.z and uepsx.y.z should be the upper and lower limits of 
      !the expanded box. if the molecule is smaller than this then 
      !reduce leps and upeps accordingly
      !note leps/ueps not yet defined..

      !2011-05-10 Converted to coord derived type
      amin=limgunit(1)%min
      amax=limgunit(1)%max

      !find global limits IN GRID UNITS, both, molecule and objects, 
      !are considered
      if (nobject.gt.1) then
         do ii=2,nobject
            !2011-05-10  Using operations on coord type variables 
            !defined in module operators_on_coordinates 
            amin=min(amin,limgunit(ii)%min)
            amax=max(amax,limgunit(ii)%max)
         end do
      end if

      rmax=rdmx
      
      !do i=1,natom
      !   xmin=min(xn2(1,i),xmin) etc...
      !   xmax=max(xn2(1,i),xmax)etc...
      !   rmax=max(rad3(i),rmax)
      !end do	

      if(rionst.ne.0.)rmax=max(rmax,exrad)
      rmax=rmax*scale

      !2011-05-10  Using operations on coord type variables defined
      !in module operators_on_coordinates 
      amin=amin-rmax
      amax=amax+rmax

      !2011-05-10  Using operations on coord and int_coord 
      !type variables defined in module operators_on_coordinates 
      limeps%min=max(int(amin)-2,1)
      limeps%max=min(int(amax)+2,igrid)

      !2011-05-10  Changed to array operations
      idebmap=.true.
      iepsmp=int_coord(0,0,0)
      !point is out of any kind of object (probably in solution)
      !iepsmp(i,j,k,1)=0 no longer needed, calloc superseded malloc
      !iepsmp(i,j,k,2)=0 no longer needed, calloc superseded malloc
      !iepsmp(i,j,k,3)=0 no longer needed, calloc superseded malloc
      !idebmap(i,j,k)=.true.

      !hgrl=1./2./scale
      !if radprb is less than half of grid spacing, then use old 
      !algorithm (sri 29March 93)
      !The new algorithm should be able to handle all scales; still 
      !remains to be tested (Sri Apr 95)

      call date_and_time(DATE=day,TIME=time,VALUES=values)
      start=values(5)*3600.+values(6)*60.+values(7)+values(8)/1000.
      if(verbose)write(6,*) 'start vw surface at ' ,&
                                &time(1:2),':',time(3:4),':',time(5:6)

      if(gaussian.eq.0) then
         call setout(0.)
      else
         call setgaussian(0.)
      endif 


      call date_and_time(DATE=day,TIME=time,VALUES=values)
      finish=values(5)*3600.+values(6)*60.+values(7)+values(8)/1000.
      if(verbose)write(6,*) 'fill in re-entrant regions at ',&
                                &time(1:2),':',time(3:4),':',time(5:6)


      if(ideveloper) then
         if(verbose)write(6,*) 'time elapsed :',finish-start,' s' 
      else
         if(verbose)write(6,"(a,f10.2,a)") ' time elapsed :',finish-start,' s' 
      end if

      !2011-05-10 Array is deallocated with ordinary F95 statement
      if(allocated(limgunit)) deallocate(limgunit)

      !2011-05-17 Parameters transfered via module architecture

      if(gaussian.eq.0) then
          !print *,"in vwtms,gaussian:",gaussian
          call vwtms
      else
      
      endif

      call date_and_time(DATE=day,TIME=time,VALUES=values)
      finish=values(5)*3600.+values(6)*60.+values(7)+values(8)/1000.

      if(ideveloper) then
         if(.not.isrf.and.verbose) write(6,*)&
           & 'time to turn everything in is',finish-start,' s'
      else
         if(.not.isrf.and.verbose) write(6,"(a,f10.2,a)")&
           & ' time to turn everything in is',finish-start,' s'
      end if

      !comment out membrane stuff for a moment..
      !if(imem) call mkmem

      end subroutine epsmak
