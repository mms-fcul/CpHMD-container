      !  ATTENTION!  This file is part of setbcmod module.
      !==================================================================
      ! 2011-06-07 All parameters are accessible via qlog and pointers
      subroutine itit
      
      !!c NOTE THIS HAS BEEN ALTERED TO GIVE THE CONVERGENCE TO THE EXACT
      !!C SOLUTION, WHICH IS INPUTTED VIA OPTION 3,(AND NOW WITH GAUSS
      !!C SIEDEL IMPLEMENTED, WITH MAPPED SUBARRAYS.)
      
      !!c  finally what we`ve all been waiting for-
      !!c do the actual pb equation iteration
      !!c first the linear, with relaxation parameter 1.
      !!c then non-linear with relaxation parameter 0.6
      
      !!c this is a modified iteration routine which runs at about
      !!c three times that of the original version. this has been 
      !!c accomplished by several features which will be mentioned
      !!c at the appropiate juncture.
      !!c note that the arrays slab,old and denom no longer exist
      !!c and that epsmap is only used in setting up arrays
      !!c in their place are several arrays, the two biggest being
      !!c phimap2 and sf.
      !!c however, the array space accessed in the main inner loop
      !!c is no larger than before, i.e. memory requirements should
      !!c not be much different.
      
      !!c some notes on array sizes. there can be no more than 10000
      !!c charges, or 12000 boundary elements, or 40000 high potential salt
      !!c sites (for the nonlinear). also, for the latter the potential
      !!c in salt should not exceed 10kt
      
      
      integer, parameter :: nxran = 60, nyran = 60
      real :: rmsl(nxran),rmsn(nxran),rmaxl(nxran),rmaxn(nxran)
      
      ! 2011-06-07 Arrays and variables are accessible via qlog and
      
      character(24) :: day
      character :: symb
      character(70) :: title,iplot(nyran)
      logical :: qstopper,resdat,once,igt
      integer :: star,fin,sta1(igrid),sta2(igrid),fi1(igrid),fi2(igrid)
      real :: Green(-10:10,-10:10,-10:10),grdn(5)
      integer :: n,m,l
      !----------------------------------------------------------------------
      ! 2011-06-07 Declarations added due to IMPLICIT NONE
      real :: start,ap1,ap2,ap3,ap4,finish,grden,epsrat
      real :: om1,om2,om3,om4,res3,res4,rmsch2,rnorm2,rmxch2
      real :: rmxch,rtemp2,rtemp3,sixth,rnormch,sumdown,sumup
      real :: temp,temp1,temp2,temp3,temp4,tmp,rmsch,ymin,ymax
      real*4  timediff, timearray(2), dtime 
      
      integer :: icount2abac,icount2bbac,iadd1,iadd2,ibin
      integer :: iclr,iscl,imk,idif1x,idif2x,idif1y,idif2y
      integer :: idif1z,idif2z,icgrid,i,ii,ihgd,ihgd2,ires
      integer :: inc1za,inc1zb,inc1ya,inc1yb,inc1xa,inc1xb
      integer :: inc2za,inc2zb,inc2ya,inc2yb,inc2xa,inc2xb
      integer :: iplt,j,k,ix,iy,iz,ishift,itemp,itemp1,itemp2
      integer :: itemp3,itemp4,itest,iw,lat1,lat2,long1,long2
      integer :: npoint,nn,midg,it
      !!c-------------------------------------------------------------
      call date_and_time(DATE=day,TIME=time,VALUES=values)
      start=values(5)*3600.+values(6)*60.+values(7)+values(8)/1000.
      timediff = dtime( timearray)
      timetot=timetot+timediff
      if(verbose) write(6,*)'now iterating on ',day(1:4),'-',day(5:6),'-',& !Lin Li: add date
      &day(7:8),' at ',time(1:2),':',time(3:4),':',time(5:6)
      
      !!c allocate memory to arrays
      !!c some initialization
      
      once=.true.
      icount2abac=icount2a; icount2bbac=icount2b; epsrat=epsout/epsin 
      sixth = 1.d0/6.d0
      !!c         th120 = 1./120.
      icgrid=igrid*igrid*igrid;  ihgd=(igrid+1)/2
      if(icon2.eq.0) then
         icon1=10; icon2=1
      end if
      if(icon1.gt.nlit) icon1=nlit 
      ! 2011-06-07 Changed to array operations
      rmsl=0.; rmsn=0. ; rmaxl=0. ; rmaxn=0.
      npoint = (igrid-2)*(igrid-2)*(igrid-2)
      
      
      !!c ---------------------------------------------
      !!c MAIN SET UP ROUTINE	
      !!c ---------------------------------------------
      !!c     OLD WAY
      if(iper(1)) then
         n=0
         do iz=2,igrid-1
            iadd1=(iz-1)*igrid*igrid 
            do iy=2,igrid-1
               iadd2=(iadd1+(iy-1)*igrid +2)/2
               n=n+1; ibndx(n)=iadd2
            end do
         end do
         idif1x=(igrid-2)/2; idif2x=idif1x+1
         inc1xa=1; inc1xb=0; inc2xa=0; inc2xb=1
      end if
      !-------------------------------------------------------------
      if(iper(2)) then
         n=0
         do iz=2,igrid-1
            iadd1=(iz-1)*igrid*igrid 
            do ix=2,igrid-1
               iadd2=(iadd1+ix+1)/2
               n=n+1; ibndy(n)=iadd2
            end do
         end do
         idif1y=igrid*(igrid-2)/2; idif2y=idif1y+1
         inc1ya=(igrid/2)+1; inc1yb=inc1ya-1
         inc2ya=inc1yb;      inc2yb=inc1ya
      end if
      !-------------------------------------------------------------
      if(iper(3).or.iper(6)) then
         n=0
         do ix=2,igrid-1
            iadd1=ix+1
            do iy=2,igrid-1
               iadd2=(iadd1+(iy-1)*igrid)/2
               n=n+1;  ibndz(n)=iadd2
            end do
         end do
         idif1z=igrid*igrid*(igrid-2)/2; idif2z=idif1z+1
         inc1za=((igrid**2)/2)+1; inc1zb=inc1za
         inc2za=inc1zb;           inc2zb=inc1za
      end if
      !-------------------------------------------------------------
      !!c END OF SET UP	
      
      !!c remove qstopper file if it already exists
      inquire(file='qstop.test',exist=qstopper)
      if(qstopper) then
         open(30,file='qstop.test')
         close(30,status='delete')
         qstopper=.false.
      end if
      !!c check for resolution data
      !!c commented out, by kas, 29my 89, as syntax not vax compatible, and obseleted
      
      resdat = .false.
      !!c b++commented out by Walter (06-2001), convergence criteria given by user
      !!c       res1=0.0
      !!c       res2=0.0
      !!c e++++++++++++++++++++++++++++++++++++++++++++++++++++
      res3=0.0;  res4=0.0
      !!c     end if
      !-------------------------------------------------------------
      if(resdat) then 
         write(6,*) ' '
         write(6,*) 'linear resolution criteria are:',res1,res2,res5
         if(nnit.ne.0) then
            write(6,*) 'non-linear resolution criteria are:',res3,res4
         end if
         write(6,*) ' '
      end if
      
      if(verbose) then
         write(6,*) ' '; write(6,*) ' '
         if(gten.gt.0.0)then
            write(6,*) &
            &'  rms-change     max change    &
            &grid energy    #iterations'
         else
            write(6,*) '  rms-change     max change    &
            &   #iterations'
         endif
      end if
      !-------------------------------------------------------------
      !!c set up start and stop vectors
      sta1(2)=(igrid*igrid + igrid +4)/2; sta2(2)=sta1(2)-1
      fi1(2)=igrid*igrid - (igrid+1)/2;   fi2(2)=fi1(2)
      itemp1=igrid + 2;     itemp2=igrid*igrid -igrid -2
      do i=3,igrid-1
         sta1(i)=fi1(i-1) + itemp1; sta2(i)=fi2(i-1) + itemp1
         fi1(i)=sta1(i-1) + itemp2; fi2(i)=sta2(i-1) + itemp2
      end do
      !!c also
      lat1= (igrid-1)/2; lat2= (igrid+1)/2
      long1= (igrid*igrid - 1)/2; long2= (igrid*igrid + 1)/2
      
      ires=0
      !-------------------------------------------------------------
      if(icheb) then
         om2=1.0
         ! 2011-06-07 To get rid of GOTO statement
      else
         om2=2.0/(1.0 + sqrt(1 - spec))
         do ix=1,(icgrid+1)/2
            sf1(ix)=sf1(ix)*om2; sf2(ix)=sf2(ix)*om2
         end do
         do ix=1,icount1b 
            qval(ix)=qval(ix)*om2
         end do
         do iy=1,6
            do ix=1,icount2b
               db(iy,ix)=db(iy,ix)*om2
            end do
         end do
         sixth=sixth*om2 
      end if
      !-------------------------------------------------------------
      
      om1=1.0-om2; i=1;  iw=1
      !debug tested           write(6,*)'ITIT2--->',om2,om1
      do iz=1,igrid
         do iy=1,igrid
            do ix=1,igrid
               phimap3(iw)=phimap(ix,iy,iz)
               iw=iw+1
            end do
         end do
      end do
      !-------------------------------------------------------------
      ! 2011-06-07 DO-EXIT construct instead of IF-GOTO
      D1689: do
         do ix=1,(icgrid+1)/2
            iy=ix*2; phimap1(ix)=phimap3(iy-1)
            phimap2(ix)=phimap3(iy)
         end do
         ihgd2=ihgd-1
         !!c b+++++++++++++++++++++++++++++++++++++++
         !!c        if (.not.iper(1)) then THIS IS A PART OF A DIFFERENT pbc ASSIGNMENT
         !!c e+++++++++++++++++++++++++++++++++++++++
         !-------------------------------------------------------------
         !!c       set x boundary values 
         star=(igrid+1)/2; iy=(igrid*(igrid+1)/2)-igrid+1
         fin=(igrid*(igrid-1)-2)/2
         do ix=star,fin
            iy=iy+igrid
            bndx1(ix)=phimap1(iy)
            bndx2(ix)=phimap1(iy+ihgd2)
         end do
         
         star=(igrid+2)/2
         iy=(igrid*(igrid+2)/2) - igrid +1
         fin=(igrid*(igrid-1)-1)/2
         do ix=star,fin
            iy=iy+igrid
            bndx3(ix)=phimap2(iy)
            bndx4(ix)=phimap2(iy+ihgd2)
         end do
         !-------------------------------------------------------------
         !!c        end if THIS IS A PART OF A DIFFERENT pbc ASSIGNMENT
         !-------------------------------------------------------------
         ! 2011-06-07 DO - EXIT construct in iteration loop instead of IF - GOTO
         D1000: do
            !!c clear rms, max change
            rmsch = 0.0; rmxch = 0.00
            !!c if there is no salt then the main loop is executed without sf
            !!c saving about 15% in execution time
            !-------------------------------------------------------------
            if(rionst.gt.0.0) then
               do n = 2, igrid-1
                  star=sta1(n); fin=fi1(n)
                  do ix = star,fin
                     temp1 = phimap2(ix) + phimap2(ix-1)
                     temp2 = phimap2(ix+lat1) + phimap2(ix-lat2)
                     temp3 = phimap2(ix+long1) + phimap2(ix-long2)
                     phimap1(ix)=phimap1(ix)*om1+(temp1+&
                     & temp2+temp3)*sf1(ix)
                  end do
               end do
               !!c otherwise the main loop is as below:
            else
               do n = 2, igrid-1
                  star=sta1(n); fin=fi1(n)
                  do ix = star,fin
                     temp1 = phimap2(ix) + phimap2(ix-1)
                     temp2 = phimap2(ix+lat1) + phimap2(ix-lat2)
                     temp3 = phimap2(ix+long1) + phimap2(ix-long2)
                     ! debug tested               write(6,*) '---1000-cycle--->',n,ix,temp1,temp2,temp3
                     phimap1(ix) = phimap1(ix)*om1 +&
                     &(temp1+temp2+temp3)*sixth
                  end do
               end do
            end if
            !-------------------------------------------------------------
            !!c the above loops are about fourtimes faster than the original
            !!c loop over all grid points for several reasons, the biggest being that
            !!c we are only solving laplace's equation (unless salt is present), which
            !!c numerically much simpler, hence faster. we put all we leave out, back
            !!c in below, ending up with an equivalent calculation, but much faster.
            
            !!c b+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if (iper(1).and..false.) then
               !!c calculating first slice
               
               ix=1+igrid*lat2-igrid; ishift=igrid*lat1
               !not coming here        write(6,*) '---1000-cycle--->',ix,lat2,lat1,igrid,ishift
               if(rionst.gt.0.0) then
                  D2000: do n=1,(igrid-3)*lat2+1
                     ix=ix+igrid
                     !!c scandisco tutti i punti, anche i bordi (in questa fase del 
                     !!c ciclo ho solo spigoli xy) e li scarto a posteriori
                     itest=mod(n,igrid)
                     if(itest.eq.0..or.itest.eq.lat2) then
                        if (iper(2)) then
                           if(itest.eq.0.)then
                              temp2 = phimap2(ix+lat1-ishift)+&
                              &phimap2(ix-lat2)
                           else
                              temp2 = phimap2(ix+lat1)+&
                              &phimap2(ix-lat2+ishift)
                           end if
                        else 
                           cycle D2000
                        end if
                     else
                        temp2 = phimap1(ix-lat1)  + phimap1(ix+lat2)
                     end if
                     !!c      secondo template x+ x- (al contrario della fase successiva)
                     temp1 = phimap2(ix)      + phimap2(ix-1+lat1)
                     temp3 = phimap2(ix+long1)+ phimap2(ix-long2)
                     
                     rtemp2=phimap2(ix+lat1+lat1)+phimap2(ix-lat2+lat1)
                     rtemp3=phimap2(ix+long1+lat1)+phimap2(ix-long2+lat1)
                     
                     phimap1(ix)=.5*(phimap1(ix+lat1)+&
                     & phimap1(ix))*om1+(temp1+ .5*(temp2+temp3+&
                     & rtemp2+rtemp3))*(sf1(ix)+sf1(ix+lat1))*.5
                     !!c now updating last slice
                     phimap1(ix+lat1)=phimap1(ix)
                  end do D2000
               else
                  D2100: do n=1,(igrid-3)*lat2+1
                     ix=ix+igrid
                     !!c scandisco tutti i punti, anche i bordi (in questa fase del ciclo 
                     !!c ho solo spigoli xy) e li scarto a posteriori
                     itest=mod(n,igrid)
                     if(itest.eq.0..or.itest.eq.lat2) then
                        if (iper(2)) then
                           if(itest.eq.0.)then
                              temp2 = phimap2(ix+lat1-ishift)+&
                              &phimap2(ix-lat2)
                           else
                              temp2 = phimap2(ix+lat1)+&
                              &phimap2(ix-lat2+ishift)
                           end if
                        else 
                           cycle D2100
                        end if
                     else
                        temp2 = phimap1(ix-lat1)  + phimap1(ix+lat2)
                     end if
                     temp1 = phimap2(ix)+phimap2(ix-1+lat1)
                     temp3 = phimap2(ix+long1) + phimap2(ix-long2)
                     
                     phimap1(ix) = phimap1(ix)*om1 + &
                     & (temp1+temp2+temp3)*sixth
                     !!c now updating last slice
                     phimap1(ix+lat1)=phimap1(ix)
                  end do D2100
               end if
               
               star=(igrid+1)/2; fin=(igrid*(igrid-1)-2)/2
               iy=(igrid*(igrid+1)/2) - igrid + 1
               do ix=star,fin
                  iy=iy+igrid
                  bndx1(ix)=phimap1(iy)
                  bndx2(ix)=phimap1(iy+ihgd2)
               end do
            end if
            !!c e+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
            !!c first we add back the dielectric boundary points, by recalculating them
            !!c individually. note this is still vectorised by means of a gathering
            !!c load by the compiler.
            
            !!C$DIR NO_RECURRENCE 
            do k=1,icount2a
               ix=idpos(k)
               temp1=phimap2(ix-1)*db(1,k)+phimap2(ix)*db(2,k)
               temp2=phimap2(ix-lat2)*db(3,k)+phimap2(ix+lat1)*db(4,k)
               temp3=phimap2(ix-long2)*db(5,k)+phimap2(ix+long1)*db(6,k)
               phimap1(ix)= phimap1(ix) + temp1+temp2+temp3
            end do
            !!c Reset x boundary values altered in above loops.
            star=(igrid+1)/2; fin=(igrid*(igrid-1)-2)/2
            iy=(igrid*(igrid+1)/2) - igrid +1
            !!C$DIR NO_RECURRENCE
            do ix=star,fin
               iy=iy+igrid
               phimap1(iy)=bndx1(ix)
               phimap1(iy+ihgd2)=bndx2(ix)
               ! debug tested        write(6,*)'****1000-cycle****',ix,iy,ihgd2,bndx1(ix),bndx2(ix)
            end do
            !!c next we add back an adjustment to all the charged grid points due to
            !!c the charge assigned. the compiler directive just reassures the vector
            !!c compiler that all is well as far as recurrence is concerned, i.e. it
            !!c would think there is a recurrence below, where as in fact there is none.
            !!C$DIR NO_RECURRENCE 
            do k=1,icount1a
               temp=qval(k)
               ix=iqpos(k)
               phimap1(ix)=phimap1(ix) + temp
            end do
            !!c if periodic boundary condition option
            !!c force periodicity using wrap around update of boundary values:
            !!c 2nd slice-->last
            !!c last-1 slice-->first
            
            !----------------------------------------------------------------------
            !!c z periodicity
            if(iper(3)) then
               do iz = 1,(igrid-2)**2,2
                  temp1=ibndz(iz);    itemp1=temp1
                  temp2=temp1+idif1z; itemp2=temp2
                  temp3=temp2+inc1za; itemp3=temp3
                  temp4=temp1+inc1zb; itemp4=temp4
                  !!c           iz=1  ===       iz=64
                  phimap1(itemp1)=phimap2(itemp2)
                  !!c           iz=65   ===      iz=2
                  phimap1(itemp3)=phimap2(itemp4)
               end do
            end if
            !----------------------------------------------------------------------
            if(iper(6).and..false.) then
               sumdown=0.0 ;  sumup=0.0
               do iz = 1,(igrid-2)**2,2
                  temp1=ibndz(iz)
                  temp2=temp1+idif1z; itemp2=temp2
                  temp4=temp1+inc1zb; itemp4=temp4
                  sumup=sumup+phimap2(itemp2)
                  sumdown=sumdown+phimap2(itemp4)
               end do
               tmp=((igrid-2)**2+1)/2.
               sumup=sumup/tmp-vdrop%z
               sumdown=sumdown/tmp+vdrop%z
               do iz = 1,(igrid-2)**2,2
                  temp1=ibndz(iz); itemp1=temp1
                  temp2=temp1+idif1z
                  temp3=temp2+inc1za; itemp3=temp3
                  phimap1(itemp1)=sumup
                  phimap1(itemp3)=sumdown
                  !!c           write(6,*)sumup,sumdown
               end do
            end if
            !----------------------------------------------------------------------
            !!c y periodicity
            if(iper(2)) then
               do iy = 1,(igrid-2)**2,2
                  temp1=ibndy(iy);    itemp1=temp1
                  temp2=temp1+idif1y; itemp2=temp2
                  temp3=temp2+inc1ya; itemp3=temp3
                  temp4=temp1+inc1yb; itemp4=temp4
                  phimap1(itemp1)=phimap2(itemp2)
                  phimap1(itemp3)=phimap2(itemp4)
               end do
            end if
            !----------------------------------------------------------------------
            !!c x periodicity (old way)
            if(iper(1)) then
               do  ix = 1,(igrid-2)**2,2
                  temp1=ibndx(ix);    itemp1=temp1
                  temp2=temp1+idif1x; itemp2=temp2
                  temp3=temp2+inc1xa; itemp3=temp3
                  temp4=temp1+inc1xb; itemp4=temp4
                  phimap1(itemp1)=phimap2(itemp2)
                  phimap1(itemp3)=phimap2(itemp4)
               end do
            end if
            !----------------------------------------------------------------------
            if(icheb) then
               !----------------------- Begin body of omalt subroutine -------------
               it=2*i-1
               om3=1./(1.-om2*spec*0.25)
               if(om1.lt.1.e-6) om3=1./(1.-om2*spec*0.5)
               om4=om3/om2; om2=om3; om1=1.0-om2
               if(rionst.gt.0.0) then
                  if(mod(it,2).eq.1) then
                     do ix=1,nhgp
                        sf1(ix)=sf1(ix)*om4
                     end do
                  else
                     do ix=1,nhgp
                        sf2(ix)=sf2(ix)*om4
                     end do
                  end if
               end if
               
               do ix=1,icount1b
                  qval(ix)=qval(ix)*om4
               end do
               
               do iy=1,6
                  do ix=1,icount2b
                     db(iy,ix)=db(iy,ix)*om4
                  end do
               end do
               sixth=sixth*om4
            end if
            !-------------------------- End  body of omalt subroutine -------------
            !debug tested        write(6,*) '---omalt---',om1,sixth
            !----------------------------------------------------------------------
            !!c Next update phimap2 using the new phimap1
            if(rionst.gt.0.0) then	
               do n = 2, igrid-1
                  star=sta2(n); fin=fi2(n)
                  do ix = star,fin
                     temp1 = phimap1(ix)       + phimap1(ix+1)
                     temp2 = phimap1(ix+lat2)  + phimap1(ix-lat1)
                     temp3 = phimap1(ix+long2) + phimap1(ix-long1)
                     phimap2(ix) =phimap2(ix)*om1+&
                     & (temp1+temp2+temp3)*sf2(ix)
                  end do
               end do
               !----------------------------------------------------------------------
            else
               do n = 2, igrid-1
                  star=sta2(n); fin=fi2(n)
                  do ix = star,fin
                     temp1 = phimap1(ix)       + phimap1(ix+1)
                     temp2 = phimap1(ix+lat2)  + phimap1(ix-lat1)
                     temp3 = phimap1(ix+long2) + phimap1(ix-long1)
                     phimap2(ix) =phimap2(ix)*om1+&
                     &(temp1+temp2+temp3)*sixth
                  end do
               end do
            end if
            !!c b+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if (iper(1).and..false.) then
               !!c calculating first slice
               ix=long2-igrid; ishift=igrid*lat1
               if(rionst.gt.0.0) then
                  D2200: do n=0,long1-igrid
                     ix=ix+igrid
                     !!c scandisco tutti i punti, anche i bordi (in questa fase 
                     !!c del ciclo ho solo spigoli xy) e li scarto a posteriori
                     itest=mod(n,igrid)
                     if(itest.eq.0..or.itest.eq.lat1) then
                        if (iper(2)) then
                           if(itest.eq.0.)then
                              temp2=phimap1(ix-lat1+ishift)+&
                              &phimap1(ix+lat2)
                           else
                              temp2=phimap1(ix-lat1)+&
                              &phimap1(ix+lat2-ishift)
                           end if
                        else 
                           cycle D2200
                        end if
                     else
                        temp2 = phimap1(ix-lat1)  + phimap1(ix+lat2)
                     end if
                     !!c cosi' non tocco i bordi iy=1 o iy=igrid
                     !!c      segue secondo template x- x+
                     temp1 = phimap1(ix+lat1)  + phimap1(ix+1)
                     temp3 = phimap1(ix-long1) + phimap1(ix+long2)
                     !!c             bisogna sistemare ancora questi due, forse toglierli
                     rtemp2 = phimap1(ix+lat2+lat1)  + phimap1(ix)
                     rtemp3 = phimap1(ix+long2+lat1) + &
                     &phimap1(ix-long1+lat1)
                     
                     phimap2(ix)=.5*(phimap2(ix+lat1)+phimap2(ix))*om1+&
                     &(temp1+ .5*(temp2+temp3+rtemp2+rtemp3))*&
                     &(sf2(ix)+sf2(ix+lat1))*.5
                     !!c now updating last slice
                     phimap2(ix+lat1)=phimap2(ix)
                  end do D2200
                  !----------------------------------------------------------------------
               else
                  D2300: do n=0,long1-igrid
                     ix=ix+igrid; itest=mod(n,igrid)
                     if(itest.eq.0..or.itest.eq.lat1) then
                        if (iper(2)) then
                           if(itest.eq.0.)then
                              temp2=phimap1(ix-lat1+ishift)+&
                              &phimap1(ix+lat2)
                           else
                              temp2=phimap1(ix-lat1)+&
                              &phimap1(ix+lat2-ishift)
                           end if
                        else 
                           cycle D2300
                        end if
                     else
                        temp2 = phimap1(ix-lat1)  + phimap1(ix+lat2)
                     end if
                     temp1 = phimap1(ix+lat1)  + phimap1(ix+1)
                     temp3 = phimap1(ix-long1) + phimap1(ix+long2)
                     
                     phimap2(ix)=phimap2(ix)*om1+(temp1+temp2+temp3)*sixth
                     !!c now updating last slice
                     phimap2(ix+lat1)=phimap2(ix)
                  end do D2300
               end if
               !----------------------------------------------------------------------
               star=(igrid+2)/2; fin=(igrid*(igrid-1)-1)/2
               iy=(igrid*(igrid+2)/2) - igrid +1
               do ix=star,fin
                  iy=iy+igrid
                  bndx3(ix)=phimap2(iy)
                  bndx4(ix)=phimap2(iy+ihgd2)
               end do
            end if
            !!c e+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            do k=icount2a+1,icount2b
               ix=idpos(k)
               temp1=phimap1(ix)*db(1,k)+phimap1(ix+1)*db(2,k)
               temp2=phimap1(ix-lat1)*db(3,k)+phimap1(ix+lat2)*db(4,k)
               temp3=phimap1(ix-long1)*db(5,k)+phimap1(ix+long2)*db(6,k)
               phimap2(ix)=phimap2(ix) + temp1+temp2+temp3
            end do
            !----------------------------------------------------------------------
            !!c reset x  boundary condition
            star=(igrid+2)/2; fin=(igrid*(igrid-1)-1)/2
            iy=(igrid*(igrid+2)/2) - igrid +1
            do ix=star,fin
               iy=iy+igrid
               phimap2(iy)=bndx3(ix)
               phimap2(iy+ihgd2)=bndx4(ix)
            end do
            
            !!C$DIR NO_RECURRENCE 
            do k=icount1a+1,icount1b
               temp=qval(k)
               ix=iqpos(k)
               phimap2(ix)=phimap2(ix) + temp
            end do
            !----------------------------------------------------------------------
            !!c z periodicity
            if(iper(3)) then
               do iz = 2,(igrid-2)**2,2
                  temp1=ibndz(iz);    itemp1=temp1
                  temp2=temp1+idif2z; itemp2=temp2
                  temp3=temp2+inc2za; itemp3=temp3
                  temp4=temp1+inc2zb; itemp4=temp4
                  !!c           iz =1     ===  iz = 64
                  phimap2(itemp1)=phimap1(itemp2)
                  !!c           iz = 65  ===   iz = 2
                  phimap2(itemp3)=phimap1(itemp4)
               end do
            end if
            if(iper(6).and..false.) then
               sumdown=0.0;  sumup=0.0
               do iz = 2,(igrid-2)**2,2
                  temp1=ibndz(iz)
                  temp2=temp1+idif2z; itemp2=temp2
                  temp4=temp1+inc2zb; itemp4=temp4
                  sumup=sumup+phimap1(itemp2)
                  sumdown=sumdown+phimap1(itemp4)
               end do
               tmp=((igrid-2)**2-1)/2.
               sumup=sumup/tmp-vdrop%z
               sumdown=sumdown/tmp+vdrop%z
               do iz = 2,(igrid-2)*(igrid-2),2
                  temp1=ibndz(iz);    itemp1=temp1
                  temp2=temp1+idif2z
                  temp3=temp2+inc2za; itemp3=temp3
                  phimap2(itemp1)=sumup
                  phimap2(itemp3)=sumdown
               end do
            end if
            !----------------------------------------------------------------------
            !!c y periodicity
            if(iper(2)) then
               do iy = 2,(igrid-2)**2,2
                  temp1=ibndy(iy);    itemp1=temp1
                  temp2=temp1+idif2y; itemp2=temp2
                  temp3=temp2+inc2ya; itemp3=temp3
                  temp4=temp1+inc2yb; itemp4=temp4
                  phimap2(itemp1)=phimap1(itemp2)
                  phimap2(itemp3)=phimap1(itemp4)
               end do
            end if
            !----------------------------------------------------------------------
            !!c x periodicity (old way)
            if(iper(1)) then
               do ix = 2,(igrid-2)**2,2
                  temp1=ibndx(ix);    itemp1=temp1
                  temp2=temp1+idif2x; itemp2=temp2
                  temp3=temp2+inc2xa; itemp3=temp3
                  temp4=temp1+inc2xb; itemp4=temp4
                  phimap2(itemp1)=phimap1(itemp2)
                  phimap2(itemp3)=phimap1(itemp4)
               end do
            end if
            !----------------------------------------------------------------------
            ! 2011-06-08 Subroutine omalt is small, thus moved its body here
            if(icheb) then
               !----------------------- Begin body of omalt subroutine -------------
               it=2*i
               om3=1./(1.-om2*spec*0.25)
               if(om1.lt.1.e-6) om3=1./(1.-om2*spec*0.5)
               om4=om3/om2; om2=om3; om1=1.0-om2
               if(rionst.gt.0.0) then
                  if(mod(it,2).eq.1) then
                     do ix=1,nhgp
                        sf1(ix)=sf1(ix)*om4
                     end do
                  else
                     do ix=1,nhgp
                        sf2(ix)=sf2(ix)*om4
                     end do
                  end if
               end if
               
               do ix=1,icount1b
                  qval(ix)=qval(ix)*om4
               end do
               
               do iy=1,6
                  do ix=1,icount2b
                     db(iy,ix)=db(iy,ix)*om4
                  end do
               end do
               sixth=sixth*om4
               !-------------------------- End  body of omalt subroutine -------------
            end if
            !----------------------------------------------------------------------
            !!c we also save time by only checking convergence every ten
            !!c iterations, rather than every single iteration.
            if(mod(i,icon1).eq.(icon1-1)) then
               do ix=2,(icgrid+1)/2,icon2
                  phimap3(ix)=phimap2(ix)
               end do
            end if
            
            if(gten.gt.0.0) then
               grden=0.0
               do ix=1,icount1a
                  iy=iqpos(ix)
                  grden=grden + phimap1(iy)*gval(ix)
               end do
               do ix=icount1a+1,icount1b
                  iy=iqpos(ix)
                  grden=grden+phimap2(iy)*gval(ix)
               end do
               !!c b++++++++++++modified++to save on grdn dimension++++++++++++
               ii=mod(i,5)
               grdn(ii)=grden/2.0
               if(i.gt.10) then
                  igt=.true.
                  do ix=1,5
                     do iy=1,5
                        if(abs(grdn(iy)-grdn(ix)).gt.gten) igt=.false.
                     end do
                  end do
                  if(igt)then
                     write(6,*) (grdn(iy),iy=1,5)
                     ires=1
                  end if
                  !!c e+++++++++++++++++++++++++++++++++++++++++++++++++++++
               end if
            end if
            
            if((mod(i,icon1).eq.0).or.(ires.eq.1)) then
               rnorm2=0
               do ix=2,(icgrid+1)/2,icon2
                  temp2=phimap3(ix)-phimap2(ix)
                  rnorm2=rnorm2+temp2**2
                  rmxch=max(rmxch,abs(temp2))
               end do
               rmsch = sqrt(real(icon2)*rnorm2/npoint)
               rnormch= sqrt(rnorm2)
               rmsch2=rmsch
               rmxch2=rmxch 
               if(verbose) then

                  if(ideveloper) then
                     if(gten.gt.0.)then
                        write(6,*) rmsch2,rmxch2,grden,' at',i,' iterations'
                     else
                        write(6,*) rmsch2,rmxch2,' at',i,'iterations'
                     endif
                  else
                     if(gten.gt.0.)then
                        write(6,"(3e13.5,a,i8,a)") rmsch2,rmxch2,grden,' at',i,' iterations'
                     else
                        write(6,"(2e13.5,a,i8,a)") rmsch2,rmxch2,' at',i,' iterations'
                     endif
                  end if



               end if
               !!c b++++++++changed and to or++++++++++++++++++++
               if(rmsch.le.res1.or.rmxch.le.res2.or.rnormch.le.res5) ires=1
               !!c e+++++++++++++++++++++++++++++++++++++++++++++
               if(igraph.and.(once)) then 
                  do j=i-9,i
                     ibin = (j-1.)*(nxran-1.)/(nlit-1.) + 1
                     rmsl(ibin) = rmsch
                     rmaxl(ibin) = rmxch
                  end do
               end if
            end if
            !!c check to see if accuracy is sufficient, or if a qstop command
            !!c has been issued
            inquire(file='qstop.test',exist=qstopper)
            if(qstopper) ires=1
            !!c LOOP
            i=i+1
            if (max(gten,res1,res2,res5).lt.tol) then
               !!c     if(gten.lt.1.e-6.and.res1.lt.1.e-6.and.res2.lt.1.e-6) then
               if(i.le.nlit.and.ires.eq.0) then
                  cycle D1000
               else
                  exit
               end if
            else
               if((i.le.nlit.or.iautocon).and.ires.eq.0) then
                  cycle D1000
               else
                  exit
               end if
            end if
         end do D1000
         
         !!c	end of iteration loop
         !----------------------------------------------------------------------
         !!c       remap into phimap
         do iy=1,(icgrid-1)/2 
            ix=iy*2 
            phimap3(ix-1)=phimap1(iy)
            phimap3(ix)=phimap2(iy)
         end do
         if(once) then
            iw=1 
            do iz=1,igrid
               do iy=1,igrid
                  do ix=1,igrid
                     phimap(ix,iy,iz)=phimap3(iw)
                     iw=iw+1
                  end do
               end do
            end do
            phimap(igrid,igrid,igrid)=phimap1((icgrid+1)/2)
         else
            iw=1 
            do iz=1,igrid
               do iy=1,igrid
                  do ix=1,igrid
                     phimap(ix,iy,iz)=phimap(ix,iy,iz)-phimap3(iw)
                     iw=iw+1
                  end do
               end do
            end do
            phimap(igrid,igrid,igrid)=phimap(igrid,igrid,igrid)  &
            &  - phimap1((icgrid+1)/2)
         end if
         !!c ++da cacciare!!++++++++++++++++++++++++++++++++
         if (.false.) then
            write(*,*)"pot000:",phimap((igrid+1)/2,(igrid+1)/2,(igrid+1)/2)
         end if
         !!c ++da cacciare!!++++++++++++++++++++++++++++++++
         if (.false.) then
            write(6,*)"lo faccio"
            open(52,file='potinplane',form='formatted')
            do ix=(igrid+1)/2-10,(igrid+1)/2+10
               do iy=(igrid+1)/2-10,(igrid+1)/2+10
                  write(52,*)ix-(igrid+1)/2,iy-(igrid+1)/2,&
                  & phimap(ix,iy,(igrid+1)/2)
               end do
            end do
            close(52)
            open(52,file='Green10.dat',form='formatted')
            do ix=-10,10
               do iy=-10,10
                  do iz=-10,10
                     read(52,*)n,m,l,Green(n,m,l)
                  end do
               end do
            end do
            close(52)
            open(52,file='Green10.bin',form='unformatted')
            write(52)Green
            close(52)
            
         end if
         !!c e++++++++++++++++++++++++++++++++++++++++++++++
         !!c e++++++++++++++++++++++++++++++++++++++++++++++
         
         if(diff.and.(once)) then
            write(6,*) 'now doing uniform dielectric run....'
            once=.false.
            om3= 2.0/(1 +(3.14159265359/real(igrid)))
            sixth=sixth*om3/om2
            do i=1,icgrid+1
               phimap3(i)=phimap3(i)*epsrat 
            end do
            do i=1,icount1b
               qval(i)=qval(i)*om3/om2
            end do
            do i=1,ibc
               itemp=cgbp(i)%value2
               qval(itemp)=cgbp(i)%value1*om3
            end do
            om1=1-om3 ; icount2a=0;  icount2b=0
            nlit=int(7.8*real(igrid)/pi)
            i=1
            ! 2011-06-07 GOTO is removed due to DO-EXIT construct 
         else
            exit D1689
         end if
      end do D1689
      icount2a=icount2abac;  icount2b=icount2bbac
      
      call date_and_time(DATE=day,TIME=time,VALUES=values)
      timediff = dtime( timearray)
      timetot=timetot+timediff
      finish=values(5)*3600.+values(6)*60.+values(7)+values(8)/1000.
      if(verbose) then
         write(6,*)'finished qdiffx linear iterations at: ',&
         &time(1:2),':',time(3:4),':',time(5:6)

         if(ideveloper) then
            write(6,*)'total time elapsed so far: ',timediff,' s'
         else
            write(6,"(a,f10.2,a)")' total time elapsed so far: ',timediff,' s'
         end if

         !write(6,*)'# loops                  : ',(i-1)

         if(ideveloper) then
            write(6,*)'mean,max change (kT/e)   : ',rmsch2,rmxch2
         else
         end if

      end if
      !----------------------------------------------------------------------
      !!c  plot convergence history   
      if(igraph) then
         iclr = 1; iscl = 1;  imk = 0; iplt = 0;  symb = 'M'
         title = '    linear iteration convergence history   '
         call conplt(rmaxl,title,iclr,iscl,imk,iplt,symb,1,nlit,&
         & iplot,ymin,ymax)
         iclr = 0
         call conplt(rmsl,title,iclr,iscl,imk,iplt,symb,1,nlit,&
         & iplot,ymin,ymax)
         iscl = 0; imk = 1
         call conplt(rmaxl,title,iclr,iscl,imk,iplt,symb,1,nlit,&
         & iplot,ymin,ymax)
         symb = 'A' ;  iplt = 1
         call conplt(rmsl,title,iclr,iscl,imk,iplt,symb,1,nlit,&
         & iplot,ymin,ymax)
      end if
      !----------------------------------------------------------------------
      !!c code phimap corner, for use in transference from irises to convex
      !!c and via versa
      ap1=phimap(1,1,1)
      ap2=ap1*10000
      ap3=real(int(ap2))
      if(ap3.gt.0) then
         ap4=(ap3+0.8)/10000
      else
         ap4=(ap3-0.8)/10000
      end if
      phimap(1,1,1)=ap4
      
      if(ipoten) then
         midg = (igrid+1)/2
         do m = 1,5
            n = (igrid - 1)/4; nn = (m-1)*n + 1
            write(6,*)'phi',nn,midg
            write(6,*)(phimap(nn,midg,ii),ii=1,igrid)
         end do
      end if
      end subroutine itit
