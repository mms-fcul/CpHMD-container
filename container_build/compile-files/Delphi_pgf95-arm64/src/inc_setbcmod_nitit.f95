      !  ATTENTION!  This file is part of setbcmod module.
      !==================================================================
      ! 2011-06-09 All other parameters are accessible via qlog and pointers
      subroutine nitit(qfact)
      !	subroutine nitit(idpos,db,sf1,sf2,iqpos,qval,icount2a,icount2b,
      
      !!c NOTE THIS HAS BEEN ALTERED TO GIVE THE CONVERGENCE TO THE EXACT
      !!C SOLUTION, WHICH IS INPUTTED VIA OPTION 3,(AND NOW WITH GAUSS
      !!C SIEDEL IMPLEMENTED, WITH MAPPED SUBARRAYS.)
      
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
      
      
      integer, parameter :: nxran=60, nyran=60
      real :: rmsl(nxran),rmsn(nxran),rmaxl(nxran),rmaxn(nxran)
      
      ! 2011-06-09 Arrays and variables below are accessible via qlog
      character(24) :: day
      character :: symb
      character(70) :: title,iplot(nyran)
      logical :: qstopper,resdat
      integer :: star,fin,sta1(igrid),sta2(igrid),fi1(igrid),fi2(igrid)
      !!c b+++++++++++++++++++++++++
      integer :: icont,icountplus,itshift
      logical :: ichangeom,inew,inewfirst,istop
      real :: conv(3),omcomp,tmp,der,relparprev,qfact,deb
      real :: factor,delom,linrelpar,cost
      character(18) :: nlstr
      !--------------------------------------------------------------------------
      ! 2011-06-09 Declarations added due to IMPLICIT NONE
      real :: start,finish,fraction,sixth,om1,om2,rmsch,rmsch2
      real :: temp,temp1,temp2,temp3,temp4,rmxch,rmxch2,rnorm2
      real :: rnormch,ymin,ymax,debfct,derprec,fac1
      integer :: i,j,k,n,m,ii,iadd1,iadd2,icgrid,ibin,iclr,imk,iplt
      integer :: idif1x,idif2x,idif1y,idif2y,idif1z,idif2z
      integer :: inc1xa,inc1xb,inc2xa,inc2xb,ihgd,ihgd2
      integer :: inc1ya,inc1yb,inc2ya,inc2yb,iscl,iw,ix,iy,iz
      integer :: inc1za,inc1zb,inc2za,inc2zb,itnum,lat1,lat2
      integer :: itemp1,itemp2,itemp3,itemp4,long1,long2
      integer :: mcon1,midg,nn,npoint,ires
      real*4  timediff, timearray(2), dtime
      !--------------------------------------------------------------------------
      nlstr='                  '; icountplus=0; factor=1.0;itshift=0
      ichangeom=.true. ; inewfirst=.false.; istop=.true.; inew=.true.
      !!c e+++++++++++++++++++++++++
      debfct = epsout/(deblen*scale)**2
      cost=debfct/(2.*rionst*epsout)
      npoint = (igrid-2)**3
      !!c-------------------------------------------------------
      call date_and_time(DATE=day,TIME=time,VALUES=values)
      timediff = dtime( timearray)
      timetot=timetot+timediff
      start=values(5)*3600.+values(6)*60.+values(7)+values(8)/1000.
      if(verbose)write(6,*)'now iterating at: ',&
      &time(1:2),':',time(3:4),':',time(5:6)
      !!c-------------------------------------------------------
      !!c some initialization
      !!c b+++++++++++++++++++++++++++++++++++++++++++++++++++++
      conv=0.
      !!c e++++++++++++++++++++++++++++++++++++++++++++++++++++
      tmp=abs(chi2*chi4) ; icont=0; sixth = 1./6.
      icgrid=igrid*igrid*igrid ; ihgd=(igrid+1)/2
      if(icon2.eq.0) then
         mcon1=10; icon2=1
      end if
      itnum=0 ; fraction=0.0 
      if(icon1.gt.nlit) icon1=nlit 
      j=0
      do  iz=1,igrid
         do  iy=1,igrid
            do ix=1,igrid
               j=j+1; deb=0.
               if (idebmap(ix,iy,iz)) deb=1.
               phimap3(j)=deb
            end do
         end do
      end do
      do ix=1,nhgp
         iy=ix*2
         debmap1(ix)=phimap3(iy-1)
         debmap2(ix)=phimap3(iy)
      end do
      rmsl=0.; rmsn=0.; rmaxl=0.; rmaxn=0.
      
      !!c ---------------------------------------------
      !!c MAIN SET UP ROUTINE	
      !!c ---------------------------------------------
      if(iper(1)) then
         n=0
         do iz=2,igrid-1
            iadd1=(iz-1)*igrid*igrid 
            do iy=2,igrid-1
               iadd2=(iadd1+(iy-1)*igrid +2)/2
               n=n+1; ibndx(n)=iadd2
            end do
         end do
         idif1x=(igrid-2)/2;   idif2x=idif1x+1
         inc1xa=1; inc1xb=0; inc2xa=0; inc2xb=1
      end if
      if(iper(2)) then
         n=0
         do iz=2,igrid-1
            iadd1=(iz-1)*igrid*igrid 
            do ix=2,igrid-1
               iadd2=(iadd1+ix+1)/2
               n=n+1; ibndy(n)=iadd2
            end do
         end do
         idif1y=igrid*(igrid-2)/2 ; idif2y=idif1y+1
         inc1ya=(igrid/2)+1 ;       inc1yb=inc1ya-1
         inc2ya=inc1yb;             inc2yb=inc1ya
      end if
      if(iper(3)) then
         n=0
         do ix=2,igrid-1
            iadd1=ix+1
            do iy=2,igrid-1
               iadd2=(iadd1+(iy-1)*igrid)/2
               n=n+1;  ibndz(n)=iadd2
            end do
         end do
         idif1z=igrid*igrid*(igrid-2)/2; idif2z=idif1z+1
         inc1za=((igrid**2)/2)+1 ;       inc1zb=inc1za
         inc2za=inc1zb;                  inc2zb=inc1za
      end if
      
      !!c remove qstopper file if it already exists
      inquire(file='qstop.test',exist=qstopper)
      if(qstopper) then
         open(30,file='qstop.test')
         close(30,status='delete')
         qstopper=.false.
      end if
      !!c check for resolution data
      
      !!c set up start and stop vectors
      sta1(2)=(igrid*igrid + igrid +4)/2;  sta2(2)=sta1(2)-1
      fi1(2)=igrid**2 - (igrid+1)/2 ;       fi2(2)=fi1(2)
      itemp1=igrid + 2 ;   itemp2=igrid*igrid -igrid -2
      do i=3,igrid-1
         sta1(i)=fi1(i-1) + itemp1 ;  sta2(i)=fi2(i-1) + itemp1
         fi1(i)=sta1(i-1) + itemp2 ;   fi2(i)=sta2(i-1) + itemp2
      end do
      !!c also
      lat1= (igrid-1)/2;        lat2= (igrid+1)/2
      long1= (igrid**2 - 1)/2; long2= (igrid**2 + 1)/2
      ires=0; iw=1
      
      do iz=1,igrid
         do iy=1,igrid
            do ix=1,igrid
               iw=iw+1
               phimap3(iw)=phimap(ix,iy,iz)
            end do
         end do
      end do
      do ix=2,icgrid+1,2
         iy=ix/2 ; phimap1(iy)=phimap3(ix); phimap2(iy)=phimap3(ix+1)
      end do
      
      star=(igrid+1)/2;  fin=(igrid*(igrid-1)-2)/2
      iy=(igrid*(igrid+1)/2) - igrid + 1 ; ihgd2=ihgd-1
      do ix=star,fin
         iy=iy+igrid; bndx1(ix)=phimap1(iy)
         bndx2(ix)=phimap1(iy+ihgd2)
      end do
      
      star=(igrid+2)/2 ; fin=(igrid*(igrid-1)-1)/2
      iy=(igrid*(igrid+2)/2) - igrid +1
      do ix=star,fin
         iy=iy+igrid ; bndx3(ix)=phimap2(iy); bndx4(ix)=phimap2(iy+ihgd2)
      end do
      
      om2=2.0/(1.0 + sqrt(1 - spec))
      if (.not.imanual) ichangeom=.false.
      relparprev=om2;  linrelpar=om2

      if(ideveloper) then
         write(6,*)'linear rel. parameter = ',linrelpar
      else
         write(6,"(a,f8.4)")' linear rel. parameter = ',linrelpar
      end if

      if (imanual) then
         write(6,"(a,f8.4)")' non linear fixed rel. parameter =',relpar 
      else
         write(6,"(a,f8.4)")' non linear initial rel. parameter =',relpar
         write(6,"(a,f8.4)")' q factor',qfact
      end if
      om1=1.0-om2
      write(6,*) ' ' ;  write(6,*) ' '
      write(6,*) '  rms-change     max change         #iterations'
      
      do ix=1,(icgrid+1)/2
         sf1(ix)=sf1(ix)*om2 ;  sf2(ix)=sf2(ix)*om2
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
      i=1
      !------------------------------------------------------------------
      ! 2011-06-09 DO-EXOT construct in iteration loop instead of IF-GOTO
      D1000: do
         !!c clear rms, max change
         rmsch = 0.0 ;  rmxch = 0.0
         !!c if there is no salt then the main loop is executed without sf
         !!c saving about 15% in execution time
         if(rionst.gt.0.0) then
            do n = 2, igrid-1
               star=sta1(n);  fin=fi1(n)
               do ix = star,fin
                  temp1 = phimap2(ix)+phimap2(ix-1)
                  temp2 = phimap2(ix+lat1)+phimap2(ix-lat2)
                  temp3 = phimap2(ix+long1)+phimap2(ix-long2)
                  phimap1(ix) =phimap1(ix)*om1 +(qmap1(ix) + &
                  &temp1+temp2+temp3)*sf1(ix)
                  !!c e++++++++++++++++++++++++++++temporaneo?**************
               end do
            end do
            !!c otherwise the main loop is as below:
         else
            do n = 2, igrid-1
               star=sta1(n)
               fin=fi1(n)
               do ix = star,fin
                  temp1 = phimap2(ix)+phimap2(ix-1)
                  temp2 = phimap2(ix+lat1) + phimap2(ix-lat2)
                  temp3 = phimap2(ix+long1) + phimap2(ix-long2)
                  phimap1(ix) = phimap1(ix)*om1 + (temp1+temp2+temp3)*sixth
               end do
            end do
         end if
         !!c the above loops are about fourtimes faster than the original
         !!c loop over all grid points for several reasons, the biggest being that
         !!c we are only solving laplace's equation (unless salt is present), which
         !!c numerically much simpler, hence faster. we put all we leave out, back
         !!c in below, ending up with an equivalent calculation, but much faster.
         
         !!c first we add back the dielectric boundary points, by recalculating them
         !!c individually. note this is still vectorised by means of a gathering
         !!c load by the compiler.
         
         do k=1,icount2a
            ix=idpos(k)
            temp1=phimap2(ix-1)*db(1,k)+phimap2(ix)*db(2,k)
            temp2=phimap2(ix-lat2)*db(3,k)+phimap2(ix+lat1)*db(4,k)
            temp3=phimap2(ix-long2)*db(5,k)+phimap2(ix+long1)*db(6,k)
            phimap1(ix)= phimap1(ix) + temp1+temp2+temp3
         end do
         !!c       end if
         !!c
         !!c next we add back an adjustment to all the charged grid points due to
         !!c the charge assigned. the compiler directive just reassures the vector
         !!c compiler that all is well as far as recurrence is concerned, i.e. it
         !!c would think there is a recurrence below, where as in fact there is none.
         
         !!c Now reset boundary values altered in above loops.
         star=(igrid+1)/2 ; fin=(igrid*(igrid-1)-2)/2
         iy=(igrid*(igrid+1)/2) - igrid +1
         do ix=star,fin
            iy=iy+igrid; phimap1(iy)=bndx1(ix); phimap1(iy+ihgd2)=bndx2(ix)
         end do
         do k=1,icount1a
            temp=qval(k); ix=iqpos(k); phimap1(ix)=phimap1(ix) + temp
         end do
         !!c if periodic boundary condition option
         !!c force periodicity using wrap around update of boundary values:
         !!c 2nd slice-->last
         !!c last-1 slice-->first
         
         !!c z periodicity
         if(iper(3)) then
            do iz = 1,(igrid-2)**2,2
               temp1=ibndz(iz);    itemp1=temp1
               temp2=temp1+idif1z; itemp2=temp2
               temp3=temp2+inc1za; itemp3=temp3
               temp4=temp1+inc1zb; itemp4=temp4
               phimap1(itemp1)=phimap2(itemp2)
               phimap1(itemp3)=phimap2(itemp4)
            end do
         end if
         
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
         
         !!c x periodicity
         if(iper(1)) then
            do ix = 1,(igrid-2)**2,2
               temp1=ibndx(ix);    itemp1=temp1
               temp2=temp1+idif1x; itemp2=temp2
               temp3=temp2+inc1xa; itemp3=temp3
               temp4=temp1+inc1xb; itemp4=temp4
               phimap1(itemp1)=phimap2(itemp2)
               phimap1(itemp3)=phimap2(itemp4)
            end do
         end if
         
         !!c Next update phimap2 using the new phimap1
         
         if(rionst.gt.0.0) then	
            do n = 2, igrid-1
               star=sta2(n);  fin=fi2(n)
               do ix = star,fin
                  temp1 = phimap1(ix)+ phimap1(ix+1)
                  temp2 = phimap1(ix+lat2)+phimap1(ix-lat1)
                  temp3 = phimap1(ix+long2) + phimap1(ix-long1)
                  phimap2(ix) =phimap2(ix)*om1 + &
                  &(qmap2(ix)+temp1+temp2+temp3)*sf2(ix)
               end do
            end do
         else
            do n = 2, igrid-1
               star=sta2(n); fin=fi2(n)
               do ix = star,fin
                  temp1 = phimap1(ix)+phimap1(ix+1)
                  temp2 = phimap1(ix+lat2) + phimap1(ix-lat1)
                  temp3 = phimap1(ix+long2) + phimap1(ix-long1)
                  phimap2(ix) =phimap2(ix)*om1 + (temp1+temp2+temp3)*sixth
               end do
            end do
         end if
         do k=icount2a+1,icount2b
            ix=idpos(k)
            temp1=phimap1(ix)*db(1,k)+phimap1(ix+1)*db(2,k)
            temp2=phimap1(ix-lat1)*db(3,k)+phimap1(ix+lat2)*db(4,k)
            temp3=phimap1(ix-long1)*db(5,k)+phimap1(ix+long2)*db(6,k)
            phimap2(ix)=phimap2(ix) + temp1+temp2+temp3
         end do
         !!c       end if
         !!c reset boundary condition
         star=(igrid+2)/2; fin=(igrid*(igrid-1)-1)/2
         iy=(igrid*(igrid+2)/2) - igrid +1
         do ix=star,fin
            iy=iy+igrid; phimap2(iy)=bndx3(ix); phimap2(iy+ihgd2)=bndx4(ix)
         end do
         do k=icount1a+1,icount1b
            temp=qval(k); ix=iqpos(k); phimap2(ix)=phimap2(ix) + temp
         end do
         
         !!c z periodicity
         if(iper(3)) then
            do iz = 2,(igrid-2)**2,2
               temp1=ibndz(iz);    itemp1=temp1
               temp2=temp1+idif2z; itemp2=temp2
               temp3=temp2+inc2za; itemp3=temp3
               temp4=temp1+inc2zb; itemp4=temp4
               phimap2(itemp1)=phimap1(itemp2)
               phimap2(itemp3)=phimap1(itemp4)
            end do
         end if
         
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
         
         !!c x periodicity
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
         
         !!c we also save time by only checking convergence every ten
         !!c iterations, rather than every single iteration.
         
         !!c store phi2 in phi3 to compare against next iteration
         if(mod(i,icon1).eq.(icon1-1)) then
            do ix=2,(icgrid+1)/2,icon2
               phimap3(ix)=phimap2(ix)
            end do
         end if
         
         !!c check convergence
         if(mod(i,icon1).eq.0) then
            !!c rmsch= rms change
            !!c store in rmsch2
            !!c rmxch= max change
            !!c store in rmxch2
            rmsch2=rmsch; rmxch2=rmxch ; rnorm2=0
            ! 2011-06-09 Did not understand meaning of those two IFs
            do ix=2,(icgrid+1)/2,icon2
               temp2=phimap3(ix)-phimap2(ix)
               rnorm2=rnorm2+temp2*temp2
               rmxch=max(rmxch,abs(temp2))
            end do
            !!c b+++++++++++++++++++++++++++++++++++++++++++++++++++++++
            conv(3)=conv(2); conv(2)=conv(1); conv(1)=rmxch
            !!c e+++++++++++++++++++++++++++++++++++++++++++++++++++++++
            rmsch = sqrt(real(icon2)*rnorm2/npoint)
            rnormch=sqrt(rnorm2)/relpar
            !!c b+++++++++++++++++++++++++++++++++++
            if((rmsch.le.res1.or.rmxch.le.res2.or.&
            & rnormch.le.res5).and.itnum.gt.22) ires=1
            !!c e+++++++++++++++++++++++++++++++++++
            if(itnum.eq.0) then 

               if(ideveloper) then
                  write(6,*) rmsch,rmxch,' at ',i,'iterations'
               else
                  write(6,"(2e13.5,a,i8,a)") rmsch,rmxch,' at',i,' iterations'
               end if
               !+++++++++++++++++++++++++++++++++++'

               istop=.not.(rmxch.ge.0.22)
               if(igraph) then 
                  do j=i-9,i
                     ibin = (j-1.)*(nxran-1.)/(nlit-1.) + 1
                     rmsl(ibin) = rmsch;  rmaxl(ibin) = rmxch
                  end do
               end if
               !!c b+++++++++++++++++++++++++++++++++++++
               !!c ottimizzazione di primo passo
               inewfirst=(i.ge.nlit-mod(nlit,icon1).and.&
               &istop.and. itnum.eq.0)
               if (.not.imanual.and.inewfirst.and.qfact.gt.3.2) then
                  factor=exp(-qfact*2.1)+1.E-6; ichangeom=.true.
               end if
            end if
            !!c nonlinear part
            inew=inew.or.inewfirst
            if (itnum.ne.0.) then

               if(ideveloper) then

                  write(6,*) rmsch,rmxch,itnum,'it. ',nlstr
               else

                  write(6,"(2e13.5,i8,2a)") rmsch,rmxch,itnum,' it. ',nlstr
               end if
 

               if (.not.imanual) then
                  derprec=der;  der=(conv(1)-conv(2))/conv(1)
                  if(rmxch.lt.1.e-6) then
                     factor=1.2; ichangeom=.true.
                  else  
                     if (der.gt.0..and.(.not.inew)) then
                        icountplus=icountplus+1
                        factor=factor*(1.-der)**.99
                        !!c                 write(6,*)'factor:',factor
                        if (der.gt.0.55) then
                           ichangeom=.true.
                           if (der.ge.1.) factor=1.e-5
                        end if
                        if(der.gt.0.35.and.conv(1).gt..1)then
                           ichangeom=.true.
                           factor=(factor*.05/conv(1))**4
                        end if
                     end if
                     if ((der.gt.0.and.conv(1).gt..1).and.inew.and..not.&
                     & inewfirst) then
                        !!c               if ((der.gt.0.and.conv(1).gt..1).and..not.inew.and..not.
                        ichangeom=.true.
                        factor=min((factor*.05/conv(1)),&
                        & factor*(1.-der)**.86)
                     end if
                     if (der.le.0.) then
                        icountplus=0
                        !!c                  write(6,*)'fatto',relpar,itnum,der,derprec
                        factor=1.
                        if (itnum.gt.24.and.itnum.lt.24+.75*(nnit-24)&
                        &.and.rmxch.lt..03.and.derprec.le.0.) then
                           if(der.gt.-.2.and.derprec.gt.-.2) then
                              write(6,*) 'Trying to speed up&
                              & the convergence process'
                              factor=1.1; ichangeom=.true.
                              if(relpar.lt..2.and.der.gt.-.05.and. derprec &
                              & .gt. -.05) factor=1-45.226*(relpar-.2)
                           end if
                        end if
                     end if
                     if (icountplus.ge.2) ichangeom=.true.
                     !!c               write(6,*) 'der:',der
                  end if
                  inewfirst=.false.
               end if
               !!c e++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            end if
            !!c end of convergence check
         end if
         !!c check to see if accuracy is sufficient, or if a qstop command
         !!c has been issued
         !!ccommented out qstopper as obsoleteinquire(file='qstop.test',exist=qstopper)
         !!c          if(qstopper) ires=1
         
         !!c LOOP
         i=i+1
         !!c b++++++++++++++++++++++++++++++++++++++++++++
         if((i.le.nlit.or..not.istop).and.ires.eq.0)  cycle D1000
         !!c e++++++++++++++++++++++++++++++++++++++++++++
         !!c	end of iteration loop
         !!c first pass write header
         if(nnit.gt.0.and.itnum.eq.0) then
            write(6,*) ' '
            write(6,*) 'now for the non-linear iterations'
            write(6,*) ' '
            write(6,*) '  rms-change     max change         #iterations'
         end if
         
         !!c icon1 = ogni quanti blocchi di iterazioni verifica convergenza
         !!c nlit = quante iterazioni nel blocco 
         !!c b+++OCT 2000
         icon1=10; nlit=10
         !!c e+++++++++++++
         itnum=itnum+1
         if(itnum.gt.nnit.or.ires.eq.1.and..not.inew) exit D1000
         i=1
         !!c b++++++++OCT 2000++++++++++++++++++++++++
         inew=ichangeom
         if (ichangeom) then
            relpar=relpar*factor
            if (relpar.lt.1.E-4) then
               write(6,*)'estimation ',relpar,' 1E-4 preferred'
               relpar=1.E-4
            end if
            factor=1.

            if(ideveloper) then
               write(6,*)'                 New relaxation parameter =',&
               & relpar
            else
               write(6,"(a,f10.2)")'                  &               
                  & New relaxation parameter =',relpar
            end if

            ichangeom=.false.
            icountplus=0;  omcomp=relpar/relparprev
            relparprev=relpar;  om1=1.0-relpar
            do ix=1,(icgrid+1)/2
               sf1(ix)=sf1(ix)*omcomp
               sf2(ix)=sf2(ix)*omcomp
            end do
            do ix=1,icount1b
               qval(ix)=qval(ix)*omcomp
            end do
            do iy=1,6
               do ix=1,icount2b
                  db(iy,ix)=db(iy,ix)*omcomp
               end do
            end do
            sixth=sixth*omcomp; icont=icont+1
         end if
         
         fraction= fraction + 0.05
         if(fraction.gt.1.0) then
            fraction=1.0;  nlstr='full non-linearity'
         end if
		 
         fac1=fraction*cost
         if (tmp.lt.1.e-6.and..false.) then
            do ix=1,nhgp
               temp1=phimap1(ix)*debmap1(ix)
               temp2=phimap2(ix)*debmap2(ix)
               temp3=temp1*temp1
               temp4=temp2*temp2
               qmap1(ix)=fac1*temp3*temp1*(chi3 + temp3*chi5)
               qmap2(ix)=fac1*temp4*temp2*(chi3 + temp4*chi5)
               
               !!c b++++soglia messa per incrementare stabilita' della convergenza
               !2012-04-24 Chuan set cutoff values of qmap1 and qmap2 nonzero
               if (temp3.gt.2500) then
                  qmap1(ix)=fac1*2500.*50.*(chi3+2500.*chi5)
               end if             
               if (temp4.gt.2500) then 
                  qmap2(ix)=fac1*2500.*50.*(chi3+2500.*chi5)
               end if
               !!c e+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               !!c           qmap2(ix)=fac1*(-2.*rionst*(sinh(temp2)-temp2))
            end do
         else
            do ix=1,nhgp
               temp1=phimap1(ix)*debmap1(ix)
               temp2=phimap2(ix)*debmap2(ix)
               temp3=temp1*temp1
               temp4=temp2*temp2
               qmap1(ix)=fac1*temp3*(chi2+temp1*(chi3 + &
               &temp1*(chi4+temp1*chi5)))
               qmap2(ix)=fac1*temp4*(chi2+temp2*(chi3 + &
               &temp2*(chi4+temp2*chi5)))
               
               !!c b++++soglia messa per incrementare stabilita' della convergenza
               !2012-04-24 Chuan set cutoff values of qmap1 and qmap2 nonzero
               if (temp3.gt.2500) then
                  qmap1(ix)=fac1*2500.*(chi2+50.*(chi3+&
                  &50.*(chi4+50.*chi5)))
               end if             
               if (temp4.gt.2500) then 
                  qmap2(ix)=fac1*2500.*(chi2+50.*(chi3+&
                  &50.*(chi4+50.*chi5)))
               end if
               !!c e+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               
               !!c e+++++++++++++++++++++++++
            end do
         end if
         !!c now go do a couple of linear iterations
      end do D1000
      
      do iy=1,(icgrid-1)/2 
         ix=iy*2;  phimap3(ix)=phimap2(iy);  phimap3(ix-1)=phimap1(iy)
      end do
      iw=1 
      do iz=1,igrid
         do iy=1,igrid
            do  ix=1,igrid
               phimap(ix,iy,iz)=phimap3(iw); iw=iw+1
            end do
         end do
      end do
      phimap(igrid,igrid,igrid)=phimap1((icgrid+1)/2)
      !!c b++++++++++++++++++++++++++++++++++++++++++++
      if (relpar.lt.0.05) then
         write(6,*) 'Convergence is more reliable &
         &if relaxation parameter>0.05'
         write(6,*) 'If it is possible, it is&
         & advisable to increase it'
         write(6,*)' '
      end if
      !!c e++++++++++++++++++++++++++++++++++++++++++++++
      call date_and_time(DATE=day,TIME=time,VALUES=values)
      finish=values(5)*3600.+values(6)*60.+values(7)+values(8)/1000.
      write(6,*)' '
      write(6,*)'finished qdiffx linear iterations at: ',&
      &time(1:2),':',time(3:4),':',time(5:6)
      timediff = dtime( timearray)
      timetot=timetot+timediff

      if(ideveloper) then
         write(6,*)'time taken : ',timediff,' s'
      else
         write(6,"(a,f10.2,a)")' time taken : ',timediff,' s'
      end if


!      write(6,*)'# full non-linear loops  : ',(i-1)
!      write(6,*)'mean,max change (kT/e)   : ',rmsch,rmxch
      
      !!c  plot convergence history   
      if(igraph) then
         iclr = 1; iscl = 1; imk = 0; iplt = 0;  symb = 'M'
         title = '    linear iteration convergence history   '
         call conplt(rmaxl,title,iclr,iscl,imk,iplt,symb,1,nlit,&
         &  iplot,ymin,ymax)
         iclr = 0
         call conplt(rmsl,title,iclr,iscl,imk,iplt,symb,1,nlit,&
         &  iplot,ymin,ymax)
         iscl = 0; imk = 1
         call conplt(rmaxl,title,iclr,iscl,imk,iplt,symb,1,nlit,&
         &  iplot,ymin,ymax)
         symb = 'A'; iplt = 1
         call conplt(rmsl,title,iclr,iscl,imk,iplt,symb,1,nlit,&
         &  iplot,ymin,ymax)
      end if
      
      !!c	give some intermediate output of phi
      if(ipoten) then
         midg = (igrid+1)/2
         do m = 1,5
            n = (igrid - 1)/4;  nn = (m-1)*n + 1
            write(6,*)'phi',nn,midg
            write(6,*)(phimap(nn,midg,ii),ii=1,igrid)
         end do
      end if
      end subroutine nitit
