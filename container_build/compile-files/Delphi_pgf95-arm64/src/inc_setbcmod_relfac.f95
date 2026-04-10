      !  ATTENTION!  This file is part of setbcmod module.
      !==================================================================
      ! 2011-06-05 All parameters are accessible via qlog and pointers
      subroutine relfac
      !	subroutine relfac(idpos,db,sf1,sf2,icount2a,icount2b,spec,ibnum,
      
      ! 2011-06-07 Variables and arrays below are accessible via
      ! 2011-06-07 Some local arrays and variables
      real :: sn1(igrid),sn2(igrid),sn3(igrid)
      integer :: star,fin,sta1(igrid),sta2(igrid),fi1(igrid),fi2(igrid)
      ! 2011-06-07 Declarations added due to IMPLICIT NONE
      real :: sixth=1./6.,temp,rtemp2,rtemp3,recipr
      real :: temp1,temp2,temp3,temp4
      integer :: ix,iy,iz,icgrid,ihgd,ihgd2,i,k,iadd1,iadd2,n
      integer :: idif1z,idif2z, idif1y,idif2y,inc1ya,inc1yb
      integer :: inc1za,inc1zb,inc2za,inc2zb,inc2ya,inc2yb
      integer :: itemp1,itemp2,itemp3,itemp4,iw,lat1,lat2
      integer :: long1,long2
      
      phimap1=0.0
      phimap2=0.0
      phimap3=0.0
      icgrid=igrid**3; ihgd=(igrid+1)/2
      !!c b++++++++++++debug++++++++++++++++++++++++
      if(iper(2)) then
         n=0
         do iz=2,igrid-1
            iadd1=(iz-1)*igrid*igrid 
            do  ix=2,igrid-1
               iadd2=(iadd1+ix+1)/2
               n=n+1; ibndy(n)=iadd2
            end do
         end do
         idif1y=igrid*(igrid-2)/2; idif2y=idif1y+1
         inc1ya=(igrid/2)+1; inc1yb=inc1ya-1
         inc2ya=inc1yb; inc2yb=inc1ya
      end if
      if(iper(3)) then
         n=0
         do ix=2,igrid-1
            iadd1=ix+1
            do  iy=2,igrid-1
               iadd2=(iadd1+(iy-1)*igrid)/2
               n=n+1; ibndz(n)=iadd2
            end do
         end do
         idif1z=igrid*igrid*(igrid-2)/2; idif2z=idif1z+1
         inc1za=((igrid**2)/2)+1; inc1zb=inc1za
         inc2za=inc1zb; inc2zb=inc1za
      end if
      !!c
      !!c set up start and stop vectors
      sta1(2)=(igrid**2 + igrid +4)/2; sta2(2)=sta1(2)-1
      fi1(2)=igrid**2 - (igrid+1)/2;   fi2(2)=fi1(2)
      itemp1=igrid + 2;  itemp2=igrid**2 -igrid -2
      do i=3,igrid-1
         sta1(i)=fi1(i-1) + itemp1;  sta2(i)=fi2(i-1) + itemp1
         fi1(i)=sta1(i-1) + itemp2;   fi2(i)=sta2(i-1) + itemp2
      end do
      
      !!c also
      lat1= (igrid-1)/2; lat2= (igrid+1)/2
      long1= (igrid**2 - 1)/2; long2= (igrid**2 + 1)/2
      
      !!c set up sn array for lowest eigenstate
      i=0
      sn1(1)=0.0;  sn1(igrid)=0.0
      sn2(1)=0.0;  sn2(igrid)=0.0
      sn3(1)=0.0;  sn3(igrid)=0.0
      !debug  write(6,*) '----In relfac-9---->',igrid,size(sn1),size(sn2),size(sn3)
      do ix=2,igrid-1
         temp=3.14159265359*real(ix-1)/real(igrid-1)
         sn1(ix)=sqrt(2.0)*sin(temp)/sqrt(real(igrid-1))
         sn2(ix)=sn1(ix); sn3(ix)=sn1(ix)
      end do
      recipr=1.0/sqrt(real(igrid))
      if(iper(1)) then
         sn1=recipr
      end if
      if(iper(2)) then
         sn1=recipr
      end if
      if(iper(3)) then
         sn3=recipr
      end if
      
      iw=1
      do iz=1,igrid
         temp3=sn3(iz)
         do iy=1,igrid
            temp2=temp3*sn2(iy)
            do ix=1,igrid
               phimap3(iw)=temp2*sn1(ix)
               iw=iw+1
            end do
         end do
      end do
      temp=0.0
      do ix=2,icgrid-1,2
         iy=ix/2; phimap2(iy)=phimap3(ix)
         temp=temp + phimap3(ix)*phimap3(ix)
      end do
      
      if(rionst.gt.0.0) then
         do n = 2, igrid-1
            star=sta1(n)
            fin=fi1(n)
            do ix = star,fin
               temp1 = phimap2(ix) + phimap2(ix-1)
               temp2 = phimap2(ix+lat1) + phimap2(ix-lat2)
               temp3 = phimap2(ix+long1) + phimap2(ix-long2)
               phimap1(ix) = (temp1+temp2+temp3)*sf1(ix)
            end do
         end do
         
         !!c otherwise the main loop is as below:
      else
         
         do n = 2, igrid-1
            star=sta1(n)
            fin=fi1(n)
            do ix = star,fin
               temp1 = phimap2(ix) + phimap2(ix-1)
               temp2 = phimap2(ix+lat1) + phimap2(ix-lat2)
               temp3 = phimap2(ix+long1) + phimap2(ix-long2)
               phimap1(ix) =  (temp1+temp2+temp3)*sixth
            end do
         end do
      end if
      
      !!c b+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (iper(1)) then
         !!c calculating first slice
         
         ix=1+igrid*lat2
         if(rionst.gt.0.0) then
            do n=1,(igrid-3)*lat2+1
               temp1 = phimap2(ix)+phimap2(ix-1+lat1)
               temp2 = phimap2(ix+lat1)+phimap2(ix-lat2)
               rtemp2=phimap2(ix+lat1+lat1)+phimap2(ix-lat2+lat1)
               temp3 = phimap2(ix+long1)+phimap2(ix-long2)
               rtemp3=phimap2(ix+long1+lat1)+phimap2(ix-long2+lat1)
               phimap1(ix)=(temp1+.5*(temp2+temp3+rtemp2+rtemp3))*&
               &(sf1(ix)+sf1(ix+lat1))*.5
               !c now updating last slice
               phimap1(ix+lat1)=phimap1(ix)
               ix=ix+igrid
            end do
         else
            do n=1,(igrid-3)*lat2+1
               temp1 = phimap2(ix)+phimap2(ix-1+lat1)
               temp2 = phimap2(ix+lat1) + phimap2(ix-lat2)
               temp3 = phimap2(ix+long1) + phimap2(ix-long2)
               phimap1(ix)=(temp1+temp2+temp3)*sixth
               !!c now updating last slice
               phimap1(ix+lat1)=phimap1(ix)
               ix=ix+igrid
            end do
         end if
      end if
      !!c e+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      !!C$DIR NO_RECURRENCE 
      do k=1,icount2a
         ix=idpos(k)
         temp1=phimap2(ix-1)*db(1,k)+phimap2(ix)*db(2,k)
         temp2=phimap2(ix-lat2)*db(3,k)+phimap2(ix+lat1)*db(4,k)
         temp3=phimap2(ix-long2)*db(5,k)+phimap2(ix+long1)*db(6,k)
         phimap1(ix)= phimap1(ix) + temp1+temp2+temp3
      end do
      !!c       end if
      
      !!c Now reset boundary values altered in above loops.
      star=igrid*(igrid+1)/2
      fin=igrid*(igrid*(igrid-1)-2)/2
      !!C$DIR NO_RECURRENCE
      do ix=star,fin,igrid
         phimap1(ix+1)=0.0
         phimap1(ix+ihgd)=0.0
      end do

      temp=0.0
      do ix=1,(icgrid-1)/2
         temp=temp + phimap1(ix)*phimap3(2*ix-1)
      end do
      !!c if periodic boundary condition option
      !!c force periodicity using wrap around update of boundary values:
      !!c 2nd slice-->last
      !!c last-1 slice-->first
      
      !!c z periodicity
      if(iper(3)) then
         do iz = 1,(igrid-2)**2,2
            temp1=ibndz(iz)
            temp2=temp1+idif1z
            temp3=temp2+inc1za
            temp4=temp1+inc1zb
            itemp1=temp1; itemp2=temp2; itemp3=temp3; itemp4=temp4
            phimap1(itemp1)=phimap2(itemp2)
            phimap1(itemp3)=phimap2(itemp4)
         end do
      end if
      
      !!c y periodicity
      
      if(iper(2)) then
         do iy = 1,(igrid-2)**2,2
            temp1=ibndy(iy)
            temp2=temp1+idif1y
            temp3=temp2+inc1ya
            temp4=temp1+inc1yb
            itemp1=temp1; itemp2=temp2; itemp3=temp3; itemp4=temp4
            phimap1(itemp1)=phimap2(itemp2)
            phimap1(itemp3)=phimap2(itemp4)
         end do
      end if
      
      !!c Next update phimap3 using the new phimap1
      
      if(rionst.gt.0.0) then	
         do n = 2, igrid-1
            star=sta2(n)
            fin=fi2(n)
            do ix = star,fin
               temp1 = phimap1(ix) + phimap1(ix+1)
               temp2 = phimap1(ix+lat2) + phimap1(ix-lat1)
               temp3 = phimap1(ix+long2) + phimap1(ix-long1)
               phimap3(ix) =(temp1+temp2+temp3)*sf2(ix)
            end do
         end do
         
      else
         
         do n = 2, igrid-1
            star=sta2(n)
            fin=fi2(n)
            do ix = star,fin
               temp1 = phimap1(ix) + phimap1(ix+1)
               temp2 = phimap1(ix+lat2) + phimap1(ix-lat1)
               temp3 = phimap1(ix+long2) + phimap1(ix-long1)
               phimap3(ix) = (temp1+temp2+temp3)*sixth
            end do
         end do
      end if
      !!c b+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (iper(1)) then
         !!c calculating first slice
         
         ix=long2
         if(rionst.gt.0.0) then
            do n=1,long1-igrid-1
               ix=ix+igrid
               temp1 = phimap1(ix+lat1)  + phimap1(ix+1)
               temp2 = phimap1(ix+lat2)  + phimap1(ix-lat1)
               rtemp2 = phimap1(ix+lat2+lat1)  + phimap1(ix)
               temp3 = phimap1(ix+long2) + phimap1(ix-long1)
               rtemp3 = phimap1(ix+long2+lat1) + phimap1(ix-long1+lat1)
               phimap3(ix)=(temp1+.5*(temp2+temp3+rtemp2+rtemp3))*&
               &(sf2(ix)+sf2(ix+lat1))*.5
               !!c now updating last slice
               phimap3(ix+lat1)=phimap3(ix)
            end do
         else
            do n=1,long1-igrid-1
               ix=ix+igrid
               temp1 = phimap1(ix+lat1)  + phimap1(ix+1)
               temp2 = phimap1(ix+lat2)  + phimap1(ix-lat1)
               temp3 = phimap1(ix+long2) + phimap1(ix-long1)
               phimap3(ix)=(temp1+temp2+temp3)*sixth
               !!c now updating last slice
               phimap3(ix+lat1)=phimap3(ix)
            end do
         end if
      end if
      !!c e+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      do k=icount2a+1,icount2b
         ix=idpos(k)
         temp1=phimap1(ix)*db(1,k)+phimap1(ix+1)*db(2,k)
         temp2=phimap1(ix-lat1)*db(3,k)+phimap1(ix+lat2)*db(4,k)
         temp3=phimap1(ix-long1)*db(5,k)+phimap1(ix+long2)*db(6,k)
         phimap3(ix)=phimap3(ix) + temp1+temp2+temp3
      end do
      !!c reset boundary condition
      
      star=(igrid+2)/2
      iy=(igrid*(igrid+2)/2) - igrid +1
      fin=(igrid*(igrid-1)-1)/2
      ihgd2=ihgd-1
      !!C$DIR NO_RECURRENCE
      do ix=star,fin
         iy=iy+igrid
         phimap3(iy)=0.0
         phimap3(iy+ihgd2)=0.0
      end do
      
      !!c z periodicity
      
      if(iper(3)) then
         do iz = 2,(igrid-2)**2,2
            temp1=ibndz(iz);    temp2=temp1+idif2z
            temp3=temp2+inc2za; temp4=temp1+inc2zb
            itemp1=temp1;  itemp2=temp2; itemp3=temp3; itemp4=temp4
            phimap3(itemp1)=phimap1(itemp2)
            phimap3(itemp3)=phimap1(itemp4)
         end do
      end if
      
      !!c y periodicity
      
      if(iper(2)) then
         do iy = 2,(igrid-2)**2,2
            temp1=ibndy(iy);    temp2=temp1+idif2y
            temp3=temp2+inc2ya; temp4=temp1+inc2yb
            itemp1=temp1; itemp2=temp2; itemp3=temp3; itemp4=temp4
            phimap3(itemp1)=phimap1(itemp2)
            phimap3(itemp3)=phimap1(itemp4)
         end do
      end if
      
      temp=0.0
      do ix=1,(icgrid-1)/2
         temp=temp + (phimap3(ix)*phimap2(ix))
      end do
      spec=(2.0*temp)
      !!c following needed as spec exceeds 1.0 occasionally in focussing calculations
      !!c (SS May 8, 1998)
      if(spec.gt.1.0)spec=0.995
      if(verbose) then
         write(6,*) ' '

         if(ideveloper) then
            write(6,"(a,f15.10)")' gauss-seidel spectral radius is',spec
         else
            write(6,"(a,f10.4)")' gauss-seidel spectral radius is',spec
         end if
         write(6,*) ' '
      end if
      end subroutine relfac
