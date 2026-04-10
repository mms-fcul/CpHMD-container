!#####################################################################
!ATTENTION!  This file is part of crgarrmod module.
!            Do not compile separately!
!
!2011-05-30 Parameters transfered via qlog and pointers modules
!         dimension xn2(3,natom),xn1(3,natom)
!         integer extracrg
!         dimension chrgv2(natom+extracrg)
!
!         integer nqgrdtonqass(natom+extracrg),ncrgmx,realsiz
!         real atmeps(natom+extracrg),atmcrg(4,1),chgpos(3,ncrgmx)
!         dimension iatmmed(natom+nobject),iepsmp(igrid,igrid,igrid,3)
!         real medeps(0:nmedia),rx,ry,rz,cmid(3)
!         real rad3(natom)
!         integer imed,jx,jy,jz,epsdim
!         logical ipassed,verbose
!         character*15 atinf(natom)
!
!         dimension xo(3),chrgv4(natom),cqplus(3),cqmin(3)
!         integer crgatn(natom+extracrg)


!#####################################################################

      subroutine crgarr
     
      !include "pointer.h"


      integer :: epsdim
      type(coord) :: xo,rxyz
      type(int_coord) :: jxyz

      !2011-05-30 Declarations added due to IMPLICIT NONE
      integer :: ic1,ic2,i,ix,iiii,ii,imed,jx,jy,jz
      real :: chrg,rgrid
      logical :: ipassed

      !2011-05-30 First determined ic1=nqass - number of assigned 
      !charges in order to avoid re-sizing of arrays atmcrg,chgpos and
      !crgatn. 
      
      !Part1 (easy) : from molecules
      ic1=0
      do i=1,Natom
         if(abs(delphipdb(i)%chrgv4).gt.1.e-6) ic1=ic1+1
      end do

      !2011-05-30 Part 2 (more difficult) : from charge distributions
      !           running distrTOpoint without assigning values to
      !           arrays not allocated yet
      if (ndistr.gt.0) call distrTOpoint(ic1,.false.)

      nqass=ic1

      !atmcrg contains grid positions of all charges AND the charge in 
      !the 4th field atmeps6 contains 6*epsilon/epkt as a function of 
      !ic2-th charge internal
      !atmeps6  to the grid NO LONGER USED
      !nqgrdtonqass maps ic2 to ic1
      !atmeps contains epsilon/epkt as a funcion of ic1-th general 
      !charge

      allocate(atmcrg(nqass),chgpos(nqass))
      allocate(crgatn(nqass),atmeps(nqass))

      epsdim=nobject+natom+2
   
      !find charge moments for dipole approximation
      qnet=0.0;  qplus=0.0; qmin=0.0 
      cqplus=coord(0.,0.,0.); cqmin=coord(0.,0.,0.); 

      ic1=0
      do ix=1,natom
         if (abs(delphipdb(ix)%chrgv4).gt.1.e-6) then
            ic1=ic1+1
            chrg=delphipdb(ix)%chrgv4
            atmcrg(ic1)%xyz=xn2(ix)
            atmcrg(ic1)%value=chrg
            chgpos(ic1)=xn1(ix)
            crgatn(ic1)=ix
            qnet=qnet + chrg

            if (chrg.gt.0.) then
               qplus=qplus + chrg

               !2011-05-30 Using operations on coord type variables 
               !           defined in module operators_on_coordinates 
               cqplus=cqplus+(chrg*atmcrg(ic1)%xyz)
            else
               qmin=qmin + chrg
               cqmin=cqmin+(chrg*atmcrg(ic1)%xyz)
            end if
         end if
      end do
         
      if(verbose)&
            & write(6,*)"number of charges coming from molecules ",ic1

      !insert charges from charge distributions
      if (ndistr.gt.0) call distrTOpoint(ic1,.true.)

      !debug++++++++++++++++++++++++WWW
      if (.false.) then
         open(52,file='Charges.txt',form='formatted')
         do iiii=1,ic1
            write(52,*) iiii,chgpos(iiii)
         end do
         close (52)
      end if
      !end debug+++++++++++++++++++

      !assign charges for boundary conditions
      !ic1 = number of charges

      !divide by charge totals
      if (qplus.gt.1.e-6) cqplus=cqplus/qplus

      if (abs(qmin).gt.1.e-6) cqmin=cqmin/qmin

      !select those charges which will be charging the grid
      !Arrays of correct size are already allocated
      rgrid=igrid; ic2=0

      !2011-06-02 First determine correct size of arrays chgrv2 and
      !           nqgrdtongass
      do ix=1,nqass
         if ((atmcrg(ix)%xyz.vandgt.1.).and.&
            & (atmcrg(ix)%xyz.vandlt.rgrid))  ic2=ic2+1
      end do
      
      nqgrd=ic2 ; ic2=0
      allocate(chrgv2(nqgrd),nqgrdtonqass(nqgrd))

      do ix=1,nqass
         !crgatn(crg number)=atom number or natom+objectnumber or -
         !distr.number
         ii=crgatn(ix)

         if (ii.lt.0) then
            !now we have to consider charge distributions
            !in this case the distribution is not linked to any object
            !(jx,jy,jz)=coordinates of closest grid point to charge
            !(rx,ry,rz)=coordinates of the charge relatives to the 
            !current grid point

            !2011-06-02 Using operations on coord and int_coord type 
            !           variables defined  in module 
            !           operators_on_coordinates 
            jxyz=int(atmcrg(ix)%xyz+0.5)
            rxyz=atmcrg(ix)%xyz -float(jxyz)

            if (rxyz%z.gt.rxyz%x) then
               if (rxyz%z.gt.-rxyz%x) then
                  if (rxyz%z.gt.rxyz%y) then
                     if (rxyz%z.gt.-rxyz%y) then
                        imed=iepsmp(jx,jy,jz)%k
                     else
                        imed=iepsmp(jx,jy-1,jz)%j
                     end if
                  else
                     imed=iepsmp(jx,jy,jz)%j
                  end if
               else
                  if (rxyz%y.gt.rxyz%x) then
                     if (rxyz%y.gt.-rxyz%x) then
                        imed=iepsmp(jx,jy,jz)%j
                     else
                        imed=iepsmp(jx-1,jy,jz)%i
                     end if
                  else
                     imed=iepsmp(jx,jy-1,jz)%j
                  end if
               end if
            else
               if (rxyz%z.gt.-rxyz%x) then
                  if (rxyz%y.gt.rxyz%x) then
                     imed=iepsmp(jx,jy,jz)%j
                  else
                     if (rxyz%y.gt.-rxyz%x) then
                        imed=iepsmp(jx,jy,jz)%i
                     else
                        imed=iepsmp(jx,jy-1,jz)%j
                     end if
                  end if
               else
                  if (rxyz%z.gt.rxyz%y) then
                     imed=iepsmp(jx,jy-1,jz)%j
                  else
                     if (rxyz%z.gt.-rxyz%y) then
                        imed=iepsmp(jx,jy,jz)%j
                     else
                        imed=iepsmp(jx,jy,jz-1)%k
                     end if
                  end if
               end if
            end if

            imed=imed/epsdim
         else
            imed=iatmmed(ii)
         end if
      
         atmeps(ix)=medeps(imed)

         if ((atmcrg(ix)%xyz.vandgt.1.).and.&
                                & (atmcrg(ix)%xyz.vandlt.rgrid))  then
            ic2=ic2+1
            chrgv2(ic2)%xyz=atmcrg(ix)%xyz
            chrgv2(ic2)%value=atmcrg(ix)%value
            nqgrdtonqass(ic2)=ix
         end if

      end do

      ipassed=.false.
      
      do i=1,nqass
         ii=crgatn(i)
      
         if (ii.gt.0.and.ii.le.natom) then
            if (delphipdb(ii)%rad3.le.0.) then
               ipassed=.true.
               write(6,'(I4,A16,'' is charged! Radius moved from zero to''&
               &    ,f5.1)')ii,delphipdb(ii)%atinf,radpolext
               delphipdb(ii)%rad3=radpolext
            end if
         end if
      end do

      if(ipassed) write(6,*)'BE CAREFUL!! A WRONG ASSIGNMENT FOR THE &
           &RADIUS MIGHT LEAD TO INACCURATE REACTION FIELD ENERGY !!!'

      end subroutine crgarr
