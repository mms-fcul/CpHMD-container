      !  ATTENTION!  This file is part of setcrgmod module.
      !==================================================================
      ! 2011-06-02 Parameters transfred via qlog and pointers modules
      subroutine mkdbsf
      !	subroutine mkdbsf(ibnum,nsp,dbval,icount2a,icount2b,sfd,natom,
      ! 2011-06-02 Arrays are accessible via pointers module
      !!c b+++++++++++++++++++++++
      real :: denom,vecttemp(6)
      integer :: epsdim
      !!c e+++++++++++++++++++++++
      
      ! 2011-06-03 Fixed-size arrays dbval and sfd and variables
      integer deb,it(6)
      !--------------------------------------------------------------
      ! 2011-06-03 Declarations added due to IMPLICIT NONE
      real :: debfct,dbs,sixeps,sixth,sixsalt,temp
      integer :: ibnum1,ibnum2,ibnum3,idbs,icgrid,isgrid,ieps
      integer :: i,j,k,l,iv,iz,iy,ix,iw,itmp,itemp
      
      ! 2011-06-03 Arrays iepsv and phimap declared as allocatable
      allocate(iepsv(nsp),phimap3(ngp))
      
      debfct = epsout/(deblen*scale)**2
      sixeps=epsout*6.0; sixth=1.0/6.0
      sixsalt=sixth*((1/(1+debfct/sixeps))-1.0)
      isgrid=igrid*igrid
      !!c initialisation of idpos+++++++++++++++
      epsdim=natom+nobject+2
      !!c +++++++++++++++++++++++++++++++++++++
      ibnum1=0; ibnum2=nsp/2; idbs=0; dbs=0.
      
      if(gaussian.eq.1) then
         idbwrt=.false. !LinLi
      endif

      if(idbwrt) then
         open(13,file=dbnam(:dblen))
         write(13,*) "DELPHI DB FILE"
         write(13,*) "FORMAT NUMBER=1"
         write(13,*) "NUMBER OF BOUNDARY POINTS= ",ibnum
      end if
      
      do  ix=1,ibnum
         i=ibgrd(ix)%i ; j=ibgrd(ix)%j; k=ibgrd(ix)%k
         if (idirectalg.eq.0) then
            ! 2011-06-03 Changed to array operations
            it=0
            !c b+++++++++++++++++
            if((iepsmp(i,j,k)%i/epsdim).ne.0) it(1)=1
            if((iepsmp(i,j,k)%j/epsdim).ne.0) it(2)=1
            if((iepsmp(i,j,k)%k/epsdim).ne.0) it(3)=1
            if((iepsmp(i-1,j,k)%i/epsdim).ne.0) it(4)=1
            if((iepsmp(i,j-1,k)%j/epsdim).ne.0) it(5)=1
            if((iepsmp(i,j,k-1)%k/epsdim).ne.0) it(6)=1
            ieps=it(1)+it(2)+it(3)+it(4)+it(5)+it(6)
         else
            ieps=0 ; temp=0.
            if(gaussian.eq.0)then !LinLi: old style
!               print *,'LinLi:Non-Gaussian in mkdbsf'
               itmp=iepsmp(i,j,k)%i/epsdim
               temp=temp+medeps(itmp)
               vecttemp(1)=medeps(itmp)
               itmp=iepsmp(i,j,k)%j/epsdim
               temp=temp+medeps(itmp)
               vecttemp(2)=medeps(itmp)
               itmp=iepsmp(i,j,k)%k/epsdim
               temp=temp+medeps(itmp)
               vecttemp(3)=medeps(itmp)
               itmp=iepsmp(i-1,j,k)%i/epsdim
               temp=temp+medeps(itmp)
               vecttemp(4)=medeps(itmp)
               itmp=iepsmp(i,j-1,k)%j/epsdim
               temp=temp+medeps(itmp)
               vecttemp(5)=medeps(itmp)
               itmp=iepsmp(i,j,k-1)%k/epsdim
               temp=temp+medeps(itmp)
               vecttemp(6)=medeps(itmp)
            else !LinLi: this is Gaussian
!               print *,'LinLi:Gaussian in mkdbsf'
               temp=temp+gepsmp2(i,j,k)%x
               vecttemp(1)=gepsmp2(i,j,k)%x

               temp=temp+gepsmp2(i,j,k)%y
               vecttemp(2)=gepsmp2(i,j,k)%y

               temp=temp+gepsmp2(i,j,k)%z
               vecttemp(3)=gepsmp2(i,j,k)%z

               temp=temp+gepsmp2(i-1,j,k)%x
               vecttemp(4)=gepsmp2(i-1,j,k)%x

               temp=temp+gepsmp2(i,j-1,k)%y
               vecttemp(5)=gepsmp2(i,j-1,k)%y

               temp=temp+gepsmp2(i,j,k-1)%z
!               temp=temp+gepsmp(i,j,k-1)%z !mistake
               vecttemp(6)=gepsmp2(i,j,k-1)%z
            endif

            !!c e+++++++++++++++++
         end if
         deb=0
         if (idebmap(i,j,k)) deb=1
         if(deb.eq.1) idbs=idbs+1
         
         iw=isgrid*(k-1) + igrid*(j-1) + i
         iv=(iw+1)/2
         if(iw.ne.(2*iv)) then
            ibnum1=ibnum1+1 ; ibnum3=ibnum1
         else
            ibnum2=ibnum2+1 ; ibnum3=ibnum2
         end if
         
         idpos(ibnum3) = iv ;  iepsv(ibnum3) = ieps
         
         if (idirectalg.eq.0) then
            db(1,ibnum3)=dbval(it(4),ieps,deb)
            db(2,ibnum3)=dbval(it(1),ieps,deb)
            db(3,ibnum3)=dbval(it(5),ieps,deb)
            db(4,ibnum3)=dbval(it(2),ieps,deb)
            db(5,ibnum3)=dbval(it(6),ieps,deb)
            db(6,ibnum3)=dbval(it(3),ieps,deb)
         else
            !!c b+++++++++++++++
            denom=temp+deb*debfct
            if (rionst.eq.0.) then
               db(1,ibnum3)=vecttemp(4)/denom -sixth
               db(2,ibnum3)=vecttemp(1)/denom -sixth
               db(3,ibnum3)=vecttemp(5)/denom -sixth
               db(4,ibnum3)=vecttemp(2)/denom -sixth
               db(5,ibnum3)=vecttemp(6)/denom -sixth
               db(6,ibnum3)=vecttemp(3)/denom -sixth
            else
               db(1,ibnum3)=vecttemp(4)/denom
               db(2,ibnum3)=vecttemp(1)/denom
               db(3,ibnum3)=vecttemp(5)/denom
               db(4,ibnum3)=vecttemp(2)/denom
               db(5,ibnum3)=vecttemp(6)/denom
               db(6,ibnum3)=vecttemp(3)/denom
            end if
            !!c e+++++++++++++++          
         end if
         if(idbwrt) write(13,'(i2,1x,i2,1x,i2,6f12.8)') &
         & i,j,k,(db(l,ibnum3),l=1,6)
         !400         format(i2,x,i2,x,i2,6f12.8)
         dbs=dbs+db(1,ibnum3)
      end do
      if(idbwrt) close(13)
      if(verbose) write(6,*) &
      &"no. dielectric boundary points in salt = ",idbs
      
      !!c realign idpos and db,compressing to contingous space
      icount2a=ibnum1
      icount2b=icount2a+ibnum2-(nsp/2)
      
      itemp=(nsp/2)
      do  ix=icount2a+1,icount2b
         itemp=itemp+1
         idpos(ix)=idpos(itemp)
         iepsv(ix)=iepsv(itemp)
         db(1,ix)=db(1,itemp)
         db(2,ix)=db(2,itemp)
         db(3,ix)=db(3,itemp)
         db(4,ix)=db(4,itemp)
         db(5,ix)=db(5,itemp)
         db(6,ix)=db(6,itemp)
      end do
      
      !!c set saltmaps 1 and 2
      !!c NB phimap3 used as a dummy flat 65**3 array
      if(rionst.gt.0.) then
         iw=1
         do iz=1,igrid
            do iy=1,igrid
               do ix=1,igrid
                  deb=0
                  if (idebmap(ix,iy,iz)) deb=1
                  phimap3(iw)=sixth + deb*sixsalt
                  iw=iw+1
               end do
            end do
         end do
         
         iy=0
         icgrid=igrid*igrid*igrid
         sf1((icgrid+1)/2)=phimap3(icgrid)
         
         do ix=1,icgrid-2,2
            iy=iy+1
            sf1(iy)=phimap3(ix)
            sf2(iy)=phimap3(ix+1)
         end do
         
         do ix=1,icount2a
            !!c b+++++++++++++++
            if (idirectalg.ne.0) then
               sf1(idpos(ix))= 0.0
            else
               !!c e+++++++++++++++
               i=1
               if(sf1(idpos(ix)).eq.sixth) i=0
               sf1(idpos(ix))=sfd(iepsv(ix),i)
            end if
         end do
         
         do ix=icount2a+1,icount2b
            !!c b+++++++++++++++
            if (idirectalg.ne.0) then
               sf2(idpos(ix))= 0.0
            else
               !!c e+++++++++++++++
               i=1
               if(sf2(idpos(ix)).eq.sixth) i=0
               sf2(idpos(ix))=sfd(iepsv(ix),i)
            end if
         end do
      end if
      if(allocated(iepsv)) deallocate(iepsv)
      if(allocated(phimap3)) deallocate(phimap3)
      end subroutine mkdbsf
