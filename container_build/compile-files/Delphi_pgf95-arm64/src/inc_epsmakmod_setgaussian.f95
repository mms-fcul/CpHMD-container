!         ATTENTION!  This file is part of epsmakmod module.
!              Do not compile separately!
!==================================================================
! 2011-05-12  Parameters transfered via modules
         subroutine setgaussian(radprobe)
           implicit none
!	subroutine setout(natom,radprobe,exrad,scale,igrid,nobject
!     &  ,oldmid,nmedia,ionlymol,isitmd,debug,verbose)
!	include 'pointer.h'

!	dimension xn2(3,natom),rad3(natom),xn(3)
!	dimension iepsmp(igrid,igrid,igrid,3)
!!c b+++++++++++
!        integer nobject,imedia,objecttype,epsdim,kind,iac,nmedia
!        character*96 dataobject(nobject,2),strtmp,strtmp1
!        dimension iatmmed(natom+nobject)
!        real limgunit(3,2,nobject),vectz(3)
!        real radprobe,exrd,modul2,axdist
!        integer ioff(3,1),ismin(3),ismax(3)
!        dimension iepsmp(igrid,igrid,igrid,3)
!        integer nobject,imedia,objecttype,epsdim,kind,iac,nmedia
!        character*96 dataobject(nobject,2),strtmp,strtmp1
!        dimension iatmmed(natom+nobject)
!        real limgunit(3,2,nobject),vectz(3)
!        real radprobe,exrd,modul2,axdist
!        real rmid,oldmid(3),tmpvect(3),tmp,tmp1,dx,dy,dz,dist,zeta
!        real xa(3),xb(3),radius,modul,x2,y2,z2,tmpvect1(3),mod2,tmp2
!        real tmpvect2(3),xc(3),xd(3),xp(3),alpha,tan2,dot,modx,mody
!        real vectx(3),vecty(3),dx1,dx2,dx3,dis2min1,dis2min2
!        logical*1 idebmap(igrid,igrid,igrid)
!        logical isitmd,verbose   !  qlog module
! 2011-05-12 Some local variables
!           real :: sq1(-15:15),sq2(-15:15),sq3(-15:15)
!           real :: rad2a1(-15:15),rad2a2(-15:15),rad2a3(-15:15)
!           type(coord) :: sq(-15:15) , rad2aav(-15:15)
           type(coord) :: sq(-150:150),rad2aav(-150:150),vtemp2(-150:150)!Gaussian enlarged the values
           logical :: itobig,itest2,ipore,ionlymoldebug
           character(96) :: strtmp,strtmp1
! 2011-05-12 Non-standard integer variable, thus necessary
           integer :: epsdim, objecttype,iflag
! 2011-05-12 Non-standard real variable, thus necessary
           real :: modul,modul2, mod2,modx,mody,dentemp 
! 2011-05-12 Non-standard type of variable, thus necessary
           type(coord) :: xa,xb,xc,xd,xp,ddist,xyz2,dxyz,ddxyz,xn,vecdist
           type(coord) :: tmpvect1,tmpvect2,origin
           type(coord) :: vectx,vecty,vectz,rad2av,fxn,vtemp
           type(int_coord) :: ismin,ismax,idist,idist1,ixyz,itest,ixn,i123
           real        ::coeff,stepsize
           integer*8    ::longint
!!c b+++++++++++
!!c here radprb is not zero only if one wants to map the extended surf.
!!c that has to be done with radprb(1) if solute-solvent interface is
!!c concerned
!!c imedia = medium number for a object

!!c a non-zero entry in iepsmp indicates an atom # plus 1 (to properly treat
!!c re-entrant mid-points later (15 Aug 93)

!------------------------------------------------------------------------------------------
! 2011-05-27 Declarations added due to IMPLICIT NONE
           integer :: iboxt,iac,ibox,ii,igrdc,i,j,k,imedia,iv,ix,iy,iz,kind
           integer :: limmax,lim,n
           integer,dimension(1:6) ::inwater !mid point
           logical,dimension(1:6) ::ifinwater !mid point
           real :: alpha,dis2min2,dis2min1,dot,distsq,exrd,dist,rad,rad2
           real :: rad2a,rad4,radius,radmax2,rmid,rad5,radp2,radtest,radprobe
           real :: radp,tan2,temp,tmp,tmp1,tmp2,epstemp
           real :: peak,pi,den,distance2,sigmatime,radsq !sigma=sigma*radius. peak=4/(3*pi**0.5*sigma**3)

!------------------------------------------------------------------------------------------

           if(gaussian.eq.1.and.inhomo.eq.0.and.logs) then
              goto 1010
           endif

           pi=3.14159265359
           sigmatime=3.0
           peak=4/(3*sqrt(pi)*sigma**3)
           gepsmp2=coord(repsout,repsout,repsout)
           gepsmp=coord(0.0,0.0,0.0)
           if(verbose)write(6,*)&
              & 'Starting creating Van der Waals  Epsilon Map '
!!c e+++++++++++

           epsdim=natom+nobject+2
!debug           write(6,*)'-----setout-->',epsdim,size(iatmmed)
           iboxt=0;   radmax2=0.0; rmid=real((igrid+1)/2);  itest2=.false.
           if (natom.gt.0) then
! 2011-05-12 Label 691 removed
              do  ix=1,natom
! 2011-05-13  Chaged to derive-type array delphipdb (pointers module)
                 radmax2=max(radmax2,delphipdb(ix)%rad3)
              end do
!691	  continue
!!c b+++++++++++++++++++++++++++++
!!c this is probably the best way to do it,depending on which surf. is desired
              temp=max(radprobe,exrad)
!                print *,'Gaussian:in setout:temp,radprobe,exrad,radmax2',temp,radprobe,exrad,radmax2
!!c e++++++++++++++++++++++++++++++
!              radmax2=scale*(radmax2+temp) 
              radmax2=sigmatime*scale*(radmax2*sigma+temp) !Gaussian: 3 sigma plus temp. sigma=sigma*radius.Now radmax is 3 sigma.
              lim=1+radmax2
!                print *,'Gaussian:setout:radmax2,lim',radmax2,lim
              limmax = 1200 !Gaussian:original value:12
              itobig=.false.
              if(lim.gt.limmax) itobig=.true.
              igrdc=(2*lim+1)**3
              allocate(ioff(igrdc))
              if(.not.itobig) then
! 2011-05-12  removed goto 7878 statement
!              if(itobig) goto 7878
                 radtest= (radmax2 + 0.5*sqrt(3.0))**2
                 ibox=0
!                print*,'Gaussian:temp,radmax2,',temp,radmax2 
! 2011-05-12 Strange statement. May allocate or may not allocate array that
!            used later in the program irrespectively of itobig value, thus
!            moved array allocation before if condition
!                i_ioff= memalloc(i_ioff,4,3*igrdc)
                 do  ix=-lim,lim
                    do  iy=-lim,lim
                       do  iz=-lim,lim
                          idist=int_coord(ix,iy,iz);dist=real(idist.dot.idist)
                          ddist = dist + 0.25 + float(idist)
                          if((dist.lt.radtest).or.(ddist.vorlt.radtest)) then
                             ibox=ibox+1
                             ioff(ibox)=idist
                          else
                          end if
                       end do
                    end do
                 end do
              end if
! 2011-05-12 Labels replaced by end do and end if.
!! c +++++++++++++++++++++++++
           end if
!!c set interiors in OBJECTS   
           do ii=1,nobject
              ipore=.false.
              strtmp=dataobject(ii,1)
              write(6,*)' '
              if (strtmp(1:4).ne.'is a') then
                 strtmp1=dataobject(ii,2)
                 read(strtmp(8:10),'(I3)')objecttype
                 read(strtmp(12:14),'(I3)')imedia    

!!c completing iatmmed with imedia data
                 iatmmed(natom+ii)=imedia
                 if(verbose)write(6,*)&
                      & '(setout) object number',ii,' in medium',imedia
!!c check if object sits within limits of box and calculating ismin 
!!c and ismax accordingly
                 exrd=exrad*scale
                 temp=radprobe*scale+exrd
!  2011-05-13  Using operations on coord and int_coord type variables defined
!              in module operators_on_coordinates 
                 ismin= int(limgunit(ii)%min-temp-1.)
                 ismin=min(ismin,igrid)
                 ismin=max(ismin,1)
                 ismax= int(limgunit(ii)%max+temp+1.)
                 ismax=min(ismax,igrid)
                 ismax=max(ismax,1)

! 2011-05-13 Changed multiple IF's to select case
                 select case(objecttype)
                 case(1)!   if (objecttype.eq.1) then
!!c     dealing with a sphere
                    read(strtmp(16:80),*)kind,xb,radius
                    ipore=(kind.eq.2)
                    radius=radius*scale
                    radp=radius+exrd
                    rad2=radius*radius
                    radp2=radp*radp
!  2011-05-13  Using operations on coord and int_coord type variables defined
!              in module operators_on_coordinates 
                    xb=(xb-oldmid)*scale+rmid
                    do ix=ismin%i,ismax%i
                       do iy=ismin%j,ismax%j
                          do iz=ismin%k,ismax%k
                             ixyz=int_coord(ix,iy,iz)
                             xyz2=(float(ixyz)-xb)*(float(ixyz)-xb)
                             if(sum(xyz2).le.radp2) idebmap(ix,iy,iz)=ipore
                             xyz2%x=xyz2%x+0.25

                             ddist=sum(xyz2)+float(ixyz)-xb
                             
                             if(ddist%x.le.rad2) &
                                  &iepsmp(ix,iy,iz)%i=natom+1+ii+imedia*epsdim
                             if(ddist%y.le.rad2) &
                                  &iepsmp(ix,iy,iz)%j=natom+1+ii+imedia*epsdim
                             if(ddist%z.le.rad2) &
                                  &iepsmp(ix,iy,iz)%k=natom+1+ii+imedia*epsdim
                          end do
                       end do
                    end do


                 case(2) !if (objecttype.eq.2) then
!!c             dealing with a cylinder 
                    read(strtmp(16:80),*)kind,xa,xb,radius
                    ipore=(kind.eq.2)
                    read(strtmp1,'(5f8.3)')vectz,modul,modul2
!!c here we work in grid units
                    radius=radius*scale
                    rad2=radius*radius
                    radp=radius+exrd
                    radp2=radp*radp
                    modul=modul*scale
                    modul2=modul*modul
                    tmp=exrd*modul
!  2011-05-13  Using operations on coord and int_coord type variables defined
!              in module operators_on_coordinates 
                    xb=(xb-oldmid)*scale + rmid
                    vectz=vectz*scale
                    do ix=ismin%i,ismax%i
                       do iy=ismin%j,ismax%j
                          do iz=ismin%k,ismax%k
!!c               vectz=A-B; modul=|A-B|;tmpvect1=P-B; tmp1=(A-B)(P-B)
!!c               mod2=(P-B)**2; modul2=(A-B)**2, tmp=exrad*|A-B|.
!  2011-05-13  Using operations on coord and int_coord type variables defined
!              in module operators_on_coordinates 
                             ixyz=int_coord(ix,iy,iz)
                             tmpvect1=float(ixyz)-xb
                             tmp1=tmpvect1.dot.vectz
                             if ((tmp1.ge.-tmp).and.(tmp1.le.modul2+tmp)) then
                                mod2=tmpvect1.dot.tmpvect1
                                if((mod2-(tmp1/modul)**2).le.radp2)&
                                     & idebmap(ix,iy,iz)=ipore
                             end if

! 2011-05-13   x-offset
                             tmpvect1%x=tmpvect1%x+.5
                             tmp1=tmpvect1.dot.vectz
                             if ((tmp1.ge.0.0).and.(tmp1.le.modul2)) then
                                mod2=tmpvect1.dot.tmpvect1
                                if ((mod2-(tmp1/modul)**2).le.rad2) then
                                   iepsmp(ix,iy,iz)%i=natom+1+ii+imedia*epsdim
                                end if
                             end if
! 2011-05-13   y-offset
                             tmpvect1%x=tmpvect1%x-.5
                             tmpvect1%y=tmpvect1%y+.5
                             tmp1=tmpvect1.dot.vectz
                             if ((tmp1.ge.0.0).and.(tmp1.le.modul2)) then
                                mod2=tmpvect1.dot.tmpvect1
                                if ((mod2-(tmp1/modul)**2).le.rad2) then
                                   iepsmp(ix,iy,iz)%j=natom+1+ii+imedia*epsdim
                                end if
                             end if
! 2011-05-13   z-offset
                             tmpvect1%y=tmpvect1%y-.5
                             tmpvect1%z=tmpvect1%z+.5
                             tmp1=tmpvect1.dot.vectz
                             if ((tmp1.ge.0.0).and.(tmp1.le.modul2)) then
                                mod2=tmpvect1.dot.tmpvect1
                                if ((mod2-(tmp1/modul)**2).le.rad2) then
                                   iepsmp(ix,iy,iz)%k=natom+1+ii+imedia*epsdim
                                end if
                             end if
                          end do
                       end do
                    end do

                 case(3) ! if (objecttype.eq.3) then
!!c             dealing with a cone     
                    read(strtmp(16:80),*)kind,xa,xb,alpha 
                    ipore=(kind.eq.2)
!!c conversion degrees --> radiants
                    alpha=alpha*3.14159265359/180.
                    tan2=tan(alpha*.5)**2
                    read(strtmp1,'(5f8.3)')vectz,modul,modul2
!!c here we work in grid units
                    modul=modul*scale
                    modul2=modul*modul
!  2011-05-13  Using operations on coord and int_coord type variables defined
!              in module operators_on_coordinates 
                    xb=(xb-oldmid)*scale + rmid
                    vectz=vectz*scale
                    do ix=ismin%i,ismax%i
                       do iy=ismin%j,ismax%j
                          do iz=ismin%k,ismax%k
                             ixyz=int_coord(ix,iy,iz)
                             tmpvect1=(float(ixyz)-rmid)/scale+oldmid
!!cWWW	if ((ix.eq.17).and.(iy.eq.6).and.(iz.eq.13)) then
!!c            write(6,*) ismin, ismax
!!c		  write(6,*) tmpvect1
!!c	write(6,*) vectz, modul, modul2
!!c          end if
! 2011-05-13 Parameters transfered via module architecture
                             call distobj(tmpvect1,dist,vecdist,ii,exrad,.true.)
                             if (dist.le.0.0) idebmap(ix,iy,iz)=ipore
                             
!!c               vectz=A-B; tmpvect1=P-B; tmp1=(A-B)(P-B)
!!c               mod2=(P-B)**2; modul2=(A-B)**2.

! 2011-05-13   x-offset
!  2011-05-13  Using operations on coord and int_coord type variables defined
!              in module operators_on_coordinates 
                             tmpvect1=float(ixyz)-xb
                             tmpvect1%x=tmpvect1%x+0.5
                             tmp1=tmpvect1.dot.vectz

                             if ((tmp1.ge.0.0).and.(tmp1.le.modul2)) then
!!c                   mod2=(tmpvect1(1)**2+tmpvect1(2)**2+tmpvect1(3)**2)
!!c                   tmp2=(1+tan2)*tmp1*tmp1/modul2+tan2*(modul2-2*tmp1)
                                mod2=(tmpvect1.dot.tmpvect1)*modul2
!                                mod2=(tmpvect1(1)**2+tmpvect1(2)**2+tmpvect1(3)**2)*modul2
                                tmp2=(1+tan2)*tmp1*tmp1
                                if (mod2.le.tmp2) then
                                   iepsmp(ix,iy,iz)%i=natom+1+ii+imedia*epsdim
                                end if
                             end if
! 2011-05-13   y-offset
                             tmpvect1%x=tmpvect1%x-0.5
                             tmpvect1%y=tmpvect1%y+0.5
                             tmp1=tmpvect1.dot.vectz
                             if ((tmp1.ge.0.0).and.(tmp1.le.modul2)) then
!!c                   mod2=tmpvect1(1)**2+tmpvect1(2)**2+tmpvect1(3)**2
!!c                   tmp2=(1+tan2)*tmp1*tmp1/modul2+tan2*(modul2-2*tmp1)
                                mod2=(tmpvect1.dot.tmpvect1)*modul2
!                                mod2=(tmpvect1(1)**2+tmpvect1(2)**2+tmpvect1(3)**2)*modul2
                                tmp2=(1+tan2)*tmp1*tmp1
                                if (mod2.le.tmp2) then
                                   iepsmp(ix,iy,iz)%j=natom+1+ii+imedia*epsdim
                                end if
                             end if
! 2011-05-13   z-offset
                             tmpvect1%y=tmpvect1%y-0.5
                             tmpvect1%z=tmpvect1%z+0.5
                             tmp1=tmpvect1.dot.vectz
                             if ((tmp1.ge.0.0).and.(tmp1.le.modul2)) then
!!c                   mod2=tmpvect1(1)**2+tmpvect1(2)**2+tmpvect1(3)**2
!!c                   tmp2=(1+tan2)*tmp1*tmp1/modul2+tan2*(modul2-2*tmp1)
                                mod2=(tmpvect1.dot.tmpvect1)*modul2
!                                mod2=(tmpvect1(1)**2+tmpvect1(2)**2+tmpvect1(3)**2)*modul2
                                tmp2=(1+tan2)*tmp1*tmp1
                                if (mod2.le.tmp2) then
                                   iepsmp(ix,iy,iz)%k=natom+1+ii+imedia*epsdim
                                end if
                             end if
                             
                          end do
                       end do
                    end do
!            end if

                 case(4)  !          if (objecttype.eq.4) then
!!c             dealing with a parallelepiped
                    read(strtmp(16:80),*)kind,xa,xb,xc,xd 
                    ipore=(kind.eq.2)
                    read(strtmp1,'(12f8.3)')vectz,modul,vecty,mody,vectx, modx

!!c             conversion to axial symmetry points
!  2011-05-13  Using operations on coord and int_coord type variables defined
!              in module operators_on_coordinates 
                    xb=0.5*(xc+xb)
                    xb=(xb-oldmid)*scale + rmid
                    
                    modul=modul*scale
                    modul2=modul*modul
                    vectz=vectz*scale
                    modx=modx*scale
                    tmp1=modx*modx/2.
                    vectx=vectx*scale
                    mody=mody*scale
                    tmp2=mody*mody/2.
                    vecty=vecty*scale
                    do ix=ismin%i,ismax%i
                       do iy=ismin%j,ismax%j
                          do iz=ismin%k,ismax%k
!!c             vectz=A-B;vectx=C-D;vecty=(C+D)/2-B;tmp1=|C-D|/2   
!!c             modul2=(A-B)**2, tmp2=mody/2
!!c             dot=(P-B)(..-..); 
                             ixyz=int_coord(ix,iy,iz)
                             xp=(float(ixyz)-rmid)/scale + oldmid
! 2011-05-13 Parameters transfered via module architecture
                             call distobj(xp,dist,vecdist,ii,exrad,.true.)
                             if (dist.le.0.0) idebmap(ix,iy,iz)=ipore
!!c now xp=P-B;
! 2011-05-13  x-offset
                             xp=float(ixyz)-xb
                             xp%x=xp%x+0.5
                             dot=vectz.dot.xp
                             if ((dot.ge.0.0).and.(dot.le.modul2)) then
                                dot=vectx.dot.xp
                                if (abs(dot).le.tmp1) then
                                   dot=vecty.dot.xp
                                   if (abs(dot).le.tmp2) then
                                      iepsmp(ix,iy,iz)%i=natom+1+ii+imedia*epsdim
                                   end if
                                end if
                             end if
! 2011-05-13  y-offset
                             xp%x=xp%x-0.5
                             xp%y=xp%y+0.5
                             dot=vectz.dot.xp
                             if ((dot.ge.0.0).and.(dot.le.modul2)) then
                                dot=vectx.dot.xp
                                if (abs(dot).le.tmp1) then
                                   dot=vecty.dot.xp
                                   if (abs(dot).le.tmp2) then
                                      iepsmp(ix,iy,iz)%j=natom+1+ii+imedia*epsdim
                                   end if
                                end if
                             end if
! 2011-05-13  z-offset
                             xp%y=xp%y-0.5
                             xp%z=xp%z+0.5
                             dot=vectz.dot.xp
!                             xp(1)=ix-xb(1)
!                             xp(2)=iy-xb(2)
!                             xp(3)=iz+.5-xb(3)
!                             dot=vectz(1)*xp(1)+vectz(2)*xp(2)+vectz(3)*xp(3)
                             if ((dot.ge.0.0).and.(dot.le.modul2)) then
                                dot=vectx.dot.xp
!                                dot=vectx(1)*xp(1)+vectx(2)*xp(2)+vectx(3)*xp(3)
                                if (abs(dot).le.tmp1) then
                                   dot=vecty.dot.xp
!                                   dot=vecty(1)*xp(1)+vecty(2)*xp(2)+vecty(3)*xp(3)
                                   if (abs(dot).le.tmp2) then
                                      iepsmp(ix,iy,iz)%k=natom+1+ii+imedia*epsdim
                                   end if
                                end if
                             end if
                             
                          end do
                       end do
                    end do
                 end select

              end if
           end do
           
!!c end of setting in OBJECTS
!!c e++++++++++++++++++++++++

!!c set interiors in MOLECULES

           if(itest2.or.itobig) write(6,*)'setout method 1',itest2,itobig
           DoATOMS: do iv=1, natom
!!c restore values
             rad= delphipdb(iv)%rad3
              rad= sigmatime*sigma*rad ! 3 sigma

!  2011-05-13  Using operations on coord and int_coord type variables defined
!              in module operators_on_coordinates 
              xn=xn2(iv)
!              xn(1)=xn2(1,iv);  xn(2)=xn2(2,iv); xn(3)=xn2(3,iv)
! 2011-05-13 Removed GOTO statement

              if(rad.lt.1.e-6) then
!debug                 write(6,*)'------setout22:iv-->',iv,ixyz,iepsmp(ix,iy,iz),epsdim,iatmmed(iv)
                 cycle DoATOMS
              end if
!              if(rad.lt.1.e-6) goto 608
!!c scale radius to grid
!                print *,'rad,scale,rad*scale,radprobe',rad,scale,rad*scale,radprobe
              rad= rad*scale;  rad5= (rad + 0.5)**2;  radp = rad + exrad*scale
!!c b++++++++++++++++++++++++++++       
              rad= rad + radprobe*scale
!!c e++++++++++++++++++++++++++++rad2 is sigmatimed rad square. radsq is original rad square.
              rad4= (rad + 0.5)**2
              rad2= rad*rad
              radp2= radp*radp
              radsq=rad2/sigmatime**2


!!c set dielectric map
!!c check if sphere sits within limits of box
              itest2=.false.
!  2011-05-13  Using operations on coord and int_coord type variables defined
!              in module operators_on_coordinates 
              ismin=int(xn-radmax2-1.)
              ismax=int(xn+radmax2+1.)
              itest=ismin
              ismin=min(ismin,igrid) ; ismin=max(ismin,1)
              if(itest.vorne.ismin) itest2=.true.
              itest=ismax
              ismax=min(ismax,igrid) ; ismax=max(ismax,1)
              if(itest.vorne.ismax) itest2=.true.
!--------------Maxim Sep_15, 2011---------------------------
!              if(itest2) write(6,*)'-->',iv,xn,delphipdb(iv)%rad3,ismin,radmax2
!-----------------------------------------------------------


!!--------- slow method--------------------------
!              if(itest2.or.itobig) then
              if(itobig) then
!debug                 write(6,*)'==itest2====>',iv,nmedia,ionlymol
                 
!                    print *,'setout:slow method'
!                    print *,'itest2,itobig',itest2,itobig
! 2011-05-13 Seems redundant statement
!                 num=num+1
                 rad2a = rad2 - 0.25
                 do iz =  ismin%k,ismax%k
                    do iy =  ismin%j,ismax%j
                       do ix =  ismin%i,ismax%i
!  2011-05-13  Using operations on coord and int_coord type variables defined
!              in module operators_on_coordinates 
                          ixyz=int_coord(ix,iy,iz);  dxyz=float(ixyz)-xn
                          distsq=dxyz.dot.dxyz;      dxyz=dxyz+distsq

!!c b+++++++++++++++++++
                          if(dxyz%x.lt.rad2a) then	            
                             iepsmp(ix,iy,iz)%i=iv+1+iatmmed(iv)*epsdim
                          end if
                          if(dxyz%y.lt.rad2a) then
                             iepsmp(ix,iy,iz)%j=iv+1+iatmmed(iv)*epsdim
                          end if
                          if(dxyz%z.lt.rad2a) then
                             iepsmp(ix,iy,iz)%k=iv+1+iatmmed(iv)*epsdim
                          end if
!!c e+++++++++++++++++++
                          if(distsq.lt.radp2)  idebmap(ix,iy,iz) =.false.
                       end do
                    end do
                 end do
! 2011-05-13 Labels removed
              else


!!c faster method
!!c IT HAS PROBLEMS!!!! Walter (be careful before using 
!!c also in multidilectric case!!!.and..not.isitmd

                 rad2a = rad2 - 0.25
!  2011-05-13  Using operations on coord and int_coord type variables defined
!              in module operators_on_coordinates 
                 ixn=nint(xn);  fxn=float(ixn)-xn;  rad2av=rad2a-fxn !rad2av is never used
                 do  ix=-lim,lim
                    vtemp=real(ix)+fxn; sq(ix)=vtemp*vtemp; rad2aav(ix)=rad2a-vtemp
                    vtemp2(ix)=vtemp
                 end do
!!c *********************************************************************
!                print *,'Gaussian:nmedia,ionlymol:',nmedia,ionlymol
!!c b+++adjust inter-atom, different epsilon bgps+++04/2004 Walter+
                 if(nmedia.gt.1.and.ionlymol) then
!                    print *,'Gaussian:setout:fast 1 method'
                    do i=1,ibox
!                    do 2004 i=1,ibox
!  2011-05-14  Using operations on coord and int_coord type variables defined
!              in module operators_on_coordinates 
                       i123=ioff(i)
                       ixyz=ixn+i123
                       ix=ixyz%i ; iy=ixyz%j ; iz=ixyz%k
                       distsq = sq(i123%i)%x +sq(i123%j)%y + sq(i123%k)%z
                       write(6,*)'------->',i,ibox,i123,ixyz,distsq

!debug                       write(6,*)'-1->', iv,distsq,rad2aav(i123%i)%x,rad2aav(i123%j)%y,&
!                            & rad2aav(i123%k)%z
                       if(distsq.lt.rad2aav(i123%i)%x)  then
!                       if(distsq.lt.rad2a1(i1))  then
                          iac=mod(iepsmp(ix,iy,iz)%i,epsdim)-1	      
!                          iac=mod(iepsmp(ix,iy,iz,1),epsdim)-1	      
                          if(iac.eq.-1.or.iac.gt.natom) then
!!c occhio! non so cosa fa con i pori!!
                             iepsmp(ix,iy,iz)%i=iv+1+iatmmed(iv)*epsdim
                          else
!  2011-05-14  Using operations on coord and int_coord type variables defined
!              in module operators_on_coordinates 
                             ddxyz=float(ixyz)-xn ; ddxyz%x= ddxyz%x+0.5
                             dis2min1=ddxyz.dot.ddxyz - rad2
!                             dx1=ix+0.5-xn(1)
!                             dx2=iy-xn(2)
!                             dx3=iz-xn(3)
!                             dis2min1=dx1*dx1+dx2*dx2+dx3*dx3-rad2
   
                             ddxyz=float(ixyz)-xn2(iac) ; ddxyz%x= ddxyz%x+0.5
                             dis2min2=ddxyz.dot.ddxyz -(delphipdb(iac)%rad3*scale)**2 
!                             dx1=ix-xn2(1,iac)+0.5
!                             dx2=iy-xn2(2,iac)
!                             dx3=iz-xn2(3,iac)
!                             dis2min2=dx1*dx1+dx2*dx2+dx3*dx3-(rad3(iac)*scale)**2

                             if (dis2min2.gt.dis2min1) iac=iv
                             iepsmp(ix,iy,iz)%i=iac+1+iatmmed(iac)*epsdim
                          endif
!!c         write(6,*) "atomo:",iv," mezzo:",iepsmp(ix,iy,iz,1,2),"asse 1"
                       end if
                       if(distsq.lt.rad2aav(i123%j)%y) then
!                       if(distsq.lt.rad2a2(i2)) then
                          iac=mod(iepsmp(ix,iy,iz)%j,epsdim)-1	      
!                          iac=mod(iepsmp(ix,iy,iz,2),epsdim)-1	      
                          if(iac.eq.-1.or.iac.gt.natom) then
!!c occhio! non so cosa fa con i pori!!
                             iepsmp(ix,iy,iz)%j=iv+1+iatmmed(iv)*epsdim
!                             iepsmp(ix,iy,iz,2)=iv+1+iatmmed(iv)*epsdim
                          else
!  2011-05-14  Using operations on coord and int_coord type variables defined
!              in module operators_on_coordinates 
                             ddxyz=float(ixyz)-xn ; ddxyz%y= ddxyz%y+0.5
                             dis2min1=ddxyz.dot.ddxyz - rad2
!                             dx1=ix-xn(1)
!                             dx2=iy+0.5-xn(2)
!                             dx3=iz-xn(3)
!                             dis2min1=dx1*dx1+dx2*dx2+dx3*dx3-rad2
   
                             ddxyz=float(ixyz)-xn2(iac) ; ddxyz%y= ddxyz%y+0.5
                             dis2min2=ddxyz.dot.ddxyz -(delphipdb(iac)%rad3*scale)**2 
!                             dx1=ix-xn2(1,iac)
!                             dx2=iy+0.5-xn2(2,iac)
!                             dx3=iz-xn2(3,iac)
!                             dis2min2=dx1*dx1+dx2*dx2+dx3*dx3-(rad3(iac)*scale)**2
                             if (dis2min2.gt.dis2min1) iac=iv
                             iepsmp(ix,iy,iz)%j=iac+1+iatmmed(iac)*epsdim
!                             iepsmp(ix,iy,iz,2)=iac+1+iatmmed(iac)*epsdim
                          endif
!!c         write(6,*) "atomo:",iv," mezzo:",iepsmp(ix,iy,iz,2,2),"asse 2"
                       end if

                       if(distsq.lt.rad2aav(i123%k)%z) then
!                       if(distsq.lt.rad2a3(i3)) then
                          iac=mod(iepsmp(ix,iy,iz)%k,epsdim)-1	      
!                          iac=mod(iepsmp(ix,iy,iz,3),epsdim)-1	      
                          if(iac.eq.-1.or.iac.gt.natom) then
!!c occhio! non so cosa fa con i pori!!
                             iepsmp(ix,iy,iz)%k=iv+1+iatmmed(iv)*epsdim
!                             iepsmp(ix,iy,iz,3)=iv+1+iatmmed(iv)*epsdim
                          else
!  2011-05-14  Using operations on coord and int_coord type variables defined
!              in module operators_on_coordinates 
                             ddxyz=float(ixyz)-xn ; ddxyz%z= ddxyz%z+0.5
                             dis2min1=ddxyz.dot.ddxyz - rad2 
!                             dx1=ix-xn(1)
!                             dx2=iy-xn(2)
!                             dx3=iz+0.5-xn(3)
!                             dis2min1=dx1*dx1+dx2*dx2+dx3*dx3-rad2
 
                             ddxyz=float(ixyz)-xn2(iac) ; ddxyz%z= ddxyz%z+0.5
                             dis2min2=ddxyz.dot.ddxyz -(delphipdb(iac)%rad3*scale)**2 
!                             dx1=ix-xn2(1,iac)
!                             dx2=iy-xn2(2,iac)
!                             dx3=iz+0.5-xn2(3,iac)
!                             dis2min2=dx1*dx1+dx2*dx2+dx3*dx3-(rad3(iac)*scale)**2
                             if (dis2min2.gt.dis2min1) iac=iv
                             iepsmp(ix,iy,iz)%k=iac+1+iatmmed(iac)*epsdim
!                             iepsmp(ix,iy,iz,3)=iac+1+iatmmed(iac)*epsdim
                          endif
!!c         write(6,*) "atomo:",iv," mezzo:",iepsmp(ix,iy,iz,3,2),"asse 3"
                       end if

                       if(distsq.lt.radp2)   idebmap(ix,iy,iz)=.false.
                    end do
!2004                   continue
                 else 
!!c ****************************************************




!                    print *,'Gaussian:setout:fast 2 method'
!!C$DIR NO_RECURRENCE
!                    write(6,*)'Gaussian:iv,nmedia,ionlymol,ibox:',iv,nmedia,ionlymol,ibox
                    do  i=1,ibox

        
!                    do 9024 i=1,ibox
!  2011-05-14  Using operations on coord and int_coord type variables defined
!              in module operators_on_coordinates 
                       i123=ioff(i); ixyz=ixn+i123
                       ix=ixyz%i ; iy=ixyz%j ; iz=ixyz%k
!Gaussian:Avoid reaching grids out side (1,igrid)
                       if(ix.lt.igrid.and.ix.gt.1.and.iy.lt.igrid.and.&
                        & iy.gt.1.and.iz.lt.igrid.and.iz.gt.1) then 
                       distsq = sq(i123%i)%x +sq(i123%j)%y + sq(i123%k)%z
!                       i1= ioff(1,i); i2= ioff(2,i);  i3= ioff(3,i)
!                       ix=ixn1+ i1; iy=iyn1+ i2;   iz=izn1+ i3
!                       distsq = sq1(i1) +sq2(i2) + sq3(i3)
!!c b+++++++++++++++++++
!debug                       write(6,*)'-2->', iv,distsq,rad2aav(i123%i)%x,rad2aav(i123%j)%y,&
!                            & rad2aav(i123%k)%z
!                       write(6,*)'-2a->', i,ibox,i123,ixn
!                       write(6,*)'-2b->', i,ibox,sq(i123%i)%x,sq(i123%j)%y,sq(i123%k)%z
!                       print *,'Gaussian:distsq,rad2aav(i123%i)%x:',distsq,rad2aav(i123%i)%x
 
                       if(distsq.lt.rad2aav(i123%i)%x)  then
!                          iepsmp(ix,iy,iz)%i=iv+1+iatmmed(iv)*epsdim
!                          gepsmp(ix,iy,iz)%x=repsin
                          distance2=distsq+0.25+vtemp2(i123%i)%x
!                          den=peak*exp(-(distance2/(sigma**2*radsq)))
                          den=exp(-(distance2/(sigma**2*radsq)))
                          gepsmp(ix,iy,iz)%x=1-(1-gepsmp(ix,iy,iz)%x)*(1-den)

!               gepsmp2(ix,iy,iz)%x=gepsmp(ix,iy,iz)%x*repsin+(1-gepsmp(ix,iy,iz)%x)*repsout
!                       if(distsq.lt.rad2a1(i1))  then
!                          iepsmp(ix,iy,iz,1)=iv+1+iatmmed(iv)*epsdim
!         write(6,*) "atomo:",iv," mezzo:",iepsmp(ix,iy,iz)%i,"asse 1"
!!c         write(6,*) "atomo:",iv," mezzo:",iepsmp(ix,iy,iz,1,2),"asse 1"
                       end if
 
                       if(distsq.lt.rad2aav(i123%j)%y) then
!                          iepsmp(ix,iy,iz)%j=iv+1+iatmmed(iv)*epsdim
!                          gepsmp(ix,iy,iz)%y=repsin
                          distance2=distsq+0.25+vtemp2(i123%j)%y
                          den=exp(-(distance2/(sigma**2*radsq)))

!                           if(ix.eq.100.and.iy.eq.100.and.iz.eq.100) &
!                          den=peak*exp(-(distance2/(sigma**2*radsq)))
                          gepsmp(ix,iy,iz)%y=1-(1-gepsmp(ix,iy,iz)%y)*(1-den)
!               gepsmp2(ix,iy,iz)%y=gepsmp(ix,iy,iz)%y*repsin+(1-gepsmp(ix,iy,iz)%y)*repsout

!                          iepsmp(ix,iy,iz,2)=iv+1+iatmmed(iv)*epsdim
!!c         write(6,*) "atomo:",iv," mezzo:",iepsmp(ix,iy,iz,2,2),"asse 2"
                       end if
 
                       if(distsq.lt.rad2aav(i123%k)%z) then
!                          iepsmp(ix,iy,iz)%k=iv+1+iatmmed(iv)*epsdim
!                          gepsmp(ix,iy,iz)%z=repsin
                          distance2=distsq+0.25+vtemp2(i123%k)%z
                          den=exp(-(distance2/(sigma**2*radsq))) !Gaussian: make density at center is 1.
!                          den=peak*exp(-(distance2/(sigma**2*radsq)))
                          gepsmp(ix,iy,iz)%z=1-(1-gepsmp(ix,iy,iz)%z)*(1-den)
!               gepsmp2(ix,iy,iz)%z=gepsmp(ix,iy,iz)%z*repsin+(1-gepsmp(ix,iy,iz)%z)*repsout
!                       if(distsq.lt.rad2a1(i1))  then
!
!                       if(distsq.lt.rad2a3(i3)) then
!                          iepsmp(ix,iy,iz,3)=iv+1+iatmmed(iv)*epsdim
!!c         write(6,*) "atomo:",iv," mezzo:",iepsmp(ix,iy,iz,3,2),"asse 3"
                       end if

!!c e++++++++++++++++++++
!                       if(distsq.lt.radp2)   idebmap(ix,iy,iz)=.false.
!Gaussian change method of generating idebmap

                    endif  !for detecting outside cube
                    end do
!!9024	      continue
                 end if
              end if
! 2011-05-13 Removed label 608
! 608           continue
!        write(6,*)'------setout22:iv-->',iv,ixyz,iepsmp(ix,iy,iz),epsdim,iatmmed(iv)

!!c end do of atoms

           end do DoATOMS


!           cutoff=0.9
!           print *,"cutoff,sigma:",cutoff,sigma
           if(.false.)then !cutoff: make the top flat
              do ix=1,igrid
                 do iy=1,igrid
                    do iz=1,igrid
                       gepsmp(ix,iy,iz)%x=gepsmp(ix,iy,iz)%x/cutoff
                       gepsmp(ix,iy,iz)%y=gepsmp(ix,iy,iz)%y/cutoff
                       gepsmp(ix,iy,iz)%z=gepsmp(ix,iy,iz)%z/cutoff
                       if(gepsmp(ix,iy,iz)%x.gt.1.0) gepsmp(ix,iy,iz)%x=1.0
                       if(gepsmp(ix,iy,iz)%y.gt.1.0) gepsmp(ix,iy,iz)%y=1.0
                       if(gepsmp(ix,iy,iz)%z.gt.1.0) gepsmp(ix,iy,iz)%z=1.0
                    enddo
                 enddo
              enddo
           endif

1010       continue
!           print *,"debug: srfcut,inhomo:",srfcut,inhomo
           do ix=1,igrid
              do iy=1,igrid
                 do iz=1,igrid
                   gepsmp2(ix,iy,iz)%x=gepsmp(ix,iy,iz)%x*repsin+(1-gepsmp(ix,iy,iz)%x)*repsout
                   gepsmp2(ix,iy,iz)%y=gepsmp(ix,iy,iz)%y*repsin+(1-gepsmp(ix,iy,iz)%y)*repsout
                   gepsmp2(ix,iy,iz)%z=gepsmp(ix,iy,iz)%z*repsin+(1-gepsmp(ix,iy,iz)%z)*repsout
!################### for set epsout in protein larger than in water:##########
                   if(gepsmp(ix,iy,iz)%x.lt.0.02)then
                       gepsmp2(ix,iy,iz)%x=80.0
                   endif
                   if(gepsmp(ix,iy,iz)%y.lt.0.02)then
                       gepsmp2(ix,iy,iz)%y=80.0
                   endif
                   if(gepsmp(ix,iy,iz)%z.lt.0.02)then
                       gepsmp2(ix,iy,iz)%z=80.0
                   endif
!################### end for this epsout in protein different than in water######
                 enddo
              enddo
           enddo


           if(inhomo.eq.1)then !reduce epsilon out side protein
!              epstemp=repsin*(1-srfcut)+repsout*srfcut
              epstemp=srfcut
              if(epstemp.lt.repsin) then
                 print*,'srfcut is lower than epsin.'
                 stop
              endif
!              print*,"epstemp",epstemp
              do ix=1,igrid
                 do iy=1,igrid
                    do iz=1,igrid
                      if(gepsmp2(ix,iy,iz)%x.gt.epstemp)then
                         gepsmp2(ix,iy,iz)%x=repsin
                      end if
                      if(gepsmp2(ix,iy,iz)%y.gt.epstemp)then
                         gepsmp2(ix,iy,iz)%y=repsin
                      end if
                      if(gepsmp2(ix,iy,iz)%z.gt.epstemp)then
                         gepsmp2(ix,iy,iz)%z=repsin
                      end if
                    enddo
                 enddo
              enddo
           endif
 
           ibnum=0;longint=0;dentemp=0.1
!           ifinwater=.true.
           do i=2,igrid
              do j=2,igrid
                do k=2,igrid
                   if(gepsmp(i,j,k)%x.gt.dentemp.or.gepsmp(i,j,k)%y.gt.dentemp &
                   & .or.gepsmp(i,j,k)%z.gt.dentemp.or.gepsmp(i-1,j,k)%x.gt.dentemp &
                   & .or.gepsmp(i,j-1,k)%y.gt.dentemp.or.gepsmp(i,j,k-1)%z.gt.dentemp)then
                   
                         longint=longint+1
                         idebmap(i,j,k)=.false. !gaussian change method of generating idebmap
                  endif
                enddo
              enddo
           enddo
           ibnum=longint
           allocate(ibgrd(ibnum))

           n=0
           do i=2,igrid
              do j=2,igrid
                do k=2,igrid

                   if(gepsmp(i,j,k)%x.gt.dentemp.or.gepsmp(i,j,k)%y.gt.dentemp &
                   & .or.gepsmp(i,j,k)%z.gt.dentemp.or.gepsmp(i-1,j,k)%x.gt.dentemp &
                   & .or.gepsmp(i,j-1,k)%y.gt.dentemp.or.gepsmp(i,j,k-1)%z.gt.dentemp)then


                     n=n+1
                     ibgrd(n)%i=i
                     ibgrd(n)%j=j
                     ibgrd(n)%k=k
                  endif
                enddo
              enddo
           enddo
 
! Array deallocation is made by ordinary F95 statement
           if(allocated(ioff)) deallocate(ioff)
!           i_ioff= memalloc(i_ioff,0,0)
!!c b+++++debug+++++++++++++++++++++++++++++++
           if (.false.) then !Gaussian
              open(52,file='linlieps_k',form='formatted')
              write (52,*)igrid,igrid*igrid
              do ix=1,igrid
                 do iy=1,igrid
                    do iz=1,igrid
                    write (52,*)ix,iy,iz,iepsmp(ix,iy,iz)%k
                    end do
                 end do
              end do
              close (52)

           end if


           if (.false.) then !debug: output density and eps
               open(14,file="cube_density", form='formatted')
              
              write(14,*)'qdiffxs4 with an improved surfacing routine'
              write(14,*) 'Gaussian cube format phimap'
              coeff=0.5291772108
              stepsize=1.0/scale
              origin=(oldmid-stepsize*(igrid-1)/2)/coeff
!             do i=1,3
!                origin(i)=(oldmid(i)-stepsize*(igrid-1)/2)/coeff
!             enddo
              write(14,'(i5,3f12.6)') 1, origin
!              write(14,'(i5,3f12.6)') 1, (origin(i),i=1,3)
              write(14,'(i5,3f12.6)') igrid, stepsize/coeff,0.0,0.0
              write(14,'(i5,3f12.6)') igrid, 0.0,stepsize/coeff,0.0
              write(14,'(i5,3f12.6)') igrid, 0.0,0.0,stepsize/coeff
              write(14,'(i5,4f12.6)') 1,0.0,0.0,0.0,0.0
              
              do ix = 1,igrid
                 do iy = 1,igrid
!                    do iz=1,igrid
                    write(14,'(6E13.5)')(gepsmp(ix,iy,iz)%z,iz=1,igrid)
!              write(14,*),ix,iy,iz,gepsmp(ix,iy,iz)%x,gepsmp(ix,iy,iz)%y,gepsmp(ix,iy,iz)%z

                 end do
              end do
              close(14)

              open(14,file="cube_eps", form='formatted')
              
              write(14,*)'qdiffxs4 with an improved surfacing routine'
              write(14,*) 'Gaussian cube format phimap'
              coeff=0.5291772108
              stepsize=1.0/scale
              origin=(oldmid-stepsize*(igrid-1)/2)/coeff
!             do i=1,3
!                origin(i)=(oldmid(i)-stepsize*(igrid-1)/2)/coeff
!             enddo
              write(14,'(i5,3f12.6)') 1, origin
              write(14,'(i5,3f12.6)') igrid, stepsize/coeff,0.0,0.0
              write(14,'(i5,3f12.6)') igrid, 0.0,stepsize/coeff,0.0
              write(14,'(i5,3f12.6)') igrid, 0.0,0.0,stepsize/coeff
              write(14,'(i5,4f12.6)') 1,0.0,0.0,0.0,0.0
              
              do ix = 1,igrid
                 do iy = 1,igrid
                    write(14,'(6E13.5)')(gepsmp2(ix,iy,iz)%z,iz=1,igrid)
                 end do
              end do
              close(14)
           endif
!!c e++++++++++++++++++++++++++++++++++++++++++
           if(verbose)write(6,*)'Ending creating Van der Waals  Epsilon Map '

         return
         end subroutine setgaussian

