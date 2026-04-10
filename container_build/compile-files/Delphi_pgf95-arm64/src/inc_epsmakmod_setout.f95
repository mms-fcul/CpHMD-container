!#####################################################################
!ATTENTION! This file is part of epsmakmod module.
!           Do not compile separately!
!
! 2011-05-12  Parameters transfered via modules
!#####################################################################

      subroutine setout(radprobe)
      
      implicit none

      !2011-05-12 Some local variables
      type(coord) :: sq(-15:15) , rad2aav(-15:15)
      logical :: itobig,itest2,ipore,ionlymoldebug
      character(96) :: strtmp,strtmp1
      !2011-05-12 Non-standard integer variable, thus necessary
      integer :: epsdim, objecttype
      !2011-05-12 Non-standard real variable, thus necessary
      real :: modul,modul2, mod2,modx,mody         
      !2011-05-12 Non-standard type of variable, thus necessary
      type(coord) :: xa,xb,xc,xd,xp,ddist,xyz2,dxyz,ddxyz,xn,vecdist
      type(coord) :: tmpvect1,tmpvect2
      type(coord) :: vectx,vecty,vectz,rad2av,fxn,vtemp
      type(int_coord) :: ismin,ismax,idist,idist1,ixyz,itest,ixn,i123
      !here radprb is not zero only if one wants to map the extended 
      !surf. that has to be done with radprb(1) if solute-solvent 
      !interface is concerned
      !imedia = medium number for a object
      !a non-zero entry in iepsmp indicates an atom # plus 1 (to 
      !properly treat re-entrant mid-points later (15 Aug 93)
      !2011-05-27 Declarations added due to IMPLICIT NONE
      integer :: iboxt,iac,ibox,ii,igrdc,i,imedia,iv,ix,iy,iz,kind
      integer :: limmax,lim
      real :: alpha,dis2min2,dis2min1,dot,distsq,exrd,dist,rad,rad2
      real :: rad2a,rad4,radius,radmax2,rmid,rad5,radp2,radtest,radprobe
      real :: radp,tan2,temp,tmp,tmp1,tmp2

      if(verbose)write(6,*)&
                     & 'Starting creating Van der Waals  Epsilon Map '

      epsdim=natom+nobject+2
      iboxt=0; radmax2=0.0; rmid=real((igrid+1)/2); itest2=.false.

      if (natom.gt.0) then
         !2011-05-12 Label 691 removed
         do ix=1,natom
            !2011-05-13 Chaged to derive-type array delphipdb (pointers 
            !module)
            radmax2=max(radmax2,delphipdb(ix)%rad3)
         end do

         !this is probably the best way to do it,depending on which 
         !surf. is desired
         temp=max(radprobe,exrad)

         radmax2=scale*(radmax2+temp);  lim=1+radmax2
         limmax = 12
         itobig=.false.
         if(lim.gt.limmax) itobig=.true.
         igrdc=(2*lim+1)**3
         allocate(ioff(igrdc))

         if (.not.itobig) then
            !2011-05-12 removed goto 7878 statement
            radtest= (radmax2 + 0.5*sqrt(3.0))**2
            ibox=0
                 
            !2011-05-12 Strange statement. May allocate or may not 
            !allocate array that used later in the program 
            !irrespectively of itobig value, thus moved array 
            !allocation before if condition
            do ix=-lim,lim
               do iy=-lim,lim
                  do iz=-lim,lim
                     !2011-05-13 Replaced by faster operation
                     !2011-05-13  Using operations on coord and 
                     !int_coord type variables defined in module 
                     !operators_on_coordinates 
                     idist=int_coord(ix,iy,iz)
                     dist=real(idist.dot.idist)
                     ddist = dist + 0.25 + float(idist)
                     
                     if ((dist.lt.radtest).or.(ddist.vorlt.radtest))then
                        ibox=ibox+1
                        ioff(ibox)=idist
                     end if
                  end do
               end do
            end do
         end if
      end if

      !set interiors in OBJECTS   
      do ii=1,nobject
         ipore=.false.
         strtmp=dataobject(ii,1)
         write(6,*)' '
      
         if (strtmp(1:4).ne.'is a') then
            strtmp1=dataobject(ii,2)
            read(strtmp(8:10),'(I3)')objecttype
            read(strtmp(12:14),'(I3)')imedia     
          
            !completing iatmmed with imedia data
            iatmmed(natom+ii)=imedia
            if(verbose)write(6,*)&
                     & '(setout) object number',ii,' in medium',imedia

            !check if object sits within limits of box and calculating 
            !ismin and ismax accordingly
            exrd=exrad*scale
            temp=radprobe*scale+exrd
 
            !2011-05-13  Using operations on coord and int_coord type 
            !variables defined in module operators_on_coordinates 
            ismin=int(limgunit(ii)%min-temp-1.)
            ismin=min(ismin,igrid)
            ismin=max(ismin,1)

            ismax=int(limgunit(ii)%max+temp+1.)
            ismax=min(ismax,igrid)
            ismax=max(ismax,1)
                 
            !2011-05-13 Changed multiple IF's to select case
            select case(objecttype)
            case(1) !if (objecttype.eq.1) then
                    !dealing with a sphere
               read(strtmp(16:80),*)kind,xb,radius
               ipore=(kind.eq.2)
               radius=radius*scale
               radp=radius+exrd
               rad2=radius*radius
               radp2=radp*radp

               !2011-05-13  Using operations on coord and int_coord 
               !type variables defined in module 
               !operators_on_coordinates 
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
                    !dealing with a cylinder 
               read(strtmp(16:80),*)kind,xa,xb,radius
               ipore=(kind.eq.2)
               read(strtmp1,'(5f8.3)')vectz,modul,modul2

               !here we work in grid units
               radius=radius*scale
               rad2=radius*radius
               radp=radius+exrd
               radp2=radp*radp
               modul=modul*scale
               modul2=modul*modul
               tmp=exrd*modul

               !2011-05-13  Using operations on coord and int_coord 
               !type variables defined in module 
               !operators_on_coordinates 
               xb=(xb-oldmid)*scale + rmid
               vectz=vectz*scale

               do ix=ismin%i,ismax%i
                  do iy=ismin%j,ismax%j
                     do iz=ismin%k,ismax%k
                        !vectz=A-B; modul=|A-B|;tmpvect1=P-B; 
                        !tmp1=(A-B)(P-B)
                        !mod2=(P-B)**2; modul2=(A-B)**2, 
                        !tmp=exrad*|A-B|.
                        !2011-05-13  Using operations on coord and 
                        !int_coord type variables defined in module 
                        !operators_on_coordinates 
                        ixyz=int_coord(ix,iy,iz)
                        tmpvect1=float(ixyz)-xb
                        tmp1=tmpvect1.dot.vectz
                        
                        if((tmp1.ge.-tmp).and.(tmp1.le.modul2+tmp))then
                           mod2=tmpvect1.dot.tmpvect1
                           if((mod2-(tmp1/modul)**2).le.radp2)&
                                             & idebmap(ix,iy,iz)=ipore
                        end if

                        !x-offset
                        tmpvect1%x=tmpvect1%x+.5
                        tmp1=tmpvect1.dot.vectz
                        if ((tmp1.ge.0.0).and.(tmp1.le.modul2)) then
                           mod2=tmpvect1.dot.tmpvect1
                           if ((mod2-(tmp1/modul)**2).le.rad2) then
                              iepsmp(ix,iy,iz)%i=natom+1+ii+&
                                                        &imedia*epsdim
                           end if
                        end if

                        !y-offset
                        tmpvect1%x=tmpvect1%x-.5
                        tmpvect1%y=tmpvect1%y+.5
                        tmp1=tmpvect1.dot.vectz
                        if ((tmp1.ge.0.0).and.(tmp1.le.modul2)) then
                           mod2=tmpvect1.dot.tmpvect1
                           if ((mod2-(tmp1/modul)**2).le.rad2) then
                              iepsmp(ix,iy,iz)%j=natom+1+ii+&
                                                        &imedia*epsdim
                           end if
                        end if

                        !z-offset
                        tmpvect1%y=tmpvect1%y-.5
                        tmpvect1%z=tmpvect1%z+.5
                        tmp1=tmpvect1.dot.vectz
                        if ((tmp1.ge.0.0).and.(tmp1.le.modul2)) then
                           mod2=tmpvect1.dot.tmpvect1
                           if ((mod2-(tmp1/modul)**2).le.rad2) then
                              iepsmp(ix,iy,iz)%k=natom+1+ii+&
                                                        &imedia*epsdim
                           end if
                        end if
                     end do
                  end do
               end do
            case(3) !if (objecttype.eq.3) then
                    !dealing with a cone     
               read(strtmp(16:80),*)kind,xa,xb,alpha 
               ipore=(kind.eq.2)
 
               !conversion degrees --> radiants
               alpha=alpha*3.14159265359/180.
               tan2=tan(alpha*.5)**2
               read(strtmp1,'(5f8.3)')vectz,modul,modul2
 
               !here we work in grid units
               modul=modul*scale
               modul2=modul*modul

               !2011-05-13  Using operations on coord and int_coord 
               !type variables defined in module 
               !operators_on_coordinates 
               xb=(xb-oldmid)*scale + rmid
               vectz=vectz*scale

               do ix=ismin%i,ismax%i
                  do iy=ismin%j,ismax%j
                     do iz=ismin%k,ismax%k
                        ixyz=int_coord(ix,iy,iz)
                        tmpvect1=(float(ixyz)-rmid)/scale+oldmid

                        !2011-05-13 Parameters transfered via module 
                        !architecture
                        call distobj(tmpvect1,dist,vecdist,ii,exrad,&
                                    &.true.)

                        if (dist.le.0.0) idebmap(ix,iy,iz)=ipore
                             
                        !vectz=A-B; tmpvect1=P-B; tmp1=(A-B)(P-B)
                        !mod2=(P-B)**2; modul2=(A-B)**2.

                        !x-offset
                        !2011-05-13  Using operations on coord and 
                        !int_coord type variables defined in module 
                        !operators_on_coordinates 
                        tmpvect1=float(ixyz)-xb
                        tmpvect1%x=tmpvect1%x+0.5
                        tmp1=tmpvect1.dot.vectz

                        if ((tmp1.ge.0.0).and.(tmp1.le.modul2)) then
                           mod2=(tmpvect1.dot.tmpvect1)*modul2
                           tmp2=(1+tan2)*tmp1*tmp1
                           if (mod2.le.tmp2) then
                              iepsmp(ix,iy,iz)%i=natom+1+ii+&
                                                        &imedia*epsdim
                           end if
                        end if

                        !y-offset
                        tmpvect1%x=tmpvect1%x-0.5
                        tmpvect1%y=tmpvect1%y+0.5
                        tmp1=tmpvect1.dot.vectz

                        if ((tmp1.ge.0.0).and.(tmp1.le.modul2)) then
                           mod2=(tmpvect1.dot.tmpvect1)*modul2
                           tmp2=(1+tan2)*tmp1*tmp1
                           if (mod2.le.tmp2) then
                              iepsmp(ix,iy,iz)%j=natom+1+ii+&
                                                        &imedia*epsdim
                           end if
                        end if
 
                        !z-offset
                        tmpvect1%y=tmpvect1%y-0.5
                        tmpvect1%z=tmpvect1%z+0.5
                        tmp1=tmpvect1.dot.vectz
                        if ((tmp1.ge.0.0).and.(tmp1.le.modul2)) then
                           mod2=(tmpvect1.dot.tmpvect1)*modul2
                           tmp2=(1+tan2)*tmp1*tmp1
                           if (mod2.le.tmp2) then
                              iepsmp(ix,iy,iz)%k=natom+1+ii+&
                                                        &imedia*epsdim
                           end if
                        end if
                           
                     end do
                  end do
               end do
            case(4) !if (objecttype.eq.4) then
                    !dealing with a parallelepiped
               read(strtmp(16:80),*)kind,xa,xb,xc,xd 
               ipore=(kind.eq.2)
               read(strtmp1,'(12f8.3)')vectz,modul,vecty,mody,vectx,modx
 
               !conversion to axial symmetry points
               !2011-05-13  Using operations on coord and int_coord 
               !type variables defined in module 
               !operators_on_coordinates 
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
                        !vectz=A-B;vectx=C-D;vecty=(C+D)/2-B;
                        !tmp1=|C-D|/2   
                        !modul2=(A-B)**2, tmp2=mody/2
                        !dot=(P-B)(..-..); 
                        ixyz=int_coord(ix,iy,iz)
                        xp=(float(ixyz)-rmid)/scale + oldmid

                        !2011-05-13 Parameters transfered via module 
                        !architecture
                        call distobj(xp,dist,vecdist,ii,exrad,.true.)

                        if (dist.le.0.0) idebmap(ix,iy,iz)=ipore
            
                        !now xp=P-B;
                        !x-offset
                        xp=float(ixyz)-xb
                        xp%x=xp%x+0.5
                        dot=vectz.dot.xp

                        if ((dot.ge.0.0).and.(dot.le.modul2)) then
                           dot=vectx.dot.xp
                           if (abs(dot).le.tmp1) then
                              dot=vecty.dot.xp
                              if (abs(dot).le.tmp2) then
                                 iepsmp(ix,iy,iz)%i=natom+1+ii+&
                                                        &imedia*epsdim
                              end if
                           end if
                        end if

                        !y-offset
                        xp%x=xp%x-0.5
                        xp%y=xp%y+0.5
                        dot=vectz.dot.xp

                        if ((dot.ge.0.0).and.(dot.le.modul2)) then
                           dot=vectx.dot.xp
                           if (abs(dot).le.tmp1) then
                              dot=vecty.dot.xp
                              if (abs(dot).le.tmp2) then
                                 iepsmp(ix,iy,iz)%j=natom+1+ii+&
                                                        &imedia*epsdim
                              end if
                           end if
                        end if

                        !z-offset
                        xp%y=xp%y-0.5
                        xp%z=xp%z+0.5
                        dot=vectz.dot.xp

                        if ((dot.ge.0.0).and.(dot.le.modul2)) then
                           dot=vectx.dot.xp
                           if (abs(dot).le.tmp1) then
                              dot=vecty.dot.xp
                              if (abs(dot).le.tmp2) then
                                 iepsmp(ix,iy,iz)%k=natom+1+ii+&
                                                        &imedia*epsdim
                              end if
                           end if
                        end if
                            
                     end do
                  end do
               end do
            end select
         end if
      end do !end of setting in OBJECTS

      !set interiors in MOLECULES
      if(itest2.or.itobig) write(6,*)'setout method 1',itest2,itobig
     
      DoATOMS: do iv=1,natom
         !restore values
         rad= delphipdb(iv)%rad3

         !Using operations on coord and int_coord type variables defined
         !in module operators_on_coordinates 
         xn=xn2(iv)

         !2011-05-13 Removed GOTO statement
         if (rad.lt.1.e-6) then
            cycle DoATOMS
         end if

         !scale radius to grid
         rad=rad*scale; rad5=(rad+0.5)**2; radp=rad+exrad*scale
         rad=rad+radprobe*scale
         rad4=(rad+0.5)**2; rad2=rad*rad; radp2=radp*radp

         !set dielectric map
         !check if sphere sits within limits of box
         itest2=.false.

         !2011-05-13 Using operations on coord and int_coord type 
         !variables defined in module operators_on_coordinates 
         ismin=int(xn-radmax2-1.)
         ismax=int(xn+radmax2+1.)
         itest=ismin
         ismin=min(ismin,igrid) ; ismin=max(ismin,1)
         if(itest.vorne.ismin) itest2=.true.
         itest=ismax
         ismax=min(ismax,igrid) ; ismax=max(ismax,1)
         if(itest.vorne.ismax) itest2=.true.

         if (itest2.or.itobig) then !slow method
            !2011-05-13 Seems redundant statement
            rad2a = rad2 - 0.25

            do iz=ismin%k,ismax%k
               do iy=ismin%j,ismax%j
                  do ix=ismin%i,ismax%i
                     !2011-05-13  Using operations on coord and 
                     !int_coord type variables defined in module 
                     !operators_on_coordinates 
                     ixyz=int_coord(ix,iy,iz); dxyz=float(ixyz)-xn
                     distsq=dxyz.dot.dxyz;     dxyz=dxyz+distsq

                     if (dxyz%x.lt.rad2a) then	            
                        iepsmp(ix,iy,iz)%i=iv+1+iatmmed(iv)*epsdim
                     end if

                     if (dxyz%y.lt.rad2a) then
                        iepsmp(ix,iy,iz)%j=iv+1+iatmmed(iv)*epsdim
                     end if

                     if (dxyz%z.lt.rad2a) then
                        iepsmp(ix,iy,iz)%k=iv+1+iatmmed(iv)*epsdim
                     end if

                     if(distsq.lt.radp2) idebmap(ix,iy,iz) =.false.
                  end do
               end do
            end do
         else !faster method
            !IT HAS PROBLEMS!!!! Walter (be careful before using 
            !also in multidilectric case!!!.and..not.isitmd
            rad2a=rad2-0.25

            !2011-05-13  Using operations on coord and int_coord type 
            !variables defined in module operators_on_coordinates 
            ixn=nint(xn); fxn=float(ixn)-xn; rad2av=rad2a-fxn

            do ix=-lim,lim
               vtemp=real(ix)+fxn
               sq(ix)=vtemp*vtemp
               rad2aav(ix)=rad2a-vtemp
            end do

            !adjust inter-atom, different epsilon bgps+++04/2004 Walter
            if (nmedia.gt.1.and.ionlymol) then
               do i=1,ibox
                  !2011-05-14  Using operations on coord and int_coord 
                  !type variables defined in module 
                  !operators_on_coordinates 
                  i123=ioff(i)
                  ixyz=ixn+i123
                  ix=ixyz%i; iy=ixyz%j; iz=ixyz%k
                  distsq=sq(i123%i)%x+sq(i123%j)%y+sq(i123%k)%z

                  if (distsq.lt.rad2aav(i123%i)%x)  then
                     iac=mod(iepsmp(ix,iy,iz)%i,epsdim)-1	      

                     if (iac.eq.-1.or.iac.gt.natom) then
                        !occhio! non so cosa fa con i pori!!
                        iepsmp(ix,iy,iz)%i=iv+1+iatmmed(iv)*epsdim
                     else
                        !2011-05-14  Using operations on coord and 
                        !int_coord type variables defined in module 
                        !operators_on_coordinates 
                        ddxyz=float(ixyz)-xn; ddxyz%x= ddxyz%x+0.5
                        dis2min1=ddxyz.dot.ddxyz-rad2  

                        ddxyz=float(ixyz)-xn2(iac); ddxyz%x= ddxyz%x+0.5
                        dis2min2=ddxyz.dot.ddxyz-&
                                       &(delphipdb(iac)%rad3*scale)**2 

                        if (dis2min2.gt.dis2min1) iac=iv
                        iepsmp(ix,iy,iz)%i=iac+1+iatmmed(iac)*epsdim
                     end if
                  end if

                  if (distsq.lt.rad2aav(i123%j)%y) then
                     iac=mod(iepsmp(ix,iy,iz)%j,epsdim)-1	      
                     if (iac.eq.-1.or.iac.gt.natom) then
                        !occhio! non so cosa fa con i pori!!
                        iepsmp(ix,iy,iz)%j=iv+1+iatmmed(iv)*epsdim
                     else
                        !2011-05-14  Using operations on coord and 
                        !int_coord type variables defined in module 
                        !operators_on_coordinates 
                        ddxyz=float(ixyz)-xn; ddxyz%y= ddxyz%y+0.5
                        dis2min1=ddxyz.dot.ddxyz-rad2
  
                        ddxyz=float(ixyz)-xn2(iac); ddxyz%y= ddxyz%y+0.5
                        dis2min2=ddxyz.dot.ddxyz&
                                      &-(delphipdb(iac)%rad3*scale)**2 

                        if (dis2min2.gt.dis2min1) iac=iv
                        iepsmp(ix,iy,iz)%j=iac+1+iatmmed(iac)*epsdim
                     end if
                  end if

                  if (distsq.lt.rad2aav(i123%k)%z) then
                     iac=mod(iepsmp(ix,iy,iz)%k,epsdim)-1	      
                     if (iac.eq.-1.or.iac.gt.natom) then
                        !occhio! non so cosa fa con i pori!!
                        iepsmp(ix,iy,iz)%k=iv+1+iatmmed(iv)*epsdim
                     else
                        !2011-05-14  Using operations on coord and 
                        !int_coord type variables defined in module 
                        !operators_on_coordinates 
                        ddxyz=float(ixyz)-xn; ddxyz%z= ddxyz%z+0.5
                        dis2min1=ddxyz.dot.ddxyz-rad2 
 
                        ddxyz=float(ixyz)-xn2(iac); ddxyz%z=ddxyz%z+0.5
                        dis2min2=ddxyz.dot.ddxyz&
                                      &-(delphipdb(iac)%rad3*scale)**2 

                        if (dis2min2.gt.dis2min1) iac=iv
                        iepsmp(ix,iy,iz)%k=iac+1+iatmmed(iac)*epsdim
                     end if
                  end if

                  if(distsq.lt.radp2) idebmap(ix,iy,iz)=.false.
               end do
            else 
               do i=1,ibox
                  !2011-05-14  Using operations on coord and int_coord 
                  !type variables defined in module 
                  !operators_on_coordinates 
                  i123=ioff(i); ixyz=ixn+i123
                  ix=ixyz%i; iy=ixyz%j; iz=ixyz%k
                  distsq = sq(i123%i)%x +sq(i123%j)%y + sq(i123%k)%z

                  if (distsq.lt.rad2aav(i123%i)%x)  then
                     iepsmp(ix,iy,iz)%i=iv+1+iatmmed(iv)*epsdim
                  end if

                  if (distsq.lt.rad2aav(i123%j)%y) then
                     iepsmp(ix,iy,iz)%j=iv+1+iatmmed(iv)*epsdim
                  end if

                  if (distsq.lt.rad2aav(i123%k)%z) then
                     iepsmp(ix,iy,iz)%k=iv+1+iatmmed(iv)*epsdim
                  end if

                  if (distsq.lt.radp2) idebmap(ix,iy,iz)=.false.
               end do
            end if
         end if 
      end do DoATOMS !end do of atoms

      !Array deallocation is made by ordinary F95 statement
      if(allocated(ioff)) deallocate(ioff)

      if (debug) then
         open(52,file='iepsmapnewfirst',form='formatted')
         write (52,*)igrid,igrid*igrid
         do ix=1,igrid
            do iy=1,igrid
               iz=(igrid+1)/2
               write (52,*)ix,iy,iz,iepsmp(ix,iy,iz)%k
            end do
         end do
         close (52)
      end if

      if(verbose)write(6,*)'Ending creating Van der Waals  Epsilon Map '

      return

      end subroutine setout

