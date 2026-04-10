!#####################################################################
!ATTENTION! This file is part of crgarrmod module.
!           Do not compile separately!
!2011-06-01 All other parameters transfered via qlog and pointers 
!           modules
!           Additional logical parameter ass needed when subroutine is 
!           run only to determine size of allocatable arrays where 
!            assigned charges are stored 
!#####################################################################

      subroutine distrTOpoint(ic1,ass)

      !include "pointer.h"

      character(80) :: strtmp
      
      !linkobj=number of object which is thought to contain the charge
      !        integer crgatn(1),natom,linkobj,crgatnval

      !2011-06-01 Changed to coord type variables
      type(coord) :: xb,xa,xu,xv,xw,tmpvect,xc,xd,xp,xm,tmpv

      !ccharge = total charge in distribution, as a fraction of e

      !2011-06-01 pi is declared as a parameter in qlog module 
      !           pi=3.1416
      logical :: ass

      !2011-05-30 Declarations added due  to IMPLICIT NONE    
      real :: temp,tmp,tmp1,tmp2,modul,h,aw,alpha,crgassigned
      real :: rmid,deltau,deltaw,dtetavol,drhovol,du,dttvol
      real :: delta,fractcharge,fractcharge1,fract,nrof,radius
      real :: ro,rho,teta,tt,tu,tv,ttmax,vcurr,ucurr,dv,deltav
      real :: ccharge,fnt
      integer :: crgatnval,ic1,linkobj,ncell,i,nro,nfi,nteta,ntetamax
      integer :: nt,ntt,ntv,itv,ntu,itu

      h=1./scale; rmid=real((igrid+1)/2)
      if(ass)write(6,*) "number of charge distrib. present ",ndistr 

      do i=1,ndistr
         strtmp=datadistr(i)
         read(strtmp(14:16),*)linkobj
         crgatnval=-i
         if (linkobj.gt.0) crgatnval=natom+linkobj
         read(strtmp(18:23),*)ccharge

         !2011-06-01 Multiple IFs replaced by SELECT CASE construct
         select case(strtmp(8:10))
         case('  1') !if (strtmp(8:10).eq.'  1') then
                     !here we have a spherical distribution
            read(strtmp(27:80),*)xb,radius
               
            if (ass) then
               qnet=qnet+ccharge
                  
               if (ccharge.ne.0.0) then
                  if (ccharge.gt.0.0) then
                     qplus=qplus+ccharge
               
                     !2011-06-01 Using operations on coord type 
                     !variables defined
                     !in module operators_on_coordinates 
                     cqplus=cqplus+(ccharge*xb)
                  else
                     qmin=qmin+ccharge
                     cqmin=cqmin+(ccharge*xb)
                  end if
               end if
            end if

            !if (strtmp(11:12).eq.' v') then
            !   uniform in the volume
            !   tmp=scale*scale*scale*radius*radius*radius
            !   fractcharge=3*charge/(4*pi*tmp)
            !   do nro=1,radius*scale
            !      do nteta=1,pi*nro
            !         tmp=nro*sin(nteta/nro)
            !         do nfi=1,2*pi*tmp
            !            ic1=ic1+1
            !            chgpos(1,ic1)=xb(1)+tmp*cos(nfi/tmp)/scale 
            !            chgpos(2,ic1)=xb(2)+tmp*sin(nfi/tmp)/scale
            !            chgpos(3,ic1)=xb(3)+nro*cos(nteta/nro)/scale
            !            atmcrg(1,ic1)=&
            !                   &(chgpos(1,ic1)-cmid(1))*scale+rmid
            !            atmcrg(2,ic1)=&
            !                   &(chgpos(2,ic1)-cmid(2))*scale+rmid
            !            atmcrg(3,ic1)=&
            !                   &(chgpos(3,ic1)-cmid(3))*scale+rmid
            !            atmcrg(4,ic1)=fractcharge
            !            crgatn(ic1)=crgatnval
            !         end do
            !      end do
            !   end do
            !end if
               
            !the new version,improved, needs less charges to fill 
            !the object
            !2011-06-01 Multiple IFs replaced by SELECT CASE construct 
            select case(strtmp(11:12))
            case(' v') !if (strtmp(11:12).eq.' v') then
               !uniform in the volume
               crgassigned=0.0
               tmp=scale*scale*scale*radius*radius*radius
               fractcharge1=3*ccharge/(4*pi*tmp)
               tmp=radius*scale
        
               do nro=0,int(sqrt(2*tmp)-1.0)
                  ncell=ic1
                        
                  !tmp1=current radius*scale=1/deltateta
                  tmp1=(tmp-.5*(nro+1)*(nro+1))
                      
                  fractcharge=fractcharge1*(nro+1)*&
                                         &  2*sin(1./(2*tmp1))*(tmp1+&
                                         &  (nro+1)*(nro+1)/(12*tmp1))
               
                  write(6,*)"distrib: ro:",tmp1/scale
                      
                  do nteta=0,int(pi*tmp1-.5)
                     !tmp2=1/deltaphi
                     tmp2=tmp1*sin((.5+nteta)/tmp1)
                         
                     do nfi=0,int(2*pi*tmp2-.5)
                        ic1=ic1+1
                              
                        !2011-06-01 Changed to coord and 
                        !grid_value variable types
                        if (ass) then
                           chgpos(ic1)%x=xb%x+&
                                      &(tmp2*cos((nfi+.5)/tmp2)/scale)
                           chgpos(ic1)%y=xb%y+&
                                      &(tmp2*sin((nfi+.5)/tmp2)/scale)
                           chgpos(ic1)%z=xb%z+&
                                      &(tmp1*cos((nteta+.5)/tmp1)/scale)

                           atmcrg(ic1)%xyz=((chgpos(ic1)-&
                                                    &cmid)*scale+rmid)
                           atmcrg(ic1)%value=crgatnval

                           crgatn(ic1)=crgatnval
                        end if
                     end do
                  end do
                         
                  crgassigned=crgassigned+fractcharge*(ic1-ncell)
               end do
                   
               if (abs(ccharge-crgassigned).gt.1.e-6) then
                  ic1=ic1+1
                         
                  if (ass) then
                     chgpos(ic1)=xb; atmcrg(ic1)%value=crgatnval
                     atmcrg(ic1)%xyz=((chgpos(ic1)-cmid)*scale+rmid)
                     atmcrg(ic1)%value=ccharge-crgassigned
                     crgatn(ic1)=crgatnval
                     if(fractcharge*ccharge.lt.0.)&
                                     &write(6,*)"WARNING fc*charge<0!"
                     write(6,*)"charge not uniformly assigned:",&
                                                  &ccharge-crgassigned
                  end if
               end if
            case(' s')  !if (strtmp(11:12).eq.' s') then
               !uniform on the surface
               tmp=scale*scale*radius*radius
               fractcharge=ccharge/(4*pi*tmp)
               ro=radius*scale
                   
               do nteta=1,int(pi*ro)
                  tmp=ro*sin(nteta/ro)
                  do nfi=1,int(2*pi*tmp)
                     ic1=ic1+1
                     if (ass) then
                        chgpos(ic1)%x=xb%x+(tmp*cos(nfi/tmp)/scale)
                        chgpos(ic1)%y=xb%y+(tmp*sin(nfi/tmp)/scale)
                        chgpos(ic1)%z=xb%z+(radius*cos(nteta/ro))
                        atmcrg(ic1)%xyz=(chgpos(ic1)-cmid)*scale+rmid
                        atmcrg(ic1)%value=fractcharge
                        crgatn(ic1)=crgatnval
                     end if
                  end do
               end do
            end select
         case('  2') !if (strtmp(8:10).eq.'  2') then
            !here we have a cylindric distribution
            read(strtmp(27:80),*)xa,xb,radius
            call basisortho(xa,xb,xu,xv,xw,modul)
         
            if (ass) then
               qnet=qnet+ccharge
      
               if (ccharge.ne.0.0) then
                  if (ccharge.gt.0.0) then
                     qplus=qplus+ccharge
                     cqplus=cqplus+((0.5*ccharge)*(xa+xb))
                  else
                     qmin=qmin+ccharge
                     cqmin=cqmin+((0.5*ccharge)*(xa+xb))
                  end if
               end if
            end if
         
            select case(strtmp(11:12))
            case(' v') !if (strtmp(11:12).eq.' v') then
                
               !uniform in the volume
               !tmp=scale*scale*scale*radius*radius*modul
               !fractcharge=charge/(pi*tmp)
               !tmp1=radius*scale-.5
               !deltanro=tmp1/int(tmp1)
               !do nro=tmp1,0,-deltanro
               !   do nteta=1,2*pi*nro
               !      tmp=nro*sin(float(nteta)/nro)
               !      tmp1=nro*cos(float(nteta)/nro)
               !      do nt=.5,scale*modul-.5
               !         fractcharge=fract/
               !         ic1=ic1+1
               !         chgpos(1,ic1)=xb(1)+(tmp1*xu(1)&
               !                     &+tmp*xv(1)+nt*xw(1))*h
               !         chgpos(2,ic1)=xb(2)+(tmp1*xu(2)&
               !                     &+tmp*xv(2)+nt*xw(2))*h
               !         chgpos(3,ic1)=xb(3)+(tmp*xv(3)+nt*xw(3))*h
               !         atmcrg(1,ic1)=(chgpos(1,ic1)-&
               !                                  &cmid(1))*scale+rmid
               !         atmcrg(2,ic1)=(chgpos(2,ic1)-&
               !                                  &cmid(2))*scale+rmid
               !         atmcrg(3,ic1)=(chgpos(3,ic1)-&
               !                                  &cmid(3))*scale+rmid
               !         atmcrg(4,ic1)=fractcharge
               !         crgatn(ic1)=crgatnval
               !      end do
               !   end do
               !end do
               !write(6,*)nro,nteta,nt

               !the new version,improved, needs less charges to fill 
               !the object
               xp=0.5*(xa+xb)

               tmp1=radius*radius*modul
               fract=ccharge/(pi*tmp1)
               rho=radius-.5*h; drhovol=h
               ttmax=.5*(modul-h); tt=ttmax; dttvol=h

               !2011-06-01 GOTO replaced by DO-IF-EXIT construct
               D100: do
                  ntetamax=2*pi*rho*scale
        
                  if (ntetamax.ne.0) then
                     do nteta=1,ntetamax
                        dtetavol=2.*pi/ntetamax
                        teta=nteta*dtetavol
                        tmp1=rho*cos(teta)
                        tmp2=rho*sin(teta)
                        fractcharge=fract*rho*drhovol*dtetavol*dttvol
                        if (tt.eq.0.) fractcharge=2.*fractcharge 

                        ic1=ic1+1

                        !2011-06-01 Changed to coord and 
                        !grid_value variable types 
                        !(xu%z =0)     
                        if (ass) then
                           chgpos(ic1)=xp+(tmp1*xu+tmp2*xv)+tt*xw
                           atmcrg(ic1)%xyz=(chgpos(ic1)-&
                                                     &cmid)*scale+rmid
                           atmcrg(ic1)%value=fractcharge
                           crgatn(ic1)=crgatnval
                        end if
      
                        if (tt.gt.0.) then
                           ic1=ic1+1
                  
                           if (ass) then
                              chgpos(ic1)=chgpos(ic1-1)-2.*tt*xw
                              atmcrg(ic1)%xyz=(chgpos(ic1)-&
                                                     &cmid)*scale+rmid
                              atmcrg(ic1)%value=fractcharge
                              crgatn(ic1)=crgatnval
                           end if
                        end if               
                     end do
                  else
                     if (tt.gt.0.) then
                        ic1=ic1+1
                                  
                        if (ass) then
                           chgpos(ic1)=xp-tt*xw
                           atmcrg(ic1)%xyz=(chgpos(ic1)-&
                                                     &cmid)*scale+rmid
                           atmcrg(ic1)%value=&
                                      &pi*fract*drhovol*drhovol*dttvol
                           crgatn(ic1)=crgatnval
                        end if
                     else
                        dttvol=2.*dttvol
                     end if
                      
                     ic1=ic1+1
                         
                     if (ass) then
                        chgpos(ic1)=xp+tt*xw
                        atmcrg(ic1)%xyz=(chgpos(ic1)-cmid)*scale+rmid
                        atmcrg(ic1)%value=&
                                      &pi*fract*drhovol*drhovol*dttvol
                        crgatn(ic1)=crgatnval
                     end if
                  end if
                     
                  tt=tt-.5*dttvol
                  dttvol=min(drhovol,dttvol+.5*h)
                  tt=tt-.5*dttvol
            
                  if (tt.gt.-.5*dttvol) then
                     if (tt.lt..5*dttvol) then
                        dttvol=tt+.5*dttvol
                        tt=0.0
                     end if
              
                     cycle D100
                  end if

                  tt=ttmax; dttvol=h; rho=rho-.5*drhovol
                  drhovol=drhovol+.5*h; rho=rho-.5*drhovol
            
                  if (rho.le.-.5*drhovol) then
                     exit D100
                  else
                     if (rho.lt..5*drhovol) then
                        drhovol=rho+.5*drhovol
                        rho=0.0
                     end if
                  end if
               end do D100
            case(' s') !if (strtmp(11:12).eq.' s') then
               !uniform on the LATERAL surface
               tmp=scale*scale*radius*modul 
               fractcharge=ccharge/(2*pi*tmp)
               ro=radius*scale

               do nteta=1,int(2*pi*ro)
                  tmp=ro*sin(nteta/ro)
                  tmp1=ro*cos(nteta/ro)
                     
                  do nt=1,int(modul*scale)
                     ic1=ic1+1
 
                     !2011-06-01  xu%z (xu(3) before) is 
                     !always zero (from basisortho subroutine)
                     if (ass) then
                        chgpos(ic1)=xb+((tmp1*xu)+(tmp*xv)&
                                                &+(real(nt)*xw))/scale
                        atmcrg(ic1)%xyz=(chgpos(ic1)-cmid)*scale+rmid
                        atmcrg(ic1)%value=fractcharge
                        crgatn(ic1)=crgatnval
                     end if
                  end do
               end do
            end select
         case('  3') !if (strtmp(8:10).eq.'  3') then
            !here we have a conic distribution
            read(strtmp(27:80),*)xa,xb,alpha 
            alpha=alpha*pi/180
                 
            call basisortho(xa,xb,xu,xv,xw,modul)
              
            radius=tan(alpha)*modul

            if (ass) then
               qnet=qnet+ccharge
               write(6,*)"center of charge for cone not yet done!Only net charge"
   
               if (.false..and.ccharge.ne.0.0) then
                  if (ccharge.gt.0.0) then
                     !!c è da verificare assegnazione qcharge!!!
                     qplus=qplus+ccharge
                     cqplus=cqplus+(ccharge/3.)*((2.*xa)+xb)
                  else
                     qmin=qmin+ccharge
                     cqmin=cqmin+(ccharge/3.)*((2.*xa)+xb)
                  end if
               end if
            end if

            select case(strtmp(11:12))
            case(' v') !if (strtmp(11:12).eq.' v') then
               !uniform in the volume
               tmp=(scale**3)*radius*radius*modul
               fractcharge=3*ccharge/(pi*tmp)
         
               do nro=1,int(radius*scale)
                  do nteta=1,int(2*pi*nro)
                     tmp=nro*sin(real(nteta)/nro)
                     tmp1=nro*cos(real(nteta)/nro)
                     do nt=1,int(modul*(scale-nro/radius))
                        ic1=ic1+1
                           
                        if (ass) then
                           chgpos(ic1)=xb+((tmp1*xu)+(tmp*xv)&
                                                &+(real(nt)*xw))/scale
                           atmcrg(ic1)%xyz=(chgpos(ic1)-&
                                                     &cmid)*scale+rmid
                           atmcrg(ic1)%value=fractcharge
                           crgatn(ic1)=crgatnval
                        end if
                     end do
                  end do
               end do
            case(' s') !if (strtmp(11:12).eq.' s') then
               !uniform on the LATERAL surface
               tmp=scale*scale*radius*modul
               fractcharge=ccharge*cos(alpha)/(pi*tmp)
         
               do nro=1,int(radius*scale/sin(alpha))
                  nrof=nro*sin(alpha)
                  do nteta=1,int(2*pi*nrof)
                     tmp=nrof*sin(real(nteta)/nrof)
                     tmp1=nrof*cos(real(nteta)/nrof)
                     nt=modul*(scale-nrof/radius)
                     ic1=ic1+1

                     if (ass) then
                        chgpos(ic1)=xb+((tmp1*xu)+(tmp*xv)+&
                                                 &(real(nt)*xw))/scale
                        atmcrg(ic1)%xyz=(chgpos(ic1)-cmid)*scale+rmid
                        atmcrg(ic1)%value=fractcharge
                        crgatn(ic1)=crgatnval
                     end if
                  end do
               end do
            end select
         case('  4') !if (strtmp(8:10).eq.'  4') then
            !here we have a parallelepiped distribution
            read(strtmp(27:80),*)xa,xb,xc,xd 

            !conversion to axial symmetry points
            !nuova terna
            xu=xb-xa; tmp1=sqrt(xu.dot.xu)
            xv=xc-xa; tmp2=sqrt(xv.dot.xv)
            xw=xd-xa; tmp=sqrt(xw.dot.xw)

            !ora trovo il centro xm 
            xm=(xd+xc+xb-xa)*0.5

            !assegno momenti di carica
            if (ass) then
               qnet=qnet+ccharge
      
               if (ccharge.ne.0.0) then
                  if (ccharge.gt.0.0) then
                     qplus=qplus+ccharge
                     cqplus=cqplus+(ccharge*xm)
                  else
                     qmin=qmin+ccharge
                     cqmin=cqmin+ccharge*xm
                  end if
               end if
            end if

            !ordering by lenght |xw| < |xu| < |xv|
            if (tmp2.lt.tmp1) call swap(xv,xu,tmp2,tmp1)
            if (tmp1.lt.tmp)  call swap(xu,xw,tmp1,tmp)
            if (tmp2.lt.tmp1) call swap(xv,xu,tmp2,tmp1)   

            !from now on  xu,xv are versors, xw is a vector
            xu=xu/tmp1; xv=xv/tmp2

            crgassigned=0.0
            fractcharge=ccharge/(tmp*tmp1*tmp2)
            ntt=int(sqrt(scale*tmp))
         
            do nt=1,ntt
               fnt=(1./sqrt(scale*tmp))+(nt-1)*(1./sqrt(scale*tmp))
               deltaw=fnt*sqrt(tmp/scale)                   
               deltau=h; ucurr=(tmp1+h)/2; du=deltau
      
               D1000: do
                  deltau=min(deltau,deltaw)
                  du=(du+deltau)/2; ucurr=ucurr-du
                  if (ucurr.lt.0.0) exit D1000
             
                  deltav=h; vcurr=(tmp2+h)/2; dv=deltav

                  D1010: do
                     deltav=min(deltav,deltau)
                     dv=(dv+deltav)/2; vcurr=vcurr-dv
                     if (vcurr.lt.0.0) exit D1010  
                     aw=.5*(1-fnt*fnt)
                     fractcharge1=fractcharge*deltau*deltav*deltaw

                     ic1=ic1+1
                     crgassigned=crgassigned+fractcharge1
      
                     if (ass) then
                        chgpos(ic1)=xm+(ucurr*xu)+(vcurr*xv)+(aw*xw)
                        atmcrg(ic1)%xyz=(chgpos(ic1)-cmid)*scale+rmid
                        atmcrg(ic1)%value=fractcharge1
                        crgatn(ic1)=crgatnval
                     end if

                     if (aw.gt.1.e-6) then
                        ic1=ic1+1
                        crgassigned=crgassigned+fractcharge1
                  
                        if (ass) then
                           chgpos(ic1)=xm+ucurr*xu+vcurr*xv-(aw*xw)
                           atmcrg(ic1)%xyz=&
                                        &(chgpos(ic1)-cmid)*scale+rmid
                           atmcrg(ic1)%value=fractcharge1
                           crgatn(ic1)=crgatnval
                        end if
                     end if
            
                     if (ucurr.gt.1.e-6) then
                        ic1=ic1+1
                        crgassigned=crgassigned+fractcharge1
                  
                        if (ass) then
                           chgpos(ic1)=&
                                     &xm-(ucurr*xu)+(vcurr*xv)+(aw*xw)
                           atmcrg(ic1)%xyz=&
                                        &(chgpos(ic1)-cmid)*scale+rmid
                           atmcrg(ic1)%value=fractcharge1
                           crgatn(ic1)=crgatnval
                        end if
                            
                        if (aw.gt.1.e-6) then
                           ic1=ic1+1
                           crgassigned=crgassigned+fractcharge1
                  
                           if (ass) then
                              chgpos(ic1)=xm-(ucurr*xu)&
                                                  &+(vcurr*xv)-(aw*xw)
                              atmcrg(ic1)%xyz=&
                                        &(chgpos(ic1)-cmid)*scale+rmid
                              atmcrg(ic1)%value=fractcharge1
                              crgatn(ic1)=crgatnval
                           end if
                        end if
                     end if

                     if (vcurr.gt.1.e-6) then    
                        ic1=ic1+1
                        crgassigned=crgassigned+fractcharge1
                        
                        if (ass) then
                           chgpos(ic1)=&
                                     &xm+(ucurr*xu)-(vcurr*xv)+(aw*xw)
                           atmcrg(ic1)%xyz=&
                                        &(chgpos(ic1)-cmid)*scale+rmid
                           atmcrg(ic1)%value=fractcharge1
                           crgatn(ic1)=crgatnval
                        end if
                            
                        if (aw.gt.1.e-6) then
                           ic1=ic1+1
                           crgassigned=crgassigned+fractcharge1
                         
                           if (ass) then
                              chgpos(ic1)=xm+(ucurr*xu)-&
                                                   &(vcurr*xv)-(aw*xw)
                              atmcrg(ic1)%xyz=&
                                        &(chgpos(ic1)-cmid)*scale+rmid
                              atmcrg(ic1)%value=fractcharge1
                              crgatn(ic1)=crgatnval
                           end if
                        end if
                            
                        if (ucurr.gt.1.e-6) then
                           ic1=ic1+1
                           crgassigned=crgassigned+fractcharge1
                  
                           if (ass) then
                              chgpos(ic1)=xm-(ucurr*xu)-&
                                                   &(vcurr*xv)+(aw*xw)
                              atmcrg(ic1)%xyz=&
                                        &(chgpos(ic1)-cmid)*scale+rmid
                              atmcrg(ic1)%value=fractcharge1
                              crgatn(ic1)=crgatnval
                           end if
                         
                           if (aw.gt.1.e-6) then
                              ic1=ic1+1
                              crgassigned=crgassigned+fractcharge1
                                     
                              if (ass) then
                                 chgpos(ic1)=xm-(ucurr*xu)-&
                                                   &(vcurr*xv)-(aw*xw)
                                 atmcrg(ic1)%xyz=&
                                        &(chgpos(ic1)-cmid)*scale+rmid
                                 atmcrg(ic1)%value=fractcharge1
                                 crgatn(ic1)=crgatnval
                              end if
                           end if
                        end if
                     end if
                 
                     dv=deltav
                     deltav=deltav+h
                  end do D1010

                  du=deltau
                  deltau=deltau+h
               end do D1000
            end do

            if (abs(ccharge-crgassigned).gt.1.e-6) then
               ic1=ic1+1
            
               if (ass) then
                  chgpos(ic1)=xm
                  atmcrg(ic1)%xyz=(chgpos(ic1)-cmid)*scale+rmid
                  atmcrg(ic1)%value=ccharge-crgassigned

                  if(fractcharge*ccharge.lt.0.)&
                                     &write(6,*)"WARNING fc*charge<0!"
                  write(6,*)"charge not uniformly assigned:",&
                                                  &ccharge-crgassigned
               end if
            end if
         case('  8') !if (strtmp(8:10).eq.'  8') then
            !here we have a discrete distribution of point charge
            read(strtmp(27:80),*)xb
            ic1=ic1+1
         
            if (ass) then
               qnet=qnet+ccharge
      
               if (ccharge.ne.0.0) then
                  if (ccharge.gt.0.0) then
                     qplus=qplus+ccharge
                     cqplus=cqplus+(ccharge*xb)
                  else
                     qmin=qmin+ccharge
                     cqmin=cqmin+ccharge*xb
                  end if
               end if

               chgpos(ic1)=xb               
               atmcrg(ic1)%xyz=(chgpos(ic1)-cmid)*scale+rmid
               atmcrg(ic1)%value=ccharge
               crgatn(ic1)=crgatnval
            end if
         case('  9') !if (strtmp(8:10).eq.'  9') then
                     !here we have a linear  distribution
            read(strtmp(27:80),*)xa,xb
      
            if (ass) then
               qnet=qnet+ccharge
      
               if (ccharge.ne.0.0) then
                  if (ccharge.gt.0.0) then
                     qplus=qplus+ccharge
                     cqplus=cqplus+(0.5*ccharge)*(xa+xb)
                  else
                     qmin=qmin+ccharge
                     cqmin=cqmin+(0.5*ccharge)*(xa+xb)
                  end if
               end if
            end if
            
            tmpvect=xa-xb; modul=tmpvect.dot.tmpvect
              
            tmp=1.0/(scale*modul)
            fractcharge=ccharge*tmp
            ntt=int(scale*modul)
         
            do nt=1,ntt
               fnt=tmp+(nt-1)*tmp
               ic1=ic1+1
   
               if (ass) then
                  chgpos(ic1)=xb+(fnt*tmpvect)               
                  atmcrg(ic1)%xyz=(chgpos(ic1)-cmid)*scale+rmid
                  atmcrg(ic1)%value=fractcharge
                  crgatn(ic1)=crgatnval
               end if
            end do
         case(' 10') !if (strtmp(8:10).eq.' 10') then
                     !here we have a circular surface distribution
            !c'è qualche discrepanza!! almeno con il manuale
            read(strtmp(27:80),*)xa,xb,radius
            call basisortho(xa,xb,xu,xv,xw,modul)
      
            if (ass) then
               qnet=qnet+ccharge

               if (ccharge.ne.0.0) then
                  if (ccharge.gt.0.0) then
                     qplus=qplus+ccharge
                     cqplus=cqplus+(ccharge*xa)
                  else
                     qmin=qmin+ccharge
                     cqmin=cqmin+(ccharge*xa)
                  end if
               end if
            end if

            tmp=scale*scale*radius*radius
            fractcharge=ccharge/(pi*tmp)
                  
            do nro=1,int(radius*scale)
               do nteta=1,int(2*pi*nro)
                  tmp=nro*sin(real(nteta)/nro)
                  tmp1=nro*cos(real(nteta)/nro)
                  ic1=ic1+1
       
                  if (ass) then
                     chgpos(ic1)=xb+(((tmp1*xu)+(tmp*xv))/scale)
                     atmcrg(ic1)%xyz=(chgpos(ic1)-cmid)*scale+rmid
                     atmcrg(ic1)%value=fractcharge
                     crgatn(ic1)=crgatnval
                  end if
               end do
            end do
         case(' 11') !if (strtmp(8:10).eq.' 11') then
                     !here we have a rectangular surface distribution
            read(strtmp(27:80),*)xb,xc,xd
      
            if (ass) then
               qnet=qnet+ccharge
   
               if (ccharge.ne.0.0) then
                  if (ccharge.gt.0.0) then
                     qplus=qplus+ccharge
                     cqplus=cqplus+(0.5*ccharge)*(xc+xd)
                  else
                     qmin=qmin+ccharge
                     cqmin=cqmin+(0.5*ccharge)*(xc+xd)
                  end if
               end if
            end if
              
            xu=xc-xd; xv=2.*xb-xd-xc
            tmp1=xu.dot.xu; tmp2=xv.dot.xv

            !tmp1=(C-D)**2;tmp2=(E-D)**2
            !xu=(C-D),xv=(E-D) are not necessarily unitary vectors
            !call inner(xu,xu,tmp1)
            !call inner(xv,xv,tmp2)
              
            fractcharge=ccharge/(tmp1*tmp2*(scale**2))

            ntv=int((scale*tmp2)/2)
          
            do itv=1,ntv
               tv=1./(scale*tmp2)+(itv-1)/(scale*tmp2)
               ntu=int((scale*tmp1)/2)
               do itu=1,ntu
                  tv=1./(scale*tmp1) + (itu-1)/(scale*tmp1)
                  ic1=ic1+1
                  if (ass) then
                     chgpos(ic1)=xb+(tu*xu)+(tv*xv)
                     atmcrg(ic1)%xyz=(chgpos(ic1)-cmid)*scale+rmid
                     atmcrg(ic1)%value=fractcharge
                     crgatn(ic1)=crgatnval
                  end if
       
                  ic1=ic1+1
          
                  if (ass) then
                     chgpos(ic1)=xb-(tu*xu)+(tv*xv)
                     atmcrg(ic1)%xyz=(chgpos(ic1)-cmid)*scale+rmid
                     atmcrg(ic1)%value=fractcharge
                     crgatn(ic1)=crgatnval
                  end if
         
                  ic1=ic1+1
            
                  if (ass) then
                     chgpos(ic1)=xb+(tu*xu)-(tv*xv)
                     atmcrg(ic1)%xyz=(chgpos(ic1)-cmid)*scale+rmid
                     atmcrg(ic1)%value=fractcharge
                     crgatn(ic1)=crgatnval
                  end if
      
                  ic1=ic1+1
            
                  if (ass) then
                     chgpos(ic1)=xb-(tu*xu)-(tv*xv)
                     atmcrg(ic1)%xyz=(chgpos(ic1)-cmid)*scale+rmid
                     atmcrg(ic1)%value=fractcharge
                     crgatn(ic1)=crgatnval
                  end if
               end do
            end do
         end select
      end do
 
      end subroutine distrTOpoint
       


