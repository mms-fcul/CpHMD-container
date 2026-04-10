      !  ATTENTION!  This file is part of wrtsitmod module.
      !==================================================================
      ! 2011-06-13 All other parameters are accessible via qlog and pointers
      subroutine rforce(afield,scrg)
      !	  subroutine rforce(afield,natom,scspos,scrg,
      type(coord) :: afield(natom+nobject),sfield,vtemp,xyz
      real :: scrg(icount2b)
      ! 2011-06-13 Arrays below are declared in pointers module
      !------------------------------------------------------------------
      ! 2011-06-13 Declarations added due to IMPLICIT NONE
      integer :: i,j,ipmax=0,cont=0,iat
      real :: sc,vvtemp,fact1,trullo,dist,sdist,temp
      !------------------------------------------------------------------
      ! 2011-06-13 Changed to array operations
      afield=coord(0.,0.,0)
      
      !  2011-06-12  Using operations on coord type variables defined
      do i=1,icount2b
         sfield=coord(0.,0.,0.)
         sc=scrg(i)
         
         !!c erased some comments (Walter)
         !!c calculate total field on this surface element due to ALL charges
         do j=1,nqass
            vtemp=scspos(i)-chgpos(j)
            dist=vtemp.dot.vtemp
            sdist=sqrt(dist)
            temp=atmcrg(j)%value/(dist*sdist)
            sfield=sfield+temp*vtemp
            !!c add negative of this to the charged atom
            iat=crgatn(j)
            !!c b+++++++++++++++++++++++++++++
            if (iat.lt.0) cycle
            !!c e++++++++++++++++++++++++++++++
            if(iat.eq.0)then
               write(6,*)'problems with crgatn '
               stop
            endif
            afield(iat)=afield(iat)-((sc*temp)*vtemp)
         end do
         sfield=sfield*scrg(i)
         j=atsurf(i)
         afield(j)=afield(j)+sfield
      end do
      end subroutine rforce
      
      !!c b+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rforceeps1(afield,scrg)
      
      type(coord) :: afield(natom+nobject),sfield,vtemp,xyz
      type(int_coord) :: ixyz
      real :: scrg(icount2b)
      ! 2011-06-13 Arrays below are declared in pointers module
      !------------------------------------------------------------------
      ! 2011-06-13 Declarations added due to IMPLICIT NONE
      integer :: p,i,j,iat,cont=0,ipmax=0
      real fact,fact1,deltas,sigmap,dist,sdist
      real :: sum=0., zax=0., realds=0., rmyds=0.0
      real :: trullo,temp,sc,vvtemp
      !------------------------------------------------------------------
      !!c       era un 4 bytes
      afield=coord(0.,0.,0.)
      fact=-2.*pi*80/(79.)
      do p=1,icount2b
         sfield=coord(0.,0.,0.)
         ixyz=nint((scspos(p)-oldmid)*scale)
         xyz=(float(ixyz)/scale) + oldmid
         sc=0.5*scrg(p)
         vtemp=scspos(p)-xyz; vvtemp=vtemp.dot.scsnor(p)
         fact1=0.8/(scale*scale)-(vvtemp*vvtemp)
         deltas=Pi*fact1; sigmap=scrg(p)/deltas
         realds=realds+scrg(p)/0.0393; rmyds=rmyds+deltas
         trullo=fact*sigmap*sigmap*deltas
         if(p.eq.400) then
            write(6,*)"normale",scsnor(p)
            write(*,*)"P",vtemp
            write(*,*)scspos(p)
            write(*,*)"area",scrg(p)/0.0393,scrg(p)*scale*scale/0.0393
         end if
         j=atsurf(p)
         !!c erased some comments (Walter)
         do iat=1,natom
            if(iat.eq.j) then
               cont=cont+1; afield(iat)%x=afield(iat)%x+trullo
            end if
         end do
         !!c calculate total field on this surface element due to ALL charges
         do i=1,-nqass
            vtemp=scspos(p)-chgpos(i); dist=vtemp.dot.vtemp
            
            sdist=sqrt(dist)
            temp=atmcrg(i)%value/(dist*sdist)
            
            sfield=sfield+temp*vtemp
            !!c add negative of this to the charged atom
            iat=crgatn(i)
            !!c b+++++++++++++++++++++++++++++
            if (iat.lt.0) exit
            !!c e++++++++++++++++++++++++++++++
            if(iat.eq.0)then
               write(6,*)'problems with crgatn '
               stop
            endif
            if(iat.eq.j) then
               afield(iat)=afield(iat)-sc*temp*vtemp+fact1*scsnor(p)
            end if
         end do
         sfield=sfield*sc
         
      end do
      write(*,*)"supcalc",realds/(16*3.14159265359),"mia",rmyds/(16*3.14159265359)
      end subroutine rforceeps1
	  
      !!c ***********************************************************
      subroutine rforcenew      
      !!c  atmcrg (1..3,nqass) atomic charge pos. in grid units, like xn2(1..3,natom)
      !!c     chgpos real charge pos in Amstrong
      !!c     scspos, surface  charge position in Amstrong
      real :: polariz(natom)
      !!c       erano un 4 bytes
      type(coord) :: fxyz,ffxyz,xyzi,xyzq,xyzp,distiq,distip
      !------------------------------------------------------------------
      ! 2011-06-13 Declarations added due ti IMPLICIT NONE
      integer :: i,imed,q,p
      real :: cost1, cost2,fact,fact1,fact2,fact3,chrgv4(natom)
      real :: qq,qi,qp
      real :: riq,riq2,rip,rip2,eps
      real distiqx,distiqy,distiqz,distipx,distipy,distipz
      !------------------------------------------------------------------
      if (.not.ionlymol) then
         write(6,*)"Not yet ready to give forcefield in case of objects"
         stop
      end if
      
      do i=1,natom
         fxyz=coord(0.,0.,0.)
         qi=delphipdb(i)%chrgv4; imed=iatmmed(i); eps=medeps(imed)*epkt
         cost1=3./(eps+2.); cost2=polariz(i)*epkt*cost1;xyzi=xn1(i)
         !!c ciclo sulle cariche vere
         do q=1,nqass
            if (crgatn(q).eq.i) cycle
            qq=atmcrg(q)%value/atmeps(q)
            xyzq=chgpos(q); distiq=xyzi-xyzq ; riq2=distiq.dot.distiq
            riq=sqrt(riq2);  fact=qq/(riq*riq2); fact1=-2.*fact*fact
            !!c     termine elettroforetico
            fxyz=fxyz+(qi*fact)*distiq
            !!c     termine dielettroforetico q=p
            ffxyz=fact1*distiq
            !!c     termine dielettroforetico q<>p col simmetrico
            do  p=1,q-1
               if (crgatn(p).eq.i) cycle
               qp=atmcrg(p)%value/atmeps(p)
               xyzp=chgpos(p); distip=xyzi-xyzp; rip2=distip.dot.distip
               rip=sqrt(rip2)
               fact2=qp*fact/(rip*rip2)
               fact3=(distiq.dot.distip)/(riq2*rip2)
               ffxyz=ffxyz+fact2*((xyzp-xyzq)-&
               &(3.*fact3)*((riq2*distip).dot.(rip2*distiq)))
            end do
            fxyz=fxyz+(ffxyz*cost2)
         end do
         
         !!c ciclo sulle cariche di polarizzazione
         do q=1,icount2b
            qq=schrg(q)*epkt
            xyzq=scspos(q); distiq=xyzi-xyzq; riq2=distiq.dot.distiq
            riq=sqrt(riq2); fact=qq/(riq*riq2); fact1=-2.*fact*fact
            
            !!c     termine elettroforetico
            fxyz=fxyz+(qi*fact)*distiq
            !!c     termine dielettroforetico p=q
            ffxyz=fact1*distiq
            !!c     termine dielettroforetico q<>p col simmetrico
            !!c     prima parte: p scorre le cariche vere
            do p=1,nqass
               if (crgatn(p).eq.i) cycle
               qp=atmcrg(p)%value/atmeps(p)
               xyzp=chgpos(p); distip=xyzi-xyzp; rip2=distip.dot.distip
               rip=sqrt(rip2); fact2=qp*fact/(rip*rip2)
               fact3=(distiq.dot.distip)/(riq2*rip2)
               ffxyz=ffxyz+fact2*((xyzp-xyzq)-&
               &(3.*fact3)*((riq2*distip).dot.(rip2*distiq)))
            end do
            !!c     seconda parte, p scorre le cariche di polarizzazione
            do  p=1,q-1
               qp=schrg(p)*epkt
               xyzp=scspos(p); distip=xyzi-xyzp; rip2=distip.dot.distip
               rip=sqrt(rip2); fact2=qp*fact/(rip*rip2)
               fact3=(distiq.dot.distip)/(riq2*rip2)
               ffxyz=ffxyz+fact2*((xyzp-xyzq)-&
               &(3.*fact3)*((riq2*distip).dot.(rip2*distiq)))
            end do
            fxyz=fxyz+ffxyz*cost2
         end do
         rfield(i)=fxyz*cost1
      end do
      end subroutine rforcenew