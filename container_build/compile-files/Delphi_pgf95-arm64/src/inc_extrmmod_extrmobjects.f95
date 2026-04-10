      !  ATTENTION!  This file is part of extrmmod module.
      !==================================================================
      
      subroutine extrmobjects    
      !   2011-05-06  Declared in qlog module             
      !   2011-05-06   Dataobject and limobject declared in pointers module and
      !   2011-05-06  Are parts of derived-type array delphipdb declared in 
      !------------------------------------------------------------------------------
      !  2011-05-06   Some local variables             
      integer :: count,objecttype
      !  2011-05-06  Using new derived type variables instead of arrays
      type(coord) :: xa,xb,xc,xd,tmpvect,tmpvect1,tmpvect2
      type(coord) :: emin, emax 
      real :: alpha,sign
      real :: modul,modul2,tmp1, tmp2
      !------------------------------------------------------------------------------ 
      ! 2011-06-14 Declarations added due to IMPLICIT NONE
      integer :: ix
      real :: radius,tmp
      !-------------------------------------------------------------------
      !!c  limobject contains extreme values of each object 
      !!c  for a molecule it has extreme but without radii
      !!c  descritpion:   limobject(object number,coord,min_max)
      
      !!c find extrema of each object and  according to them 
      
      numbmol=0
      rdmx=0.01
      do count=1,nobject
         !!c extracts object type number
         strtmp=dataobject(count,1)
         if (strtmp(1:4).eq.'is a') then          
            !!c here we have a molecule, calculating also max radius
            numbmol=numbmol+1
            if(verbose) &
            & write(6,*) 'Object number',count,' is a molecule'    
            if (natom.eq.0) then
               write(6,*)'A molecule with no atoms?'
               stop
            end if
            ! 2011-05-06  Limobject is now derived-type array described
            emin=delphipdb(1)%xyz
            emax=delphipdb(1)%xyz
            rdmx=delphipdb(1)%rad3
            do ix=2,natom
               !  2011-05-06  Using operations on coord type variables  defined
               emin=min(emin,delphipdb(ix)%xyz)
               emax=max(emax,delphipdb(ix)%xyz)
               rdmx=max(rdmx,delphipdb(ix)%rad3)
            end do
            limobject(count)%min=emin
            limobject(count)%max=emax
            !  2011-05-06   Label 2103 replaced by end do
            
         else
            read(strtmp(8:10),'(I3)')objecttype
            !!c now deals with different types of objects
            !  2011-05-06 Multiple IF's replaced y select case
            
            select case(objecttype)
            case(1)        ! if (objecttype.eq.1) then
               !c! here we have a sphere
               if(verbose)write(6,*)&
               &  'Object number',count,' is a sphere.'      
               read(strtmp(20:80),*)xa,radius
               !  2011-05-06  Subroutine warning at the end of this module just
               if (radius.le.1./scale) call warning
               !  2011-05-06  Using operations on coord type variables  defined
               limobject(count)%min=xa-radius
               limobject(count)%max=xa+radius
               
            case(2)        ! if (objecttype.eq.2) then
               
               !!c here we have a cylinder
               if(verbose) then
                  write(6,*)'Object number',count,' is a cylinder.'
                  write(6,*)strtmp
               end if
               read(strtmp(20:80),*)xa,xb,radius 
               
               !  2011-05-06  Using operations on coord type variables defined
               tmpvect=xa-xb
               modul2=tmpvect.dot.tmpvect
               modul=sqrt(modul2)
               write(dataobject(count,2),'(5f8.3)')tmpvect,modul, modul2
               ! 2011-05-06   Message replaced by call to subroutine at the end of this module
               if((radius.le.1./scale).or.(modul2.le.1./(scale**2)))&
               & call warning
               tmpvect1=tmpvect*tmpvect
               tmpvect2=radius*sqrt(modul2-tmpvect1)/modul
               limobject(count)%min=xb-tmpvect2+min(0.,tmpvect)
               limobject(count)%max=xb+tmpvect2+max(0.,tmpvect)
               
               
            case(3)       ! if (objecttype.eq.3) then
               !!c here we have a cone
               if(verbose) then
                  write(6,*)'object',count,' is a cone.'
                  write(6,*)strtmp
               end if
               read(strtmp(20:80),*)xa,xb,alpha 
               alpha=tan(alpha*1.5707/180.)
               
               !  2011-05-07  Using operations on coord type variables  defined
               tmpvect=xa-xb
               modul2=tmpvect.dot.tmpvect
               
               write(dataobject(count,2),'(5f8.3)')tmpvect,&
               sqrt(modul2),modul2
               ! 2011-05-07   Message replaced by call to subroutine at the end of this module
               if (modul2.le.1./(scale**2)) call warning
               ! 2011-05-07  Using operations on coord type variables  defined
               tmpvect1=tmpvect*tmpvect
               tmpvect2=alpha*sqrt(modul2-tmpvect1)
               limobject(count)%min=min(xb,xa-tmpvect2)
               limobject(count)%max=max(xb,xa+tmpvect2)
               
            case(4)      !     if (objecttype.eq.4) then
               !!c here we have a box 
               if(verbose) write(6,*)&
               &  'object number',count,' is a box.'
               read(strtmp(20:80),*)xa,xb,xc,xd
               ! 2011-05-07  Using operations on coord type variables  defined
               tmpvect=xd-xa
               tmp=tmpvect.dot.tmpvect
               ! 2011-05-07   Message replaced by call to subroutine at the end of this module
               if (tmp.le.1./(scale**2)) call warning
               ! 2011-05-07  Using operations on coord type variables  defined
               tmpvect2=xc-xa
               tmp2=tmpvect2.dot.tmpvect2
               ! 2011-05-07   Message replaced by call to subroutine at the end of this module
               if (tmp2.le.1./(scale**2)) call warning
               ! 2011-05-07  Using operations on coord type variables  defined
               tmpvect1=xb-xa
               tmp1=tmpvect1.dot.tmpvect1
               ! 2011-05-07   Message replaced by call to subroutine at the end of this module
               if (tmp1.le.1./(scale**2)) call warning
               write(dataobject(count,2),'(12f8.3)')tmpvect,sqrt(tmp),&
               &  tmpvect2,sqrt(tmp2),tmpvect1,sqrt(tmp1)
               write(6,*)"Vertex Coordinates: ",dataobject(count,1)
               write(6,*)"modulus shortest direction:",sqrt(tmp)
               write(6,*)"modulus second   direction:",sqrt(tmp2)
               write(6,*)"modulus third    direction:",sqrt(tmp1)
               ! 2011-05-07  Using operations on coord type variables  defined
               xb=0.5*(xb+xc)
               xa=xb+tmpvect
               
               
               ! 2011-05-07  Using new operations on coord type variables  defined
               limobject(count)%min=  &
               &   submin(xa,xb,tmpvect)+minsign(0.5,tmpvect1)+&
               &   minsign(0.5,tmpvect2)
               limobject(count)%max=  &
               &   submax(xa,xb,tmpvect)+maxsign(0.5,tmpvect1)+&
               &   maxsign(0.5,tmpvect2)
               !!c now tmpvect=A-B; tmpvect1=C-D ; tmpvect2=E-C (axial notation)
               !!c axial changes only A->(B+C)/2+D-A and B->(B+C)/2
               !!c find the minimum
               ! find the maximum
            end select
            
         end if
      end do
      end subroutine extrmobjects
