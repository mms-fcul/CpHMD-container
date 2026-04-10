!#####################################################################
!ATTENTION!  This file is part of epsmakmod module.
!            Do not compile separately!
!
!this routine calculate the outgoing normal vector onto the surface
!(possibly extended by prbrad) 
!be careful, dist is an oriented distance!!supposed to be positive 
!if the point is external to the object
!in this context, the vector d is (G-H)/dist
!where G is the point and H its projection onto the surface
!(Walter, last update August 5, 2001) 
!
!2011-05-16 Parameters mainly are transfered via module architecture and
!           recursive subroutines are supported in F95 standard.
!           Parameters below may differ from call to call thus
!           are transferred traditionally via formal arguments.
!           xp - vector transfered from calling program
!           dist - calculated distance (output)
!           ii - number of object (input)
!           prbrad - radius of probe sphere (input)
!           inout - branching parameter from the caling program (input)
!
!2011-05-16  Array is assessible via module archtecture
!#####################################################################

      recursive subroutine distobj(xp,dist,dxyz,ii,prbrad,inout)
      
      implicit none

      !2011-05-16 Nobject, count are not used in this subroutine
      integer :: objecttype,ii,sign
      logical :: inout,imovepoint
      character(96) :: strtmp,strtmp1
      !2011-05-16 All real (3) arrays are converted into coord type 
      !           variable described in pointers module
      type(coord) :: xa,xb,xc,xd,xp,vectx,vecty,vectz,dxyz
      !2011-05-16 Non-implicit real variables, thus explicit 
      !           declaration           
      real :: modul,modul2,modx,mody
      !2011-05-27 Declarations added due to IMPLICIT NONE
      real :: dist,prbrad,calpha,alpha,dot,czeta,fact,extent,prx,pry,prz
      real :: radius,pntsh,talpha,tang,tmp,tmp1,tmp2,tet,xxx
      integer :: i

      !inout is option to calculate only dist and not dx,dy,dz
      strtmp=dataobject(ii,1)
      strtmp1=dataobject(ii,2)
      read(strtmp(8:10),'(I3)')objecttype

      extent=0.0
 
      !2011-05-16 Multiple If's replaced by select case statement
      select case(objecttype)
      case(1) !if (objecttype.eq.1) then
              !sphere
         read(strtmp(20:80),*)xb,radius
     
         !2011-05-16  Using operations on coord and int_coord type 
         !variables defined in module operators_on_coordinates 
         dxyz=xp-xb; tmp= dxyz.dot.dxyz

         if (tmp.eq.0.) then 
            !point exactly in the center
            dist=-radius-prbrad; dxyz%x=-1.; return
         end if

         tmp1=sqrt(tmp)
         dist=tmp1-radius-extent
         if (.not.inout) dxyz=dxyz/tmp1
      case(2) !if (objecttype.eq.2) then
              !cylinder
         read(strtmp(20:80),*)xa,xb,radius
         read(strtmp1,'(5f8.3)')vectz,modul,modul2

         !2011-05-16  Using operations on coord and int_coord type 
         !variables defined in module operators_on_coordinates 
         dxyz=xp-xb; dot=dxyz.dot.vectz
         
         !dot=(P-B)(A-B)
         !dot=dx*vectz(1)+dy*vectz(2)+dz*vectz(3)
         zeta=dot/modul

         !tmp=|P-B|**2
         tmp=dxyz.dot.dxyz

         !axdist is |P-B|sin(teta)
         tmp1=tmp-zeta*zeta

         if (tmp1.lt.0.0.and.tmp1.gt.-1.e-3) then
            axdist=0.0
         else
            axdist=sqrt(tmp1)
         end if

         dist=-extent-zeta

         if (.not.inout) then
            if (dist.ge.0.) then
               !inferiore
               if (axdist.gt.radius+extent) then
                  !laterale inferiore
                  tmp1=(radius+extent)/axdist
                  tmp2=(extent+zeta*tmp1)/modul

                  !2011-05-16  Using operations on coord and int_coord 
                  !type variables defined in module 
                  !operators_on_coordinates 
                  dxyz=(dxyz*(1-tmp1))+(vectz*tmp2)
                  dist=sqrt(dxyz.dot.dxyz) ; dxyz=dxyz/dist
               else
                  !assiale inferiore
                  dxyz=(-vectz)/modul
               end if
            else
               dist=zeta-(modul+extent)
               if (dist.ge.0.) then
                  !superiore
                  if (axdist.gt.radius+extent) then
                     !laterale superiore
                     tmp1=(radius+extent)/axdist
                     tmp2=(-extent+zeta*tmp1)/modul-1.
                     dxyz=(dxyz*(1-tmp1))+(vectz*tmp2)
                     dist=sqrt(dxyz.dot.dxyz) ; dxyz=dxyz/dist
                  else
                     !assiale superiore
                     dxyz=vectz/modul
                  end if
               else
                  !z-centrale
                  tmp1=axdist-radius-extent
                  
                  if (tmp1.lt.dist.or.tmp1.lt.-extent-zeta) then
                     !vicini alle basi da dentro
                     if (-dist.le.zeta+extent) then
                        !base superiore
                        dxyz=vectz/modul
                     else
                        !base inferiore
                        dist=-zeta-extent
                        dxyz=(-vectz)/modul
                     end if
                  else
                     !vicino ai lati (laterale), da dentro o da fuori
                     if (axdist.eq.0.) then
                        !point exactly onto the axis
                        dist=-radius-prbrad ; return
                     end if

                     dist=tmp1
                     tmp2=zeta/modul                   
                     dxyz=(dxyz-(vectz*tmp2))/axdist
                  end if
               end if
            end if
         else
            if (dist.ge.0.) then
               !inferiore, assiale inf. gia' a posto
               if (axdist.gt.radius+extent) then
                  !laterale inferiore
                  tmp1=(radius+extent)/axdist
                  tmp2=(extent+zeta*tmp1)/modul

                  !2011-05-16  Using operations on coord and int_coord 
                  !type variables defined in module 
                  !operators_on_coordinates 
                  dxyz=(dxyz*(1-tmp1))+(vectz*tmp2)
                  dist=sqrt(dxyz.dot.dxyz)
               end if
            else
               dist=zeta-(modul+extent)
               if (dist.ge.0.) then
                  !superiore
                  if (axdist.gt.radius+extent) then
                     !laterale superiore, assiale gia' a posto
                     tmp1=(radius+extent)/axdist
                     tmp2=(-extent+zeta*tmp1)/modul-1.
                     dxyz=(dxyz*(1-tmp1))+(vectz*tmp2)
                     dist=sqrt(dxyz.dot.dxyz)
                  end if
               else
                  !z-centrale
                  tmp1=axdist-radius-extent
                  if (tmp1.lt.dist.or.tmp1.lt.-extent-zeta) then
                     !vicini alle basi da dentro, base sup. gia' a posto
                     if (-dist.gt.zeta+extent) dist=-zeta-extent
                  !base inferiore
                  else
                     !vicino ai lati, da dentro o da fuori
                     dist=tmp1
                  end if
               end if
            end if
         end if
      !else di object type 2
      case(3) !if (objecttype.eq.3) then
         !cone, xb is the tip
         read(strtmp(20:80),*)xa,xb,alpha

         !conversion degrees --> radiants
         alpha=alpha*Pi/180
         talpha=tan(alpha*.5)
         read(strtmp1,'(5f8.3)')vectz,modul,modul2
         radius=modul*talpha
         calpha=cos(alpha*.5)
         tang=tan((Pi-alpha)/4.)

         !2011-05-16  Using operations on coord and int_coord type 
         !variables defined in module operators_on_coordinates 
         dxyz=xp-xb
         dot=dxyz.dot.vectz

         !dot=(P-B)(A-B)
         !dot=dx*vectz(1)+dy*vectz(2)+dz*vectz(3)
         zeta=dot/modul

         !czeta= complement to zeta
         czeta=modul-zeta
 
         !tmp=|P-B|**2
         tmp=dxyz.dot.dxyz

         !axdist is |P-B|sin(teta)
         axdist=sqrt(tmp-zeta*zeta)
         xxx=axdist-radius
         
         if (.not.inout) then
            if (xxx.le.0.) then
               if (czeta.lt.0.) then
                  !exactly below the basis (zona 1)
                  dist=-czeta-prbrad; dxyz=vectz/modul
                  return
               else if(czeta.le.tang*(radius-axdist)) then
                  !internal, close to the basis (zona 5)
                  dist=-czeta-prbrad; dxyz=vectz/modul
                  return
               else if(axdist.le.zeta*talpha) then
                  if (axdist.gt.0.) then
                     !internal, close to lateral (zona 6)
                     dist=calpha*(axdist-zeta*talpha)-prbrad
                     tmp1=-calpha/axdist
                     tmp2=calpha*(talpha+zeta/axdist)/modul
                     dxyz=(dxyz*tmp1)+(vectz*tmp2)
                     return
                  else
                     if (zeta.gt.0.) then
                        !exactly on the axis
                        fact=zeta*talpha*calpha*calpha
                        tmp1=modul2-vectz%x*vectz%x

                        if (tmp1.gt.0.)then
                           !vectz not aligned with x axis 
                           tmp1=1./sqrt(tmp1)
                           tmp2=fact*(talpha-vectz%x*tmp1)/modul
                           dxyz=tmp2*vectz 
                           dxyz%x=dxyz%x+fact*modul*tmp1
                           dist=-zeta*talpha*calpha-prbrad
                           return
                        else
                           !vectz aligned with x axis 
                           tmp2=fact*talpha/modul
                           dxyz=vectz*(fact*talpha); dxyz%y=dxyz%y+fact
                           dist=-zeta*talpha*calpha - prbrad
                           return
                        end if
                     else
                        !exactly on tip (apex)
                        dxyz=vectz/modul; dist=-prbrad; return
                     end if
                  end if
               end if
            else if(czeta.lt.0.) then
               !below basis, laterally (part of zona 2)
               tmp1=radius/axdist
  
               !here axdist > 0 for sure, dist probably too
               tmp2=(zeta*tmp1)/modul-1.
               dxyz=(dxyz*(1.-tmp1))+(vectz*tmp2)
               dist=dxyz.dot.dxyz; dxyz=dxyz/dist; dist=dist-prbrad 
               return
            end if

            !tmp1 = |T-B|
            tmp1=calpha*(zeta+axdist*talpha)

            if (tmp1.le.0.) then
               !below the tip (zona 4)
               dist=sqrt(tmp)
               dxyz=dxyz/dist
               dist=dist-prbrad ; return
            else if(tmp1.ge.modul/calpha) then
               !beyond basis, external (remaining of zone 2)
               tmp1=radius/axdist
   
               !here axdist > 0 for sure, dist probably too
               tmp2=(zeta*tmp1)/modul-1.
               dxyz=(dxyz*(1.-tmp1))+(vectz*tmp2)
               dist=sqrt(dxyz.dot.dxyz) ; dxyz=dxyz/dist
               dist=dist-prbrad; return
            else if(axdist/zeta.gt.talpha) then
               !lateral external (zone 3)
               tmp1=calpha/axdist
               tmp2=calpha*(-talpha-zeta/axdist)/modul
               dxyz=(dxyz*tmp1)+(vectz*tmp2)
               dist=calpha*(axdist-zeta*talpha)-prbrad
               return
            end if
         else
            if (xxx.le.0.) then
               if (czeta.le.tang*(radius-axdist)) then
                  !exactly below the basis (zona 1)
                  !internal, close to the basis (zona 5)
                  dist=-czeta-prbrad ; return
               else if(axdist.le.zeta*talpha) then
                  if (axdist.gt.0.) then
                     !internal, close to lateral (zona 6)
                     dist=calpha*(axdist-zeta*talpha)-prbrad; return
                  else
                     if (zeta.ge.0.) then
                        !exactly on the axis
                        dist=-zeta*talpha*calpha - prbrad ; return
                     end if
                  end if
               end if
            else if(czeta.lt.0.) then
               !below basis, laterally (part of zona 2)
               tmp1=radius/axdist

               !here axdist > 0 for sure, dist probably too
               tmp2=(zeta*tmp1)/modul-1.
               dxyz=(dxyz*(1.-tmp1))+(vectz*tmp2)
               dist=sqrt(dxyz.dot.dxyz) -prbrad
               return
            end if

            !tmp1 = |T-B|
            tmp1=calpha*(zeta+axdist*talpha)

            if (tmp1.le.0.) then
               !below the tip (zona 4)
               dist=sqrt(tmp) - prbrad ; return
            else if(tmp1.ge.modul/calpha) then
               !beyond basis, external (remaining of zone 2)
               tmp1=radius/axdist
        
               !here axdist > 0 for sure, dist probably too
               tmp2=(zeta*tmp1)/modul-1.
               dxyz=(dxyz*(1.-tmp1))+(vectz*tmp2)
               dist=sqrt(dxyz.dot.dxyz); dxyz=dxyz/dist
               dist=dist-prbrad; return
            else if(axdist/zeta.gt.talpha) then
               !lateral external (zone 3)
               dist=calpha*(axdist-zeta*talpha) -prbrad; return
            end if
         end if
      !else di objecttype 3
      case(4) !if(objecttype.eq.4) then
              !box
         read(strtmp(20:80),*)xa,xb,xc,xd
         read(strtmp1,'(12f8.3)')vectz,modul,vecty,mody,vectx, modx

         !normalize vectors
         !2011-05-16  Using operations on coord and int_coord type 
         !variables defined in module operators_on_coordinates 
         vectx=vectx/modx ; vecty=vecty/mody ; vectz=vectz/modul

         !calulate projections
         dxyz=xp-xa

         !dot=(P-A)(D-A) new notation
         prz=dxyz.dot.vectz

         !dot=(P-A)(C-A) new notation
         pry=dxyz.dot.vecty

         !dot=(P-A)(B-A) new notation
         prx=dxyz.dot.vectx

         !calculate components of the projection of P onto the extended 
         !surface
         !change temporarily modulus values only for convenience
         modx=modx+extent; mody=mody+extent; modul=modul+extent
         sign=1
         
         if (prx.lt.-extent) then 
            prx=-extent; sign=0
         else
            if (prx.gt.modx) then
               prx=modx; sign=0
            end if
         end if
 
         if (pry.lt.-extent) then
            pry=-extent; sign=0
         else
            if (pry.gt.mody) then
               pry=mody; sign=0
            end if
         end if
         
         if (prz.lt.-extent) then
            prz=-extent; sign=0
         else
            if (prz.gt.modul) then
               prz=modul; sign=0
            end if
         end if

         !now if sign=1 we are inside the box so dist will eventually 
         !be negative
         if (sign.eq.1) then
            i=1; dist=prx+extent
            if (prx+extent.gt.modx-prx) then
               i=-1; dist=modx-prx
            end if
            
            if (pry+extent.lt.dist) then
               i=2; dist=pry+extent
            end if

            if (mody-pry.lt.dist) then
               i=-2; dist=mody-pry
            end if

            if (prz+extent.lt.dist) then
               i=3; dist=prz+extent
            end if

            if (modul-prz.lt.dist) then
               i=-3; dist=modul-prz
            end if

            dist=-dist

            !i moduli li lascio sbagliati perche' non li riuso
            !mody=mody-extent
            !modul=modul-extent

            !2011-05-16 GOTO() statement is replaced by select case 
            if (.not.inout) then
               select case(i+4)
               case(1)
                  dxyz=vectz
               case(2)
                  dxyz=vecty
               case(3)
                  dxyz=vectx
               case(4)
                  write(6,*)'problems in distobject'
                  stop
               case(5)
                  dxyz=-vectx
               case(6)
                  dxyz=-vecty
               case(7)
                  dxyz=-vectz
               case default
                  write(6,*)'problems in distobject'
                  stop
               end select
            end if
         else    
            !i moduli li lascio sbagliati perche' non li riuso
            !modx=modx-extent
            !mody=mody-extent
            !modul=modul-extent

            !2011-05-16  Using operations on coord and int_coord type 
            !variables defined in module operators_on_coordinates 
            dxyz=((prx*vectx)+(pry*vecty)+(prz*vectz))-dxyz
            dist=sqrt(dxyz.dot.dxyz)

            if (.not.inout) then
               if (dist.eq.0.) write(6,*)'Wrong step in distobj!'
               dxyz=(-dxyz)/dist
            end if
         end if
      end select

      dist=dist-prbrad;  return

      !bgp lies on extended surface or center or axdist=0, correcting 
      !2011-05-16 How this part is accessed?
      pntsh=0.
      if (imovepoint) then
         pntsh=1.e-5
         if (vectz%y.eq.0.) tet=1.5707963268
      else 
         tet=atan(vectz%z/vectz%y)
         xp%y=xp%y+pntsh*sin(tet)
         xp%z=xp%z+pntsh*cos(tet)
      end if

      call distobj(xp,dist,dxyz,ii,prbrad+1.e-5-pntsh,inout)

      if (.not.imovepoint) dist=0.0
           
      dist=dist-prbrad

      return
 
      end subroutine distobj
