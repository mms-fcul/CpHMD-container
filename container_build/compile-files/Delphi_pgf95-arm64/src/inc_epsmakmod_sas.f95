!#####################################################################
!ATTENTION! This file is part of epsmakmod module.
!           Do not compile separately!
!
!2011-05-26 Other parameters transfered via module architecture
!#####################################################################

      subroutine sas(crd,nacc)

      implicit none

      !coi contains circle of intersection data for pairs when
      !objects are involved
      type ccoi
         type(coord) :: xyz !vector applied in the center of the sphere 
                            !and pointing to the center of the coi
         real :: rad        !coi radius
         integer :: is      !number of the intersections between an 
                            !atom and an object
      end type ccoi

      integer, parameter :: nver=520,nedge=1040
      !2011-05-26 Arrays declared in pointers module and allocated in 
      !vwtms sub
      type(coord) :: crd(natom), ver(nver),tv
      integer :: edgv(2,nedge),edg(nedge),oti(nver),st(nedge)
      !a third field in pls takes into account the correspondence 
      !nprobj-npr
      !2011-05-26 Array declared in pointers module and allocated here
      !           integer*4 pls(3,1)
      !2011-05-26 Removed variables declared in qlog module and 
      !redeclared some variables as coord type
      integer :: objecttype,itmp,jprec,kk,nprtobj,nprobj,inter
      integer :: ii,jj,nside,dim,dim1,kind
      character(96) :: strtmp,strtmp1
      real :: side,dist
      type(coord) :: omin,omax,vectx,xmin,xmax,xa,xb, xc, xd, xq
      type(coord) :: vectz, vecty, xyzm,x123, dxyz, dx123, tij123
      type(coord) :: rmv(3),cf123, dy123
      type(int_coord) :: ix123, ic123
      type(ccoi), allocatable :: coi(:),coitemp(:)
      type(int_coord), allocatable :: plstemp(:)
      type(coord), allocatable :: expostemp(:)
      real :: tmp,tmp1
      real :: radius,modul,mod2
      real :: alpha,dot
      real :: dx,dy,dz,seno,cose,modx,mody,modul2
      real :: tmp3,rad2
      !2011-05-26 Declaration due IMPLICIT NONE
      integer :: nacc,nacct,i,j,k,ie,ia2,ie2,ie1,ilvl,ia1,iv,iv1,iv2
      integer :: ip,liml,limu,ne,nbv,nlvl,npr,nprp,nprx,nst,nprt,nvo
      integer :: nxa,nvi, nv
      real :: ctf,ctf2,cst,dctf,d2,del,dx1,dx2,dx3,h,ds2,pre,rad, radj
      real :: rij,rdn,rv1,rv2,rvmg,rsx2,dctf2,sm1,sm2,snt,tmp2,tta,vmg
      real :: cbln,csp,dmg,rdx2,tm

      ver=coord(0.,0.,0.); edgv=0.
      nacct=0; nprt=0; nprtobj=0; nlvl=5; nvi=12

      radpmax=max(radprb(1),radprb(2))
      write(6,*)'radpmax,radprb ',radpmax,radprb(1),radprb(2)

      cbln=2.*(rdmx+radpmax)

      call cubedata(1.0,cbln)

      dim=(lcb+1)*(mcb+1)*(ncb+1)

      !2011-05-26 Arrays are allocated with ordinary ALLOCATE F95 
      !           statement
      !           WARNING! Allocation may be incorrect since memalloc 
      !           keeps old values of the array in memory, if exist, 
      !           while ALLOCATE statement kills all old values. Arrays 
      !           are allocated again later in sas subroutines and in 
      !           vwtms subroutine 
      allocate(cbn1(dim),cbn2(dim))

      dim1=27
      if ((nobject-numbmol).gt.0) dim1=max(dim,27)

      allocate(cbal(dim1*(natom+nobject-numbmol)))

      call cube(crd,radprb(1),cbn1,cbn2)

      !generate a template of vertices....
      tta=2.*pi/nvi
      do i=1,nvi
         rdn=(i-1)*tta
         ver(i)=coord(cos(rdn),sin(rdn),0.)
         j=i+1
         if(i.eq.nvi)j=1
         edgv(1,i)=i; edgv(2,i)=j
      end do

      nv=nvi; ne=nvi; ie1=1; ie2=0

      do ilvl=1,nlvl
         ie1=ie2+1; ie2=ne
         do ie=ie1,ie2
            iv1=edgv(1,ie); iv2=edgv(2,ie)

            !2011-05-26  Using operations on coord and int_coord type 
            !variables defined in module operators_on_coordinates 
            xyzm=ver(iv1)+ver(iv2); vmg=sqrt(xyzm.dot.xyzm)

            nv=nv+1
            ver(nv)=xyzm/vmg

            ne=ne+1; edg(ie)=ne; edgv(1,ne)=iv1; edgv(2,ne)=nv
            ne=ne+1; edgv(1,ne)=nv; edgv(2,ne)=iv2
         end do
      end do
      
      ne=ie2
      edg(ie1:ne)=-1
      write(6,*)'# of vertices = ',nv,' # of edges = ',ne
      ast=1

      nacc=0;  npr=0;  nprobj=0;  nprp=0
      !it finds pairs.......
      do i=1,natom
         rad=r0(i); rad2=r02(i)
         if(delphipdb(i)%rad3.eq.0.) cycle
         x123=crd(i); ix123=int((x123-xyzo)*cbai)
         liml=cbn1(ix123%i+1+(lcb+1)*ix123%j+(lcb+1)*(mcb+1)*ix123%k)
         limu=cbn2(ix123%i+1+(lcb+1)*ix123%j+(lcb+1)*(mcb+1)*ix123%k)

         !2011-06-17 Resizing of arrays keeping old values intact
         if ((npr+limu-liml+1).gt.nprt) then
            nprt=nprt+5000
            if (allocated(pls)) then
               allocate(plstemp(nprt-5000)) ; plstemp=pls
               deallocate(pls); allocate(pls(nprt))
               pls(1:nprt-5000)=plstemp; deallocate(plstemp)
            else 
               allocate(pls(nprt))
            end if
         end if

         if ((nprobj+limu-liml+1).gt.nprtobj) then
            nprtobj=nprtobj+1000
            if (allocated(coi)) then
               allocate(coitemp(nprtobj-1000)); coitemp=coi
               deallocate(coi); allocate(coi(nprtobj))
               coi(1:nprtobj-1000)=coitemp; deallocate(coitemp)
            else 
               allocate(coi(nprtobj))
            end if
         end if

         jprec=0

         do jj=liml,limu
            j=cbal(jj)
            if (j.le.natom) then
               radj=r0(j)
               if (delphipdb(j)%rad3.gt.0..and.j.gt.i) then
                  ctf=rad+radj; ctf2=ctf*ctf
                  dctf=abs(rad-radj); dctf2=dctf*dctf
                  dx123=crd(j)-x123 ; d2=dx123.dot.dx123
                  del=ctf2-d2

                  if (del.gt.0.01.and.d2.gt.dctf2) then
                     npr=npr+1
                     pls(npr)=int_coord(i,j,0)
                  end if
               end if
            else
               if (j.ne.jprec) then
                  !it finds out if there is intersection between i and 
                  !kk and it generates suitable parameters 
                  !kk= objectnumber
                  kk=j-natom 
                  strtmp=dataobject(kk,1)
                  strtmp1=dataobject(kk,2)
                  inter=0
                  xq=crd(i)
                  read(strtmp(8:10),'(I3)')objecttype

                  select case(objecttype)
                  case(1) !if (objecttype.eq.1) then
                          !dealing with a sphere
                     call distobj(xq,dist,dxyz,kk,radprb(1),.false.)

                     if (abs(dist).lt.r0(i)) then
                        inter=inter+1; npr=npr+1; nprobj=nprobj+1
                        pls(npr)=int_coord(i,j,nprobj)
                        coi(nprobj)=ccoi((-dxyz)*dist,&
                                          &sqrt(r02(i)-dist**2),inter)
                     end if
                  case(2) !if (objecttype.eq.2) then
                          !dealing with a cylinder
                     call distobj(xq,dist,dxyz,kk,radprb(1),.true.)
               
                     if (abs(dist).lt.r0(i)) then
                        read(strtmp(20:80),*)xa,xb,radius
                        read(strtmp1,'(5f8.3)')vectz,modul,modul2
                        tmp=zeta+radprb(1)
                          
                        if (abs(tmp).lt.r0(i))then !side B
                           inter=inter+1; npr=npr+1; nprobj=nprobj+1
                           pls(npr)=int_coord(i,j,nprobj)
                           coi(nprobj)=ccoi(vectz*(tmp/modul),&
                                           &sqrt(r02(i)-tmp**2),inter)
                        end if
                          
                        tmp=zeta-radprb(1)
                          
                        if (abs(tmp).lt.r0(i)) then !side A
                           inter=inter+1; npr=npr+1; nprobj=nprobj+1
                           pls(npr)=int_coord(i,j,nprobj)
                           coi(nprobj)=ccoi(vectz*(tmp/modul),&
                                           &sqrt(r02(i)-tmp**2),inter)
                        end if

                        tmp=axdist-radius-radprb(1)

                        if (abs(tmp).lt.r0(i)) then !lateral,closest, 
                                                    !approximating 
                                                    !with planes
                           if (axdist.eq.0.) then
                              write(6,*)'cannot use planar approximation',i,j
                           else
                              inter=inter+1; npr=npr+1; nprobj=nprobj+1
                              pls(npr)=int_coord(i,j,nprobj)
                              coi(nprobj)=ccoi((tmp/axdist)*&
                                         &(xq-xb-vectz*(zeta/modul)),&
                                         &sqrt(r02(i)-tmp**2),inter)
                           end if
                        end if
                          
                        tmp=axdist+radius+radprb(1)
                          
                        if (abs(tmp).lt.r0(i)) then !lateral,farthest, 
                                                    !approximating 
                                                    !with planes
                           write(6,*)'planar approx very poor, sphere too big',i,j
                           if (axdist.eq.0.) then
                              write(6,*)'cannot use planar approximation',i,j
                           else
                              inter=inter+1;npr=npr+1;nprobj=nprobj+1
                              pls(npr)=int_coord(i,j,nprobj)
                              coi(nprobj)=ccoi((tmp/axdist)*&
                                         &(xq-xb-vectz*(zeta/modul)),&
                                         &sqrt(r02(i)-tmp**2),inter)
                           end if
                        end if
                     end if
                  case(3) !if (objecttype.eq.3) then
                          !dealing with a cone
                     call distobj(xq,dist,dxyz,kk,radprb(1),.true.)

                     if (abs(dist).lt.r0(i)) then
                        read(strtmp(20:80),*)xa,xb,alpha
                        alpha=alpha*pi/180
                        seno=sin(alpha); cose=cos(alpha)

                        read(strtmp1,'(5f8.3)')vectz,modul,modul2
                        tmp=zeta+radprb(1)

                        if (abs(tmp).lt.r0(i)) then !side B
                           inter=inter+1; npr=npr+1; nprobj=nprobj+1
                           pls(npr)=int_coord(i,j,nprobj)
                           coi(nprobj)=ccoi(vectz*(tmp/modul),&
                                           &sqrt(r02(i)-tmp**2),inter)
                        end if

                        tmp=axdist*cose-radprb(1)-seno*(modul-zeta)
 
                        if (abs(tmp).lt.r0(i)) then !lateral,closest, 
                                                    !approximating 
                                                    !with planes
                           if (axdist.eq.0.) then
                              write(6,*)'no planar approx since sphere too large',i,j
                           else
                              inter=inter+1;npr=npr+1;nprobj=nprobj+1
                              tmp1=tmp*cose/axdist
                              tv=tmp1*(xb-xq+vectz*((zeta-&
                                            &axdist*seno/cose)/modul))
                              pls(npr)=int_coord(i,j,nprobj)
                              coi(nprobj)=ccoi(tv,&
                                           &sqrt(r02(i)-tmp**2),inter)
                           end if
                        end if

                        tmp=seno*(modul-zeta)+cose*axdist+radprb(1)

                        if (tmp.lt.r0(i)) then
                           if (zeta.gt.modul+radprb(1)&
                                               &.or.axdist.eq.0.) then
                              write(6,*)'cannot use planar approx in this pos',i,j
                           else !lateral,farthest, approximating with 
                                !planes
                              inter=inter+1;npr=npr+1;nprobj=nprobj+1
                              pls(npr)=int_coord(i,j,nprobj)
                              tv=(tmp/axdist)*(-cose*(xq-xb)+&
                               &((cose*zeta+axdist*seno)/modul)*vectz)
                              coi(nprobj)=ccoi(tv,&
                                           &sqrt(r02(i)-tmp**2),inter)
                           end if
                        end if

                        !beyond the tip, close to the axis
                        if (zeta.gt.modul+radprb(1)&
                                          &.and.axdist.lt.1.0e-6) then
                           inter=inter+1; npr=npr+1; nprobj=nprobj+1
                           tmp=(dist-r0(i))*(1.-cose)/modul
                           pls(npr)=int_coord(i,j,nprobj)
                           coi(nprobj)=ccoi(tmp*vectz,r0(i)*seno,inter)
                        end if
                     end if
                  case(4) !if (objecttype.eq.4) then
                          !dealing with a parallelepiped
                     read(strtmp(20:80),*)xa,xb,xc,xd
                     read(strtmp1,'(12f6.2)')vectz,modul,vecty,mody,&
                                                           &vectx,modx

                     !now using the newest notation for box vertices
                     !new notation:vectx=B-A,vecty=C-A,vectz=D-A;
                     !             xq=P-A;
                     !chiamate da togliere!!!!!!!!!!!!!!!!!!!!!!!!!!
                     xq=xq-xa

                     !working on z interval
                     dot=vectz.dot.xq 
                     tmp=dot/modul; tmp1=radprb(1)+modul

                     if (tmp.gt.-radprb(1)-rad.and.tmp.lt.rad+tmp1) then
                        tmp2=tmp1-tmp; rdx2=rad2-tmp2**2
                        tmp3=radprb(1)+tmp; rsx2=rad2-tmp3**2

                        if (rdx2.gt.0.0) then
                           inter=inter+1; npr=npr+1; nprobj=nprobj+1
                           tmp2=tmp2/modul
                           pls(npr)=int_coord(i,j,nprobj)
                           coi(nprobj)=ccoi(vectz*tmp2,sqrt(rdx2),inter)
                        end if
                          
                        if (rsx2.gt.0.0) then
                           inter=inter+1; npr=npr+1; nprobj=nprobj+1
                           tmp3=-tmp3/modul
                           pls(npr)=int_coord(i,j,nprobj)
                           coi(nprobj)=ccoi(vectz*tmp3,sqrt(rsx2),inter)
                        end if

                        !working on y interval
                        dot=vecty.dot.xq 
                        tmp=dot/mody; tmp1=radprb(1)+mody
                        tmp2=tmp1-tmp
                        rdx2=rad2-tmp2**2; tmp3=radprb(1)+tmp
                        rsx2=rad2-tmp3**2
                          
                        if (rdx2.gt.0.0) then
                           inter=inter+1; npr=npr+1; nprobj=nprobj+1
                           pls(npr)=int_coord(i,j,nprobj)
                           coi(nprobj)=ccoi(vecty*(tmp2/mody),&
                                                    &sqrt(rdx2),inter)
                        end if
                         
                        if (rsx2.gt.0.0) then
                           inter=inter+1; npr=npr+1; nprobj=nprobj+1
                           pls(npr)=int_coord(i,j,nprobj)
                           coi(nprobj)=ccoi(vecty*(-tmp3/mody),&
                                                    &sqrt(rsx2),inter)
                        end if

                        !working on x interval
                        dot=vectx.dot.xq
                        tmp=dot/modx;tmp1=radprb(1)+modx;tmp2=tmp1-tmp
                        rdx2=rad2-tmp2**2; tmp3=radprb(1)+tmp
                        rsx2=rad2-tmp3**2

                        if (rdx2.gt.0.0) then
                           inter=inter+1; npr=npr+1; nprobj=nprobj+1
                           pls(npr)=int_coord(i,j,nprobj)
                           coi(nprobj)=ccoi(vectx*(tmp2/modx),&
                                                    &sqrt(rdx2),inter)
                        end if
                          
                        if (rsx2.gt.0.0) then
                           inter=inter+1; npr=npr+1; nprobj=nprobj+1
                           pls(npr)=int_coord(i,j,nprobj)
                           coi(nprobj)=ccoi(vectx*(-tmp3/modx),&
                                                    &sqrt(rsx2),inter)
                        end if
                     end if
                  end select
               end if 
            end if !end of j.le.natom

            jprec=j

         end do

         if (npr.eq.nprp) then
            ast(i)=0
         end if
           
         nprp=npr
           
      end do

      !2011-05-26 Added deallocation in order to avoid confusion
      !with later allocation of the same arrays
      if(allocated(cbn1)) deallocate(cbn1)
      if(allocated(cbn2)) deallocate(cbn2)
      if(allocated(cbal)) deallocate(cbal)
      write(6,*)'# of pairs = ',npr

      !2011-05-26 Temporarly removed time calculations
      cbln=rdmx+radpmax

      call cubedata(2.0,cbln)

      dim=(lcb+1)*(mcb+1)*(ncb+1)
      allocate(cbn1(dim))
      allocate(cbn2(dim))

      dim1=27
      if ((nobject-numbmol).gt.0) dim1=max(dim,27)

      allocate(cbal(dim1*(natom+nobject-numbmol)))

      call cube(crd,radprb(1),cbn1,cbn2)

      nprx=0
      do ip=1,npr
         i=pls(ip)%i; j=pls(ip)%j

         if (j.le.natom) then
            dx123=crd(j)-crd(i); d2=dx123.dot.dx123; dmg=sqrt(d2)
            pre=1.+(r02(i)-r02(j))/d2
            tij123=crd(i)+((0.5*pre)*dx123)
            rij=0.5*sqrt((r0(i)+r0(j))**2-d2)&
                                       &*sqrt(d2-(r0(i)-r0(j))**2)/dmg
         else
            nprobj=pls(ip)%k; rij=coi(nprobj)%rad

            !pay attention, here dx has a different meaning from 
            !previous one
            dx123=coi(nprobj)%xyz; d2=dx123.dot.dx123
            tij123=crd(i)+dx123
            dmg=sqrt(d2)
         end if

         dx1=dx123%x; dx2=dx123%y; dx3=dx123%z
         rvmg=sqrt(dx1*dx1+dx2*dx2)

         if (rvmg.gt.1.0e-8) then
            rv1=-dx2/rvmg; rv2=dx1/rvmg; cst=dx3/dmg

            snt=sqrt(1.-cst*cst)
            !snt=rvmg/dmg !doesn't lead to any improved performance

            csp=1.0-cst; tm=csp*rv1; sm1=snt*rv1;  sm2=snt*rv2
            rmv(1)=coord(tm*rv1+cst,tm*rv2,sm2)
            rmv(2)=coord(tm*rv2,csp*rv2*rv2+cst,-sm1)
            rmv(3)=coord(-sm2,sm1,cst)
         else
            rmv(1)=coord(1.,0.,0.)
            rmv(2)=coord(0.,1.,0.)
            rmv(3)=coord(0.,0.,1.)
         end if
           
         nvo=0; nbv=0

         !assign memory to expos if needed
         !2011-06-17 Re-sizing array keeping old value
         if ((nacc+nv).gt.nacct) then
            nacct=nacct+1000
            if (allocated(expos)) then
               allocate(expostemp(nacct-1000)); expostemp=expos
               deallocate(expos); allocate(expos(nacct))
               expos(1:nacct-1000)=expostemp; deallocate(expostemp)
            else 
               allocate(expos(nacct))         
            end if
         end if

         D10: do iv=1,nvi
            !+rm(7)*ver(3,iv) has been removed because it is always 
            !zero
            !2011-05-26 Using operations on coord and int_coord type 
            !variables defined in module operators_on_coordinates 
            cf123%x=rmv(1)%x*ver(iv)%x+rmv(1)%y*ver(iv)%y
            cf123%y=rmv(2)%x*ver(iv)%x+rmv(2)%y*ver(iv)%y
            cf123%z=rmv(3)%x*ver(iv)%x+rmv(3)%y*ver(iv)%y
            cf123=tij123+(cf123*rij)

            if (j.gt.natom) then
               inter=coi(nprobj)%is

               if (inter.gt.1) then
                  !if inter>1 we are close to tips in object, thus, 
                  !false vertices might have been previously 
                  !generated, so now, if a vertex is outside the 
                  !object it is fictiously checked for occlusion by 
                  !oti but discarded afterwards
                  !kk= objectnumber
                  kk=j-natom; xq=cf123

                  call distobj(xq,dist,dxyz,kk,radprb(1),.true.)

                  if (dist.gt.5.0e-4) then 
                     oti(iv)=j; cycle D10
                  end if
               end if
            end if

            ic123=int((cf123-xyzo)*cbai)

            liml=cbn1(ic123%i+1+(lcb+1)*ic123%j+&
                                             &(lcb+1)*(mcb+1)*ic123%k)
            limu=cbn2(ic123%i+1+(lcb+1)*ic123%j&
                                            &+(lcb+1)*(mcb+1)*ic123%k)

            D05: do ii=liml,limu
               k=cbal(ii)

               if (k.gt.natom) then
                  oti(iv)=k; cycle D05
               end if

               dy123=crd(k)-cf123; ds2=dy123.dot.dy123

               if (ds2.lt.rs2(k)) then
                  oti(iv)=k; cycle D10
               end if
            end do D05
              
            nvo=nvo+1; nacc=nacc+1; expos(nacc)=cf123; oti(iv)=0

         end do D10
 
         nst=0
 
         if (nlvl.gt.0) then
            do ie=nvi,1,-1
               ia1=oti(edgv(1,ie))
               ia2=oti(edgv(2,ie))
               if(ia1.gt.0.and.ia1.eq.ia2) cycle
               nst=nst+1; st(nst)=ie
            end do
         end if
	
         if (nst.gt.0) then
            D030: do
               ie=st(nst); nst=nst-1 
               ia1=oti(edgv(1,ie)); ia2=oti(edgv(2,ie))

               if ((ia1.gt.natom).or.(ia2.gt.natom)) then
                  if(nst.gt.0) cycle D030
                  exit D030
               end if

               iv=ie+nvi

               !rm(7)*ver(3,iv) has been removed because it is always 
               !zero
               cf123=coord(rmv(1).dot.ver(iv),rmv(2).dot.ver(iv),&
                                                  &rmv(3).dot.ver(iv))
               cf123=tij123+(cf123*rij)

               if (ia1.ne.0) then
                  dy123=crd(ia1)-cf123; ds2=dy123.dot.dy123

                  if (ds2.lt.rs2(ia1)) then
                     oti(iv)=ia1
                     
                     if (edg(ie).gt.0) then
                        nst=nst+1; st(nst)=edg(ie)+1
                     end if
                     
                     if(nst.gt.0) cycle D030
                     exit D030
                  end if
               end if

               if (ia2.ne.0) then
                  dy123=crd(ia2)-cf123; ds2=dy123.dot.dy123

                  if (ds2.lt.rs2(ia2)) then
                     oti(iv)=ia2

                     if (edg(ie).gt.0) then
                        nst=nst+1; st(nst)=edg(ie)
                     end if
                     
                     if(nst.gt.0) cycle D030
                     exit D030
                  end if
               end if

               ic123=int((cf123-xyzo)*cbai)

               liml=cbn1(ic123%i+1+(lcb+1)*ic123%j+&
                                             &(lcb+1)*(mcb+1)*ic123%k)
               limu=cbn2(ic123%i+1+(lcb+1)*ic123%j+&
                                             &(lcb+1)*(mcb+1)*ic123%k)

               D055: do ii=liml,limu
                  k=cbal(ii)

                  if (k.gt.natom) then
                     oti(iv)=k; cycle D055
                  end if

                  dy123=crd(k)-cf123; ds2=dy123.dot.dy123

                  if (ds2.lt.rs2(k)) then
                     oti(iv)=k

                     if (edg(ie).gt.0) then
                        nst=nst+1; st(nst)=edg(ie)+1
                        nst=nst+1; st(nst)=edg(ie)
                     end if
                         
                     if(nst.gt.0) cycle D030
                     exit D030
                  end if
               end do D055

               nvo=nvo+1; nacc=nacc+1; expos(nacc)=cf123

               oti(iv)=0
               if (edg(ie).gt.0) then
                  if (edg(edg(ie)+1).gt.0.or.ia2.gt.0) then
                     nst=nst+1; st(nst)=edg(ie)+1
                  end if

                  if (edg(edg(ie)).gt.0.or.ia1.gt.0) then
                     nst=nst+1; st(nst)=edg(ie)
                  end if
               end if

               if(nst.le.0) exit D030
            end do D030
         end if

         if (nvo.gt.0) then
            !considering pairs also where one 'partner' is an object
            nprx=nprx+1; ast(i)=0
            if (j.le.natom) ast(j)=0
         end if
           
      end do
      
      if(allocated(pls))  deallocate(pls)
      if(allocated(coi))  deallocate(coi)
      if(allocated(cbn1)) deallocate(cbn1)
      if(allocated(cbn2)) deallocate(cbn2)
      if(allocated(cbal)) deallocate(cbal)

      nxa=0
      do i=1,natom
         if(ast(i).eq.0)nxa=nxa+1
      end do

      if (nobject-numbmol.gt.1) then
         write(6,*)'now calculating object-object exposed vertices'
         jj=0
         do ii=1,nobject
            strtmp=dataobject(ii,1)
            read(strtmp(16:18),*)kind

            if (strtmp(1:4).ne.'is a'.and.kind.ne.2) then
               if (jj.eq.0) then
                  omin=limobject(ii)%min
                  omax=limobject(ii)%max
               end if

               jj=1
               omin=min(omin,limobject(ii)%min)
               omax=min(omax,limobject(ii)%max)
            end if
         end do

         !make cubedata
         !nside = number of initial subdivisions
         nside=4
         !for objects, only water probes involved

         !OBS: it is advisable to improve the sideinter estimation as a 
         !function of smallest object dimensions and of radprobe as 
         !well, here probably is too small!!
         tmp=max(omax%x-omin%x,omax%y-omin%y)
         tmp=max(tmp,omax%z-omin%z)+2*radprb(1)

         side=tmp/real(nside)
         xmin=omin-(radprb(1)-0.5*side)
         xmax=omax+(radprb(1)-0.5*side)
         xmax=omin+side*(0.5+float(int(0.999+((xmax-xmin)/side))))

         h=1./scale
         sideinter=max(radprb(1)/4.,0.05); sidemin=sideinter/4
         write(6,*)'Generating vertices between objects:'
         write(6,*)'Finite difference grid spacing:',h
         write(6,*)'Initial side:',side
         write(6,*)'Intermediate side:',sideinter
         write(6,*)'Minimum side:',sidemin, xmin,xmax
         nacct=nacc+(side/sidemin)**3
         write(6,*)'threshold for number of exposed vertices:',nacct

         !2011-05-26 Array is already allocated
         !2011-05-26 Other parameters are transfered via qlog module
         call objvertices(nacc,xmin,xmax,side,h, (numbmol.gt.0))
      end if

      write(6,*)'# pairs analyzed (atom-atom and atom-object)= ',npr
      write(6,*)'# exposed pairs (atom-atom and atom-object)= ',nprx
      write(6,*)'no. arc points  = ',nacc
      write(6,*)'no. surface atoms  = ',nxa,' nbur = ',natom-nxa
      !if(nacc.gt.exmax)stop 'nacc limit exceeded'
      
      end subroutine sas
