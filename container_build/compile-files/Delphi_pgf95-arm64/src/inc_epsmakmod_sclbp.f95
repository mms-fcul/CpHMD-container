!#####################################################################
!ATTENTION! This file is part of epsmakmod module.
!           Do not compile separately!
!
!program to reposition the boundary grid points on the molecular surface
!(S. Sridharan         May 1994)
!
!2011-05-19 Other parameters to subroutine are transferred via qlog
!           and pointers modules
!           Parameter iibnum is now local because it can take different
!           values when sclbp is called from msrf and vwtms
!2011-05-19 Arrays declared in pointers module allocated in vwtms 
!           subroutine
!#####################################################################

      subroutine sclbp(vcrd,vnr,iibnum,iab1,iab2)

      implicit none

      integer :: iab1(0:lcb1,0:mcb1,0:ncb1),iab2(0:lcb1,0:mcb1,0:ncb1)
      type(coord) :: vnr(iibnum), vcrd(iibnum)
      integer :: nbra(1000)
      logical :: out,outcb(-2:2,-2:2,-2:2)
      !2011-05-19 Leaving only variables not declared in qlog and
      !           pointers modules
      integer :: iaprec,objecttype,nearest,even,dim1,kind
      integer :: dim,prevmed,med
      integer :: epsdim,imezzo(6),iente(6),ix,iy,iz
      type(int_coord) :: ixyz
      integer iac1,iac2
      logical lga,lgd,iflag,precedenza,vicinanza,flag
      character(96) :: strtmp
      real disev,disodd
      type(coord) :: vectx,vnor,s123
      real tmp(6),radpmax,dst,zeta
      type(coord) :: modul,vectz,mod2,x1xyz,xg, u123,xxyyzz
      real radius,temp
      type(coord) vecty,xq,dot,dxyz,dixyz,dx123, dr
      real alpha,tan2,dx,dy,dz
      real dix,diy,diz,hgs,ds2min1,ds2min2
      type(int_coord) :: lmncb1,it,jxyz
      !2011-05-27 Declarations added due to IMPLICIT NONE
      integer :: iibnum,iac,ia,i,iacl,iii,ii,jjx,jjy,jjz,jx,jy,jz
      integer :: jzi,jyi,jxi,liml,limu,kk,nnbr,ncbp
      real :: x1,cbln,del,dis2,dis,dist,dmn,dcr,ctf,dmxx,dmn1,dmn2
      real :: cba,ds2,dsr,rmn,rdist,dmx

      !2011-05-19 Array allocated by ordinary F95 allocate statement
      allocate(internal(nobject))
      epsdim=natom+nobject+2; iflag=.false.; iac1=0; iac2=0 

      iall=0
    
      !hgs= half grid spacing
      hgs=1./(2.*scale)

      !2011-05-17 Changed to array operations
      outcb=.true.; outcb(-1:1,-1:1,-1:1)=.false. 

      !convertion from grid to real coordinates(can also use routine 
      !gtoc)
      x1=1.0/scale

      !2011-05-19  Using operations on coord and int_coord type 
      !variables defined in module operators_on_coordinates 
      x1xyz=oldmid-(0.5*x1*(igrid+1))
      radpmax=max(radprb(1),radprb(2))

      if (extot.eq.0.and.radpmax.gt.0.0.and.(nobject.gt.1&
                                                &.or.natom.gt.1)) then
         !find extrema
         !here one should consider the global system (Walter)
         write(6,*)'Scaling routine in action!'

         !2011-05-19  Using operations on coord and int_coord type 
         !variables defined in module operators_on_coordinates 
         cmin=coord(6000.,6000.,6000.);cmax=coord(-6000.,-6000.,-6000.)
         do ii=1,nobject
            cmin=min(cmin,limobject(ii)%min)
            cmax=max(cmax,limobject(ii)%max)
         end do

         !2011-05-19 All other parameters are transfered via modules 
         !qlog and pointers
         call sas(xn1,extot)
      end if

      del=radpmax; del=max(del,1./(2.*scale)); cbln=rdmx+del

      call cubedata(2.0,cbln)

      dim=(lcb+1)*(mcb+1)*(ncb+1)
      allocate(cbn1(dim),cbn2(dim))

      dim1=27
      if ((nobject-numbmol).gt.0) dim1=max(dim,27)
      allocate(cbal(dim1*(natom+nobject-numbmol)))

      !2011-05-19 All parameters transfered via modules qlog and 
      !pointers
      call cube(xn1,radprb(1),cbn1,cbn2)

      ncbp=0

      D500: do i=1,iibnum
         !+per trattare molecole con diversa epsilon++01/2002+
         if (iibnum.ne.ibnumsurf.and.numbmol.gt.1) then
            !2011-05-19 Converted to int_coord derived type
            ixyz=int(vcrd(i))
            ix=ixyz%i; iy=ixyz%j; iz=ixyz%k 

            iflag=.false.

            iente(1)=mod(iepsmp(ix,iy,iz)%i,epsdim)
            iente(2)=mod(iepsmp(ix,iy,iz)%j,epsdim)
            iente(3)=mod(iepsmp(ix,iy,iz)%k,epsdim)
            iente(4)=mod(iepsmp(ix-1,iy,iz)%i,epsdim)
            iente(5)=mod(iepsmp(ix,iy-1,iz)%j,epsdim)
            iente(6)=mod(iepsmp(ix,iy,iz-1)%i,epsdim)
 
            imezzo(1)=iepsmp(ix,iy,iz)%i/epsdim
            imezzo(2)=iepsmp(ix,iy,iz)%j/epsdim
            imezzo(3)=iepsmp(ix,iy,iz)%k/epsdim
            imezzo(4)=iepsmp(ix-1,iy,iz)%i/epsdim
            imezzo(5)=iepsmp(ix,iy-1,iz)%j/epsdim
            imezzo(6)=iepsmp(ix,iy,iz-1)%k/epsdim

            !guardo se ho due molecole con diversa epsilon nel 
            !punto,interno
            if(imezzo(1).ne.imezzo(6).and.imezzo(1)*imezzo(6).ne.0) &
                 &iflag=(iente(1).le.natom+1.and.iente(6).le.natom+1)

            !iflag sarà vero se il bgp è interno ma non con un oggetto
            do ii=2,6
               if(imezzo(ii).ne.imezzo(ii-1).and.&
                                      &imezzo(ii)*imezzo(ii-1).ne.0) &
                 &iflag=iflag.or.&
                    &(iente(ii).le.natom+1.and.iente(ii-1).le.natom+1)
            end do
         end if

         !2011-05-19  Using operations on coord and int_coord type 
         !variables defined in module operators_on_coordinates 
         xg=vcrd(i)*x1+x1xyz

         !find the closest surface atom to the gridpoint
         it=int((xg-xyzo)*cbai)

         dmn=100.; ds2min1=1000.; ds2min2=1000.; prevmed=0; iac=0
         nnbr=0
         lmncb=int_coord(lcb,mcb,ncb)

         if ((it.vorlt.0).or.(it.vorgt.lmncb)) then
            !if the bgp is outside the cube, probably it is due to some 
            !object
            do ii=1,nobject
               strtmp=dataobject(ii,1)
               read(strtmp(16:18),*)kind
               if (strtmp(1:4).ne.'is a'.and.kind.ne.2) then
                  if ((xg.vandle.(limobject(ii)%max+x1)).and.&
                              (xg.vandgt.(limobject(ii)%min-x1))) then
                     nnbr=nnbr+1
                     nbra(nnbr)=ii+natom
                     liml=1
                     limu=0
                  end if
               end if
            end do
            
            if(liml.ne.1.or.limu.ne.0) write(6,*)'bgp close to nothing'
         else 
            !2011-05-19 Changed 1d array to 3d array as in cube 
            !subroutine
            liml=cbn1(it%i+1+(lcb+1)*it%j+(lcb+1)*(mcb+1)*it%k)
            limu=cbn2(it%i+1+(lcb+1)*it%j+(lcb+1)*(mcb+1)*it%k)
         end if
           
         iaprec=0
         
         do kk=liml,limu
            ia=cbal(kk)

            if (ia.le.natom) then
               !b+aggiunto iflag per salvare comunque in atsurf valore 
               !del + vicino (01/02)
               !non sono sicurissimo perche' non ricordo esattmente che 
               !fa poi con prevmed...
               if (iflag) then
                  dx123=xg-xn1(ia)
                  dis2=(dx123.dot.dx123)-&
                                &delphipdb(ia)%rad3*delphipdb(ia)%rad3

                  !dis2, and ds2min are distance**2 from center - 
                  !radius**2 so they can be <0
                  if (dis2.lt.ds2min1) then
                     iac2=iac1; ds2min2=ds2min1; iac1=ia; ds2min1=dis2
                  else if (dis2.le.ds2min2) then 
                     iac2=ia; ds2min2=dis2
                  end if
               else
                  if (ast(ia).eq.0) then
                     nnbr=nnbr+1; nbra(nnbr)=ia
                  end if
               end if
            else
               if (ia.ne.iaprec) then
                  iaprec=ia; nnbr=nnbr+1; nbra(nnbr)=ia
               end if
            end if
         end do

         if (iflag) then
            atsurf(i)=iac1

            if (iac1*iac2.eq.0.or.iac1.eq.iac2) then
               write(6,*)'Problems in Scaling multidielectric Boundary Grid Points'
               stop
            end if
            
            atndx(i)=-1

            !looks like atndx is used to build Delunay surface, so 
            !excluding these bgps
            !2011-05-19 Using operations on coord and int_coord type 
            !variables defined in module operators_on_coordinates 
            dx123=xn1(iac2)-xn1(iac1)
            temp=dx123.dot.dx123
            temp=0.5*(ds2min2-ds2min1)/temp
            vcrd(i)=xg+(temp*dx123); vnr(i)=coord(0.,0.,0.)
            cycle D500
         else
            do ii=1,nnbr
               ia=nbra(ii)
               med=iatmmed(ia)
               lgd=(med.ne.prevmed)

               if (ia.gt.natom) then
                  iii=ia-natom 
                  xq=xg

                  !try to find closest VdW surface, better if it is 
                  !buried
                  !internal is used for the object to which surface the 
                  !bgp is closer
                  call distobj(xq,dist,dixyz,iii,0.0,.false.)

                  precedenza=ia.gt.iac.and.(iac.gt.natom.or.iac.eq.0)
                  vicinanza=abs(dist).lt.abs(dmn) 
                  lga=(precedenza.and.(vicinanza.or.dist.lt.0.)).or.&
                                            &(vicinanza.and.dmn.gt.0.)

                  if ((dist.lt.dmn.and..not.lgd).or.(lga.and.lgd)) then
                     dmn=dist; iac=ia; prevmed=med
                     dr=dixyz*(dist-radprb(1)); vnor=dixyz
                  end if

                  internal(iii)=(dist.lt.0.0)
               else
                  dx123=xg-xn1(ia) 
                  dis=sqrt(dx123.dot.dx123)-delphipdb(ia)%rad3
                  precedenza=ia.gt.iac.or.iac.gt.natom
                  vicinanza=abs(dis).lt.abs(dmn)
                  lga=(precedenza.and.(vicinanza.or.dis.lt.0.)).or.&
                                            &(vicinanza.and.dmn.gt.0.)

                  if ((dis.lt.dmn.and..not.lgd).or.(lga.and.lgd)) then
                     prevmed=med
                     dmn=dis; iac=ia
                  end if
               end if
            end do

            atsurf(i)=iac

         end if

         if (iac.eq.0.and.iac1.eq.0) then
            write(6,*)'no close atom or object for boundary point ',i
            stop
         end if

         !if iac is an object dr has alredy been calculated
         !and HAS a DIFFERENT value!!!!!!
         if (iac.le.natom) then
            dr=xg-xn1(iac)
         end if

         dsr=sqrt(dr.dot.dr)
         out=.true.

         if (radpmax.gt.0.0)then
            !u should have the same value as previous one
            if (iac.le.natom) then
               u123=xn1(iac)+(((r0(iac)*dr)/dsr))
            else
               u123=xg-dr
            end if

            it=int((u123-xyzo)*cbai)
            nnbr=0

            !2011-05-19 Changed 1d to 3d array as in cube subroutine
            liml=cbn1(it%i+1+(lcb+1)*it%j+(lcb+1)*(mcb+1)*it%k)
            limu=cbn2(it%i+1+(lcb+1)*it%j+(lcb+1)*(mcb+1)*it%k)

            do kk=liml,limu
               ia=cbal(kk)
               if (ia.le.natom) then
                  dx123=u123-xn1(ia); ds2=dx123.dot.dx123
                  if(ds2.lt.rs2(ia))out=.false.
               else
                  if (ia.ne.iac.and.(.not.internal(ia-natom))) then
                     !I want to know if u is within the shell 
                     !sorrounding the object
                     xq=u123

                     call distobj(xq,dist,dixyz,ia-natom,0.,.true.)

                     if (dist.gt.0.0.and.dist.lt.radprb(1)-1.e-6) &
                                                         & out=.false.
                  end if
               end if
            end do
         end if

         if (out) then
            ncbp=ncbp+1
            if (iac.le.natom) then
               vcrd(i)=xn1(iac)+(dr*(delphipdb(iac)%rad3/dsr))
               vnr(i)=dr/dsr
            else
               vcrd(i)=(xg-(radprb(1)*vnor))-dr
               vnr(i)=vnor
            end if

            atndx(i)=iac
         else
            atndx(i)=0
         end if
      end do D500

      if(allocated(cbn1)) deallocate(cbn1)
      if(allocated(cbn2)) deallocate(cbn2)
      if(allocated(cbal)) deallocate(cbal)

      !scale the re-entrant points with respect to expos
      !if radprb = 0.0 we are done.
      if (radpmax.gt.0.0) then
         iall=0
         cba=1./grdi

         D700: do i=1,iibnum
            !b+++mol. con diversa eps ++++01/02+++++++++++++++
            if(atndx(i).eq.-1) cycle D700

            if (atndx(i).eq.0) then
               s123=vcrd(i)*x1+x1xyz
       
               !2011-05-19 mn(x,y,z) and grdi were assigned values in 
               !INDVER subroutine
               !now coord type variable mnxyz is declared in pointers 
               !module and real grdi declared and thus accessible in 
               !qlog module
               xxyyzz=(s123-mnxyz)*grdi; jxyz=int(xxyyzz)
               jx=jxyz%i; jy=jxyz%j; jz=jxyz%k

               dxyz=xxyyzz-float(jxyz)

               dmn1=min(dxyz%x,dxyz%y,dxyz%z)
               dmx=max(dxyz%x,dxyz%y,dxyz%z)
               dmn2=1.0-dmx; dcr=min(dmn1,dmn2); ctf=cba*(1+dcr)

               ctf=ctf*ctf;  iacl=0;  rmn=100.
               do jjx=jx-1,jx+1
                  do jjy=jy-1,jy+1
                     do jjz=jz-1,jz+1
                        do ii=iab1(jjx,jjy,jjz),iab2(jjx,jjy,jjz)
                           iac= icume(ii)
                           dist=(s123-expos(iac)).dot.(s123-expos(iac))

                           if (dist.lt.rmn) then
                              rmn=dist
                              iacl=iac
                           end if
                        end do
                     end do
                  end do
               end do
        
               if (.not.(iacl.gt.0.and.rmn.lt.ctf)) then
                  do jxi=-2,2
                     do jyi=-2,2
                        do jzi=-2,2
                           if (outcb(jxi,jyi,jzi)) then
                              jjx=jx+jxi
                              if (jjx.ge.0.and.jjx.le.lcb1) then
                                 jjy=jy+jyi
                                 if (jjy.ge.0.and.jjy.le.mcb1) then
                                    jjz=jz+jzi
                                    if (jjz.ge.0.and.jjz.le.ncb1) then
                                       do ii=iab1(jjx,jjy,jjz),&
                                                    &iab2(jjx,jjy,jjz)
                                          iac= icume(ii)
                                          dist=(s123-expos(iac))&
                                               &.dot.(s123-expos(iac))
                                          if (dist.lt.rmn) then
                                             rmn=dist; iacl=iac
                                          end if
                                       end do
                                    end if
                                 end if
                              end if
                           end if
                        end do
                     end do
                  end do
                   
                  if (iacl.le.0) then
                     iall=iall+1
                     do iac=1,extot
                        dist=(s123-expos(iac)).dot.(s123-expos(iac))
                        if (dist.lt.rmn) then
                           rmn=dist; iacl=iac
                        end if
                     end do
                  end if
               end if

               dxyz=s123-expos(iacl); rdist=sqrt(dxyz.dot.dxyz)

               if (rdist.eq.0) then
                  dist=0.0
               else
                  !if inside any object, radprb(2)...
                  dst=0.; flag=.true.
                  
                  D400: do ii=1,nobject
                     strtmp=dataobject(ii,1)
                     read(strtmp(16:18),*)kind
                  
                     if (strtmp(1:4).ne.'is a'.and.kind.ne.2) then
                        xq=s123

                        call distobj(xq,dist,dixyz,ii,0.,.true.)

                        !assuming that if the VdW point is half grid 
                        !space into an object that means that this 
                        !belongs to an atom buried in the object
                        if (dst.lt.-hgs) then 
                           dist=radprb(2)/rdist
                           flag=.false.; exit D400
                        end if
                     end if
                  end do D400

                  if(flag) dist=radprb(1)/rdist
               end if

               vcrd(i)=expos(iacl)+(dxyz*dist)

               if (rdist.gt.1.0e-8) then
                  vnr(i)=(-dxyz)/rdist
               else
                  write(6,*)'bdp close to arcp ',i,rdist
               end if
            end if
         end do D700
      end if
      
      if (allocated(internal)) deallocate(internal)

      !do i=1,ibnum
      !   write(2,'(i8,3f12.6,i5)') &
      !                      &i,vcrd(1,i),vcrd(1,i),vcrd(1,i),atndx(i)
      !end do
      !close (2)
      !Varrebbe la pena di capire % of ... cosa significa.
      !write(6,*)'% of boundary points contacting solvent = ',&
      !                                 &float(ncbp)/float(ibnum)*100.

      end subroutine sclbp
