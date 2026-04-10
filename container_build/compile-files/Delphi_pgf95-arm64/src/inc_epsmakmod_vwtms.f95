!#####################################################################
!ATTENTION! This file is part of epsmakmod module.
!           Do not compile separately!
!
!subroutine to take a van der Waals epsmap and expand it into a
!molecular volume map..
!procedure is to first generate a list of accessible points
!whose collective volume may contain any midpoint in the box.
!next send these points to be mapped into indexing arrays...
!then take the vanderwaals epsmap and grow it out to the molecular
!epsmap, checking against the set of accessible points..
!
!can't do the vw surface by projecting midpoints to accessible surface 
!and testing if that point is in or outside the accessible volume 
!(20 Oct 92)
!(but can do the contact part of MS that way, Nov 93)
!
!program to expand the vw surface into Richards' molecular surface
!(S. Sridharan	May 1994)
!
!2011-05-17 Parameters transferred via modules qlog and pointers
!
!integer cbn1(1),cbn2(1),cbal(1)  !2011-05-16 Arrays declared
!integer iab1(1),iab2(1),icume(1) !in pointers module
!pointer (i_egrid,egrid)  2011-05-16 transfered to pointers.f95
!#####################################################################

      subroutine vwtms

      logical :: outcb(-2:2,-2:2,-2:2)
      character(80) :: line
      integer :: nbra(1000)
      !2011-05-17 bndeps is declared in pointers module           
      integer :: eps(6),nt
      integer :: itmp(6),dim1,cont
      integer iaprec,objecttype,imap(4,6),dim,sign,epsdim,isign
      integer :: kind,eps2(6),j3
      logical :: remov
      character(96) :: strtmp
      type(coord) :: itemp, rtemp,xq,dxyz,dr123,dx123,u123
      type(coord) :: goff(6),xg,x1xyz,s123,xxyyzz
      type(int_coord) :: ixyz,it123,ixyz2,jxyz
      type(int_coord), allocatable :: ibgrd_temp(:)
      logical out,nbe(0:6),intb,exists,flag
      !2011-05-27 Declarations added  due to IMPLICIT NONE           
      real :: offf,cbln,cba,del,dis,dmn,dist,ds2,off,dsr,r0a
      real :: s1,s2,s3,x1,radp,prbrd12,prbrd22,rsm,rsm1,prbrd2
      integer :: ia,i,iac,ibgp,ii,iext,iarv,iacv,iiii,imedia
      integer :: iiord,ios,it1,it2,it3,iv,ix2,iy2,iz2,ix,iy,iz
      integer :: iii,iord,j,jj,jjj,jx,jy,jz,limu,liml,kk,m,k,n
      integer :: mr,mpr,n2,nacc,natm,nmmt,ndv,ncms,ncav,nnn,nnbr
      integer :: nn,ntt,nmt,n1,Nar,iqqq

      !2011-05-17 Arrays allocated by ordinary F95 allocate statement
      allocate(ibnd(ibmx),r0(natom),r02(natom),rs2(natom))
      allocate(ast(natom))
      allocate(bndeps(igrid,igrid,igrid,2))

      offf=(igrid+1.)/2.
      epsdim=natom+nobject+2

      !imap maps from midpoint position to iepsmap entry positions
      !2011-05-17 Changed to array operations
      imap=0

      imap(1,4)=-1
      imap(2,5)=-1
      imap(3,6)=-1
      imap(4,1)=1
      imap(4,2)=2
      imap(4,3)=3
      imap(4,4)=1
      imap(4,5)=2
      imap(4,6)=3

      !2011-05-17 Changed to array operations
      outcb=.true.; outcb(-1:1,-1:1,-1:1)=.false. 
      nbe=.false.; nbe(1:5)=.true.

      goff=coord(0.,0.,0.) 
      off=0.5/scale

      !2011-05-17 Changed to coord derived variable type
      goff(1)%x=off; goff(2)%y=off; goff(3)%z=off
      goff(4)%x=-off; goff(5)%y=-off; goff(6)%z=-off

      radpmax=max(radprb(1),radprb(2)) !transferred via qlog module

      !convertion from grid to real coordinates(can also use routine 
      !gtoc)
      !2011-05-17  Using operations on coord and int_coord type 
      !variables defined in module operators_on_coordinates 
      x1=1.0/scale
      x1xyz=oldmid-(0.5*x1*(igrid+1))

      !find extrema

      !find global extrema
      cmin=coord(6000,6000,6000) ; cmax=coord(-6000,-6000,-6000)
      do ii=1,nobject
         cmin=min(cmin,limobject(ii)%min)
         cmax=max(cmax,limobject(ii)%max)
      end do

      !find vanderwaals boundary
      n=0; nn=0; nmt=0; nmmt=0

      !NB change limits to those of the molecule.
      !set for iepsmp NOT equal to unity
      do k=limeps%min%k+1,limeps%max%k-1
         do j=limeps%min%j+1,limeps%max%j-1
            do i=limeps%min%i+1,limeps%max%i-1
               !one distinguishes between internal,external,
               !internal bgp and external bgp
               iext=0
               ibgp=0

               !2011-05-17 Changed to iepsmp to int_coord derived type
               itmp(1)=iabs(iepsmp(i,j,k)%i)/epsdim
               itmp(2)=iabs(iepsmp(i,j,k)%j)/epsdim
               itmp(3)=iabs(iepsmp(i,j,k)%k)/epsdim
               itmp(4)=iabs(iepsmp(i-1,j,k)%i)/epsdim
               itmp(5)=iabs(iepsmp(i,j-1,k)%j)/epsdim
               itmp(6)=iabs(iepsmp(i,j,k-1)%k)/epsdim

               if(itmp(1).eq.0) iext=1
               if(itmp(1).ne.itmp(6)) ibgp=1
         
               do cont=2,6
                  if(itmp(cont).eq.0) iext=1
                  if(itmp(cont).ne.itmp(cont-1)) ibgp=1
               end do

               !assignement of right values to bndeps according to the 
               !point nature
               !from now ibnum is the total number of internal and 
               !external boundary grid points
               if (ibgp.gt.0) then
                  n=n+1
                  bndeps(i,j,k,1)=n
                  bndeps(i,j,k,2)=iext
                  if (iext.gt.0) nn=nn+1
                  ibnd(n)=int_coord(i,j,k)
               else
                  bndeps(i,j,k,1)=0
                  bndeps(i,j,k,2)=0
               end if

               if (debug) then
                  nt=0
                  !passed from mod to dim, I need the medium
                  if((iepsmp(i,j,k)%i/epsdim).gt.0)nt=nt+1
                  if((iepsmp(i,j,k)%j/epsdim).gt.0)nt=nt+1
                  if((iepsmp(i,j,k)%k/epsdim).gt.0)nt=nt+1
                  if((iepsmp(i-1,j,k)%i/epsdim).gt.0)nt=nt+1
                  if((iepsmp(i,j-1,k)%j/epsdim).gt.0)nt=nt+1
                  if((iepsmp(i,j,k-1)%k/epsdim).gt.0)nt=nt+1

                  if (nbe(nt).neqv.((iext.eq.1).and.(ibgp.eq.1)))then
                     write(6,*)'PROBLEMS1 ',i,j,k
                   
                     !2011-05-17 converts grid to real coordinates 
                     !using operations on coord and int_coord type 
                     !variables defined in module 
                     !operators_on_coordinates
                     itemp=float(int_coord(i,j,k))
                     rtemp=((itemp-offf)/scale)+oldmid

                     write(6,*)rtemp
                     write(6,*)iepsmp(i,j,k),iepsmp(i-1,j,k)%i,&
                                  &iepsmp(i,j-1,k)%j,iepsmp(i,j,k-1)%k
                  end if
               end if
            end do
         end do
      end do

      ibnum=n
      ibnumsurf=nn
      nn=0

      if (verbose) then
         write(6,*)'boundary points facing continuum solvent= ',&
                   &ibnumsurf
         write(6,*)'total number of boundary points before elab.= ',&
                   &ibnum
      end if

      if (ibnum.gt.ibmx) then
         write(6,*)'ibnum= ',ibnum,' is greater than ibmx = ',ibmx
         write(6,*)'increase ibmx in vwtms.f'
         stop
      end if

      if (debug) then
         open(52,file='bgpPrima.log',form='formatted')
         
         !grid coordinates
         write(52,*)igrid,ibnum
         do iiii=1,ibnum
            write(52,*)iiii,ibnd(iiii)
         end do
         close (52)
      end if !end debug

      !2011-05-17 Arrays allocated by ordinary F95 allocate statement
      if (radpmax.lt.1.e-6) then
         allocate(ibgrd(ibnum))

         !2011-05-17 Changed to array operations, but keeping some
         !assignment in a cycle due to array size mismatch
         ast=0
         do i=1,ibnum
            ibgrd(i)=ibnd(i)
         end do
      else
         if (.not.ionlymol.and.radprb(1).ne.radprb(2)) then
            do i=1,natom
               !2011-05-17  Using operations on coord and int_coord 
               !type variables defined in module 
               !operators_on_coordinates 
               xq=xn1(i); r0a=delphipdb(i)%rad3+radprb(1)

               do ii=1,nobject
                  strtmp=dataobject(ii,1)
                  read(strtmp(16:18),*)kind
                  if (strtmp(1:4).ne.'is a'.and.kind.ne.2) then
                     if (((xq-delphipdb(i)%rad3).vandlt.&
                        &limobject(ii)%max).and.((xq+&
                        &delphipdb(i)%rad3).vandgt.&
                        &limobject(ii)%min)) then
                      
                        call distobj(xq,dist,dxyz,ii,0.,.true.)
 
                        !only full buried atoms have the proberadius 
                        !changed to radprb(2)
                        if (dist.lt.-delphipdb(i)%rad3) &
                                     & r0a=delphipdb(i)%rad3+radprb(2)
                     end if
                  end if
               end do
               
               r0(i)=r0a; r02(i)=r0a*r0a; rs2(i)=0.99999*r02(i)
            end do
         else
            do i=1,natom
               r0a=delphipdb(i)%rad3+radprb(1)
               r0(i)=r0a; r02(i)=r0a*r0a; rs2(i)=0.99999*r02(i)
            end do
         end if

         if(iacs) write(6,'(a21)') " opening surface file"
         if(iacs) open(40,file='hsurf2.dat')

         !make a list of accessible points..,expos. all scaling of grid
         !points will be done to thses points..
         prbrd12=radprb(1)*radprb(1)
         prbrd22=radprb(2)*radprb(2)

         !calculate an identity for this conformation
         rsm=0.0
         if (ionlymol) then
            do i=1,natom
               rsm=rsm+delphipdb(i)%rad3*sum(abs(xn1(i)))
            end do
         end if

         inquire(file='ARCDAT',exist=exists)
         flag=.true.
         if (exists.and.ionlymol) then
            open(1,file='ARCDAT',form='unformatted',status='old',&
                &iostat=ios)
            read(1)natm,radp,nacc,rsm1

            !maybe it could be improved taking into account objects 
            !(Walter 1999)
            if (natm.eq.natom.and.radp.eq.radprb(1).and.&
                                       & abs(rsm1-rsm).lt.1.0e-8) then
               if(verbose)write(6,*)'reading accessible surface arcs data from file  ARCDAT'
               extot=nacc
               allocate(expos(extot))

               !2011-05-19 Subroutine arcio is too short and called 
               !only from this subroutine, thus no meaning to be 
               !separate piece of code
               read(1)ast; read(1) expos
               if(verbose) write(6,*)'no. of arc points read = ',nacc
               close (1)
               flag=.false.
            else
               close (1)
            end if
         end if

         if (flag) then
            call sas(xn1,extot)

            if (extot.gt.0.and.ionlymol) then
               !maybe can be improved by taking into account objects
               write(6,*)&
                    &'writing accessible surface arcs data to  ARCDAT'
               open(1,file='ARCDAT',form='unformatted',iostat=ios)
               
               if (ios.ne.0) then
                  write(6,*) 'error opening ARCDAT file'; stop
               end if

               write(1)natom,radprb(1),extot,rsm             

               if (natom.eq.0) then
                  write(1)0; write(1) expos
               else
                  write(1)ast; write(1)expos
               end if

               close(1)
            end if

            if (debug) then
               open(52,file='Vertices.txt',form='formatted')
               do iiii=1,extot
                  write (52,*) expos(iiii)
               end do
               close (52)
            end if !end debug

         end if

         del=1./scale
         del=max(del,radpmax)
         cbln=rdmx+del

         call cubedata(2.0,cbln)

         dim=(lcb+1)*(mcb+1)*(ncb+1)

         !2011-05-17 Array allocation is done by ordinary F95 ALLOCATE 
         !statement Array allocated as 3D as in the cube subroutine
         allocate(cbn1(dim),cbn2(dim))

         dim1=27
         if ((nobject-numbmol).gt.0) dim1=max(dim,27)
         allocate(cbal(dim1*(natom+nobject-numbmol)))

         call cube(xn1,radprb(1),cbn1,cbn2)

         !link the accessible points into iab1 and iab2
         call indverdata(radpmax,scale)

         cba=1./grdi; NNN=(lcb1+1)*(mcb1+1)*(ncb1+1)
         allocate(iab1(NNN),iab2(NNN))
         allocate(icume(extot))

         if(verbose)&
               &write(6,"(a,f8.3)")' grid for indexing accessible points = ',cba

         call indver(extot,iab1,iab2)

         !write out surface data
         if (iacs) then
            line= ' '; line(1:6)='ATOM  '; line(14:14)='O' 
            line(18:20)='SP '
            
            do i=1,extot
               xg=expos(i)
               iv=1
       
               !2011-05-24 Subroutine watput is small and called only 
               !once, thus the code transfered into calling subroutine
               write(line(7:11),'(i5)')i; write(line(24:26),'(i3)')iv

               !2011-05-24 In original coordinates from array xo were 
               !written into th line.
               !Logic of this piece, however, requires coordinates 
               !from xg to be written
               write(line(31:54),'(3(f8.3))')xg
               write(40,'(a80)') line
            end do
            close (40)
         end if

         !now start the expansion
         !m1= the number of boundary points removed from list
         ncav=0
         n1=1; n2=ibnum

         !m= number of new boundary elements..
         mpr=100000; ndv=0

         D100: do
            m=0; mr=0
            
            do i=n1,n2
               ixyz=ibnd(i); ix=ixyz%i; iy=ixyz%j; iz=ixyz%k

               !considering both internal and external b.g.p.
               if (bndeps(ix,iy,iz,1).ne.0) then

                  !still has to be considered what is external and 
                  !what internal!!!!!WWW
                  !remov is true if it is an internal midpoint close 
                  !to an interface where a molecule is present 
                  !(expansion has to take place also in objects)
                  remov=.false.

                  !tengo il mod perche' deve prendere solo punti in 
                  !atomi
                  eps(1)=mod(iepsmp(ix,iy,iz)%i,epsdim)
                  eps(2)=mod(iepsmp(ix,iy,iz)%j,epsdim)
                  eps(3)=mod(iepsmp(ix,iy,iz)%k,epsdim)
                  eps(4)=mod(iepsmp(ix-1,iy,iz)%i,epsdim)
                  eps(5)=mod(iepsmp(ix,iy-1,iz)%j,epsdim)
                  eps(6)=mod(iepsmp(ix,iy,iz-1)%k,epsdim)

                  remov=((eps(1).gt.1.and.eps(1).le.natom+1).or.&
                                 &(eps(2).gt.1.and.eps(2).le.natom+1))
                  remov=((eps(3).gt.1.and.eps(3).le.natom+1).or.&
                        &(eps(4).gt.1.and.eps(4).le.natom+1)).or.remov
                  remov=((eps(5).gt.1.and.eps(5).le.natom+1).or.&
                        &(eps(6).gt.1.and.eps(6).le.natom+1)).or.remov

                  !da farsi solo se pores eps2 contiene il mezzo
                  eps2(1)=(iepsmp(ix,iy,iz)%i/epsdim)
                  eps2(2)=(iepsmp(ix,iy,iz)%j/epsdim)
                  eps2(3)=(iepsmp(ix,iy,iz)%k/epsdim)
                  eps2(4)=(iepsmp(ix-1,iy,iz)%i/epsdim)
                  eps2(5)=(iepsmp(ix,iy-1,iz)%j/epsdim)
                  eps2(6)=(iepsmp(ix,iy,iz-1)%k/epsdim)

                  !cWWW there is still an issue in case there are both 
                  !molecules and objects: since parent object of 
                  !reentrant points is only known in sclbp, filling 
                  !reentrant regions due to molecules in objects might 
                  !fail
                  remov=remov.and.(numbmol.gt.0)

                  !2011-05-17 Using operations on coord and int_coord 
                  !type variables defined in module 
                  !operators_on_coordinates 
                  xg=(x1*float(ixyz))+x1xyz

                  D200: do j=1,6
                     !essere in poro ==> eps2=0 and eps >0
                     if (eps(j).eq.0.or.(remov.and.eps(j).gt.natom+1)&
                        &.or.(eps2(j).eq.0.and.eps(j).gt.0)) then
                        prbrd2=prbrd22
                        if (eps(j).eq.0.or.eps2(j).eq.0) prbrd2=prbrd12

                        !add midpoint offset to grid point..
                        s123=xg+goff(j)

                        !determine if this virgin midpoint is in or out
                        !2011-05-18 mn(x,y,z) and grdi were 
                        !assigned values in INDVER subroutine
                        !now coord type variable mnxyz is declared 
                        !in pointers module and real grdi declared 
                        !and thus accessible in qlog module
                        xxyyzz=(s123-mnxyz)*grdi; jxyz=int(xxyyzz)
                        jx=jxyz%i; jy=jxyz%j; jz=jxyz%k

                        !2011-05-18 Indexes lcb1, mcb1, ncb1 are 
                        !transfred via qlog module and are set in 
                        !INDVER subroutine
                        if ((jxyz.vorle.0).or.(jxyz.vorge.lmncb1)) then
                           write(6,*)'midpoint out of cube' 
                           write(6,'(2i5,3f8.3,3i6,3i8)')&
                                               &i,j,xxyyzz,jxyz,lmncb1
                           write(6,*)iepsmp(ix,iy,iz)%i
                           write(6,*)iepsmp(ix,iy,iz)%j 
                           write(6,*)iepsmp(ix,iy,iz)%k
                           write(6,*)iepsmp(ix-1,iy,iz)%i
                           write(6,*)iepsmp(ix,iy-1,iz)%j
                           write(6,*)iepsmp(ix,iy,iz-1)%k
                        end if
                       
                        dmn=1000.
                        iacv=0

                        !2011-05-18 Repeating piece of the code now 
                        !is in the separate file
                        include 'inc_epsmakmod_vwtms.inc' 
                           
                        !-1,0,0
                        jx=jx-1

                        include 'inc_epsmakmod_vwtms.inc' 
                           
                        !1,0,0
                        jx=jx+2

                        include 'inc_epsmakmod_vwtms.inc' 
                           
                        !0,-1,0
                        jx=jx-1
                        jy=jy-1

                        include 'inc_epsmakmod_vwtms.inc' 
                           
                        !0,1,0
                        jy=jy+2

                        include 'inc_epsmakmod_vwtms.inc' 
                           
                        !0,0,-1
                        jy=jy-1
                        jz=jz-1

                        include 'inc_epsmakmod_vwtms.inc' 
                           
                        !0,0,1
                        jz=jz+2

                        include 'inc_epsmakmod_vwtms.inc' 
                          
                        !nn=2
                        !1,0,1
                        jx=jx+1

                        include 'inc_epsmakmod_vwtms.inc'

                        !-1,0,1
                        jx=jx-2

                        include 'inc_epsmakmod_vwtms.inc'

                        !0,1,1
                        jx=jx+1
                        jy=jy+1

                        include 'inc_epsmakmod_vwtms.inc'
     
                        !0,-1,1
                        jy=jy-2

                        include 'inc_epsmakmod_vwtms.inc'

                        !-1,-1,0
                        jz=jz-1
                        jx=jx-1

                        include 'inc_epsmakmod_vwtms.inc'

                        !1,-1,0
                        jx=jx+2

                        include 'inc_epsmakmod_vwtms.inc'

                        !1,1,0
                        jy=jy+2

                        include 'inc_epsmakmod_vwtms.inc'

                        !-1,1,0
                        jx=jx-2

                        include 'inc_epsmakmod_vwtms.inc'

                        !-1,0,-1
                        jz=jz-1
                        jy=jy-1

                        include 'inc_epsmakmod_vwtms.inc'

                        !1,0,-1
                        jx=jx+2

                        include 'inc_epsmakmod_vwtms.inc'

                        !0,1,-1
                        jx=jx-1
                        jy=jy+1

                        include 'inc_epsmakmod_vwtms.inc'
 
                        !0,-1,-1
                        jy=jy-2

                        include 'inc_epsmakmod_vwtms.inc'

                        !nn=3
                        !-1,-1,-1
                        jx=jx-1

                        include 'inc_epsmakmod_vwtms.inc'

                        !1,-1,-1
                        jx=jx+2

                        include 'inc_epsmakmod_vwtms.inc'

                        !1,1,-1
                        jy=jy+2
                         
                        include 'inc_epsmakmod_vwtms.inc' 
 
                        !-1,1,-1
                        jx=jx-2
                         
                        include 'inc_epsmakmod_vwtms.inc'

                        !-1,1,1
                        jz=jz+2
                         
                        include 'inc_epsmakmod_vwtms.inc'

                        !1,1,1
                        jx=jx+2
                         
                        include 'inc_epsmakmod_vwtms.inc'
 
                        !1,-1,1
                        jy=jy-2
                         
                        include 'inc_epsmakmod_vwtms.inc'

                        !-1,-1,1
                        jx=jx-2

                        include 'inc_epsmakmod_vwtms.inc'

                        !it might be in the contact region; find 
                        !the closest atom surface
                        !2011-05-18  Using operations on coord and 
                        !int_coord type variables defined
                        !in module operators_on_coordinates 
                        it123=int((s123-xyzo)*cbai)

                        dmn=100.; iac=0; nnbr=0
                        lmncb=int_coord(lcb,mcb,ncb)

                        if((it123.vorlt.0).or.(it123.vorgt.lmncb))then
                           !if the bgp is outside the cube, 
                           !probably it is due to some object
                           do ii=nobject,1,-1
                              strtmp=dataobject(ii,1)
                              read(strtmp(16:18),*)kind
                              if(strtmp(1:4).ne.'is a'.and.kind.ne.2) &
                                 &then
                                 if ((s123.vandle.(limobject(ii)%max+&
                                    &x1)).and.(s123.vandgt.&
                                    &(limobject(ii)%min-x1))) then
                                    nnbr=nnbr+1; nbra(nnbr)=ii+natom
                                    liml=0; limu=0
                                 end if
                              end if
                           end do

                           if(liml.ne.0.or.limu.ne.0) &
                                   &write(6,*)'a bgp close to nothing'
                        else
                           !2011-05-19 Changed 1d array to 3d array as 
                           !in cube subroutine
                           liml=cbn1(it123%i+1+(lcb+1)*it123%j+&
                                             &(lcb+1)*(mcb+1)*it123%k)
                           limu=cbn2(it123%i+1+(lcb+1)*it123%j+&
                                             &(lcb+1)*(mcb+1)*it123%k)
                        end if

                        iaprec=0
                        DOKK: do kk=liml,limu
                           ia=cbal(kk)
                           if (ia.eq.0) write(6,*)'problems with cube'
                        
                           if (ia.le.natom.and.ia.gt.0) then
                              if (ast(ia).eq.0) then
                                 nnbr=nnbr+1
                                 nbra(nnbr)=ia
                              end if
                           else
                              if (ia.ne.iaprec.and.eps(j).eq.0) then
                                 !different from sclbp, I have to 
                                 !consider the object only if the
                                 !midpoint is external O SE C'E' UN 
                                 !PORO E VEDERE IN SEGUITO
                                 iaprec=ia

                                 !assuming any object is not buried
                                 nnbr=nnbr+1
                                 nbra(nnbr)=ia
                              end if
                           end if
                        end do DOKK

                        DOII: do ii=1,nnbr
                           ia=nbra(ii)

                           if (ia.gt.natom) then
                              iii=ia-natom
                              !2011-05-18  Using operations on coord 
                              !and int_coord type variables defined
                              !in module operators_on_coordinates 
                              xq=s123

                              !check if bgp is inside VdW of object
                              call distobj(xq,dist,dxyz,iii,0.0,&
                                          &.false.)

                              !an object can compete with an atom for 
                              !a midpoint only if this midpoint 
                              !is out of the object itself
                              if (dist.ge.0..and.dist.lt.dmn) then
                                 dmn=dist; iac=ia 
                                 dr123=(-dxyz)*(radprb(1)-dist)
                              end if
                           else
                              dx123=s123-xn1(ia); ds2=dx123.dot.dx123 
                              dis=sqrt(ds2)-delphipdb(ia)%rad3

                              !dis= distance to atom surface
                              if (dis.lt.dmn)then
                                 dmn=dis
                                 iac=ia
                              end if
                           end if
                        end do DOII

                        if (iac.eq.0) then
                           if (debug) then
                              write(6,*)"bgp:",i,&
                                  &" might be a cavity point",ix,iy,iz
                              write(6,*)"midpoint",j&
                                         &," in position [A]",s1,s2,s3
                              write(6,*)"it1:",it1," it2:",it2,&
                                                          &" it3:",it3
                           end if

                           ncav=ncav+1

                           !possibly a cavity point
                        else
                           !check to see if it is in the contact 
                           !region of that atom or object by 
                           !projecting it to the atom's acc surface 
                           !and checking against the acc volumes of 
                           !nearby atoms
                           if (iac.le.natom) then
                              dr123=s123-xn1(iac)
                              dsr=sqrt(dr123.dot.dr123)
                              u123=xn1(iac)+((r0(iac)*dr123)/dsr)
                           else
                              u123=s123-dr123
                           end if

                           it123=int((u123-xyzo)*cbai)

                           liml=cbn1(it123%i+1+(lcb+1)*it123%j+&
                                             &(lcb+1)*(mcb+1)*it123%k)
                           limu=cbn2(it123%i+1+(lcb+1)*it123%j+&
                                             &(lcb+1)*(mcb+1)*it123%k)

                           DLIM: do kk=liml,limu
                              ia=cbal(kk)
                              if (ia.le.natom) then
                                 dx123=u123-xn1(ia)
                                 ds2=dx123.dot.dx123
                                 flag=.true.

                                 if (ds2.lt.rs2(ia)) then
                                    flag=.false.; exit DLIM
                                 end if
                              else
                                 if (ia.ne.iac.and.eps(j).eq.0) then
                                    xq=u123

                                    call distobj(xq,dist,dxyz,&
                                           &ia-natom,radprb(1),.true.)

                                    flag=.true.
                                    if (dist.lt.-1.e-6) then
                                       flag=.false.; exit DLIM
                                    end if

                                    !oriented distance from extended 
                                    !object surface if negative => 
                                    !reentrant region
                                 end if
                              end if
                           end do DLIM

                           !it is in the contact region. flag the 
                           !midpoint so it is not checked again
                           !iac is atom number...NOT increased by 1

                           !2011-05-18 To get rid of goto 201 
                           !statements above
                           if (flag) then
                              eps(j)=-iac; eps2(j)=-iac; cycle D200
                           end if
                        end if

                        eps(j)=1 !eps = 1 means cavity or reentrant

                        !remap iepsmp
                        if (iac.eq.0) then 
                           !this is an assumption, still to deeply 
                           !understand meaning of cavity here and to 
                           !improve this choice!WWW
                           if (ia.gt.0) then
                              imedia=iatmmed(ia)
                           else
                              write(6,*)&
                               &'assigning arbitrary epsilon in cavity'
                              imedia=iatmmed(1)
                           end if
                        else
                           imedia=iatmmed(iac)
                        end if
                        
                        select case(imap(4,j))
                        case(1)
                           iepsmp(ix+imap(1,j),iy+imap(2,j),&
                                 &iz+imap(3,j))%i=eps(j)+imedia*epsdim
                        case(2)
                           iepsmp(ix+imap(1,j),iy+imap(2,j),&
                                 &iz+imap(3,j))%j=eps(j)+imedia*epsdim
                        case(3)
                           iepsmp(ix+imap(1,j),iy+imap(2,j),&
                                 &iz+imap(3,j))%k=eps(j)+imedia*epsdim
                        case default
                           write(6,*)'?????'
                        end select

                        eps2(j)=imedia

                        !not assigning the owner but only the medium, 
                        !the former will be assigned in the scale 
                        !routine
        
                        !check to see if the nearest neighbour status 
                        !has been changed..
                        ix2=ix; iy2=iy; iz2=iz 

                        !if the nearest neighbour is a box boundary 
                        !point then skip this since box boundary 
                        !points can not also be dielctric boundary 
                        !points
                        !2011-05-18 Multiple IFs replaced by SELECT 
                        !CASE
                        select case(j)
                        case(1)
                           ix2=ix+1; if(ix2.eq.igrid) cycle D200
                        case(2)
                           iy2=iy+1; if(iy2.eq.igrid) cycle D200
                        case(3)
                           iz2=iz+1; if(iz2.eq.igrid) cycle D200
                        case(4)
                           ix2=ix-1; if(ix2.eq.1) cycle D200
                        case(5)
                           iy2=iy-1; if(iy2.eq.1) cycle D200
                        case(6)
                           iz2=iz-1; if(iz2.eq.1) cycle D200
                        end select

                        !once again one distinguishes between 
                        !internal,external,internal bgp and external 
                        !bgp
                        iext=0
                        ibgp=0

                        !2011-05-18 Changed to i nt_coord derived type
                        itmp(1)=iabs(iepsmp(ix2,iy2,iz2)%i)/epsdim
                        itmp(2)=iabs(iepsmp(ix2,iy2,iz2)%j)/epsdim
                        itmp(3)=iabs(iepsmp(ix2,iy2,iz2)%k)/epsdim
                        itmp(4)=iabs(iepsmp(ix2-1,iy2,iz2)%i)/epsdim
                        itmp(5)=iabs(iepsmp(ix2,iy2-1,iz2)%j)/epsdim
                        itmp(6)=iabs(iepsmp(ix2,iy2,iz2-1)%k)/epsdim

                        if(itmp(1).eq.0) iext=1
                        if(itmp(1).ne.itmp(6)) ibgp=1

                        do cont=2,6
                           if(itmp(cont).eq.0) iext=1
                           if(itmp(cont).ne.itmp(cont-1)) ibgp=1
                        end do

                        if (debug) then
                           nt=0
                           if((iepsmp(ix2,iy2,iz2)%i/epsdim).gt.0)   &
                                                              &nt=nt+1
                           if((iepsmp(ix2,iy2,iz2)%j/epsdim).gt.0)   &
                                                              &nt=nt+1
                           if((iepsmp(ix2,iy2,iz2)%k/epsdim).gt.0)   &
                                                              &nt=nt+1
                           if((iepsmp(ix2-1,iy2,iz2)%i/epsdim).gt.0) &
                                                              &nt=nt+1
                           if((iepsmp(ix2,iy2-1,iz2)%j/epsdim).gt.0) &
                                                              &nt=nt+1
                           if((iepsmp(ix2,iy2,iz2-1)%k/epsdim).gt.0) &
                                                              &nt=nt+1
                           if(nbe(nt).neqv.(ibgp.eq.1.and.iext.eq.1))&
                              & then
                              write(6,*)'PROBLEMS3',ix2,iy2,iz2
                           end if
                        end if !end+debugging
                        
                        if ((ibgp.eq.0).and.&
                                   &(bndeps(ix2,iy2,iz2,1).ne.0)) then
                           !reset bndeps for that point (i.e. remove 
                           !bgp flag).
                           !a bgp become internal
                           ibnumsurf=ibnumsurf-bndeps(ix2,iy2,iz2,2)
                           bndeps(ix2,iy2,iz2,1)=0
                           bndeps(ix2,iy2,iz2,2)=0
                           mr=mr+1
                        else
                           if (ibgp.eq.1.and.iext.eq.0.and.&
                                      &bndeps(ix2,iy2,iz2,2).eq.1)then
                              !an ext  bgp is turned into an internal 
                              !bgp
                              ibnumsurf=ibnumsurf-1
                              bndeps(ix2,iy2,iz2,2)=0
                           end if
                        end if

                        if (ibgp.eq.1.and.&
                                     &bndeps(ix2,iy2,iz2,1).eq.0) then
                           !create a new boundary point..
                           m=m+1
                           bndeps(ix2,iy2,iz2,1)=n2+m
                           ibnd(n2+m)=int_coord(ix2,iy2,iz2)
                           bndeps(ix2,iy2,iz2,2)=iext
                           ibnumsurf=ibnumsurf+bndeps(ix2,iy2,iz2,2)
                        end if
                     end if

                     !now jump to the next midpoint of the same grid 
                     !point
                  end do D200

                  !remap iepsmp in case there have been changes..
                  !(that is some  0's became -1's)
                  !in other words: midpoint must remain external to 
                  !objects
                  do jj=1,6
                     !in this way I can deal with eps(jj)<0
                     isign=1

                     !iord=owner of the midpoint jj before change or 
                     !after eps=1
                     iiord=comp(iepsmp(ix+imap(1,jj),iy+imap(2,jj),&
                                           &iz+imap(3,jj)),imap(4,jj))
                     iord=mod(iiord,epsdim)

                     !the last changed for sure has not iord<0
                     !there can be iord<0 due to nearest neighbors 
                     !already changed
                     if (iord.lt.0) cycle

                     !if it has changed at previous step, dont change 
                     !anymore
                     if (eps(jj).lt.0) then
                        isign=-1
                        if (iord.eq.0) iord=1
                     end if
                     
                     jjj=iabs(iiord)/epsdim

                     select case(imap(4,jj))
                     case(1)
                        iepsmp(ix+imap(1,jj),iy+imap(2,jj),&
                             &iz+imap(3,jj))%i=isign*(iord+jjj*epsdim)
                     case(2)
                        iepsmp(ix+imap(1,jj),iy+imap(2,jj),&
                             &iz+imap(3,jj))%j=isign*(iord+jjj*epsdim)
                     case(3)
                        iepsmp(ix+imap(1,jj),iy+imap(2,jj),&
                             &iz+imap(3,jj))%k=isign*(iord+jjj*epsdim)
                     case default
                         write(6,*)'?????'
                     end select

                     !left iord definition with mod since if it is <> 
                     !0, it keeps its identity
                  end do

                  !at this point one still can trace what changes have 
                  !been made check to see if this is still a boundary 
                  !point
                  !once again one distinguishes between 
                  !internal,external,internal bgp and external bgp
                  iext=0
                  ibgp=0

                  itmp(1)=iabs(iepsmp(ix,iy,iz)%i)/epsdim
                  itmp(2)=iabs(iepsmp(ix,iy,iz)%j)/epsdim
                  itmp(3)=iabs(iepsmp(ix,iy,iz)%k)/epsdim
                  itmp(4)=iabs(iepsmp(ix-1,iy,iz)%i)/epsdim
                  itmp(5)=iabs(iepsmp(ix,iy-1,iz)%j)/epsdim
                  itmp(6)=iabs(iepsmp(ix,iy,iz-1)%k)/epsdim

                  if(itmp(1).eq.0) iext=1
                  if(itmp(1).ne.itmp(6)) ibgp=1

                  do cont=2,6
                     if(itmp(cont).eq.0) iext=1
                     if(itmp(cont).ne.itmp(cont-1)) ibgp=1
                  end do

                  if (debug) then	
                     nt=0
                     if((iabs(iepsmp(ix,iy,iz)%i)/epsdim).gt.0)nt=nt+1
                     if((iabs(iepsmp(ix,iy,iz)%j)/epsdim).gt.0)nt=nt+1
                     if((iabs(iepsmp(ix,iy,iz)%k)/epsdim).gt.0)nt=nt+1
                     if((iabs(iepsmp(ix-1,iy,iz)%i)/epsdim).gt.0) &
                                                              &nt=nt+1
                     if((iabs(iepsmp(ix,iy-1,iz)%j)/epsdim).gt.0) &
                                                              &nt=nt+1
                     if((iabs(iepsmp(ix,iy,iz-1)%k)/epsdim).gt.0) &
                                                              &nt=nt+1

                     if (nbe(nt).neqv.(ibgp.eq.1.and.iext.eq.1)) then
                        write(6,*)'PROBLEMS4',ix,iy,iz
                        write(6,*)&
                           &'epsdim=',epsdim,'ibgp=',ibgp,'iext=',iext
                        write(6,*)'itmp',itmp
                        write(6,*)iepsmp(ix,iy,iz)%i
                        write(6,*)iepsmp(ix,iy,iz)%j
                        write(6,*)iepsmp(ix,iy,iz)%k
                        write(6,*)iepsmp(ix-1,iy,iz)%i
                        write(6,*)iepsmp(ix,iy-1,iz)%j
                        write(6,*)iepsmp(ix,iy,iz-1)%k
                     end if
                  end if

                  !if not now a boundary element change bndeps
                  if ((iext.eq.0).or.(ibgp.eq.0)) then
                     ibnumsurf=ibnumsurf-bndeps(ix,iy,iz,2)
                     if(ibgp.eq.1) bndeps(ix,iy,iz,2)=iext
                     
                     if (ibgp.eq.0) then
                        bndeps(ix,iy,iz,1)=0
                        bndeps(ix,iy,iz,2)=0
                        mr=mr+1
                        if(iext.eq.1)&
                          &write(6,*)'!!!born a new external point!!!'
                     end if
                  end if

               end if !if end for whether bndeps is nonzero
            end do !next boundary point FINISH

            n1=n2+1
            n2=n2+m
            if(verbose) write(6,*)'bgp added m=',m,&
                                              &' bgp removed  mr =',mr

            if (m.gt.mpr) then
               ndv=ndv+1
               if (ndv.gt.20) then !Lin Li: the value used to be 2,
                                   !        sometimes not enough
                  write(6,*)'surface iteration did not converge'
                  stop
               end if
            else
               ndv=0
            end if

            !2011-05-18 Replaced goto 100 statement
            if(m.le.0) exit D100
         end do D100

         if (n2.gt.ibmx) then
            write(6,*)'ibnd upper bound ',n2,' exceeds ibmx'
            stop
         endif

         !2011-05-18 Memory cleanup is done with DEALLOCATE F95 
         !statement
         if(allocated(cbn1)) deallocate(cbn1)
         if(allocated(cbn2)) deallocate(cbn2)
         if(allocated(cbal)) deallocate(cbal)

         if (verbose) then
            write(6,*)&
              &'no. cavity mid-points inaccessible to solvent = ',ncav
         end if

         !consolidate the list, removing dead boundary points, adding 
         !new ones..
         j=0
         ncms=0

         !2011-05-19  Array is re-sized keeping old values in the 
         !memory. 
         if (allocated(ibgrd)) then
            Nar=size(ibgrd)
            if (Nar.lt.ibmx) then
               allocate(ibgrd_temp(Nar)) ; ibgrd_temp=ibgrd
               deallocate(ibgrd); allocate(ibgrd(ibmx))
               ibgrd(1:Nar)=ibgrd_temp; deallocate(ibgrd_temp)
            end if
         else
            allocate(ibgrd(ibmx))
         end if

         do i=1,n2
            ixyz=ibnd(i)
            ix=ixyz%i; iy=ixyz%j; iz=ixyz%k

            if (bndeps(ix,iy,iz,1).ne.0) then
               j=j+1
               bndeps(ix,iy,iz,1)=j

               !2011-05-19 Precaution not to exceed array size (see 
               !above comment)
               if (j.le.ibmx) then
                  ibgrd(j)=ixyz
               else
                  write(6,*) 'j=',j,' is larger than ibmx= ',ibmx,&
                             ', thus stopped...'
                  stop
               end if
            end if

            do jj=1,6
               !2011-05-17  Using operations on coord and int_coord 
               !type variables defined in module 
               !operators_on_coordinates 
               ntt=comp(iepsmp(ix+imap(1,jj),iy+imap(2,jj),&
                                           &iz+imap(3,jj)),imap(4,jj))
               nt=mod(ntt,epsdim)

               if (nt.lt.0) then
                  ntt=-ntt

                  select case(imap(4,jj))
                  case(1)
                     iepsmp(ix+imap(1,jj),iy+imap(2,jj),&
                                                 &iz+imap(3,jj))%i=ntt
                  case(2)
                     iepsmp(ix+imap(1,jj),iy+imap(2,jj),&
                                                 &iz+imap(3,jj))%j=ntt
                  case(3)
                     iepsmp(ix+imap(1,jj),iy+imap(2,jj),&
                                                 &iz+imap(3,jj))%k=ntt
                  case default
                     write(6,*)'?????'
                  end select

                  if (nt.eq.-1) then
                     ntt=ntt-1

                     select case(imap(4,jj))
                     case(1); iepsmp(ix+imap(1,jj),iy+imap(2,jj),&
                                                 &iz+imap(3,jj))%i=ntt
                     case(2); iepsmp(ix+imap(1,jj),iy+imap(2,jj),&
                                                 &iz+imap(3,jj))%j=ntt
                     case(3); iepsmp(ix+imap(1,jj),iy+imap(2,jj),&
                                                 &iz+imap(3,jj))%k=ntt
                     case default;  write(6,*)'?????'
                     end select
                  end if
               end if

               if (debug) then
                  ntt=comp(iepsmp(ix+imap(1,jj),iy+imap(2,jj),&
                                           &iz+imap(3,jj)),imap(4,jj))
                  nt=mod(ntt,epsdim)
                  jjj=ntt/epsdim
                  if (nt.eq.0.and.jjj.gt.0) then
                     write(6,*)'PROBLEMS 5',ix,iy,iz,jj
                  end if
               end if
            end do
         end do

         if (j.gt.ibmx)then
            write(6,*) 'no. ms points exceeds ibmx'; stop
         end if
         
         ibnum=j
      
         if (verbose) then 
            write(6,*)'after surface elaboration ibnum= ',ibnum
            write(6,*)'    and               ibnumsurf= ',ibnumsurf
         end if

      end if

      if (debug) then
         open(52,file='bgpDurante.log',form='formatted')
         write(52,*)igrid,ibnum
         do iiii=1,ibnum
            write(52,*)iiii,ibnd(iiii)
            write(52,*)iiii,ibnd(3*iiii-2),ibnd(3*iiii-1),ibnd(3*iiii)
         end do
         close (52)
      endif

      if(allocated(bndeps)) deallocate(bndeps)
      if(allocated(ibnd)) deallocate(ibnd)

      !scale bondary grid point positions relative to acc data
      if (isolv.and.(irea.or.logs.or.lognl.or.isen.or.isch)) then
         if(verbose)write(6,*)'scaling boundary grid points'

         !2011-05-19 Arrays allocated by ordinary F95 allocate 
         !statement
         allocate(scspos(ibnum))

         do j=1,ibnum
            !2011-05-19  Using operations on coord and int_coord type 
            !variables defined in module operators_on_coordinates 
            scspos(j)=float(ibgrd(j))
         end do

         allocate(scsnor(ibnum),atsurf(ibnum),atndx(ibnum))

         !2011-05-19 Other parameters to subroutine are transfered via 
         !module architecture (declared in pointers and qlog modules)
         call sclbp(scspos,scsnor,ibnum,iab1,iab2)

         if (verbose) then
            write(6,*) iall,&
                   & ' points had to be assigned by global comparison'
         end if
 
         if (.not.isite.and.allocated(scsnor)) deallocate(scsnor)
      end if

      if (isrf.and..not.ivertnorm)then
         if (ionlymol) then
            !2011-05-19 Parameters transfered via module architecture
            allocate(egrid(igrid,igrid,igrid))

            call msrf

            if(allocated(egrid) )deallocate(egrid)
         else
            write(6,*)'msrf routine cannot be run'
            write(6,*)'because there are also geometric objects'
         end if
      end if

      if (.not.isitsf.and..not.isite.and.&
                                   &.not.(isch.and.scrgfrm.ne.0)) then
         if(allocated(atndx)) deallocate(atndx)
         if(allocated(atsurf)) deallocate(atsurf)
      end if

      if(allocated(iab1))  deallocate(iab1)
      if(allocated(iab2))  deallocate(iab2)
      if(allocated(icume)) deallocate(icume)
      if(allocated(r0))    deallocate(r0)
      if(allocated(r02))   deallocate(r02)
      if(allocated(rs2))   deallocate(rs2)
      if(allocated(ast))   deallocate(ast)

      end subroutine vwtms
