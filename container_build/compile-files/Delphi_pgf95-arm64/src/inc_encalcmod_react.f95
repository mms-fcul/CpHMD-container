!#####################################################################
!ATTENTION! This file is part of encalcmod module.
!           Do not compile separately!
!
!2011-06-11 All other parameters are accessible via qlog and pointers
!           modules
!2011-06-11 Arrays below are declared in pointers module
!           real phimap(igrid,igrid,igrid)
!           dimension sqs(1),ibgrd(3,ibnum),cqs(ncrgmx)
!2011-06-11 Array mv seems is not used
!           dimension mv(6)
!           dimension cgrid(1),xo(3),xn(3),sen(1)
!           dimension spdiv(1),gchrg(icount1b)
!           dimension spot(1),scspos(3,ibnum)
!           integer gchrgp(3,icount1b)
!gchrgp= (i,j,k)  for the icount1b charged grid points
!#####################################################################

      subroutine react(ergs,ergas,iisitpot)

      real :: r2,Qd(5)
      logical ido
      type(coord) :: xo,xn,pxyz
      real :: radius,ergsrea,cost
      character(80) :: line
      integer :: jj,qq
      real :: dista,tmp(ibnum*iisitpot),tmp1
      type(int_extrema) :: lim
      type(int_coord) :: ixyz
      !Declarations added due to IMPLICIT NONE
      real*8 :: ergs,ergas,en,en1,ent
      real :: epsdim,dist,dist4,en2,en3,fact,prod,ptemp
      real :: schrgj,temp,temp1,temp2,temp3,spt1,sixth
      integer :: iisitpot,ii,i,j,ix,iy,iz,imed,ibnum

      epsdim=nobject+natom+2; ibnum=icount2b

      allocate(sitephi(5,npotenziali+nvincfit),schrg(ibnum))
      allocate(sqs(ibnum),spdiv(ibnum),sen(ibnum),spot(ibnum))
      allocate(cqs(nqass))

      !fact=six divide by 4 pi, multiplied by epkt, scale
      !deleted epsval in order to have usual srf charge definition
      fact=0.9549296586/(2.0*scale*epkt)

      sixth=1.0/6.0 ;  en1=0.0;  ergsrea=0.0

      !calculate surface charge

      !goff = (igrid + 1.)/2.
      !eps2=1.
      !eps1=80
      !co=1.*561/eps2
      !d=3
      !epsdim=5
      !
      !write(6,*)'faccio..',co
      !ix=goff
      !iy=goff
      !xx=(ix-goff)/scale +oldmid(1)
      !yy=(iy-goff)/scale +oldmid(2)
      !
      !do 9000 iz=1,igrid
      !   zz=(iz-goff)/scale +oldmid(3)
      !   dista=sqrt(xx**2+yy**2+(zz+d)**2)
      !   if (zz.lt.0) then
      !      phi=1./dista
      !      phi=phi+(eps2-eps1)/((eps1+eps2)*sqrt(xx**2+yy**2+
      !          &(zz-d)**2))
      !   else
      !      phi=2.*eps2/((eps1+eps2)*dista)
      !   end if
      !   write(6,*)'point:',iz,zz
      !   write(6,*)'calculated=',phimap(ix,iy,iz)
      !   write(6,*)'analytical:',co*phi
      !   write(6,*)''
      !9000 continue

      do i=1,ibnum
         ix=ibgrd(i)%i; iy=ibgrd(i)%j;  iz=ibgrd(i)%k
         temp1=phimap(ix+1,iy,iz)+phimap(ix-1,iy,iz)
         temp2=phimap(ix,iy+1,iz)+phimap(ix,iy-1,iz)
         temp3=phimap(ix,iy,iz+1)+phimap(ix,iy,iz-1)
         spdiv(i)=phimap(ix,iy,iz)-(temp1+temp2+temp3)*sixth
      end do

      if (ibc.ne.0)then 
         !this part removes fixed charge assigned to boundary 
         !2011-06-11 Subroutine rdiv is small, thus moved its body here
         allocate(cgrid(igrid,igrid,igrid))

         !--------------- Begin body of rdiv ----------------
         do i=1,ibnum
            cgrid(ibgrd(i)%i,ibgrd(i)%j,ibgrd(i)%k)=spdiv(i)
         end do

         do i=1,ibc
            ix=cgbp(i)%xyz%x;iy=cgbp(i)%xyz%y; iz=cgbp(i)%xyz%z
            cgrid(ix,iy,iz)=cgrid(ix,iy,iz)-cgbp(i)%value1

            !NO LONGERthe self reaction energy close to surface is not 
            !counted twice
         end do

         do i=1,ibnum
            spdiv(i)=cgrid(ibgrd(i)%i,ibgrd(i)%j,ibgrd(i)%k)
         end do
         !-------------- End body of rdiv ------------------

         deallocate(cgrid)
      end if

      !spdiv is equal to the 'charge' as would appear on the grid
      !to replace the boundary elements
      en1=0.0
      do i=1,ibnum
         schrg(i)=spdiv(i)*fact; en1=en1+(schrg(i))
      end do

      !attempt to spread charge about a bit..
      !ispread=.true.

      !generate pseudo distances for surface points
      if (irea.or.logs.or.lognl.or.isen.or.isch) then
         do i=1,ibnum
            sqs(i)=0.5*(scspos(i).dot.scspos(i))
         end do

         !this is not formally correct, in principle isitpot 
         !doesnt require dipole calculation
         if (isitpot) then
            if (nvincfit.ge.3) then
               !dipole calculation
               pxyz=coord(0.0,0.0,0.0)
               tmp1=abs(medeps(0)*epkt-1.0)

               do i=1,ibnum
                  !sto assumendo che non ci siano oggetti di mezzo se
                  !faccio questo calcolo
                  tmp(i)=1.0
                  if (atndx(i).ne.-1.and.tmp1.gt.tol) then
                     imed=iatmmed(atsurf(i))
                     tmp(i)=1.0+medeps(imed)*tmp1/(medeps(imed)-&
                            &medeps(0))
                  end if
	
                  !non sono sicuro che schrg(j) sia proprio la carica
                  pxyz=pxyz+((schrg(i)*tmp(i))*scspos(i))
               end do

               sitephi(4,npotenziali+1)=pxyz%x
               sitephi(4,npotenziali+2)=pxyz%y
               sitephi(4,npotenziali+3)=pxyz%z
            end if

            if (nvincfit.ge.8) then
               !quadrupole calculation
               Qd=0.

               do i=1,ibnum
                  r2=2.0*sqs(i)
                  Qd(1)=Qd(1)+schrg(i)*(3.0D0*scspos(i)%x**2-r2)*tmp(i)
                  Qd(2)=Qd(2)+schrg(i)*(3.0D0*scspos(i)%y**2-r2)*tmp(i)
                  Qd(3)=Qd(3)+&
                        &schrg(i)*3.0D0*scspos(i)%x*scspos(i)%y*tmp(i)
                  Qd(4)=Qd(4)+&
                        &schrg(i)*3.0D0*scspos(i)%x*scspos(i)%z*tmp(i)
                  Qd(5)=Qd(5)+&
                        &schrg(i)*3.0D0*scspos(i)%y*scspos(i)%z*tmp(i)
               end do

               sitephi(4,npotenziali+4)=Qd(1)
               sitephi(4,npotenziali+5)=Qd(2)
               sitephi(4,npotenziali+6)=Qd(3)
               sitephi(4,npotenziali+7)=Qd(4)
               sitephi(4,npotenziali+8)=Qd(5)
            end if
         end if

         if (logs.or.lognl) then
            !calculate back energies
            en=0.0 ; en2=0.0; en3=0.0; ent=0.0

            !calculating interaction energy between each real charge 
            !and the surrounding polarization charges
            if (nmedia.gt.0) then          
               do i=1,nqass
                  radius=radpolext
                  ii=crgatn(i)
                  if(ii.gt.0.and.ii.le.natom) radius=delphipdb(ii)%rad3
                  cost=(1./(epkt*atmeps(i))-1.) 
                  if (atmeps(i).le.0.) write(6,*)'atmeps error',i,ii

                  if (radius.le.0.) then
                     write(6,*)'charged atom number',ii,&
                             &'radius changed from zero to ',radpolext
                     write(6,*)'BE CAREFUL!! REACTION FIELD ENERGY &
                                         & MIGHT BE OVERESTIMATED!!!!'
                     radius=radpolext
                  end if
          
                  ergsrea=ergsrea+0.5*atmcrg(i)%value*&
                                         & atmcrg(i)%value*cost/radius
               end do

               ergsrea=ergsrea*epkt

               if(ideveloper) then
                  write(6,*) 'self-reaction field energy :    '&
                          &,ergsrea,' kt'
               else
                  write(6,"(a,f20.4,a)") ' self-reaction field energy      :'&
                          &,ergsrea,' kt'
               end if

               if(inrgwrt) write(42,'(a30,f10.4,a4)') 'self-reaction&
                                       &  field energy:',ergsrea,' kt'
            end if

            if (ibufz) then
               do i=1,nqass
                  cqs(i)=(chgpos(i).dot.chgpos(i))/2.
               end do

               en=0.0

               lim%min=2+bufz%min ; lim%max=igrid-1-bufz%max 

               do j=1,ibnum
                  ixyz=ibgrd(i)
                  ido=.true.
                  if((ixyz.vorlt.lim%min).or.(ixyz.vorgt.lim%max)) &
                                                         & ido=.false.

                  if (ido) then
                     dist4=sqs(j); schrgj=schrg(j)/sqrt(2.0)
                     en1=0.0
                     ent=ent+schrg(j)
                        
                     do i=1,nqass
                        prod=scspos(j).dot.chgpos(i)
                        dist=dist4+cqs(i)-prod
                        temp=atmcrg(i)%value/sqrt(dist)
                        en1=en1+temp
                     end do

                     en=en+en1*schrgj; sen(j)=schrgj*en1
                  end if
               end do

               !deleted epsval
               ergs=en*epkt/2.0
                  
               write(6,'(f10.4)')'tot s.charge, no epsin carrying :',ent
               write(6,*) 'corrected reaction field energy:      ',ergs, ' kt'
               write(6,*)'total reaction field energy    :  ',ergsrea+ergs, 'kt'
                  
               if (inrgwrt) then
                  write(42,'(a30,f10.4,a4)') 'corrected reaction field energy:',ergs,' kt'
                  write(42,'(a30,f10.4,a4)') 'total reaction field energy:',ergs+ergsrea,' kt'
               end if
            else !if(.not.ibufz) then
               do i=1,ibnum
                  ptemp=0.0
                  do j=1,nqass
                     xo=scspos(i)-chgpos(j); dist=sqrt(xo.dot.xo)
                     ptemp=ptemp+atmcrg(j)%value/dist
                  end do

                  spot(i)=ptemp
               end do
 
               en=0.0; ent=0.0
                  
               do i=1,ibnum
                  temp=spot(i)*schrg(i); ent=ent+schrg(i); en=en+temp
               end do
 
               !deleted epsval
               ergs=en*epkt/2.0

               write(6,'(a,f20.4)')" total s.charge,no epsin carrying:",ent
               if (ibnum.eq.0.and.igrid.le.5) then
                  write(*,*),"Midpoints are out side the cube and delphi cannot determine the molecular surface."
                  write(*,*),"Please enlarge the gsize or decrease the perfil value."
               end if
                  
               if(ideveloper) then
                  write(6,*)'corrected reaction field energy: ',ergs," kt"
                  write(6,*)'total reaction field energy :   ',ergsrea+ergs,' kt'
               else
                  write(6,"(a,f20.4,a)")' corrected reaction field energy :',ergs," kt"
                  write(6,"(a,f20.4,a)")' total reaction field energy     :',ergsrea+ergs,' kt'
               end if


               if (inrgwrt) then
                  write(42,'(a35,f10.4,a4)') 'corrected reaction field energy:',ergs,' kt'
                  write(42,'(a30,f10.4,a4)') 'total reaction field energy:',ergs+ergsrea,' kt'
               end if
            end if
               
         end if !end of solvation calculation..

         if (isch) then
            write(6,'(2a28)')"writing surface charge file:",&
                                &trim(scrgnam)
            if(scrgfrm.eq.0) open(41,file=trim(scrgnam))

            if ((scrgfrm.eq.1).or.(scrgfrm.eq.2)) then 
               open(41,file=trim(scrgnam))
               write(41,'(a17)') "DELPHI FORMAT PDB"
               write(41,'(a15)') "FORMAT NUMBER=2"
               write(41,*)"      bgp# atom SC   res#      pos","                    scrg          surf energy"
            end if

            !adjusted for spread charges
            if (ibufz) then
               lim%min=2+bufz%min ; lim%max=igrid-1-bufz%max 
            end if

            D450: do i=1,ibnum
               ixyz=ibgrd(i)

               if (ibufz) then
                  ido=.true.
                  if((ixyz.vorlt.lim%min).or.(ixyz.vorgt.lim%max)) &
                                                          &ido=.false.

                  if(.not.ido) cycle D450
               end if
                  
               xo=scspos(i); en1=schrg(i)
                    
               !deleted epsval
               spt1=spot(i)*en1*epkt/2.0

               if(scrgfrm.eq.0)write(41,'(i5,3f8.3,f8.3)') i,xo,en1

               if ((scrgfrm.eq.1).or.(scrgfrm.eq.2)) then
                  jj=atsurf(i)

                  !qq= residue sequence number
                  read(delphipdb(jj)%atinf(12:15),'(I4)')qq

                  !2011-06-11 Subroutine watpte is small, thus moved 
                  !its body here
                  !Begin body of watpte.f
                  line=" "
                  line(1:6)="ATOM  "
                  line(18:20)="SC "

                  write(line(7:11),'(i5)') i
                  write(line(23:26),'(i4)') qq
                  write(line(31:54),'(3f8.3)')xo
                  write(line(56:67),'(g12.5)') en1
                  write(line(69:80),'(g12.5)') spt1
                  !End body of watpte.f

                  write(line(12:16),'(I5)')jj
                  write(41,'(a)') trim(line)
               end if
            end do D450

            close(41)
         end if

         !if isen, calculate surface energy positions and values
         if (isen) then
            !sen will now not be correct since schrg was changed
            write(6,'(a40)') "writing surface energy file: surfen.dat"
               
            open(41,file="surfen.dat")
            en2=0.0
            do i=1,ibnum
               xo=scspos(i)
               en1=sen(i)/2.0
               write(41,'(i5,3f8.3,f8.3)') i,xo,en1
            end do
            close(41)
         end if

      end if
         
      deallocate(sqs,spdiv,sen,spot,cqs)

      end subroutine react
