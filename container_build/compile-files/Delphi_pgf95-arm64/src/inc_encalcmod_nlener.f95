!#####################################################################
!ATTENTION!  This file is part of encalcmod module.
!            Do not compile separately!
!#####################################################################
      
      subroutine nlener(ergnl,igridout)

      !2011-06-11 Declarations added due to IMPLICIT NONE
      double precision :: ergnl,ergsolv,ergosm
      integer :: izero,i,j,k,n,saved
      type(coord) :: txyz,gxyz
      real :: dist
      integer :: z1p,z1m,z2p,z2m,igridout
      real :: carica,espp1,espp2,espm1,espm2,cs1,cs2,p1,p2,m1,m2
      real :: carica2,cutedgesi,cutedgesj,cutedgesk,goff
      real :: dhi1,dhi2,dhi3,dhi4,dhi5,tmpcar,phi,c

      !Secondly, scan all the grid points and calculate ergnl
      saved=0 ; ergsolv=0.0; ergosm=0.0;  n=0
      goff = (igrid + 1.)/2.
      gxyz=(-goff/scale) + oldmid
      c=scale*scale*scale
      dhi5=-chi5/6.; dhi4=-chi4/5.; dhi3=-chi3/4.
      dhi2=-chi2/3.; dhi1=-chi1/2.
      if (logions)allocate(sout(igrid*igrid*igrid))

      do k=1,igrid
         cutedgesk=1.
         if (k.eq.1.or.k.eq.igrid) cutedgesk=.5
         do j=1,igrid
            cutedgesj=cutedgesk
            if (j.eq.1.or.j.eq.igrid) cutedgesj=cutedgesk*.5
            do i=1,igrid
               cutedgesi=cutedgesj
               if (i.eq.1.or.i.eq.igrid) cutedgesi=cutedgesj*.5
 
               if (idebmap(i,j,k)) then 
                  phi=phimap(i,j,k)

                  !esp..=exp(..)-1
                  !p1=-z1p*phi
                  !espp1=sinh(p1)+2.*sinh(p1/2)**2
                  !m1=z1m*phi
                  !espm1=sinh(m1)+2.*sinh(m1/2)**2
                  !
                  !carica2=0.0
                  !if (cs2.gt.0.) then
                  !   p2=-z2p*phi
                  !   espp2=sinh(p2)+2.*sinh(p2/2)**2
                  !   m2=z2m*phi
                  !   espm2=sinh(m2)+2.*sinh(m2/2)**2
                  !   ergosm=ergosm-cs2*(z2m*espp2+z2p*espm2)
                  !   carica2=-2.*cs2*z2p*sinh(sum2*phi)*exp(diff2*phi)
                  !end if
                  !
                  !add the osmotic pressure term ergosm in kt units
                  !ergosm=ergosm-cs1*(z1m*espp1+z1p*espm1)
                  !ergosmb=ergosmb-&
                  !         &2.*cs1*exp(-.5*z1p*phi)*sinh(-.5*z1p*phi)
                  !if (ergosmb.ne.ergosm) &
                  !                   &write(6,*)'azzz',ergosmb,ergosm
                  !carica=carica2-&
                  !          &2.*cs1*z1p*sinh(sum1*phi)*exp(diff1*phi)

                  !Horner scheme for charge and osmotic term
                  tmpcar=chi5  *phi+chi4
                  tmpcar=tmpcar*phi+chi3
                  tmpcar=tmpcar*phi+chi2
                  tmpcar=tmpcar*phi+chi1
                  carica=cutedgesi*tmpcar*phi/c

                  tmpcar=dhi5  *phi+dhi4
                  tmpcar=tmpcar*phi+dhi3
                  tmpcar=tmpcar*phi+dhi2
                  tmpcar=tmpcar*phi+dhi1
                  ergosm=ergosm+cutedgesi*tmpcar*phi*phi      
 
                  if (logions) then
                     !if the gp is in solution and the contribution is 
                     !higher than a threshold, then put this 
                     !information in a list
                     n=n+1
                     sout(n)%xyz=float(int_coord(i,j,k))+gxyz
                     sout(n)%value=carica !lilin
                  end if
  
                  ergsolv=ergsolv-carica*phi
  
               end if
            end do
         end do
      end do
      
      ergsolv=ergsolv*.5*.0006023
      ergosm=-ergosm*.0006023/c

      igridout=n

      !2011-06-11 Array resizing seems to be unnecessary (array is 
      !           deallocated soon in the subroutine clbnonl
      !if (logions) then
      !   i_sout=memalloc(i_sout,realsiz,4*igridout)
      !end if


      if(ideveloper) then
         write(6,*)'rho*phi/2 term in solution      :',-ergsolv,'kt'
         write(6,*)'osmotic pressure term           :',-ergosm,'kt'
      else
         write(6,"(a,f20.4,a)")' rho*phi/2 term in solution      :',-ergsolv,' kt'
         write(6,"(a,f20.4,a)")' osmotic pressure term           :',-ergosm,' kt'
      end if

      ergnl=ergnl+ergosm+ergsolv
      
      end subroutine nlener
