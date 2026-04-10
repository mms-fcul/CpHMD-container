!#####################################################################
!ATTENTION! This file is part of epsmakmod module.
!           Do not compile separately!
!
!program to create the accessible surface arcs object-object
!Walter March 3, 2000
!
!2011-05-27 Recursive subroutines are supported by default in F95
!2011-05-27 Removed redundant declarations already made in qlog and
!           pointers modules and changed some arrays to coord type 
!           variables
!#####################################################################
      recursive subroutine objvertices(nacc,xmin,xmax,side,h,lookmol)
      
      implicit none

      type(int_coord) :: etrm,ixyz
      integer :: objecttype,itmp,nt,numsup
      integer :: epsdim,kind
      logical :: lookmol,intersdentro,inacqua
      character(96) :: strtmp
      real :: side,h
      type(coord) :: dxyz,xp,xyz
      real :: ro,sro,dx,dy,dz,tmp,dist
      type(coord) :: cubcrd,vmin,vmax,xmin,xmax
      type(grid_value) :: distv(3)
      real :: coeff
      real :: beta,delta,dist1,dist2
      !2011-05-27 Declaration added due to IMPLICIT NONE
      integer :: nacc,ii,i,ix,iy,iz,j,k
      real :: goff,dot,x,y,z

      goff=(igrid+1.)/2.
      epsdim=natom+nobject+2

      !2011-05-27 Using operations on coord and int_coord type 
      !variables defined in module operators_on_coordinates 
      etrm=1+int(0.5+((xmax-xmin)/side))

      !now look into the little cubes
      !oss. xmin = firts centers vector
      x=xmin%x
      D030: do i=1,etrm%i
         y=xmin%y
         D020: do j=1,etrm%j
            z=xmin%z
            D010: do k=1,etrm%k                 
               !now (x,y,z) is the little cube center
               if (side.le.h.and.lookmol) then
                  !in order to see if this cube is inside a molecule
                  !find the closest midpoint (here vmin/vmax used by 
                  !chance)
                  lookmol=.false.
                  vmin=coord(x,y,z)

                  !2011-05-27 Subroutine ctog is too small, thus 
                  !transfered its body here
                  
                  !Begin ctog body
                  vmax=(vmin-oldmid)*scale + goff
                  !End ctog body

                  ixyz=int(vmax+0.5)
     
                  !now vmax is in the h-sized cube system
                  dxyz=vmax-float(ixyz)

                  sro=sqrt(dxyz%x*dxyz%x+dxyz%y*dxyz%y)
                  ro=sqrt(dxyz.dot.dxyz)

                  !reasonably approximated algorithm, being the exact 
                  !one too expensive
                  if (sro.eq.0.) then
                     if (dxyz%z.ge.0.) then
                         nt=iepsmp(ix,iy,iz)%k
                     else
                         nt=iepsmp(ix,iy,iz-1)%k
                     end if
                  else
                     if (abs(dxyz%z/ro).lt.0.5) then
                        if(dxyz%x/sro.lt.-0.707) nt=iepsmp(ix-1,iy,iz)%i
                        if(dxyz%y/sro.lt.-0.707) nt=iepsmp(ix,iy-1,iz)%j
                        if(dxyz%x/sro.ge.0.707) nt=iepsmp(ix,iy,iz)%i
                        if(dxyz%y/sro.ge.0.707) nt=iepsmp(ix,iy,iz)%j
                     else
                        nt=iepsmp(ix,iy,iz-1)%k
                        if(dxyz%z.gt.0.) nt=iepsmp(ix,iy,iz)%k
                     end if
                  end if
 
                  nt=mod(nt,epsdim)

                  !modification pores: here it is expecting to find 
                  !atoms but it might see also pores?
                  if (nt.le.natom+1.and.nt.gt.0) cycle D010
               end if

               !calculate if the object surf close to the cube are more 
               !or less than 2
               numsup=0

               !for objects, only water probes involved
               tmp=radprb(1)+0.5*side
               xyz=coord(x,y,z)
               D300: do  ii=1,nobject
                  !find first three objects intersecting circumscribed 
                  !sphere to the current cube
                  strtmp=dataobject(ii,1)
                  read(strtmp(16:18),*)kind
                  if (kind.eq.3) write(6,*)'PROBLEMSobj'

                  if (strtmp(1:4).ne.'is a'.and.kind.ne.2) then
                     !2011-05-27 Using operations on coord and 
                     !int_coord type variables defined in module 
                     !operators_on_coordinates 
                     if((xyz.vorlt.limobject(ii)%min).or.&
                                  & (xyz.vorgt.limobject(ii)%max).or.&
                                  & numsup.gt.2 ) cycle D300 
                     xp=xyz

                     call distobj(xp,dist,dxyz,ii,radprb(1),.false.)

                     if (abs(dist).lt.0.86603*side) then
                        !look if it is within the circumscribed sphere
                        numsup=numsup+1
                        distv(numsup)%xyz=dxyz; distv(numsup)%value=dist
                     end if
                  end if
               end do D300

               D400: do
                  if((numsup.gt.2.and.side.gt.sidemin).or.(numsup.eq.2 &
                     &  .and.side.gt.sideinter)) then
                     !if next side is one fourth of the current one, 
                     !then  1/2-1/8
                     tmp=0.25*side
                     vmin=xyz-tmp; vmax=xyz+tmp

                     call objvertices(nacc,vmin,vmax,.5*side,h,lookmol)
                  else  
                     if (numsup.eq.2.and.side.le.sideinter) then
                        !calculate the closest intersection between the 
                        !two surfaces
                        dot=distv(1)%xyz.dot.distv(2)%xyz
                        tmp=1.-dot*dot
                        if(tmp.le.tol) exit D400
                        coeff=1./tmp
                        dist1=distv(1)%value; dist2=distv(2)%value
                        delta=(dist2*dot-dist1)*coeff
                        beta=(dist2-dist1*dot)*coeff
                        intersdentro=(4*(delta**2+&
                              &beta**2-2*delta*beta*dot).le.3*side**2)
                        if (.not.intersdentro) exit D400

                        !inacqua means out of the  SAS but I ignore it
                        inacqua=(dist1.ge.0..and.dist2.ge.0.)

                        if (side.gt.sidemin) cycle D400

                        dxyz=(delta*distv(1)%xyz)-(beta*distv(2)%xyz)

                        !check whether it is a false positive vertex
                        tmp=.5*side+tol
                        if ((abs(dxyz).vandle.tmp).and.&
                                       &(abs(delta*beta).gt.tol)) then

                           !if it is inside, then let's write it!
                           nacc=nacc+1
                           expos(nacc)=dxyz+xyz
                        end if
                     end if
                  end if
               end do D400

               !in any other case skip this cube
               !NOTE: I am skipping also in case there are 3 surfaces 
               !and side< sidemin!!!
               z=z+side
            end do D010

            y=y+side
         end do D020

         x=x+side
      end do D030

      end subroutine objvertices
