!#####################################################################
!2011-05-18 Piece of repeating code included in the
!           vwtmc subroutine via INCLUDE statement
!           DO NOT COMPILE SEPARATELY!
!2011-05-18 Array changed from 1D to 3D
!#####################################################################
      liml=iab1(jx+1+(lcb1+1)*jy+(lcb1+1)*(mcb1+1)*jz)
      limu=iab2(jx+1+(lcb1+1)*jy+(lcb1+1)*(mcb1+1)*jz)

      do ii=liml,limu
         iarv= icume(ii)

         !2011-05-17  Using operations on coord and int_coord type 
         !variables defined in module operators_on_coordinates 
         dist=(s123-expos(iarv)).dot.(s123-expos(iarv))

         if (dist.lt.prbrd2) then
            eps(j)=-1
            eps2(j)=-1
            cycle D200
         else if (dist.lt.dmn) then
            iacv=iarv
            dmn=dist
         end if
      end do
