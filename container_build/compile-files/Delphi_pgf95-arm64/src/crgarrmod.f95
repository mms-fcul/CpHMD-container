!#####################################################################
!2011-05-30 Created by PK from subroutine crgarr and subroutines within.
!#####################################################################
      module crgarrmod
      
      use qlog
      use pointers
      use operators_on_coordinates
      implicit none
      
      contains
!====================================================================>
      include 'inc_crgarrmod_crgarr.f95'       !subroutine crgarr
!====================================================================>
      include 'inc_crgarrmod_distrTOpoint.f95' !subroutine distrTOpoint2
!====================================================================>
      include 'inc_crgarrmod_basisortho.f95'   !subroutines basisortho
                                               !and swap
!====================================================================>
      end module crgarrmod
