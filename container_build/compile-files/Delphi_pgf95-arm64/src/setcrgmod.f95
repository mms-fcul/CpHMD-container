      ! 2011-06-02 Created by PK from files mkdbsf.f, setfcrg.f,
      module setcrgmod
      use qlog
      use pointers
      use operators_on_coordinates
      implicit none
      contains
      !============================================================================>
      include 'inc_setcrgmod_mkdbsf.f95'        ! subroutine crgarr
      !============================================================================>
      include 'inc_setcrgmod_setfcrg.f95'       ! subroutine distrTOpoint2
      !============================================================================>
      include 'inc_setcrgmod_setcrg.f95'        ! subroutines basisortho
      !============================================================================>
      include 'inc_setcrgmod_wrteps.f95'        ! subroutines basisortho
      !============================================================================>
      end module setcrgmod
