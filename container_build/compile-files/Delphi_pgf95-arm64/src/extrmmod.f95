      module extrmmod
      use pointers
      use qlog
      use operators_on_coordinates
      use misc
      implicit none
      character(96) :: strtmp
      
      contains
      
      !============================================================================>
      include 'inc_extrmmod_extrmobjects.f95'   ! subroutine extrmobjects
      !============================================================================>
      include 'inc_extrmmod_extrm.f95'          ! subroutines extrm,off,warning
      !============================================================================>
      include 'inc_extrmmod_wrtprm.f95'         ! subroutines wrtprm,wrt
      !============================================================================>
      end module extrmmod
