      ! 2011-06-02 Created by PK from files setbc4.f, phintp4.f,
      module setbcmod
      use qlog
      use pointers
      use operators_on_coordinates
      use misc2                       ! for subroutine phintp   
      implicit none
      contains
      !============================================================================>
      include 'inc_setbcmod_setbc.f95'   ! subroutines setbc 
      !============================================================================>
      include 'inc_setbcmod_relfac.f95'  ! subroutine relfac
      !============================================================================>
      include 'inc_setbcmod_itit.f95'    ! subroutine itit
      !============================================================================>
      include 'inc_setbcmod_conplt.f95'  ! subroutine conplt
      !============================================================================>
      include 'inc_setbcmod_nitit.f95'   ! subroutine nitit
      !============================================================================>
      end module setbcmod
