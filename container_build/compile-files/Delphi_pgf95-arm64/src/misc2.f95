!#####################################################################
! 2011-06-12 Created by PK for subroutines phintp and ass,rfind
!            called from more than one module
!           
!#####################################################################
      module misc2
      use qlog
      use pointers
      use operators_on_coordinates
      use misc
      implicit none
      
      contains
!====================================================================>
      include 'inc_misc2_phintp.f95' !subroutine phintp
!====================================================================>
      include 'inc_misc2_ass.f95'    !subroutine ass (combined radass
                                     !and crgass) and rfind
!====================================================================>

      end module misc2
