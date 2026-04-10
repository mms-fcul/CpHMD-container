!#####################################################################
!2011-06-10 Created by PK from files encalc4.f,react2.f,clb.f,nlener.f
!#####################################################################
      module encalcmod
      use qlog
      use pointers
      use operators_on_coordinates
      implicit none
      
      contains
!====================================================================>
      include 'inc_encalcmod_encalc.f95'  ! subroutine encalc
!====================================================================>
      include 'inc_encalcmod_react.f95'   ! subroutine react
!====================================================================>
      include 'inc_encalcmod_clb.f95'     ! subroutines clb,clbmedia,
                                          ! clbnonl,clbtot
!====================================================================>
      include 'inc_encalcmod_nlener.f95'  ! subroutine nlener
!====================================================================>
      end module encalcmod
