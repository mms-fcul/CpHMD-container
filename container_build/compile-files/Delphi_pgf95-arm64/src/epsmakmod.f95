!#####################################################################
! 2011-05-10-27 Created by PK from subroutine epsmak and subroutines 
!               within.
!#####################################################################
      module epsmakmod
      
      use qlog
      use pointers
      use operators_on_coordinates
      implicit none
 
      contains
!====================================================================>
      include 'inc_epsmakmod_epsmak.f95'      !subroutine epsmak
!====================================================================>
      include 'inc_epsmakmod_setout.f95'      !subroutine setout
      include 'inc_epsmakmod_setgaussian.f95'      !subroutine setout
!====================================================================>
      include 'inc_epsmakmod_vwtms.f95'       !subroutine vwtms
!====================================================================>
      include 'inc_epsmakmod_sclbp.f95'       !subroutine sclbp
!====================================================================>
      include 'inc_epsmakmod_msrf.f95'        !subroutine msrf
!====================================================================>
      include 'inc_epsmakmod_indver.f95'      !subroutines indver and
                                              !indverdata
!====================================================================>
      include 'inc_epsmakmod_mkvtl.f95'       !subroutine mkvtl
!====================================================================>
      include 'inc_epsmakmod_fxhl.f95'        !subroutine fxhl
!====================================================================>
      include 'inc_epsmakmod_ex.f95'          !subroutine ex
!====================================================================>
      include 'inc_epsmakmod_sas.f95'         !subroutine sas
!====================================================================>
      include 'inc_epsmakmod_cube.f95'        !subroutines cube and 
                                              !cubedata
!====================================================================>
      include 'inc_epsmakmod_objvertices.f95' !subroutine objvertices
!====================================================================>
      include 'inc_epsmakmod_distobj.f95'     !subroutine distobj
!====================================================================>
      end module epsmakmod
