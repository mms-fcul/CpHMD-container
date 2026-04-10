      !  2011-04-15 Created by PK from files qinttot.F
      
      module qqintmod
      use qlog 
      use pointers
      use misc
      implicit none
      contains   
      !============================================================================>
      include 'inc_qqintmod_qqint.f95'   ! subroutine qqint
      !============================================================================>
      include 'inc_qqintmod_defprm.f95'  ! subroutine defprm
      !============================================================================>
      include 'inc_qqintmod_rdflnm.f95'  ! subroutine rdflnm
      !============================================================================>
      include 'inc_qqintmod_funcint.f95' ! subroutine funcint
      !============================================================================>
      include 'inc_qqintmod_ppi.f95'      ! subroutine ppi
      !============================================================================>
      include 'inc_qqintmod_prm1.f95'    ! subroutine prm1
      !============================================================================>
      include 'inc_qqintmod_rdprm.f95'   ! subroutine rdprm
      !============================================================================>
      include 'inc_qqintmod_statint.f95' ! subroutine statint
      !============================================================================>
      include 'inc_qqintmod_yesno.f95'   ! subroutines yesno and spuchar
      !============================================================================>
      include 'inc_qqintmod_createpdb.f95'   ! subroutines yesno and spuchar
      !============================================================================>
      end module qqintmod
