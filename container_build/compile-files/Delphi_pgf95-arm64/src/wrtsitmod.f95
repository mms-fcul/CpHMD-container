      ! 2011-06-12 Created by PK from files wrtsit4.f,rforce.f,ts.f,phintp4.f
      module wrtsitmod
      use qlog
      use pointers
      use operators_on_coordinates
      use misc                   ! for subroutines elb, up, ichash
      use misc2                  ! for subroutines phintp and ass
      implicit none
      contains
      !============================================================================>
      include 'inc_wrtsitmod_wrtsit.f95'  ! subroutine wrtsit
      !============================================================================>
      include 'inc_wrtsitmod_rforce.f95'  ! subroutines rforce,rforceeps1,
      ! rforcenew (not used?)
      !============================================================================>
      include 'inc_wrtsitmod_tops.f95'    ! subroutine tops
      !============================================================================>
      include 'inc_wrtsitmod_debtp.f95'   ! subroutine debtp
      !============================================================================>
      include 'inc_wrtsitmod_wrtphi.f95'  ! subroutines wrtphi,expand,
      ! phicon,submaps,rdphimapo,
      !  wrtphimap
      !============================================================================>
      
      end module wrtsitmod
