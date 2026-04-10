      !   Created 2011-04-28 by PK from setrc4.f, getatm2.f, radass4.f,
      !    crgass4.f, chkcrg.f, form2.f
      module setrcmod
      use pointers
      use qlog
      use misc
      use misc2
      implicit none
      !  2011-04-28   Varialbles used only in this module
      integer, private :: idfrm
      character(80), private :: line,filnam
      character(6), private :: head,str6
      logical, private :: ifrm,iatinf,iatrad,iatcrg
      ! 2011-04-28, declared in pointers module
      contains
      !============================================================================>
      include 'inc_setrcmod_setrc.f95'   ! subroutine setrc
      !============================================================================>
      include 'inc_setrcmod_getatm.f95'   ! subroutine getatm
      !============================================================================>
      include 'inc_setrcmod_chkcrg.f95'   ! subroutine chkcrg
      !============================================================================>
      end module setrcmod
