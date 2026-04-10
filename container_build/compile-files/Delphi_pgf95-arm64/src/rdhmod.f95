      ! module rdh created by PK from rdhrad.f and rdhcrg.f
      module rdhmod
      use qlog        ! 2011-04-23 instead of include 'qlog.h'
      use pointers    ! 2011-04-23 instead of include 'qdiffpar4.h
      use misc        ! 2011-04-28 contains various small subs
      implicit none
      contains
      !============================================================================>
      include 'inc_rdhmod_rdhrad.f95'   ! subroutine rdhrad
      !============================================================================>
      include 'inc_rdhmod_rdhcrg.f95'   ! subroutine rdhcrg
      !============================================================================>
      include 'inc_rdhmod_ent.f95'      ! subroutines ent1 and mess1
      !============================================================================>
      end module rdhmod
      
