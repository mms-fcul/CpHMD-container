      !  ATTENTION!  This file is part of setrcmod module.
      subroutine chkcrg
      character*80 fmt
      character*4 err,strng,resnme,resnum,rsnm
      !-------------------------------------------------------------------
      ! 2011-06-14 Declarations added due to IMPLICIT NONE
      real :: crgsm,crg,eror
      integer :: i,itmp
      !-------------------------------------------------------------------
      crgsm=0.0
      rsnm='    '
      resnme='    '
      fmt='(''!!! WARNING: '',2a4,'' has a net charge of '',f8.4)'
      
      ! b+++++++++++++++++++++++++++++++++++++++++++++++
      if (isitsf) then
         do i=1,natom
            !  2011-05-05   Values are taken from derived-type array 
            crg=delphipdb(i)%chrgv4
            resnum=delphipdb(i)%atinf(12:15)
            strng=delphipdb(i)%atinf(7:10)
            if(resnum.ne.rsnm.or.strng.ne.resnme)then
               if(abs(crgsm).gt.1.0e-4)then
                  eror=abs(crgsm)-1.0
                  ! 2011-05-05  Label 10 to format statement removed
                  if(abs(eror).gt.1.0e-4) &
                  &write(6,trim(fmt))resnme,rsnm,crgsm
               endif
               resnme=strng
               rsnm=resnum
               crgsm=crg
            else
               crgsm=crgsm+crg
            endif
            ! specific section to find maximum residue number
            ! 2011-05-05  Goes here when isitsf = .true.
            read(resnum,'(i4)')itmp
            if(itmp.gt.resnummax) resnummax=itmp
         end do
      else
         ! e++++++++++++++++++++++++++++++++++++++++++++
         ! 2011-05-05 Seemingly redundant piece of code
         do i=1,natom
            crg=delphipdb(i)%chrgv4
            resnum=delphipdb(i)%atinf(12:15)
            strng=delphipdb(i)%atinf(7:10)
            if(resnum.ne.rsnm.or.strng.ne.resnme)then
               if(abs(crgsm).gt.1.0e-4)then
                  eror=abs(crgsm)-1.0
                  if(abs(eror).gt.1.0e-4) &
                  &write(6,trim(fmt))resnme,rsnm,crgsm
               endif
               resnme=strng
               rsnm=resnum
               crgsm=crg
            else
               crgsm=crgsm+crg
            endif
         end do
      endif
      
      if(abs(crgsm).gt.1.0e-4)then
         eror=abs(crgsm)-1.0
         if(abs(eror).gt.1.0e-4) &
         &    write(6,trim(fmt))resnme,rsnm,crgsm
      endif
      end subroutine chkcrg