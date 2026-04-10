      !  ATTENTION!  This file is part of setcrgmod module.
      !==================================================================
      subroutine wrteps(imaxwrd,epsdim)
      !!c	kim sharp/9 feb 88
      !!c	reformats epsilon array epsmap
      !!c	to compact array (MKG format)
      !!c	old format: epsmap(65,65,65,3) real*4
      !!c	new format: neps(5,65,65,3) integer*2
      !!c	where first index 1-65 now compressed
      !!c	into 1-5 plus offset into 16 bit words
      !!c	compact format also contains oldmid, the
      !!c	center of the protein in real coordinates
      !!c	compaction is effected by storing
      !!c	real eps (which take values of 0. and 1.)
      !!c	as bits in a 16 bit word
      !!c	access is via pointers idimx and ioffset
      !!c	thus x arrary indices of reps 0-15 -> word 1
      !!c	16-31 -> word 2 etc
      	
      
      !!c	 parameter  (imaxwrd = igrid/16 + 1)
      ! 2011-06-05 Arrays are declared in pointers module and
      !	integer iepsmp(igrid,igrid,igrid,3),epsdim
      !	logical*1 idebmap(igrid,igrid,igrid)
      
      !!c attenzione!!! le dimensioni non sono quelle che sembrano! int*2!!!
      !	dimension neps(imaxwrd,igrid,igrid,3)	
      !	dimension keps(imaxwrd,igrid,igrid)	
      
      !!c	compact fine epsilon map
      ! 2011-06-05 Probable bug, should be igrid instead of ngrid
      integer :: idimx(igrid),ip2(igrid)
      !!c	array of pointers to words
      integer :: ioffset(igrid)		
      !!c	array of pointers to bit offsets
      	integer ::  i,j
      	character*80 filnam
      
      !!c	integer epsmap(65,65,65)
      !---------------------------------------------------------------------
      ! 2011-06-05 Declarations added due to IMPLICIT NONE
      integer :: epsdim,imaxwrd,ix,iy,iz,k1,kmap,imap
      type(int_coord) :: j123
      !!c---------------------------------------------------------------------
      ! 2011-06-05 Arrays are declared in pointers module 
      allocate(neps(imaxwrd,igrid,igrid),keps(imaxwrd,igrid,igrid)) 
      
      write(6,*)' setting up pointers...'
      do ix = 1, igrid
         idimx(ix) = ix/16 + 1; ioffset(ix) = mod(ix,16)
         ip2(ix)=2**(mod(ix,16))
      end do
      ! 2011-06-05 Changed to array operations
      ip2(15:igrid:16)=-2**15
      write(6,*)' clearing bits...'
      neps=int_coord(0,0,0) ; keps=0
      write(6,*)' generating compact fine epsilon array...'
      do iz=1,igrid
         do  ix=1,igrid
            i=idimx(ix)
            do  iy=1,igrid
               j123=int_coord(0,0,0)
               k1=0
               !!c divide solvente da non solvente
               if((iepsmp(ix,iy,iz)%i/epsdim).ne.0) j123%i=ip2(ix)
               if((iepsmp(ix,iy,iz)%j/epsdim).ne.0) j123%j=ip2(ix)
               if((iepsmp(ix,iy,iz)%k/epsdim).ne.0) j123%k=ip2(ix)
               if(idebmap(ix,iy,iz)) k1=ip2(ix)
               ! 2011-06-05 Another bug? SHould it be ix instead of i?
               neps(i,iy,iz)= neps(i,iy,iz)+j123 !Lin chaged it from ix to i. 2012-06-12
               keps(i,iy,iz)=keps(i,iy,iz)+k1 !Lin changed it from ix to i
            end do
         end do
      end do
      
      kmap = 1
      
      write(6,*)' writing to compact epsilon file'
      open(17,file=trim(epsnam),form='unformatted')
      filnam = ' '
      ! 2011-06-05 Did not understand the meaning of this statement
      inquire(17,name = filnam)
      write(6,*)' '
      write(6,*)'dielectric map written to file'
      write(6,*)filnam
      write(6,*)' '
      imap = 0
      write (17) kmap, scale, oldmid
      write (17) neps
      write (17) keps
      close (17)
      
      !!c b++++++++++++++++++++++++++++++++++++++++++
      if (.false.) then
         write(6,*)' PER STEFANO'
         open(17,file='StefanoDebmap.bin',form='unformatted')
         write (17) idebmap
         close (17)
      end if
      
      if (debug) then
         open(52,file='iepsmapnew',form='unformatted')
         write (52) kmap, scale, oldmid
         write (52) iepsmp
         close (52)
      end if
      !!c e++++++++++++++++++++++++++++++++++++++++++
      deallocate(neps,keps)
      end subroutine wrteps