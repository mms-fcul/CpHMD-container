!Created 2011-04-28 by PK. Contains subs called from more than one 

      module misc
      implicit none

      contains

      !------------------------------------------------------------------
      subroutine up(txt,len) !convert character string to upper case

      character(len) :: txt
      character(len) :: temp
      !2011-04-24  DATA statements removed
      character(26) :: ualpha="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      character(26) :: lalpha="abcdefghijklmnopqrstuvwxyz"
      integer :: len,i,match

      temp=" "

      do i=1,len
         if ((txt(i:i).ge.'a').and.(txt(i:i).le.'z')) then
            match=index(lalpha,txt(i:i))
            temp(i:i)=ualpha(match:match)
         else
            temp(i:i) = txt(i:i)
         end if
      end do

      txt = temp

      end subroutine up

      !------------------------------------------------------------------
      subroutine elb(txt,len)

      !eliminate leading blanks from a character string
      !removed GOTO statements
      integer :: len
      character(len) :: txt
      character(len) :: temp
      integer :: i,nonb

      temp=" "

      do i=1,len
         if (txt(i:i).ne.' ') then
            nonb = i
            exit
         end if
      end do

      temp=txt(nonb:len)
      txt=temp

      end subroutine elb

      !------------------------------------------------------------------
      !2011-04-24 added parameter Nrec, number of records in charge or 
      !           radii table, needed for correct hash table 
      !           normalization if number of records in charge and 
      !           radii file differs
      !           This function is called from rdh and setrc modules
      !
      !produce hash number from atom and residue name,residue number and
      !chain name for charge assignment
      !------------------------------------------------------------------   
      integer function ichash(atxt,rtxt,ntxt,ctxt,Nrec)

      character(6) :: atxt
      character(3) :: rtxt
      character(4) ::  ntxt
      character :: ctxt
      integer :: i,j,n=1,Nrec
      character*38 :: string="* 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"

      do i=1,3
         j=index(string,rtxt(i:i))
         n=5*n+j
      end do

      do i=1,6
         j=index(string,atxt(i:i))
         n=5*n+j
      end do

      do i=1,4
         j=index(string,ntxt(i:i))
         n=5*n+j
      end do

      do i=1,1
           j=index(string,ctxt(i:i))
           n=5*n+j
      end do

      n=iabs(n)
      ichash=mod(n,Nrec)+1 !see comment at the function start

      end function ichash

      !------------------------------------------------------------------
      subroutine form(fname,nam1,ifrm)

      character(80) :: fname,line,asci
      integer :: len
      logical :: ifrm
      integer :: ias,i,ios,nam1
      asci(1:40)="1234567890 .-+#,$asdfghjklzxcvbnmqwertyu"
      asci(41:80)="iopASDFGHJKLZXCVBNMQWERTYUIOP)(}{][/    "

      !2011-04-30 nam1 is obsolete used trim instead, 204 format removed
      open(9,file=fname(:nam1),form='formatted')

      !2011-04-30  Used iostat, instead of err- and end=
      read(9,'(a)',iostat=ios)line           
      close(9)

      if (ios.eq.0) then !line read successfully
         ias=0

         do i=1,80
            if (index(asci,line(i:i)).eq.0) then
               ias=ias+1
            end if
         end do

         if (ias.gt.10) then
            ifrm=.false.
         else
            ifrm=.true.
         end if

         return      
      else if (ios.lt.0)  then !2011-04-30 instead of err=500
         write(6,*) "unexpected error reading pdb file"
         write(6,*) "assuming formatted file!"
         ifrm=.true.
         return !2011-04-30  instead of goto 1000
      else !2011-04-30  instead of end=600 
         write(6,*) "unexpected end of pdb file!"
         write(6,*) "assuming formatted file!"
         ifrm=.true.
         return !2011-04-30 instead goto  1000
      end if

      end subroutine form

      end module misc
