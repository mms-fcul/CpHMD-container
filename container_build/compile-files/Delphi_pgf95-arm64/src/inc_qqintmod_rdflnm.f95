      !  ATTENTION!  This file is part of qqinttmod module.
      !==================================================================
      subroutine rdflnm(j1,j2,direc,fnam,flen)
      character(len=100) :: direc,calph= ' '       
      integer :: un(10),unum = 0 ,flen,six,five    
      character(len=10) :: cnumb='1234567890'      
      character(len=80) :: fnam                    
      !-------------------------------------------------------
      ! 2011-06-14 Declarations added due to IMPLICIT NONE
      integer :: i,j,k,m,n,j1,j2
      !-------------------------------------------------------
      
      calph(1:40)="ABCDEFGHIJKLMNOPQRSTUVWXYZ.:_-+=!@#$^123"
      calph(41:80)="4567890abcdefghijklmnopqrstuvwxyz|\/?><;"
      flen=0
      fnam="fort." 
      if(j1.ne.0) then
         j=j1 ; n=1 ; k=j+4
         do
            k=k+1
            if(index(cnumb,direc(k:k)).eq.0) exit 
            un(n)=index(cnumb,direc(k:k)); n=n+1
         end do
         
         do i=1,n-1
            fnam(5+i:5+i)=cnumb(un(i):un(i))
         end do
         flen=4+n
      end if
      !!c j2=position of first letter in file descriptor, i.e. the letter
      !!c f of "file"
      if(j2.ne.0) then
         j=j2;  k=j+5
         if(index(calph,direc(k:k)).ne.0) k=k-1
         n=1 ;  m=k+1
         do 
            k=k+1
            if(index(calph,direc(k:k)).eq.0) exit
            n=n+1
         end do
         
         fnam(1:n-1)=direc(m:j+4+n)
         flen=n-1
      end if
      
      end subroutine rdflnm    ! changed to F95 syntax on 2011-04-15
