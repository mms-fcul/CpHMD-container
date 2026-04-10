      !  ATTENTION!  This file is part of qqinttmod module.
      !==================================================================
      !  Subroutine createpdb, from file creapdb.c
      subroutine createpdb(prepsout,pnumbmol)  
      
      integer :: mediacount=0, link,kind
      ! kind=0 dielectric; kind=1 metal; kind=2 pore,kind=3 fictious,end of pore 
      real :: medeps(0:80)
      character(10) :: filename
      character(160) :: line=" "
      character(7) :: str
      character :: ans,choice,c
      logical ex
      !-------------------------------------------------------------------
      ! 2011-06-14 Declarations added due to IMPLICIT NONE
      real :: prepsout
      integer :: pnumbmol,number,n,count,distrtype,flag,irstat
      integer :: n1,n2,tmp,objecttype,count1
      !-------------------------------------------------------------------
      pnumbmol=0
      medeps(mediacount)=prepsout
      write(6,*)"External Dielectric Constant ",medeps(mediacount)
      
      write(6,'(a)',ADVANCE='NO')&
      &"Do you want to overwrite file fort.13, if present? [y/n]:"
      read(5,*) ans
      
      if (ans.ne.'y') return
      !  2011-04-19
      !  Strictly speaking in fortran this statement is not necessary since
      !  file fort.13 by default is associated with the stream Nr.13
      open(13,file='fort.13',status='unknown',iostat=irstat)
      write(13,*)
      count=1
      ! cycle on different objects
      do while(ans.eq.'y'.or.ans.eq.'Y') 
         flag=0; kind=0
         ! writes down the line: OBJECT #object #objecttype #media */
         write(13,'(a,i3)',advance='no') "OBJECT ",count
         
         write(6,*)"Input objecttype number : "
         write(6,*)"(molecule from pdb file [0], "
         write(6,*)"sphere [1], cylinder [2], "
         write(6,'(a)',advance='no')"cone [3], box [4]) "
         read(5,*) objecttype
         write(13,'(i4)',advance='no') objecttype
         
         tmp=-1
         write(6,'(a)',advance='no')&
         &("Input internal dielectric constant : ")
         read(5,*) repsin
         do count1=1,mediacount
            if (repsin.eq.medeps(count1)) tmp=count1
         end do
         if(tmp.eq.-1) then
            mediacount=mediacount+1
            medeps(mediacount)=repsin
            tmp=mediacount  
         end if
         write(13,'(i3,f8.3)') tmp,repsin   !  writes medianumber and epsilon
         line=" "
         select case(objecttype)
         case (0)
            pnumbmol=pnumbmol+1
            write(6,'(a,i4)',advance='no')"insert filename of object : ",count     
            read(5,'(a)') filename
            inquire(file=trim(filename),EXIST=ex)
            if(ex) then
               open(10,file=trim(filename),status='old')
               do while(n.eq.0)
                  read(10,'(a)',iostat=n) line
                  if(line(1:4).eq.'ATOM'.or.line(1:6).eq.'HETATM')&
                  & write(13,'(a)') trim(line)
                  close(13)
               end do
            else
               write(6,*) "can't open this file"
            end if
         case default
            write(6,'(a,i3,a)',advance='no') &
            &"Insert data of object",count," (spaced by commas):"
            read(5,'(a)') line
            write(13,'(a)',advance='no') "DATA  "
            write(13,'(3i4,1x,a)') objecttype,tmp,kind,trim(line)
            if (kind.eq.2.and.objecttype.eq.2) then
               count=count+1; number=20
               write(13,'(a)',advance='no') "OBJECT"
               write(13,'(3i4,1x,f8.3)')count,number,tmp,repsin
               n1=10; n2=3
               write(13,'(a)',advance='no') "DATA  "
               write(13,'(3i4,1x,a)')n1,tmp,n2,trim(line)
               !    2011-04-19 Transferred form creapdb.c, but is it meaningful to repeat?
               count=count+1; number=20
               write(13,'(a)',advance='no') "OBJECT"
               write(13,'(3i4,1x,f8.3)')count,number,tmp,repsin
               n1=10; n2=3
               write(13,'(a)',advance='no') "DATA  "
               write(13,'(3i4,1x,a)')n1,tmp,n2,trim(line)
            end if
         end select
         write(6,'(a)',advance='no')"Are you going to insert a new object?(y/n): "
         read(5,*)ans
         count=count+1
      end do
      write(6,'(a)',advance='no')&
      &"Are you going to insert charge distributions?(y/n): "
      read(5,*) ans
      if(ans.eq.'y'.or.ans.eq.'Y') then
         
         count=1
         ! cycle on different distributions 
         do while((ans.eq.'y').or.ans.eq.'Y')
            
            !write down the line: CRGDST #distribution              */
            write(13,'(a)',advance='no') "CRGDST "; write(13,'(i3)') count
            
            write(6,*)"Input distrtype number : "
            write(6,*)"Shape: sphere [1], cylinder [2], cone [3], box [4], "
            write(6,*)"point charge [8], segment [9], disk [10], &
            &rectangular plate [11] "
            read(5,*) distrtype
            write(6,'(a)',advance='no')&
            &"Is it a volumic (v) or surfacial (s) charge? : "
            read(5,*)choice
            write(6,'(a)',advance='no')&
            &"Input objectnumber to which this distribution is linked &
            &(0 if is free) : "
            read(5,*)link
            write(6,'(a)',advance='no')&
            &"Input total charge in this distribution : "
            read(5,*) charge
            
            line=" "
            
            write(6,'(a,i3)',advance='no')&
            &"Insert data of distribution ",count,":"
            read(5,*) line
            write(13,'(a)',advance='no')"DATA   "
            write(13,'(i3,1x,a1,1x,i3,1x,f8.3)') &
            &distrtype,choice,link,charge,line
            
            write(6,'(a)',advance='no')&
            &"Are you going to insert a new distribution?(y/n): " 
            read(5,*)ans
            count=count+1;
         end do
      end if
      ! write number of different media at the beginning of file*/
      rewind(13);
      write(13,'(a)',advance='no')"MEDIA  "
      write(13,'(i3)') mediacount
      close(13);
      end subroutine createpdb
