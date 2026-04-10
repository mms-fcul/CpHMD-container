      !  ATTENTION!  This file is part of qqinttmod module.
      !==================================================================
      subroutine rdprm
      character(80) line,asci,enc,filnam    
      character(10) :: cnumb="0123456789"   
      character(400) :: mline               
      integer, dimension(20) :: lin1,itype  
      logical :: iext                       
      real ::  cz1,cz2,z1p,z1m,z2p,z2m      
      integer :: ifn,ias,iblen,i,ifillen,ios,j1,j2,j3,j4,ilen
      integer :: n
      real :: dfact
      !-------------------------------------------------------------------
      !!c read in run parameters
      prmnam(prmlen+1:len(prmnam))=' '
      asci="1234567890 .-+#,$asdfghjklzxcvbnmqwertyuio"
      asci=trim(asci)//"pASDFGHJKLZXCVBNMQWERTYUIOP)(}{][/"

      inquire(file=prmnam(:prmlen),exist=iext)
      if(.not.iext) then
         write(6,*) "parameter file ",trim(prmnam),&
         &" does not seem to exist."
         write(6,*)'Therefore stopping'
         stop
      end if
      
      write(6,*) "opening parameter file ",trim(prmnam)
      open(10,file=prmnam(:prmlen))
      
      !!c set current file length to prmlen
      !!c set current file to prmnam
      ifn=10; ifillen=prmlen
      filnam(:ifillen)=prmnam(:prmlen)
      
      !!c goto here if a new file has just been opened
      !!c and reset line number to 1 (only used by type 1 files)
      
      D10: do
         itype(ifn)=0; lin1(ifn)=1
         
         !!c goto here to read a new line from an already opened file
         
         D20 : do
            line(1:80)=' '
            read(ifn,'(a80)',iostat=ios) line
            if(ios.ne.0) then
               if(ios.gt.0) write(6,*) "trouble at mill"
               close(ifn)
               if(ifn.gt.10) then
                  ifn=ifn-1; cycle D20
               else
                  exit D10
               end if
            end if
            !!c added 9/9/92 to fixed qincludes being not commented out!
            if(line(1:1).eq."!") line(2:80)=" "
            
            !!c do below if the type of the prm file is undetermined
            !!c normally the first line, unless that was a qinclude
            !!c statement
            
            if(itype(ifn).eq.0) then
               !!c check as to whether the first character of the file is a number
               !!c	write(6,*) "checking for type"
               if(index(cnumb,line(1:1)).eq.0) then
                  j1=index(line,"QINCLUDE")
                  j2=j1+index(line,"qinclude")
                  if(j2.ne.0) then
                     !!c	write(6,*) "there is an include file"
                     !!c there is a qinclude statement before the current file type is
                     !!c determined. is this a command passer or a file opener?
                     j4=len_trim(line)  
                     j3=index(asci,line(j2+9:j2+9))
                     if(j3.eq.0) then
                        !!c the qinclude statement does not open a file, but passes a
                        !!c command, therefore pass that command
                        call ppi(line(j2+10:j4-2),j4-2-j2-9)
                     else      ! 2011-04-17 instead of if(j3.ne.0) then
                        !	write(6,*) "it opens a file"
                        ! the include statement opens a file, therefore increment ifn
                        ! go back to 10 if it exists
                        filnam=line(j2+9:j4-1)
                        ifillen=j4-j2-9
                        if(filnam.eq." ") then
                           ! 2011-04-17 instead of call system("echo $HOME > tmp101")
                           call getenv('HOME',uhomenam) 
                           uhomelen=len_trim(uhomenam)  
                           filnam=uhomenam(:uhomelen)//"/qpref.prm"
                           ifillen=uhomelen+10
                           write(6,*) "opening default file qpref.prm"
                        end if
                        inquire(file=filnam(:ifillen),exist=iext)
                        if(iext) then
                           ifn=ifn+1
                           write(6,*) "opening file: ",filnam(:ifillen)
                           open(ifn,file=filnam(:ifillen))
                           cycle D10
                        end if
                        if(.not.iext) then
                           write(6,*) "the file ",filnam(:ifillen)," specified in"
                           write(6,*) &
                           &"a qinclude statement does not exist. continuing"
                        end if
                     end if
                     !!c end of dealing with qinclude
                     
                  else ! 2011-04-17  instead of if(j2.eq.0) then
                     
                     !!c no qinclude, therefore define this as a type 2 parameter file
                     !!c deal with possible multiply lines, deal with errors or
                     !!c errroneous end of files
                     itype(ifn)=2;  iblen=1
                     D30: do
                        ilen=len_trim(line)   
                        !	write(6,*) "passing first line",line(:ilen)
                        if(line(ilen:ilen).eq."&") then
                           mline(iblen:iblen+ilen-2)=line(:ilen-1)
                           iblen=iblen+ilen-1
                           read(ifn,'(a80)',iostat=ios) line
                           if(ios.ne.0) exit D30
                           cycle D30  !   goto 30
                        else
                           mline(iblen:iblen+ilen-1)=line(:ilen)
                           exit D30
                        end if
                     end do D30
                     if(ios.ne.0) then
                        if(ios.gt.0) then
                           write(6,*) "error reading file :",trim(filnam)
                           write(6,*) "continuing.."
                        else
                           write(6,*) "error, continuation slash on last line of"
                           write(6,*) "file ",trim(filnam)
                           write(6,*) "untoward things may happen"
                        end if
                        close(ifn)
                        if(ifn.gt.10) then
                           ifn=ifn-1; cycle D20
                        else
                           exit D10
                        end if
                        !!                               goto 100
                     end if
                     call ppi(mline,ilen)
                  end if
                  
               else      !    if(index(cnumb,line(1:1)).ne.0) then
                  !!c the firstcharacter is a number , therefore define this is as a type 1
                  !!c without a qinclude on line 1
                  itype(ifn)=1
                  call prm1(line,lin1(ifn))
                  lin1(ifn)=lin1(ifn)+1
               end if
               
               !!c leave, having determined ttype unless a qinclude was encountered
               cycle D20
            end if
            
            !!c type is already determined
            !!c need to check for qincludes! different for type 1 and type 2
            
            j1=index(line,"QINCLUDE")
            j2=j1+index(line,"qinclude")
            if(j2.ne.0) then
               !!c	write(6,*) "include statement found"
               !!c there is a qinclude statement
               !!c determined. is this a command passer or a file opener?
               j4=len_trim(line)     !  2011-04-17 instead of call namleb(line,j4)
               j3=index(asci,line(j2+9:j2+9))
               if(j3.eq.0) then
                  !!c the qinclude statement does not open a file, but passes a
                  !!c command, therefore pass that command
                  call ppi(line(j2+10:j4-2),j4-2-j2-9)
               else     !2011-04-17 instead of if(j3.ne.0) then
                  !!c the include statement opens a file, therefore increment ifn
                  !!c go back to 10 if it exists
                  filnam=line(j2+9:j4-1)
                  ifillen=j4-j2-9
                  if(filnam.eq." ") then
                     call getenv('HOME',uhomenam) 
                     uhomelen=len_trim(uhomenam)  
                     filnam=uhomenam(:uhomelen)//"/qpref.prm"
                     ifillen=uhomelen+10
                     write(6,*) "opening default file qpref.prm"
                  end if
                  inquire(file=filnam(:ifillen),exist=iext)
                  
                  if(iext) then
                     ifn=ifn+1
                     write(6,*) "opening file: ",filnam(:ifillen)
                     open(ifn,file=filnam(:ifillen))
                     cycle D10
                  else      !  2011-04-17 instead of if(.not.iext) then
                     write(6,*) "the file ",filnam(:ifillen)," specified in"
                     write(6,*) "a qinclude statement does not exist. continuing"
                  end if
                  
                  !!c end of opening fil "if",j3
               end if
               cycle D20
            end if
            
            if(itype(ifn).eq.1) then
               call prm1(line,lin1(ifn))
               lin1(ifn)=lin1(ifn)+1
            end if
            
            if(itype(ifn).eq.2) then
               iblen=1
               D31: do
                  ilen=len_trim(line)   
                  if(line(ilen:ilen).eq."&") then
                     mline(iblen:iblen+ilen-2)=line(:ilen-1)
                     iblen=iblen+ilen-1
                     read(ifn,'(a80)',iostat=ios) line
                     if(ios.ne.0) exit D31
                     cycle D31
                  else
                     mline(iblen:iblen+ilen-1)=line(:ilen)
                     exit D31
                  end if
               end do D31
               if(ios.ne.0) then
                  if(ios.gt.0) then
                     write(6,*) "error reading file :"
                     write(6,*) trim(filnam)
                     write(6,*) "closing file and continuing.."
                  else
                     write(6,*) "error, continuation slash on last line of"
                     write(6,*) "file ",trim(filnam)
                     write(6,*) "which is now closed. untoward things may happen"
                  end if
                  close(ifn)
                  if(ifn.gt.10) then
                     ifn=ifn-1; cycle D20
                  else
                     exit D10
                  end if
               end if
               call ppi(mline,ilen)
            end if
            
            !!c read next line
            cycle D20
            
            !!c deal with end of file
            !500                write(6,*) "trouble at mill"
         end do D20
      end do D10
      
      if((repsin.lt.0.).or.(repsout.lt.0.)) then
         repsin=abs(repsin); repsout=abs(repsout); diff=.true.
      end if
      
      !!c        dfact = 3.047*sqrt(repsout/80.)
      dfact=0.01990076478*sqrt(temperature*repsout)
      !!c b++++++++++++++++++++++
      !!c define a number that indicates whether or not there is some salt 
      z1p=ival(1);  z1m=ival(2)
      z2p=ival2(1); z2m=ival2(2)
      !!c now cz1 and cz2 are concentration of positive ion !!!
      cz1=conc(1)*z1m;  cz2=conc(2)*z2m
      rionst=(cz1*z1p*(z1p+z1m)+cz2*z2p*(z2p+z2m))/2.
      
      !!c coefficients in Taylor series of the charge concentration
      !!c  apart from n! (order >=1)
      !2012-04-24 chuan Correct coefficients in Taylor series. NOT in compact 
      !
      chi1=-2.*rionst
      chi2=cz1*z1p*z1p**2-cz1*z1m*z1m**2+cz2*z2p*z2p**2-&
      &cz2*z2m*z2m**2
      chi2=chi2/2.
      chi3=cz1*z1p*z1p**3+cz1*z1m*z1m**3+cz2*z2p*z2p**3+&
      &cz2*z2m*z2m**3
      chi3=-chi3/6.
      chi4=cz1*z1p*z1p**4-cz1*z1m*z1m**4+cz2*z2p*z2p**4-&
      &cz2*z2m*z2m**4
      chi4=chi4/24.
      chi5=cz1*z1p*z1p**5+cz1*z1m*z1m**5+cz2*z2p*z2p**5+&
      &cz2*z2m*z2m**5
      chi5=-chi5/120.
      
      !!c convert ionic strength to debye length
      
      if(rionst.gt.1.e-6) then
         deblen = dfact/sqrt(rionst) 
         if (nnit.gt.0) lognl=.true.
      else
         logions=.false.; deblen = 1.e6
      end if
      
      !!c test for unformatted pdb and frc files
      if(.not.ipdbrd) then
         open(13,file=pdbnam(:pdblen),form='formatted')
         read(13,'(a80)',iostat=n)line     
         ias=0
         do  i=1,80
            if(index(asci,line(i:i)).eq.0) ias=ias+1
         end do ! instead of 600 continue
         if(ias.gt.10) ipdbrd=.true.
         close(13)
      end if
      
      if(ifrcwrt) then
         open(15,form='formatted')
         read(13,'(a80)',iostat=n)line     
         ias=0
         do i=1,80
            if(index(asci,line(i:i)).eq.0) ias=ias+1
         end do                            
         if(ias.gt.10) ifrcrd=.true.
         close(15)
      end if
      
      !!c epkt assignment as a function of temperature
      epkt=167100.9162872952/temperature
      !!c set epsin and epsout (=epkt adjusted dielectrics such that
      !!c all distances are in angstroms, charges in e)
      
      epsin = repsin/epkt; epsout = repsout/epkt
      end subroutine rdprm
