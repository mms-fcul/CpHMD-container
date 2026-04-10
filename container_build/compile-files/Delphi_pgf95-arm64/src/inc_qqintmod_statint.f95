      !  ATTENTION!  This file is part of qqinttmod module.
      !==================================================================
      subroutine statint(line,lineln)
      character(100) ::  line
      character(30)  :: type=" "                                
      character(20)  :: quant=" "                               
      character(26) :: calph1="abcdefghijklmnopqrstuvwxyz"      
      character(26) :: calph2="ABCDEFGHIJKLMNOPQRSTUVWXYZ"      
      character(2) :: stat2(100)                                
      character(2) :: typ2                                      
      character(10) :: stat10(100)                                
      character(10) :: typ10
      character :: chr                                          
      integer :: lineln,typenum,statnum                         
      logical :: itf,icon                                       
      integer :: i,j,ieq,ieq1,ieq2
      !-------------------------------------------------------------------
      !!c have to decide whether to use decode or internal read
      !!c the latter works better (and is easier!) for real quantities
      !!c because if a period is absent from the number it assumes it
      !!c is a very small number. this must be (?) an error. anyway
      !!c read SEEMS to work so far.
      !!c number of different statements
      
      !!c data list of two letter codes
      stat2=' '; stat10=' '
      stat2(1)="GS"
      stat2(2)="SC"
      stat2(3)="PF"
      stat2(4)="ID"
      stat2(5)="ED"
      stat2(6)="PR"
      stat2(7)="IR"
      stat2(8)="IS"
      stat2(9)="BC"
      stat2(10)="LI"
      stat2(11)="NI"
      stat2(12)="MD"
      stat2(13)="FC"
      stat2(14)="LP"
      stat2(15)="LG"
      stat2(16)="CI"
      stat2(17)="CF"
      stat2(18)="PX"
      stat2(19)="PY"
      stat2(20)="PZ"
      stat2(21)="AC"
      stat2(22)="XU"
      stat2(23)="GC"
      stat2(24)="RF"
      stat2(25)="CI"
      stat2(26)="SP"
      stat2(27)="CS"
      !!c b+++++++++++++++++++++++
      stat2(29)="RL"
      stat2(30)="RR"
      stat2(31)="S2"
      stat2(32)="R2"
      stat2(33)="+1"
      stat2(34)="-1"
      stat2(35)="+2"
      stat2(36)="-2"
      stat2(37)="MC"
      stat2(38)="XC"
      stat2(39)="NC"
      stat2(40)="VX"
      stat2(41)="VY"
      stat2(42)="VZ"
      stat2(43)="AD"
      stat2(44)="TE"
      !!c e+++++++++++++++++++++++
      stat2(45)="CF"
      stat2(46)="SG"
      stat2(47)="IH"
      stat2(48)="SF"
      stat2(49)="GA"
      
      do i=1,100
         if(stat2(i).ne." ") statnum=i
      end do
      !!c data list of six letter codes
      
      stat10(1)="GSIZE"
      stat10(2)="SCALE"
      stat10(3)="PERFIL"
      stat10(4)="INDI"
      stat10(5)="EXDI"
      stat10(6)="PRBRAD"
      stat10(7)="IONRAD"
      stat10(8)="SALT"
      stat10(9)="BNDCON"
      stat10(10)="LINIT"
      stat10(11)="NONIT"
      stat10(12)="MEMDAT"
      stat10(13)="FCRG"
      stat10(14)="LOGPOT"
      stat10(15)="LOGGRP"
      stat10(16)="CONINT"
      stat10(17)="CONFRA"
      stat10(18)="PBX"
      stat10(19)="PBY"
      stat10(20)="PBZ"
      stat10(21)="AUTOC"
      stat10(22)="EXITUN"
      stat10(23)="GRDCON"
      stat10(24)="RELFAC"
      stat10(25)="CHEBIT"
      stat10(26)="SOLVPB"
      stat10(27)="CLCSRF"
      stat10(28)="PHICON"
      !!c b+++++++++++++++++++++++
      stat10(29)="RADPOL"
      stat10(30)="RELPAR"
      stat10(31)="SALT2"
      stat10(32)="RADPR2"
      stat10(33)="VAL+1"
      stat10(34)="VAL-1"
      stat10(35)="VAL+2"
      stat10(36)="VAL-2"
      stat10(37)="RMSC"
      stat10(38)="MAXC"
      stat10(39)="NORMC"
      stat10(40)="VDROPX"
      stat10(41)="VDROPY"
      stat10(42)="VDROPZ"
      stat10(43)="ATPODS"
      stat10(44)="TEMPER"
      stat10(45)="CUTOFF"
      stat10(46)="SIGMA"
      stat10(47)="INHOMO"
      stat10(48)="SRFCUT"
      stat10(49)="GAUSSIAN"
      !!c e+++++++++++++++++++++++
      
      !!c convert type to uppercase
      
      do i=1,lineln
         j=index(calph1,line(i:i))
         if(j.ne.0) line(i:i)=calph2(j:j)
      end do
      type=' ';quant=' '
      ieq=index(line(1:lineln),"=")
      ieq1=ieq-1 ; ieq2=lineln-ieq
      type(1:ieq1)=line(1:ieq1)
      quant(1:ieq2)=line(ieq+1:lineln)
      typenum=0
      
      !!c two and six letter codes are matched exactly
      !!c otherwise the statement type is compared to a list of longer names
      !!c where the first letter is compared, then the rest, tocut down on work
      
      if(ieq1.eq.2) then
         typ2=type(1:2)
         do i=1,statnum
            if(typ2.eq.stat2(i)) typenum=i
         end do
      end if
      
      if((ieq1.le.10).and.(ieq1.gt.2)) then
         typ10=type(1:10)
         do i=1,statnum
            if(typ10.eq.stat10(i)) typenum=i
         end do
      end if
      
      if(ieq1.gt.6) then
         select case(type(1:1)) ! 2011-04-17 instead of multiply IFs
         case ("A")             ! 2011-04-17 instead of if(chr.eq."A") then
            if(type(1:ieq1).eq."AUTOMATICCONVERGENCE") typenum=21
            if(type(1:ieq1).eq."AUTOCONVERGENCE") typenum=21
            if(type(1:ieq1).eq."AUTOCON") typenum=21
            if(type(1:ieq1).eq."ATOMPOTDIST") typenum=43
         case("B")              ! 2011-04-17 instead of if(chr.eq."B") then
            if(type(1:ieq1).eq."BOXFILL") typenum=3
            if(type(1:ieq1).eq."BOUNDARYCONDITION") typenum=9
            if(type(1:ieq1).eq."BOUNDARYCONDITIONS") typenum=9
         case("C")           ! 2011-04-17 instead of if(chr.eq."C") then
            if(type(1:ieq1).eq."CONVERGENCEINTERVAL") typenum=16
            if(type(1:ieq1).eq."CONVERGENCEFRACTION") typenum=17
         case("D")  
            
         case("E")
            if(type(1:ieq1).eq."EXTERIORDIELECTRIC") typenum=5
            if(type(1:ieq1).eq."EXTERNALDIELECTRIC") typenum=5
            if(type(1:ieq1).eq."EXITUNIFORMDIELECTRIC") typenum=22
         case("F")
            if(type(1:ieq1).eq."FANCYCHARGE") typenum=13
         case("G")
            if(type(1:ieq1).eq."GRIDSIZE") typenum=1
            if(type(1:ieq1).eq."GRIDCONVERGENCE") typenum=23
         case("H")
         case("I")
            if(type(1:ieq1).eq."IONRADIUS") typenum=7
            if(type(1:ieq1).eq."IONICSTRENGTH") write(6,*)&
            &'the option IONICSTRENGTH is no longer &
            &available, try SALT or SALT1'
            ! no longer available  typenum=8
            if(type(1:ieq1).eq."INTERIORDIELECTRIC") typenum=4
            if(type(1:ieq1).eq."ITERATIONS") typenum=10
            if(type(1:ieq1).eq."ITERATION") typenum=10
         case("J")
            
         case("K")
            
         case("L")
            if(type(1:ieq1).eq."LINEARITERATION") typenum=10
            if(type(1:ieq1).eq."LINEARITERATIONS") typenum=10
            if(type(1:ieq1).eq."LOGFILEPOTENTIALS") typenum=14
            if(type(1:ieq1).eq."LOGFILECONVERGENCE") typenum=15
         case("M")
            if(type(1:ieq1).eq."MEMBRANEDATA") typenum=12
            if(type(1:ieq1).eq."MAXCONVERGENCE") typenum=38
         case("N")
            if(type(1:ieq1).eq."NONLINEARITERATION") typenum=11
            if(type(1:ieq1).eq."NONLINEARITERATIONS") typenum=11
            if(type(1:ieq1).eq."NORMCONVERGENCE") typenum=39
         case("O")
         case("P")
            if(type(1:ieq1).eq."PERIODICBOUNDARYX") typenum=18
            if(type(1:ieq1).eq."PERIODICBOUNDARYY") typenum=19
            if(type(1:ieq1).eq."PERIODICBOUNDARYZ") typenum=20
            if(type(1:ieq1).eq."PROBERADIUS") typenum=6
            if(type(1:ieq1).eq."PERCENTBOXFILL") typenum=3
            if(type(1:ieq1).eq."PERCENTFILL") typenum=3
         case("Q")
         case("R")
            if(type(1:ieq1).eq."RELAXATIONFACTOR") typenum=24
            ! b+++++++++++++++++++++++
            if(type(1:ieq1).eq."RADPOLEXT") typenum=29
            if(type(1:ieq1).eq."RELPAR") typenum=30
            if(type(1:ieq1).eq."RMSCONVERGENCE") typenum=37
            ! 2011-04-17 I believe those two statements are misplaced
            ! e++++++++++++++++++++++
         case("S")
            if(type(1:ieq1).eq."SALTCONC") typenum=8
            if(type(1:ieq1).eq."SALTCONCENTRATION") typenum=8
            if(type(1:ieq1).eq."SPHERICALCHARGEDISTRIBUTION") typenum=13
            
         case("T")
            if(type(1:ieq1).eq."TEMPERATURE") typenum=44
         case("U")
         case("V")
         case("W")
         case("X")
         case("Y")
         case("Z")
            
         end select
      end if
      !-----------------------------------------------------------------------------------------
      select case(typenum)     ! 2011-04-17  instead of multiple IFs
      case default             ! 2011-04-17 instead of if(typenum.eq.0) then
         write(6,*) "!!!!!!!!!!!!!"
         write(6,*) "the statement"
         write(6,*) line(1:lineln)
         write(6,*) "could not be interpreted"
         write(6,*) "!!!!!!!!!!!!!"
      case(1)
         read(quant,*) igrid
         !!c	write(6,*) "grid size. ok ",igrid
      case(2)
         read(quant,*) scale
         !!c	write(6,*) "scale. ok",scale
      case(3)
         read(quant,*) perfil
         !!c	write(6,*) "perfil ok",perfil
      case(4)
         read(quant,*) repsin
            if(ideveloper) then
               write(6,*) "epsin. ok",repsin
            else 
               write(6,"(a,f8.2)") " epsin. ok",repsin
            end if
      case(5)
         read(quant,*) repsout
         !!c	write(6,*) "epsout ok",repsout
      case(6)
         !!c b++++++++++++++++++++++++++++++
         read(quant,*) radprb(1)
         !!c	read(quant,'(2f4.2)') radprb(1)
         !!c e++++++++++++++++++++++++++++++
      case(7)
         read(quant,*) exrad
         !!c b++++++++++++++++++++++++++++++
      case(8)
         read(quant,*) conc(1)
         !!c e++++++++++++++++++++++++++++++
      case(9)
         read(quant,*) ibctyp
      case(10)
         read(quant,*) nlit
         iautocon=.false.
         !!c	write(6,*) "nlit ok",nlit
      case(11)
         read(quant,*) nnit
         !!c b++++++++++++++++++++++++++++++++++++++++
         if (nnit.le.20) then
            write(6,*)'At least 30 nonlinear iterations'
            nnit=30
         end if
         !!c e++++++++++++++++++++++++++++++++++++++++
         !!c	write(6,*) "nnit ok",nnit
         
      case(12)
         call yesno(type,ieq1,quant,ieq2,itf)
         imem=itf
         !!c	write(6,*) "imem = ",imem
      case(13)   
         call yesno(type,ieq1,quant,ieq2,itf)
         isph=itf
         !!c	write(6,*) "isph = ",isph
      case(14)
         call yesno(type,ieq1,quant,ieq2,itf)
         ipoten=itf
      case(15)
         call yesno(type,ieq1,quant,ieq2,itf)
         igraph=itf
      case(16)
         read(quant,*) icon1
         !!c	write(6,*) "icon1=",icon1
         
      case(17)
         read(quant,*) icon2
         !!c	write(6,*) "icon2=",icon2
         
      case(18)
         call yesno(type,ieq1,quant,ieq2,itf)
         iper(1)=itf
         !!c	write(6,*) "iper(1)=",iper(1)
         
      case(19)
         call yesno(type,ieq1,quant,ieq2,itf)
         iper(2)=itf
         !!c	write(6,*) "iper(2)=",iper(2)
         
      case(20)
         call yesno(type,ieq1,quant,ieq2,itf)
         iper(3)=itf
         !!c	write(6,*) "iper(3)=",iper(3)
         
      case(21)
         call yesno(type,ieq1,quant,ieq2,itf)
         iautocon=itf
         !!c	write(6,*) "automatic convergence= ",itf
         
      case(22)
         call yesno(type,ieq1,quant,ieq2,itf)
         iexun=itf
         
      case(23)
         read(quant,*)gten
         iautocon=.true.
         
      case(24)
         read(quant,*)uspec
         iuspec=.true.
         
      case(25)
         call yesno(type,ieq1,quant,ieq2,itf)
         icheb=itf
         
      case(26)
         call yesno(type,ieq1,quant,ieq2,itf)
         isolv=itf
         
      case(27)
         call yesno(type,ieq1,quant,ieq2,itf)
         isrf=itf
         
      case(28)
         call yesno(type,ieq1,quant,ieq2,itf)
         iconc=itf
         !!c b++++++++++++++++++++++++++++++
      case(29)
         read(quant,*) radpolext
         ! write(6,*)'radpolext= ',radpolext
      case(30)
         read(quant,*) relpar
         imanual=.true.
      case(31)
         read(quant,*) conc(2)
         !!c now having different probes if the molecule faces the water or not
      case(32)
         read(quant,*) radprb(2)
      case(33)
         read(quant,*) ival(1)
      case(34)
         read(quant,*) ival(2)
      case(35)
         read(quant,*) ival2(1)
      case(36)
         read(quant,*) ival2(2)
      case(37)
         read(quant,*) res1
      case(38)
         read(quant,*) res2
      case(39)
         read(quant,*) res5
      case(40)
         iper(4)=.true.
         read(quant,*) vdrop%x
      case(41)
         iper(5)=.true.
         read(quant,*) vdrop%y
      case(42)
         iper(6)=.true.
         read(quant,*) vdrop%z
      case(43)
         read(quant,*) atompotdist
      case(44)
         read(quant,*) temperature !  conversion from Celsius to absolute
         temperature=temperature+273.15
      case(45)
         read(quant,*) cutoff
      case(46)
         read(quant,*) sigma
      case(47)
         read(quant,*) inhomo
      case(48)
         read(quant,*) srfcut
      case(49)
         read(quant,*) gaussian


      end select
      !!c e+++++++++++++++++++++++++++++
      end subroutine statint
