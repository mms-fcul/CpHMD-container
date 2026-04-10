      !  ATTENTION!  This file is part of qqinttmod module.
      !==================================================================
      subroutine prm1(line,il)
      character(80) :: line   ! changed to F95 syntax on 2011-04-16
      ! 2011-06-14 Declarations added due to IMPLICIT NONE
      integer :: il,i,j,k
      
      if(line.eq." ") return !  2011-04-16 instead of goto 1000
	  
      select case(il)  ! 2011-04-16 instead of computed GOTO statement
         
      case(1)
         read(line,*)igrid
         
      case(2)
         j=index(line,"g")
         i=index(line,"G")
         k=max(i,j)
         if(k.eq.0) then
            read(line,*) perfil
            !	scale=0.
         else
            read(line(:k-1),*) scale
            !	perfil=0.
         end if
      case(3)
         read(line,*)offset
      case(4)
         read(line,*)repsin,repsout
         
         ! b++++++++++++++++++not tested actually+++
      case(5)
         read(line,*)conc(1),conc(2)
      case(6)
         
         read(line,*)exrad,radprb(1),radprb(2),radpolext,relpar,&
         &              atompotdist
         ! da testare, manca imanual!!!
         ! e++++++++++++++++++
         
      case(7)
         read(line,*)ibctyp
      case(8)
         read(line,*)iper
      case(9)
         if((line(1:1).eq.'a').or.(line(1:1).eq.'A')) then
            iautocon=.true. ;  gten=0. ;  nlit=0
         else
            if(index(line,".").eq.0) then
               read(line,*)nlit ; iautocon=.false.; gten=0.
            else
               read(line,*)gten ;  iautocon=.true. ;  nlit=0
            end if
         end if
      case(10)
         read(line,*)nnit
         ! b++++++++++++++++++++++++++++++++++++++++
         if (nnit.le.20) then
            write(6,*)'At least 30 nonlinear iterations'
            nnit=30
         end if
         ! e++++++++++++++++++++++++++++++++++++++++
      case(11)
         read(line,*)iconc,ibios
      case(12)
         read(line,*)isite
      case(13)
         read(line,*)iatout
      case(14)
         toplbl=line(:60)
      case(15)
         read(line,*)isph
      case(16)
         read(line,*)ipdbwrt
      case(17)
         read(line,*)ifrcwrt
      case(18)
         if((index(line,"g").ne.0).or.(index(line,"G").ne.0))logg=.true.
         if((index(line,"s").ne.0).or.(index(line,"S").ne.0))logs=.true.
         if((index(line,"c").ne.0).or.(index(line,"C").ne.0))logc=.true.
         if((index(line,"a").ne.0).or.(index(line,"A").ne.0))loga=.true.
         ! b++++++++++++++++++++++++w Oct 2000
         if((index(line,"ion").ne.0).or.(index(line,"ION").ne.0))&
         &    logions=.true.
         ! e+++++++++++++++++++++++++++++++++++
      case(19)
         read(line,*)igraph,ipoten,icon1,icon2
      case(20)
         read(line,*)imem
      case(21)
         read(line,*)phiwrt
      case(22)
         read(line,*)iacs,isen,isch
      end select
      end subroutine prm1 ! changed to F95 syntax on 2011-04-16
