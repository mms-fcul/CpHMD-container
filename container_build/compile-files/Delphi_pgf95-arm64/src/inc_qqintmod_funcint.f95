      !  ATTENTION!  This file is part of qqinttmod module.
      !==================================================================
      subroutine funcint(line,lineln)
      character(100) :: line                               
      character(100) :: calph= ' ', direc=' ', direc2=' ' 
      integer :: lineln,frmlen                             
      character(10) :: cnumb = '1234567890'                
      character(10) :: func(10),type=' ', frmnam           
      character(26) :: calph1="abcdefghijklmnopqrstuvwxyz" 
      character(26) :: calph2="ABCDEFGHIJKLMNOPQRSTUVWXYZ" 
      logical :: ij3=.false. ,ifrm=.false.                 
      logical :: ipr,iflag,ij4                             
      !-------------------------------------------------------------------
      ! 2011-06-14 Declarations added due to IMPLICIT NONE
      integer :: i,j,k,l,ieq1,ieq2,ipos,itype,j0,j1,j2
      integer :: k1,k2,k3,k4
      !-------------------------------------------------------------------
      line(lineln+1:100)=' '
      calph(1:40)='ABCDEFGHIJKLMNOPQRSTUVWXYZ.,:_-+=!@#$^123'
      calph(41:56)="4567890)(|\/?><;"
      do i=1,lineln
         if(line(i:i).eq."(") exit    
      end do
      ieq1=i-1; ieq2=lineln-ieq1;  type(1:ieq1)=line(1:ieq1)
      
      ! convert type to uppercase
      do i=1,ieq1
         j=index(calph1,type(i:i))
         if(j.ne.0) type(i:i)=calph2(j:j)
      end do
      direc(1:ieq2)=line(ieq1+1:lineln) 
      direc2(1:ieq2)=line(ieq1+1:lineln)
      
      !!c convert direc to uppercase, not direc2
      do i=1,ieq2
         j=index(calph1,direc(i:i))
         if(j.ne.0) direc(i:i)=calph2(j:j)
      end do
      ! changed 2011-04-15 to SELECT CASE statement from multiple IFs.
      select case(type(1:ieq1))
      case('CENTER')    ;   itype=1
      case('CENT')      ;   itype=1
      case('ACENTER')   ;   itype=2
      case('ACENT')     ;   itype=2
      case('IN')        ;   itype=3
      case('READ')      ;   itype=3
      case('OUT')       ;   itype=4
      case('WRITE')     ;   itype=4
      case('ENERGY')    ;   itype=5
      case('SITE')      ;   itype=6
      case('BUFFZ')     ;   itype=7
      case('QPREF')     ;   itype=8
      case('INSOBJ')    ;   itype=9
      case default      ;   itype=0
      end select
      !!c see if the function contains a format statement
      !!c use direc2 to maintain case sensitivity
      k1=index(direc(1:ieq2),"FRM=")
      k2=index(direc(1:ieq2),"FORM=")
      k3=index(direc(1:ieq2),"FORMAT=")
      if((k1+k2+k3).ne.0) then
         if(k1.ne.0) k=k1+4
         if(k2.ne.0) k=k2+5
         if(k3.ne.0) k=k3+7
         !    Changed 2011-04-15 to DO EXIT loop from GOTO statements      
         	       k4=k+1        
         do
            if(k4.gt.len_trim(direc)) exit
            if((direc(k4:k4).eq.",").or.(direc(k4:k4).eq.")")) exit
            k4=k4+1
         end do
         k4=k4-1
         if(index(calph,direc(k:k)).eq.0) then
            k=k+1; k4=k4-1
         end if
         frmlen=k4+1-k
         frmnam=direc2(k:k4)
         call up(frmnam,k4-k+1)
         ifrm=.true.
      end if
      
      j1=index(direc,"UNIT=")
      j2=index(direc,"FILE=")
      if((j1+j2).ne.0) ij3=.true.
      ! 2011-04-16 changed from multiple IFs to SELECT CASE  
      select case(itype)  
         
      case(1)              !  2011-04-16 changed from IF(itype.eq.0)
         
         j0=0
         if((direc(1:ieq2).eq."()").or.(direc(1:ieq2).eq."(0,0,0)")) then
            offset=coord(0.,0.,0.) ;  j0=1
         end if
         
         if(ij3) call rdflnm(j1,j2,direc2,centnam,centlen)
         
         if((j1+j2).ne.0) then
            offset=coord(999.,0.,0.)
         end if
         
         if(index(direc,"AN=1").ne.0) then
            offset=coord(999.,999.,0.)
         end if
         
         if((j0+j1+j2).eq.0) then
            k1=2
            o302: do                  
               k2=k1;  ipr=.false.
               do i=k2,ieq2
                  if(direc(i:i).eq.".") ipr=.true.
                  if((direc(i:i).eq.",").or.(direc(i:i).eq.")"))then
                     if(.not.ipr) then
                        do l=ieq2+1,i+1,-1
                           direc(l:l)=direc(l-1:l-1)
                           direc2(l:l)=direc2(l-1:l-1)
                        end do
                        direc(i:i)="." ; direc2(i:i)="."
                        ieq2=ieq2+1 ;    k1=i+2
                        if(k1.eq.ieq2) exit o302 
                        cycle o302               
                     else
                        ipr=.false. 
                     end if
                  end if
               end do
               exit o302
            end do o302
            ! 2011-07-19 Turned out formatted read here is wrong!
            read(direc(2:ieq2-1),*) offset 
         end if
         
      case(2)   ! 2011-04-16  changed from IF(itype.eq.2)
         
         k1=2
         o300: do
            k2=k1 ;  ipr=.false.
            do i=k2,ieq2
               if(direc(i:i).eq.".") ipr=.true.
               if((direc(i:i).eq.",").or.(direc(i:i).eq.")"))then
                  if(.not.ipr) then
                     do l=ieq2+1,i+1,-1
                        direc(l:l)=direc(l-1:l-1)
                        direc2(l:l)=direc2(l-1:l-1)
                     end do
                     direc(i:i)="." ; direc2(i:i)="."
                     ieq2=ieq2+1 ;    k1=i+2
                     if(k1.eq.ieq2) exit o300 
                     cycle o300               
                  else
                     ipr=.false.
                  end if
               end if
            end do
            exit o300
         end do o300
         ! 2011-07-19 Turned out formatted read here is wrong!
         read(direc(2:ieq2-1),*)acent
         iacent=.true.
      case(3)   !  2011-04-16  changed from IF(itype.eq.3)
         
         if(index(direc(1:ieq2),"(SIZ").ne.0) then
            sizfrm=0
            if(ij3) call rdflnm(j1,j2,direc2,siznam,sizlen)
            siznam(sizlen+1:len(siznam))=' '
            if(ifrm) then
               if(sizfrm.eq.0) then
                  write(6,*) "unknown size file format ",trim(frmnam)
               end if
            end if
         end if
         
         if(index(direc(1:ieq2),"(CRG").ne.0) then
            crgfrm=0
            if(ij3) call rdflnm(j1,j2,direc2,crgnam,crglen)
            if(ifrm) then
               if(crgfrm.eq.0) then
                  write(6,*) "unknown charge file format ",trim(frmnam)
               end if
            end if
         end if
         
         if(index(direc(1:ieq2),"(PDB").ne.0.or.index(direc(1:ieq2),"(MODPDB").ne.0) then
            pdbfrm=0
            if(ij3) call rdflnm(j1,j2,direc2,pdbnam,pdblen)
            if(ifrm) then
               if(frmnam(:frmlen).eq."UN") pdbfrm=1
               if(frmnam(:frmlen).eq."MOD") pdbfrm=2
               if(frmnam(:frmlen).eq."PQR") pdbfrm=3
               if(pdbfrm.gt.0) ipdbrd=.true.
               if(pdbfrm.eq.0) then
                  write(6,*) "unknown pdb file format",trim(frmnam)
               end if
            end if
         end if
         
         !##############Lin Li:modpdb4 option:
         if(index(direc(1:ieq2),"(MODPDB4").ne.0) then
            pdbfrm=40
            if(ij3) call rdflnm(j1,j2,direc2,pdbnam,pdblen)
            if(ifrm) then
               if(frmnam(:frmlen).eq."UN") pdbfrm=41
               if(frmnam(:frmlen).eq."MOD") pdbfrm=42
               if(frmnam(:frmlen).eq."PQR") pdbfrm=43
               if(pdbfrm.gt.40) ipdbrd=.true.
               if(pdbfrm.eq.40) then
                  write(6,*) "unknown pdb file format",trim(frmnam)
               end if
            end if
         end if
         
         if(index(direc(1:ieq2),"(FRC").ne.0) then
            frcifrm=0
            if(ij3) call rdflnm(j1,j2,direc2,frcinam,frcilen)
            if(frcinam(:frclen).eq."self") iself=.true.
            if(frcinam(:frclen).eq."SELF") iself=.true.
            if(ifrm) then
               if(frcifrm.eq.0) then
                  write(6,*) "unknown input frc file format for ",&
                  &trim(frmnam)
               end if
            end if
         end if
         
         if(index(direc(1:ieq2),"(PHI").ne.0) then
            phiifrm=0
            if(ij3) call rdflnm(j1,j2,direc2,phiinam,phiilen)
            if(ifrm) then
               if(phiifrm.eq.0) then
                  write(6,*) "unknown input phimap format",trim(frmnam)
               end if
            end if
         end if
      case(4)     !   2011-04-16   changed from IF(itype.eq.4)
         iflag=.true.
         if(index(direc(1:ieq2),"OFF)").ne.0) iflag=.false.
         ij4=iflag.and.ij3
         
         if(index(direc(1:ieq2),"(PHI").ne.0) then
            phiwrt=iflag ;  phifrm=0
            if(ij4) call rdflnm(j1,j2,direc2,phinam,philen)
            if(ifrm) then
               if(frmnam(:frmlen).eq."BIOSYM") phifrm=1
               if(frmnam(:frmlen).eq."GRASP")  phifrm=2
               ! b+++++++++++++++++++++++++++++++++++++++++++
               if(frmnam(:frmlen).eq."CCP4")   phifrm=3
               if(frmnam(:frmlen).eq."DIFF")   phifrm=4
               ! 2011-07-23 New piece of code from Maxim
               if(frmnam(:frmlen).eq."CUBE")   phifrm=5
               ! e+++++++++++++++++++++++++++++++++++++++++++
               if(phifrm.eq.1) ibios=.true.
               if(phifrm.eq.0) then
                  write(6,*) "unknown phimap format: ",trim(frmnam)
               end if
            end if
         end if
         
         if(index(direc(1:ieq2),"(SRF").ne.0) then
            isrf=.true. ;  ibem=.false.
            if(ij4) call rdflnm(j1,j2,direc2,srfnam,srflen)
            if(ifrm) then
               if(frmnam(:frmlen).eq."BEM")ibem=.true.
            end if
         end if
         
         if(index(direc(1:ieq2),"(FRC").ne.0) then
            isite=iflag ;  frcfrm=0
            if(ij4) call rdflnm(j1,j2,direc2,frcnam,frclen)
            if(ifrm) then
               if(frmnam(:frmlen).eq."RC") frcfrm=1
               if(frmnam(:frmlen).eq."R") frcfrm=2
               if(frmnam(:frmlen).eq."UN") frcfrm=3
               if((frcfrm.eq.1).or.(frcfrm.eq.2)) irea=.true.
               if((frcfrm.gt.3).or.(frcfrm.lt.1)) then
                  write(6,*) "unknown frc format: ",trim(frmnam)
               end if
            end if
         end if
         
         if(index(direc(1:ieq2),"(EPS").ne.0) then
            epswrt=iflag ;  epsfrm=0
            if(ij4) call rdflnm(j1,j2,direc2,epsnam,epslen)
            if(ifrm) then
               if(epsfrm.eq.0) then
                  write(6,*) "unknown eps format: ",trim(frmnam)
               end if
            end if
         end if
         
         if(index(direc(1:ieq2),"(MODPDB").ne.0) then
            iatout=iflag ;   mpdbfrm=0
            if(ij4) call rdflnm(j1,j2,direc2,mpdbnam,mpdblen)
            if(ifrm) then
               ! ++++++++++WWW+++++++++++
               if(frmnam(:frmlen).eq."PQR") mpdbfrm=1
               ! ++++++++++WWW+++++++++++
               if(mpdbfrm.eq.0) then
                  write(6,*) "unknown modified pdb format: ",&
                  & trim(frmnam)
               end if
            end if
         end if
         
         !#################Lin Li: for 4 digits precession:
         
         if(index(direc(1:ieq2),"(MODPDB4").ne.0) then
            iatout=iflag ;   mpdbfrm=40
            if(ij4) call rdflnm(j1,j2,direc2,mpdbnam,mpdblen)
            if(ifrm) then
               ! ++++++++++WWW+++++++++++
               if(frmnam(:frmlen).eq."PQR") mpdbfrm=41
               ! ++++++++++WWW+++++++++++
               if(mpdbfrm.eq.40) then
                  write(6,*) "unknown modified pdb format: ",&
                  & trim(frmnam)
               end if
            end if
         end if
         
         if(index(direc(1:ieq2),"(UNPDB").ne.0) then
            ipdbwrt=iflag ;  updbfrm=0
            if(ij4) call rdflnm(j1,j2,direc2,updbnam,updblen)
            if(ifrm) then
               if(updbfrm.eq.0) then
                  write(6,*) "unknown unformatted pdb format: ",&
                  &trim(frmnam)
               end if
            end if
         end if
         
         if(index(direc(1:ieq2),"(UNFRC").ne.0) then
            ifrcwrt=iflag ;  ufrcfrm=0
            if(ij4) call rdflnm(j1,j2,direc2,ufrcnam,ufrclen)
            if(ifrm) then
               if(ufrcfrm.eq.0) then
                  write(6,*) "unknown unformatted frc format: ",&
                  &trim(frmnam)
               end if
            end if
         end if
         
         if(index(direc(1:ieq2),"(ENERGY").ne.0) then
            inrgwrt=iflag ;  nrgfrm=0
            if(ij4) call rdflnm(j1,j2,direc2,nrgnam,nrglen)
            if(ifrm) then
               if(nrgfrm.eq.0) then
                  write(6,*) "unknown energy file format: ",&
                  &trim(frmnam)
               end if
            end if
         end if
         
         if(index(direc(1:ieq2),"(GCRG").ne.0)   iwgcrg=.true.
         if(index(direc(1:ieq2),"(HSURF").ne.0)  iacs=.true.
         if(index(direc(1:ieq2),"(DB").ne.0)     idbwrt=.true.
         if(index(direc(1:ieq2),"(SURFEN").ne.0) isen=.true.
         if(index(direc(1:ieq2),"(SCRG").ne.0) then
            isch=iflag ;  scrgfrm=0
            if(ij4) call rdflnm(j1,j2,direc2,scrgnam,scrglen)
            if(ifrm) then
               if(frmnam(:frmlen).eq."PDB") scrgfrm=1
               if(frmnam(:frmlen).eq."PDBA")&
               & write(6,*)"option PDBA no longer supported!"
               if((scrgfrm.gt.2).or.(scrgfrm.lt.1)) &
               & write(6,*)"unknown surface charge format: ",&
               &scrgnam(:scrglen)
            end if
         end if
         
      case(5)  !  changed 2011-04-16  from IF(itype.eq.5) 
         
         direc(ieq2:ieq2)=","
         if(index(direc(1:ieq2),"GRID,").ne.0) logg=.true.
         
         if(index(direc(1:ieq2),"S,").ne.0)then
            ipos=index(direc(1:ieq2),"S,")
            if(direc(ipos-1:ipos-1).ne."A")then
               logs=.true.
            else
               logas=.true.
            end if
         end if
         
         if(index(direc(1:ieq2),"G,").ne.0)then
            ipos=index(direc(1:ieq2),"G,")
            if(direc(ipos-1:ipos-1).ne."A")then
               logg=.true.
            else
               loga=.true.
            end if
         end if
         
         !!c b+++++++++++++++++++++++++++++++++++++w 2000
         if(index(direc(1:ieq2),"ION,").ne.0) logions=.true.
         if(index(direc(1:ieq2),"IONIC_C").ne.0) logions=.true.
         !!c e+++++++++++++++++++++++++++++++++++++++++++
         if(index(direc(1:ieq2),"SOLVATION").ne.0) logs=.true.
         if(index(direc(1:ieq2),"SOL,").ne.0) logs=.true.
         if(index(direc(1:ieq2),"C,").ne.0) logc=.true.
         if(index(direc(1:ieq2),"COU,").ne.0) logc=.true.
         if(index(direc(1:ieq2),"COULOMBIC").ne.0) logc=.true.
         if(index(direc(1:ieq2),"AS,").ne.0) logas=.true.
         if(index(direc(1:ieq2),"ANASURF,").ne.0) logas=.true.
         if(index(direc(1:ieq2),"ANALYTICSURFACE").ne.0) logas=.true.
         if(index(direc(1:ieq2),"AG,").ne.0) loga=.true.
         if(index(direc(1:ieq2),"ANAGRID,").ne.0) loga=.true.
         if(index(direc(1:ieq2),"ANALYTICGRID").ne.0) loga=.true.
         
      case(6)  ! 2011-04-16   changed from IF(itype.eq.6)
         
         !!c empty bracket means reset all quantities to false
         if(ieq2.eq.2) then
            isita=.false.; isitx=.false.; isitq=.false.; isitp=.false.
            !!c b+++++++++++++++++++++++
            isitap=.false.; isitmd=.false.; isittf=.false.
            isitsf=.false.; isitpot=.false.; isitdeb=.false.
            !!c e+++++++++++++++++++++++
            isiti=.false.; isitr=.false.; isitc=.false.
            isitf=.false.; isitt=.false.
         end if
         
         direc(ieq2:ieq2)=","
         if(index(direc(1:ieq2),"ATOM,").ne.0) isita=.true.
         if(index(direc(1:ieq2),"CHARGE,").ne.0) isitq=.true.
         if(index(direc(1:ieq2),"POTENTIAL,").ne.0) isitp=.true.
         !!c b+++++++++++++++++++++++++++++
         if(index(direc(1:ieq2),"ATOMICPOT,").ne.0) isitap=.true.
         if(index(direc(1:ieq2),"DEBYEFRACTION,").ne.0) isitdeb=.true.
         !!c e+++++++++++++++++++++++++++++
         if(index(direc(1:ieq2),"FIELD,").ne.0) isitf=.true.
         if(index(direc(1:ieq2),"REACTION,").ne.0) isitr=.true.
         if(index(direc(1:ieq2),"COULOMB,").ne.0) isitc=.true.
         if(index(direc(1:ieq2),"COORDINATES,").ne.0) isitx=.true.
         if(index(direc(1:ieq2),"SALT,").ne.0) isiti=.true.
         if(index(direc(1:ieq2),"TOTAL,").ne.0) isitt=.true.
         
         !!c NB can only do last letters like this if there is no cooincindence with
         !!c those above
         if(index(direc(1:ieq2),"A,").ne.0) isita=.true.
         if(index(direc(1:ieq2),"Q,").ne.0) isitq=.true.
         if(index(direc(1:ieq2),"P,").ne.0) isitp=.true.
         if(index(direc(1:ieq2),"R,").ne.0) isitr=.true.
         if(index(direc(1:ieq2),"C,").ne.0) isitc=.true.
         if(index(direc(1:ieq2),"X,").ne.0) isitx=.true.
         if(index(direc(1:ieq2),"I,").ne.0) isiti=.true.
         if(index(direc(1:ieq2),"T,").ne.0) isitt=.true.
         
         !!c b++++++++++++++++++
         if(index(direc(1:ieq2),"ATPO,").ne.0) isitap=.true.
         if(index(direc(1:ieq2),"DEB,").ne.0) isitdeb=.true.
         !!c e++++++++++++++++++
         !!c extra ones for fields
         if(index(direc(1:ieq2),"F,").ne.0) then
            if(index(direc(1:ieq2),"RF,").ne.0) isitrf=.true.
            if(index(direc(1:ieq2),"CF,").ne.0) isitcf=.true.
            !!c b++++++++++++++++++
            if(index(direc(1:ieq2),"MDF").ne.0) isitmd=.true.
            if(index(direc(1:ieq2),"SF").ne.0) isitsf=.true.
            if(index(direc(1:ieq2),"TF").ne.0) isittf=.true.
            if(.not.(isitrf.or.isitcf.or.isitmd.or.isittf.or.isitsf))&
            &  isitf=.true.
            !!c	  if((.not.isitrf).and.(.not.isitcf).and.(.not.isitmd)&
            !!c                &.and.(.not.isittf)) isitf=.true.
            !!c e++++++++++++++++++
         end if
         
         if(isitr) irea=.true.
         if(isitrf) irea=.true.
         !!c b++++++++++++++++++
         if(isitmd.or.isittf.or.isitsf) irea=.true.
         if(isitsf) iself=.true.
         !!c e++++++++++++++++++
         if(isitt) irea=.true.
         
      case(7)          !  2011-04-16    changed from IF(itype.eq.7)
         
         read(direc(2:ieq2-1),'(6i3)') bufz
         !!c 333            format(6i3)
         ibufz=.true.
         
      case(8)            !  2011-04-16    changed from IF(itype.eq.8)
         
         !!c b++++++++++++++++++++++++++++++++++++++++++++++++++++
      case(9)            !  2011-04-16    changed from IF(itype.eq.9)
         icreapdb=.true.
      case default   !    2011-04-16 changed from if (itype.eq.0)
         write(6,*) "!!!!!!!!!!!!!!!!!"
         write(6,*) "The function specifier ",type(1:ieq1)," is"
         write(6,*) "not recognised. Therefore the function will ",&
         &  "not be processed"
         write(6,*) "!!!!!!!!!!!!!!!!!"
      end select
      !!c e++++++++++++++++++++++++++++++++++++++++++++++++++
      end subroutine funcint ! changed 2011-04-15 to F95 syntax