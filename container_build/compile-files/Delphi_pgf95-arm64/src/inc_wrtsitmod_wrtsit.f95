      !  ATTENTION!  This file is part of wrtsitmod module.
      !==================================================================
      ! 2011-06-11 All other parameters are accessible via qlog and pointers
      subroutine wrtsit(iisitsf)
      
      !!c write a file containing site potentials and/or fields and/or atom
      !!c information. adapted april 92 to enable flexible output/input
      !!c possible fields include atom information, coordinates, charge, potential
      !!c salt concentration, reaction potentials, coulombic potential, and fields
      !!c the default of which is coordinates, charges, potentials and fields.
      
      !!c nqass = number of assigned charges, icount2b= number of boundary elements
      !!c scspos= position in angstroms of induced surface charges
      !!c xn1= positions of all atoms in angstroms, natm=number of atoms
      
      ! 2011-06-12 Arrays and variables below are accessible via
      !!c b+++++++++++++++++++++++
      character(4) :: resnum
      logical :: isitmp,isitmp1,residsf(resnummax),atmsf(natom*iisitsf)
      character(80) :: filnam
      character(65) :: datum
      character(80) :: line,vrow,oline
      character(15) :: otemp,f200
      character(16) :: atdes
      character(6) :: head
      character(24) :: crdstr
      character(5) :: atnum
      type(coord) :: xo,xn,xu,trgelf,xu2,xo2,xt,exyz,cxyz
      type(coord) :: fxyz, fl,fu,vtemp,rxyz
      type(object_min_max) :: bedge
      real :: scomp(30),sold(30)
      logical :: ifrm,ifrm2,iext,iqass,ofrm
      integer :: noradu,iresnum,iisitsf,jtmp,ncrgs,ncrg
      integer :: i,j,k,inum,nnatom,ios,ii,jj,idist,it
      real :: aphi,radu,rads,debyefraction,sixth,goff,phirt,phict
      real :: chrgv,crgs,eps,etot,ff,finish,fn,dist,phiac,phiat
      real :: phic,phias,phir,phit,phirtt,qphiv,pot,phiv,phii
      real :: tcrgs,vphi,temp,sdist,phiact
      !---------------------------------------------------------------------------
      f200='(f10.4)'
      goff=(igrid+1.)/2.
      residsf=.false.
      !!c e+++++++++++++++++++++++
      !!c      call leggigreen
      allocate(rfield(ibnum))
      !!c initialize some parameters
      atdes=" "; phirt=0; phict=0; sixth=1.0/6.0
      iqass=.true.; ofrm=.true.
      !  2011-06-12  Using operations on coord type variables defined
      bedge%min=oldmid-(0.5*(igrid-1)/scale)
      bedge%max=oldmid+(0.5*(igrid-1)/scale)
      
      !!c ifrm is true if ANY of the flags for special output have been set
      ifrm=isita.or.isitq.or.isitp.or.isitf.or.isitr.or.isitt
      ifrm=ifrm.or.isitc.or.isitx.or.isiti.or.isitrf.or.isitcf
      !!c b+++++isitap=site atomic potential ++++ isitdeb= site debyemap fraction
      ifrm=ifrm.or.isitap.or.isittf.or.isitdeb
      !!c +++++++isitsf=site surface charge and electric field at site's 
      !!          Solvent Accessible Surface
      isitmp1=ifrm
      ifrm=ifrm.or.isitsf
      !!c +++++++isitmd=site reaction and coulombic fields for molecular dynamics
      if (isitmd.or.isitpot) then
         if (isitmd) allocate(atmforce(natom))
         frcfrm=-1 ; ofrm=.false.
      end if
      !!c e++++++++++++++++++++++++++++++++++++++++
      if((.not.ifrm).and.((frcfrm.eq.0).or.(frcfrm.eq.3))) then
         !!c set default to standard frc file
         isitx=.true.; isitq=.true.; isitf=.true.; isitp=.true.
      end if
      select case(frcfrm)
      case(1)               !   if(frcfrm.eq.1) then
         isitx=.true.; isitq=.true.; isitr=.true.; isitc=.true.
      case(2)               !   if(frcfrm.eq.2) then
         isitx=.true.; isitq=.true.;  isitr=.true.
      case(3)               !   if(frcfrm.eq.3) then
         ofrm=.false.
      end select
      
      vrow=' ';  datum=' ';  j=1;  k=1
      if(isita) then
         datum(k:k+4)="ATOM "
         vrow(j:j+14)="ATOM DESCRIPTOR"
         j=j+20; k=k+5
      end if
      if(isitx) then
         vrow(j+4:j+27)="ATOM COORDINATES (X,Y,Z)"
         datum(k:k+11)="COORDINATES "
         k=k+12; j=j+30
      end if
      if(isitq) then
         vrow(j+3:j+8)="CHARGE"
         datum(k:k+6)="CHARGE "
         k=k+7; j=j+10
      end if
      if(isitp) then
         vrow(j+2:j+9)="GRID PT."
         datum(k:k+10)="POTENTIALS "
         k=k+11; j=j+10
      end if
      if(isiti) then
         vrow(j+1:j+8)="SALT CON"
         datum(k:k+4)="SALT "
         j=j+10; k=k+5
      end if
      if(j.gt.80) then
         isitr=.false.;  isitc=.false.
         !!c b+++++++++++++++++++++++++++++++++++++++
         isitap=.false.; isitdeb=.false.; isittf=.false.; isitsf=.false.
         !!c e+++++++++++++++++++++++++++++++++++++++
         isitf=.false.; isitrf=.false.; isitcf=.false.; isitt=.false.
      end if
      if(isitr) then
         vrow(j:j+9)=" REAC. PT."
         datum(k:k+8)="REACTION "
         k=k+9; j=j+10
      end if
      if(j.gt.80) then
         isitc=.false.
         !!c b+++++++++++++++++++++++++++++++++++++++
         isitap=.false.; isitdeb=.false.; isittf=.false.; isitsf=.false.
         !!c e+++++++++++++++++++++++++++++++++++++++
         isitf=.false.; isitrf=.false.; isitcf=.false.; isitt=.false.
      end if
      if(isitc) then
         vrow(j:j+9)=" COUL. POT"
         datum(k:k+9)="COULOMBIC "
         k=k+10;  j=j+10
      end if
      !!c b+++++++++++++++++++++++++++++++++++++++
      if(j.gt.80) then
         isitap=.false.;  isitdeb=.false.;  isitf=.false.
         isitrf=.false.;  isitcf=.false.;   isittf=.false.
         isitt=.false.;   isitsf=.false.
      end if
      if(isitap) then
         vrow(j+2:j+9)="ATOM PT."
         datum(k:k+10)="ATOMIC PT. "
         k=k+11;  j=j+10
      end if
      if(j.gt.80) then
         isitdeb=.false.; isitf=.false.;  isitrf=.false.
         isitcf=.false.;  isittf=.false.;  isitt=.false.
         isitsf=.false.
      end if
      if(isitdeb) then
         vrow(j+3:j+13)="DEBFRACTION"
         datum(k:k+11)="DEBFRACTION "
         k=k+12; j=j+14
      end if
      !!c e+++++++++++++++++++++++++++++++++++++++
      if(j.gt.60) then
         isitf=.false.;  isitrf=.false.; isitcf=.false.
         !!c b+++++++++++++++++++++++++++++++++++++++
         isittf=.false.;  isitsf=.false.
         !!c e+++++++++++++++++++++++++++++++++++++++
         isitt=.false.
      end if
      
      if(isitf) then
         vrow(j+4:j+28)="GRID FIELDS: (Ex, Ey, Ez)"
         datum(k:k+5)="FIELDS "
         j=j+30
      end if
      if(j.gt.60) then
         isitrf=.false.; isitcf=.false.
         !!c b+++++++++++++++++++++++++++++++++++++++
         isittf=.false.; isitsf=.false.
         !!c e+++++++++++++++++++++++++++++++++++++++
         isitt=.false.
      end if
      
      if(isitrf) then
         vrow(j+4:j+28)="REAC. FORCE: (Rx, Ry, Rz)"
         datum(k:k+5)="RFORCE "
         j=j+30
      end if
      if(j.gt.60) then
         isitcf=.false.
         !!c b+++++++++++++++++++++++++++++++++++++++
         isittf=.false.; isitsf=.false.
         !!c e+++++++++++++++++++++++++++++++++++++++
         isitt=.false.
      end if
      
      if(isitcf) then
         vrow(j+4:j+28)="COUL. FORCE: (Cx, Cy, Cz)"
         datum(k:k+5)="CFORCE "
         j=j+30
      end if
      if(j.gt.60) then
         !!c b+++++++++++++++++++++++++++++++++++++++
         isittf=.false.
         isitsf=.false.
         !!c e+++++++++++++++++++++++++++++++++++++++
         isitt=.false.
      end if
      
      if(isittf) then
         vrow(j+4:j+28)="TOTAL FORCE: (Tx, Ty, Tz)"
         datum(k:k+5)="TFORCE "
         j=j+30
      end if
      if(j.gt.70) then
         isitt=.false.
      end if
      
      if(isitt) then
         vrow(j+4:j+9)=" TOTAL"
         datum(k:k+5)="TOTAL "
         j=j+10
      end if
      
      !!c b+++++++++++++++++++++++++++++++++++++++
      if(j.gt.50) isitsf=.false.
      if(isitsf) then
         vrow(j+4:j+68)="sCharge,    x          y       z       &
         &surf.E°n,surf. E[kT/(qA)]"
         datum(k:k+34)="SCh, x, y, z, surf En, surf. E"
         j=j+50
      end if
      !!c e+++++++++++++++++++++++++++++++++++++++
      
      !!c if site potentials required and unformatted read/write, skip
      !c during formatted frc file read/write can write unformatted frc.pdb
      pot=0.0 ; exyz=coord(0.,0.,0.); cxyz=coord(0.,0.,0.)
      if (.not.(isitmd.or.isitpot)) then
         !!c open files, write to log file..
         write(6,*)'  '
         write(6,*)'writing potentials at given sites...'
         write(6,*)'  '
      end if
      
      if(iself) then
         write(6,*) "using the current pdb file"
         ifrm2=.true.;  iqass=.false.
      else
         if (.not.isitpot) then
            inquire(file=trim(frcinam),exist=iext)
            if(.not.iext) then
               write(6,*) "the input frc file ",trim(frcinam),&
               & " does      not exist"
               write(6,*) "exiting..."
            else
               write(6,*)'coordinates, etc for potential &
               &output read from file ',trim(frcinam)
               write(6,*)'  '
            end if
            ! 2011-06-13 Argument frcilen is obsolete, using trim function instead                 
            call form(frcinam,frcilen,ifrm2)
         end if
      end if
      
      !!c if unformatted may not contain all the info needed for all options, i.e
      !!c atom info
      if((.not.ifrm2).and.(isita)) then
         write(6,*) "atom info flag turned off cos this &
         &unformatted file does"
         write(6,*) "not contain atom info"
         isita=.false.; iqass=.false.
      end if
      
      if(.not.ifrm2) iqass=.false.
      
      if(.not.iself.and..not.isitpot) then
         if(ifrm2) then
            open(15,file=trim(frcinam),iostat=ios)
         else
            open(15,file=trim(frcinam), &
            & form="unformatted",iostat=ios)
         end if
         if(ios.ne.0) then
            write(6,*) "error reading the frc input file ",&
            &trim(frcinam)
            if(allocated(rfield)) deallocate(rfield)
            return
         end if
      end if
      
      !!c b+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(isitsf) then
         !!c isitsf assumes ifrm2=.true.
         inquire(file=frcinam(:frcilen),exist=iext)
         if(.not.iext) then
            write(6,*) "the input frc file ",frcinam(:frcilen),&
            &" does not exist"
            write(6,*) "exiting..."
         else
            write(6,*)'coordinates, etc for potential &
            &output read from file'
            write(6,*)frcinam(:frcilen)
            write(6,*)'  '
         end if
         open(15,file=frcinam(:frcilen),iostat=ios)
         if(ios.ne.0) then
            write(6,*) "error reading the frc input file ",&
            &trim(frcinam)
            if(allocated(rfield)) deallocate(rfield)
            return
         end if
         
         D302: do
            read(15,'(a)',iostat=ios)line
            if(ios.ne.0) exit D302
            head = line(1:6)
            if(head.ne.' ') call up(head,6)
            if((head.ne.'ATOM  ').and.(head.ne.'HETATM')) cycle D302
            resnum= line(24:27)
            read(resnum,'(i4)')iresnum
            residsf(iresnum)=.true.
         end do D302
         close(15)
      end if
      !!c e++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(ofrm) open(16,file=trim(frcnam))
      if(.not.(ofrm.or.isitmd.or.isitpot)) then
         open(16,file=trim(frcnam),form="unformatted")
      end if
      if (.not.isitmd.and..not.isitpot) then 
         write(6,*)'potentials written to file ',trim(frcnam)
         write(6,*)'  '
      end if
      
      !!c write header to file...(10 lines)
      
      if(ofrm) then
         write(16,*)'DELPHI SITE POTENTIAL FILE'
         write(16,*)'grid size,percent fill:',igrid,perfil
         write(16,*)'outer diel. and first one assigned :',&
         &repsout, medeps(1)*epkt
         !!c b+++++++++++++++ mi piacerebbe ma non compatibile con vecchi scripts
         !!c         do i = 1,nmedia
         !!c           write(16,*)'dielectric in medium nr. ',i,':',medeps(i)*epkt
         !!c         end do
         !!c e+++++++++++++++
         write(16,*)'ionic strength (M):',rionst
         write(16,*)'ion excl., probe radii:',exrad,radprb(1),radprb(2)
         write(16,*)'linear, nolinear iterations:',nlit,nnit
         write(16,*)'boundary condition:',ibctyp
         write(16,*)'Data Output: ',datum
         write(16,*)'title: ',toplbl
         write(16,*)'    '
         write(16,*)'    '
         write(16,'(a)')vrow
      end if
      
      if((.not.ofrm).and.(.not.(isitmd.or.isitpot))) then
         line=" "; line="DELPHI FRC FILE"; write(16) line
         line=" "; line="FORMAT NUMBER=1"; write(16) line
         line=" "; line(1:5)="DATA="; line(6:70)=datum
         write(16) line
         line=" "; write(line,'(i5,3f10.4)') igrid,perfil,repsout,rionst
         write(16) line
         !!c b+++++++++++++++
         do i = 1,nmedia
            write(6,*)'dielectric in medium nr. ',i,':',medeps(i)*epkt
         end do
         !!c e+++++++++++++++
         line=" "
         write(line,'(3f10.4,3i5)') exrad,radprb,nlit,nnit,ibctyp
         write(16) line
      end if
      
      !!c write 2 formatting lines
      !!c read atom coordinate file
      etot = 0.;  nnatom = 0; chrgv=0.0
      
      if((.not.iself).and.(isitrf.or.isitmd.or.isittf)) then
         write(6,*)"WARNING! Cannot calculate reaction forces without"
         write(6,*) "WARNING! Using internal (self) coordinates"
         isitrf=.false.; isittf=.false.; isitmd=.false.
      end if
      
      if(isitrf.or.isitmd.or.isittf) then
         if (nmedia.eq.1.and.abs(medeps(1)*epkt-1.).lt.1.e-6)then
            !!c b+++++++++++++++++++++++++++++++++++++++
            !!c	    call rforce(rfield,natm,scspos,schrg,atsurf,
            !!c     &    icount2b,atmcrg,xn1,chgpos,nqass,crgatn,nobject)
            ! 2011-06-12 All other parameters are accessible via qlog and pointers 
            call rforceeps1(rfield,schrg)
            !!c e+++++++++++++++++++++++++++++++++++++++
         else
            call rforce(rfield,schrg)
         end if
      end if
      !!c b+++++++++++++++++++++++++++++++++++++++
      if (isitpot) then
         do ii=1,npotenziali
            ! 2011-06-12 Something unclear with sitephi array (see react sub)
            xo=coord(sitephi(1,ii),sitephi(2,ii),sitephi(3,ii))
            xn=((xo-oldmid)*scale)+goff
            call phintp(xn,sitephi(4,ii))
         end do
      else
         !!c b+++++++++++++++++++++++++++++++++++++++
         D304: do
            !!c beginning of the big loop on natom
            if(iself) then
               if(nnatom.eq.natom) exit D304
               nnatom=nnatom+1
               xo=xn1(nnatom)
               chrgv=delphipdb(nnatom)%chrgv4
               !!c b+++++++++++++++++++++++++++++++++++++++
               radu=delphipdb(nnatom)%rad3*scale
               !!c e+++++++++++++++++++++++++++++++++++++++
               atm = delphipdb(nnatom)%atinf(1:4)
               res = delphipdb(nnatom)%atinf(7:9)
               rnum = delphipdb(nnatom)%atinf(12:15)
               chn = delphipdb(nnatom)%atinf(11:11)
            else
               if(ifrm2) then
                  read(15,'(a)',iostat=ios)line
                  if(ios.ne.0) exit D304
                  head = line(1:6)
                  call up(head,6)
                  if((head.ne.'ATOM  ').and.(head.ne.'HETATM')) cycle D304
                  nnatom = nnatom + 1
                  crdstr = line(31:54); atnum  = line(7:11)
                  !!c positions, atom number 
                  !!c non ho l'informazione concernente il raggio
                  read(crdstr,'(3f8.3)')xo
                  read(atnum,'(i5)')inum
               else
                  read(15,iostat=ios) xo,radu,chrgv
                  if(ios.ne.0) exit D304
                  nnatom=nnatom+1
               end if
            end if
            !!c end of atom read..
            
            isitmp=(isitq.and.iqass).or.isitap.or.isitp
            if((isita.or.isitmp).and..not.iself) then
               atm = line(12:16);   res = line(18:20)
               rnum = line(23:26);  chn = line(22:22)
               if(atm.ne.' ') then
                  call up(atm,6);     call elb(atm,6)
               end if
               if(res.ne.' ') then
                  call up(res,3);     call elb(res,3)
               end if
               if(rnum.ne.' ') then
                  call up(rnum,4);    call elb(rnum,4)
               end if
               if(atm.ne.' ') then
                  call up(chn,1);     call elb(chn,1)
               end if
            end if
            !!c scale atoms to grid space
            xn=((xo-oldmid)*scale)+goff
            if(isita) then
               atdes(1:4)=atm(1:4);    atdes(6:8)=res(1:3)
               atdes(12:15)=rnum(1:4); atdes(10:10)=chn(1:1)
            end if
            !!c assignchargeto atom, searching for decreasingly specific specification
            !!c note if no charge record found, is assumed to be 0.0
            
            !!c b+++++++++++++++++++++++++++++++++++++++++++++++++
            if((.not.iself.and.ifrm2).and.isitmp)then
               chrgv=0.0
               call ass(chrgv,charge,chash,Ncmax)
               if(isitap) call ass(radu,radii,rhash,Nrmax)
               radu=radu*scale
            end if
            
            if(isitsf) then
               read(rnum,'(i4)')iresnum
               atmsf(nnatom)=.false.
               if(residsf(iresnum)) atmsf(nnatom)=.true.
            end if
            !!c e+++++++++++++++++++++++++++++++++++++++++++++++++
            if(isitap.and.abs(chrgv).ge.1.e-6) then
               rads=min(radu,atompotdist*scale)
               xt=xn ; xt%x=xt%x+rads
               call phintp(xt,vphi)
               aphi=vphi
               xt=xn ; xt%x=xt%x-rads
               call phintp(xt,vphi)
               aphi=aphi+vphi
               xt=xn ; xt%y=xt%y+rads
               call phintp(xt,vphi)
               aphi=aphi+vphi
               xt=xn ; xt%y=xt%y-rads
               call phintp(xt,vphi)
               aphi=aphi+vphi
               xt=xn ; xt%z=xt%z+rads
               call phintp(xt,vphi)
               aphi=aphi+vphi
               xt=xn ; xt%z=xt%z-rads
               call phintp(xt,vphi)
               aphi=(aphi+vphi)/6.
            end if
            
            if(isitp.or.isiti.or.(isitap.and.abs(chrgv).lt.1.e-6))then
               call phintp(xn,vphi)
               if(isitap.and.abs(chrgv).lt.1.e-6) aphi=vphi
               if(isitp) then
                  qphiv=chrgv*vphi;  etot=etot+qphiv; phiv=vphi
               end if
               
               if(isiti) then
                  !!c NB we have changed the iconc action so that the phimap has NOT been
                  !!c converted to salt concentrations. therefore 
                  write(6,*)'WRTSIT:these salt concentrations'
                  write(6,*)'do NOT have the benefit of idebmap (as yet)'
                  if(nnit.ne.0) then
                     !!c b+++++++++++++++++++++
                     temp = vphi*chi5+chi4; temp = vphi*temp+chi3
                     temp = vphi*temp+chi2; temp = chi1+temp*vphi
                     phii = vphi*temp
                  else
                     phii= -rionst*2.0*vphi
                     !!c e++++++++++to check+++++++++++++
                  end if
               end if
               !!c end if isitp or isiti, salt and or potentials
            end if
            !!c b+++++++++++++++++++++
            if(isitdeb) then
               !!c       it calculates the fraction of closest grid points that are in solution
               write(6,*)'Calculating Debye Fraction'
               call debtp(xn,debyefraction)
            end if
            !!c e+++++++++++++++++++++++
            
            if(isitf) then
               xn%x=xn%x+1.          ! xn(1) = xn(1) + 1.
               call phintp(xn,fu%x)
               xn%x=xn%x-2.          ! xn(1) = xn(1) - 2.
               call phintp(xn,fl%x)
               xn%x=xn%x+1.; xn%y=xn%y+1.
               call phintp(xn,fu%y)
               xn%y=xn%y-2.          ! xn(2) = xn(2) - 2.
               call phintp(xn,fl%y)
               xn%y=xn%y+1. ; xn%z=xn%z+1.
               call phintp(xn,fu%z)
               xn%z=xn%z-2.          ! xn(3) = xn(3) - 2.
               call phintp(xn,fl%z)
               xn%z=xn%z+1.          ! xn(3) = xn(3) + 1.
               !!c b+++++++++++Walter OCt 2000
               !!c the electric field is opposite the potential gradient
               !!c so I change the sign
               fxyz=(fl-fu)*(0.5*scale)
               !!c e+++++++++++++++++++++++++++
            end if
            
            if(isitt) then
               !!c check if this point is within the box.
               it=0
               if((xo.vorlt.bedge%min).or.(xo.vorgt.bedge%max)) it=1
               
               if(it.eq.0) then
                  xo2=((xo-oldmid)*scale)+goff
                  !!c first find reaction field from surface elements inside of the box..
                  phir=0.0 ; phias=0.0; ncrgs=0; tcrgs=0.0; sold=0.
                  
                  do i=1,icount2b
                     vtemp=xo-scspos(i); dist=sqrt(vtemp.dot.vtemp)
                     !!c find analytic potential from this induced charge..=phias
                     ncrgs=ncrgs+1; tcrgs=tcrgs+schrg(i)
                     phirtt=schrg(i)/dist
                     !!c medeps either epsin contain the 561.0 factor....
                     phirtt=phirtt*epkt ; phir=phir+phirtt
                     xu2=float(ibgrd(i))
                     crgs=schrg(i)
                     !!c ++++++++++1 took place of repsin because eps is no more included 
                     !!c       in schrg , surface charge
                     call tops(xu2,xo2,crgs,1.,scale,phiat,trgelf,1)
                     phiat=phiat*2.0; phias=phias+phiat
                     idist=int(dist)+1
                     sold(idist)=sold(idist)+phiat-phirtt
                     !!c                    write(6,*) phias
                  end do
                  temp=0.0
                  write(6,*) "Writing sold(1:30) and temp "
                  do i=1,30
                     temp=temp+sold(i)
                     write(6,*) sold(i),temp
                  end do
                  write(6,*) " "
                  
                  !!c next find the colombic potential for that site from charges within the box
                  phic=0.0; phiac=0.0; ncrg=0
                  do i=1,nqass
                     it=0
                     if((chgpos(i).vorlt.bedge%min).or.&
                     &(chgpos(i).vorgt.bedge%max)) it=1
                     if(it.eq.0) then
                        ncrg=ncrg+1
                        vtemp=xo-chgpos(i); dist=sqrt(vtemp.dot.vtemp)
                        
                        if(dist.lt.5.0) then
                           if(dist.gt.1.e-6) then
                              temp=atmcrg(i)%value/dist
                              !!c b+++++++++++++++++++++++++++++++++++++
                              phic=phic + temp/atmeps(i)
                              !!c e+++++++++++++++++++++++++++++++++++++
                           end if
                           !!c find analytic potential from this real charge..=phiac
                           xu=chgpos(i); crgs=atmcrg(i)%value
                           xu2=((xu-oldmid)*scale)+goff
                           eps=atmeps(i)*epkt
                           call tops(xu2,xo2,crgs,eps,scale,&
                           &phiact,trgelf,1)
                           phiac=phiac+phiact
                        end if
                     end if
                  end do
                  !!c medeps, either epsin contain the 561.0 factor....
                  phiac=phiac*2.0
                  
                  !!c find the grid potentials..
                  call phintp(xn,phiv)
                  open(7,file="extra.dat")
                  write(7,*) phic,phir,phiv,phias,phiac,ncrg,ncrgs,tcrgs
                  close(7)
                  phit=phic+phir+phiv-phias-phiac
               else
                  phit=0.0
               end if
               
               !!c phit contains the total corrected potential
               
            end if
            
            if(isitr) then
               scomp=0.; sold=0.
               phir=0.
               do i=1,icount2b	
                  vtemp=xo-scspos(i); dist=sqrt(vtemp.dot.vtemp)
                  idist=int(dist)+1
                  if(idist.le.30)sold(idist)=sold(idist)+&
                  &(epkt*schrg(i)/dist)
                  phir=phir + schrg(i)/dist
               end do
               !!c medeps either epsin contains the 561.0 factor....
               phir=phir*epkt 
               do i=1,30
                  if(i.eq.1) scomp(i)=sold(i)
                  if(i.ne.1) scomp(i)=scomp(i-1)+sold(i)
               end do
               phirt=phirt+phir*chrgv
            end if
            
            if(isitrf.or.isitmd.or.isittf)then 
               !!c b++++++++++++++++++++++++++++++++++++++++++++
               !!c medeps either epsin contains the 561.0 factor....
               rxyz=rfield(nnatom)*epkt
               !!c e+++++++++++++++++++++++++++++++++++++++++++++
            end if
            
            if(isitcf.or.isitmd.or.isittf) then
               cxyz=coord(0.,0.,0.)
               if(abs(chrgv).gt.1.e-6) then
                  do i=1,nqass
                     vtemp=xo-chgpos(i); dist=vtemp.dot.vtemp
                     if(dist.gt.1.e-6) then
                        sdist=sqrt(dist)*dist
                        !!c b+++++++++++++++++++++
                        temp=atmcrg(i)%value/(atmeps(i)*sdist)
                        !!c e+++++++++++++++++++++
                        cxyz=cxyz+(vtemp*temp)
                     end if
                  end do
                  !!c atmeps and medeps and epsin contain the 561.0 factor....
                  cxyz=cxyz*chrgv
               end if
            end if
            
            if(isitc) then
               phic=0
               do i=1,nqass
                  vtemp=xo-chgpos(i); dist=vtemp.dot.vtemp
                  if(dist.gt.1.e-6) then
                     sdist=sqrt(dist)
                     temp=atmcrg(i)%value/sdist
                     !!c b++++++++++++++++++++++++++++
                     phic=phic + temp/atmeps(i)
                     !!c e++++++++++++++++++++++++++++
                  end if
               end do
               !!c atmeps and medeps and epsin contain the 561.0 factor....
               phict=phict+phic*chrgv
            end if
            
            !!c write out calculated/assigned charges
            oline=' '
            j=1
            if(isita) then
               oline(j:j+15)=atdes(1:16)
               j=j+20
            end if
            !!c NB need otemp cos can not write into a substring apparently
            !!c NB otemp needs to be at least 15 long to avoid an error!!
            if(isitx) then
               write(otemp,f200) xo%x
               oline(j:j+9)=otemp
               j=j+10
               write(otemp,f200) xo%y
               oline(j:j+9)=otemp
               j=j+10
               write(otemp,f200) xo%z
               oline(j:j+9)=otemp
               j=j+10
            end if
            
            if(isitq) then
               write(otemp,f200) chrgv
               oline(j:j+9)=otemp
               j=j+10
            end if
            
            if(isitp) then
               write(otemp,f200) phiv
               !!c               write(6,*) phiv,phi
               oline(j:j+9)=otemp
               j=j+10
            end if
            
            if(isiti) then
               write(otemp,f200) phii
               oline(j:j+9)=otemp
               j=j+10
            end if
            
            if(isitr) then
               write(otemp,f200) phir
               oline(j:j+9)=otemp
               j=j+10
            end if
            
            if(isitc) then
               write(otemp,f200) phic
               oline(j:j+9)=otemp
               j=j+10
            end if
            !!c b+++++++++++++++++++++++++++++
            if(isitap) then
               write(otemp,f200) aphi
               oline(j:j+9)=otemp
               j=j+10
            end if
            if(isitdeb) then
               write(otemp,f200) debyefraction
               oline(j:j+9)=otemp
               j=j+10
            end if
            !!c e+++++++++++++++++++++++++++++
            if(isitf) then
               write(otemp,f200) fxyz%x
               oline(j:j+9)=otemp
               j=j+10
               write(otemp,f200) fxyz%y
               oline(j:j+9)=otemp
               j=j+10
               write(otemp,f200) fxyz%z
               oline(j:j+9)=otemp
               j=j+10
            end if
            
            if(isitrf) then
               write(otemp,f200) rxyz%x
               oline(j:j+9)=otemp
               j=j+10
               write(otemp,f200) rxyz%y
               oline(j:j+9)=otemp
               j=j+10
               write(otemp,f200) rxyz%z
               oline(j:j+9)=otemp
               j=j+10
            end if
            
            if(isitcf) then
               write(otemp,f200) cxyz%x
               oline(j:j+9)=otemp
               j=j+10
               write(otemp,f200) cxyz%y
               oline(j:j+9)=otemp
               j=j+10
               write(otemp,f200) cxyz%z
               oline(j:j+9)=otemp
               j=j+10
            end if
            
            !!c b+++++++++++++++++++++++++++++
            if(isittf) then
               vtemp=rxyz+cxyz
               write(otemp,f200) vtemp%x
               oline(j:j+9)=otemp
               j=j+10
               write(otemp,f200) vtemp%y
               oline(j:j+9)=otemp
               j=j+10
               write(otemp,f200) vtemp%z
               oline(j:j+9)=otemp
               j=j+10
            end if
            if(isitmd) then
               vtemp=rxyz+cxyz
               atmforce(nnatom)=vtemp
               write(6,*)'atom:',nnatom,'rx=',rxyz%x,'cx=',&
               &cxyz%x,'tx=',vtemp%x
               write(6,*)'atom:',nnatom,'ry=',rxyz%y,'cy=',&
               &cxyz%y,'ty=',vtemp%y
               write(6,*)'atom:',nnatom,'rz=',rxyz%z,'cz=',&
               &cxyz%z,'tz=',vtemp%z
            end if
            !!c e+++++++++++++++++++++++++++++
            
            if(isitt) then
               write(otemp,f200) phit
               oline(j:j+9)=otemp
               j=j+10
            end if
            
            if(ofrm.and.isitmp1) write(16,'(a80)') oline
            if(.not.ofrm.and.isitmp1) then
               if(isita) write(16) atdes
               if(isitx) write(16) xo
               if(isitq) write(16) chrgv
               if(isitp) write(16) phiv
               if(isiti) write(16) phii
               if(isitr) write(16) phir
               if(isitc) write(16) phic
               !!c b++++++++++++++++++++++++++++++
               if(isitap) write(16) aphi
               !!c e++++++++++++++++++++++++++++++
               if(isitf) write(16) fxyz
               if(isitrf) write(16) rxyz
               if(isitcf) write(16) cxyz
               !!c b++++++++++++++++++++++++++++++
               if(isittf) write(16) rxyz+cxyz
               !!c e++++++++++++++++++++++++++++++
            end if
            
            !!c	end of file
         end do D304
      end if
      
      !!c b++++++++++++++++++++++++++++++++++++++++++++++
      if(isitsf) then
         do jj=1,ibnum
            i=atsurf(jj)
            if(atmsf(i).and.(atndx(jj).gt.0)) then
               !!c if the bgp belongs to the interesting site 
               !!c attention: using always radprb(1), in some case might be inappropriate
               xo=scspos(jj)+(radprb(1)*scsnor(jj))
               xn=((xo-oldmid)*scale)+goff
               xn%x=xn%x+1.          ! xn(1) = xn(1) + 1.
               call phintp(xn,fu%x)
               xn%x=xn%x-2.          ! xn(1) = xn(1) - 2.
               call phintp(xn,fl%x)
               xn%x=xn%x+1.; xn%y=xn%y+1.
               call phintp(xn,fu%y)
               xn%y=xn%y-2.          ! xn(2) = xn(2) - 2.
               call phintp(xn,fl%y)
               xn%y=xn%y+1. ; xn%z=xn%z+1.
               call phintp(xn,fu%z)
               xn%z=xn%z-2.          ! xn(3) = xn(3) - 2.
               call phintp(xn,fl%z)
               xn%z=xn%z+1.          ! xn(3) = xn(3) + 1.
               
               fxyz=(fl-fu)
               fn=0.5*scale*(fxyz.dot.scsnor(jj))
               ff=0.5*scale*sqrt(fxyz.dot.fxyz)
               
               if(ofrm) then
                  jtmp=j
                  write(otemp,f200) schrg(jj)
                  oline(jtmp:jtmp+9)=otemp
                  jtmp=jtmp+10
                  
                  write(otemp,f200) xo%x
                  oline(jtmp:jtmp+9)=otemp
                  jtmp=jtmp+10
                  write(otemp,f200) xo%y
                  oline(jtmp:jtmp+9)=otemp
                  jtmp=jtmp+10
                  write(otemp,f200) xo%z
                  oline(jtmp:jtmp+9)=otemp
                  jtmp=jtmp+10
                  
                  write(otemp,f200) fn
                  oline(jtmp:jtmp+9)=otemp
                  jtmp=jtmp+10
                  
                  write(otemp,f200) ff
                  oline(jtmp:jtmp+9)=otemp
                  jtmp=jtmp+10
                  
                  write(16,'(a80)') oline
               end if
               if(.not.ofrm) write(16) schrg(jj),fn
               
            endif
         end do
      end if
      if(allocated(atndx)) deallocate(atndx)
      if(allocated(atsurf)) deallocate(atsurf)
      !!c e++++++++++++++++++++++++++++++++++++++++++++++
      if(.not.iself) close(15)
      if(verbose) then
         write(6,*)'   '
         write(6,*)'number of atom coordinates read  : ',nnatom
         write(6,*)'   '
      end if
      etot = etot/2.
      if(ofrm) then
         if(frcfrm.eq.0)  then
            write(16,*)'total energy = ',etot,' kt'
            if(isitr) write(16,*)"corrected reaction field energy= ",&
            & phirt/2," kt"
            if(isitap)write(16,*)'Atomic potential for charged &
            &atoms is averaged over a spherical surface of &
            &less than',atompotdist,'A'
         end if
         if(frcfrm.eq.1) then
            write(16,*) "corrected reaction field energy= ",phirt/2," kt"
            write(16,*) "total coulombic energy     = ",phict/2," kt"
            if(isitap)write(16,*)'Atomic potential for charged &
            &atoms is averaged over a spherical surface &
            &of less than',atompotdist,'A'
         end if
         if(frcfrm.eq.2) then
            write(16,*) "corrected reaction field energy= ",phirt/2," kt"
            if(isitap)write(16,*)'Atomic potential for charged &
            &atoms is averaged over a spherical surface of &
            &less than',atompotdist,'A'
         end if
      end if
      close(16)
      
      !!c end of formatted frc read/write and unformatted frc write
      !!c END of unformatted frc.pdb read and frc write
      
      finish=0.
      if(verbose)write(6,*) 'frc stuff now done at',finish
      if(allocated(rfield)) deallocate(rfield)
      
      end subroutine wrtsit