!#####################################################################
!here the actual assignment from the hash tables is made
!2005-05-02   Other subroutine arguments are passed via modules
!2011-05-03   Parameters for subroutine ass (combined radass and crgass)
!             val   - output value of radius/charge
!             param - derived-type (parameter_file_record) array for 
!                     radius/charge records
!             hash  - derived-type (hash_table) array for radius/charge 
!                     hash table
!             Nmax  - size of derived-type arrays 
!2011-05-03   New parameters param and hash transfers either radius or 
!             charge parameters and hash arrays 
!2011-05-03   atm, res, rnum, chn are accessible via pointers module
!#####################################################################

      subroutine ass(val,param,hash,Nmax)

      !include "qdiffpar4.h" !removed due to module architecture

      type(parameter_file_record) :: param(Nmax)
      type(hash_table) :: hash(Nmax)
             
      character(3) :: tres='   '
      character(6) :: tatm='      '
      character(6) :: satm='*     ' !wildcard for atom
      character(1) :: achn=' '
      character(3) :: ares='   '
      character(4) :: arnum='    '

      !2011-05-03 Variants of records searched are combined from radass 
      !           and crgass subroutines.

      !2011-06-12 Declarations added due to IMPLICIT NONE
      real :: val
      integer :: Nmax,n,ifind

      !2011-06-16 Changed to rfindnew subroutine to avoid hash tables
      !Search for record (atm ; res ; rnum ; ch) 
      call rfindnew(atm,res,rnum,chn,ifind,n,param,Nmax)
      !call rfind(atm,res,rnum,chn,ifind,n,param,hash,Nmax)
      !call rfind(atm,res,rnum,chn,ifind,n)
      
      if (ifind.eq.0) then
         !Search for record (atm ; res ; rnum ; ' ') 
         call rfindnew(atm,res,rnum,achn,ifind,n,param,Nmax)
         !call rfind(atm,res,rnum,achn,ifind,n,param,hash,Nmax)
         
         if (ifind.eq.0) then
            !Search for record (atm ; res ; ' ' ; ' ') 
            call rfindnew(atm,res,arnum,achn,ifind,n,param,Nmax)
            !call rfin(atm,res,arnum,achn,ifind,n,param,hash,Nmax)

            if (ifind.eq.0) then
               !Search for record (atm ; ' ' ; ' ' ; ' ') 
               call rfindnew(atm,ares,arnum,achn,ifind,n,param,Nmax)
               !call rfind(atm,ares,arnum,achn,ifind,n,param,hash,Nmax)
               
               if (ifind.eq.0) then
                  tatm = atm(1:1)//'     '

                  !Search for record (generic_atom_name;res; ' ';' ')
                  call rfindnew(tatm,ares,arnum,achn,ifind,n,param,Nmax)
                  !call rfind(tatm,ares,arnum,achn,ifind,n,param,&
                  !           &hash,Nmax)
               
                  if (ifind.eq.0) then
                     ! Search for record (* ; ' ' ; ' ' ; ch) 
                     call rfindnew(satm,ares,arnum,chn,ifind,n,&
                                   &param,Nmax)
                     !call rfind(satm,ares,arnum,chn,ifind,n,param,&
                     !           &hash,Nmax) 
                     
                     if (ifind.eq.0) then
                        !Search for record (atm ; res ; ' ' ; ch) 
                        call rfindnew(atm,res,arnum,chn,ifind,n,&
                                      &param,Nmax) 
                        !call rfind(atm,res,arnum,chn,ifind,n,param,&
                        !           &hash,Nmax) 

                        if (ifind.eq.0) then
                           !Search for record (atm ; ' ' ; rnum ; ch) 
                           call rfindnew(atm,ares,rnum,chn,ifind,n,&
                                         &param,Nmax) 
                           !call rfind(atm,ares,rnum,chn,ifind,n,&
                           !           &param,hash,Nmax) 
                           
                           if (ifind.eq.0) then
                              !Search for record (atm;' ';' ';ch) 
                              call rfindnew(atm,ares,arnum,chn,ifind,&
                                            &n,param,Nmax) 
                              !call rfind(atm,ares,arnum,chn,ifind,&
                              !           &n,param,hash,Nmax) 
                              
                              if (ifind.eq.0) then
                                 !Search for record (atm;' ';rnum;' ') 
                                 call rfindnew(atm,ares,rnum,achn,&
                                               &ifind,n,param,Nmax) 
                                 !call rfind(atm,ares,rnum,achn,&
                                 !           &ifind,n,param,hash,Nmax) 
                              end if
                           end if
                        end if
                     end if
                  end if
               end if
            end if
         end if
      end if
      
      if (ifind.eq.0) then
         val=0.0; norecord=1;n=1
      else
         val = param(n)%value; norecord=0
      end if

      end subroutine ass

      !--------------------------------------------------------------------------
      !2011-05-03  Subroutines rfind and cfind are now combined
      !2011-05-03  Parameters added to combined subroutines
      !            param - derived-type (parameter_file_record) array 
      !                    for radius/charge records
      !            hash  - derived-type (hash_table) array for radius/
      !                    charge hash table
      !            Nmax  - size of derived-type arrays 
      subroutine rfind(tatm,tres,trnum,tchn,ifind,n,param,hash,Nmax)
      
      !Local variables specific only for this subroutine
      type(parameter_file_record) :: param(Nmax)
      type(hash_table) :: hash(Nmax)

      character(3) :: tres
      character(6) :: tatm
      character(1) :: tchn
      character(4) :: trnum

      !2011-06-12 Declarations added due to IMPLICIT NONE
      real :: val
      integer :: n,itemp,ifind,Nmax

      !find entry nres in hash table and check match with tres
      !include 'qdiffpar4.h' !removed due to module architecture
      !2011-05-02 Nrmax argument is added for ichash function (misc 
      !           module)

      n = ichash(tatm,tres,trnum,tchn,Nmax)

      !check for match
      ifind = 0
      
      do !2011-05-02 instead of "100 continue"
         !while no match and not at end of link list
         !2011-05-02  Changed to derived-type array for radii hash table
         !            (pointers module)
         write(6,*) '!atm=',tatm,'  res=',tres,'  rnum=',trnum,&
                    &'chn=',tchn
         write(6,*)' n=',n,' Record in hash: ',hash(n)
         write(6,*)' n=',n,' Record in param: ',param(hash(n)%inumb)
         
         if (hash(n)%inumb.eq.0) then !no match
            ifind = 0
            exit
         end if

         !2011-05-02  Changed to derive-type arrays for radii records
         !            and hash table
         itemp=hash(n)%inumb
         if ((tres.eq.param(itemp)%rnam)&
            &.and.(tatm.eq.param(itemp)%atnam)&
            &.and.(trnum.eq.param(itemp)%rnum)&
            &.and.(tchn.eq.param(itemp)%chn))  then
            n = hash(n)%inumb; ifind = 1; exit
         else
            if (hash(n)%ilink.ne.0) then
               n = hash(n)%ilink
            else
               ifind = 0; exit
            end if
         end if
      end do !2011-05-02 instead of go to 100
      
      write(6,*) '------------------------'
      
      end subroutine rfind

      !--------------------------------------------------------------- 
      !2011-06-16 Created new version of rfind subroutine to avoid
      !           hash tables 
      !--------------------------------------------------------------- 
      subroutine rfindnew(tatm,tres,trnum,tchn,ifind,n,param,Nmax)
      
      character(3) :: tres
      character(6) :: tatm
      character(1) :: tchn
      character(4) :: trnum
      
      type(parameter_file_record) :: param(Nmax)

      !2011-06-16 Declarations added due to IMPLICIT NONE
      integer :: n,i,ifind,Nmax

      ifind=0
      do i=1,Nmax
         if ((tres.eq.param(i)%rnam)&
            &.and.(tatm.eq.param(i)%atnam)&
            &.and.(trnum.eq.param(i)%rnum)&
            &.and.(tchn.eq.param(i)%chn))  then
            n=i;ifind=1; exit
         end if
      end do
      
      end subroutine rfindnew
