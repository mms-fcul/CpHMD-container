!#####################################################################
!2011-05-06 created by PK in order to define
!           operations with coord and int_coord 
!           derived-type variables
!#####################################################################

      module operators_on_coordinates

      use pointers, ONLY : coord,int_coord

      !2011-05-06  Expands - operator for coord
      interface operator(-)
         module procedure vminus, vminus2, vminus3, neg
 
         !2011-05-10  Expands - operator for int_coord
         module procedure ivminus, ivminus2, ivminus3, ineg
      end interface operator(-)

      !2011-05-06  Expands + operator for coord
      interface operator(+)
         module procedure vplus, vplus2, vplus3

         !2011-05-10  Expands + operator for int_coord
         module procedure ivplus, ivplus2, ivplus3
      end interface operator(+)

      !2011-05-06  Expands * operator for coord
      interface operator(*)
         module procedure vmult, vmult2, vmult3
      
         !2011-05-06  Expands * operator for int_coord
         module procedure ivmult, ivmult2, ivmult3
      end interface operator(*)

      !2011-05-06  Expands / operator for coord 
      interface operator(/)
         module procedure vdiv
         module procedure ivdiv,ivdiv2
      end interface operator(/)

      !2011-05-06  Scalar vector multiplication
      interface operator(.dot.)
         module procedure vdot
         module procedure ivdot
      end interface operator(.dot.)

      !2011-05-20  Vector vector multiplication
      interface operator(.x.)
         module procedure vx
         module procedure ivx
      end interface operator(.x.)

      !2011-05-10  Combines logical OR and LT
      interface operator(.vorlt.)
         module procedure vorlt, vorlt2
         module procedure ivorlt, ivorlt2
      end interface operator(.vorlt.)

      !2011-05-17  Combines logical AND and LT
      interface operator(.vandlt.)
         module procedure vandlt, vandlt2
         module procedure ivandlt, ivandlt2
      end interface operator(.vandlt.)

      !2011-05-13  Combines logical OR and LE
      interface operator(.vorle.)
         module procedure vorle, vorle2
         module procedure ivorle, ivorle2
      end interface operator(.vorle.)

      !2011-05-17  Combines logical AND and LE
      interface operator(.vandle.)
         module procedure vandle, vandle2
         module procedure ivandle,ivandle2
      end interface operator(.vandle.)

      !2011-05-10  Combines logical OR and GT
      interface operator(.vorgt.)
         module procedure vorgt, vorgt2
         module procedure ivorgt, ivorgt2
      end interface operator(.vorgt.)

      !2011-05-10  Combines logical OR and GE
      interface operator(.vorge.)
         module procedure vorge, vorge2
         module procedure ivorge, ivorge2
      end interface operator(.vorge.)

      !2011-05-10  Combines logical AND and GT
      interface operator(.vandgt.)
         module procedure vandgt, vandgt2
         module procedure ivandgt, ivandgt2
      end interface operator(.vandgt.)

      !2011-05-10  Combines logical AND and GE
      interface operator(.vandge.)
         module procedure vandge, vandge2
         module procedure ivandge, ivandge2
      end interface operator(.vandge.)

      !2011-05-13  Combines logical OR and NE
      interface operator(.vorne.)
         module procedure vorne, vorne2
         module procedure ivorne, ivorne2
      end interface operator(.vorne.)

      !2011-05-06  Expands sqrt intrinsic 
      interface sqrt
         module procedure vsqrt
      end interface sqrt

      !2011-05-17  Expands sum intrinsic 
      interface sum
         module procedure vsum
         module procedure ivsum
      end interface sum

      !2011-05-10  Expands int intrinsic 
      interface int
         module procedure ivint
      end interface int

      !2011-05-13  Expands nint intrinsic 
      interface nint
         module procedure ivnint
      end interface nint

      !2011-05-06  Expands float intrinsic 
      interface float
         module procedure vfloat
      end interface float

      !2011-05-06  Expands abs intrinsic 
      interface abs
         module procedure vabs
         module procedure ivabs
      end interface abs

      !2011-05-06  Expands min intrinsic 
      interface min
         module procedure vmin, vmin2, vmin3, vmin4
         module procedure ivmin, ivmin2, ivmin3, ivmin4
      end interface

      !2011-05-06  Expands max intrinsic 
      interface max
         module procedure vmax, vmax2, vmax3, vmax4
         module procedure ivmax, ivmax2, ivmax3, ivmax4
      end interface max

      !2011-05-07  New minsign operation
      interface minsign
         module procedure vmsign, vmsign2
      end interface minsign

      !2011-05-07  New maxsign operation
      interface maxsign
         module procedure vpsign, vpsign2
      end interface maxsign

      !2011-05-07  New submin operation
      interface submin
         module procedure vsubmin
      end interface submin

      !2011-05-07  New submax operation
      interface submax
         module procedure vsubmax
      end interface submax

      !2011-05-18  New comp operation
      interface comp
         module procedure vcomp, ivcomp
         !module procedure vcomp2, ivcomp2
      end interface comp

      !2011-05-13  Expands assignment 
      !interface assignment(=)
      !   module procedure vassign
      !   module procedure viassign
      !end interface assignment(=)

      contains
      !--------------------------------------------------------------
      !2011-05-06  Substracts two variables of coord derived type    
      function vminus(a,b)
         type(coord), intent(in) :: a,b
         type(coord) :: vminus
         vminus%x=a%x-b%x
         vminus%y=a%y-b%y
         vminus%z=a%z-b%z
      end function vminus

      !--------------------------------------------------------------
      !2011-05-10  Substracts two variables of int_coord derived type
      function ivminus(a,b)
         type(int_coord), intent(in) :: a,b
         type(int_coord) :: ivminus
         ivminus%i=a%i-b%i
         ivminus%j=a%j-b%j
         ivminus%k=a%k-b%k
      end function ivminus

      !--------------------------------------------------------------
      !2011-05-06  Substracts real variable from coord derived type    
      function vminus2(a,b)
         type(coord), intent(in) :: a
         real, intent(in) :: b
         type(coord) :: vminus2
         vminus2%x=a%x-b
         vminus2%y=a%y-b
         vminus2%z=a%z-b
      end function vminus2

      !--------------------------------------------------------------
      !2011-05-06  Substracts integer variable from int_coord derived 
      !            type
      function ivminus2(a,b)
         type(int_coord), intent(in) :: a
         integer, intent(in) :: b
         type(int_coord) :: ivminus2
         ivminus2%i=a%i-b
         ivminus2%j=a%j-b
         ivminus2%k=a%k-b
      end function ivminus2

      !--------------------------------------------------------------
      !2011-05-06  Substracts  coord derived type from real variable   
      function vminus3(a,b)
         type(coord), intent(in) :: b
         real, intent(in) :: a
         type(coord) :: vminus3
         vminus3%x=a-b%x
         vminus3%y=a-b%y
         vminus3%z=a-b%z
      end function vminus3

      !--------------------------------------------------------------
      !2011-05-10  Substracts int_coord derived type from integer 
      !            variable
      function ivminus3(a,b)
         type(int_coord), intent(in) :: b
         integer, intent(in) :: a
         type(int_coord) :: ivminus3
         ivminus3%i=a-b%i
         ivminus3%j=a-b%j
         ivminus3%k=a-b%k
      end function ivminus3

      !--------------------------------------------------------------


      !2011-05-06  Negation of variable of coord derived type    
      function neg(a)
         type(coord), intent(in) :: a
         type(coord) :: neg
         neg%x=-a%x
         neg%y=-a%y
         neg%z=-a%z
      end function neg

      !--------------------------------------------------------------
      !2011-05-10  Negation of variable of int_coord derived type    
      function ineg(a)
         type(int_coord), intent(in) :: a
         type(int_coord) :: ineg
         ineg%i=-a%i
         ineg%j=-a%j
         ineg%k=-a%k
      end function ineg

      !--------------------------------------------------------------
      !2011-05-06  Adds two variables of coord derived type    
      function vplus(a,b)
         type(coord), intent(in) :: a,b
         type(coord) :: vplus
         vplus%x=a%x+b%x
         vplus%y=a%y+b%y
         vplus%z=a%z+b%z
      end function vplus

      !--------------------------------------------------------------
      !2011-05-10  Adds two variables of int_coord derived type    
      function ivplus(a,b)
         type(int_coord), intent(in) :: a,b
         type(int_coord) :: ivplus
         ivplus%i=a%i+b%i
         ivplus%j=a%j+b%j
         ivplus%k=a%k+b%k
      end function ivplus

      !--------------------------------------------------------------
      !2011-05-06  Adds real variable to coord derived type    
      function vplus2(a,b)
         type(coord), intent(in) :: a
         real, intent(in) :: b
         type(coord) :: vplus2
         vplus2%x=a%x+b
         vplus2%y=a%y+b
         vplus2%z=a%z+b
      end function vplus2

      !--------------------------------------------------------------
      !2011-05-10  Adds integer variable to int_coord derived type    
      function ivplus2(a,b)
         type(int_coord), intent(in) :: a
         integer, intent(in) :: b
         type(int_coord) :: ivplus2
         ivplus2%i=a%i+b
         ivplus2%j=a%j+b
         ivplus2%k=a%k+b
      end function ivplus2

      !--------------------------------------------------------------
      !2011-05-06  Adds  coord derived type to real variable
      function vplus3(a,b)
         type(coord), intent(in) :: b
         real, intent(in) :: a
         type(coord) :: vplus3
         vplus3%x=a+b%x
         vplus3%y=a+b%y
         vplus3%z=a+b%z
      end function vplus3

      !--------------------------------------------------------------
      !2011-05-10  Adds  int_coord derived type to integer variable
      function ivplus3(a,b)
         type(int_coord), intent(in) :: b
         integer, intent(in) :: a
         type(int_coord) :: ivplus3
         ivplus3%i=a+b%i
         ivplus3%j=a+b%j
         ivplus3%k=a+b%k
      end function ivplus3

      !--------------------------------------------------------------
      !2011-05-06  Multiplies two variables coord derived type    
      function vmult(a,b)
         type(coord), intent(in) :: a,b
         type(coord) :: vmult
         vmult%x=a%x*b%x
         vmult%y=a%y*b%y
         vmult%z=a%z*b%z
      end function vmult

      !--------------------------------------------------------------
      !2011-05-06  Multiplies two variables int_coord derived type    
      function ivmult(a,b)
         type(int_coord), intent(in) :: a,b
         type(int_coord) :: ivmult
         ivmult%i=a%i*b%i
         ivmult%j=a%j*b%j
         ivmult%k=a%k*b%k
      end function ivmult

      !--------------------------------------------------------------
      !2011-05-06  Multiplies real variable to coord derived type    
      function vmult2(a,b)
         type(coord), intent(in) :: a
         real, intent(in) :: b
         type(coord) :: vmult2
         vmult2%x=a%x*b
         vmult2%y=a%y*b
         vmult2%z=a%z*b
      end function vmult2

      !--------------------------------------------------------------
      !2011-05-10  Multiplies integer variable to int_coord derived 
      !            type
      function ivmult2(a,b)
         type(int_coord), intent(in) :: a
         integer, intent(in) :: b
         type(int_coord) :: ivmult2
         ivmult2%i=a%i*b
         ivmult2%j=a%j*b
         ivmult2%k=a%k*b
      end function ivmult2

      !--------------------------------------------------------------
      !2011-05-10  Multiplies coord derived type to real variable   
      function vmult3(a,b)
         type(coord), intent(in) :: b
         real, intent(in) :: a
         type(coord) :: vmult3
         vmult3%x=a*b%x
         vmult3%y=a*b%y
         vmult3%z=a*b%z
      end function vmult3

      !--------------------------------------------------------------
      !2011-05-10  Multiplies int_coord derived type to integer 
      !            variable
      function ivmult3(a,b)
         type(int_coord), intent(in) :: b
         integer, intent(in) :: a
         type(int_coord) :: ivmult3
         ivmult3%i=a*b%i
         ivmult3%j=a*b%j
         ivmult3%k=a*b%k
      end function ivmult3

      !--------------------------------------------------------------
      !2011-05-06 Scalar vector product  for coord  
      function vdot(a,b)
         type(coord), intent(in) :: a,b
         real :: vdot
         vdot=a%x*b%x+a%y*b%y+a%z*b%z
      end function vdot

      !--------------------------------------------------------------
      !2011-05-10 Scalar vector product for int_coord
      function ivdot(a,b)
         type(int_coord), intent(in) :: a,b
         integer :: vdot
         ivdot=a%i*b%i+a%j*b%j+a%k*b%k
      end function ivdot

      !--------------------------------------------------------------
      !2011-05-20 Vector product  for coord   type
      function vx(a,b)
         type(coord), intent(in) :: a,b
         type(coord) :: vx
         vx%x=a%y*b%z-a%z*b%y
         vx%y=a%z*b%x-a%x*b%z
         vx%z=a%x*b%y-a%y*b%x
      end function vx

      !--------------------------------------------------------------
      !2011-05-20 Vector product  for coord   type
      function ivx(a,b)
         type(int_coord), intent(in) :: a,b
         type(int_coord) :: ivx
         ivx%i=a%j*b%k-a%k*b%j
         ivx%j=a%k*b%i-a%i*b%k
         ivx%k=a%i*b%j-a%j*b%i
      end function ivx

      !--------------------------------------------------------------
      !2011-05-06  Square root of variable of coord derived type    
      function vsqrt(a)
         type(coord), intent(in) :: a
         type(coord) vsqrt
         vsqrt%x=sqrt(a%x)
         vsqrt%y=sqrt(a%y)
         vsqrt%z=sqrt(a%z)
      end function vsqrt

      !--------------------------------------------------------------
      !2011-05-06  Absolute value of variable of coord derived type    
      function vabs(a)
         type(coord), intent(in) :: a
         type(coord) vabs
         vabs%x=abs(a%x)
         vabs%y=abs(a%y)
         vabs%z=abs(a%z)
      end function vabs

      !--------------------------------------------------------------
      !2011-05-10  Converts coord type into int_coord type
      function ivint(a)
         type(coord), intent(in) :: a
         type(int_coord) ivint
         ivint%i=int(a%x)
         ivint%j=int(a%y)
         ivint%k=int(a%z)
      end function ivint

      !--------------------------------------------------------------
      !2011-05-13  Converts coord type into int_coord type with nint
      function ivnint(a)
         type(coord), intent(in) :: a
         type(int_coord) ivnint
         ivnint%i=nint(a%x)
         ivnint%j=nint(a%y)
         ivnint%k=nint(a%z)
      end function ivnint

      !--------------------------------------------------------------
      !2011-05-13  Converts int_coord type into coord type
      function vfloat(a)
         type(int_coord), intent(in) :: a
         type(coord) :: vfloat
         vfloat%x = real(a%i)
         vfloat%y = real(a%j)
         vfloat%z = real(a%k)
      end function vfloat

      !--------------------------------------------------------------
      !2011-05-10  Absolute value of variable of int_coord derived type
      function ivabs(a)
         type(int_coord), intent(in) :: a
         type(int_coord) :: ivabs
         ivabs%i=abs(a%i)
         ivabs%j=abs(a%j)
         ivabs%k=abs(a%k)
      end function ivabs

      !--------------------------------------------------------------
      !2011-05-06  
      function vmin(b,a)
         type(coord), intent(in) :: a
         real :: b
         type(coord) vmin
         vmin%x=min(b,a%x)
         vmin%y=min(b,a%y)
         vmin%z=min(b,a%z)
      end function vmin

      !--------------------------------------------------------------
      !2011-05-10 
      function ivmin(b,a)
         type(int_coord), intent(in) :: a
         integer :: b
         type(int_coord) ivmin
         ivmin%i=min(b,a%i)
         ivmin%j=min(b,a%j)
         ivmin%k=min(b,a%k)
      end function ivmin

      !--------------------------------------------------------------
      !2011-05-06  
      function vmin2(b,a)
         type(coord), intent(in) :: b
         real :: a
         type(coord) vmin2
         vmin2%x=min(a,b%x)
         vmin2%y=min(a,b%y)
         vmin2%z=min(a,b%z)
      end function vmin2

      !--------------------------------------------------------------
      !2011-05-06  
      function ivmin2(b,a)
         type(int_coord), intent(in) :: b
         integer :: a
         type(int_coord) ivmin2
         ivmin2%i=min(a,b%i)
         ivmin2%j=min(a,b%j)
         ivmin2%k=min(a,b%k)
      end function ivmin2

      !--------------------------------------------------------------
      !2011-05-07  
      function vmin3(a,b)
         type(coord), intent(in) :: a,b
         type(coord) vmin3
         vmin3%x=min(a%x,b%x)
         vmin3%y=min(a%y,b%y)
         vmin3%z=min(a%z,b%z)
      end function vmin3

      !--------------------------------------------------------------
      !2011-05-10  
      function ivmin3(a,b)
         type(int_coord), intent(in) :: a,b
         type(int_coord) ivmin3
         ivmin3%i=min(a%i,b%i)
         ivmin3%j=min(a%j,b%j)
         ivmin3%k=min(a%k,b%k)
      end function ivmin3

      !--------------------------------------------------------------
      !2011-05-24  
      function vmin4(a)
         type(coord), intent(in) :: a
         real :: vmin4
         vmin4=min(a%x,a%y,a%z)
      end function vmin4

      !--------------------------------------------------------------
      !2011-05-24  
      function ivmin4(a)
         type(int_coord), intent(in) :: a
         integer :: ivmin4
         ivmin4=min(a%i,a%j,a%k)
      end function ivmin4

      !--------------------------------------------------------------
      !2011-05-06  
      function vmax(b,a)
         type(coord), intent(in) :: a
         real, intent(in) :: b
         type(coord) vmax
         vmax%x=max(b,a%x)
         vmax%y=max(b,a%y)
         vmax%z=max(b,a%z)
      end function vmax

      !--------------------------------------------------------------
      !2011-05-10 
      function ivmax(b,a)
         type(int_coord), intent(in) :: a
         integer, intent(in) :: b
         type(int_coord) ivmax
         ivmax%i=max(b,a%i)
         ivmax%j=max(b,a%j)
         ivmax%k=max(b,a%k)
      end function ivmax

      !--------------------------------------------------------------
      !2011-05-06  
      function vmax2(b,a)
         type(coord), intent(in) :: b
         real :: a
         type(coord) vmax2
         vmax2%x=max(a,b%x)
         vmax2%y=max(a,b%y)
         vmax2%z=max(a,b%z)
      end function vmax2

      !--------------------------------------------------------------
      !2011-05-06  
      function ivmax2(b,a)
         type(int_coord), intent(in) :: b
         integer :: a
         type(int_coord) ivmax2
         ivmax2%i=max(a,b%i)
         ivmax2%j=max(a,b%j)
         ivmax2%k=max(a,b%k)
      end function ivmax2

      !--------------------------------------------------------------
      !2011-05-07  
      function vmax3(a,b)
         type(coord), intent(in) :: a,b
         type(coord) vmax3
         vmax3%x=max(a%x,b%x)
         vmax3%y=max(a%y,b%y)
         vmax3%z=max(a%z,b%z)
      end function vmax3

      !--------------------------------------------------------------
      !2011-05-10  
      function ivmax3(a,b)
         type(int_coord), intent(in) :: a,b
         type(int_coord) ivmax3
         ivmax3%i=max(a%i,b%i)
         ivmax3%j=max(a%j,b%j)
         ivmax3%k=max(a%k,b%k)
      end function ivmax3

      !--------------------------------------------------------------
      !2011-05-24  
      function vmax4(a)
         type(coord), intent(in) :: a
         real :: vmax4
         vmax4=max(a%x,a%y,a%z)
      end function vmax4

      !--------------------------------------------------------------
      !2011-05-24  
      function ivmax4(a)
         type(int_coord), intent(in) :: a
         integer :: ivmax4
         ivmax4=max(a%i,a%j,a%k)
      end function ivmax4

      !--------------------------------------------------------------
      !2011-05-06  Divides coord derived type by real    
      function vdiv(a,b)
         type(coord), intent(in) :: a
         real, intent(in) :: b
         type(coord) :: vdiv
         vdiv%x=a%x/b
         vdiv%y=a%y/b
         vdiv%z=a%z/b
      end function vdiv

      !--------------------------------------------------------------
      !2011-05-10  Divides int_coord derived type by real    
      function ivdiv(a,b)
         type(int_coord), intent(in) :: a
         real, intent(in) :: b
         type(coord) :: ivdiv
         ivdiv%x=real(a%i)/b
         ivdiv%y=real(a%j)/b
         ivdiv%z=real(a%k)/b
      end function ivdiv

      !--------------------------------------------------------------
      !2011-05-20  Divides int_coord derived type by integer    
      function ivdiv2(a,b)
         type(int_coord), intent(in) :: a
         integer, intent(in) :: b
         type(int_coord) :: ivdiv2
         ivdiv2%i=a%i/b
         ivdiv2%j=a%j/b
         ivdiv2%k=a%k/b
      end function ivdiv2

      !--------------------------------------------------------------
      !2011-05-07
      function vmsign(b,a)
         type(coord),intent(in) :: a 
         real,intent(in) :: b
         type(coord) :: vmsign
         if(a%x.gt.0.) then
            vmsign%x=-b*a%x
         else
            vmsign%x=b*a%x
         end if
         if(a%y.gt.0.) then
            vmsign%y=-b*a%y
         else
            vmsign%y=b*a%y
         end if
         if(a%z.gt.0.) then
            vmsign%z=-b*a%z
         else
            vmsign%z=b*a%z
         end if
      end function vmsign

      !--------------------------------------------------------------
      !2011-05-07
      function vmsign2(b,a)
         type(coord),intent(in) :: b
         real,intent(in) :: a
         type(coord) :: vmsign2
         if(b%x.gt.0.) then
            vmsign2%x=-a*b%x
         else
            vmsign2%x=a*b%x
         end if
         if(b%y.gt.0.) then
            vmsign2%y=-a*b%y
         else
            vmsign2%y=a*b%y
         end if
         if(b%z.gt.0.) then
            vmsign2%z=-a*b%z
         else
            vmsign2%z=a*b%z
         end if
      end function vmsign2

      !--------------------------------------------------------------
      !2011-05-07
      function vpsign(b,a)
         type(coord),intent(in) :: a
         real,intent(in) :: b
         type(coord) :: vpsign
         if(a%x.lt.0.) then
            vpsign%x=-b*a%x
         else
            vpsign%x=b*a%x
         end if
         if(a%y.lt.0.) then
            vpsign%y=-b*a%y
         else
            vpsign%y=b*a%y
         end if
         if(a%z.lt.0.) then
            vpsign%z=-b*a%z
         else
            vpsign%z=b*a%z
         end if
      end function vpsign

      !--------------------------------------------------------------
      !2011-05-07
      function vpsign2(b,a)

         type(coord), intent(in) :: b
         real,intent(in) :: a
         type(coord) :: vpsign2
         if(b%x.lt.0.) then
            vpsign2%x=-a*b%x
         else
            vpsign2%x=a*b%x
         end if
         if(b%y.lt.0.) then
            vpsign2%y=-a*b%y
         else
            vpsign2%y=a*b%y
         end if
         if(b%z.lt.0.) then
            vpsign2%z=-a*b%z
         else
            vpsign2%z=a*b%z
         end if
      end function vpsign2

      !--------------------------------------------------------------
      !2011-05-07
      function vsubmin(a,b,c)
         type(coord), intent(in) :: a,b,c
         type(coord) :: vsubmin
         if(c%x.gt.0.) then
            vsubmin%x=b%x
         else
            vsubmin%x=a%x
         end if
          if(c%y.gt.0.) then
            vsubmin%y=b%y
         else
            vsubmin%y=a%y
         end if
          if(c%z.gt.0.) then
            vsubmin%z=b%z
         else
            vsubmin%z=a%z
         end if
      end function vsubmin

      !--------------------------------------------------------------
      !2011-05-07
      function vsubmax(a,b,c)
         type(coord), intent(in) :: a,b,c
         type(coord) :: vsubmax
         if(c%x.lt.0.) then
            vsubmax%x=b%x
         else
            vsubmax%x=a%x
         end if
          if(c%y.lt.0.) then
            vsubmax%y=b%y
         else
            vsubmax%y=a%y
         end if
          if(c%z.lt.0.) then
            vsubmax%z=b%z
         else
            vsubmax%z=a%z
         end if
      end function vsubmax

      !--------------------------------------------------------------
      !2011-05-10
      function vorlt(a,b)
         type(coord), intent(in) :: a,b
         logical vorlt
         if((a%x.lt.b%x).or.(a%y.lt.b%y).or.(a%z.lt.b%z)) then
            vorlt=.true.
         else
            vorlt=.false.
         end if
      end function vorlt

      !--------------------------------------------------------------
      !2011-05-17
      function ivorlt(a,b)
         type(int_coord), intent(in) :: a,b
         logical ivorlt
         if((a%i.lt.b%i).or.(a%j.lt.b%j).or.(a%k.lt.b%k)) then
            ivorlt=.true.
         else
            ivorlt=.false.
         end if
      end function ivorlt

      !--------------------------------------------------------------
      !2011-05-10
      function vorlt2(a,b)
         type(coord), intent(in) :: a
         real, intent(in) :: b
         logical vorlt2
         if((a%x.lt.b).or.(a%y.lt.b).or.(a%z.lt.b)) then
            vorlt2=.true.
         else
            vorlt2=.false.
         end if
      end function vorlt2

      !--------------------------------------------------------------
      !2011-05-17
      function ivorlt2(a,b)
         type(int_coord), intent(in) :: a
         integer, intent(in) :: b
         logical ivorlt2
         if((a%i.lt.b).or.(a%j.lt.b).or.(a%k.lt.b)) then
            ivorlt2=.true.
         else
            ivorlt2=.false.
         end if
      end function ivorlt2

      !--------------------------------------------------------------
      !2011-05-17
      function vandlt(a,b)
         type(coord), intent(in) :: a,b
         logical vandlt
         if((a%x.lt.b%x).and.(a%y.lt.b%y).and.(a%z.lt.b%z)) then
            vandlt=.true.
         else
            vandlt=.false.
         end if
      end function vandlt

      !--------------------------------------------------------------
      !2011-05-17
      function ivandlt(a,b)
         type(int_coord), intent(in) :: a,b
         logical ivandlt
         if((a%i.lt.b%i).and.(a%j.lt.b%j).and.(a%k.lt.b%k)) then
            ivandlt=.true.
         else
            ivandlt=.false.
         end if
      end function ivandlt

      !--------------------------------------------------------------
      !2011-05-17
      function vandlt2(a,b)
         type(coord), intent(in) :: a
         real, intent(in) :: b
         logical vandlt2
         if((a%x.lt.b).and.(a%y.lt.b).and.(a%z.lt.b)) then
            vandlt2=.true.
         else
            vandlt2=.false.
         end if
      end function vandlt2

      !--------------------------------------------------------------
      !2011-05-17
      function ivandlt2(a,b)
         type(int_coord), intent(in) :: a
         integer, intent(in) :: b
         logical ivandlt2
         if((a%i.lt.b).and.(a%j.lt.b).and.(a%k.lt.b)) then
            ivandlt2=.true.
         else
            ivandlt2=.false.
         end if
      end function ivandlt2

      !--------------------------------------------------------------
      !2011-05-13
      function vorle(a,b)
         type(coord), intent(in) :: a,b
         logical vorle
         if((a%x.le.b%x).or.(a%y.le.b%y).or.(a%z.le.b%z)) then
            vorle=.true.
         else
            vorle=.false.
         end if
      end function vorle

      !--------------------------------------------------------------
      !2011-05-13
      function ivorle(a,b)
         type(int_coord), intent(in) :: a,b
         logical ivorle

         if((a%i.le.b%i).or.(a%j.le.b%j).or.(a%k.le.b%k)) then
            ivorle=.true.
         else
            ivorle=.false.
         end if
      end function ivorle
      
      !--------------------------------------------------------------
      !2011-05-10
      function vorle2(a,b)
         type(coord), intent(in) :: a
         real, intent(in) :: b
         logical vorle2
   
         if((a%x.le.b).or.(a%y.le.b).or.(a%z.le.b)) then
            vorle2=.true.
         else
            vorle2=.false.
         end if
      end function vorle2

      !--------------------------------------------------------------
      !2011-05-10
      function ivorle2(a,b)
         type(int_coord), intent(in) :: a
         integer, intent(in) :: b
         logical :: ivorle2
    
         if((a%i.le.b).or.(a%j.le.b).or.(a%k.le.b)) then
            ivorle2=.true.
         else
            ivorle2=.false.
         end if
      end function ivorle2

      !--------------------------------------------------------------
      !2011-05-17
      function vandle(a,b)
         type(coord), intent(in) :: a,b
         logical vandle
  
         if((a%x.le.b%x).and.(a%y.le.b%y).and.(a%z.le.b%z)) then
            vandle=.true.
         else
            vandle=.false.
         end if
      end function vandle

      !--------------------------------------------------------------
      !2011-05-17
      function ivandle(a,b)
         type(int_coord), intent(in) :: a,b
         logical ivandle
        
         if((a%i.le.b%i).and.(a%j.le.b%j).and.(a%k.le.b%k)) then
            ivandle=.true.
         else
            ivandle=.false.
         end if
      end function ivandle

      !--------------------------------------------------------------
      !2011-05-17
      function vandle2(a,b)
         type(coord), intent(in) :: a
         real, intent(in) :: b
         logical vandle2
    
         if((a%x.le.b).and.(a%y.le.b).and.(a%z.le.b)) then
            vandle2=.true.
         else
            vandle2=.false.
         end if
      end function vandle2

      !--------------------------------------------------------------
      !2011-05-17
      function ivandle2(a,b)
         type(int_coord), intent(in) :: a
         integer, intent(in) :: b
         logical ivandle2
        
         if((a%i.le.b).and.(a%j.le.b).and.(a%k.le.b)) then
            ivandle2=.true.
         else
            ivandle2=.false.
         end if
      end function ivandle2

      !--------------------------------------------------------------
      !2011-05-10
      function vorgt(a,b)
         type(coord), intent(in) :: a,b
         logical vorgt
     
         if((a%x.gt.b%x).or.(a%y.gt.b%y).or.(a%z.gt.b%z)) then
            vorgt=.true.
         else
            vorgt=.false.
         end if
      end function vorgt

      !--------------------------------------------------------------
      !2011-05-10
      function ivorgt(a,b)
         type(int_coord), intent(in) :: a,b
         logical ivorgt
      
         if((a%i.gt.b%i).or.(a%j.gt.b%j).or.(a%k.gt.b%k)) then
            ivorgt=.true.
         else
            ivorgt=.false.
         end if
      end function ivorgt

      !--------------------------------------------------------------
      !2011-05-10
      function vorgt2(a,b)
         type(coord), intent(in) :: a
         real, intent(in) :: b
         logical vorgt2
  
         if((a%x.gt.b).or.(a%y.gt.b).or.(a%z.gt.b)) then
            vorgt2=.true.
         else
            vorgt2=.false.
         end if
      end function vorgt2

      !--------------------------------------------------------------
      !2011-05-10
      function ivorgt2(a,b)
         type(int_coord), intent(in) :: a
         integer, intent(in) :: b
         logical ivorgt2
   
         if((a%i.gt.b).or.(a%j.gt.b).or.(a%k.gt.b)) then
            ivorgt2=.true.
         else
            ivorgt2=.false.
         end if
      end function ivorgt2

      !--------------------------------------------------------------
      !2011-05-10
      function vorge(a,b)
         type(coord), intent(in) :: a,b
         logical vorge
        
         if((a%x.ge.b%x).or.(a%y.ge.b%y).or.(a%z.ge.b%z)) then
            vorge=.true.
         else
            vorge=.false.
         end if
      end function vorge

      !--------------------------------------------------------------
      !2011-05-10
      function ivorge(a,b)
         type(int_coord), intent(in) :: a,b
         logical ivorge
    
         if((a%i.ge.b%i).or.(a%j.ge.b%j).or.(a%k.ge.b%k)) then
            ivorge=.true.
         else
            ivorge=.false.
         end if
      end function ivorge

      !--------------------------------------------------------------
      !2011-05-10
      function vorge2(a,b)
         type(coord), intent(in) :: a
         real, intent(in) :: b
         logical vorge2
    
         if((a%x.ge.b).or.(a%y.ge.b).or.(a%z.ge.b)) then
            vorge2=.true.
         else
            vorge2=.false.
         end if
      end function vorge2

      !--------------------------------------------------------------
      !2011-05-10
      function ivorge2(a,b)
         type(int_coord), intent(in) :: a
         integer, intent(in) :: b
         logical ivorge2
    
         if((a%i.ge.b).or.(a%j.ge.b).or.(a%k.ge.b)) then
            ivorge2=.true.
         else
            ivorge2=.false.
         end if
      end function ivorge2

      !--------------------------------------------------------------
      !2011-05-17
      function vandgt(a,b)
         type(coord), intent(in) :: a,b
         logical :: vandgt
        
         if((a%x.gt.b%x).and.(a%y.gt.b%y).and.(a%z.gt.b%z)) then
            vandgt=.true.
         else
            vandgt=.false.
         end if
      end function vandgt

      !--------------------------------------------------------------
      !2011-05-17
      function ivandgt(a,b)
         type(int_coord), intent(in) :: a,b
         logical :: ivandgt
        
         if((a%i.gt.b%i).and.(a%j.gt.b%j).and.(a%k.gt.b%k)) then
            ivandgt=.true.
         else
            ivandgt=.false.
         end if
      end function ivandgt

      !--------------------------------------------------------------
      !2011-05-17
      function vandgt2(a,b)
         type(coord), intent(in) :: a
         real, intent(in) :: b
         logical :: vandgt2
    
         if((a%x.gt.b).and.(a%y.gt.b).and.(a%z.gt.b)) then
            vandgt2=.true.
         else
            vandgt2=.false.
         end if
      end function vandgt2

      !--------------------------------------------------------------
      !2011-05-17
      function ivandgt2(a,b)
         type(int_coord), intent(in) :: a
         integer, intent(in) :: b
         logical :: ivandgt2
        
         if((a%i.gt.b).and.(a%j.gt.b).and.(a%k.gt.b)) then
            ivandgt2=.true.
         else
            ivandgt2=.false.
         end if
      end function ivandgt2

      !--------------------------------------------------------------
      !2011-05-17
      function vandge(a,b)
         type(coord), intent(in) :: a,b
         logical :: vandge
        
         if((a%x.ge.b%x).and.(a%y.ge.b%y).and.(a%z.ge.b%z)) then
            vandge=.true.
         else
            vandge=.false.
         end if
      end function vandge

      !--------------------------------------------------------------
      !2011-05-17
      function ivandge(a,b)
         type(int_coord), intent(in) :: a,b
         logical :: ivandge
        
         if((a%i.ge.b%i).and.(a%j.ge.b%j).and.(a%k.ge.b%k)) then
            ivandge=.true.
         else
            ivandge=.false.
         end if
      end function ivandge

      !--------------------------------------------------------------
      !2011-05-17
      function vandge2(a,b)
         type(coord), intent(in) :: a
         real, intent(in) :: b
         logical :: vandge2
    
         if((a%x.ge.b).and.(a%y.ge.b).and.(a%z.ge.b)) then
            vandge2=.true.
         else
            vandge2=.false.
         end if
      end function vandge2

      !--------------------------------------------------------------
      !2011-05-17
      function ivandge2(a,b)
         type(int_coord), intent(in) :: a
         integer, intent(in) :: b
         logical :: ivandge2
     
         if((a%i.ge.b).and.(a%j.ge.b).and.(a%k.ge.b)) then
            ivandge2=.true.
         else
            ivandge2=.false.
         end if
      end function ivandge2

      !--------------------------------------------------------------
      !2011-05-13
      function vorne(a,b)
         type(coord), intent(in) :: a
         real, intent(in) :: b
         logical vorne
    
         if((a%x.ne.b).or.(a%y.ne.b).or.(a%z.ne.b)) then
            vorne=.true.
         else
            vorne=.false.
         end if
      end function vorne

      !--------------------------------------------------------------
      !2011-05-13
      function vorne2(a,b)
         type(coord), intent(in) :: a,b
         logical vorne2
        
         if((a%x.ne.b%x).or.(a%y.ne.b%y).or.(a%z.ne.b%z)) then
            vorne2=.true.
         else
            vorne2=.false.
         end if
      end function vorne2

      !--------------------------------------------------------------
      !2011-05-13
      function ivorne(a,b)
         type(int_coord), intent(in) :: a
         integer, intent(in) :: b
         logical ivorne

         if((a%i.ne.b).or.(a%j.ne.b).or.(a%k.ne.b)) then
            ivorne=.true.
         else
            ivorne=.false.
         end if
      end function ivorne

      !--------------------------------------------------------------
      !2011-05-13
      function ivorne2(a,b)
         type(int_coord), intent(in) :: a,b
         logical ivorne2
        
         if((a%i.ne.b%i).or.(a%j.ne.b%j).or.(a%k.ne.b%k)) then
            ivorne2=.true.
         else
            ivorne2=.false.
         end if
      end function ivorne2

      !--------------------------------------------------------------
      !2011-05-13
      function vsum(a)
         type(coord), intent(in) :: a
         real :: vsum
         vsum=a%x+a%y+a%z
      end function vsum

      !--------------------------------------------------------------
      !2011-05-17
      function ivsum(a)
         type(int_coord), intent(in) :: a
         integer :: ivsum
         ivsum=a%i+a%j+a%k
      end function ivsum

      !--------------------------------------------------------------
      !2011-05-13
      !subroutine vassign(a,b)
      !   type(coord),intent(out) :: a
      !   real, intent(in) :: b
      !   a%x=b ; a%y=b ; a%z=b
      !end subroutine vassign

      !--------------------------------------------------------------
      !2011-05-13
      !subroutine viassign(a,b)
      !   type(int_coord),intent(out) :: a
      !   integer, intent(in) :: b
      !   a%i=b ; a%j=b ; a%k=b
      !end subroutine viassign

      !--------------------------------------------------------------
      !2011-05-18
      function vcomp(a,b)
         type(coord), intent(in) :: a
         integer, intent(in) :: b
         real vcomp
        
         select case(b)
            case(1)
               vcomp=a%x
            case(2)
               vcomp=a%y
            case(3)
               vcomp=a%z
            case default
               vcomp=0.0
         end select
      end function vcomp

      !--------------------------------------------------------------
      !2011-05-18
      function ivcomp(a,b)
         type(int_coord), intent(in) :: a
         integer, intent(in) :: b
         integer ivcomp
    
         select case(b)
            case(1)
               ivcomp=a%i
            case(2)
               ivcomp=a%j
            case(3)
               ivcomp=a%k
            case default
               ivcomp=0
         end select
      end function ivcomp

      !--------------------------------------------------------------
      !2011-05-18
      !function vcomp2(a,b)
      !   real, intent(in) :: a
      !   integer, intent(in) :: b
      !   type(coord) vcomp2
      !   select case(b)
      !      case(1)
      !         vcomp2%x=a
      !      case(2)
      !         vcomp2%y=a
      !      case(3)
      !         vcomp2%z=a
      !   end select
      !end function vcomp2

      !--------------------------------------------------------------
      !2011-05-18
      !function ivcomp2(a,b)
      !   integer, intent(in) :: a
      !   integer, intent(in) :: b
      !   type(int_coord) :: ivcomp2
      !   select case(b)
      !      case(1)
      !         ivcomp2%i=a
      !      case(2)
      !         ivcomp2%j=a
      !      case(3)
      !         ivcomp2%k=a
      !   end select
      !end function ivcomp2


      end module operators_on_coordinates
