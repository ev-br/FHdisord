!-----------------------------------------------------------
!  adding/dropping a row and a column via N**2 operations 
!   +  full N**3 inverse 
!   +  changing a single row
!  using BLAS / LAPACK
!
!
!  this version [det_n2___.f, w/triple underscore] : 
!   1. fixed the SIGN of the inv_m1
!   2. rank-two updates for both _p2 & _m2
!   3. change the contents of a row & column [det only]
!
!-----------------------------------------------------------
	module det_n2
	implicit none; save

!------------------------------------
! Implementation issues:
!   1) inv_p1: need to set a vector = 0; 
!       Currnet implementation uses dscal; 
!       Will it be any better to use dcopy from a pre-set attay of zeros?
!   2) inv_p1 & inv_m1: outer product of two vectors;
!        Current implementation uses dger; Would dgemm() do better?
!FIXED   3) inv_m1: row & column to be dropped are swapped w/the last ones,
!           and the former are then being dropped. 
!           The SIGN of the determinant might get wrong, and it is NOT
!           being taken care of. 
!------------------------------------

	public :: full_inv, det_p1, inv_p1, det_m1, inv_m1, det_r, inv_r
!	public :: det_p2

	public :: det_p2_, inv_p2, det_m2, inv_m2, det_rc

	private

	real*8, external :: ddot



!*****************************************************************
	contains



!-------------------------------
!--- Det: add a row & a column
!-------------------------------
	real*8 function  det_p1(pm,lda,a,u,v,z,s)
	integer :: lda, pm       ! leading dimension and actual size of A
	real*8  :: a(lda,lda)
      real*8  :: v(lda),u(lda),z(lda),s
!
!  Input: a(lda,lda), s, v(lda), u(lda)
!         all arrays are of the size pm<lda 
!
!  Output: z(lda)
!          det. ratio @ det_p1

!  z=a1 DOT u
	call dgemv('N',pm,pm,1.d0,a,lda,u,1,0.d0,z,1)


!  \lambda = v DOT z; \rho = s - \lambda; det = 1/\rho
	det_p1 = s - ddot(pm,v,1,z,1)

	end function det_p1


!-------------------------------
!--- Add row & column: update the inverse 
!-------------------------------
	subroutine  inv_p1(pm,lda,a,det,v,w,z)
	integer :: lda, pm       ! leading dimension and actual size of A
	real*8  :: a(lda,lda),det,v(lda),w(lda),z(lda) 

	real*8 :: rho
!
!  Input: a(lda,lda), v(lda)
!         z(lda), det --- as set by det_p1
!         
!          all arrays are of the size pm<lda 
!
!  Output: a(lda,lda) ==> (pm+1) * (pm+1) inverse
!          pm = pm( @ input) + 1

	rho = 1.d0/det	

!  w = transp(a1)*v
	call dgemv('T',pm,pm,1.d0,a,lda,v,1,0.d0,w,1)

! b^{-1} 
	call dger(pm,pm,rho,z,1,w,1,a,lda)

! last row
       call dscal(pm,0.d0,a(pm+1,1),lda) ! set =0 first
       call daxpy(pm, -rho, w,1,a(pm+1,1),lda)

! last column
       call dscal(pm,0.d0,a(1,pm+1),1) ! set =0 first
       call daxpy(pm, -rho, z,1,a(1,pm+1),1)

! (pm+1,pm+1)
	a(pm+1,pm+1)=rho

! update pm
	pm=pm+1

	end subroutine inv_p1



!-------------------------------
!--- Det: drop r-th row & c-th column
!-------------------------------
      real*8 function  det_m1(pm,lda,a,r,c)
      integer :: pm,lda       ! actual size & leading dimension of A
      real*8  :: a(lda,lda)    
	integer:: c,r
	integer :: s
!
!  Input:  a(lda,lda), pm
!
!  Output: det. ratio:  det BIG / det SMALL (=1/rho)
!

	s=1; if(c/=pm)s=-s; if(r/=pm)s=-s
	det_m1 = 1.d0*s/a(c,r)

      end function det_m1


!-------------------------------
!--- Inv: drop r-th row and c-th column
!-------------------------------
      subroutine  inv_m1(pm,lda,a,r,c)
      integer :: lda, pm       ! leading dimension and actual size of A
      real*8  :: a(lda,lda)
      integer :: c,r

	real*8 :: rh

!  Input: a(lda,lda) pm*pm
!         r, c --- row and  column to be dropped
!
!  Output: a(lda,lda) ==> (pm-1) * (pm-1)
!          pm = pm( @ input ) - 1
!
!   How:  swaps to-be-dropped and last row & cols, 
!         and then drops the former ones
!
!   The latter is done using: 
!       a(pm,pm)         is   \rho
!       a(1:pm,pm)       is  -\rho * z  -- last column
!       a(pm,1:pm)       is  -\rho * w  -- last row
!       a(1:pm-1,1:pm-1) is  a^{-1} + \rho z*w 

! swap c-th and last row (c-th column of A <==> c-th row of A^{-1})
	if(c/=pm) call dswap(pm,a(c,1),lda,a(pm,1),lda)
      
! swap r-th and last column
	if(r/=pm) call dswap(pm,a(1,r),1,a(1,pm),1)
	
!-------------------------------- drop the last row & column
	rh = -1.d0/a(pm,pm)

! stripe out the outer product z*w
	call dger(pm-1,pm-1,rh,a(1,pm),1,a(pm,1),lda,a,lda)

! update pm
	pm = pm - 1


      end subroutine inv_m1


!-------------------------------
!--- Det: change a single row
!-------------------------------
	real*8 function  det_r(pm,lda,a,r,v)
	integer :: lda, pm       ! leading dimension and actual size of A
	real*8  :: a(lda,lda), v(lda)
      integer :: r
!
!  Input: a(lda,lda) --- inverse matrix
!         v(lda)     --- a row to be added
!         r          --- a row number 
!           all arrays are of the size pm<lda 
!
!  Output:
!           det. ratio @ det_r

! \lambda =  ( v DOT A_r )
      det_r = 1.d0 + ddot(pm,v,1,a(1,r),1)

	end function det_r


!-------------------------------
!--- inv: change a single row
!-------------------------------
	subroutine  inv_r(pm,lda,a,r,det,v,w,z)
	integer :: lda, pm       ! leading dimension and actual size of A
	real*8  :: a(lda,lda), v(lda),w(lda),z(lda),det
      integer :: r

      real*8  :: rho
!
!  Input: a(lda,lda) --- inverse matrix
!         v(lda)     --- a row to be added
!         w(lda),z(lda)     --- working arrays
!         r          --- row number 
!         det        --- det ratio, as set by det_r
!           all arrays are of the size pm<lda 
!
!  Output:
!         a(lda,lda) contains an updated inverse matrix


      rho = -1.d0/det

! z_i = A_{i,r}
      call dcopy(pm,a(1,r),1,z,1)

! w = v DOT A 
      call dgemv('T',pm,pm,1.d0,a,lda,v,1,0.d0,w,1)


! A+ \rho* z DOT w^T
      call dger(pm,pm,rho,z,1,w,1,a,lda)

! one cannot get rid of z() array since if one otherwise plugs {a(1,r),1}
! directly into dger(...) instead of {z,1}, it all goes nuts.
! probably, dger() spoils the z array via blocking or the like.

	end subroutine inv_r


!-------------------------------
!--- det: add two rows & columns
!-------------------------------
	real*8 function  det_p2_(pm,lda,a,u2,v2,s2,z2) 
	integer :: lda, pm       ! leading dimension and actual size of A
	real*8  :: a(lda,lda)
        real*8  :: v2(lda,2),u2(lda,2),s2(2,2),z2(lda,2)

!
!  Input: a(lda,lda), v2(lda,2), u2(lda,2), s2(2,2)
!         all arrays are of the size pm<lda 
!
!  Output:    det. ratio @ det_p2
!             s2 (lda,2) --- original s2 is LOST
!	      z2(lda,2)
!

! z2 = a^{-1} . u2
	call dgemm('N','N',pm,2,pm,1.d0,a,lda,u2,lda,0.d0,z2,lda)

! s2 = s2 - v^T . z
	call dgemm('T','N',2,2,pm,-1.d0,v2,lda,z2,lda,1.d0,s2,2)


! return det|c2|
	det_p2_ = s2(1,1)*s2(2,2) - s2(1,2)*s2(2,1)

	end function det_p2_



!-------------------------------
!--- inv: add two rows & columns
!-------------------------------
	subroutine inv_p2(pm,lda,a,det,c2,z2,v2,w2,x2)
	integer :: lda, pm       ! leading dimension and actual size of A
	real*8  :: det
	real*8  :: a(lda,lda)
        real*8  :: w2(lda,2),z2(lda,2),v2(lda,2),x2(lda,2),c2(2,2)
	real*8  :: cc(2,2)

!
!  Input: a(lda,lda), 
!         z2(lda,2), c2(2,2), det -- as set by det_p2_
!         v2(lda,2)
!         w2(lda,2),x2(lda,2) -- working arrays
!         all arrays are of the size pm<lda 
!
!  Output: a(lda,lda) ==> (pm+2) * (pm+2) inverse
!          pm = pm( @ input) + 2


	!print*,'-------- inv_p2'; print*


! invert c2
	cc(1,1) =  c2(2,2); cc(2,2) =  c2(1,1)
	cc(1,2) = -c2(1,2); cc(2,1) = -c2(2,1)
	cc = cc/det

	a(pm+1,pm+1)=cc(1,1); a(pm+2,pm+1)=cc(2,1);  ! a corner
	a(pm+1,pm+2)=cc(1,2); a(pm+2,pm+2)=cc(2,2);

! last 2 cols : x=z DOT cc, & q = -z
        call dgemm('N','N',pm,2,2,1.d0,z2,lda,cc,2,0.0,x2,lda ) 
        call dscal(pm,0.d0,a(1,pm+1),1)                 ! set =0 first
        call dscal(pm,0.d0,a(1,pm+2),1) 
        call daxpy(pm, -1.d0, x2,1,a(1,pm+1),1)         ! copy:
        call daxpy(pm, -1.d0, x2(1,2),1,a(1,pm+2),1)    ! [ x2 & x2(1,2) are the starting addr of the 1st and 2nd clmn, resp ]

	!print*,'-------- inv_p2 1'; 

! now need w= ( A^{-1} )^T . v
	call dgemm('T','N',pm,2,pm,1.d0,a,lda,v2,lda,0.d0,w2,lda)

	!print*,'-------- inv_p2 2'; 

! P = A^{-1} + x. w^T
	call dgemm('N','T',pm,pm,2,1.d0,x2,lda,w2,lda,1.d0,a,lda)

	!print*,'-------- inv_p2 3'; 

! last 2 rows :  z = w . cc^T   (can overwite z now) 
	call dgemm('N','T',pm,2,2,1.d0,w2,lda,cc,2,0.d0,z2,lda)

	!print*,'-------- inv_p2 4'; 


        call dscal(pm,0.d0,a(pm+1,1),lda)             ! set =0 first
        call dscal(pm,0.d0,a(pm+2,1),lda)  
        call daxpy(pm, -1.d0, z2,1,a(pm+1,1),lda)     ! copy
        call daxpy(pm, -1.d0, z2(1,2),1,a(pm+2,1),lda)

! update pm
	pm = pm + 2


	end subroutine inv_p2



!-------------------------------
!--- det: remove two rows & columns
!-------------------------------
	real*8 function  det_m2(pm,lda,a,r1,r2,c1,c2) 
	integer :: lda, pm       ! leading dimension and actual size of A
	real*8  :: a(lda,lda)
	integer :: c1,c2,r1,r2

	integer :: cc1,cc2,rr1,rr2
	real*8  :: det

!
!  Input: a(lda,lda), 
!         c1,c2,r1,r2 -- rows & cols to be dropped
!
!  Output:    det. ratio @ det_p2
!

	cc1=min(c1,c2); cc2=max(c1,c2)
	rr1=min(r1,r2); rr2=max(r1,r2)    ! sort

	det = a(cc1,rr1)*a(cc2,rr2) - a(cc1,rr2)*a(cc2,rr1)

	det_m2=1.d0/det


	end function det_m2


!-------------------------------
!--- inv: remove two rows & columns
!-------------------------------
	subroutine inv_m2(pm,lda,a,det,r1,r2,c1,c2,z2,w2)
	integer :: lda, pm       ! leading dimension and actual size of A
	real*8  :: det
	real*8  :: a(lda,lda), z2(lda,2),w2(lda,2)
	integer :: c1,c2,r1,r2
	
	integer :: cc1,cc2,rr1,rr2,     j
	real*8  :: ss(2,2)
!
!  Input: a(lda,lda)
!         det         -- as set by det_m2 
!         c1,c2,r1,r2 -- rows & cols to be dropped
!         z2(lda,2), w2(lda,2) -- working arrays
!
!  Output: a(lda,lda) ==> (pm-2) * (pm-2) inverse
!          pm = pm( @ input) - 2
!

	cc1=min(c1,c2); cc2=max(c1,c2)
	rr1=min(r1,r2); rr2=max(r1,r2)    ! sort

! swap cc2-nd and last row (c-th column of A <==> c-th row of A^{-1})
	if(cc2/=pm)   call dswap(pm,a(cc2,1),lda,a(pm,1),lda)	
	if(cc1/=pm-1) call dswap(pm,a(cc1,1),lda,a(pm-1,1),lda)	

! swap rr2-nd and last column
	if(rr2/=pm)   call dswap(pm,a(1,rr2),1,a(1,pm),1)
	if(rr1/=pm-1) call dswap(pm,a(1,rr1),1,a(1,pm-1),1)

! invert ss
	ss(1,1) =  a(pm,pm);   ss(2,2) =  a(pm-1,pm-1)
	ss(1,2) = -a(pm-1,pm); ss(2,1) = -a(pm,pm-1)
	ss = ss*det

! copy the last rows
	call dcopy(pm-2,a(pm-1,1),lda,z2(1,1),1)
	call dcopy(pm-2,a(pm,1),  lda,z2(1,2),1)

	call dgemm('N','T',pm-2,2,2,1.d0,z2,lda,ss,2,0.d0,w2,lda)

! copy the last columns (into z)
	call dcopy(pm-2,a(1,pm-1),1,z2(1,1),1)
	call dcopy(pm-2,a(1,pm),  1,z2(1,2),1)

! stripe out an outer product z.w^{T}
	call dgemm( 'N','T',pm-2,pm-2,2,-1.d0,z2,lda,w2,lda,1.d0,a,lda )

! update pm
	pm = pm - 2


	end subroutine inv_m2



!-------------------------------
!--- Det: change a single row
!-------------------------------
	real*8 function  det_rc(pm,lda,a,r,v,c,z,w)
	integer :: lda, pm       ! leading dimension and actual size of A
	real*8  :: a(lda,lda), v(lda), z(lda),w(lda)
	integer :: r,c

	real*8 :: la, la1,ttt
!
!  Input: a(lda,lda) --- inverse matrix
!         v(lda)     --- a row to be added
!         r          --- the row number 
!         z(lda)     --- a column to be added
!         c          --- the column number  
!         w(lda)     --- a working array
!
!           all arrays are of the size pm<lda 
!
!  Output:
!           det. ratio @ det_r

! \lambda =  ( v DOT A_r )
        la = 1.d0 + ddot(pm,v,1,a(1,r),1)

! \lambda_bar = ( A_c DOT z )
	la1 = 1.d0 + ddot(pm,z,1,a(c,1),lda)

	ttt=0.d0
	if( a(r,c) /= 0 )then

!  w = transp(a1)*v
	call dgemv('T',pm,pm,1.d0,a,lda,v,1,0.d0,w,1)

! ttt = (w DOT z )  * a(r,c)
	ttt = ddot(pm,w,1,z,1)* a(r,c)

	endif

	det_rc = la*la1  -ttt

	end function det_rc




!-------------------------------
!--- inv & det of a matrix (honest N**3)
!-------------------------------
	real*8 function  full_inv(pm,lda,a)
	integer :: lda, pm  ! leading dimension and actual size of A
	real*8  :: a(lda,lda)

	integer, allocatable :: ipiv(:)
	real*8, allocatable :: work(:)
	integer :: info, i,lwork,icnt
!
!  Input: A(lda,lda), pm<lda 
!
!  Output: A contains the inverse
!          full_inv contains the determinant

	lwork=lda
	allocate(ipiv(1:pm), work(1:lwork))

	call dgetrf(pm,pm,a,lda,ipiv, info) 

        full_inv = 1d0
         icnt=0
         do i=1,pm
           full_inv = full_inv * a(i,i)
           if (ipiv(i).ne.i) then
             icnt = icnt+1
           endif
         enddo
         if (mod(icnt,2).eq.1) full_inv = -full_inv



	call dgetri(pm,a,lda,ipiv,work,lwork,info )
	if(info/=0)print*,'dgetri info = ',info

	deallocate( work,ipiv )

	end function full_inv


	end module det_n2
