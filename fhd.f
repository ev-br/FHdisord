!
! Determinant diagrammatic MC code for the Fermi-Hubbard model
! at finite-temperature, 
!  using the updating strategy a la O. Goulko & M. Wingate, arxiv:0910.3909
!  
!  for details see
!  http://montecarlo.csi.cuny.edu/umass/fermion.html
!
!  In this code, the worm stuff is commented out, using the !w clause; Large worm-related 
!  functions and subroutines are deleted; these can be copy-pasted back from the original 
!  code from the reference above
!
!
!  THIS one is a rename of ./mess/fhd_jul18.f
!

!*******************************************
! TODO:
!  DONE 1.  get rid of infile etc --- make it a la EVK & Kris
!  DONE 2.  add blocking analysis to prnt()
!  DONE 3.  updates a la Goulko&Wingate
!  DONE 4.  worm measurement: N**2 correlated
!  DONE 5.  two sorts of measurements, cheap & expensive
!       5.  cleanup
!       6.  filling the matrix @ updates: a la add() or a la leap_add(), which one's faster?
!*******************************************



!---------------------
!--- Variables
!---------------------
	module vrbls
	implicit none

	real*8,  parameter :: nul=0.0d0, un1=1.0d0, un2 = 2.0d0
	real*8   :: pi
	integer, parameter :: i_hu = huge(1)

	integer, parameter :: nkink_m = 512  ! max # of kinks on one site
	integer, parameter :: nmnm_max = 8192 ! max total # of kinks

!------------------------------------------------------
!     Phys. parameters
!------------------------------------------------------

      integer, parameter :: d=3, dd=6       ! d is the global dimension; 
	                                      ! if you need d/=3, change it in here and in the parameter file  
      integer :: N(1:d)     ! Array for Numbers of sites per dimensions
      integer :: N2(1:d)
      integer :: Nsite      ! Total number of sites

	real*8 :: U, U_in    ! MINUS interaction strength, initial for thermalization
	real*8 :: beta       ! inverse temperature
	real*8 :: mu         ! chem. potential
	real*8, allocatable  :: musite(:)   ! site-dependent chem.potential

        real*8  ::   disord    ! disorder strength: \mu = \mu_av + (2r-1)*disord

	real*8 :: ef, kf     ! Fermi energy & momentum

!w        real*8 :: eta        ! GF- vs Z-sector weight

	real*8,parameter :: gf_eta=1.038d0   ! 1 + GF anomalous dimension, U(1) universality class

!--- Green function array
	integer :: Ntab                   ! Actual Number of sites per dimension for tabulation
	integer :: mtau

	real*8, allocatable :: GR_DAT(:,:,:)	! Green function array
	real*8, allocatable :: GRD_DAT(:,:,:)     ! Green function derivative
 
	real*8 :: bmt, bmt1     ! a shorthand for beta/mtau, and its inverse

	real*8 :: g0000         ! green function w/ all arg = 0, non-interacting density

!------------------------------------------------------
!     Sites and Associations
!------------------------------------------------------

      integer, allocatable :: ass(:,:)
      integer, allocatable :: back(:)
      integer, allocatable :: x(:,:)  ! coordinates x(1:d,site)


!     The array ass(...) specifies the nearest-neighbor sites in both positive and
!     negative i-direction. For example, in a 2-D lattice, like the one below,
!     the associations to the io site are:  
!     (here +(-)1 = positive (negative) x-direction, 
!     +(-)2 = positive (negative) y-direction)
      
!     ass(+1,i0)=i1;    ass(-1,i0)=i3;    ass(+2,i0)=i4;    ass(-2,i0)=i2
      
!                            ...    
!                             |  
!                     ... -- i2 -- ...
!                       |     |     |
!                ...-- i3 -- i0 -- i1 --.... 
!                       |     |     |
!                     ... -- i4 -- ... 
!                             |  
!                            ...
!           
!     * The array back(..) specifies the 'opposite direction', e.g. in 3D
!           back(1) = 4, back(4)=1 (positive x is 1, negative x is 4), and so forth
! 


!------------------------------------------------------
!     Update related
!------------------------------------------------------

	real*8 :: bun           !  = beta*U*Nsite 
!w	real*8 :: udt           !  = U*ddt_same

!w	integer :: n_ca, n_le             ! # of sites within the cre/ann, leap cube
!w	integer, allocatable :: uzli(:)   ! site list for cube()
!w	logical, allocatable :: yes(:)    ! flags for cube()

!w	integer :: rx_ca, rx_le            ! x- and \tau- radii for cre/ann, leaps
!w	real*8  :: rt_ca, rt_le

!w	integer, allocatable :: uzli_cube_ca(:,:), uzli_cube_le(:,:) ! site lists 

! mega-leaps
!w	integer :: nt_same             ! # of \tau points for leaps
!w	real*8  :: ddt_same            ! corresp. interval length
!w	integer :: n_ev_cut            ! # of events after the cutoff
!w	integer, allocatable :: dx_le(:,:)
!w	real*8,  allocatable :: dt_le(:), w_le(:)
!w	real*8 :: w_norm            ! normalization factor for w_le(:)

! nearest kink:
!w	integer, allocatable :: name_best(:)
!w	real*8, allocatable :: dist_best(:), distance(:) 


! leaps_o
        real*8 :: udt_o ;    ! U*dt_o, should be of order 1
        real*8 :: dt_o  ;    ! used in leaps_o


!------------------------------------------------------
!     Configuration
!------------------------------------------------------

! The data structure is purposedly redundant, for in all cases where 
! there is a memory/speed tradeoff, it is solved in favor of speed.
!
! The interaction vertices (a.k.a 'kinks')  live on the sites of a d-dimensional lattice
! and on the [0,\beta] interval of the \tau-continuum. The kinks are labelled by
! the 'names', e.g. unique numbers. The configuration consists of:
!
! 1. Overall list of kinks' names, managed by the name manager, see GetName() and DropName() subroutines
!    Given the kink name, one has access to the kink's site, temporal position (a.k.a 'time', also 'tau'),
!    and the position in TheMatrix, see below. 
!
! 2. Each site hosts a list of names of kinks which belong to this site. The order of kinks in these
!    lists is arbitrary.
!
! 3. The configuration weight matrix is not stored itself; one always deals with its inverse, 
!    a.k.a. TheMatrix :). 
!
! 4. ira/masha are not kinks. They are stored separately, see below.
!
	integer,   allocatable :: kname(:,:)   ! list of kinks's names @ a given site: kname(j,site)
	integer*2, allocatable :: nkink(:)     ! # of kinks @ site: nkink(site)

	integer, allocatable   :: ksite(:)      ! ksite(name) => site of a kink 'name'
	real*8, allocatable    :: ktau(:)       ! ktau(name) => tau of a kink 'name'

	integer, allocatable   :: row(:), clmn(:) ! row(name), clmn(name) => position of a 'name' in TheMatrix
	integer, allocatable   :: nm_row(:),nm_clmn(:)  ! nm_row(row) => name of the kink associated with the row

! NameManager data
        integer                        :: nmnm         ! total number of kinks
        INTEGER, DIMENSION(-2:nmnm_max) :: namelist, numnam

!w        logical :: present             ! =.true. is ira & masha are present in the configuration
!w        integer, parameter :: ira=-2, masha=-1

! TheMatrix
	integer             :: pm,lda      ! actual size & leading dimension
	real*8, allocatable :: matr(:,:)   ! TheMatrix, inverse
	real*8, allocatable :: m_u(:), m_v(:), m_w(:),m_z(:)       ! working arrays for rank-1 updates
	real*8              :: m_s

	real*8, allocatable :: m_u2(:,:),m_v2(:,:),m_s2(:,:),m_c2(:,:)  ! for rank-2 updates

!-------------------------------------------------------
!     Measurement related
!-------------------------------------------------------

!
!   Errorbar analysis is performed via blocking method. The maximum number of blocks is
!   fixed, once this number is reached, the blocks are collated so that one has twice
!   as few blocks of twice the size. 
!
	integer, parameter :: b_n_max=500   ! max # of blocks
	integer  :: Z_b                     ! block size

	integer :: b_n, i_b    ! # of filled blocks  & current nbr. of measurements 
	                       ! in the block into the (b_n+1)-th block

        integer :: i_b1, b_n1        ! same as i_b, but for diag_gim estimator 
        integer :: Z_b1

! Partition function
	real*8 :: Z            

! Diagonal estimators:
	real*8  :: PE, KE, ndens     ! PE, KE and density estimators 
	real*8  :: PE_stat(b_n_max), KE_stat(b_n_max)        ! files
	real*8  :: ndens_stat(b_n_max)

        real*8  :: gim, gim_stat(b_n_max)

	real*8  :: g_uu(0:100),g_ud(0:100)  ! Density-densiy correltaors, a.k.a 'dance-dance'


! Integrated ira-masha correlator:
!w	real*8  :: im               ! estimator
!w	real*8  :: im_stat(b_n_max) ! file
 
! Service histograms:
!w	integer, parameter :: im_t_max=15
!w	real*8  :: im_t(-im_t_max:im_t_max)          ! t(ira-masha) histogram

	integer, parameter :: nmnm_sup = 7500        ! kink # distribution
	real*8, dimension(0:nmnm_sup) :: nmnm_distr

	integer, parameter :: det_distr_max=25       ! det. distribution
	real*8 :: det_distr(-det_distr_max:det_distr_max)


!------------------------------------------------------
!      Flow control and output
!------------------------------------------------------

!
!  Since configuration weight determinant ratios (a.k.a.'det') are
!  recalculated incrementally, the roundoff error buildup is possible
!  Also, if the current configuration is close to the nodal surface,
!  the det-s might get extravagantly large/small, which also might lead to
!  error acumulation. In order to deal with these, determinants are recalculated 
!  from scratch (which costs N**3 operations!) (i) every so often (once per N**4 
!  proves to be OK), and (ii) every time unusually large det pops out (the latter
!  signals that the pre-update configuration was close to a nodal surface)
!
!
	real*8 :: tolerance    ! determinant recalculation tolerance level
	real*8 :: recalc_rate  ! recalc. rate counter

! Update manager probabilities
	real*8, dimension(1:10) :: prob

      real*8 :: step, step_p, step_w, step_m, step_m1, step_t, step_r
!     step    counts MC steps
!     step_p  the number of steps between printouts
!     step_w  the number of steps between writing to disk
!     step_m  the number of steps between measurements
!     step_t  the number of steps for thermolization
!     step_r  the number of steps for det. recalculation

!  Counters for printing, writing, measuring, thermalization, det recalculation
      real*8 :: i_p, i_w, i_m,i_m1, i_t, i_r

      real*8  :: n_sw   ! number of sweeps for thermalization, 1 sweep = \beta*U*Nsite updates

! Address/accept counters	for various updates
	real*8 :: c_a_v, a_a_v, c_d_v, a_d_v
!w	real*8 :: c_ann, a_ann, c_cre, a_cre
!w	real*8 :: c_l_a, a_l_a, c_l_d, a_l_d
	real*8 :: c_lo_a, a_lo_a, c_lo_d, a_lo_d
!w	real*8 :: c_r, a_r

	real*8 :: i_backsc     ! 'backscattering' counter: leap_add immediately followed by a leap_drop

! Kink number output
      integer :: nm_min=0, nm_max=i_hu
      real*8  :: nm_av=0

! 'status' variable holds the name of the update; proves useful for debugging
	character*6 :: status, prevstatus

! Output channel #:
	integer, parameter :: OUT_UNIT=22        ! For MPI code, writing to the stdout is useless, 
											 ! thus all the output goes to file(s)

! rndm() seed, see the 'mumbers' module
	integer :: r_ij, r_kl
	integer, allocatable :: ar_ij(:), ar_kl(:)


! Walltime constraints
	character*8  :: date               ! "how much on your watch?" -- "Six o'cloch"
	character*10 :: time
	character*5  :: zone

	integer :: tvalues(8)
      
	real*8 :: time_prev, time_curr, time_elapsed, time_limit
	integer :: hr_curr, hr_prev


! MPIzation related, see below:
	character*50 :: myparfname, fnamesuffix, outputfname, rndmstr
!	character*50 :: str,istr, anfname                     ! service vrbls


!!!	integer :: site1,site2, x1(d),x2(d)
!!!	real*8 :: t1,t2


	end module vrbls




!===============================================
!=== Main block :)
!===============================================
	program MAIN
	use vrbls; use det_n2
	use mumbers; !use tree
	implicit none 

	logical :: acpt                   ! accepted update flag
	integer :: ddd, i,j,i1, ans1, ans2, ans3 , dummy
	real*8  :: r, dp, det , dt
        character*6 :: cvoid


! ------- pi pi pi pi ------------
	pi = 4.d0*atan(1.d0)

!--------- Set the clock --------
	call date_and_time(date, time, zone, tvalues)
	time_curr = tvalues(5)*3600.d0 + tvalues(6)*60.d0 + tvalues(7)  ! seconds 
	hr_curr = tvalues(5)

	time_prev=0.d0; time_elapsed = 0.d0; hr_prev = 0


!--------- Get the par-file name from the command line --------------
        call  GETARG(1 , myparfname)

        ! determine the suffix for output files
        fnamesuffix = myparfname(4:)

        outputfname = 'out'//trim(fnamesuffix)   ! send output here


!------- Reading parameter file ------------------

        open(OUT_UNIT,file=outputfname,position='append')      ! to be open until the MC starts
        write(OUT_UNIT,*) ' reading ', myparfname,' ...'

        open( 1, file=trim(myparfname) )
      
        read(1,*) ddd        ! Global dimension
	   if(ddd/=d)then; print*,'d/=3'; call mystop; endif

        Nsite=1
        do i=1,d 
           read(1,*) N(i)  ! Number of sites per given dimension
           Nsite=Nsite*N(i)
	   N2(i)=N(i)/2
        end do
      
	read(1,*) ans1  ! 0 if new configuration, 1 if old one
	read(1,*) ans2  ! 0 if new statistics,    1 if old one
	read(1,*) ans3  ! 0 if new rndm() seed, 1 if read one
	
	read(1,*) mu            ! Chemical potential
        read(1,*) disord        ! disorder strength
	read(1,*) U,U_in        ! - U, interaction, & initial for thermalization
	read(1,*) beta          ! Inverse temperature
	    bun = beta*U*Nsite  ! Shorthand
	    step_r = bun*bun    ! step for determinant recalculation
		
!w	read(1,*) eta          ! G- vs. Z-sector weight
	read(1,*) mtau         ! # of points in \tau for tabulation of GFs
						   ! Again, for inlining reasons it's necessary to have static variables here
!	   if(mtau>mt_max)then; print*,'mtau > mtau_max'; call mystop; 
!	   endif               
	   bmt = beta/mtau; bmt1 = 1.d0/bmt                   ! shorthands
	read(1,*) tolerance    ! det recalculation tolerance
        read(1,*) udt_o        ! used in leap_o updates
            dt_o = udt_o/U
!w	read(1,*) rx_ca, rt_ca ! cre/ann x- and \tau- radii
!w	    if(2*rt_ca>beta)then; print*,'rt_ca > beta /2'; call mystop; 
!w		endif
!w	read(1,*) rx_le, rt_le ! cre/ann x- and \tau- radii
!w	    if(2*rt_le>beta)then; print*,'rt_le > beta /2'; call mystop; 
!w		endif 
!w	read(1,*) nt_same              ! number of points for tabulation of weights for leaps
!w	    ddt_same = rt_le/nt_same   ! shorthands
!w	    udt=U*ddt_same
      read(1,*) n_sw                 ! number of sweeps for thermalization
      read(1,*) step_p               ! step for printing
      read(1,*) step_w               ! step for writing to disk
      read(1,*) step_m               ! step for measuring
      read(1,*) step_m1              ! step for measuring expensive updates
	read(1,*) time_limit           ! time limit [hrs]
 	    time_limit = time_limit*3600   ! seconds      
	read(1,*) cvoid                ! horizontal separator :)
! Update manager probabilities: 
      read(1,*) prob(1); prob(2)=prob(1)+prob(1)                ! add/drop
      read(1,*) dp; prob(3)=prob(2)+dp;  prob(4)=prob(3)+dp     ! leap_o_add/drop
!w      read(1,*) dp; prob(5)=prob(4)+dp;  prob(6)=prob(5)+dp     ! cre/ann
!w      read(1,*) dp; prob(7)=prob(6)+dp;  prob(8)=prob(7)+dp     ! leap_add/drop
!w      read(1,*) dp; prob(9)=prob(8)+dp                          ! hop
	if(abs(prob(4)-1.d0)>1.0d-10)then
	   print*,'Update manager probs. normalization problem...'
	   call mystop
	endif

!--- Random number generator:
!
!  Number of seeds provided MUST be equal to the groupsize (one seed per process)
!  Otherwise the program crashes with mysterious error message
!
!
       if(ans3==0)then     ! init rndm() state

          read(1,*)cvoid


          read(1,*) r_ij,r_kl
          call init_rndm(r_ij,r_kl)

        else                ! read rndm() state

          rndmstr='rndm'//trim(fnamesuffix)//'.dat'
          call read_rndm(rndmstr)

        endif



      close(1)

	write(OUT_UNIT,*)'....done'


! leaps_o : FIXME
!        udt_o = 0.7d0 ; 
!        dt_o = udt_o/U ; 



!---------- DONE reading, allocating memory:

!--- Lattice
	allocate (ass(dd,Nsite),back(dd),x(1:d,1:Nsite))
      CALL ASSA

!--- Green functions
	allocate( GR_DAT(0:mtau+1,1:Nsite,1:Nsite), GRD_DAT(0:mtau+1,1:Nsite,1:Nsite) )

!--- site-dependent chemical potential
	allocate( musite(1:Nsite) )

!--- Configuration
	allocate( nkink(1:Nsite), kname(nkink_m,1:Nsite) )  ! site data
	allocate( ksite(-2:nmnm_max),ktau(-2:nmnm_max) )    ! global lists

	allocate( row(-2:nmnm_max),clmn(-2:nmnm_max) )      ! c/r <-> name links
	
	lda=128; allocate(nm_row(lda),nm_clmn(lda))         ! Allocate 128*128 matrix first. 
														! more memory will be allocated later if necessary	
	allocate(matr(1:lda,1:lda),m_u(1:lda),m_v(1:lda))   ! TheMatrix
	allocate(m_w(1:lda),m_z(1:lda))

	allocate(m_u2(lda,2),m_v2(2,lda),m_s2(2,2),m_c2(lda,2) )
 

! Name Manager is initialized in either init_cnf or rd_cnf

! Site lists that contain sites within cubes with given edge lengths around given site
! it proves convenient to have TWO sets of lists: one for leap-type updates, one for cre/ann updates 
!w	allocate( uzli(1:Nsite),yes(1:Nsite) )
!w	n_ca=cube(1,rx_ca,uzli); allocate(uzli_cube_ca(1:n_ca,1:Nsite))
!w	 call tab_cube(rx_ca,n_ca,uzli_cube_ca)
!w	n_le=cube(1,rx_le,uzli); allocate(uzli_cube_le(1:n_le,1:Nsite))
!w	 call tab_cube(rx_le,n_le,uzli_cube_le)
!w	deallocate(uzli,yes)

! Working arrays for the closest neighbor selection in leap_drop.
!w	allocate( name_best(1:n_le),dist_best(1:n_le) )
!w	allocate( distance(1:nkink_m) )

!------------ DONE with memory allocation, proceed to tabulations

!--- Tabulate Green functions
	write(OUT_UNIT,*)'tabulating GFs...'  
	Ntab=maxval(N(:)/2+1)                          ! the longest jump w/PBC is N/2
!	if(ntab>ntab_max)then; write(OUT_UNIT,*)'Ntab_max is too small'; 
!	call mystop; endif
!	if(mtau>mt_max)then; write(OUT_UNIT,*)'mt_max is too small'; 
!	call mystop; endif
	call TABULATE 
	write(OUT_UNIT,*)'...done'

!+++++++++++++++++++++++++++++++++++++++ @#$   FOR TESTING PURPOSES ONLY
!	site1=1; x1=x(:,site1)
!	site2=site1; !site2=ass(1,site2); site2=ass(3,site2); 
!	x2=x(:,site2)
!	t1=0;
!	open(1,file='g0_fh.txt') 
!	do j=-9,9
!	   t2= beta*j/10.d0
!	   write(1,*),t2,GREENFUN(x1,t1,x2,t2)
!	   print*,j,t2,GREENFUN(x1,t1,x2,t2)
!	enddo
!	close(1)
!	pause
!+++++++++++++++++++++++++++++++++++++++


	
!--- Tabulate the weights for leap_add/drop updates
!w	call tab_leaps

!--- Structure output
      nm_av=0; nm_min=i_hu; nm_max=0
      nmnm_distr=0; det_distr=0.d0; 


! 'how much on your watch?' a.k.a initialize the clock
	time_prev = time_curr; hr_prev = hr_curr

	call date_and_time(date, time, zone, tvalues)
	time_curr = tvalues(5)*3600.d0 + tvalues(6)*60.d0 + tvalues(7)  ! seconds 
	hr_curr = tvalues(5)

	dt = time_curr - time_prev
	if( hr_curr < hr_prev )dt = dt + 24*3600.d0   ! take care of the runs that span across midnight 


	write(OUT_UNIT,*)' '
	write(OUT_UNIT,*)' '
	write(OUT_UNIT,*)' init. done: ',dt,' sec'	
	
!--- close output
	close(OUT_UNIT)


!----------- DONE with initializations



!=========== Read / init conf & thermalize =========

      if (ans1 == 1) then
	   call RD_CNF         ! read the condiguration
	   call therm2         ! thermalize if requested (via worm scheme)
      else
         call INIT_CNF       ! init the configuration 
	   call therm1         ! thermalize (via diagonal scheme)
      end if

! nullify counters
      i_p=0.d0; i_w=0.d0; i_m=0.d0; i_m1=0.d0; step=0.d0; 
      c_a_v=0.d0; a_a_v=0.d0; c_d_v=0.d0; a_d_v=0.d0
!w      c_cre=0.d0; a_cre=0.d0; c_ann=0.d0; a_ann=0.d0
!w      c_l_a=0.d0; a_l_a=0.d0; c_l_d=0.d0; a_l_d=0.d0
      c_lo_a=0.d0; a_lo_a=0.d0; c_lo_d=0.d0; a_lo_d=0.d0
!w      c_r=0.d0; a_r=0.d0
      recalc_rate=0.d0

! Construct TheMatrix
	det=recalc_matrix(pm,lda,matr)


!============== Read / init  statistics ==============
      if (ans2 == 1) then
         call RD_STAT
      else
         call INIT_STAT
      end if


!==================== main MC loop  ===============
	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'  '
	write(OUT_UNIT,*)' Start MC loop'
	close(OUT_UNIT)

	do;

         step=step+un1                          
         i_p=i_p+un1; i_w=i_w+un1; i_m=i_m+un1; i_m1=i_m1+1.d0;  i_r=i_r+1.d0


!------------------------ update --------------------------
	   r = rndm()
	   if     (r<prob(1))then;    call add(acpt,det)
	   else if(r<prob(2))then;    call drop(acpt,det)
	   else if(r<prob(3))then;    call leap_o_add(acpt,det)
	   else if(r<prob(4))then;    call leap_o_drop(acpt,det)
!w	   else if(r<prob(5))then;    call create(acpt,det)
!w	   else if(r<prob(6))then;    call annihilate(acpt,det)
!w	   else if(r<prob(7))then;    call leap_add(acpt,det)
!w	   else if(r<prob(8))then;    call leap_drop(acpt,det)
!w	   else if(r<=prob(9))then;   call hop(acpt,det)
	   else; print*,'am I nuts or what???'; call mystop
	   endif

!------------- recalculate if necessary -------------
	   if( acpt ) then

! det. distribution
!            i = floor(log10(abs(det)))
!	      if(abs(i)<det_distr_max) det_distr(i) = det_distr(i)+1.d0

! recalculate
	      if(abs(det)>tolerance) det=recalc_matrix(pm,lda,matr)

	   endif

	   if (i_r == step_r) then; i_r=0; 
	        det=recalc_matrix(pm,lda,matr); 
	   endif
!-----------------------------------------------------


         nm_av=nm_av+(un1*nmnm)/step_p        
         IF(nmnm>nm_max)nm_max=nmnm
         IF(nmnm<nm_min)nm_min=nmnm

	 if(nmnm<nmnm_sup) nmnm_distr(nmnm) = nmnm_distr(nmnm) + 1.d0


         if (i_m == step_m)  then; i_m=nul; call MEASURE; end if
         if (i_m1== step_m1) then; i_m1=nul; call MEASURE_1; endif
         if (i_p == step_p)  then; i_p=nul; call PRNT;   end if  
         if (i_w == step_w)  then;  i_w=nul; call WRT;     
	   end if 

	enddo


	contains  

!*************************  subroutines below ONLY ************************

!-----------------
!--- Add a kink
!-----------------
	subroutine add(acpt,det)
	logical :: acpt             ! a flag to be returned to the main loop
	real*8  :: det              ! det value to be returned to the main loop

	real*8  :: ratio
	integer :: name, site,j,nk, xnew(d), xi(d),xm(d),xv(d),vova, sv
	real*8  :: tnew, ti, tm, tv



   !     if(nmnm>0)return


	c_a_v = c_a_v + un1; acpt=.false.
	prevstatus=status; status='_add__' 

	if(pm==lda)then; call resize_matrix(-1); endif


!------- check
	if(nmnm>=nmnm_max)then;
	   open(OUT_UNIT,file=outputfname,position='append')
	   write(OUT_UNIT,*)'add: nmnm > nmnm_max', nmnm,nmnm_max, step
	   close(OUT_UNIT)
	   call mystop
	endif
!---------- end check

	site=Nsite*rndm()+1.d0
	if(site>Nsite)site=Nsite    ! select a site to play on
	xnew(:) = x(:,site) 
	tnew = beta*rndm()                    ! select tau		

!------------- determinant
	if(pm==0)then; det=g0000  !ratio = g0000**2
	else

!w! ira-masha 
!w      if(present)then
!w          xm=x(:,ksite(masha)); tm=ktau(masha)
!w          m_u(row(masha)) = GREENFUN(site,tnew,ksite(masha),tm)
!w          xi=x(:,ksite(ira)); ti=ktau(ira)
!W          m_v(clmn(ira)) = GREENFUN(ksite(ira),ti,site,tnew)
!w      endif

! calcualte the det ratio
	do j=1,nmnm; vova=namelist(j);
	   sv=ksite(vova) 
	   xv=x(:,sv); tv=ktau(vova); 
	   m_v(clmn(vova)) = GREENFUN(sv,tv,site,tnew)
	   m_u(row(vova))  = GREENFUN(site,tnew,sv,tv)
	enddo
	m_s = g0000
	det = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)   ! det ratio itself

	endif  ! pm==0
!---------------------------

	ratio = det*det
	ratio = ratio*bun/(nmnm+1)


! Metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif

	acpt=.true.

! update configuration 

!	call CheckGetName;    ! this slows things down a lot 

	call GetName(name)                               ! get a name for a new kink 
	nk = nkink(site)
	nkink(site) = nk+1; kname(nk+1,site) = name      ! on-site info
	ksite(name)=site; ktau(name)=tnew                ! global list

	if(pm==0)then                                    ! update TheMatrix
	   pm=1; matr(pm,pm)=1.d0/g0000
	else
	   call inv_p1(pm,lda,matr,det,m_v,m_w,m_z)        
	endif

	clmn(name) = pm; nm_clmn(pm) = name              ! update row/column <-> kink name links
	row(name)  = pm; nm_row(pm)  = name; 

	a_a_v = a_a_v + un1

	end subroutine add



!-----------------
!--- Drop a kink
!-----------------
	subroutine drop(acpt,det)
	logical :: acpt
	real*8  :: det
!
! This is complementary to subroutine add above
!
	real*8 :: ratio 
	integer :: site,j,nk, name, num,r,c, vova

	c_d_v = c_d_v + un1;   acpt=.false.
	prevstatus=status; status='_drop_'

	if(nmnm==0)return    ! nothing to drop yet :)
 
	num=nmnm*rndm()+1.d0                           ! play a kink to be dropped
	if(num>nmnm)num=nmnm; name=namelist(num); 
	site=ksite(name); nk = nkink(site)
	
!---------- determinant
	if(pm==1)then; det=1.d0/matr(pm,pm)
	else
	    r = row(name); c=clmn(name)
	    det = det_m1(pm,lda,matr,r,c)
	endif
!----------------------

	ratio = det*det*bun/nmnm
	ratio = un1/ratio

! Metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif

	acpt=.true.

! update configuration 
	if(pm==1)then; pm=0; nkink(site)=0
	else

	   do j=1,nk; if(kname(j,site)==name)exit   ! find name
	   enddo

	   kname(j,site) = kname(nk,site); nkink(site) = nk - 1    ! on-site

	   vova = nm_row(pm); row(vova) = r; nm_row(r) = vova      ! matrix links
	   vova = nm_clmn(pm); clmn(vova) = c; nm_clmn(c) = vova
	   call inv_m1(pm,lda,matr,r,c)                            ! TheMatrix

	endif

!	call CheckDropName(name);
	call DropName(name)

	a_d_v = a_d_v + un1

	end subroutine drop






!----------------------------
!--- Leap masha & add a kink : a la Goulko & Wingate
!----------------------------
	subroutine leap_o_add(acpt,det)
	logical :: acpt
	real*8  :: det

        integer :: num, site, lesha, nk, vova, name, j, sv
        real*8  :: tau, tlesha, dist_name, dist_vova, ratio, tv

!        if(nmnm>0)return

	if(pm==lda)then; call resize_matrix(-1); endif

!------- check --------------------------------------
	if(nmnm>=nmnm_max)then;                 ! too many kinks
	   open(OUT_UNIT,file=outputfname,position='append')
	   write(OUT_UNIT,*)'leap_add: nmnm > nmnm_max', nmnm,nmnm_max,step
	   close(OUT_UNIT)
	   call mystop
	endif
!---------- end check ----------------------------------

	c_lo_a = c_lo_a + 1.d0


! play lesha: adding will be in its neighborhood
        num=nmnm*rndm()+1.d0    
	if(num>nmnm)num=nmnm; lesha=namelist(num); 
	site=ksite(lesha); 
        nk = nkink(site)
 
!---------- check ----------------
	if(nkink(site)==nkink_m)then           ! too many kinks on the particular site
	   open(OUT_UNIT,file=outputfname,position='append')
	   write(OUT_UNIT,*)'leap_add: nk(lesha)= ',nkink(site), step
	   close(OUT_UNIT)
	   call mystop
	endif
!-----------------------------



! play where to add
        tlesha= ktau(lesha) ;
        tau = tlesha + dt_o*rndm() ;   ! forward in time only
        if(tau>beta)tau=tau-beta


! check detailed balance: the new one must be the nearest to lesha
        if(nk>1)then   
            vova = find_nearest(lesha);
            dist_vova = time_distance(tlesha,ktau(vova))
            dist_name = time_distance(tlesha,tau)
            if( dist_name > dist_vova )return
        endif
!---------------------
 

!------------- determinant
	if(pm==0)then; det=g0000 
	else

! calcualte the det ratio
	  do j=1,nmnm; vova=namelist(j);
	     sv=ksite(vova) ;   tv=ktau(vova) 
	     m_v(clmn(vova)) = GREENFUN(sv,tv,site,tau)
	     m_u(row(vova))  = GREENFUN(site,tau,sv,tv)
	  enddo
	  m_s = g0000
	  det = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)   ! det ratio itself

	endif  ! pm==0
!---------------------------


	ratio = det*det*udt_o*nmnm/(nmnm+1.d0)

! Metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif

	acpt=.true.

! update configuration 
!	call CheckGetName;    ! this slows things down a lot 

	call GetName(name)                               ! get a name for a new kink 
	nkink(site) = nk+1; kname(nk+1,site) = name      ! on-site info
	ksite(name)=site; ktau(name)=tau               ! global list

	if(pm==0)then                                    ! update TheMatrix
	   pm=1; matr(pm,pm)=1.d0/g0000
	else
	   call inv_p1(pm,lda,matr,det,m_v,m_w,m_z)        
	endif

	clmn(name) = pm; nm_clmn(pm) = name              ! update row/column <-> kink name links
	row(name)  = pm; nm_row(pm)  = name; 

	a_lo_a = a_lo_a + 1.d0



!        print*,'leap_add: ', step
!        print*,'lesha = ', lesha, site, tlesha
!        print*,'name  = ', name, ksite(name), ktau(name)
!        print*, '@site : ', kname(1:nk+1,site)
!        print*,'++++++++++++++'; print*
!
        call try_leap_o_drop(acpt, det, ratio, lesha)


	end subroutine leap_o_add


!-----------------
!--- Leap masha & drop a kink : after Goulko & Wingate
!-----------------
	subroutine try_leap_o_drop(acpt,det, ratio_, lesha)
	logical :: acpt
	real*8  :: det, ratio_

        integer :: num, lesha, name, vova, nk, site, r,c, j
        real*8  :: tlesha, ratio


	if(nmnm==0)return          ! nothing to drop yet :)

! play lesha [the nearest one will be dropped]
!        num=nmnm*rndm()+1.d0    
!	if(num>nmnm)num=nmnm; lesha=namelist(num); 
	site=ksite(lesha); 
        tlesha= ktau(lesha) ; 
        nk = nkink(site)

!------------------- find the name of the kink nearest to lesha
        if(nk<2) return    ! can't do it
        name= find_nearest(lesha)
        if( time_distance(tlesha, ktau(name))>dt_o )return   ! detailed balance
!------------------------------------


 
!---------- determinant
	r = row(name); c=clmn(name)
	det = det_m1(pm,lda,matr,r,c)
!------------------------

	ratio =  det*det*udt_o*(nmnm-1)/nmnm
	!ratio = 1.d0/ratio

        if( abs(ratio-ratio_)>1d-10 )then

            print* ; print*
            print*,"try_drop: ", step
            print*,'ratio = ', ratio
            print*,'ratio_= ', ratio_

            print*
            print*,'lesha = ', lesha, ksite(lesha), ktau(lesha)
            print*, 'name =  ', name, ksite(name), ktau(name)
            print*, '@site : ', kname(1:nk,site)
         
           pause
        endif

	end subroutine try_leap_o_drop



!-----------------
!--- Leap masha & drop a kink : after Goulko & Wingate
!-----------------
	subroutine leap_o_drop(acpt,det)
	logical :: acpt
	real*8  :: det

        integer :: num, lesha, name, vova, nk, site, r,c, j
        real*8  :: tlesha, ratio

	status='leap_d'; acpt=.false.; 

	if(nmnm==0)return          ! nothing to drop yet :)

	c_lo_d = c_lo_d + 1.d0

! play lesha [the nearest one will be dropped]
        num=nmnm*rndm()+1.d0    
	if(num>nmnm)num=nmnm; lesha=namelist(num); 
	site=ksite(lesha); 
        tlesha= ktau(lesha) ; 
        nk = nkink(site)

!------------------- find the name of the kink nearest to lesha
        if(nk<2) return    ! can't do it
        name= find_nearest(lesha)
        if( time_distance(tlesha, ktau(name))>dt_o )return   ! detailed balance
!------------------------------------


!---------- determinant
	r = row(name); c=clmn(name)
	det = det_m1(pm,lda,matr,r,c)
!------------------------

	ratio =  det*det*udt_o*(nmnm-1)/nmnm
	ratio = 1.d0/ratio


!------------ Metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif

	acpt=.true.

!---------- update configuration 

! update configuration 
	if(pm==1)then; pm=0; nkink(site)=0
	else

	   do j=1,nk; if(kname(j,site)==name)exit   ! find name
	   enddo

	   kname(j,site) = kname(nk,site); nkink(site) = nk - 1    ! on-site

	   vova = nm_row(pm); row(vova) = r; nm_row(r) = vova      ! matrix links
	   vova = nm_clmn(pm); clmn(vova) = c; nm_clmn(c) = vova
	   call inv_m1(pm,lda,matr,r,c)                            ! TheMatrix

	endif

!	call CheckDropName(name);
	call DropName(name)

	a_lo_d = a_lo_d + 1.d0


	end subroutine leap_o_drop



!-----------------------------
!--- distance from 1st arg to 2nd arg
!-----------------------------
        real*8 function time_distance(t0,t1)
        real*8  :: t0,t1
        real*8 :: dt
!c	dt = abs(t1-t0); dt = min(dt, beta-dt)
!c	if( dt<=rt_le )then

	dt=t1-t0; if(dt<0)dt=dt+beta
        time_distance=dt ; 

	end function time_distance


!-----------------------------
!--- find the kink nearest to arg
!-----------------------------
        integer function find_nearest(lesha)
        integer, intent(in) :: lesha
        integer :: site, nk, j, vova, name
        real*8 :: dist, dist_best, t_lesha        

        site = ksite(lesha)
        nk=nkink(site);  if(nk<2)then; print*,' find_nearest called @site= ',site,' w/ nk= ', nk ; stop; endif
        t_lesha = ktau(lesha)

        dist_best=1d200; name=-101 ; ! this will crash the calling routine if things go nuts
        do j=1,nk
            vova=kname(j,site);
            if( vova /= lesha )then
               dist = time_distance( t_lesha, ktau(vova) ); 
               if( dist<dist_best )then
                   name=vova;
                   dist_best = dist
               endif
            endif
        enddo

        find_nearest= name 

	end function find_nearest







!------------------
!--- Measurements
!------------------
	subroutine measure

	integer :: j
	real*8  :: eta1

!w	if(present)then       !   GF sector
!w
!w	    j = floor((ktau(masha)-ktau(ira))*im_t_max/beta)    
!w	    im_t(j) = im_t(j)+1.d0                       ! t(ira-masha) histogram

!w            im = im + 1.d0                               ! integrated correlator

!w	else                  !   Z sector

	    i_b = i_b + 1 ; 	Z = Z + un1
	    PE = PE + 1.d0*nmnm/beta               ! PE
	    KE = KE + diag_KE()                    ! KE
	    ndens = ndens + diag_dens()            ! number density 
!!	    call dance_dance                       ! dance-dance correlators


!w	endif  ! present


!-------------- is a current block full?
	if( i_b == Z_b )then    ! wrap it up
	        b_n = b_n + 1; i_b=0
		PE_stat(b_n) = PE; PE = 0.d0
		KE_stat(b_n) = KE; KE = 0.d0
		ndens_stat(b_n) = ndens; ndens = 0.d0


!w! normalize [see create()] & rescale via *L^{gf_eta}/(beta*Nsite)**2
!w		eta1 = n_ca*2.d0*rt_ca/eta
!w	        eta1 = eta1 * N(1)**(gf_eta)/beta/Nsite
!w		im_stat(b_n) = im*eta1; im = 0.d0      ! normalize via eta1*Z

		if(b_n == b_n_max)then                  
!w                    call collate( im_stat, b_n)
		    call collate( PE_stat, b_n )
                    call collate( KE_stat, b_n )
                    call collate( ndens_stat, b_n )
                    b_n = b_n/2; Z_b = Z_b * 2.d0
		endif
	endif



	end subroutine measure


!------------------
!--- Measurements: expensive
!------------------
	subroutine measure_1

        i_b1 = i_b1+1
        gim = gim + diag_gim();

!-------------- is a current block full?
        if( i_b1==Z_b1 )then  ! wrap it up too
              b_n1=b_n1 + 1; i_b1=0
              gim_stat(b_n1) = gim;  gim=0.d0

              if(b_n1==b_n_max)then  
                   call collate( gim_stat, b_n1 )
                   b_n1 = b_n1/2; Z_b1 = Z_b1*2.d0
              endif


        endif
        

        end subroutine measure_1



!--------------------------------
!--- measure density: insert an extra kink @ random place
!--------------------------------
	real*8 function diag_dens()

	real*8 :: det,ti, tn
	integer :: site, xn(d), xi(d), j, vova 

! select where to insert a kink
	site=Nsite*rndm()+1.d0; if(site>Nsite)site=Nsite
	xn(:)=x(:,site); tn=beta*rndm()
	
!------------- determinant
	if(nmnm==0)then; det=g0000  
	else

	  do j=1,pm; vova=nm_row(j)               ! fill a column
	     ti=ktau(vova)
	     m_u(j)  = GREENFUN(site,tn,ksite(vova),ti)
	  enddo

	  do j=1,pm; vova=nm_clmn(j)             ! fill a row
	     ti=ktau(vova)
	     m_v(j) = GREENFUN(ksite(vova),ti,site,tn)
	  enddo

	  m_s = g0000

	  det = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)   ! det ratio itself

	endif  ! nmnm==0

! estimator
	diag_dens = det*2.d0   ! a factor of 2 reflects the spin summation

	end function diag_dens



!--------------------------------
!--- measure density: remove 1/2 kink (nmnm times per call)
!--------------------------------
	real*8 function diag_dens_m1_()

        integer :: num, name, r,c,j
        real*8  :: det, sum, dummy


        if(nmnm==0)then; dummy=0.d0;
        else; 

           sum=0.d0
           do num=1,nmnm
              name=namelist(num); 
              r=row(name); c=clmn(name)
              if(pm==1)then; det=1.d0/matr(pm,pm);  
	      else
	         det = det_m1(pm,lda,matr,r,c)
	      endif
              sum = sum + 1.d0/det
           enddo 
           sum = sum/nmnm  

           dummy = 2.d0*sum*nmnm/(U*beta*Nsite)

           diag_dens_m1_=dummy

        endif


	end function diag_dens_m1_


!--------------------------------
!--- measure density: remove 1/2 kink  (just one kink per call)
!--------------------------------
	real*8 function diag_dens_m1_oneshot()

        integer :: num, name, r,c,j
        real*8  :: det, sum


        if(nmnm==0)then; diag_dens_m1_oneshot=0.d0;
        else; 

	   num=nmnm*rndm()+1.d0 ; if(num>nmnm)num=nmnm; 
           name=namelist(num); 

           r=row(name); c=clmn(name)
           if(pm==1)then; det=1.d0/matr(pm,pm);   ! print*,'det=',det; pause
	   else
	         det = det_m1(pm,lda,matr,r,c)
	   endif

           det = 1.d0/ det

           diag_dens_m1_oneshot = 2.d0*det*nmnm/(U*beta*Nsite)
        endif


	end function diag_dens_m1_oneshot


!--------------------------------
!--- measure g_im: remove 2 1/2-kinks (all pairs of kinks per call)
!--------------------------------
	real*8 function diag_gim()

        integer :: name1,name2, r1,c1,r2,c2,i,j
        real*8  :: det1,det2, sum, dummy


        if(nmnm<2)then; sum=0.d0;
        else; 

           sum=0.d0
           do i=1,nmnm; name1=namelist(i)
              do j=1,nmnm; name2=namelist(j)
                if(name1/=name2)then

                  r1=row(name1); 
                  c2=clmn(name2)

                  det1 = 1.d0/det_m1(pm,lda,matr,r1,c2)
                  det2 = det1
                  
                  sum = sum + det1*det2  

                endif 
              enddo
           enddo

           !sum = sum/ nmnm*(nmnm-1)                   ! nmnm*(nmnm-1) measurements were made
           !sum=sum*nmnm*(nmnm-1)/(U*beta*Nsite)**2    ! normalization
           sum = sum/(U*beta*Nsite)**2
        endif

        diag_gim = sum* N(1)**(gf_eta)

	end function diag_gim



!--------------------------------
!--- measure g_im: remove 2 1/2-kinks (one shot per call)
!--------------------------------
	real*8 function diag_gim_oneshot()

        integer :: num, name1,name2, r1,c1,r2,c2,j
        real*8  :: det1,det2, sum, dummy


        if(nmnm<2)then; diag_gim_oneshot=0.d0;
        else; 

	   num=nmnm*rndm()+1.d0 ; if(num>nmnm)num=nmnm; 
           name1=namelist(num); 

111        continue
           num=nmnm*rndm()+1.d0; if(num>nmnm)num=nmnm
           name2=namelist(num)
           if(name2==name1) goto 111 ; 

           r1=row(name1); c1=clmn(name1);
           r2=row(name2); c2=clmn(name2);  

	   det1 = 1.d0/det_m1(pm,lda,matr,r1,c2)
           det2 = det1  !1.d0/det_m1(pm,lda,matr,r2,c1)

           dummy= det1*det2*nmnm*(nmnm-1)/(U*beta*Nsite)**2

           diag_gim_oneshot = dummy* N(1)**(gf_eta)

        endif


	end function diag_gim_oneshot





!-----------------------------------------------
!--- measure kinetic energy from the Z-sector
!-----------------------------------------------
	real*8 function diag_KE()

	integer :: site1, site2, x1(d),x2(d), vova, xv(d), j, i
	real*8  :: t, det,tv
!
!  This estmator is for the tight-binding dispersion ONLY (n.n. sites)
!
	site1=Nsite*rndm()+1.d0; if(site1>Nsite)site1=Nsite; 

	i=dd*rndm()+1; if(i>dd)i=dd; site2 = ass(i,site1)
		
	x1=x(:,site1); x2=x(:,site2)
	t=beta*rndm()
	

!------------- determinant
!	m_s = GREENFUN(x2,t,x1,t)        ! a corner 
	m_s = GREENFUN(site2,t,site1,t)

	if(pm==0) then; det = m_s        
	else

	  do j=1,pm; vova=nm_row(j)               ! fill a column
	     tv=ktau(vova)
	     m_u(j)  = GREENFUN(site2,t,ksite(vova),tv)
	  enddo

	  do j=1,pm; vova=nm_clmn(j)             ! fill a row
	     tv=ktau(vova)
	     m_v(j) = GREENFUN(ksite(vova),tv,site1,t)
	  enddo

	  det = det_p1(pm,lda,matr,m_u,m_v,m_z,m_s)   ! det ratio itself

	endif

! return value
	diag_KE = dd*Nsite*det  ! d*Nsite is the # of bonds,       
				! a factor of 2 is due to the spin summation ---- NO !!!!
				! ATTENTION : this is the energy PER COMPONENT
				! NEED to multiply by 2 to get total kin. energy
				! (reason: Ekin(one component) = 2t*det --- )


	end function diag_KE


!================  DONE with measurements; various service functions below 



!------------------------
!--- Collate the array
!------------------------
      subroutine collate(arr,n)
      real*8, dimension(*) :: arr
      integer :: n, Zb        ! array length & # of measurements per array entry

      integer :: i

      do i=1,n/2
          arr(i)=arr(2*i)+arr(2*i-1)
      enddo

      end subroutine collate


!-------------------------------
!--- Analyze block statistics: average and dispersion
!-------------------------------
	subroutine bSTAT(arr,n,Zb,av,err)
	real*8, dimension(*) :: arr
	integer              :: n, Zb
	real*8               :: av, err

	real*8 :: av2, dif

	av  = sum( arr(1:n) )/Zb/n
	av2 = sum( arr(1:n)**2 )/Zb/Zb/n

				!av2 = av2 + (arr(j)/Zb)**2

	dif = av2 - av*av; 	if(dif<0.d0)dif=0.d0

	err = sqrt( dif ) / sqrt(1.d0*n)


	end subroutine bSTAT


!-------------------------------
!--- Merge blocks & emit av +/- err
!-------------------------------
        subroutine mrg(arr,n,Zb)
        integer, intent(in)              :: n, Zb
        real*8, dimension(1:n), intent(in) :: arr

        real*8  :: av, err, arr1(1:n)
        integer :: i, n1, zb1


        zb1 = zb;       arr1(1:n) = arr(1:n); n1=n

! assuming that the unit 1 is already open by the calling routine
        write(1,*),'-----------'   

        do;

! emit
      call bSTAT(arr1,n1,zb1,av,err)
      write(1,777) av, err,n1
 777  format(4x,g12.5,4x,' +/- ',g12.5,8x,I3)

! enough?
        if(n1<3)exit

! merge
        n1=INT(n1/2); zb1=zb1*2
        do i=1,n1
            arr1(i) =  arr1(2*i-1) + arr1(2*i)
        enddo

        enddo

        write(1,*),'------------'; write(1,*)' ' ;


        end subroutine mrg






!------------------------------------
!--- Print out and check the runtime
!------------------------------------
	subroutine prnt
	integer :: i,j, name
	real*8 :: xxx, yyy, dt 

	real*8 :: PE_av, PE_err, KE_av, KE_err, im_av, im_err
	real*8 :: dens_av, dens_err

	logical :: lastone


! 'how much on your watch?'
	time_prev = time_curr; hr_prev = hr_curr

	call date_and_time(date, time, zone, tvalues)
	time_curr = tvalues(5)*3600.d0 + tvalues(6)*60.d0 + tvalues(7)  ! seconds 
	hr_curr = tvalues(5)

	dt = time_curr - time_prev
	if( hr_curr < hr_prev )dt = dt + 24*3600.d0   ! across the midnight

	time_elapsed = time_elapsed + dt

	lastone=.false.
	if( time_elapsed > time_limit )  lastone=.true.   ! time to wrap up
!-------------------------------------------

	open(1,file=outputfname,position='append')

	write(1,*)''
        write(1,*)'-------------------------', time_elapsed/3600,' hrs'
	
	if(i_t<step_t)then
	  write(1,*)' thermalization: ',1.d2*i_t/step_t,' % done'
	  write(1,*)' current U = ', bun/beta/Nsite
	  write(1,*)'  '
	endif

! is there enough statistics?
	if(b_n>3)then

        do i=1,d
        write(1,fmt=700)i,N(i)
	end do
 700    format (8x,'N(',i1,') =', i4)

        write(1,*)' '
	write(1,*)' 1/T  = ', beta, '  -U  = ', bun/beta/Nsite
	write(1,*)' \mu  = ', mu, ' nnn = ', g0000
	write(1,*)'  '

        write(1,*)' MC step (mln) =', step/1.d6,'  Z(mln) = ',Z/1.d6
	write(1,*)' Z_b(mln) = ',Z_b/1d6,' b_n = ', b_n
	write(1,*)' Z_b1(mln)= ',Z_b1/1d6,' b_n1= ', b_n1
        write(1,*)' '


	write(1,fmt=771)nmnm,nm_max,nm_min,nm_av
 771  format(' nmnm-> ',I5,' now [ max = ',I5,' min = ',I5,  ' average = ',G11.5,' ]')


!--- pot. energy -----------------------------------
!	call bSTAT(PE_stat(1:b_n),b_n,Z_b,PE_av,PE_err)
!	write(1,*)'  '
!      write(1,fmt=701) PE_av, PE_err
! 701  format(8x,'PE =',g12.5,4x,' +/- ',g10.3)

        write(1,*)'doing PE: '
        call mrg( PE_stat(1:b_n), b_n, Z_b )


      
!--- kin. energy -----------------------------------
        write(1,*)'doing KE: '
        call mrg( KE_stat(1:b_n), b_n, Z_b )

      

!--- number density ---------------------------------
        write(1,*)'doing dens: '
        call mrg( ndens_stat(1:b_n), b_n, Z_b )


!--- gim
       write(1,*)'doing diag_gim:'
       call mrg( gim_stat(1:b_n1), b_n1, Z_b1 )

! Fermi energy & momentum for the non-interacting gas of the same density:
!	kf = (3.d0*pi*pi*dens_av)**(1.d0/3.d0)
!	ef = kf**2 !/2.d0                          ! m=1
!
!	write(1,fmt=777)ef, kf
! 777  format(8x,'E_F =',g12.5,4x,' k_F = ',g12.5 )

!--- integrated correlator ---------------------------
!w	call bSTAT(im_stat(1:b_n),b_n,Z_b,im_av,im_err)
!w	write(1,*)' '
!w      write(1,fmt=703) im_av, im_err
!w 703  format(8x,'g_im(w=0,k=0) =',g12.5,4x,' +/- ',g12.5)


!-----------------------------------------------------
	write(1,*)
	write(1,*)' address/accept:'
	write(1,*)' add/drop   ',a_a_v/(c_a_v+.01),' / ',a_d_v/(c_d_v+.01)
	write(1,*)' leap_o a/d ',a_lo_a/(c_lo_a+.01),' / ',a_lo_d/(c_lo_d+.01)
!w	write(1,*)' cre/ann    ',a_cre/(c_cre+.01),' / ',a_ann/(c_ann+.01)
!w	write(1,*)' leap a/d   ',a_l_a/(c_l_a+.01),' / ',a_l_d/(c_l_d+.01)
!w	write(1,*)' hop        ',a_r/(c_r+.001)
	write(1,*)' recalc ', recalc_rate/(step + 0.0001)   !w,' / '  ,  i_backsc/(a_l_a + .0001)

!---- nullify counters
      nm_max=0; nm_min=i_hu; nm_av=nul

!=============================  writeouts: various service distributions

!w! write t(ira-masha) distribution
!w	if(im_t(0)>0)then;	yyy = PE_av/U/Nsite/im_t(0)
!w	else; yyy=1.d0
!w	endif
!w
!w	open(2,file='ct3d'//trim(fnamesuffix)//'.dat')
!w	do j=-im_t_max,im_t_max-1
!w	  write(2,*)beta*(j+0.5d0)/im_t_max, im_t(j)*yyy
!w	enddo
!w  	write(2,*)beta, im_t(im_t_max)*yyy  ! the last point
!w	close(2)


! write kink number distrr
      xxx=sum(nmnm_distr(1:nmnm_sup))
      if(xxx==0.d0)xxx=1.d0

      open(2,file='nmnm'//trim(fnamesuffix)//'.dat')
	do j=0,nmnm_sup; write(2,*)j,nmnm_distr(j)/xxx
      enddo
	close(2)

!! write det distr
!	xxx=sum(det_distr(-det_distr_max:det_distr_max))
!	if(xxx==0.d0)xxx=1.d0
!
!	open(2,file='det'//trim(fnamesuffix)//'.dat')
!	do j=-det_distr_max,det_distr_max
!	   write(2,*)j, det_distr(j)/xxx
!	enddo
!	close(2)



!---------  uncomment this if you want to see the configuration -----------

! write kink visualization
!	open(2,file=trim(workdirpath)//'viz'//trim(fnamesuffix)//'.dat')
!	do i=1,nmnm; name=namelist(i)
!		write(2,*)x(:,ksite(name)),ktau(name)
!	enddo
!	close(2)
!
!	if(present)then
!
!	open(2,file=trim(workdirpath)//'viz_i'//trim(fnamesuffix)//'.dat')
!		write(2,*)x(:,ksite(ira)), ktau(ira)
!	close(2)
!
!	open(2,file=trim(workdirpath)//'viz_m'//trim(fnamesuffix)//'.dat')
!	    write(2,*)x(:,ksite(masha)),ktau(masha)
!	close(2)
!
!	endif
!---------------------------------------------------------------------------

	endif   ! b_n>3

	close(1)


! time to wrap up? --- write everything, release allocated memory and wrap up.
	if(lastone)then; 
	    open(1,file=outputfname,position='append')
		write(1,*)'Time to wrap up, dude..'; 
		close(1)
	    call wrt
		call mystop
	endif


	end subroutine prnt
	


!-------------------
!--- Write data
!-------------------
	subroutine wrt
	integer :: j, name

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'write .....'

! configuration
	open(1,file='cnf'//trim(fnamesuffix)//'.dat')
	    write(1,*)beta
!w            write(1,*)present
!w            if(present)then
!w                write(1,*)ksite(ira),ktau(ira),row(ira),clmn(ira)
!w                write(1,*)ksite(masha), ktau(masha),row(masha),clmn(masha)
!w            endif
	  write(1,*)nmnm,pm
	  write(1,*)nm_row(1:pm)
	  write(1,*)nm_clmn(1:pm)
	  do j=1,nmnm
	     name=namelist(j)
	     write(1,*)name, ksite(name),ktau(name), row(name),clmn(name)
	  enddo
	  write(1,*)namelist(nmnm+1:nmnm_max)
	close(1)

! rndm state
	rndmstr = 'rndm'//trim(fnamesuffix)//'.dat'
	call save_rndm(rndmstr)
!
! There is a switch in the parfile as to read 
! the generator state or to seed it anew.
!

! statistics
	open(1,file='stat'//trim(fnamesuffix)//'.dat')
	  write(1,*)d,N
	  write(1,*)beta, U, mu
	  write(1,*)step, Z
	  write(1,*)Z_b, b_n, i_b
	  write(1,*)PE
	  write(1,*)PE_stat(1:b_n)
	  write(1,*)KE
	  write(1,*)KE_stat(1:b_n)
	  write(1,*)ndens
	  write(1,*)ndens_stat(1:b_n)
!w	  write(1,*)im
!w	  write(1,*)im_stat(1:b_n)
          write(1,*)Z_b1, b_n1, i_b1
          write(1,*)gim
          write(1,*)gim_stat(1:b_n1) 
	  write(1,*)'------'
	  write(1,*)g_uu(0:N2(1))
	  write(1,*)g_ud(0:N2(1))
	close(1)

	write(OUT_UNIT,*)'writing done!'
	close(OUT_UNIT)

	end subroutine wrt


!---------------------
!--- Read statistics
!---------------------
	subroutine rd_stat
	integer :: ddd
	real*8  :: dummy
	character*6 :: cvoid


	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'reading stat .....'

	open(1,file='stat'//trim(fnamesuffix)//'.dat')
	  read(1,*)ddd,N
	      if(ddd/=d) call mystop
	  read(1,*)dummy, dummy, dummy
	  read(1,*)step, Z
	  read(1,*)Z_b, b_n, i_b
	  read(1,*)PE
	  read(1,*)PE_stat(1:b_n)
	  read(1,*)KE
	  read(1,*)KE_stat(1:b_n)
	  read(1,*)ndens
	  read(1,*)ndens_stat(1:b_n)
!w	  read(1,*)im
!w	  read(1,*)im_stat(1:b_n)
          read(1,*)Z_b1, b_n1, i_b1
          read(1,*)gim
          read(1,*)gim_stat(1:b_n1)
	  read(1,*)cvoid
	  read(1,*)g_uu(0:N2(1))
	  read(1,*)g_ud(0:N2(1))
	close(1)

	det_distr=0.d0; nmnm_distr=0.d0

	write(OUT_UNIT,*)'... done'
	close(OUT_UNIT)

	end subroutine rd_stat


!---------------------
!--- Init statistics
!---------------------
	subroutine init_stat

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'init stat .....'

	Z_b = ceiling(beta*U*Nsite) + 1    ! The block size for statistics 
	Z_b = Z_b / 8                      ! [something to start with]
	i_b = 0; b_n = 0

	Z=0.d0
!w	im = 0.d0; im_stat=0.d0
	PE = 0.d0; PE_stat=0.d0
	KE = 0.d0; KE_stat=0.d0
	ndens = 0.d0; ndens_stat=0.d0
	g_uu=0.d0; g_ud=0.d0

        Z_b1=Z_b
        i_b1=0; b_n1=0
        gim=0.d0; gim_stat=0.d0        
!---------------------------------

!w	im_t=0.d0; 
        det_distr=0.d0; nmnm_distr=0.d0; 

	write(OUT_UNIT,*)'... done'
	close(OUT_UNIT)


	end subroutine init_stat


!--------------------------
!--- Read configuration
!--------------------------
	subroutine rd_cnf
	integer :: j,name, site
	real*8  :: bbb, f
	character*50 :: rndmstr

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'reading conf.....'
	
! configuration
!w	namelist(ira)=ira; namelist(masha)=masha
        do j=1,nmnm_max; namelist(j)=j; numnam(j)=j
        enddo; nmnm=0; pm=0

	nkink=0; kname=0

	open(1,file='cnf'//trim(fnamesuffix)//'.dat')
	    read(1,*)bbb;        f = beta / bbb    
					! If the configuration if for some \beta_1, not equal to \beta
					! then all the times are just 'stretched' by a factor of f.
					! This feature is useful for re-thermalization of configurations
					! with close \beta-s
!w          read(1,*)present
!w         if(present)then
!w              read(1,*)ksite(ira),ktau(ira),row(ira),clmn(ira)
!w            read(1,*)ksite(masha), ktau(masha),row(masha),clmn(masha)
!w           endif
	  read(1,*)nmnm,pm
	     if(nmnm>nmnm_max)then
	        print*,'rd_cnf: nmnm= ',nmnm,' while nmnm_max= ',nmnm_max
			call mystop 
	     endif
!================== allocate enough memory to fit TheMatrix
	    deallocate(matr,m_u,m_v,m_w,m_z,nm_row,nm_clmn)    ! deallocate first
	    deallocate(m_u2,m_v2,m_s2,m_c2)
		lda=(int(pm/64)+1)*64+128
		allocate(nm_row(lda),nm_clmn(lda))
		allocate(matr(lda,lda))
		allocate(m_u(lda),m_v(lda),m_w(lda),m_z(lda))
	    allocate(m_u2(lda,2),m_v2(2,lda),m_s2(2,2),m_c2(lda,2) )
!=========================================
	  read(1,*)nm_row(1:pm)
	  read(1,*)nm_clmn(1:pm)
	  do j=1,nmnm
	      read(1,*)name, site,ktau(name), row(name), clmn(name)
	      namelist(j)=name; numnam(name)=j
	      ktau(name) = ktau(name) * f          ! rescale \tau
	         ksite(name)=site
	         nkink(site)=nkink(site)+1
	         kname(nkink(site),site)=name
	  enddo
	  read(1,*)namelist(nmnm+1:nmnm_max)
	close(1)

! calculate TheMatrix
	bbb = recalc_matrix(pm,lda,matr)
	

	write(OUT_UNIT,*)'... done'
	close(OUT_UNIT)

	end subroutine rd_cnf


!--------------------------
!--- Init configuration
!--------------------------
	subroutine init_cnf
	integer :: j

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'init conf.....'

! configuration
	nkink=0; kname=0
	ksite=0; ktau=0.d0
	row=0; clmn=0; pm=0
	nm_row=0; nm_clmn=0

! name manager
!w	namelist(ira)=ira; namelist(masha)=masha	
      do j=1,nmnm_max
          namelist(j)=j; numnam(j)=j
      enddo; nmnm=0; 

! no ira, masha @ the beginning
!w      present=.false.

	write(OUT_UNIT,*)'... done'
	close(OUT_UNIT)

	end subroutine init_cnf



!----------------------------------------------
!---  Green Function (former selector + spline interp.)
!----------------------------------------------
      real*8 function GREENFUN(site1,tau1,site2,tau2)
      implicit none
      integer :: site1,site2,j, sgn
      double precision :: tau, tau1, tau2, dt, gre

      integer :: nx, ny, nz, nta  !, ntb
      double precision :: tta,ttb,ga,gb,c, gra,grb   !,p


! prepare \tau
      tau=tau1-tau2
      dt=tau; sgn=1

	if(tau < 1.d-14)then; dt=beta+tau; sgn=-1; endif
! Explanation: G(t=0) must be understood as G(t-> -0) = -G(t=\beta)
! A long way to accomplish this is below, commented out. A short way is above :).
!
!     if (abs(tau) < 1.d-14) then; dt = beta; sgn=-1
!     else if (tau > 0.d0) then; dt=tau; sgn=1
!	else; dt=beta+tau; sgn=-1
!	end if
!----------------------------------------



! prepare coords
!	j=1; nx = abs(x1(j)-x2(j)); nx = min(nx,N(j)-nx)
!	j=2; ny = abs(x1(j)-x2(j)); ny = min(ny,N(j)-ny)
!	j=3; nz = abs(x1(j)-x2(j)); nz = min(nz,N(j)-nz)

!----------------------------------- spline

	nta=dt*bmt1 !*p

      tta=dt-nta*bmt 
	ttb=tta - bmt     !dt-ntb*(beta/mtau) 
!cccccccccccccccccccccccccccccccccccccc
      
	ga=GR_DAT(nta,site1,site2)
	gb=GR_DAT(nta+1,site1,site2)

	gra=GRD_DAT(nta,site1,site2)
	grb=GRD_DAT(nta+1,site1,site2)

      c=(ga-gb)*bmt1

      gre=(c+gra)*ttb + (c+grb)*tta
      gre=gre*tta*ttb*bmt1 + gb*tta-ga*ttb
      gre=gre*bmt1


	GREENFUN = gre*sgn

      end function GREENFUN


!-------------------------------------
!     Tabulates Green function and its time derivate at positive tau.
!-------------------------------------
      subroutine TABULATE
      implicit none
	real*8, allocatable :: ham(:,:)
	integer :: site,site1,j
	real*8 :: factor, ww,ttt,term, gamma, expet(0:mtau)
	integer :: nt

	! lapack stuff
	character*1 :: jobz,uplo
	integer     :: ldh, lwork,info
	real*8, allocatable  :: work(:), eps(:)

!	integer :: x1(d),x2(d) , site0

	print*,'start with TABULATE'

! check the array sizes
!	if(Nsite>ntab_max)then
!	    print*,' ******************************'
!	    print*,' TABULATE: need to have ntab_max >= Nsite '
!	    print*,' currently Nsite = ',Nsite, 'and ntab_max = ', ntab_max
!	    call mystop
!	endif



!--------- site-dependent chemical potential 
	do site=1,Nsite
	   musite(site)=mu + (2.d0*rndm() - 1.d0)*disord
	enddo

!!!        do site=1,Nsite
!!!                print*, site, musite(site)
!!!        enddo
!!!        pause


! build the hamiltonian
	allocate(ham(1:Nsite,1:Nsite)) ; 
	ham=0.d0
	do site=1,Nsite
	   ham(site,site)=ham(site,site)-musite(site)
	   do j=1,dd; site1=ass(j,site); 
	      ham(site,site1)=ham(site,site1)-1.d0  
	   enddo
	enddo;	

!  compute eigenvalues; for LAPACK parameters and arguments, see
!  http://www.netlib.org/lapack/double/dsyev.f
!  SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

	jobz='V'  ! compute eigenvectorz
	uplo='U'  ! use upper diag of ham(:,:) --- doesn't matter really
	ldh=Nsite
	lwork=12*Nsite
	allocate( work(lwork), eps(Nsite) )

! query the optimal workspace size
	call dsyev(jobz,uplo,Nsite,ham,ldh,eps,work,-1,info)
	lwork=work(1)
	deallocate(work); allocate(work(lwork))

! diagonalize
	call dsyev(jobz,uplo,Nsite,ham,ldh,eps,work,lwork,info)

	if(info/=0)then; print*,'*** dsyev returns info = ', info
			 print*,'*** check the TABULATE routine'
			 call mystop
	endif


!	print*,'optimal lwork =', work(1),lwork
!
!	print*,'eigenvalues:'
!	do j=1,Nsite
!		print*,eps(j)
!	enddo
!
!	print*; print*,'----------- eigenvector # '
!	do j=1,Nsite;
!		!write(1,*)j,ham(j,6)
!		print*,j,ham(j,1)
!	enddo
!	print*; print*; print*

!------------- have the spectrum, proceed to the GFs
	GR_DAT=0.d0; GRD_DAT=0.d0	

	do j=1,Nsite

	  gamma=-eps(j)*beta
	  gamma=exp(gamma)+1.d0
	  gamma=-1.d0/gamma

	  ww = exp(-eps(j)*bmt) ! bmt=beta/mtau
	  do nt=0,mtau; expet(nt)=ww**nt
	  enddo

          do site=1,Nsite ; do site1=1,Nsite
	    factor = ham(site,j)*ham(site1,j)
            do nt=0,mtau
               term = factor*expet(nt)*gamma    !/Nsite
               GR_DAT(nt,site,site1) = GR_DAT(nt,site,site1) + term
               GRD_DAT(nt,site,site1) = GRD_DAT(nt,site,site1) -eps(j)*term
            enddo

          enddo; enddo ! site, site1

	enddo   ! j: eigenvalues
!------------------------


! fill in fictitious nt=mtau+1, see GREENFUN for explanation
	GR_DAT(mtau+1,:,:)=0.d0; GRD_DAT(mtau+1,:,:)=0.d0

! fill g0000
	site = 1; ttt=0.d0
	g0000 = GREENFUN(site,ttt,site,ttt)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!	site0=1; site1=site0
!	do j=1,N(1)
!	  print*,site0,site1,GR_DAT(100,site0,site1),GRD_DAT(100,site0,site1)
!	  site1=ass(1,site1); site1=ass(2,site1)
!	enddo

	print*,'done with TABULATE'

      end subroutine TABULATE





!----------------------------------------------
!---  Green Function, spline interpolation : UNIFORM conventions
!----------------------------------------------
!      real*8 function GREENFUN_uniform(x1,tau1,x2,tau2)
!      implicit none
!	integer :: x1(1:d), x2(1:d),j, sgn
!      double precision :: tau, tau1, tau2, dt, gre
!
!	integer :: nx, ny, nz, nta  !, ntb
!	double precision :: tta,ttb,ga,gb,c, gra,grb   !,p
!
!! prepare \tau
!      tau=tau1-tau2
!      dt=tau; sgn=1
!
!	if(tau < 1.d-14)then; dt=beta+tau; sgn=-1; endif
!! G(t=0) must be understood as G(t-> -0) = -G(t=\beta)
!! A long way to accomplish this is below, commented out. A short way is above :).
!!
!!     if (abs(tau) < 1.d-14) then; dt = beta; sgn=-1
!!     else if (tau > 0.d0) then; dt=tau; sgn=1
!!	else; dt=beta+tau; sgn=-1
!!	end if
!!----------------------------------------
!
!
!
!! prepare coords, don't forger about PBC
!	j=1; nx = abs(x1(j)-x2(j)); nx = min(nx,N(j)-nx)
!	j=2; ny = abs(x1(j)-x2(j)); ny = min(ny,N(j)-ny)
!	j=3; nz = abs(x1(j)-x2(j)); nz = min(nz,N(j)-nz)
!
!!----------------------------------- spline
!
!	nta=dt*bmt1    ! Recall, bmt=beta/mtau, bmt1=1/bmt 
!
!      tta=dt-nta*bmt 
!	ttb=tta - bmt   
!
!
!cccccccccccccccccccccccccccccccccccccc
!      
!	ga=GR_DAT(nta,nx,ny,nz)
!	gb=GR_DAT(nta+1,nx,ny,nz)
!
!	gra=GRD_DAT(nta,nx,ny,nz)
!	grb=GRD_DAT(nta+1,nx,ny,nz)
!
!      c=(ga-gb)*bmt1
!
!      gre=(c+gra)*ttb + (c+grb)*tta
!      gre=gre*tta*ttb*bmt1 + gb*tta-ga*ttb
!      gre=gre*bmt1
!
!
!	GREENFUN = gre*sgn
!
!      end function GREENFUN_uniform



!-------------------------------------
!     Tabulates Green function and its time derivate at positive tau: UNIFORM case (momentum's conserved)
!-------------------------------------
!      subroutine TABULATE_uniform
!      implicit none
!!
!! These tabulated values will be used in greenfun() for the spline interpolation. 
!!
!      integer :: mx, my, mz, i
!      integer :: nx, ny, nz,nt
!      double precision :: tau
!      double precision :: phase, eps, gamma !, xi
!      double precision :: pp(d), s1, s2, term
!      double precision :: co(-1000:1000,d)
!
!	integer :: xxx(1:d)
!	real*8 :: ttt
!
!	real*8 :: kx,ky,kz ,  expet(0:mtau), ww
!
!
!	do i=1,d; pp(i)=4.d0*asin(1.d0)/N(i); enddo
!
! 
!!------------------ Cosines  ------------------
!      do i=1,d
!	  do mx=-N(i)/2+1, N(i)/2;                ! The single-particle dispersion is tight-binding,
!	          co(mx,i)=-2.d0*cos(pp(i)*mx);   ! thus the spectrum is cosine.
!	end do
!	enddo
!!----------------------------------------------
!
!! nullify'em
!	GR_DAT=0.d0; GRD_DAT=0.d0	
!
!! sum over momentums 1st                 ! This weird loop sequence is UNAVOIDABLE
!	do mz=-N(3)/2+1, N(3)/2            ! in order to ensure the data locality:
!	do my=-N(2)/2+1, N(2)/2;           ! otherwise tabulation for L>10 takes
!     do mx=-N(1)/2+1, N(1)/2;           ! hours and hours.
!
!! spectrum
!        eps=co(mx,1)+co(my,2)+co(mz,3)-mu
!
!	  gamma=-eps*beta
!	  gamma=exp(gamma)+1.d0
!	  gamma=-1.d0/gamma
!
!	  ww = exp(-eps*bmt)             ! bmt=beta/mtau
!	  do nt=0,mtau; expet(nt)=ww**nt
!	  enddo
!
!! coordinates 
!          do nz = 0, Ntab-1; do ny = 0, Ntab-1;  do nx = 0, Ntab-1
!	       phase=pp(1)*nx*mx+pp(2)*ny*my+pp(3)*nz*mz
!	
!	       do nt=0,mtau
!
!		     term = cos(phase)*expet(nt)*gamma/Nsite
!
!		     GR_DAT(nt,nx,ny,nz) = GR_DAT(nt,nx,ny,nz) + term
!		     GRD_DAT(nt,nx,ny,nz) = GRD_DAT(nt,nx,ny,nz) -eps*term
!
!		enddo; enddo; enddo  ! coordinates
!	  enddo                  ! tau
!	enddo; enddo; enddo      ! momentum
!!------------------------
!
!
!! fill in fictitious nt=mtau+1, see GREENFUN for explanation
!	GR_DAT(mtau+1,:,:,:)=0.d0; GRD_DAT(mtau+1,:,:,:)=0.d0
!
!! fill g0000
!	xxx = 0; ttt=0.d0
!	g0000 = GREENFUN(xxx,ttt,xxx,ttt)
!
!! temp @#$ : ---------------------------------THIS IS FOR TESTING ONLY
!!	print*,"******** G(t=0,x,0,0) "
!!	do j=0,N(1)-1
!!	  print*,j,GR_DAT(0,j,0,0),GRD_DAT(0,j,0,0)
!!	enddo
!!	stop
!
!     end subroutine TABULATE_uniform


!-------------------------------
!--- Increase LDA (matrix storage): 
!      if lnew == -1 then just add 512
!-------------------------------
	subroutine resize_matrix(lnew)
	integer :: lnew

	integer :: lda_new 
	real*8, allocatable :: matr_new(:,:)
	real*8, allocatable :: nm_r_new(:), nm_c_new(:)

	if(lnew==-1)then; lda_new=lda+512
	else; lda_new = lnew
	endif

	allocate(matr_new(lda_new,lda_new))
	allocate(nm_r_new(lda_new),nm_c_new(lda_new))

! save TheMatrix as is
	matr_new(1:pm,1:pm)=matr(1:pm,1:pm)
	nm_r_new(1:pm)=nm_row(1:pm); nm_c_new(1:pm)=nm_clmn(1:pm)

! Resize
	deallocate(matr); deallocate(m_u,m_v,m_w,m_z,m_u2,m_v2,m_s2,m_c2)
	deallocate(nm_row,nm_clmn)

	lda=lda_new

	allocate(matr(lda,lda),m_u(lda),m_v(lda))
	allocate(m_w(lda),m_z(lda))

	allocate(m_u2(lda,2),m_v2(2,lda),m_s2(2,2),m_c2(lda,2) )

	allocate(nm_row(lda),nm_clmn(lda))


! Restore
	matr(1:pm,1:pm)=matr_new(1:pm,1:pm)
	nm_row(1:pm)=nm_r_new(1:pm); nm_clmn(1:pm)=nm_c_new(1:pm)

! Cleanup
	deallocate(matr_new,nm_r_new,nm_c_new)


	end subroutine resize_matrix


!-------------------------------
!--- Recalculate the inverse from scratch (N**3 operations)
!-------------------------------
	real*8 function recalc_matrix(pm,lda,a)
	integer :: pm,lda     ! actual size & leading dimension of A
        real*8  :: a(lda,lda)

	integer :: i,j, vova, lesha, xi(d), xj(d)
	real*8  :: ti,tj, det_big

        recalc_matrix=0.d0

	if(pm<2)return   ! no use to recalculate

	recalc_rate = recalc_rate + 1.d0

! build the matrix
	do j=1,pm; do i=1,pm
	  vova=nm_row(i); lesha=nm_clmn(j)
	  !xi=x(:,ksite(vova)); xj=x(:,ksite(lesha))
	  ti=ktau(vova); tj=ktau(lesha)
	  a(i,j)=GREENFUN(ksite(lesha),tj,ksite(vova),ti)
	enddo; enddo

! invert
	det_big = full_inv(pm,lda,a)
	
	recalc_matrix = det_big       ! return the absolute value of the determinant

	end function recalc_matrix


!--------------------------------------------
!--- Arranges associations between sites
!--------------------------------------------
      subroutine ASSA  
 
      integer :: site, site1, i, i1, i2 !, i3 
      integer :: ic(d), NN(d+1) 
!     ic(i) is the i-coordinate of a site, 
!     the relation between site indexes and coordinates reads:
!     site = 1 + (ic(1)-1) + N(1)*(ic(2)-1) +...
!               + N(1)*...*N(d-1)*(ic(d)-1)

      NN(1)=1; do i=2,d+1; NN(i)=NN(i-1)*N(i-1); enddo
       
      do i=1,d; back(i)=i+d; back(i+d)=i; enddo 
      
      ic=1; ic(1)=0
      DO site=1, Nsite

!------------ Coordinates for site ----------
         i1=1 
         DO
         if (ic(i1) < N(i1)) then
            ic(i1)=ic(i1)+1
            DO i2=1,i1-1; ic(i2)=1; END DO
            EXIT
         else; i1=i1+1;  end if 

         END DO

!-------------------------------------------------------


         DO i=1, d
            if (ic(i) < N(i)) then
               site1=site+NN(i)
            else
               site1=site+NN(i)-NN(i+1)
            end if

            ass(i,site)=site1
            ass(back(i),site1)=site
	      x(:,site)=ic(:)
                 
         END DO
      END DO
      
      end subroutine ASSA










!----------------------
!---- Getname function for the Name Manager
!----------------------
      SUBROUTINE GetName(nick)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: nick

      IF(nmnm<nmnm_max)THEN; 
	  nmnm=nmnm+1
        nick=namelist(nmnm); numnam(nick)=nmnm
      ELSE
	    ! has been moved into CheckGetName()
      ENDIF   

      END SUBROUTINE GetName

!------------------------------
!--- Check for GetName: nmnm>nmnm_max
!------------------------------
	subroutine CheckGetName

	if(nmnm>=nmnm_max)then; PRINT*,'GetName-> list is over!'
	   call mystop
	endif

	end subroutine CheckGetName

       
!----------------------
!---- DropName function for Name Manager
!----------------------
      SUBROUTINE DropName(nick)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nick
      INTEGER :: nm, lastname

      nm=numnam(nick)
      IF(nm<nmnm)THEN
         lastname=namelist(nmnm); namelist(nm)=lastname
         numnam(lastname)=nm; namelist(nmnm)=nick
         numnam(nick)=nmnm; nmnm=nmnm-1
      ELSE IF(nm==nmnm)THEN
		 nmnm=nmnm-1
      ELSE
	   ! has been moved into the CheckDropName()
      ENDIF
      RETURN 

      END SUBROUTINE DropName 


!------------------------------
!--- Check whether the name to be dropped exists
!------------------------------
	subroutine CheckDropName(nick)
	integer :: nick

	if(numnam(nick)>nmnm)then; 
	PRINT*,'DropName-> No such name:',nick; call mystop
	endif

	end subroutine CheckDropName



!---------------------------------
!--- Calculate the det from scratch & compare
!---------------------------------
	subroutine check_recalc

	real*8, allocatable :: a(:,:)
	real*8 :: det_1
	integer :: i,j

	if(pm<2)return
	
      allocate(a(lda,lda))
      det_1 = recalc_matrix(pm,lda,a)
	
      do i=1,pm; do j=1,pm
	    if(abs(a(i,j)-matr(i,j))>1d-8)then 
	       print*; print*,'check_recalc failed'
	       print*,i,j
		   print*,matr(i,j),a(i,j),abs(a(i,j)-matr(i,j))
	       print*,status,step,pm
	       print*,det_1
	       call mystop
	    endif
	enddo; enddo


	deallocate(a)


	end subroutine check_recalc


!---------------------------------
!--- Thermalize from scratch:
!      - diag updates only
!      - tune U
!---------------------------------
	subroutine therm1
	logical :: acpt
	real*8 :: r, pex, bun_in, bun_fin, det

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)' '
	write(OUT_UNIT,*)'Thermalization [diag] will be done by ', n_sw,' sweeps'
	close(OUT_UNIT)

	bun_in = beta*U_in*Nsite; bun_fin=bun

	step_t = n_sw*bun
	i_p=0.d0; i_w=0.d0; i_t=0.d0; step=0.d0

	do;

         step=step+un1                          
         i_p=i_p+un1; i_w=i_w+un1; i_t=i_t+1.d0
	   if(step>=step_t)exit

	   pex=step_t/i_t-1.d0

	   if (pex < 0.05) then; bun=bun_fin
	   else
            bun = bun_fin + (bun_in - bun_fin)*exp(-1.d0/pex)
	   end if

!--- Update: diagonal add/drop ONLY
	   r = rndm()
	   if     (r<0.5)then;    call add(acpt,det)
	   else if(r<=1.d0)then;  call drop(acpt,det)
	   endif

!------------- recalculate if necessary -------------
	   if( acpt ) then

	      if(abs(det)>tolerance) det=recalc_matrix(pm,lda,matr)

	   endif

	   if (i_r == step_r) then; i_r=0; 
	        det=recalc_matrix(pm,lda,matr); 
	   endif
!-------------------------------------------------

         if (i_p  == step_p)  then; i_p=nul; call PRNT;    end if  

	enddo

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'Thermalization done by ', step_t/1d6,	' mln steps'
	write(OUT_UNIT,*)'  '
	close(OUT_UNIT)

	if(n_sw>0) call wrt   

	end subroutine therm1


!---------------------------------
!--- Thermalize via worm scheme
!---------------------------------
	subroutine therm2
	logical :: acpt
	real*8 :: r, det


	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)' '
	write(OUT_UNIT,*)'Thermalization [worm] will be done by ',	 n_sw,' sweeps'
	close(OUT_UNIT)

	step_t = n_sw*bun
	i_p=0.d0; i_w=0.d0; i_t=0.d0; step=0.d0

	do;

         step=step+un1                          
         i_p=i_p+un1; i_w=i_w+un1; i_t=i_t+1.d0; i_r = i_r + 1.d0
	   if(step>=step_t)exit


!--- Update :)
	   r = rndm()
	   if     (r<prob(1))then;    call add(acpt,det)
	   else if(r<prob(2))then;    call drop(acpt,det)
	   else if(r<prob(3))then;    call leap_o_add(acpt,det)
	   else if(r<prob(4))then;    call leap_o_drop(acpt,det)
!w	   else if(r<prob(5))then;    call create(acpt,det)
!w	   else if(r<prob(6))then;    call annihilate(acpt,det)
!w	   else if(r<prob(7))then;    call leap_add(acpt,det)
!w	   else if(r<prob(8))then;    call leap_drop(acpt,det)
!w	   else if(r<=prob(9))then;   call hop(acpt,det)
	   else; print*,'am I nuts or what???'; call mystop
	   endif


!------------- recalculate if necessary -------------
	   if( acpt ) then

	      if(abs(det)>tolerance) det=recalc_matrix(pm,lda,matr)

	   endif

	   if (i_r == step_r) then; i_r=0; 
	        det=recalc_matrix(pm,lda,matr); 
	   endif
!-------------------------------------------------

         if (i_p  == step_p)  then; i_p=nul; call PRNT;    end if  

	enddo

	open(OUT_UNIT,file=outputfname,position='append')
	write(OUT_UNIT,*)'Thermalization done by ', step_t/1d6,	' mln steps'
	write(OUT_UNIT,*)'  '

	if(n_sw>0) call wrt   

	end subroutine therm2


!--------------------------
!--- clean up & stop 
!--------------------------
	subroutine mystop

	if(allocated(musite))deallocate(musite)
	if(allocated(ass))deallocate(ass)
	if(allocated(back))deallocate(back); 
	if(allocated(x))deallocate(x)
!w	if(allocated(uzli))deallocate(uzli) 
!w	if(allocated(yes))deallocate(yes)
!w	if(allocated(uzli_cube_le))deallocate(uzli_cube_le)
!w	if(allocated(uzli_cube_ca))deallocate(uzli_cube_ca)
!w	if(allocated(dx_le))deallocate(dx_le)
!w	if(allocated(dt_le))deallocate(dt_le)
!w	if(allocated(w_le))deallocate(w_le)
!w	if(allocated(name_best))deallocate(name_best)
!w	if(allocated(dist_best))deallocate(dist_best)
!w	if(allocated(distance))deallocate(distance)
!w	if(allocated(dx_le))deallocate(dx_le)
	if(allocated(nkink))deallocate(nkink)
	if(allocated(kname))deallocate(kname)
	if(allocated(ksite))deallocate(ksite)
	if(allocated(ktau))deallocate(ktau)
	if(allocated(row))deallocate(row)
	if(allocated(clmn))deallocate(clmn)
	if(allocated(nm_row))deallocate(nm_row)
	if(allocated(nm_clmn))deallocate(nm_clmn)
	if(allocated(matr))deallocate(matr)
	if(allocated(m_u))deallocate(m_u)
	if(allocated(m_v))deallocate(m_v)
	if(allocated(m_w))deallocate(m_w)
	if(allocated(m_z))deallocate(m_z)
	if(allocated(m_u2))deallocate(m_u2)
	if(allocated(m_v2))deallocate(m_v2)
	if(allocated(m_s2))deallocate(m_s2)
	if(allocated(m_c2))deallocate(m_c2)


!!!	call MPI_FINALIZE( ierr )
	stop

	end subroutine mystop


	end program MAIN
