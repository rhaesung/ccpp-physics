!>\file cires_ugwpv1_initialize.F90
!!

!===============================
! cu-cires ugwp-scheme
!  initialization of selected 
!  init gw-solvers (1,2,3,4)
!  init gw-source specifications
!  init gw-background dissipation
!==============================
!
! Part-0   specifications of common constants, limiters and "criiical" values
! 
!
!> This module contains common constants, limiters and "critical" values in      
    module ugwp_common
!
     use mo_ugwpv1_kind, only : wp

     implicit none
     
      real(kind=wp)   ::  pi, pi2, pih, rad_to_deg, deg_to_rad
      real(kind=wp)   ::  arad, p0s     
      real(kind=wp)   ::  grav, grav2, rgrav, rgrav2
      real(kind=wp)   ::  cpd,  rd, rv, fv            
      real(kind=wp)   ::  rdi,  rcpd, rcpd2  
      
      real(kind=wp)   ::  gor,  gr2,  grcp, gocp, rcpdl, grav2cpd
      real(kind=wp)   ::  bnv2min, bnv2max
      real(kind=wp)   ::  dw2min,  velmin, minvel        
      real(kind=wp)   ::  omega1, omega2,   omega3   
      real(kind=wp)   ::  hpscale, rhp, rhp2, rh4, rhp4, khp, hpskm
      real(kind=wp)   ::  mkzmin, mkz2min,  mkzmax, mkz2max, cdmin    
      real(kind=wp)   ::  rcpdt      
           
!      real(kind=wp), parameter ::  grav2 = grav + grav
!      real(kind=wp), parameter ::  rgrav = 1.0_wp/grav, rgrav2= rgrav*rgrav           
!      real(kind=wp), parameter ::  rdi  = 1.0_wp / rd, rcpd = 1._wp/cpd, rcpd2 = 0.5_wp/cpd  
!      real(kind=wp), parameter ::  gor  = grav/rd, rcpdt = 1._wp/(cp*dtp)
      
!      real(kind=wp), parameter ::  gr2  = grav*gor
!      real(kind=wp), parameter ::  grcp = grav*rcpd, gocp = grcp
!      real(kind=wp), parameter ::  rcpdl    = cpd*rgrav                          ! 1/[g/cp]  == cp/g
!      real(kind=wp), parameter ::  grav2cpd = grav*grcp                          !  g*(g/cp)= g^2/cp
!      real(kind=wp), parameter ::  pi2 = 2._wp*pi, pih = .5_wp*pi  
!      real(kind=wp), parameter ::  rad_to_deg=180.0_wp/pi, deg_to_rad=pi/180.0_wp
               
!      real(kind=wp), parameter ::  bnv2min = (pi2/1800._wp)*(pi2/1800._wp)
!      real(kind=wp), parameter ::  bnv2max = (pi2/30._wp)*(pi2/30._wp)      
!      real(kind=wp), parameter :: dw2min=1.0_wp,  velmin=sqrt(dw2min), minvel = 0.5_wp     
!      real(kind=wp), parameter :: omega1 = pi2/86400._wp, omega2 = 2._wp*omega1, omega3 = 3._wp*omega1     
   
!      real(kind=wp), parameter :: hpscale= 7000._wp, rhp=1._wp/hpscale, rhp2=.5_wp*rhp, rh4 = 0.25_wp*rhp 
!      real(kind=wp), parameter :: mkzmin = pi2/80.0e3_wp, mkz2min = mkzmin*mkzmin
!      real(kind=wp), parameter :: mkzmax = pi2/500._wp,   mkz2max = mkzmax*mkzmax 
!      real(kind=wp), parameter :: cdmin  = 2.e-2_wp/mkzmax 
!      real(kind=wp), parameter ::  pi = 4._wp*atan(1.0_wp),
!      real(kind=wp), parameter ::  grav =9.81_wp, cpd = 1004._wp
!      real(kind=wp), parameter ::  rd = 287.0_wp , rv =461.5_wp  
!      real(kind=wp), parameter ::  fv   = rv/rd - 1.0_wp    
!      real(kind=wp), parameter ::  arad = 6370.e3_wp
        
      contains

      subroutine init_nazdir(naz,  xaz,  yaz)
      
      implicit none
      
      integer :: naz
      real(kind=wp), dimension(naz) :: xaz,  yaz
      integer :: idir
      real(kind=wp)    :: phic, drad
      
      drad  = pi2/real(naz, kind=wp)
      if (naz.ne.4) then     
        do idir =1, naz
         Phic = drad*(real(idir, kind=wp)-1.0_wp)
         xaz(idir) = cos(Phic)
         yaz(idir) = sin(Phic)
        enddo
      else 
!                                  if (naz.eq.4) then
          xaz(1) = 1.0_wp     !E
          yaz(1) = 0.0_wp
          xaz(2) = 0.0_wp     
          yaz(2) = 1.0_wp     !N
          xaz(3) =-1.0_wp     !W
          yaz(3) = 0.0_wp
          xaz(4) = 0.0_wp
          yaz(4) =-1.0_wp     !S
      endif      
      end  subroutine init_nazdir     
!
!
!===================================================
!
!Part-1 init =>   wave dissipation + RFriction
!
!===================================================
     subroutine init_global_gwdis(levs, zkm, pmb, kvg, ktg, krad, kion, me, master)
                                    
     implicit none
     integer , intent(in)                 :: me, master
     integer , intent(in)                 :: levs
     real(kind=wp), intent(in)                      :: zkm(levs), pmb(levs)    ! in km-Pa
     real(kind=wp), intent(out), dimension(levs+1) :: kvg, ktg, krad, kion
!
!locals + data
!
     integer :: k
     real(kind=wp), parameter :: vusurf  = 2.0e-5_wp
     real(kind=wp), parameter :: musurf  = vusurf/1.95_wp
     real(kind=wp), parameter :: hpmol   = 7.0_wp
!
     real(kind=wp), parameter :: kzmin   = 0.1_wp
     real(kind=wp), parameter :: kturbo  = 100.0_wp
     real(kind=wp), parameter :: zturbo  = 130.0_wp
     real(kind=wp), parameter :: zturw   = 30.0_wp
     real(kind=wp), parameter :: inv_pra = 3.0_wp  !kt/kv =inv_pr
!
     real(kind=wp), parameter :: alpha   = 1.0_wp/86400.0_wp/15.0_wp   ! height variable see Zhu-1993 from 60-days => 6 days
     real(kind=wp)  ::  pa_alp  = 750.0_wp                             ! super-RF parameters from FV3-dycore GFSv17/16 sett
     real(kind=wp)  ::  tau_alp = 10.0_wp                              ! days   (750 Pa /10days)
!     
     real(kind=wp), parameter :: kdrag   = 1.0_wp/86400.0_wp/30.0_wp   !parametrization for WAM ion drag as e-density function
     real(kind=wp), parameter :: zdrag   = 100.0_wp
     real(kind=wp), parameter :: zgrow   = 50.0_wp
!
     real(kind=wp) ::   vumol, mumol, keddy, ion_drag
     real(kind=wp) ::   rf_fv3, rtau_fv3, ptop, pih_dlog       
!
     real(kind=wp) ::   ae1 ,ae2
!     
      
     ptop = pmb(levs)  
     rtau_fv3 = 1.0_wp/86400.0_wp/tau_alp
     pih_dlog = pih/log(pa_alp/ptop)

     do k=1, levs
      ae1 =    zkm(k)/hpmol
      vumol    = vusurf*exp(ae1)
      mumol    = musurf*exp(ae1)
      ae2 = -((zkm(k)-zturbo) /zturw)**2
      keddy    = kturbo*exp(ae2)

      kvg(k)   = vumol + keddy
      ktg(k)   = mumol + keddy*inv_pra

      krad(k)  = alpha
!
      ion_drag = kdrag
!
      kion(k)  = ion_drag!      
! add  Rayleigh_Super of FV3  for pmb < pa_alp
!
      if (pmb(k) .le. pa_alp) then
       rf_fv3=rtau_fv3*sin(pih_dlog*log(pa_alp/pmb(k)))**2
       krad(k) = krad(k) + rf_fv3 
       kion(k) = kion(k) + rf_fv3   

      endif

!      write(6,132) zkm(k), kvg(k), kvg(k)*(6.28_wp/5000._wp)**2,  kion(k) 
     enddo

      k= levs+1
      kion(k) = kion(k-1)
      krad(k) = krad(k-1)
      kvg(k)  =  kvg(k-1)
      ktg(k)  =  ktg(k-1)
      
!      if (me == master) then
!          write(6, * ) '  zkm(k), kvg(k), kvg(k)*(6.28_wp/5000._wp)**2,  kion(k) ... init_global_gwdis'     
!         do k=1, levs, 1
!          write(6,132) zkm(k), kvg(k), kvg(k)*(6.28_wp/5000._wp)**2,  kion(k), pmb(k) 
!        enddo      
!      endif
!
! 132  format( 2x, F8.3,' dis-scales:', 4(2x, E10.3))    
                                                          
     end subroutine init_global_gwdis

  end module ugwp_common
!
! ========================================================================
! Part 2 - sources
!      wave  sources
! ========================================================================
!
!    ugwp_oro_init
!
!=========================================================================
     module ugwp_oro_init
     use mo_ugwpv1_kind, only : wp
     use ugwp_common,    only : bnv2min, grav, grcp, fv, cpd, pi
     use ugwp_common,    only : mkzmin, mkz2min
     implicit none
!  
! constants and "crirtical" values to run oro-mtb_gw physics
!
!
      real(kind=wp),      parameter :: hncrit=9000.0_wp   ! max value in meters for elvmax
      real(kind=wp),      parameter :: hminmt=50.0_wp     ! min mtn height (*j*)
      real(kind=wp),      parameter :: sigfac=4.0_wp      ! mb3a expt test for elvmax factor
      real(kind=wp),      parameter :: hpmax=2500.0_wp
      real(kind=wp),      parameter :: hpmin=25.0_wp      
!
!
      real(kind=wp),      parameter :: minwnd=1.0_wp      ! min wind component (*j*)
      real(kind=wp),      parameter :: dpmin=5000.0_wp    ! minimum thickness of the reference layer in pa

      
      character(len=8)  ::  strver  = 'gfs_2018'
      character(len=8)  ::  strbase = 'gfs_2018'
      
      real(kind=wp), parameter    ::  rimin=-10.0_wp, ric=0.25_wp
      
      real(kind=wp), parameter    ::  frmax=10.0_wp, frc =1.0_wp, frmin =0.01_wp
      real(kind=wp), parameter    ::  ce=0.8_wp,     ceofrc=ce/frc, cg=0.5_wp
      real(kind=wp), parameter    ::  gmax=1.0_wp,   veleps=1.0_wp, factop=0.5_wp!
      real(kind=wp), parameter    ::  efmin=0.5_wp,  efmax=10.0_wp
      
      real(kind=wp), parameter    ::  rlolev=50000.0_wp
      integer,       parameter    ::  mdir = 8
      real(kind=wp), parameter    ::  fdir=real(mdir, kind=wp)/(8.0_wp*atan(1.0_wp))
      real(kind=wp), parameter    ::  zpgeo=2.0_wp*atan(1.0_wp)
      
      integer nwdir(mdir)
      data nwdir/6,7,5,8,2,3,1,4/
      save nwdir
      
      real(kind=wp), parameter    :: odmin = 0.1_wp, odmax = 10.0_wp
      real(kind=wp), parameter    :: fcrit_sm  = 0.7_wp, fcrit_sm2  = fcrit_sm * fcrit_sm
      real(kind=wp), parameter    :: fcrit_gfs = 0.7_wp, fcrit_v1 = 0.7_wp
      real(kind=wp), parameter    :: fcrit_mtb = 0.7_wp

      real(kind=wp),  parameter   :: zbr_pi  = zpgeo
      real(kind=wp),  parameter   :: zbr_ifs = zpgeo
 
!  
 
      real(kind=wp),    parameter ::   kxoro=6.28e-3_wp/200.0_wp    !
      real(kind=wp),    parameter ::   coro = 0.0_wp
      integer,          parameter ::   nridge=2
      real(kind=wp),    parameter ::   sigma_std=1.0_wp/100.0_wp, gamm_std=1.0_wp       
 
      real(kind=wp)    ::  cdmb                      ! scale factors for mtb
      real(kind=wp)    ::  cleff                     ! scale factors for orogw
      
      integer ::  nworo                      ! number of waves
      integer ::  nazoro                     ! number of azimuths
      integer ::  nstoro                     ! flag for stochastic launch above SG-peak

      
!------------------------------------------------------------------------------
!    small-scale orography parameters  for TOFD of Beljaars et al., 2004, QJRMS
!    SA-option can be controlled by  Integral limits of fluxes
!    in B2004: klow = 0.003 1/m ~ 2km and kinf ~ 6.28/10/(Z1)~< 1 km => meters
!    these limits can change strength of TOFD... choice of k0tr ~1/10 km (10km ~dx of C768)
!    kmax = kdis_pbl
!------------------------------------------------------------------------------
     real(kind=wp), parameter      :: kmax = 6.28_wp/(10.0_wp*25.0_wp)        ! max k-tofd
     real(kind=wp), parameter      :: k1tr = 6.28_wp/(2100.0_wp)              ! max k-transition from -1.9/slope to -2.8/slope
     real(kind=wp), parameter      :: kflt = 6.28_wp/(18.0e3_wp)              !
     real(kind=wp), parameter      :: k0tr = 6.28_wp/(10.0e3_wp)              ! min k-tofd     
     real(kind=wp), parameter      :: nk1tr = 2.8_wp
     real(kind=wp), parameter      :: nk0tr = 1.9_wp
     real(kind=wp), parameter      :: a1_tofd =  kflt ** nk1tr * 1.0e3_wp
     real(kind=wp), parameter      :: a2_tofd =  k1tr ** (nk0tr-nk1tr)    
     real(kind=wp), parameter      :: fix_tofd = 2.0_wp * 0.005_wp * 12.0_wp * 0.6_wp  !value= 0.072
!     
! B2004 scheme is based on the empirical vertical profile of the tofd divergence: 
!       Ax_tofd(Z)=exp(-[Z/ze_tofd]^3/2) / Z^1.2.....
! TOFD-flux/TMS-flux must dissipate due to PBL-diffusion with spectral damping
! Here we can enhance TOFD-impact by selecting k0tr and kmax limits
!     as functions of resolution and PBL-dissipation
!                            
     integer, parameter   :: n_tofd = 2                                          ! depth of SSO for TOFD compared with Zpbl
     real(kind=wp), parameter      :: const_tofd = 0.0759_wp                     ! alpha*beta*Cmd*Ccorr*2.109 = 12.*1.*0.005*0.6*2.109 = 0.0759
     real(kind=wp), parameter      :: ze_tofd    = 1500.0_wp                     ! BJ's z-decay in meters, 1.5 km
     real(kind=wp), parameter      :: a12_tofd   = 0.0002662_wp*0.005363_wp      ! BJ's k-spect const for sigf2 * a1*a2*exp(-[z/zdec]**1.5]
     real(kind=wp), parameter      :: ztop_tofd  = 3.0_wp*ze_tofd                ! no TOFD > this height 4.5 km
!------------------------------------------------------------------------------
!

      contains
!
      subroutine init_oro_gws(nwaves, nazdir, nstoch, effac, &
                              lonr, kxw )
!
!
      integer :: nwaves, nazdir, nstoch
      integer :: lonr                          
      real(kind=wp)    :: cdmbX
      real(kind=wp)    :: kxw
      real(kind=wp)    :: effac        ! it is analog of cdmbgwd(2) for GWs, off for now
!-----------------------------! GFS-setup for cdmb & cleff
! cdmb =  4.0 * (192.0/IMX)  
! cleff = 0.5E-5 / SQRT(IMX/192.0)  = 0.5E-5*SQRT(192./IMX)
!
      real(kind=wp), parameter :: lonr_refmb =  4.0_wp * 192.0_wp
      real(kind=wp), parameter :: lonr_refgw =  192.0_wp
      real(kind=wp), parameter :: cleff_ref  =  0.5e-5_wp           !       1256 km = 10 * 125 km ???
      
! copy  to "ugwp_oro_init"  =>  nwaves, nazdir, nstoch
 
      nworo  =  nwaves
      nazoro =  nazdir
      nstoro =  nstoch

      cdmbX = lonr_refmb/real(lonr, kind=wp)
      
      cdmb  = cdmbX
      cleff = cleff_ref * sqrt(lonr_refgw/real(lonr, kind=wp))   !* effac
!
      end subroutine init_oro_gws
!

    end module ugwp_oro_init
! =========================================================================
!
!    ugwp_conv_init
!
!=========================================================================
    module ugwp_conv_init
    
     use mo_ugwpv1_kind, only : wp

     
     implicit none
      real(kind=wp)    ::  eff_con                    ! scale factors for conv GWs
      integer ::  nwcon                      ! number of waves
      integer ::  nazcon                     ! number of azimuths
      integer ::  nstcon                     ! flag for stochastic choice of launch level above Conv-cloud
      real(kind=wp)    ::  con_dlength
      real(kind=wp)    ::  con_cldf

      real(kind=wp), parameter    :: cmin  =  5.0_wp  !2.5
      real(kind=wp), parameter    :: cmax  = 95.0_wp  !82.5
      real(kind=wp), parameter    :: cmid  = 22.5_wp
      real(kind=wp), parameter    :: cwid  = cmid
      real(kind=wp), parameter    :: bns   = 2.0e-2_wp, bns2 = bns*bns, bns4=bns2*bns2
      real(kind=wp), parameter    :: mstar = 6.28e-3_wp/2.0_wp  !  2km
      real(kind=wp)               :: dc

      real(kind=wp), allocatable  :: ch_conv(:),  spf_conv(:)
      real(kind=wp), allocatable  :: xaz_conv(:), yaz_conv(:)
     contains
!
     subroutine init_conv_gws(nwaves, nazdir, nstoch, effac, lonr, kxw)
!
     use ugwp_common,  only : pi2, arad, init_nazdir

     implicit none
     
 
      integer :: nwaves, nazdir, nstoch
      integer :: lonr
!      
! ccpp
!      
      
      real(kind=wp)    :: kxw,  effac
      real(kind=wp)    :: work1 = 0.5_wp
      real(kind=wp)    :: chk, tn4, snorm
      integer :: k

      nwcon    = nwaves
      nazcon   = nazdir
      nstcon   = nstoch
      eff_con  = effac

      con_dlength = pi2*arad/real(lonr, kind=wp)      
!
! allocate & define spectra in "selected direction": "dc" "ch(nwaves)"
!
       if (.not. allocated(ch_conv))  allocate (ch_conv(nwaves))
       if (.not. allocated(spf_conv)) allocate (spf_conv(nwaves))
       if (.not. allocated(xaz_conv)) allocate (xaz_conv(nazdir))
       if (.not. allocated(yaz_conv)) allocate (yaz_conv(nazdir))

      dc = (cmax-cmin)/real(nwaves-1, kind=wp)
!
! we may use different spectral "shapes"
! for example FVS-93 "Desabeius"
!  E(s=1, t=3,m, w, k) ~ m^s/(m*^4 + m^4) ~ m^-3 saturated tail
!
      do k = 1,nwaves
         chk         = cmin + real(k-1, kind=wp)*dc
         tn4         = (mstar*chk)**4
         ch_conv(k)  =  chk
         spf_conv(k) =  bns4*chk/(bns4+tn4)
      enddo

      snorm = sum(spf_conv)
      spf_conv = spf_conv/snorm*1.5_wp
 
      call init_nazdir(nazdir,  xaz_conv,  yaz_conv)
     end subroutine init_conv_gws


    end module ugwp_conv_init
!=========================================================================
!
!    ugwp_fjet_init
!
!=========================================================================

   module ugwp_fjet_init
   use mo_ugwpv1_kind, only : wp   

   
   
      implicit none
      real(kind=wp)    ::  eff_fj                      ! scale factors for conv GWs
      integer ::  nwfj                        ! number of waves
      integer ::  nazfj                       ! number of azimuths
      integer ::  nstfj                       ! flag for stochastic choice of launch level above Conv-cloud
!
      real(kind=wp), parameter    ::  fjet_trig=0.0_wp    ! if ( abs(frgf) > fjet_trig ) launch GW-packet

      
      real(kind=wp), parameter    :: cmin =  2.5_wp
      real(kind=wp), parameter    :: cmax = 67.5_wp
      real(kind=wp)               :: dc
      real(kind=wp), allocatable  :: ch_fjet(:) , spf_fjet(:)
      real(kind=wp), allocatable  :: xaz_fjet(:), yaz_fjet(:)
     contains
     
     subroutine init_fjet_gws(nwaves, nazdir, nstoch, effac,lonr, kxw) 

     use ugwp_common,  only : pi2, arad, init_nazdir

     implicit none

      integer :: nwaves, nazdir, nstoch
      integer :: lonr
      real(kind=wp)    :: kxw,  effac , chk
 
      integer :: k

      nwfj   =  nwaves
      nazfj  =  nazdir
      nstfj  =  nstoch
      eff_fj =  effac

       if (.not. allocated(ch_fjet))  allocate (ch_fjet(nwaves))
       if (.not. allocated(spf_fjet)) allocate (spf_fjet(nwaves))
       if (.not. allocated(xaz_fjet)) allocate (xaz_fjet(nazdir))
       if (.not. allocated(yaz_fjet)) allocate (yaz_fjet(nazdir))
 
      dc = (cmax-cmin)/real(nwaves-1, kind=wp)
      do k = 1,nwaves
         chk         = cmin + real(k-1, kind=wp)*dc
         ch_fjet(k)  =  chk
         spf_fjet(k) =  1.0_wp
      enddo
      call init_nazdir(nazdir,  xaz_fjet,  yaz_fjet)

     end subroutine init_fjet_gws

    end module ugwp_fjet_init
!
!=========================================================================
!
!
     module ugwp_okw_init
!=========================================================================
       use mo_ugwpv1_kind, only : wp      
     
       implicit none

      real(kind=wp)    ::  eff_okw                      ! scale factors for conv GWs
      integer ::  nwokw                        ! number of waves
      integer ::  nazokw                       ! number of azimuths
      integer ::  nstokw                       ! flag for stochastic choice of launch level above Conv-cloud
!
      real(kind=wp), parameter    ::  okw_trig=0.0_wp      ! if ( abs(okwp) > okw_trig ) launch GW-packet

      real(kind=wp), parameter    :: cmin =  2.5_wp
      real(kind=wp), parameter    :: cmax = 67.5_wp
      real(kind=wp)               :: dc
      real(kind=wp), allocatable  :: ch_okwp(:),   spf_okwp(:)
      real(kind=wp), allocatable  :: xaz_okwp(:),  yaz_okwp(:)

     contains
!

                             
     subroutine init_okw_gws(nwaves, nazdir, nstoch, effac, lonr, kxw)
     use ugwp_common,  only : pi2, arad, init_nazdir

     implicit none

      integer :: nwaves, nazdir, nstoch
      integer :: lonr
      real(kind=wp)    :: kxw,  effac , chk
 
      integer :: k

      nwokw   =  nwaves
      nazokw  =  nazdir
      nstokw  =  nstoch
      eff_okw =  effac

       if (.not. allocated(ch_okwp))  allocate (ch_okwp(nwaves))
       if (.not. allocated(spf_okwp)) allocate (spf_okwp(nwaves))
       if (.not. allocated(xaz_okwp)) allocate (xaz_okwp(nazdir))
       if (.not. allocated(yaz_okwp)) allocate (yaz_okwp(nazdir))
      dc = (cmax-cmin)/real(nwaves-1, kind=wp)
      do k = 1,nwaves
         chk =  cmin + real(k-1, kind=wp)*dc
         ch_okwp(k) = chk
         spf_okwp(k) = 1.0_wp
      enddo
  
      call init_nazdir(nazdir,  xaz_okwp,  yaz_okwp)
!
     end subroutine init_okw_gws
 
     end module ugwp_okw_init

!=============================== end of GW  sources
!
!  init specific  gw-solvers (1,2,3,4)
!
!===============================
!  Part -3  init  wave solvers
!===============================

  module ugwp_lsatdis_init
     use mo_ugwpv1_kind, only : wp    
     implicit none

      integer  :: nwav, nazd
      integer  :: nst
      real(kind=wp)     :: eff
      integer, parameter  :: incdim = 4, iazdim = 4
!
     contains

     subroutine initsolv_lsatdis(me, master,  nwaves, nazdir, nstoch, effac, do_physb, kxw)

     implicit none
!
     integer  :: me, master
     integer  :: nwaves, nazdir
     integer  :: nstoch
     real(kind=wp)     :: effac
     logical  :: do_physb
     real(kind=wp)     :: kxw
!
!locals: define azimuths and Ch(nwaves) - domain when physics-based soureces
!                                         are not actibve
!
     integer :: inc, jk, jl, iazi, i, j, k

     if( nwaves == 0 .or. nstoch == 1 ) then
!                                redefine from the default
       nwav = incdim
       nazd = iazdim
       nst  = 0
       eff  = 1.0_wp
     else
!                                from input_nml multi-wave spectra
       nwav = nwaves
       nazd = nazdir
       nst  = nstoch
       eff  = effac
     endif
!
     end subroutine initsolv_lsatdis
!
  end module ugwp_lsatdis_init
!
!
  module ugwp_wmsdis_init
  
     use mo_ugwpv1_kind, only : wp    
     use ugwp_common,    only : arad, pi, pi2, hpscale, rhp, rhp2, rh4, omega2 
     use ugwp_common,    only : bnv2max,  bnv2min, minvel 
     use ugwp_common,    only : mkzmin,  mkz2min, mkzmax, mkz2max, ucrit => cdmin
     
    implicit none

      real(kind=wp),     parameter   :: maxdudt = 250.0e-5_wp, maxdtdt=15.0e-2_wp
      real(kind=wp),     parameter   :: dked_min =0.01_wp,     dked_max=250.0_wp
   
      real(kind=wp),     parameter   :: gptwo=2.0_wp
 
      real(kind=wp) ,     parameter  :: bnfix  =  6.28_wp/300.0_wp, bnfix2= bnfix * bnfix
      real(kind=wp) ,     parameter  :: bnfix4 =  bnfix2 * bnfix2  
      real(kind=wp) ,     parameter  :: bnfix3 =  bnfix2 * bnfix        
!
! make parameter list that will be passed to SOLVER
! 
      integer  , parameter  :: iazidim=4        ! number of azimuths
      integer  , parameter  :: incdim=25        ! number of discrete cx - spectral elements in launch spectrum
 
      real(kind=wp) ,     parameter  :: zcimin = 2.5_wp
      real(kind=wp) ,     parameter  :: zcimax = 125.0_wp
      real(kind=wp) ,     parameter  :: zgam   = 0.25_wp
!
! Verical spectra
!
      real(kind=wp) ,     parameter  :: pind_wd = 5.0_wp/3.0_wp
      real(kind=wp) ,     parameter  :: sind_kz = 1.0_wp
      real(kind=wp) ,     parameter  :: tind_kz = 3.0_wp    
      real(kind=wp) ,     parameter  :: stind_kz = sind_kz + tind_kz  
!
! copies from kmob_ugwp namelist
!      
      real(kind=wp)                  :: nslope              ! the GW sprctral slope at small-m             
      real(kind=wp)                  :: lzstar
      real(kind=wp)                  :: lzmin
      real(kind=wp)                  :: lzmax      
      real(kind=wp)                  :: lhmet
      real(kind=wp)                  :: tamp_mpa            !amplitude for GEOS-5/MERRA-2
      real(kind=wp)                  :: tau_min             ! min of GW MF 0.25 mPa                         
      integer                        :: ilaunch
      real(kind=wp)                  :: gw_eff

      real(kind=wp)                  :: v_kxw, rv_kxw, v_kxw2 
     
     
 
!===========================================================================
      integer  :: nwav, nazd, nst
      real(kind=wp)     :: eff
 
      real(kind=wp)                :: zaz_fct, zms
      real(kind=wp), allocatable   :: zci(:), zci4(:), zci3(:),zci2(:), zdci(:)
      real(kind=wp), allocatable   :: zcosang(:), zsinang(:)
      real(kind=wp), allocatable   :: lzmet(:), czmet(:), mkzmet(:), dczmet(:), dmkz(:) 

!
! GW-eddy constants for wave-mode dissipation by background and stability of
!         "final" flow after application of GW-effects
! 
      real(kind=wp), parameter :: iPr_pt = 0.5_wp
      real(kind=wp), parameter :: lturb = 30.0_wp, sc2  = lturb*lturb      ! stable on 80-km TL lmix ~ 500 met.
      real(kind=wp), parameter :: ulturb=150.0_wp, sc2u = ulturb* ulturb   ! unstable
      real(kind=wp), parameter :: ric =0.25_wp
      real(kind=wp), parameter :: rimin = -10.0_wp, prmin = 0.25_wp
      real(kind=wp), parameter :: prmax = 4.0_wp
!
     contains
!============================================================================
     subroutine initsolv_wmsdis(me, master,  nwaves, nazdir, nstoch, effac, do_physb, kxw, version)
 
!        call initsolv_wmsdis(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
!          knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw,version)
 
     implicit none
!
!input-control for solvers:
!      nwaves, nazdir, nstoch, effac, do_physb, kxw
!
!
     integer, intent(in)  :: me, master, nwaves, nazdir, nstoch
     integer, intent(in)  :: version 
        
     real(kind=wp), intent(in)     :: effac, kxw
     logical,  intent(in) :: do_physb
      
!
!locals
!
      real(kind=wp)     :: dlzmet   
      real(kind=wp)     :: cstar,rcstar, nslope3, fnorm, zcin     
            
      integer  :: inc, jk, jl, iazi
!
      real(kind=wp) :: zang, zang1, znorm
      real(kind=wp) :: zx1, zx2, ztx, zdx, zxran, zxmin, zxmax, zx, zpexp
      real(kind=wp) :: fpc, fpc_dc
      real(kind=wp) :: ae1,ae2
     if( nwaves == 0) then
!
!     redefine from the deafault
!
       nwav   = incdim
       nazd   = iazidim
       nst    = 0
       eff    = 1.0_wp      
       gw_eff = eff
     else
!
! from input.nml
!
       nwav   = nwaves
       nazd   = nazdir
       nst    = nstoch
       gw_eff = effac
     endif  

      
      v_kxw  = kxw ; v_kxw2 = v_kxw*v_kxw 
      rv_kxw = 1.0_wp/v_kxw

      allocate ( zci(nwav),  zci4(nwav), zci3(nwav),zci2(nwav), zdci(nwav)  )
      allocate ( zcosang(nazd), zsinang(nazd) )
      allocate (lzmet(nwav), czmet(nwav), mkzmet(nwav), dczmet(nwav), dmkz(nwav) )   

!      if (me == master) then
!          print *, 'ugwp_v1/v0: init_gw_wmsdis_control '
!  
!          print *, 'ugwp_v1/v0: WMS_DIS launch layer ',    ilaunch
!          print *, 'ugwp_v1/v0: WMS_DIS tot_mflux in mpa', tamp_mpa*1000.
!          print *, 'ugwp_v1/v0: WMS_DIS lhmet in km   ' ,   lhmet*1.e-3            
!        endif

       zpexp = gptwo * 0.5_wp                    ! gptwo=2 , zpexp = 1.

!
!      set up azimuth directions and some trig factors
!
!
       zang = pi2 / real(nazd, kind=wp)

! get normalization factor to ensure that the same amount of momentum
! flux is directed (n,s,e,w) no mater how many azimuths are selected.
!
       znorm = 0.0_wp
       do iazi=1, nazd
         zang1         = real(iazi-1, kind=wp)*zang
         zcosang(iazi) = cos(zang1)
         zsinang(iazi) = sin(zang1)
         znorm         = znorm + abs(zcosang(iazi))
       enddo
!      zaz_fct = 1.0_wp
       zaz_fct = 2.0_wp / znorm            ! correction factor for azimuthal sums

!       define coordinate transform for "Ch"    ....x = 1/c stretching transform
!       ----------------------------------------------- 
! 
! x=1/Cphase transform
! Scinocca 2003.    x = 1/c stretching transform
!
          zxmax = 1.0_wp / zcimin
          zxmin = 1.0_wp / zcimax
          zxran = zxmax - zxmin
          zdx   = zxran / real(nwav-1, kind=wp)                            ! dkz
!
          ae1=zxran/zgam
          zx1   = zxran/(exp(ae1)-1.0_wp )                    ! zgam =1./4.
          zx2   = zxmin - zx1

!
! computations for zci =1/zx, stretching "accuracy" is not "accurate" spectra transform
! it represents additional "empirical" redistribution of "spectral" mode in C-space
!                               
         zms = pi2 / lzstar

        do inc=1, nwav
          ztx = real(inc-1, kind=wp)*zdx+zxmin
          ae1 = (ztx-zxmin)/zgam
          zx  = zx1*exp(ae1)+zx2                            !eq.(29-30),Scinocca-2003
          zci(inc)  = 1.0_wp /zx                                             !
          zdci(inc) = zci(inc)**2*(zx1/zgam)*exp(ae1)*zdx   !
          zci4(inc) = (zms*zci(inc))**4
          zci2(inc) = (zms*zci(inc))**2
          zci3(inc) = (zms*zci(inc))**3
        enddo
!
!
! alternatuve lzmax-lzmin without  x=1/c transform
!
!
        if (version == 1) then 

        dlzmet = (lzmax-lzmin)/ real(nwav-1, kind=wp)    
         do inc=1, nwav
          lzmet(inc) = lzmin + real(inc-1, kind=wp)*dlzmet
          mkzmet(inc) = pi2/lzmet(inc)
          zci(inc) =lzmet(inc)/(pi2/bnfix)
          zci4(inc) = (zms*zci(inc))**4
          zci2(inc) = (zms*zci(inc))**2
          zci3(inc) = (zms*zci(inc))**3
          
         enddo

         zdx = (zci(nwav)-zci(1))/ real(nwav-1, kind=wp)  
         do inc=1, nwav                  
          zdci(inc) = zdx     
         enddo   
         
         cstar = bnfix/zms
         rcstar = 1.0_wp/cstar
        ENDIF                               !   if (version == 1) then    
         
        RETURN           
!===================  Diag prints after return  ====================
         if (me == master) then
           print *
           print *,  'ugwp_v0: zcimin=' , zcimin
           print *,  'ugwp_v0: zcimax=' , zcimax
           print *,  'ugwp_v0: zgam=   ', zgam                
           print *
           
!          print *, ' ugwp_v1  nslope=', nslope
           print *
           print *,  'ugwp_v1: zcimin/zci=' , maxval(zci)
           print *,  'ugwp_v1: zcimax/zci=' , minval(zci)
           print *,  'ugwp_v1: cd_crit=', ucrit            
           print *,  'ugwp_v1: launch_level',  ilaunch
           print *, ' ugwp_v1  lzstar=', lzstar
           print *, ' ugwp_v1  nslope=', nslope

           print *
           nslope3=nslope+3.0_wp
         do inc=1, nwav           
           zcin =zci(inc)*rcstar
           fpc =  rcstar*(zcin*zcin)/(1.0_wp + zcin**nslope3)                    
           fpc_dc = fpc * zdci(inc)
           write(6,111)  inc, zci(inc), zdci(inc),ucrit, fpc, fpc_dc, 6.28e-3_wp/bnfix*zci(inc)
          enddo
         endif
         

         
 111     format( 'wms-zci', i4, 7 (3x, F8.3))

     end subroutine initsolv_wmsdis
!
!
  end module ugwp_wmsdis_init
