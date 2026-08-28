_pbl_wp!> \file satmedmfvdifq.F

!> This file contains the CCPP-compliant SATMEDMF scheme (updated version) which
!! computes subgrid vertical turbulence mixing using scale-aware TKE-based moist
!! eddy-diffusion mass-flux (TKE-EDMF) parameterization (by Jongil Han).

!! if(tte_edmf=.true.), the TKE-EDMF parameterization becomes
!! TTE(total turbulent energy)-based moist (TTE-EDMF) parameterization
!!
      module satmedmfvdifq
      use mfpbltq_mod
      use tridi_satmedmf_mod
      use mfscuq_mod
      !PCC_CANOPY_utilities
      use canopy_utils_mod

      contains

!> \defgroup module_satmedmfvdifq GFS TKE-EDMF PBL Module
!! This file contains the CCPP-compliant SATMEDMF scheme (updated version) which
!! computes subgrid vertical turbulence mixing using scale-aware TKE-based moist
!! eddy-diffusion mass-flux (TKE-EDMF) parameterization (by Jongil Han).
!> @{
!! \brief This subroutine contains all of the logic for the
!! scale-aware TKE-based moist eddy-diffusion mass-flux (TKE-EDMF, updated version) scheme.
!! For local turbulence mixing, a TKE closure model is used.
!! Updated version of satmedmfvdif.f (May 2019) to have better low level
!! inversion, to reduce the cold bias in lower troposphere,
!! and to reduce the negative wind speed bias in upper troposphere
!!
!! Incorporate the LES-based changes for TC simulation
!! (Chen et al.,2022 \cite Chen_2022)
!! with additional improvements on MF working with Cu schemes
!! Xiaomin Chen, 5/2/2022
!!
!! Incorporate the TTE-EDMF; if (tte_edmf=.true.),
!!  TKE-EDMF scheme becomes TTE-EDMF scheme and the variable 'te'
!!  is read as TTE; if (tte_edmf=.false.), the variable 'te' is
!!  read as TKE, 5/22/2025
!!
!!
!> \section arg_table_satmedmfvdifq_init Argument Table
!! \htmlinclude satmedmfvdifq_init.html
!!
      subroutine satmedmfvdifq_init (satmedmf,                          &
     &                               isatmedmf,isatmedmf_vdifq,         &
     &                               errmsg,errflg)

      logical, intent(in   ) :: satmedmf
      integer, intent(in) :: isatmedmf,isatmedmf_vdifq

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

! Consistency checks
      if (.not. satmedmf) then
        write(errmsg,fmt='(*(a))') 'Logic error: satmedmf = .false.'
        errflg = 1
        return
      end if

      if (.not. isatmedmf==isatmedmf_vdifq) then
        write(errmsg,fmt='(*(a))') 'Logic error: satmedmfvdif is ',
     &                   'called, but isatmedmf/=isatmedmf_vdifq.'
        errflg = 1
        return
      end if

      end subroutine satmedmfvdifq_init

!> \section arg_table_satmedmfvdifq_run Argument Table
!! \htmlinclude satmedmfvdifq_run.html
!!
!!\section gen_satmedmfvdifq GFS satmedmfvdifq General Algorithm
!! satmedmfvdifq_run() computes subgrid vertical turbulence mixing
!! using the scale-aware TKE-based moist eddy-diffusion mass-flux (EDMF) parameterization of
!! Han and Bretherton (2019) \cite Han_2019 .
!! -# The local turbulent mixing is represented by an eddy-diffusivity scheme which 
!! is a function of a prognostic TKE.
!! -# For the convective boundary layer, nonlocal transport by large eddies 
!! (mfpbltq.f), is represented using a mass flux approach (Siebesma et al.(2007) \cite Siebesma_2007 ).
!! -# A mass-flux approach is also used to represent the stratocumulus-top-induced turbulence
!! (mfscuq.f).
!! \section detail_satmedmfvidfq GFS satmedmfvdifq Detailed Algorithm
      subroutine satmedmfvdifq_run(im,km,ntrac,ntcw,ntrw,               &
     &     ntiw,ntke,grav,pi,rd,cp,rv,hvap,hfus,fv,eps,epsm1,              &
!The following three variables are for SA-3D-TKE
     &     def_1,def_2,def_3,sa3dtke,dku3d_h,dku3d_e,                   &
     &     rtg,u1,v1,t1,q1,usfco,vsfco,use_oceanuv,                     &
     &     swh,hlw,xmu,garea,zvfun,sigmaf,                              &
     &     psk,rbsoil,zorl,u10m,v10m,fm,fh,                             &
     &     tsea,heat,evap,stress,spd1,kpbl,                             &
     &     prsi,del,prsl,prslk,phii,phil,delt,tte_edmf,                 &
     &     dspheat,dusfc,dvsfc,dtsfc,dqsfc,hpbl,dkt,dku,tkeh,           &
     &     kinver,xkzm_m,xkzm_h,xkzm_s,dspfac,bl_upfr,bl_dnfr,          &
     &     rlmx,elmx,sfc_rlm,tc_pbl,use_lpt,                            &
     &     do_canopy, cplaqm, claie, cfch, cfrt, cclu, cpopu,           &
     &     ntqv,dtend,dtidx,index_of_temperature,index_of_x_wind,       &
     &     index_of_y_wind,index_of_process_pbl,gen_tend,ldiag3d,       &
     &     ten_t, ten_u, ten_v,                                         &
     &     errmsg,errflg)

!
      use machine  , only : kind_phys
      use funcphys , only : fpvs
      use mo_pbl_kind, only: pbl_wp
!
      implicit none
!
!----------------------------------------------------------------------
      integer, intent(in)  :: im, km, ntrac, ntcw, ntrw, ntiw,          &
     &                        ntke, ntqv
      integer, intent(in)  :: sfc_rlm
      integer, intent(in)  :: tc_pbl
      integer, intent(in)  :: use_lpt
      integer, intent(in)  :: kinver(:)
      integer, intent(out) :: kpbl(:)
      logical, intent(in)  :: gen_tend,ldiag3d
!
      real(kind=kind_phys), intent(in) :: grav,pi,rd,cp,rv,hvap,hfus,fv,   &
     &                                    eps,epsm1
      real(kind=kind_phys), intent(in) :: delt, xkzm_m, xkzm_h, xkzm_s
      real(kind=kind_phys), intent(in) :: dspfac, bl_upfr, bl_dnfr
      real(kind=kind_phys), intent(in) :: rlmx, elmx
!PCC CANOPY------------------------------------
      logical, intent(in)  :: do_canopy, cplaqm
      real(kind=kind_phys), optional, intent(in)  :: claie(:), cfch(:), &
     &                     cfrt(:), cclu(:), cpopu(:)
      !----------------------------------------------
      real(kind=kind_phys), intent(inout) :: rtg(:,:,:),    tkeh(:,:)
      real(kind=kind_phys), intent(in) ::                               &
     &                     u1(:,:),     v1(:,:),                        &
     &                     usfco(:),    vsfco(:),                       &
     &                     t1(:,:),     q1(:,:,:),                      &
!The following two variables are for SA-3D-TKE
     &                     def_1(:,:),  def_2(:,:), def_3(:,:),         &
     &                     swh(:,:),    hlw(:,:),                       &
     &                     xmu(:),       garea(:),                      &
     &                     zvfun(:),     sigmaf(:),                     &
     &                     psk(:),       rbsoil(:),                     &
     &                     zorl(:),      tsea(:),                       &
     &                     u10m(:),      v10m(:),                       &
     &                     fm(:),        fh(:),                         &
     &                     evap(:),      heat(:),                       &
     &                     stress(:),    spd1(:),                       &
     &                     prsi(:,:), del(:,:),                         &
     &                     prsl(:,:),   prslk(:,:),                     &
     &                     phii(:,:), phil(:,:)
      real(kind=kind_phys), intent(inout), dimension(:,:,:), optional ::&
     &     dtend
      integer, intent(in) :: dtidx(:,:), index_of_temperature,          &
     &       index_of_x_wind, index_of_y_wind, index_of_process_pbl
      logical, intent(in) :: use_oceanuv
      real(kind=kind_phys), intent(out) ::                              &
     &                     dusfc(:),     dvsfc(:),                      &
     &                     dtsfc(:),     dqsfc(:),                      &
     &                     hpbl(:)
      real(kind=kind_phys), intent(out) ::                              &
     &                     dkt(:,:), dku(:,:)

!
      logical, intent(in)  :: sa3dtke !flag for SA-3D-TKE scheme
!
!  flag for tke dissipative heating
      logical, intent(in)  :: dspheat
      real(kind=kind_phys), intent(out) :: ten_t(:,:),                  &
     &                                     ten_u(:,:), ten_v(:,:)
!  flag for TTE-EDMF scheme
      logical, intent(in)  :: tte_edmf
!     
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!For passing dku to the dyn_core (SA-3D-TKE scheme)
       real(kind=kind_phys), intent(out) ::
     &            dku3d_h(:,:),dku3d_e(:,:)

!
!----------------------------------------------------------------------
!***
!***  local variables
      real(kind=pbl_wp) spd1_m
!***
      integer i,is,k,n,ndt,km1,kmpbl,kmscu,ntrac1,idtend
      integer kps,kbx,kmx
      integer lcld(im),kcld(im),krad(im),mrad(im)
      integer kx1(im), kb1(im), kpblx(im)
!
      real(kind=pbl_wp) u1_loc(im,km), v1_loc(im,km),
     &                  t1_loc(im,km), q1_loc(im,km,ntrac),
     &                  hpbl_loc(im), au_loc(im,km),
     &                  f1_loc(im,km), f2_loc(im,km*max(1,ntrac-1))

      real(kind=pbl_wp) te(im,km),   tei(im,km-1), tke(im,km),
     &                     tteh(im,km), tesq(im,km-1),e2(im,0:km)
!
      real(kind=pbl_wp) theta(im,km),thvx(im,km),  thlvx(im,km),
     &                     qlx(im,km),  thetae(im,km),thlx(im,km),
     &                     slx(im,km),  svx(im,km),   qtx(im,km),
     &                     tvx(im,km),  pix(im,km),   radx(im,km-1),
     &                     dkq(im,km-1)
!
      real(kind=pbl_wp) plyr(im,km), rhly(im,km),  cfly(im,km),
     &                     qstl(im,km)
!
      real(kind=pbl_wp) dtdz1(im), gdx(im),
     &                     phih(im),  phim(im),    phihs(im),
     &                     phims(im), prn(im,km-1),
     &                     rbdn(im),  rbup(im),    thermal(im),
     &                     ustar(im), wstar(im),   hpblx(im),
     &                     ust3(im),  wst3(im),    rho_a(im),
     &                     z0(im),    crb(im),     tkemean(im),
     &                     hgamt(im), hgamq(im),
     &                     wscale(im),vpert(im),   thvs(im),
     &                     zol(im),   sflux(im),   ris(im),
     &                     sumx(im),  tx1(im),     tx2(im)
!
      real(kind=pbl_wp) radmin(im)
!
      real(kind=pbl_wp) zi(im,km+1),  zl(im,km),   zm(im,km),
     &                     xkzo(im,km-1),xkzmo(im,km-1),
     &                     xkzm_hx(im),  xkzm_mx(im), tkmnz(im,km-1),
     &                     rdzt(im,km-1),rlmnz(im,km),
     &                     al(im,km-1),  ad(im,km),   au(im,km-1),
     &                     f1(im,km),    f2(im,km*(ntrac-1))
!
      real(kind=pbl_wp) elm(im,km),   ele(im,km),
     &                     ckz(im,km),   chz(im,km),
     &                     diss(im,km-1),prod(im,km-1),
     &                     bf(im,km-1),  shr2(im,km-1), wush(im,km),
     &                     xlamue(im,km-1), xlamde(im,km-1),
     &                     gotvx(im,km), rlam(im,km-1)
!
!   variables for updrafts (thermals)
!
      real(kind=pbl_wp) tcko(im,km),  qcko(im,km,ntrac),
     &                     ucko(im,km),  vcko(im,km),
     &                     buou(im,km),  xmf(im,km)
!
!   variables for stratocumulus-top induced downdrafts
!
      real(kind=pbl_wp) tcdo(im,km),  qcdo(im,km,ntrac),
     &                     ucdo(im,km),  vcdo(im,km),
     &                     buod(im,km),  xmfd(im,km)
!
!  variables for Total Variation Diminishing (TVD) flux-limiter scheme
!
      real(kind=pbl_wp) e_half(im,km-1), e_diff(im,0:km-1),
     &                     q_half(im,km-1,ntrac-1),
     &                     qh(im,km-1,ntrac-1),
     &                     q_diff(im,0:km-1,ntrac-1)
      real(kind=pbl_wp) rrkp, phkp
      real(kind=pbl_wp) tsumn(im), tsump(im), rtnp(im)
      real(kind=pbl_wp) sfcpbl(im), vez0fun(im)
!
      logical  pblflg(im), sfcflg(im), flg(im)
      logical  scuflg(im), pcnvflg(im)
      logical  mlenflg
!
!  pcnvflg: true for unstable pbl
!
      real(kind=pbl_wp) aphi16,  aphi5,
     &                     wfac,    cfac,
     &                     gamcrt,  gamcrq, sfcfrac,
!    &                     conq,    cont,   conw,
     &                     dsdz2,   dsdzt,  dkmax,
     &                     dsig,    dt2,    dtodsd,
     &                     dtodsu,  g,      factor, dz,
     &                     gocp,    gravi,  zol1,   zolcru,
     &                     buop,    shrp,   dtn,
     &                     prnum,   prmax,  prmin,  prtke,
     &                     prscu,   pr0,    ri,
     &                     dw2,     dw2min, zk,
     &                     elmfac,  elefac, dspmax,
     &                     alp,     clwt,   cql,
     &                     f0,      robn,   crbmin, crbmax,
     &                     es,      qs,     value,  onemrh,
     &                     cfh,     gamma,  elocp,  el2orc,
     &                     epsi,    beta,   chx,    cqx,
     &                     rdt,     rdz,    qmin,   qlmin,
     &                     rimin,   rbcr,   rbint,  tdzmin,
     &                     rlmn,    rlmn0,  rlmn1,  rlmn2,
     &                     qtend,
     &                     zfac,    zfmin,  vk,     spdk2,
     &                     tkmin,   tkbmx,  disste, xkgdx,
     &                     xkinv1,  xkinv2,
     &                     zlup,    zldn,   cs0,    csmf,
     &                     tem,     tem1,   tem2,   tem3,
     &                     ptem,    ptem0,  ptem1,  ptem2
!
!The following variables are for SA-3D-TKE
      integer kk
      real(kind=pbl_wp) thetal(im,km),dku_les(im,km),dkt_les(im,km),
     &                     elmh(im,km),ele_les(im,km),pftke(im),
     &                     dkq_les(im,km),pfl(im),pfdx(im),
     &                     dku_h(im,km),dkq_h(im,km),
     &                     elmhfac,elmhmx,ckh,elm_les,
     &                     cpl1,cpl2,cpl3,cpl4,cpl5,cpl6,
     &                     cptke1,cptke2,cptke3
      integer ktkemax(im)
      real(kind=pbl_wp) tkemax(im),scl(im)
      real(kind=pbl_wp) sclmax,sclmin,dkmaxles
! end of SA-3D-TKE variables
!
      real(kind=pbl_wp) slfac
!
      real(kind=pbl_wp) vegflo, vegfup, z0lo, z0up, vc0, zc0
!
      real(kind=pbl_wp) ck0, ck1, ch0, ch1, ce0, rchck
!
!
      real(kind=pbl_wp) epotte
!
      real(kind=pbl_wp) qlcr, zstblmax, hcrinv
!
      real(kind=pbl_wp) h1

      real(kind=pbl_wp) bfac, mffac

      real(kind=pbl_wp) qice(im,km),qliq(im,km)

!PCC_CANOPY------------------------------------
      integer COUNTCAN,KCAN
      integer kount 
      real(kind=pbl_wp) FCH,       MOL,   HOL, TLCAN,
     &                     SIGMACAN,  RRCAN, BBCAN,
     &                     AACAN,     ZCAN,  ZFL, BOTCAN,
     &                     ZZ_INT,
     &                     EDDYVEST1, EDDYVEST_INT

      ! in canopy eddy diffusivity [ m**2/s ]
      real(kind=pbl_wp), allocatable  :: EDDYVESTX ( : )
      ! in canopy layer [m]
      real(kind=pbl_wp), allocatable  :: ZCANX     ( : )
      real(kind=pbl_wp), allocatable  :: ZOOOX     ( : )
      ! Declare local maximum canopy layers
      integer, parameter :: MAXCAN = 1000
      integer, parameter :: mvt   = 30 ! use 30 instead of 27
      !Based on MODIS IGBP 20 Category Dataset
      real(kind=pbl_wp) :: fch_table(mvt)         !< top of canopy (m)
       data   ( fch_table(i),i=1,mvt)   /
     &          20.0_pbl_wp, 20.0_pbl_wp, 18.0_pbl_wp, 16.0_pbl_wp,
     &          16.0_pbl_wp, 1.10_pbl_wp, 1.10_pbl_wp, 13.0_pbl_wp,
     &          10.0_pbl_wp, 1.00_pbl_wp, 5.00_pbl_wp, 2.00_pbl_wp,
     &          15.0_pbl_wp, 1.50_pbl_wp, 0.00_pbl_wp, 0.00_pbl_wp,
     &          0.00_pbl_wp, 4.00_pbl_wp, 2.00_pbl_wp, 0.50_pbl_wp,
     &          0.00_pbl_wp, 0.00_pbl_wp, 0.00_pbl_wp, 0.00_pbl_wp,
     &          0.00_pbl_wp, 0.00_pbl_wp, 0.00_pbl_wp, 0.00_pbl_wp,
     &          0.00_pbl_wp, 0.00_pbl_wp /
!----------------------------------------------
      parameter(bfac=100.0_pbl_wp)
      parameter(wfac=7.0_pbl_wp)
      parameter(gamcrt=3.0_pbl_wp,gamcrq=0.0_pbl_wp,sfcfrac=0.1_pbl_wp)
      parameter(vk=0.4_pbl_wp,rimin=-100.0_pbl_wp,slfac=0.1_pbl_wp)
      parameter(rbcr=0.25_pbl_wp,zolcru=-0.02_pbl_wp,
     &          tdzmin=1.e-3_pbl_wp)
      parameter(rlmn=30.0_pbl_wp,rlmn0=5.0_pbl_wp,rlmn1=5.0_pbl_wp,
     &          rlmn2=10.0_pbl_wp)
      parameter(prmin=0.25_pbl_wp)
      parameter(pr0=1.0_pbl_wp,prtke=1.0_pbl_wp)
      parameter(f0=1.e-4_pbl_wp,crbmin=0.15_pbl_wp,crbmax=0.35_pbl_wp)
      parameter(tkmin=1.e-9_pbl_wp,tkbmx=0.2_pbl_wp,dspmax=10.0_pbl_wp)
      parameter(qmin=1.e-8_pbl_wp,qlmin=1.e-12_pbl_wp,
     &          zfmin=1.e-8_pbl_wp)
      parameter(aphi5=5.0_pbl_wp,aphi16=16.0_pbl_wp)
      parameter(elmfac=1.0_pbl_wp,elefac=1.0_pbl_wp,cql=100.0_pbl_wp)
      parameter(dw2min=1.e-4_pbl_wp,dkmax=1000.0_pbl_wp,
     &          xkgdx=1000.0_pbl_wp)
      parameter(qlcr=3.5e-5_pbl_wp,zstblmax=2500.0_pbl_wp)
      parameter(xkinv1=0.15_pbl_wp,xkinv2=0.3_pbl_wp)
      parameter(h1=0.33333333_pbl_wp,hcrinv=250.0_pbl_wp)
      parameter(vegflo=0.1_pbl_wp,vegfup=1.0_pbl_wp,z0lo=0.1_pbl_wp,
     &          z0up=1.0_pbl_wp)
      parameter(vc0=1.0_pbl_wp,zc0=1.0_pbl_wp)
      parameter(cs0=0.4_pbl_wp,csmf=0.5_pbl_wp)
      parameter(rchck=1.5_pbl_wp,ndt=20)
      !The following variables are for SA-3D-TKE
      parameter(cpl1=0.280_pbl_wp,cpl2=0.870_pbl_wp,cpl3=0.913_pbl_wp)
      parameter(cpl4=0.153_pbl_wp,cpl5=0.278_pbl_wp,cpl6=0.720_pbl_wp)
      parameter(cptke1=0.07_pbl_wp,cptke2=0.142_pbl_wp,cptke3=0.071_pbl_wp)
      parameter(dkmaxles=300.0_pbl_wp,sclmin=500.0_pbl_wp,
     &          sclmax=2500.0_pbl_wp)
      parameter(elmhfac=1.5_pbl_wp,elmhmx=1000.0_pbl_wp,ckh=0.4_pbl_wp)

!PCC_CANOPY------------------------------------
      if (do_canopy) then
        if(.not.allocated(EDDYVESTX))
     &    allocate( EDDYVESTX    ( MAXCAN ) )
        if(.not.allocated(ZCANX))
     &    allocate( ZCANX    ( MAXCAN ) )
        if(.not.allocated(ZOOOX)) 
     &    allocate( ZOOOX    ( MAXCAN ) )
      endif
!----------------------------------------------
      if (tc_pbl == 0) then
        ck0 = 0.4_pbl_wp
        ch0 = 0.4_pbl_wp
        ce0 = 0.4_pbl_wp
      else if (tc_pbl == 1) then
        ck0 = 0.55_pbl_wp
        ch0 = 0.55_pbl_wp
        ce0 = 0.12_pbl_wp
      endif

      if(tte_edmf) then
        cfac  = 3.0_pbl_wp
        prmax = 6.0_pbl_wp
        prscu = 0.4_pbl_wp
        ck1   = 0.16_pbl_wp
        ch1   = 0.16_pbl_wp
      else
        cfac  = 4.5_pbl_wp
        prmax = 4.0_pbl_wp
        prscu = 0.67_pbl_wp
        ck1   = 0.15_pbl_wp
        ch1   = 0.15_pbl_wp
      endif

      gravi = 1.0_pbl_wp / real(grav, kind=pbl_wp)
      g = real(grav, kind=pbl_wp)
      gocp = g / real(cp, kind=pbl_wp)
      elocp  = real(hvap, kind=pbl_wp) / real(cp, kind=pbl_wp)
      el2orc = real(hvap, kind=pbl_wp) * real(hvap, kind=pbl_wp)
     &       / (real(rv, kind=pbl_wp) * real(cp, kind=pbl_wp))
!
!************************************************************************
! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!> ## Compute preliminary variables from input arguments
      dt2  = real(delt, kind=pbl_wp)
      rdt = 1.0_pbl_wp / dt2
!
! the code is written assuming ntke=ntrac
! if ntrac > ntke, the code needs to be modified
!
      ntrac1 = ntrac - 1
      km1 = km - 1
      kmpbl = km / 2
      kmscu = km / 2
!>  - Compute physical height of the layer centers and interfaces from
!! the geopotential height (\p zi and \p zl)
      do k=1,km
        do i=1,im
          zi(i,k) = real(phii(i,k), kind=pbl_wp) * gravi
          zl(i,k) = real(phil(i,k), kind=pbl_wp) * gravi
          xmf(i,k) = 0.0_pbl_wp
          xmfd(i,k) = 0.0_pbl_wp
          buou(i,k) = 0.0_pbl_wp
          buod(i,k) = 0.0_pbl_wp
          wush(i,k) = 0.0_pbl_wp
          ckz(i,k) = ck1
          chz(i,k) = ch1
          rlmnz(i,k) = rlmn0
        enddo
      enddo
      do i=1,im
        zi(i,km+1) = real(phii(i,km+1), kind=pbl_wp) * gravi
      enddo
      do k=1,km
        do i=1,im
          zm(i,k) = zi(i,k+1)
        enddo
      enddo
!>  - Compute horizontal grid size (\p gdx)
      do i=1,im
        gdx(i) = sqrt(real(garea(i), kind=pbl_wp))
      enddo
!>  - Initialize tke value at vertical layer centers and interfaces
!! from tracer (\p tke and \p tkeh)
      do k=1,km
        do i=1,im
          te(i,k) = max(real(q1(i,k,ntke), kind=pbl_wp), tkmin)
          tkeh(i,k) = 0.0_pbl_wp
          tteh(i,k) = 0.0_pbl_wp
        enddo
      enddo
      if(tte_edmf) then
        do k=1,km1
          do i=1,im
            tteh(i,k) = 0.5_pbl_wp * (te(i,k) + te(i,k+1))
          enddo
        enddo
      else
        do k = 1, km
          do i = 1, im
            tke(i,k) = te(i,k)
          enddo
        enddo
        do k=1,km1
          do i=1,im
            tkeh(i,k) = 0.5_pbl_wp * (tke(i,k) + tke(i,k+1))
          enddo
        enddo
      endif
!>  - Compute reciprocal of \f$ \Delta z \f$ (rdzt)
      do k = 1,km1
        do i=1,im
          rdzt(i,k) = 1.0_pbl_wp / (zl(i,k+1) - zl(i,k))
          prn(i,k)  = pr0
        enddo
      enddo
!
!>  - Compute reciprocal of pressure (tx1, tx2)

!>  - Compute minimum turbulent mixing length (rlmnz)

!>  - Compute background vertical diffusivities for scalars and momentum (xkzo and xkzmo)

!>  - set background diffusivities with xkzm_h & xkzm_m for gdx >= xkgdx and
!!    as a function of horizontal grid size for gdx < xkgdx
!!    \n  xkzm_hx = xkzm_h * (gdx / xkgdx)
!!    \n  xkzm_mx = xkzm_m * (gdx / xkgdx)
!
      do i=1,im
        kx1(i) = 1
        tx1(i) = 1.0_pbl_wp / real(prsi(i,1)
        tx2(i) = tx1(i)
        if(gdx(i) >= xkgdx) then
          xkzm_hx(i) = real(xkzm_h, kind=pbl_wp)
          xkzm_mx(i) = real(xkzm_m, kind=pbl_wp)
        else
          tem  = gdx(i) / xkgdx
          xkzm_hx(i) = real(xkzm_h, kind=pbl_wp) * tem
          xkzm_mx(i) = real(xkzm_m, kind=pbl_wp) * tem
        endif
      enddo
      do k = 1,km1
        do i=1,im
          xkzo(i,k)  = 0.0_pbl_wp
          xkzmo(i,k) = 0.0_pbl_wp
          if (k < kinver(i)) then
!                               minimum turbulent mixing length
            ptem      = real(prsi(i,k+1), kind=pbl_wp) * tx1(i)
            tem1      = 1.0_pbl_wp - ptem
            tem2      = tem1 * tem1 * 2.5_pbl_wp
            tem2      = min(1.0_pbl_wp, exp(-tem2))
            rlmnz(i,k)= rlmn * tem2
            rlmnz(i,k)= max(rlmnz(i,k), rlmn0)
!                               vertical background diffusivity
            tem2      = tem1 * tem1 * 10.0_pbl_wp
            tem2      = min(1.0_pbl_wp, exp(-tem2))
            xkzo(i,k) = xkzm_hx(i) * tem2
!                               vertical background diffusivity for
!                               momentum
            if (ptem >= real(xkzm_s, kind=pbl_wp)) then
              xkzmo(i,k) = xkzm_mx(i)
              kx1(i)     = k + 1
            else
              if (k == kx1(i) .and. k > 1) tx2(i) = 1.0_pbl_wp
     &           / real(prsi(i,k), kind=pbl_wp)
              tem1 = 1.0_pbl_wp - real(prsi(i,k+1), kind=pbl_wp)
     &             * tx2(i)
              tem1 = tem1 * tem1 * 5.0_pbl_wp
              xkzmo(i,k) = xkzm_mx(i) * min(1.0_pbl_wp, exp(-tem1))
            endif
          endif
        enddo
      enddo
!
!>  - Some output variables and logical flags are initialized
      do i = 1,im
         z0(i)    = 0.01_pbl_wp * real(zorl(i), kind=pbl_wp)
         rho_a(i) = real(prsl(i,1), kind=pbl_wp) / (real(rd,
     &              kind=pbl_wp) * real(t1(i,1), kind=pbl_wp)
     &            * (1.0_pbl_wp + real(fv, kind=pbl_wp)
     &            * max(real(q1(i,1,1), kind=pbl_wp), qmin)))
         dusfc(i) = 0.0_pbl_wp
         dvsfc(i) = 0.0_pbl_wp
         dtsfc(i) = 0.0_pbl_wp
         dqsfc(i) = 0.0_pbl_wp
         kpbl(i) = 1
         hpbl(i) = 0.0_pbl_wp
         kpblx(i) = 1
         hpblx(i) = 0.0_pbl_wp
         pfl(i)=1.0_pbl_wp
         pftke(i)=1.0_pbl_wp
         pblflg(i)= .true.
         sfcflg(i)= .true.
         if(rbsoil(i) > 0.0_pbl_wp) sfcflg(i) = .false.
         pcnvflg(i)= .false.
         scuflg(i)= .true.
         if(scuflg(i)) then
           radmin(i)= 0.0_pbl_wp
           mrad(i)  = km1
           krad(i)  = 1
           lcld(i)  = km1
           kcld(i)  = km1
         endif
      enddo
!
!>  - Compute a function for green vegetation fraction and surface roughness.
!! Entrainment rate in updraft is a function of vegetation fraction and surface
!! roughness length
!
      do i = 1,im
        tem = (real(sigmaf(i), kind=pbl_wp) - vegflo)
     &      / (vegfup - vegflo)
        tem = min(max(tem, 0.0_pbl_wp), 1.0_pbl_wp)
        tem1 = sqrt(tem)
        ptem = (z0(i) - z0lo) / (z0up - z0lo)
        ptem = min(max(ptem, 0.0_pbl_wp), 1.0_pbl_wp)
        vez0fun(i) = (1.0_pbl_wp + vc0 * tem1) * (1.0_pbl_wp + zc0
     &             * ptem)
      enddo
!
!>  - Compute \f$\theta\f$(theta), and \f$q_l\f$(qlx), \f$\theta_e\f$(thetae),
!! \f$\theta_v\f$(thvx),\f$\theta_{l,v}\f$ (thlvx) including ice water
      do k=1,km
        do i=1,im
          pix(i,k)   = real(psk(i), kind=pbl_wp) / real(prslk(i,k),
     &                 kind=pbl_wp)
          theta(i,k) = real(t1(i,k), kind=pbl_wp) * pix(i,k)
          qice(i,k)  = 0.0_pbl_wp
          qliq(i,k)  = 0.0_pbl_wp
          if(ntiw > 0) then
            tem = max(real(q1(i,k,ntcw), kind=pbl_wp), qlmin)
            qliq(i,k) = tem
            if(sa3dtke) then
              tem1 = max(real(q1(i,k,ntiw) + q1(i,k,5) + q1(i,k,6),
     &               kind=pbl_wp), qlmin) !for SA-3D-TKE
              qice(i,k) = tem1
            else
              tem1 = max(real(q1(i,k,ntiw), kind=pbl_wp), qlmin)
              qice(i,k) = tem1
            endif
            qlx(i,k) = tem + tem1
            ptem = real(hvap, kind=pbl_wp) * tem + real(hvap + hfus,
     &             kind=pbl_wp) * tem1
            slx(i,k) = real(cp * t1(i,k) + phil(i,k), kind=pbl_wp)
     &               - ptem
          else
            qlx(i,k) = max(real(q1(i,k,ntcw), kind=pbl_wp), qlmin)
            slx(i,k) = real(cp * t1(i,k) + phil(i,k), kind=pbl_wp)
     &               - real(hvap, kind=pbl_wp) * qlx(i,k)
            qliq(i,k) = qlx(i,k)
          endif
          tem2       = 1.0_pbl_wp + real(fv, kind=pbl_wp)
     &              * max(real(q1(i,k,1), kind=pbl_wp), qmin) - qlx(i,k)
          thvx(i,k)  = theta(i,k) * tem2
          tvx(i,k)   = real(t1(i,k), kind=pbl_wp) * tem2
          qtx(i,k)  = max(real(q1(i,k,1), kind=pbl_wp), qmin) + qlx(i,k)
          thlx(i,k)  = theta(i,k) - pix(i,k)*elocp*qlx(i,k)
          thlvx(i,k) = thlx(i,k) * (1.0_pbl_wp + real(fv, kind=pbl_wp)
     &               * qtx(i,k))
          svx(i,k)   = real(cp, kind=pbl_wp) * tvx(i,k)
          ptem1      = elocp * pix(i,k) * max(real(q1(i,k,1),
     &                 kind=pbl_wp), qmin)
          thetae(i,k)= theta(i,k) +  ptem1
!         gotvx(i,k) = g / tvx(i,k)
          gotvx(i,k) = g / thvx(i,k)
        enddo
      enddo
!
!>  - Compute an empirical cloud fraction based on
!! Xu and Randall (1996) \cite xu_and_randall_1996
      do k = 1, km
        do i = 1, im
          plyr(i,k) = 0.01_pbl_wp * real(prsl(i,k), kind=pbl_wp)   ! pa to mb (hpa)
!  --- ...  compute relative humidity
          es = 0.01_pbl_wp * real(fpvs(real(t1(i,k), kind=kind_phys)),
     &         kind=pbl_wp)       ! fpvs in pa
          qs  = max(qmin, real(eps, kind=pbl_wp) * es / (plyr(i,k)
     &        + real(epsm1, kind=pbl_wp) * es))
          rhly(i,k) = max(0.0_pbl_wp, min(1.0_pbl_wp, max(qmin,
     &                real(q1(i,k,1), kind=pbl_wp))/qs))
          qstl(i,k) = qs
        enddo
      enddo
!
      do k = 1, km
        do i = 1, im
          cfly(i,k) = 0.0_pbl_wp
          clwt = 1.0e-6_pbl_wp * (plyr(i,k) * 0.001_pbl_wp)
          if (qlx(i,k) > clwt) then
            onemrh = max(1.e-10_pbl_wp, 1.0_pbl_wp-rhly(i,k))
            tem1   = min(max((onemrh * qstl(i,k))**0.49_pbl_wp,
     &               0.0001_pbl_wp), 1.0_pbl_wp)
            tem1   = cql / tem1
            value = max(min(tem1 * qlx(i,k), 50.0_pbl_wp), 0.0_pbl_wp)
            tem2   = sqrt(sqrt(rhly(i,k)))
            cfly(i,k) = min(max(tem2 * (1.0_pbl_wp - exp(-value)),
     &                  0.0_pbl_wp), 1.0_pbl_wp)
          endif
        enddo
      enddo
!
!>  - Compute buoyancy modified by clouds
!
      do k = 1, km1
        do i = 1, im
          tem  = 0.5_pbl_wp * (svx(i,k) + svx(i,k+1))
          tem1 = 0.5_pbl_wp * real(t1(i,k) + t1(i,k+1), kind=pbl_wp)
          tem2 = 0.5_pbl_wp * (qstl(i,k) + qstl(i,k+1))
          cfh  = min(cfly(i,k+1),0.5_pbl_wp*(cfly(i,k)+cfly(i,k+1)))
          alp  = g / tem
          gamma = el2orc * tem2 / (tem1**2)
          epsi  = tem1 / elocp
          beta  = (1.0_pbl_wp + gamma * epsi * (1.0_pbl_wp + real(fv, kind=pbl_wp)))
     &          / (1.0_pbl_wp + gamma)
          chx   = cfh * alp * beta + (1.0_pbl_wp - cfh) * alp
          cqx   = cfh * alp * real(hvap, kind=pbl_wp) * (beta - epsi)
          cqx   = cqx + (1.0_pbl_wp - cfh) * real(fv, kind=pbl_wp) * g
          ptem1 = (slx(i,k+1) - slx(i,k)) * rdzt(i,k)
          ptem2 = (qtx(i,k+1) - qtx(i,k)) * rdzt(i,k)
          bf(i,k) = chx * ptem1 + cqx * ptem2
        enddo
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>  - Initialize diffusion coefficients to 0 and calculate the total
!! radiative heating rate (dku, dkt, radx)
      do k=1,km
        do i=1,im
          dku(i,k)  = 0.0_pbl_wp
          dkt(i,k)  = 0.0_pbl_wp
        enddo
      enddo
      do k=1,km1
        do i=1,im
          dkq(i,k)  = 0.0_pbl_wp
          tem       = zi(i,k+1) - zi(i,k)
          radx(i,k) = tem * real(swh(i,k) * xmu(i) + hlw(i,k),
     &                kind=pbl_wp)
        enddo
      enddo
!>  - Compute stable/unstable PBL flag (pblflg) based on the total
!! surface energy flux (\e false if the total surface energy flux
!! is into the surface)
      do i = 1,im
         sflux(i)  = real(heat(i), kind=pbl_wp) + real(evap(i)
     &             * fv, kind=pbl_wp) * theta(i,1)
         if(.not.sfcflg(i) .or. sflux(i) <= 0.0_pbl_wp)
     &     pblflg(i)=.false.
      enddo
!
!>  ## Calculate the PBL height
!!  The calculation of the boundary layer height follows Troen and Mahrt (1986) \cite troen_and_mahrt_1986 section 3. The approach is to find the level in the column where a modified bulk Richardson number exceeds a critical value.
!!  - Compute critical bulk Richardson number (\f$Rb_{cr}\f$) (crb)
!!  - For the unstable PBL, crb is a constant (0.25)
!!  - For the stable boundary layer (SBL), \f$Rb_{cr}\f$ varies
!!     with the surface Rossby number, \f$R_{0}\f$, as given by
!!     Vickers and Mahrt (2004) \cite Vickers_2004
!!     \f[
!!      Rb_{cr}=0.16(10^{-7}R_{0})^{-0.18}
!!     \f]
!!     \f[
!!      R_{0}=\frac{U_{10}}{f_{0}z_{0}}
!!     \f]
!!     where \f$U_{10}\f$ is the wind speed at 10m above the ground surface,
!!     \f$f_0\f$ the Coriolis parameter, and \f$z_{0}\f$ the surface roughness
!!     length. To avoid too much variation, we restrict \f$Rb_{cr}\f$ to vary
!!     within the range of 0.15~0.35
      do i = 1,im
        thvs(i) = real(tsea(i), kind=pbl_wp) * (1.0_pbl_wp + real(fv,
     &            kind=pbl_wp) *max(real(q1(i,1,1), kind=pbl_wp), qmin))
        if(pblflg(i)) then
!         thermal(i) = thvx(i,1)
          thermal(i) = thlvx(i,1)
          crb(i) = rbcr
        else
          thermal(i) = thvs(i)
          tem = sqrt(real(u10m(i)**2 + v10m(i)**2, kind=pbl_wp))
          tem = max(tem, 1.0_pbl_wp)
          robn = tem / (real(f0, kind=pbl_wp) * real(z0(i),
     &           kind=pbl_wp))
          tem1 = 1.e-7_pbl_wp * robn
          crb(i) = 0.16_pbl_wp * (tem1 ** (-0.18_pbl_wp))
          crb(i) = max(min(crb(i), crbmax), crbmin)
        endif
      enddo
!>  - Compute \f$\frac{\Delta t}{\Delta z}\f$ , \f$u_*\f$
      do i=1,im
         dtdz1(i)  = dt2 / (zi(i,2) - zi(i,1))
      enddo
!
      do i=1,im
         ustar(i) = sqrt(real(stress(i), kind=pbl_wp))
      enddo
!
!>  - Compute buoyancy \f$\frac{\partial \theta_v}{\partial z}\f$ (bf)
!! and the wind shear squared (shr2)
!
      do k = 1, km1
      do i = 1, im
         rdz  = rdzt(i,k)
!        bf(i,k) = gotvx(i,k)*(thvx(i,k+1)-thvx(i,k))*rdz
         dw2  = real(u1(i,k) - u1(i,k+1), kind=pbl_wp)**2
     &        + real(v1(i,k) - v1(i,k+1), kind=pbl_wp)**2
         shr2(i,k) = max(dw2,dw2min)*rdz*rdz
      enddo
      enddo
!
! Find first quess pbl height based on bulk richardson number (mrf pbl scheme)
!   and also for diagnostic purpose
!
      do i=1,im
         flg(i) = .false.
         rbup(i) = real(rbsoil(i), kind=pbl_wp)
      enddo
!>  - Given the thermal's properties and the critical Richardson number, 
!! a loop is executed to find the first level above the surface (kpblx) where 
!! the modified Richardson number is greater than the critical Richardson
!! number, using equation 10a from Troen and Mahrt (1996) \cite troen_and_mahrt_1986 
!! (also equation 8 from Hong and Pan (1996) \cite hong_and_pan_1996):
      do k = 1, kmpbl
      do i = 1, im
        if(.not.flg(i)) then
          rbdn(i) = rbup(i)
          if (tc_pbl == 0) then
            spdk2   = max(real(u1(i,k)**2 + v1(i,k)**2, kind=pbl_wp),
     &                1.0_pbl_wp)      
!           rbup(i) = (thvx(i,k)-thermal(i))*
!    &                (g*zl(i,k)/thvx(i,1))/spdk2
            rbup(i) = (thlvx(i,k) - thermal(i)) * (g * zl(i,k)
     &                / thlvx(i,1)) / spdk2

          else if (tc_pbl == 1) then
            spdk2   = real(u1(i,k) - u1(i,1), kind=pbl_wp)**2
     &              + real(v1(i,k) - v1(i,1), kind=pbl_wp)**2
     &              + bfac * ustar(i)**2
            rbup(i) = (thlvx(i,k) - thermal(i)) * (g * (zl(i,k)
     &              - zl(i,1)) / thlvx(i,1)) / spdk2
          endif
          kpblx(i) = k
          flg(i)  = rbup(i) > crb(i)
        endif
      enddo
      enddo
!>  - Once the level is found, some linear interpolation is performed to find
!! the exact height of the boundary layer top (where \f$R_{i} > Rb_{cr}\f$)
!! and the PBL height (hpbl and kpbl) and the PBL top index are saved.
      do i = 1,im
        if(kpblx(i) > 1) then
          k = kpblx(i)
          if(rbdn(i) >= crb(i)) then
            rbint = 0.0_pbl_wp
          elseif(rbup(i) <= crb(i)) then
            rbint = 1.0_pbl_wp
          else
            rbint = (crb(i)-rbdn(i))/(rbup(i)-rbdn(i))
          endif
          hpblx(i) = zl(i,k-1) + rbint * (zl(i,k) - zl(i,k-1))
          if(hpblx(i) < zi(i,kpblx(i)))
     &       kpblx(i) = kpblx(i) - 1
        else
          hpblx(i) = zl(i,1)
          kpblx(i) = 1
        endif
        hpbl(i) = hpblx(i)
        kpbl(i) = kpblx(i)
        if(kpbl(i) <= 1) pblflg(i)=.false.
      enddo
!
! update thermal at a level of slfac*hpbl for unstable pbl
!
      do i=1,im
        sfcpbl(i) = slfac * hpbl(i)
        kb1(i) = 1
        flg(i)  = .false.
        if(pblflg(i)) then
          flg(i)  = .true.
        endif
      enddo
      do k = 2, kmpbl
        do i=1,im
          if (flg(i) .and. zl(i,k) <= sfcpbl(i)) then
            kb1(i) = k
          else
            flg(i) = .false.
          endif
        enddo
      enddo
      do i=1,im
        if(pblflg(i)) kb1(i)=min(kb1(i),kpbl(i))
      enddo
!
! re-compute pbl height with the updated thermal
!
      do i=1,im
        flg(i)  = .true.
        if(pblflg(i) .and. kb1(i) > 1) then
          flg(i)  = .false.
          rbup(i) = real(rbsoil(i), kind=pbl_wp)
!         thermal(i) = thvx(i,kb1(i))
          thermal(i) = thlvx(i,kb1(i))
          kpblx(i) = kb1(i)
          hpblx(i) = zl(i,kb1(i))
        endif
      enddo
      do k = 2, kmpbl
      do i = 1, im
        if(.not.flg(i) .and. k > kb1(i)) then
          rbdn(i) = rbup(i)
          if (tc_pbl == 0) then
            spdk2   = max(real(u1(i,k)**2 + v1(i,k)**2, kind=pbl_wp),
     &                1.0_pbl_wp)      
!           rbup(i) = (thvx(i,k)-thermal(i))*
!    &                (g*zl(i,k)/thvx(i,1))/spdk2
            rbup(i) = (thlvx(i,k) - thermal(i)) * (g * zl(i,k)
     &                / thlvx(i,1)) / spdk2
          else if (tc_pbl == 1) then
            spdk2   = max(real((u1(i,k)-u1(i,1))**2+(v1(i,k)-v1(i,1))**2
     &              , kind=pbl_wp),1.0_pbl_wp) + bfac*ustar(i)**2
            rbup(i) = (thlvx(i,k) - thermal(i)) * (g * (zl(i,k) &
     &              - zl(i,1)) / thlvx(i,1)) / spdk2
          endif
          kpblx(i) = k
          flg(i)  = rbup(i) > crb(i)
        endif
      enddo
      enddo
      do i = 1,im
        if(pblflg(i) .and. kb1(i) > 1) then
          k = kpblx(i)
          if(rbdn(i) >= crb(i)) then
            rbint = 0.0_pbl_wp
          elseif(rbup(i) <= crb(i)) then
            rbint = 1.0_pbl_wp
          else
            rbint = (crb(i)-rbdn(i))/(rbup(i)-rbdn(i))
          endif
          hpblx(i) = zl(i,k-1) + rbint * (zl(i,k) - zl(i,k-1))
          if(hpblx(i) < zi(i,kpblx(i)))
     &       kpblx(i) = kpblx(i) - 1
          hpbl(i) = hpblx(i)
          kpbl(i) = kpblx(i)
        endif
      enddo
!
      if(.not.tte_edmf) then
!
!>  - Compute mean tke within pbl for TKE-EDMF
!
      do i = 1, im
        sumx(i) = 0.0_pbl_wp
        tkemean(i) = 0.0_pbl_wp
      enddo
      do k = 1, kmpbl
      do i = 1, im
        if(k < kpbl(i)) then
          dz = zi(i,k+1) - zi(i,k)   
          tkemean(i) = tkemean(i) + real(tke(i,k), kind=pbl_wp) * dz
          sumx(i) = sumx(i) + dz
        endif
      enddo
      enddo
      do i = 1, im
        if(tkemean(i) > 0.0_pbl_wp .and. sumx(i) > 0.0_pbl_wp) then
          tkemean(i) = tkemean(i) / sumx(i)
        endif
      enddo
!
      endif
!
!>  - Compute wind shear term as a sink term for updraft and downdraft
!!  velocity
!
      kps = max(kmpbl, kmscu)
      do k = 2, kps
      do i = 1, im
        dz = zi(i,k+1) - zi(i,k)
        tem = (0.5_pbl_wp * real(u1(i,k-1) - u1(i,k+1), kind=pbl_wp)
     &      / dz)**2
        tem1 = tem + (0.5_pbl_wp * real(v1(i,k-1) - v1(i,k+1),
     &         kind=pbl_wp) / dz)**2
        wush(i,k) = csmf * sqrt(tem1)
      enddo
      enddo
!
!> ## Compute Monin-Obukhov similarity parameters
!!  - Calculate the Monin-Obukhov nondimensional stability paramter, commonly
!!    referred to as \f$\zeta\f$ using the following equation from Businger et al.(1971) \cite businger_et_al_1971
!!    (eqn 28):
!!    \f[
!!    \zeta = Ri_{sfc}\frac{F_m^2}{F_h} = \frac{z}{L}
!!    \f]
!!    where \f$F_m\f$ and \f$F_h\f$ are surface Monin-Obukhov stability functions calculated in sfc_diff.f and 
!!    \f$L\f$ is the Obukhov length.
      do i=1,im
         zol(i) = max(real(rbsoil(i), kind=pbl_wp) * real(fm(i),
     &            kind=pbl_wp) * real(fm(i), kind=pbl_wp)
     &          / real(fh(i), kind=pbl_wp), rimin)
         if(sfcflg(i)) then
           zol(i) = min(zol(i),-zfmin)
         else
           zol(i) = max(zol(i),zfmin)
         endif
!>  - Calculate the nondimensional gradients of momentum and temperature (\f$\phi_m\f$ (phim) and \f$\phi_h\f$(phih)) are calculated using
!! eqns 5 and 6 from Hong and Pan (1996) \cite hong_and_pan_1996 depending on the surface layer stability:
!!   - For the unstable and neutral conditions:
!!     \f[
!!     \phi_m=(1-16\frac{0.1h}{L})^{-1/4}
!!     \phi_h=(1-16\frac{0.1h}{L})^{-1/2}
!!     \f]
!!   - For the stable regime
!!     \f[
!!     \phi_m=\phi_t=(1+5\frac{0.1h}{L})
!!     \f]
         zol1 = zol(i)*sfcfrac*hpbl(i)/zl(i,1)
         if(sfcflg(i)) then
           tem     = 1.0_pbl_wp / (1.0_pbl_wp - aphi16*zol1)
           phih(i) = sqrt(tem)
           phim(i) = sqrt(phih(i))
           tem1    = 1.0_pbl_wp / (1.0_pbl_wp - aphi16*zol(i))
           phihs(i) = sqrt(tem1)
           phims(i) = sqrt(phihs(i))
         else
           phim(i) = 1.0_pbl_wp + aphi5*zol1
           phih(i) = phim(i)
           phims(i) = 1.0_pbl_wp + aphi5*zol(i)
           phihs(i) = phims(i)
         endif
      enddo
!
!>  - The \f$z/L\f$ (zol) is used as the stability criterion for the PBL.Currently,
!!    strong unstable (convective) PBL for \f$z/L < -0.02\f$ and weakly and moderately
!!    unstable PBL for \f$0>z/L>-0.02\f$
!>  - Compute the velocity scale \f$w_s\f$ (wscale) (eqn 22 of Han et al. 2019). It
!!    is represented by the value scaled at the top of the surface layer:
!!    \f[
!!    w_s=(u_*^3+7\alpha\kappa w_*^3)^{1/3}
!!    \f]
!!    where \f$u_*\f$ (ustar) is the surface friction velocity,\f$\alpha\f$ is the ratio
!!    of the surface layer height to the PBL height (specified as sfcfrac =0.1), 
!!    \f$\kappa =0.4\f$ is the von Karman constant, and \f$w_*\f$ is the convective velocity
!!    scale defined as eqn23 of Han et al.(2019):
!!    \f[
!!    w_{*}=[(g/T)\overline{(w'\theta_v^{'})}_0h]^{1/3}
!!    \f]
      do i=1,im
        if(pblflg(i)) then
          if(zol(i) < zolcru) then
            pcnvflg(i) = .true.
          endif
          wst3(i) = gotvx(i,1)*sflux(i)*hpbl(i)
          wstar(i)= wst3(i)**h1
          ust3(i) = ustar(i)**3.0_pbl_wp
          wscale(i)=(ust3(i)+wfac*vk*wst3(i)*sfcfrac)**h1
          ptem = ustar(i)/aphi5
          wscale(i) = max(wscale(i),ptem)
        endif
      enddo
!
!>  ## The counter-gradient terms for temperature and humidity are calculated.
!! -  Equation 4 of Hong and Pan (1996) \cite hong_and_pan_1996 and are used to calculate the "scaled virtual temperature excess near the surface" (equation 9 in Hong and Pan (1996) \cite hong_and_pan_1996) for use in the mass-flux algorithm.
!
      do i = 1,im
         if(pcnvflg(i)) then
           hgamt(i) = real(heat(i), kind=pbl_wp)/wscale(i)
           hgamq(i) = real(evap(i), kind=pbl_wp)/wscale(i)
           vpert(i) = hgamt(i) + hgamq(i) *real(fv, kind=pbl_wp)
     &              * theta(i,1)
           vpert(i) = max(vpert(i),0.0_pbl_wp)
           tem = min(cfac*vpert(i),gamcrt)
           thermal(i)= thermal(i) + tem
         endif
      enddo
!
!  enhance the pbl height by considering the thermal excess
!     (overshoot pbl top)
!
      do i=1,im
         flg(i)  = .true.
         if(pcnvflg(i)) then
           flg(i)  = .false.
           rbup(i) = real(rbsoil(i), kind=pbl_wp)
         endif
      enddo
      do k = 2, kmpbl
      do i = 1, im
        if(.not.flg(i) .and. k > kb1(i)) then
          rbdn(i) = rbup(i)
          if (tc_pbl == 0) then
            spdk2   = max(real(u1(i,k)**2+v1(i,k)**2, kind=pbl_wp),
     &                1.0_pbl_wp)      
            rbup(i) = (thlvx(i,k) - thermal(i)) * (g * zl(i,k)
     &                / thlvx(i,1)) / spdk2
          else if (tc_pbl == 1) then
            spdk2   = max(real((u1(i,k) - u1(i,1))**2 + (v1(i,k)
     &              - v1(i,1))**2, kind=pbl_wp),1.0_pbl_wp)
     &              + bfac *ustar(i)**2
            rbup(i) = (thlvx(i,k) - thermal(i)) * (g * (zl(i,k) &
     &            - zl(i,1)) / thlvx(i,1)) /spdk2
          endif
          kpbl(i) = k
          flg(i)  = rbup(i) > crb(i)
        endif
      enddo
      enddo
      do i = 1,im
        if(pcnvflg(i)) then
           k = kpbl(i)
           if(rbdn(i) >= crb(i)) then
             rbint = 0.0_pbl_wp
           elseif(rbup(i) <= crb(i)) then
             rbint = 1.0_pbl_wp
           else
             rbint = (crb(i)-rbdn(i))/(rbup(i)-rbdn(i))
           endif
           hpbl(i) = zl(i,k-1) + rbint * (zl(i,k) - zl(i,k-1))
           if(hpbl(i) < zi(i,kpbl(i))) then
             kpbl(i) = kpbl(i) - 1
           endif
           if(kpbl(i) <= 1) then
              pcnvflg(i) = .false.
              pblflg(i) = .false.
           endif
        endif
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  compute tke using tte & ri for TTE-EDMF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      if(tte_edmf) then
!
      do i = 1, im
        tem = phims(i) * phims(i)
        ris(i) = zol(i) * phihs(i) / tem
        ris(i) = max(ris(i), rimin)
      enddo
      do k = 1, km1
        do i = 1, im
          ptem = sfcfrac*hpbl(i)
          if (zl(i,k) <= ptem) then
            ri = ris(i)
          else
            if(k == 1) then
              tem = gotvx(i,1) * (thlvx(i,1)-thvs(i))
              tem1 = tem / zl(i,1)
              tem1 = 0.5_pbl_wp * (tem1 + bf(i,1))
              ptem = max(real(u1(i,1)**2 + v1(i,1)**2, kind=pbl_wp),
     &               1.0_pbl_wp)
              ptem1 = ptem / (zl(i,1) * zl(i,1))
              ptem1 = 0.5_pbl_wp * (ptem1 + shr2(i,1))
              ri = max(tem1/ptem1, rimin)
            else
              tem1 = 0.5_pbl_wp * (bf(i,k-1) + bf(i,k))
              ptem1 = 0.5_pbl_wp * (shr2(i,k-1) + shr2(i,k))
              ri = max(tem1/ptem1, rimin)
            endif
          endif
          if(ri < 0.0_pbl_wp) then
            tem = 2.0_pbl_wp* ri - pr0
            epotte = ri / tem
          else
            tem = pr0 + 3.0_pbl_wp * ri
            epotte = ri / tem
          endif
          tke(i,k) = te(i,k) * (1.0_pbl_wp - epotte)
        enddo
      enddo
      do i=1,im
        tke(i,km) = tke(i,km1)
      enddo
!
!>  - Compute mean tke within pbl for TTE-EDMF
!
      do i = 1, im
        sumx(i) = 0.0_pbl_wp
        tkemean(i) = 0.0_pbl_wp
      enddo
      do k = 1, kmpbl
      do i = 1, im
        if(k < kpbl(i)) then
          dz = zi(i,k+1) - zi(i,k)     
          tkemean(i) = tkemean(i) + tke(i,k) * dz
          sumx(i) = sumx(i) + dz
        endif
      enddo
      enddo
      do i = 1, im
        if(tkemean(i) > 0.0_pbl_wp .and. sumx(i) > 0.0_pbl_wp) then
          tkemean(i) = tkemean(i) / sumx(i)
        endif
      enddo
!
      endif                      ! end of if(tte_edmf)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  look for stratocumulus
!> ## Determine whether stratocumulus layers exist and compute quantities
!!  - Starting at the PBL top and going downward, if the level is less than 2.5 km
!!    and \f$q_l\geq q_{lcr}\f$ then set kcld = k (find the cloud top index in the PBL.
!!    If no cloud water above the threshold is hound, \e scuflg is set to F. 
      do i=1,im
         flg(i)  = scuflg(i)
      enddo
      do k = 1, km1
        do i=1,im
          if(flg(i).and.zl(i,k) >= zstblmax) then
             lcld(i)=k
             flg(i)=.false.
          endif
      enddo
      enddo
      do i = 1, im
        flg(i)=scuflg(i)
      enddo
      do k = kmscu,1,-1
      do i = 1, im
        if(flg(i) .and. k <= lcld(i)) then
          if(qlx(i,k) >= qlcr) then
             kcld(i)=k
             flg(i)=.false.
          endif
        endif
      enddo
      enddo
      do i = 1, im
        if(scuflg(i) .and. kcld(i)==km1) scuflg(i)=.false.
      enddo
!>  - Starting at the PBL top and going downward, if the level is less
!!    than the cloud top, find the level of the minimum radiative heating
!!    rate wihin the cloud. If the level of the minimum is the lowest model 
!!    level or the minimum radiative heating rate is positive, then set
!!    scuflg to F.
      do i = 1, im
        flg(i)=scuflg(i)
      enddo
      do k = kmscu,1,-1
      do i = 1, im
        if(flg(i) .and. k <= kcld(i)) then
          if(qlx(i,k) >= qlcr) then
            if(radx(i,k) < radmin(i)) then
              radmin(i)=radx(i,k)
              krad(i)=k
            endif
          else
            flg(i)=.false.
          endif
        endif
      enddo
      enddo
      do i = 1, im
        if(scuflg(i) .and. krad(i) <= 1) scuflg(i)=.false.
        if(scuflg(i) .and. radmin(i)>=0.0_pbl_wp) scuflg(i)=.false.
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> ## Compute components for mass flux mixing by large thermals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  - If the PBL is convective, the updraft properties are initialized
!!    to be the same as the state variables.
      do i = 1, im
        hpbl_loc(i) = real(hpbl(i), kind=pbl_wp)
      enddo
      do k = 1, km
        do i = 1, im
          t1_loc(i,k) = real(t1(i,k), kind=pbl_wp)
          u1_loc(i,k) = real(u1(i,k), kind=pbl_wp)
          v1_loc(i,k) = real(v1(i,k), kind=pbl_wp)
          if(pcnvflg(i)) then
            tcko(i,k) = t1_loc(i,k)
            ucko(i,k) = u1_loc(i,k)
            vcko(i,k) = v1_loc(i,k)
          endif
          if(scuflg(i)) then
            tcdo(i,k) = t1_loc(i,k)
            ucdo(i,k) = u1_loc(i,k)
            vcdo(i,k) = v1_loc(i,k)
          endif
        enddo
      enddo
      do n = 1, ntrac1
      do k = 1, km
        do i = 1, im
          q1_loc(i,k,n) = real(q1(i,k,n), kind=pbl_wp)
          if(pcnvflg(i)) then
            qcko(i,k,n) = q1_loc(i,k,n)
          endif
          if(scuflg(i)) then
            qcdo(i,k,n) = q1_loc(i,k,n)
          endif
        enddo
      enddo
      enddo
!>  - Call mfpbltq(), which is an EDMF parameterization (Siebesma et al.(2007) \cite Siebesma_2007)
!!    to take into account nonlocal transport by large eddies. For details of the mfpbltq subroutine, step into its documentation ::mfpbltq
      call mfpbltq(im,im,km,kmpbl,ntcw,ntrac1,dt2,
     &    pcnvflg,zl,zm,q1_loc,t1_loc,u1_loc,v1_loc,plyr,pix,thlx,thvx,
     &    gdx,hpbl_loc,kpbl,vpert,buou,wush,tkemean,vez0fun,xmf,
     &    tcko,qcko,ucko,vcko,xlamue,real(bl_upfr, kind=pbl_wp))
!>  - Call mfscuq(), which is a new mass-flux parameterization for
!!    stratocumulus-top-induced turbulence mixing. For details of the mfscuq subroutine, step into its documentation ::mfscuq
      call mfscuq(im,im,km,kmscu,ntcw,ntrac1,dt2,
     &    scuflg,zl,zm,q1_loc,t1_loc,u1_loc,v1_loc,plyr,pix,
     &    thlx,thvx,thlvx,gdx,thetae,
     &    krad,mrad,radmin,buod,wush,tkemean,vez0fun,xmfd,
     &    tcdo,qcdo,ucdo,vcdo,xlamde,real(bl_dnfr, kind=pbl_wp))

      do i = 1, im
        if(pcnvflg(i)) then
          hpbl(i) = real(hpbl_loc(i), kind=kind_phys)
        endif
      enddo
      if (tc_pbl == 1) then
!>  - unify mass fluxes with Cu
        do i=1,im
          if(zol(i) > -0.5_pbl_wp) then
            do k = 1, km
              xmf(i,k) = 0.0_pbl_wp
            end do
          end if
        end do
!>  - taper off MF in high-wind conditions
        do i = 1,im
          tem = sqrt(real(u10m(i), kind=pbl_wp)**2+real(v10m(i),
     &          kind=pbl_wp)**2)           
          mffac = (1.0_pbl_wp - MIN(MAX(tem - 20.0_pbl_wp, 0.0_pbl_wp),
     &             10.0_pbl_wp)/10.0_pbl_wp)
          do k = 1, km
            xmf(i,k) = xmf(i,k)*mffac
            xmfd(i,k) = xmfd(i,k)*mffac
          enddo
        enddo
      endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> ## Compute Prandtl number \f$P_r\f$ (prn) and exchange coefficient varying with height
      do k = 1, kmpbl
        do i = 1, im
          if(k < kpbl(i)) then
            tem = phih(i)/phim(i)
            ptem = sfcfrac*real(hpbl(i), kind=pbl_wp)
            tem1 = max(zi(i,k+1)-ptem, 0.0_pbl_wp)
            tem2 = tem1 / (real(hpbl(i), kind=pbl_wp) - ptem)
            if(pcnvflg(i)) then
              tem = min(tem, pr0)
              prn(i,k) = tem + (pr0 - tem) * tem2
            else
              tem = max(tem, pr0)
              prn(i,k) = tem
            endif
            prn(i,k) = min(prn(i,k),prmax)
            prn(i,k) = max(prn(i,k),prmin)
!
            ckz(i,k) = ck0 + (ck1 - ck0) * tem2
            ckz(i,k) = max(min(ckz(i,k), ck0), ck1)
            chz(i,k) = ch0 + (ch1 - ch0) * tem2
            chz(i,k) = max(min(chz(i,k), ch0), ch1)
!
          endif
        enddo
      enddo
!
! Above a threshold height (hcrinv), the background vertical
! diffusivities & mixing length
!    in the inversion layers are set to much smaller values (xkinv1 &
!    rlmn1)
!
! Below the threshold height (hcrinv), the background vertical
! diffusivities & mixing length
!    in the inversion layers are increased with increasing roughness
!    length & vegetation fraction
!
      do k = 1,km1
        do i=1,im
          if(zi(i,k+1) > hcrinv) then
            tem1 = tvx(i,k+1)-tvx(i,k)
            if(tem1 >= 0.0_pbl_wp) then
              xkzo(i,k)  = min(xkzo(i,k), xkinv1)
              xkzmo(i,k) = min(xkzmo(i,k), xkinv1)
              rlmnz(i,k) = min(rlmnz(i,k), rlmn1)
            endif
          else
            tem1 = tvx(i,k+1)-tvx(i,k)
            if(tem1 > 0.0_pbl_wp) then
              ptem = xkzo(i,k) * real(zvfun(i), kind=pbl_wp)      
              xkzo(i,k) = min(max(ptem, xkinv2), xkzo(i,k))
              ptem = xkzmo(i,k) * real(zvfun(i), kind=pbl_wp)
              xkzmo(i,k) = min(max(ptem, xkinv2), xkzmo(i,k))
              ptem = rlmnz(i,k) * real(zvfun(i), kind=pbl_wp)
              rlmnz(i,k) = min(max(ptem, rlmn2), rlmnz(i,k))
            endif
          endif
        enddo
      enddo
      do k = 2,km1
        do i=1,im
          rlmnz(i,k) = 0.5_pbl_wp * (rlmnz(i,k-1) + rlmnz(i,k))
        enddo
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> ## Compute an asymtotic mixing length
!
      do k = 1, km1
        do i = 1, im
          zlup = 0.0_pbl_wp
          mlenflg = .true.
          e2(i,k) = max(2.0_pbl_wp*tke(i,k), 0.001_pbl_wp)
          do n = k, km1
            if(mlenflg) then
              dz = zl(i,n+1) - zl(i,n)
              tem1 = 2.0_pbl_wp*gotvx(i,n+1)*(thvx(i,k)-thvx(i,n+1))
              tem2 = cs0*sqrt(e2(i,n))*sqrt(shr2(i,n))
              e2(i,n+1) = e2(i,n) + (tem1 - tem2) * dz
              zlup = zlup + dz
              if(e2(i,n+1) < 0.0_pbl_wp) then
                ptem = e2(i,n+1) / (e2(i,n+1) - e2(i,n))
                zlup = zlup - ptem * dz
                zlup = max(zlup, 0.0_pbl_wp)
                mlenflg = .false.
              endif
            endif
          enddo
          zldn = 0.0_pbl_wp
          mlenflg = .true.
          do n = k, 1, -1
            if(mlenflg) then
              if(n == 1) then
                dz = zl(i,1)
                tem = real(tsea(i), kind=pbl_wp) * (1.0_pbl_wp
     &              + real(fv, kind=pbl_wp) * max(real(q1(i,1,1),
     &                kind=pbl_wp), qmin))
                tem1 = 2.0_pbl_wp*gotvx(i,n)*(tem-thvx(i,k))
                tem2 = ustar(i)*phims(i)/(vk*dz)
                tem2 = cs0*sqrt(e2(i,n))*tem2
                e2(i,n-1) = e2(i,n) + (tem1 - tem2) * dz
              else
                dz = zl(i,n) - zl(i,n-1)
                tem1 = 2.0_pbl_wp*gotvx(i,n-1)*(thvx(i,n-1)-thvx(i,k))
                tem2 = cs0*sqrt(e2(i,n))*sqrt(shr2(i,n-1))
                e2(i,n-1) = e2(i,n) + (tem1 - tem2) * dz
              endif
              zldn = zldn + dz
              if(e2(i,n-1) < 0.0_pbl_wp) then
                ptem = e2(i,n-1) / (e2(i,n-1) - e2(i,n))
                zldn = zldn - ptem * dz
                zldn = max(zldn, 0.0_pbl_wp)
                mlenflg = .false.
              endif
            endif
          enddo
!
          tem = 0.5_pbl_wp * (zi(i,k+1)-zi(i,k))
          tem1 = min(tem, rlmnz(i,k))
!>  - Following Bougeault and Lacarrere(1989), the characteristic length
!! scale (\f$l_2\f$) (eqn 10 in Han et al.(2019) \cite Han_2019) is given by:
!!\f[
!!   l_2=min(l_{up},l_{down})
!!\f]
!! and dissipation length scale \f$l_d\f$ is given by:
!!\f[
!!   l_d=(l_{up}l_{down})^{1/2}
!!\f]
!!    where \f$l_{up}\f$ and \f$l_{down}\f$ are the distances that a parcel
!! having an initial TKE can travel upward and downward before being stopped
!! by buoyancy effects.
!
! Following Rodier et. al (2017), environmental wind shear effect on
!   mixing length was included.
!
          ptem2 = min(zlup,zldn)
          rlam(i,k) = elmfac * ptem2
          rlam(i,k) = max(rlam(i,k), tem1)
          rlam(i,k) = min(rlam(i,k), real(rlmx, kind=pbl_wp))
!
          ptem2 = sqrt(zlup*zldn)
          ele(i,k) = elefac * ptem2
          ele(i,k) = max(ele(i,k), tem1)
          elmh(i,k)= elmhfac * ele(i,k)
          ele(i,k) = min(ele(i,k), real(elmx, kind=pbl_wp))
          elmh(i,k)= min(elmh(i,k), elmhmx)
!
        enddo
      enddo
!>  - Compute the surface layer length scale (\f$l_1\f$) following
!! Nakanishi (2001) \cite Nakanish_2001 (eqn 9 of Han et al.(2019) \cite Han_2019)
      do k = 1, km1
        do i = 1, im
          tem = vk * zl(i,k)
          if (zol(i) < 0.0_pbl_wp) then
            ptem = 1.0_pbl_wp - 100.0_pbl_wp * zol(i)
            ptem1 = ptem**0.2_pbl_wp
            zk = tem * ptem1
          elseif (zol(i) >= 1.0_pbl_wp) then
            zk = tem / 3.7_pbl_wp
          else
            ptem = 1.0_pbl_wp + 2.7_pbl_wp * zol(i)
            zk = tem / ptem
          endif

          if (tc_pbl == 0) then
            elm(i,k) = zk*rlam(i,k)/(rlam(i,k)+zk)
!>  - If sfc_rlm=1, use zk for elm within surface layer
            if ( sfc_rlm == 1 ) then
              if ( sfcflg(i) .and.
     &             zl(i,k) < min(100.0_pbl_wp,real(hpbl(i), kind=pbl_wp)
     &            *0.05_pbl_wp) ) elm(i,k)=zk
            endif
          else if (tc_pbl == 1) then
            ! new blending method for mixing length
            elm(i,k) = sqrt( 1.0_pbl_wp / ( 1.0_pbl_wp / (zk**2)
     &               + 1.0_pbl_wp / (rlam(i,k)**2) ) )
          endif

!
          if(k == 1) elm(i,k)=zk
!
          dz = zi(i,k+1) - zi(i,k)
          tem = max(gdx(i),dz)
          elm(i,k) = min(elm(i,k), tem)

          if (tc_pbl == 0) then
            ele(i,k) = min(ele(i,k), tem)
          else if (tc_pbl == 1) then
            ele(i,k) = elm(i,k)
          endif
!
        enddo
      enddo
      do i = 1, im
        elm(i,km) = elm(i,km1)
        ele(i,km) = ele(i,km1)
        elmh(i,km)= elmh(i,km1)
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> ## Compute eddy diffusivities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      if(tte_edmf) then
!
      do k = 1, km1
        do i = 1, im
           ptem = sfcfrac*real(hpbl(i), kind=pbl_wp)
           if (zi(i,k+1) <= ptem) then
             ri = ris(i)
           else
             ri = max(bf(i,k)/shr2(i,k),rimin)
           endif
           if(ri < 0.0_pbl_wp) then
             tem = 2.0_pbl_wp * ri - pr0
             epotte = ri / tem
           else
             tem = pr0 + 3.0_pbl_wp * ri
             epotte = ri / tem
           endif
           tkeh(i,k) = tteh(i,k) * (1.0_pbl_wp - epotte)
           tesq(i,k) = tkeh(i,k) / sqrt(tteh(i,k))
        enddo
      enddo
!
      else
!
      do k = 1, km1
        do i = 1, im
           tesq(i,k) = sqrt(tkeh(i,k))
        enddo
      enddo
!
      endif
!
      do k = 1, km1
        do i = 1, im
           tem = 0.5_pbl_wp * (elm(i,k) + elm(i,k+1))
           tem = tem * tesq(i,k)
           ri = max(bf(i,k)/shr2(i,k),rimin)
           if(k < kpbl(i)) then
             if(pcnvflg(i)) then
               dku(i,k) = ckz(i,k) * tem
               dkt(i,k) = dku(i,k) / prn(i,k)
             else
               if(ri < 0.0_pbl_wp) then ! unstable regime
                 dku(i,k) = ckz(i,k) * tem
                 dkt(i,k) = dku(i,k) / prn(i,k)
               else             ! stable regime
                 dkt(i,k) = chz(i,k) * tem
                 dku(i,k) = dkt(i,k) * prn(i,k)
               endif
             endif
           else
             if(ri < 0.0_pbl_wp) then ! unstable regime
               dku(i,k) = ck1 * tem
               dkt(i,k) = rchck * dku(i,k)
             else             ! stable regime
               dkt(i,k) = ch1 * tem
               prnum = pr0 + 2.1_pbl_wp * ri
               prnum = min(prnum,prmax)
               dku(i,k) = dkt(i,k) * prnum
             endif
           endif
!
           if(scuflg(i)) then
             if(k >= mrad(i) .and. k < krad(i)) then
               if(tte_edmf) then
                 tem1 = ck0 * tem
               else
                 tem1 = ckz(i,k) * tem
               endif
               ptem1 = tem1 / prscu
               dku(i,k) = max(dku(i,k), tem1)
               dkt(i,k) = max(dkt(i,k), ptem1)
             endif
           endif
!
           dkq(i,k) = prtke * dkt(i,k)
!
           dkt(i,k) = min(dkt(i,k),dkmax)
           dkt(i,k) = max(dkt(i,k),xkzo(i,k))
           dkq(i,k) = min(dkq(i,k),dkmax)
           dkq(i,k) = max(dkq(i,k),xkzo(i,k))
           dku(i,k) = min(dku(i,k),dkmax)
           dku(i,k) = max(dku(i,k),xkzmo(i,k))
!
        enddo
      enddo
!
!The following is for SA-3D-TKE
      if(sa3dtke) then
! 1. compute LES component of km, kh, and kq (Deardorff 1980)
! calculate thetal
        do k=1,km
          do i=1,im
          pix(i,k)   = real(psk(i), kind=pbl_wp)
     &               / real(prslk(i,k), kind=pbl_wp)
          theta(i,k) = t1_loc(i,k) * pix(i,k)
          tem=theta(i,k)/t1_loc(i,k)
          if(ntiw > 0) then
            tem1=max(q1_loc(i,k,ntcw),qlmin)+
     &           max(q1_loc(i,k,ntiw)+q1_loc(i,k,5)+q1_loc(i,k,6),qlmin)
            thetal(i,k)=theta(i,k)-(real(hvap, kind=pbl_wp)+real(hfus,
     &                  kind=pbl_wp))/real(cp, kind=pbl_wp)
     &                 *tem*tem1
          else
            tem1=max(q1_loc(i,k,ntcw),qlmin)
            thetal(i,k)=theta(i,k)-real(hvap, kind=pbl_wp)/real(cp,
     &                  kind=pbl_wp)*tem*tem1
          endif
          enddo
        enddo

        do k=1,km
          do i=1,im
            dku_les(i,k)  = 0.0_pbl_wp
            dkt_les(i,k)  = 0.0_pbl_wp
            dkq_les(i,k)  = 0.0_pbl_wp
          enddo
        enddo
!
!  eddy diffusivities at model interface (zm level) in LES scale
!
        do k = 1, km1
          do i = 1, im
            tem=gotvx(i,k)*(thetal(i,k+1)-thetal(i,k))*rdzt(i,k)
            dz = zl(i,k+1) - zl(i,k)
            tem1=(real(garea(i), kind=pbl_wp)*dz)**h1
! calculate LES mixing length
            if(tem > 0.0_pbl_wp) then
              elm_les=0.76_pbl_wp*sqrt(tke(i,k))/sqrt(tem)
              elm_les=min(elm_les,tem1)
            else
              elm_les=tem1
            endif
! calculate km, kh, and kq for LES
            dku_les(i,k)=0.1_pbl_wp*elm_les*sqrt(tkeh(i,k))
            dkt_les(i,k)=(1.0_pbl_wp+2.0_pbl_wp*elm_les/tem1)
     &                  *dku_les(i,k)
            dkq_les(i,k)=dkt_les(i,k)
            dku_les(i,k) = min(dku_les(i,k),dkmaxles)
            dkt_les(i,k) = min(dkt_les(i,k),dkmaxles)
            dkq_les(i,k) = min(dkq_les(i,k),dkmaxles)
          enddo
        enddo
!
! calculate blending coefficients for km, kt, kq, and nonlocal mixing
! finding scale of large eddies from TKE
        do i=1,im
          tkemax(i) = tke(i,1)
          ktkemax(i) = 1
        enddo
        do k = 2, kmpbl
          do i = 1, im
            if(tke(i,k) > tkemax(i)) then
              tkemax(i) = tke(i,k)
              ktkemax(i) = k
            endif
          enddo
        enddo
        do i=1,im
          flg(i) = .true.
          scl(i) = 0.0_pbl_wp
          if(zl(i,ktkemax(i)) > sclmax) then
            flg(i) = .false.
            scl(i) = sclmin
          endif
        enddo
        do k = 1, kmpbl
          do i = 1, im
            if(flg(i) .and. k > ktkemax(i)) then
              scl(i) = zl(i,k)
              tem = 0.5_pbl_wp*tkemax(i)
              if(tke(i,k) < tem) flg(i) = .false.
            endif
          enddo
        enddo
        do i=1,im
          scl(i)=max(scl(i), sclmin)
          scl(i)=min(scl(i), sclmax)
          scl(i)=max(scl(i), real(hpbl(i), kind=pbl_wp))
          pfdx(i)=gdx(i)/scl(i)
        enddo
!
        do i = 1, im
!            partition function for local fluxes
          pfl(i)=cpl1*(pfdx(i)**2+cpl2*pfdx(i)**0.5_pbl_wp-cpl3)/
     &           (pfdx(i)**2+cpl4*pfdx(i)**0.5_pbl_wp+cpl5)+cpl6
          pfl(i)=min(max(pfl(i),0.0_pbl_wp),1.0_pbl_wp)
!            partition function for TKE
          pftke(i)=(pfdx(i)**2+cptke1*pfdx(i)**(2.0_pbl_wp/3.0_pbl_wp))/
     &             (pfdx(i)**2+cptke2*pfdx(i)**(2.0_pbl_wp/3.0_pbl_wp)
     &             +cptke3)
          pftke(i)=min(max(pftke(i),0.0_pbl_wp),1.0_pbl_wp)
        enddo
!
!     blending LES and MS components of vertical km,kt, and kq
!
        do k = 1,km1
          do i=1,im
            dkq(i,k)=(1.0_pbl_wp-pfl(i))*dkq_les(i,k)+pfl(i)*dkq(i,k)
            dkt(i,k)=(1.0_pbl_wp-pfl(i))*dkt_les(i,k)+pfl(i)*dkt(i,k)
            dku(i,k)=(1.0_pbl_wp-pfl(i))*dku_les(i,k)+pfl(i)*dku(i,k)
          enddo
        enddo
!
! 2. compute MS horizontal km
!
        do k = 1, km
          do i = 1, im
            dku_h(i,k)=ckh*elmh(i,k)*sqrt(tke(i,k))
            dkq_h(i,k)=dku_h(i,k)
          enddo
        enddo
!
!  eddy diffusivities at model layer (zl level) in LES scale
!
        do k = 1, km1
          do i = 1, im
            if(k > 1) then
              dz = zl(i,k+1) - zl(i,k-1)
              tem=gotvx(i,k)*(thetal(i,k+1)-thetal(i,k-1))/dz
            else
              dz = zl(i,k+1)
              tem=gotvx(i,k)*(thetal(i,k+1)-thvs(i))/dz
            endif
            dz = zi(i,k+1) - zi(i,k)
            tem1=(real(garea(i), kind=pbl_wp)*dz)**h1
! calculate LES mixing length
            if(tem > 0.0_pbl_wp) then
              elm_les=0.76_pbl_wp*sqrt(tke(i,k))/sqrt(tem)
              elm_les=min(elm_les,tem1)
            else
              elm_les=tem1
            endif
            ele_les(i,k)=elm_les
! calculate km, kh, and kq for LES
            dku_les(i,k)=0.1_pbl_wp*elm_les*sqrt(tke(i,k))
            dkq_les(i,k)=(1.0_pbl_wp+2.0_pbl_wp*elm_les/tem1)
     &                  *dku_les(i,k)
            dkq_les(i,k)=(1.0_pbl_wp+2.0_pbl_wp*elm_les/tem1)
     &                  *dku_les(i,k)
            dku_les(i,k) = min(dku_les(i,k),dkmaxles)
            dkq_les(i,k) = min(dkq_les(i,k),dkmaxles)
          enddo
        enddo
!
        do k = 1,km1
          do i=1,im
            dku_h(i,k)=(1.0_pbl_wp-pfl(i))*dku_les(i,k)+pfl(i)
     &                *dku_h(i,k)
            dkq_h(i,k)=(1.0_pbl_wp-pfl(i))*dkq_les(i,k)+pfl(i)
     &                *dkq_h(i,k)
          enddo
        enddo
        do i = 1, im
          dku_h(i,km)=dku_h(i,km1)
          dkq_h(i,km)=dkq_h(i,km1)
        enddo
!
      endif !sa3dtke

!PCC_CANOPY------------------------------------
      if (do_canopy .and. cplaqm) then

        do k = 1, km1-1

          kount=0 !IVAI
          ZOOOX(:) = 1.0_pbl_wp  ! IVAI

          do i = 1, im

            FCH = real(cfch(i), kind=pbl_wp)               !top of canopy from AQM canopy inputs
            IF (k .EQ. 1) THEN !use model layer interfaces
              KCAN = 1
            ELSE
              IF (   FCH .GT. zi(i,k)
     &         .AND. FCH .LE. zi(i,k+1) ) THEN
                KCAN = 1
              ELSE
                KCAN = 0
              END IF
            END IF

            IF (KCAN .EQ. 1) THEN !canopy inside model layer
!           Check for other Contiguous Canopy Grid Cell Conditions

!           Not a contiguous canopy cell
            IF (   real(claie(i), kind=pbl_wp) .LT. 0.1_pbl_wp
     &        .OR. FCH .LT. 0.5_pbl_wp
     &        .OR. MAX(0.0_pbl_wp, 1.0_pbl_wp - real(cfrt(i),
     &             kind=pbl_wp)) .GT. 0.75_pbl_wp
     &        .OR. real(cpopu(i), kind=pbl_wp) .GT. 10000.0_pbl_wp
     &        .OR. (EXP(-0.5_pbl_wp*real(claie(i), kind=pbl_wp)
     &            *real(cclu(i), kind=pbl_wp)) .GT. 0.45_pbl_wp
     &            .AND. FCH  .LT. 18.0_pbl_wp) ) THEN

            ELSE ! There is a contiguous forest canopy, apply correction over canopy layers

!Raupauch M. R. A Practical Lagrangian method for relating scalar
!concentrations to
! source distributions in vegetation canopies. Q. J. R. Meteor. Soc.
! (1989), 115, pp 609-632
              MOL  = zol(i)/zl(i,k)  !Monin-Obukhov Length in layer
              HOL = FCH/MOL  !local canopy stability parameter (hc/MOL)
              ZCAN = zi(i,k+1) ! Initialize each model layer top that contains canopy (m)
                             ! Integrate across total model interface
              ZFL = ZCAN     ! Set ZFL = ZCAN
              COUNTCAN = 0   ! Initialize canopy layers

              IF (k .EQ. 1) THEN  !Find bottom in each model layer
                  BOTCAN = 0.0_pbl_wp    ! 0.5 (Jan10)
              ELSE
                  BOTCAN = zi(i,k)
              END IF

              DO WHILE (ZCAN.GE.BOTCAN)
                 ! TLCAN = Lagrangian timescale
                 TLCAN = (FCH/ustar(i)) * (
     &                   (0.256_pbl_wp * (ZCAN-(0.75_pbl_wp*FCH))
     &                 / FCH ) + (0.492_pbl_wp * EXP((-0.256_pbl_wp
     &                 * ZCAN / FCH) /0.492_pbl_wp)) )
                 IF ( HOL .LT. -0.1_pbl_wp )  THEN  !STRONG UNSTABLE
                  IF ( ZCAN/FCH .GT. 1.25_pbl_wp ) THEN !SIGMACAN = Eulerian vertical velocity variance
                    SIGMACAN = 1.25_pbl_wp*ustar(i)
                  END IF
                  IF ( ZCAN/FCH .GE. 0.175_pbl_wp
     &               .AND. ZCAN/FCH .LE. 1.25_pbl_wp ) THEN
                      SIGMACAN = ustar(i) * ( 0.75_pbl_wp + (0.5_pbl_wp
     &                         * COS((real(PI, kind=pbl_wp)
     &                         / 1.06818_pbl_wp) * (1.25_pbl_wp - (ZCAN
     &                         / FCH)))) )
                  END IF
                  IF ( ZCAN/FCH .LT. 0.175_pbl_wp )  THEN
                    SIGMACAN = 0.25_pbl_wp*ustar(i)
                  END IF
                 END IF
                 IF ( HOL .GE. -0.1_pbl_wp  .AND. HOL .LT. 0.1_pbl_wp)
     &              THEN !WEAKLY UNSTABLE to NEUTRAL
                  IF ( ZCAN/FCH .GT. 1.25_pbl_wp ) THEN
                    SIGMACAN = 1.0_pbl_wp*ustar(i)
                  END IF
                  IF ( ZCAN/FCH .GE. 0.175_pbl_wp
     &               .AND. ZCAN/FCH .LE. 1.25_pbl_wp ) THEN
                    SIGMACAN = ustar(i) * ( 0.625_pbl_wp +
     &                       (0.375_pbl_wp* COS((real(PI, kind=pbl_wp)
     &                       / 1.06818_pbl_wp) * (1.25_pbl_wp
     &                       - (ZCAN/FCH)))) )
                  END IF
                  IF ( ZCAN/FCH .LT. 0.175_pbl_wp )  THEN
                    SIGMACAN = 0.25_pbl_wp*ustar(i)
                  END IF
                 END IF
                 IF ( HOL .GE.  0.1_pbl_wp  .AND.  HOL .LT. 0.9_pbl_wp )
     &              THEN !STABLE
                  IF ( ZCAN/FCH .GT. 1.25_pbl_wp ) THEN
                    SIGMACAN = 0.25_pbl_wp*(4.375_pbl_wp - (3.75_pbl_wp
     &                       * HOL)) * ustar(i)
                  END IF
                  IF ( ZCAN/FCH .GE. 0.175_pbl_wp
     &                 .AND. ZCAN/FCH .LE. 1.25_pbl_wp ) THEN
                    RRCAN=4.375_pbl_wp-(3.75_pbl_wp*HOL)
                    AACAN=(0.125_pbl_wp*RRCAN) + 0.125_pbl_wp
                    BBCAN=(0.125_pbl_wp*RRCAN) - 0.125_pbl_wp
                    SIGMACAN = ustar(i) * ( AACAN + (BBCAN
     &                       *COS((real(PI, kind=pbl_wp)/1.06818_pbl_wp)
     &                       *(1.25_pbl_wp - (ZCAN/FCH)))) )
                  END IF
                  IF ( ZCAN/FCH .LT. 0.175_pbl_wp )  THEN
                    SIGMACAN = 0.25_pbl_wp*ustar(i)
                  END IF
                 END IF
                 IF ( HOL .GE.  0.9_pbl_wp ) THEN  !VERY STABLE
                    SIGMACAN = 0.25_pbl_wp*ustar(i)
                 END IF

                 IF ( ZCAN .EQ. ZFL )  ! THEN ! Each model layer that includes canopy
     &              EDDYVEST1 = (SIGMACAN*SIGMACAN)*TLCAN ! Model-layer top diffusivity

! Average In-Canopy diffusivity gives larger canopy correction (2022)
!                IF ( ZCAN .LE. FCH )  THEN          ! Average In-Canopy   diffusivity (2022) 
! Average model-layer diffusivity gives smaller canopy correction (pre-2022) & Jan9, 2026
                 IF ( ZCAN .LE. ZFL )  THEN         ! Average model-layer diffusivity (2026)

                   COUNTCAN = COUNTCAN + 1
                   ZCANX     (COUNTCAN) = ZCAN
                   EDDYVESTX (COUNTCAN) = (SIGMACAN*SIGMACAN)*TLCAN

                 END IF

                 ZCAN = ZCAN-0.5_pbl_wp  !step down in-canopy resolution of 0.5m

              END DO !end loop on canopy layers

            IF (COUNTCAN > 0 ) THEN

              IF (COUNTCAN .EQ. 1) THEN
                 ZZ_INT= ZCANX(COUNTCAN)  ! FCH (set to FCH if COUNTCAN is 1)
              ELSE
                 ZZ_INT= IntegrateTrapezoid(
     &                  ZCANX(COUNTCAN:1:-1), ZOOOX(COUNTCAN:1:-1))
              END IF
 
              IF (COUNTCAN .EQ. 1) THEN
                 EDDYVEST_INT = EDDYVESTX(COUNTCAN)
              ELSE

                 EDDYVEST_INT = IntegrateTrapezoid((ZCANX(COUNTCAN:1:-1)
     &                                ),EDDYVESTX(COUNTCAN:1:-1)) / 
     &                                  ZZ_INT
              END IF

! Comment out to turn OFF the integrated canopy effect
!Apr17        dkt(i,k)= (dkt(i,k)/EDDYVEST1) * EDDYVEST_INT !Scale dkt to resolved eddy diffusivity
!Apr17        dkq(i,k)= (dkq(i,k)/EDDYVEST1) * EDDYVEST_INT !Scale dkq to resolved eddy diffusivity
!Apr17        dku(i,k)= (dku(i,k)/EDDYVEST1) * EDDYVEST_INT !Scale dku to resolved eddy diffusivity

            END IF ! (COUNTCAN > 0)

            END IF   ! contiguous canopy conditions

            END IF   ! (KCAN .EQ. 1) model layer(s) containing canopy

          enddo !i

          kount = kount + 1 

        enddo !k

      endif !do_canopy .and. cplaqm

!> ## Compute TKE.
!!  - Compute a minimum TKE deduced from background diffusivity for momentum.
!
      do k = 1, km1
        do i = 1, im
          if(k == 1) then
            tem = ckz(i,1)
            tem1 = 0.5_pbl_wp * xkzmo(i,1)
          else
            tem = 0.5_pbl_wp * (ckz(i,k-1) + ckz(i,k))
            tem1 = 0.5_pbl_wp * (xkzmo(i,k-1) + xkzmo(i,k))
          endif
          ptem = tem1 / (tem * elm(i,k))
          tkmnz(i,k) = ptem * ptem
          tkmnz(i,k) = min(tkmnz(i,k), tkbmx)
          tkmnz(i,k) = max(tkmnz(i,k), tkmin)
        enddo
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  - Compute buoyancy and shear productions of TKE or TTE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      do k = 1, km1
        do i = 1, im
          if (k == 1) then
            tem = -dkt(i,1) * bf(i,1)
!           if(pcnvflg(i)) then
!             ptem1 = xmf(i,1) * buou(i,1)
!           else
              ptem1 = 0.0_pbl_wp
!           endif
            if(scuflg(i) .and. mrad(i) == 1) then
              ptem2 = xmfd(i,1) * buod(i,1)
            else
              ptem2 = 0.0_pbl_wp
            endif
            tem = tem + ptem1 + ptem2
            buop = 0.5_pbl_wp * (gotvx(i,1) * sflux(i) + tem)
!
            if(sa3dtke) then
              tem = 2.0_pbl_wp * dku_h(i,1)
              tem1 = dku(i,1)*def_1(i,1)+tem*def_2(i,1)
            else
              tem1 = dku(i,1) * shr2(i,1)
            endif
!
            tem = (u1_loc(i,2)-u1_loc(i,1))*rdzt(i,1)
!           if(pcnvflg(i)) then
!             ptem = xmf(i,1) * tem
!             ptem1 = 0.5 * ptem * (u1(i,2)-ucko(i,2))
!           else
              ptem1 = 0.0_pbl_wp
!           endif
            if(scuflg(i) .and. mrad(i) == 1) then
              ptem = ucdo(i,1)+ucdo(i,2)-u1_loc(i,1)-u1_loc(i,2)
              ptem = 0.5_pbl_wp * tem * xmfd(i,1) * ptem      
              ptem = 0.5_pbl_wp * tem * xmfd(i,1) * ptem
            else
              ptem = 0.0_pbl_wp
            endif
            ptem1 = ptem1 + ptem
!
            tem = (v1_loc(i,2)-v1_loc(i,1))*rdzt(i,1)
!           if(pcnvflg(i)) then
!             ptem = xmf(i,1) * tem
!             ptem2 = 0.5 * ptem * (v1(i,2)-vcko(i,2))
!           else
              ptem2 = 0.0_pbl_wp
!           endif
            if(scuflg(i) .and. mrad(i) == 1) then
              ptem = vcdo(i,1)+vcdo(i,2)-v1_loc(i,1)-v1_loc(i,2)      
              ptem = 0.5_pbl_wp * tem * xmfd(i,1) * ptem
            else
              ptem = 0.0_pbl_wp
            endif
            ptem2 = ptem2 + ptem

            tem2 = real(stress(i), kind=pbl_wp)*ustar(i)*phims(i)
     &            /(vk*zl(i,1))
            shrp = 0.5_pbl_wp* (tem1 + ptem1 + ptem2 + tem2)
          else
            tem1 = -dkt(i,k-1) * bf(i,k-1)
            tem2 = -dkt(i,k) * bf(i,k)
            tem  = 0.5_pbl_wp * (tem1 + tem2)
            if(pcnvflg(i) .and. k <= kpbl(i)) then
              ptem = 0.5_pbl_wp * (xmf(i,k-1) + xmf(i,k))
              ptem1 = ptem * buou(i,k)
            else
              ptem1 = 0.0_pbl_wp
            endif
            if(scuflg(i)) then
              if(k >= mrad(i) .and. k < krad(i)) then
                ptem0 = 0.5_pbl_wp * (xmfd(i,k-1) + xmfd(i,k))
                ptem2 = ptem0 * buod(i,k)
              else
                ptem2 = 0.0_pbl_wp
              endif
            else
              ptem2 = 0.0_pbl_wp
            endif
            buop = tem + ptem1 + ptem2
!
            if(sa3dtke) then
!         obtaining 3d shear production from dycore
              tem2 = 2.0_pbl_wp*dku_h(i,k)
              tem1 = dku(i,k-1)*def_1(i,k-1)
              tem2 = dku(i,k)*def_1(i,k)+tem2*def_2(i,k)
            else
              tem1 = dku(i,k-1) * shr2(i,k-1)
              tem2 = dku(i,k) * shr2(i,k)
            endif
            tem  = 0.5_pbl_wp * (tem1 + tem2)
            tem1 = (u1_loc(i,k+1)-u1_loc(i,k))*rdzt(i,k)
            tem2 = (u1_loc(i,k)-u1_loc(i,k-1))*rdzt(i,k-1)
            if(pcnvflg(i) .and. k <= kpbl(i)) then
              ptem = xmf(i,k) * tem1 + xmf(i,k-1) * tem2
              ptem1 = 0.5_pbl_wp * ptem * (u1_loc(i,k)-ucko(i,k))
            else
              ptem1 = 0.0_pbl_wp
            endif
            if(scuflg(i)) then
              if(k >= mrad(i) .and. k < krad(i)) then
                ptem0 = xmfd(i,k) * tem1 + xmfd(i,k-1) * tem2
                ptem2 = 0.5_pbl_wp * ptem0 * (ucdo(i,k)-u1_loc(i,k))
              else
                ptem2 = 0.0_pbl_wp
              endif
            else
              ptem2 = 0.0_pbl_wp
            endif
            shrp = tem + ptem1 + ptem2
            tem1 = (v1_loc(i,k+1)-v1_loc(i,k))*rdzt(i,k)
            tem2 = (v1_loc(i,k)-v1_loc(i,k-1))*rdzt(i,k-1)
            if(pcnvflg(i) .and. k <= kpbl(i)) then
              ptem = xmf(i,k) * tem1 + xmf(i,k-1) * tem2
              ptem1 = 0.5_pbl_wp * ptem * (v1_loc(i,k)-vcko(i,k))
            else
              ptem1 = 0.0_pbl_wp
            endif
            if(scuflg(i)) then
              if(k >= mrad(i) .and. k < krad(i)) then
                ptem0 = xmfd(i,k) * tem1 + xmfd(i,k-1) * tem2
                ptem2 = 0.5_pbl_wp * ptem0 * (vcdo(i,k)-v1_loc(i,k))
              else
                ptem2 = 0.0_pbl_wp
              endif
            else
              ptem2 = 0.0_pbl_wp
            endif
            shrp = shrp + ptem1 + ptem2
          endif
          if(tte_edmf) then
            if(buop > 0.0_pbl_wp) then
              prod(i,k) = 2.0_pbl_wp * buop + shrp
            else
              prod(i,k) = shrp
            endif
          else
            prod(i,k) = buop + shrp
          endif
        enddo
      enddo
!
!----------------------------------------------------------------------
!>  - First predict te due to te production & dissipation(diss)
!
      if(sa3dtke) then
!The following is for SA-3D-TKE
       dtn = dt2 / real(ndt, kind=pbl_wp)       
       do n = 1, ndt
       do k = 1,km1
        do i=1,im
           tem = sqrt(te(i,k))
!        calculating 3D TKE transport and pressure correlation
           ptem1 = ce0 / ele(i,k)
           dz = zi(i,k+1) - zi(i,k)
           tem1=(real(garea(i), kind=pbl_wp)*dz)**h1
           tem2=0.19_pbl_wp+0.51_pbl_wp*ele_les(i,k)/tem1
           ptem2= tem2 / ele_les(i,k)
           ptem=(1.0_pbl_wp-pftke(i))*ptem2+pftke(i)*ptem1
           disste = ptem * te(i,k) * tem
           tem1 = prod(i,k) + te(i,k) / dtn
           disste=max(min(disste, tem1), 0.0_pbl_wp)
           if(.not. tte_edmf) diss(i,k) = disste
!          tem=2.0*def_3(i,k)
           tem=def_3(i,k)
!          tem=min(tem,1.0)
           te(i,k) = te(i,k) + dtn * (prod(i,k)-disste+tem)
!          te(i,k) = max(te(i,k), tkmin)
           te(i,k) = max(te(i,k), tkmnz(i,k))
        enddo
       enddo
       enddo
      else
       dtn = dt2 / real(ndt, kind=pbl_wp)
       do n = 1, ndt
       do k = 1,km1
        do i=1,im
           tem = sqrt(te(i,k))
           ptem = ce0 / ele(i,k)
           disste = ptem * te(i,k) * tem
           tem1 = prod(i,k) + te(i,k) / dtn
           disste = max(min(disste, tem1), 0.0_pbl_wp)
           if(.not. tte_edmf) diss(i,k) = disste
           te(i,k) = te(i,k) + dtn * (prod(i,k)-disste)
           te(i,k) = max(te(i,k), tkmnz(i,k))
!          te(i,k) = max(te(i,k), tkmin)
        enddo
       enddo
       enddo
      endif !sa3dtke
!
! TKE dissipation for dissipative heating computation in TTE-EDMF
!
      if(tte_edmf) then
      do k = 1, km1
        do i = 1, im
          tem = sqrt(tke(i,k))
          if(sa3dtke) then
            ptem1 = ce0 / ele(i,k)
            dz = zi(i,k+1) - zi(i,k)
            tem1=(real(garea(i), kind=pbl_wp)*dz)**h1
            tem2=0.19_pbl_wp+0.51_pbl_wp*ele_les(i,k)/tem1
            ptem2= tem2 / ele_les(i,k)
            ptem=(1.0_pbl_wp-pftke(i))*ptem2+pftke(i)*ptem1
            diss(i,k) = ptem * tke(i,k) * tem
          else
            ptem = ce0 / ele(i,k)
            diss(i,k) = ptem * tke(i,k) * tem
          endif
        enddo
      enddo
      endif
!
!>  - Compute updraft & downdraft properties for TKE or TTE
!
      do k = 1, km
        do i = 1, im
          if(pcnvflg(i)) then
            qcko(i,k,ntke) = te(i,k)
          endif
          if(scuflg(i)) then
            qcdo(i,k,ntke) = te(i,k)
          endif
        enddo
      enddo
      do k = 2, kmpbl
        do i = 1, im
          if (pcnvflg(i) .and. k <= kpbl(i)) then
             dz   = zl(i,k) - zl(i,k-1)
             tem  = 0.5_pbl_wp * xlamue(i,k-1) * dz
             factor = 1.0_pbl_wp + tem
             qcko(i,k,ntke)=((1.0_pbl_wp-tem)*qcko(i,k-1,ntke)+tem*
     &                      (te(i,k)+te(i,k-1)))/factor
          endif
        enddo
      enddo
      do k = kmscu, 1, -1
        do i = 1, im
          if (scuflg(i) .and. k < krad(i)) then
            if(k >= mrad(i)) then
              dz = zl(i,k+1) - zl(i,k)
              tem  = 0.5_pbl_wp * xlamde(i,k) * dz
              factor = 1.0_pbl_wp + tem
              qcdo(i,k,ntke)=((1.0_pbl_wp-tem)*qcdo(i,k+1,ntke)+tem*
     &                       (te(i,k)+te(i,k+1)))/factor
            endif
          endif
        enddo
      enddo
!
!--------------------------------------------------------
!  compute variables for TVD flux-limiter scheme
!     on environmental subsidence and uplifting
!
      kps = max(kmpbl, kmscu)
!
! for moisture and tracers including hydrometeors
!
      do n=1,ntrac1
        do k=1,kps
          do i=1,im
            qh(i,k,n) = 0.5_pbl_wp * (q1_loc(i,k,n)+q1_loc(i,k+1,n))
          enddo
        enddo
        do k=1,kps
          do i=1,im
            q_diff(i,k,n) = q1_loc(i,k,n) - q1_loc(i,k+1,n)
          enddo
        enddo
        do i=1,im
          if(q1_loc(i,1,n) >= 0.0_pbl_wp) then
            q_diff(i,0,n) = max(0.0_pbl_wp,2.0_pbl_wp*q1_loc(i,1,n)
     &                    -q1_loc(i,2,n))-q1_loc(i,1,n)
          else
            q_diff(i,0,n) = min(0.0_pbl_wp,2.0_pbl_wp*q1_loc(i,1,n)
     &                    -q1_loc(i,2,n))-q1_loc(i,1,n)
          endif
        enddo
      enddo
!
      do n = 1, ntrac1
!
        do k = 1, kps
          do i = 1, im
            kmx = max(kpbl(i), krad(i))
            q_half(i,k,n) = qh(i,k,n)
            if((pcnvflg(i) .or. scuflg(i)) .and. k < kmx) then
              tem = 0.0_pbl_wp
              if(pcnvflg(i) .and. k < kpbl(i)) then
                tem = xmf(i,k)
              endif
              if(scuflg(i) .and.
     &            (k >= mrad(i) .and. k < krad(i))) then
                tem = tem - xmfd(i,k)
              endif
              if(tem > 0.0_pbl_wp) then
                rrkp = 0.0_pbl_wp
                if(abs(q_diff(i,k,n)) > 1.e-22_pbl_wp)
     &                 rrkp = q_diff(i,k+1,n) / q_diff(i,k,n)
                phkp = (rrkp+abs(rrkp)) / (1.0_pbl_wp+abs(rrkp))
                q_half(i,k,n) = q1_loc(i,k+1,n) +
     &                    phkp*(qh(i,k,n)-q1_loc(i,k+1,n))
              elseif (tem < 0.0_pbl_wp) then
                rrkp = 0.0_pbl_wp
                if(abs(q_diff(i,k,n)) > 1.e-22_pbl_wp)
     &                 rrkp = q_diff(i,k-1,n) / q_diff(i,k,n)
                phkp = (rrkp+abs(rrkp)) / (1.0_pbl_wp+abs(rrkp))
                q_half(i,k,n) = q1_loc(i,k,n) +
     &                     phkp*(qh(i,k,n)-q1_loc(i,k,n))
              endif
            endif
          enddo
        enddo
!
      enddo
!
! for TKE or TTE
!
      do k=1,kps
        do i=1,im
          tei(i,k) = 0.5_pbl_wp * (te(i,k)+te(i,k+1))
        enddo
      enddo

      do k=1,kps
        do i=1,im
          e_diff(i,k) = te(i,k) - te(i,k+1)
        enddo
      enddo
      do i=1,im
        if(te(i,1) >= 0.0_pbl_wp) then
          e_diff(i,0) = max(0.0_pbl_wp,2.0_pbl_wp*te(i,1)-te(i,2))-
     &                        te(i,1)
        else
          e_diff(i,0) = min(0.0_pbl_wp,2.0_pbl_wp*te(i,1)-te(i,2))-
     &                          te(i,1)
        endif
      enddo
!
      do k = 1, kps
        do i = 1, im
          kmx = max(kpbl(i), krad(i))
          e_half(i,k) = tei(i,k)
          if((pcnvflg(i) .or. scuflg(i)) .and. (k < kmx)) then
            tem = 0.0_pbl_wp
            if(pcnvflg(i) .and. k < kpbl(i)) then
              tem = xmf(i,k)
            endif
            if(scuflg(i) .and.
     &          (k >= mrad(i) .and. k < krad(i))) then
              tem = tem - xmfd(i,k)
            endif
            if(tem > 0.0_pbl_wp) then
              rrkp = 0.0_pbl_wp
              if(abs(e_diff(i,k)) > 1.e-22_pbl_wp)
     &               rrkp = e_diff(i,k+1) / e_diff(i,k)
              phkp = (rrkp+abs(rrkp)) / (1.0_pbl_wp+abs(rrkp))
              e_half(i,k) = te(i,k+1) +
     &                   phkp*(tei(i,k)-te(i,k+1))
            elseif (tem < 0.0_pbl_wp) then
              rrkp = 0.0_pbl_wp
              if(abs(e_diff(i,k)) > 1.e-22_pbl_wp)
     &               rrkp = e_diff(i,k-1) / e_diff(i,k)
              phkp = (rrkp+abs(rrkp)) / (1.0_pbl_wp+abs(rrkp))
              e_half(i,k) = te(i,k) +
     &                  phkp*(tei(i,k)-te(i,k))
            endif
          endif
        enddo
      enddo
!
!----------------------------------------------------------------------
!>  - Compute tridiagonal matrix elements for TKE or TTE
!
      do i=1,im
         ad(i,1) = 1.0_pbl_wp
         f1(i,1) = te(i,1)
      enddo
!
      do k = 1,km1
        do i=1,im
          dtodsd  = dt2/real(del(i,k), kind=pbl_wp)
          dtodsu  = dt2/real(del(i,k+1), kind=pbl_wp)
          dsig    = real(prsl(i,k), kind=pbl_wp)-real(prsl(i,k+1),
     &              kind=pbl_wp)
          rdz     = rdzt(i,k)
          tem1    = dsig * dkq(i,k) * rdz
          dsdz2   = tem1 * rdz
          au(i,k) = -dtodsd*dsdz2
          al(i,k) = -dtodsu*dsdz2
          ad(i,k) = ad(i,k)-au(i,k)
          ad(i,k+1)= 1.0_pbl_wp-al(i,k)
          tem2    = dsig * rdz
!
          if(pcnvflg(i) .and. k < kpbl(i)) then
             ptem      = 0.5_pbl_wp * tem2 * xmf(i,k)
             ptem1     = dtodsd * ptem
             ptem2     = dtodsu * ptem
             ptem      = qcko(i,k,ntke) + qcko(i,k+1,ntke)
             f1(i,k)   = f1(i,k) - ptem * ptem1
             f1(i,k+1) = te(i,k+1) + ptem * ptem2
          else
             f1(i,k+1) = te(i,k+1)
          endif
!
          if(scuflg(i)) then
            if(k >= mrad(i) .and. k < krad(i)) then
              ptem      = 0.5_pbl_wp * tem2 * xmfd(i,k)
              ptem1     = dtodsd * ptem
              ptem2     = dtodsu * ptem
              ptem      = qcdo(i,k,ntke) + qcdo(i,k+1,ntke)
              f1(i,k)   = f1(i,k) + ptem * ptem1
              f1(i,k+1) = f1(i,k+1) - ptem * ptem2
            endif
          endif
!
          kmx = max(kpbl(i), krad(i))
          if((pcnvflg(i) .or. scuflg(i)) .and. (k < kmx)) then
            ptem  = tem2 * (xmf(i,k) - xmfd(i,k))
            ptem1 = dtodsd * ptem
            ptem2 = dtodsu * ptem
            f1(i,k)   = f1(i,k) + e_half(i,k) * ptem1
            f1(i,k+1) = f1(i,k+1) - e_half(i,k) * ptem2
          endif
!
        enddo
      enddo
c
!>  - Call tridit() to solve tridiagonal problem for TKE
c
      call tridit_satmedmf(im,km,1,al,ad,au,f1,au_loc,f1_loc)
!
! Negative TKE or TTE are set to zero after borrowing it from positive
!   values within the mass-flux transport layers
!
      do i = 1,im
        tsumn(i) = 0.0_pbl_wp
        tsump(i) = 0.0_pbl_wp
        rtnp(i) = 1.0_pbl_wp
      enddo
      do k = 1,kps
        do i = 1,im
          if(pcnvflg(i) .and. scuflg(i)) then
            kbx = 1
            kmx = max(kpbl(i), krad(i))
          elseif(pcnvflg(i) .and. .not. scuflg(i)) then
            kbx = 1
            kmx = kpbl(i)
          elseif(.not. pcnvflg(i) .and. scuflg(i)) then
            kbx = mrad(i)
            kmx = krad(i)
          endif
          if((pcnvflg(i) .or. scuflg(i)) .and.
     &       (k >= kbx .and. k <= kmx)) then
            tem = f1(i,k) * real(del(i,k), kind=pbl_wp) * real(gravi,
     &            kind=pbl_wp) 
            if(f1(i,k) < 0.0_pbl_wp) tsumn(i) = tsumn(i) + tem
            if(f1(i,k) > 0.0_pbl_wp) tsump(i) = tsump(i) + tem
          endif
        enddo
      enddo
      do i = 1,im
        if(pcnvflg(i) .or. scuflg(i)) then
          if(tsump(i) > 0.0_pbl_wp .and. tsumn(i) < 0.0_pbl_wp) then
            if(tsump(i) > abs(tsumn(i))) then
              rtnp(i) = tsumn(i) / tsump(i)
            else
              rtnp(i) = tsump(i) / tsumn(i)
            endif
          endif
        endif
      enddo
      do k = 1,kps
        do i = 1,im
          if(pcnvflg(i) .and. scuflg(i)) then
            kbx = 1
            kmx = max(kpbl(i), krad(i))
          elseif(pcnvflg(i) .and. .not. scuflg(i)) then
            kbx = 1
            kmx = kpbl(i)
          elseif(.not. pcnvflg(i) .and. scuflg(i)) then
            kbx = mrad(i)
            kmx = krad(i)
          endif
          if((pcnvflg(i) .or. scuflg(i)) .and.
     &       (k >= kbx .and. k <= kmx)) then
            if(rtnp(i) < 0.0_pbl_wp) then
              if(tsump(i) > abs(tsumn(i))) then
                if(f1(i,k) < 0.0_pbl_wp) f1(i,k) = 0.0_pbl_wp
                if(f1(i,k) > 0.0_pbl_wp) f1(i,k) = (1.0_pbl_wp+rtnp(i))
     &                                            *f1(i,k)
              else
                if(f1(i,k) < 0.0_pbl_wp) f1(i,k) = (1.0_pbl_wp+rtnp(i))
     &                                            *f1(i,k)
                if(f1(i,k) > 0.0_pbl_wp) f1(i,k) = 0.0_pbl_wp
              endif
            endif
          endif
        enddo
      enddo
!
! To remove negative TKEs or TTEs which were leaked out of the mass-flux transport layers
!   by eddy diffusion or potential negative TKEs or TTEs from the diffusion scheme,
!   positive TKEs or TTEs are borrowed again now from the entire layers
!
      do i = 1,im
        tsumn(i) = 0.0_pbl_wp
        tsump(i) = 0.0_pbl_wp
        rtnp(i) = 1.0_pbl_wp
      enddo
      do k = 1,km
        do i = 1,im
          tem = f1(i,k) * real(del(i,k), kind=pbl_wp) * real(gravi,
     &          kind=pbl_wp)
          if(f1(i,k) < 0.0_pbl_wp) tsumn(i) = tsumn(i) + tem
          if(f1(i,k) > 0.0_pbl_wp) tsump(i) = tsump(i) + tem
        enddo
      enddo
      do i = 1,im
        if(tsump(i) > 0.0_pbl_wp .and. tsumn(i) < 0.0_pbl_wp) then
          if(tsump(i) > abs(tsumn(i))) then
            rtnp(i) = tsumn(i) / tsump(i)
          else
            rtnp(i) = tsump(i) / tsumn(i)
          endif
        endif
      enddo
      do k = 1,km
        do i = 1,im
          if(rtnp(i) < 0.0_pbl_wp) then
            if(tsump(i) > abs(tsumn(i))) then
              if(f1(i,k) < 0.0_pbl_wp) f1(i,k) = 0.0_pbl_wp
              if(f1(i,k) > 0.0_pbl_wp) f1(i,k) = (1.0_pbl_wp+rtnp(i))
     &                                          *f1(i,k)
            else
              if(f1(i,k) < 0.0_pbl_wp) f1(i,k) = (1.0_pbl_wp+rtnp(i))
     &                                          *f1(i,k)
              if(f1(i,k) > 0.0_pbl_wp) f1(i,k) = 0.0_pbl_wp
            endif
          endif
        enddo
      enddo
c
!>  - Recover the tendency of TKE or TTE
c
      do k = 1,km
         do i = 1,im
!           f1(i,k) = max(f1(i,k), tkmin)
            qtend = (f1(i,k)-q1_loc(i,k,ntke))*real(rdt, kind=pbl_wp)
            rtg(i,k,ntke) = rtg(i,k,ntke) + real(qtend, kind=kind_phys)
         enddo
      enddo
      if(ldiag3d) then
         idtend = dtidx(ntke+100,index_of_process_pbl)
         if(idtend>0) then
            dtend(1:im,1:km,idtend) = dtend(1:im,1:km,idtend) +
     &                                real((f1(1:im,1:km) -
     &                                q1_loc(1:im,1:km,ntke)) *
     &                                real(rdt, kind=pbl_wp),
     &                                kind=kind_phys)
         endif
      endif
c
!> ## Compute tridiagonal matrix elements for heat and moisture
c
      do i=1,im
         ad(i,1) = 1.0_pbl_wp
         f1(i,1) = t1_loc(i,1)   + dtdz1(i) * real(heat(i), kind=pbl_wp)
         f2(i,1) = q1_loc(i,1,1) + dtdz1(i) * real(evap(i), kind=pbl_wp)
      enddo
      if(ntrac1 >= 2) then
        do n = 2, ntrac1
          is = (n-1) * km
          do i = 1, im
            f2(i,1+is) = q1_loc(i,1,n)
          enddo
        enddo
      endif
c
      do k = 1,km1
        do i = 1,im
          dtodsd  = dt2/real(del(i,k), kind=pbl_wp)
          dtodsu  = dt2/real(del(i,k+1), kind=pbl_wp)
          dsig    = real(prsl(i,k), kind=pbl_wp)-real(prsl(i,k+1),
     &              kind=pbl_wp)
          rdz     = rdzt(i,k)
          tem1    = dsig * dkt(i,k) * rdz
          dsdzt   = tem1 * gocp
          if (use_lpt > 0) then
            dsdzt = dsdzt-tem1*elocp*(qliq(i,k+1)-qliq(i,k))*rdz
     &            -(1.0_pbl_wp+0.33_pbl_wp/2.5_pbl_wp)*tem1*elocp
     &            *(qice(i,k+1)-qice(i,k))*rdz
          endif
          dsdz2   = tem1 * rdz
          au(i,k) = -dtodsd*dsdz2
          al(i,k) = -dtodsu*dsdz2
          ad(i,k) = ad(i,k)-au(i,k)
          ad(i,k+1)= 1.0_pbl_wp-al(i,k)
          tem2    = dsig * rdz
!
          if(pcnvflg(i) .and. k < kpbl(i)) then
             ptem      = 0.5_pbl_wp * tem2 * xmf(i,k)
             ptem1     = dtodsd * ptem
             ptem2     = dtodsu * ptem
             tem       = t1_loc(i,k) + t1_loc(i,k+1)
             ptem      = tcko(i,k) + tcko(i,k+1)
             f1(i,k)   = f1(i,k)+dtodsd*dsdzt-(ptem-tem)*ptem1
             f1(i,k+1) = t1_loc(i,k+1)-dtodsu*dsdzt+(ptem-tem)*ptem2
             ptem      = qcko(i,k,1) + qcko(i,k+1,1)
             f2(i,k)   = f2(i,k) - ptem * ptem1
             f2(i,k+1) = q1_loc(i,k+1,1) + ptem * ptem2
          else
             f1(i,k)   = f1(i,k)+dtodsd*dsdzt
             f1(i,k+1) = t1_loc(i,k+1)-dtodsu*dsdzt
             f2(i,k+1) = q1_loc(i,k+1,1)
          endif
!
          if(scuflg(i)) then
            if(k >= mrad(i) .and. k < krad(i)) then
              ptem      = 0.5_pbl_wp * tem2 * xmfd(i,k)
              ptem1     = dtodsd * ptem
              ptem2     = dtodsu * ptem
              ptem      = tcdo(i,k) + tcdo(i,k+1)
              tem       = t1_loc(i,k) + t1_loc(i,k+1)
              f1(i,k)   = f1(i,k) + (ptem - tem) * ptem1
              f1(i,k+1) = f1(i,k+1) - (ptem - tem) * ptem2
              ptem      = qcdo(i,k,1) + qcdo(i,k+1,1)
              f2(i,k)   = f2(i,k) + ptem * ptem1
              f2(i,k+1) = f2(i,k+1) - ptem * ptem2
            endif
          endif
!
          kmx = max(kpbl(i), krad(i))
          if((pcnvflg(i) .or. scuflg(i)) .and. (k < kmx)) then
            ptem  = tem2 * (xmf(i,k) - xmfd(i,k))
            ptem1 = dtodsd * ptem
            ptem2 = dtodsu * ptem
            f2(i,k)   = f2(i,k) + q_half(i,k,1) * ptem1
            f2(i,k+1) = f2(i,k+1) - q_half(i,k,1) * ptem2
          endif
!
        enddo
      enddo
!
      if(ntrac1 >= 2) then
        do n = 2, ntrac1
          is = (n-1) * km
          do k = 1, km1
            do i = 1, im
              dtodsd = dt2/real(del(i,k), kind=pbl_wp)
              dtodsu = dt2/real(del(i,k+1), kind=pbl_wp)
              dsig  = real(prsl(i,k), kind=pbl_wp)-real(prsl(i,k+1),
     &                kind=pbl_wp)
              tem2  = dsig * rdzt(i,k)
!
              if(pcnvflg(i) .and. k < kpbl(i)) then
                ptem  = 0.5_pbl_wp * tem2 * xmf(i,k)
                ptem1 = dtodsd * ptem
                ptem2 = dtodsu * ptem
                ptem  = qcko(i,k,n) + qcko(i,k+1,n)
                f2(i,k+is) = f2(i,k+is) - ptem * ptem1
                f2(i,k+1+is)= q1_loc(i,k+1,n) + ptem * ptem2
              else
                f2(i,k+1+is) = q1_loc(i,k+1,n)
              endif
!
              if(scuflg(i)) then
                if(k >= mrad(i) .and. k < krad(i)) then
                  ptem  = 0.5_pbl_wp * tem2 * xmfd(i,k)
                  ptem1 = dtodsd * ptem
                  ptem2 = dtodsu * ptem
                  ptem  = qcdo(i,k,n) + qcdo(i,k+1,n)
                  f2(i,k+is)  = f2(i,k+is) + ptem * ptem1
                  f2(i,k+1+is)= f2(i,k+1+is) - ptem * ptem2
                endif
              endif
!
              kmx = max(kpbl(i), krad(i))
              if((pcnvflg(i) .or. scuflg(i)) .and. (k < kmx)) then
                ptem  = tem2 * (xmf(i,k) - xmfd(i,k))
                ptem1 = dtodsd * ptem
                ptem2 = dtodsu * ptem
                f2(i,k+is)   = f2(i,k+is) + q_half(i,k,n) * ptem1
                f2(i,k+1+is) = f2(i,k+1+is) - q_half(i,k,n) * ptem2
              endif
!
            enddo
          enddo
        enddo
      endif
c
!> - Call tridin() to solve tridiagonal problem for heat and moisture
c
      call tridin_satmedmf(im,km,ntrac1,al,ad,au,f1,f2,au_loc, f1_loc,
     &                     f2_loc)
!
! Negative moisture is set to zero after borrowing it from
!   positive values within the mass-flux transport layers
!
      do i = 1,im
        tsumn(i) = 0.0_pbl_wp
        tsump(i) = 0.0_pbl_wp
        rtnp(i) = 1.0_pbl_wp
      enddo
      do k = 1,kps
        do i = 1,im
          if(pcnvflg(i) .and. scuflg(i)) then
            kbx = 1
            kmx = max(kpbl(i), krad(i))
          elseif(pcnvflg(i) .and. .not. scuflg(i)) then
            kbx = 1
            kmx = kpbl(i)
          elseif(.not. pcnvflg(i) .and. scuflg(i)) then
            kbx = mrad(i)
            kmx = krad(i)
          endif
          if((pcnvflg(i) .or. scuflg(i)) .and.
     &       (k >= kbx .and. k <= kmx)) then
            tem = f2(i,k) * real(del(i,k), kind=pbl_wp) * real(gravi,
     &            kind=pbl_wp)
            if(f2(i,k) < 0.0_pbl_wp) tsumn(i) = tsumn(i) + tem
            if(f2(i,k) > 0.0_pbl_wp) tsump(i) = tsump(i) + tem
          endif
        enddo
      enddo
      do i = 1,im
        if(pcnvflg(i) .or. scuflg(i)) then
          if(tsump(i) > 0.0_pbl_wp .and. tsumn(i) < 0.0_pbl_wp) then
            if(tsump(i) > abs(tsumn(i))) then
              rtnp(i) = tsumn(i) / tsump(i)
            else
              rtnp(i) = tsump(i) / tsumn(i)
            endif
          endif
        endif
      enddo
      do k = 1,kps
        do i = 1,im
          if(pcnvflg(i) .and. scuflg(i)) then
            kbx = 1
            kmx = max(kpbl(i), krad(i))
          elseif(pcnvflg(i) .and. .not. scuflg(i)) then
            kbx = 1
            kmx = kpbl(i)
          elseif(.not. pcnvflg(i) .and. scuflg(i)) then
            kbx = mrad(i)
            kmx = krad(i)
          endif
          if((pcnvflg(i) .or. scuflg(i)) .and.
     &       (k >= kbx .and. k <= kmx)) then
            if(rtnp(i) < 0.0_pbl_wp) then
              if(tsump(i) > abs(tsumn(i))) then
                if(f2(i,k) < 0.0_pbl_wp) f2(i,k) = 0.0_pbl_wp
                if(f2(i,k) > 0.0_pbl_wp) f2(i,k) = (1.0_pbl_wp+rtnp(i))
     &                                            *f2(i,k)
              else
                if(f2(i,k) < 0.0_pbl_wp) f2(i,k) = (1.0_pbl_wp+rtnp(i))
     &                                            *f2(i,k)
                if(f2(i,k) > 0.0_pbl_wp) f2(i,k) = 0.0_pbl_wp
              endif
            endif
          endif
        enddo
      enddo
!
! To remove negative moistures which were leaked out of the mass-flux transport layers
!   by eddy diffusion or potential negative moistures from the diffusion scheme
!   especially due to downward surface latent heat flux during nighttime,
!   positive moistures are borrowed again now from the entire layers
!
      do i = 1,im
        tsumn(i) = 0.0_pbl_wp
        tsump(i) = 0.0_pbl_wp
        rtnp(i) = 1.0_pbl_wp
      enddo
      do k = 1,km
        do i = 1,im
          tem = f2(i,k) * real(del(i,k), kind=pbl_wp) * gravi
          if(f2(i,k) < 0.0_pbl_wp) tsumn(i) = tsumn(i) + tem
          if(f2(i,k) > 0.0_pbl_wp) tsump(i) = tsump(i) + tem
        enddo
      enddo
      do i = 1,im
        if(tsump(i) > 0.0_pbl_wp .and. tsumn(i) < 0.0_pbl_wp) then
          if(tsump(i) > abs(tsumn(i))) then
            rtnp(i) = tsumn(i) / tsump(i)
          else
            rtnp(i) = tsump(i) / tsumn(i)
          endif
        endif
      enddo
      do k = 1,km
        do i = 1,im
          if(rtnp(i) < 0.0_pbl_wp) then
            if(tsump(i) > abs(tsumn(i))) then
              if(f2(i,k) < 0.0_pbl_wp) f2(i,k) = 0.0_pbl_wp
              if(f2(i,k) > 0.0_pbl_wp) f2(i,k) = (1.0_pbl_wp+rtnp(i))
     &                                          *f2(i,k)
            else
              if(f2(i,k) < 0.0_pbl_wp) f2(i,k) = (1.0_pbl_wp+rtnp(i))
     &                                          *f2(i,k)
              if(f2(i,k) > 0.0_pbl_wp) f2(i,k) = 0.0_pbl_wp
            endif
          endif
        enddo
      enddo
!
! Negative hydrometeors & tracers are set to zero after
!  borrowing them from positive values within the mass-flux
!  transport layers
!
! For the negative liquid water, first borrow water from vapor
!   and then borrow it from the other layers if there is still
!   negative water
!
      if(ntrac1 >= 2) then
        is = (ntcw-1) * km
        do k = 1,kps
          do i = 1,im
            if(pcnvflg(i) .and. scuflg(i)) then
              kbx = 1
              kmx = max(kpbl(i), krad(i))
            elseif(pcnvflg(i) .and. .not. scuflg(i)) then
              kbx = 1
              kmx = kpbl(i)
            elseif(.not. pcnvflg(i) .and. scuflg(i)) then
              kbx = mrad(i)
              kmx = krad(i)
            endif
            if((pcnvflg(i) .or. scuflg(i)) .and.
     &         (k >= kbx .and. k <= kmx)) then
              if(f2(i,k+is) < 0.0_pbl_wp) then
                tem = f2(i,k) + f2(i,k+is)
                if(tem >= 0.0_pbl_wp) then
                  f2(i,k) = tem
                  f1(i,k) = f1(i,k) - elocp * f2(i,k+is)
                  f2(i,k+is) = 0.0_pbl_wp
                elseif (f2(i,k) > 0.0_pbl_wp) then
                  f2(i,k+is) = tem
                  f1(i,k) = f1(i,k) + elocp * f2(i,k)
                  f2(i,k) = 0.0_pbl_wp
                endif
              endif
            endif
          enddo
        enddo
      endif
!
! For the negative rain water, first borrow water from vapor
!   and then borrow it from the other layers if there is still
!   negative water
!
      if(ntrac1 >= 2 .and. ntrw > 0) then
        is = (ntrw-1) * km
        do k = 1,kps
          do i = 1,im
            if(pcnvflg(i) .and. scuflg(i)) then
              kbx = 1
              kmx = max(kpbl(i), krad(i))
            elseif(pcnvflg(i) .and. .not. scuflg(i)) then
              kbx = 1
              kmx = kpbl(i)
            elseif(.not. pcnvflg(i) .and. scuflg(i)) then
              kbx = mrad(i)
              kmx = krad(i)
            endif
            if((pcnvflg(i) .or. scuflg(i)) .and.
     &         (k >= kbx .and. k <= kmx)) then
              if(f2(i,k+is) < 0.0_pbl_wp) then
                tem = f2(i,k) + f2(i,k+is)
                if(tem >= 0.0_pbl_wp) then
                  f2(i,k) = tem
                  f1(i,k) = f1(i,k) - elocp * f2(i,k+is)
                  f2(i,k+is) = 0.0_pbl_wp
                elseif (f2(i,k) > 0.0_pbl_wp) then
                  f2(i,k+is) = tem
                  f1(i,k) = f1(i,k) + elocp * f2(i,k)
                  f2(i,k) = 0.0_pbl_wp
                endif
              endif
            endif
          enddo
        enddo
      endif
!
      if(ntrac1 >= 2) then
        do n = 2, ntrac1
          is = (n-1) * km
!
          do i = 1,im
            tsumn(i) = 0.0_pbl_wp
            tsump(i) = 0.0_pbl_wp
            rtnp(i) = 1.0_pbl_wp
          enddo
          do k = 1,kps
            do i = 1,im
              if(pcnvflg(i) .and. scuflg(i)) then
                kbx = 1
                kmx = max(kpbl(i), krad(i))
              elseif(pcnvflg(i) .and. .not. scuflg(i)) then
                kbx = 1
                kmx = kpbl(i)
              elseif(.not. pcnvflg(i) .and. scuflg(i)) then
                kbx = mrad(i)
                kmx = krad(i)
              endif
              if((pcnvflg(i) .or. scuflg(i)) .and.
     &           (k >= kbx .and. k <= kmx)) then
                tem = f2(i,k+is) * real(del(i,k), kind=pbl_wp) * gravi 
                if(f2(i,k+is) < 0.0_pbl_wp) tsumn(i) = tsumn(i) + tem
                if(f2(i,k+is) > 0.0_pbl_wp) tsump(i) = tsump(i) + tem
              endif
            enddo
          enddo
          do i = 1,im
            if(pcnvflg(i) .or. scuflg(i)) then
              if(tsump(i) > 0.0_pbl_wp .and. tsumn(i) < 0.0_pbl_wp) then
                if(tsump(i) > abs(tsumn(i))) then
                  rtnp(i) = tsumn(i) / tsump(i)
                else
                  rtnp(i) = tsump(i) / tsumn(i)
                endif
              endif
            endif
          enddo
          do k = 1,kps
            do i = 1,im
              if(pcnvflg(i) .and. scuflg(i)) then
                kbx = 1
                kmx = max(kpbl(i), krad(i))
              elseif(pcnvflg(i) .and. .not. scuflg(i)) then
                kbx = 1
                kmx = kpbl(i)
              elseif(.not. pcnvflg(i) .and. scuflg(i)) then
                kbx = mrad(i)
                kmx = krad(i)
              endif
              if((pcnvflg(i) .or. scuflg(i)) .and.
     &           (k >= kbx .and. k <= kmx)) then
                if(rtnp(i) < 0.0_pbl_wp) then
                  if(tsump(i) > abs(tsumn(i))) then
                    if(f2(i,k+is)<0.0_pbl_wp) f2(i,k+is)=0.0_pbl_wp
                    if(f2(i,k+is)>0.0_pbl_wp) f2(i,k+is)=(1.0_pbl_wp
     &                                              +rtnp(i))*f2(i,k+is)
                  else
                    if(f2(i,k+is)<0.0_pbl_wp) f2(i,k+is)=(1.0_pbl_wp
     &                                              +rtnp(i))*f2(i,k+is)
                    if(f2(i,k+is)>0.0_pbl_wp) f2(i,k+is)=0.0_pbl_wp
                  endif
                endif
              endif
            enddo
          enddo
!
! To remove negative hydrometeors & tracers which were leaked out of the mass-flux transport layers
!   by eddy diffusion or potential negative hydrometeors & tracers from the diffusion scheme
!   especially due to downward surface fluxes during nighttime,
!   positive hydrometeors & tracers are borrowed again now from the entire layers
!
          do i = 1,im
            tsumn(i) = 0.0_pbl_wp
            tsump(i) = 0.0_pbl_wp
            rtnp(i) = 1.0_pbl_wp
          enddo
          do k = 1,km
            do i = 1,im
              tem = f2(i,k+is) * real(del(i,k), kind=pbl_wp) * gravi
              if(f2(i,k+is) < 0.0_pbl_wp) tsumn(i) = tsumn(i) + tem
              if(f2(i,k+is) > 0.0_pbl_wp) tsump(i) = tsump(i) + tem
            enddo
          enddo
          do i = 1,im
            if(tsump(i) > 0.0_pbl_wp .and. tsumn(i) < 0.0_pbl_wp) then
              if(tsump(i) > abs(tsumn(i))) then
                rtnp(i) = tsumn(i) / tsump(i)
              else
                rtnp(i) = tsump(i) / tsumn(i)
              endif
            endif
          enddo
          do k = 1,km
            do i = 1,im
              if(rtnp(i) < 0.0_pbl_wp) then
                if(tsump(i) > abs(tsumn(i))) then
                  if(f2(i,k+is)<0.0_pbl_wp) f2(i,k+is)=0.0_pbl_wp
                  if(f2(i,k+is)>0.0_pbl_wp) f2(i,k+is)=(1.0_pbl_wp
     &               +rtnp(i))*f2(i,k+is)
                else
                  if(f2(i,k+is)<0.0_pbl_wp) f2(i,k+is)=(1.0_pbl_wp
     &               +rtnp(i))*f2(i,k+is)
                  if(f2(i,k+is)>0.0_pbl_wp) f2(i,k+is)=0.0_pbl_wp
                endif
              endif
            enddo
          enddo
!
        enddo
      endif
c
!> - Recover the tendencies of heat and moisture
c
      do  k = 1,km
         do i = 1,im
            ten_t(i,k)   = real((f1(i,k)-t1_loc(i,k))*rdt,
     &                     kind=kind_phys)
            qtend        = (f2(i,k)-q1_loc(i,k,1))*rdt
            rtg(i,k,1)   = rtg(i,k,1) + real(qtend, kind=kind_phys)
!           dtsfc(i)   = dtsfc(i)+cont*del(i,k)*ttend
!           dqsfc(i)   = dqsfc(i)+conq*del(i,k)*qtend
         enddo
      enddo
!
      do i = 1,im
        dtsfc(i) = real(rho_a(i), kind=kind_phys) * cp * heat(i)
        dqsfc(i) = real(rho_a(i), kind=kind_phys) * hvap * evap(i)
      enddo
!
      if(ldiag3d .and. .not. gen_tend) then
        idtend = dtidx(index_of_temperature,index_of_process_pbl)
        if(idtend>=1) then
          do  k = 1,km
            do i = 1,im
              dtend(i,k,idtend) = dtend(i,k,idtend) + ten_t(i,k) * delt
            enddo
          enddo
        endif
        ! Send tendencies just for QV; other tracers are below.
        idtend = dtidx(100+ntqv,index_of_process_pbl)
        if(idtend>=1) then
          do  k = 1,km
            do i = 1,im
              qtend      = (f2(i,k)-q1_loc(i,k,1))*rdt
              dtend(i,k,idtend) = dtend(i,k,idtend) + real(qtend,
     &                            kind=kind_phys) * delt
            enddo
          enddo
        endif
      endif
!
      if(ntrac1 >= 2) then
        do n = 2, ntrac1
          is = (n-1) * km
          do k = 1, km
            do i = 1, im
              qtend = (f2(i,k+is)-q1_loc(i,k,n))*rdt
              rtg(i,k,n) = rtg(i,k,n)+real(qtend, kind=kind_phys)
            enddo
          enddo
        enddo
        if(ldiag3d .and. .not. gen_tend) then
          ! Send tendencies for all tracers that were selected.
          do n = 2, ntrac1
            is = (n-1) * km
            idtend = dtidx(n+100,index_of_process_pbl)
            if(idtend>=1) then
              if(n/=ntke) then
                do k = 1, km
                  do i = 1, im
                    qtend = (f2(i,k+is)-q1_loc(i,k,n))*rdt
                    dtend(i,k,idtend) = dtend(i,k,idtend)+real(qtend,
     &                                  kind=kind_phys)*delt
                  enddo
                enddo
              endif
            endif
          enddo
        endif
      endif
!
!> ## Add TKE dissipative heating to temperature tendency
!
      if(dspheat) then
        do k = 1,km1
          do i = 1,im
!           tem = min(diss(i,k), dspmax)
!           ttend = tem / cp
            ten_t(i,k) = ten_t(i,k) + dspfac * real(diss(i,k),
     &                   kind=kind_phys) / cp
          enddo
        enddo
        if(ldiag3d .and. .not. gen_tend) then
          idtend = dtidx(index_of_temperature,index_of_process_pbl)
          if(idtend>=1) then
            do k = 1,km1
              do i = 1,im
                dtend(i,k,idtend) = dtend(i,k,idtend)+
     &                              (dspfac * real(diss(i,k),
     &                              kind=kind_phys) / cp)*delt
              enddo
            enddo
          endif
        endif
      endif
c
!> ## Compute tridiagonal matrix elements for momentum
c
      do i=1,im
         ad(i,1) = 1.0_pbl_wp + dtdz1(i) * real(stress(i), kind=pbl_wp)
     &           / real(spd1(i), kind=pbl_wp)
         f1(i,1) = u1_loc(i,1)
         f2(i,1) = v1_loc(i,1)
      enddo
c
      do k = 1,km1
        do i=1,im
          dtodsd  = dt2/real(del(i,k), kind=pbl_wp)
          dtodsu  = dt2/real(del(i,k+1), kind=pbl_wp)
          dsig    = real(prsl(i,k), kind=pbl_wp)-real(prsl(i,k+1),
     &              kind=pbl_wp)
          rdz     = rdzt(i,k)
          tem1    = dsig * dku(i,k) * rdz
          dsdz2   = tem1*rdz
          au(i,k) = -dtodsd*dsdz2
          al(i,k) = -dtodsu*dsdz2
          ad(i,k) = ad(i,k)-au(i,k)
          ad(i,k+1)= 1.0_pbl_wp-al(i,k)
          tem2    = dsig * rdz
!
          if(pcnvflg(i) .and. k < kpbl(i)) then
             ptem      = 0.5_pbl_wp * tem2 * xmf(i,k)
             ptem1     = dtodsd * ptem
             ptem2     = dtodsu * ptem
             tem       = u1_loc(i,k) + u1_loc(i,k+1)
             ptem      = ucko(i,k) + ucko(i,k+1)
             f1(i,k)   = f1(i,k) - (ptem - tem) * ptem1
             f1(i,k+1) = u1_loc(i,k+1) + (ptem - tem) * ptem2
             tem       = v1_loc(i,k) + v1_loc(i,k+1)
             ptem      = vcko(i,k) + vcko(i,k+1)
             f2(i,k)   = f2(i,k) - (ptem - tem) * ptem1
             f2(i,k+1) = v1_loc(i,k+1) + (ptem - tem) * ptem2
          else
             f1(i,k+1) = u1_loc(i,k+1)
             f2(i,k+1) = v1_loc(i,k+1)
          endif
!
          if(scuflg(i)) then
            if(k >= mrad(i) .and. k < krad(i)) then
              ptem      = 0.5_pbl_wp * tem2 * xmfd(i,k)
              ptem1     = dtodsd * ptem
              ptem2     = dtodsu * ptem
              tem       = u1_loc(i,k) + u1_loc(i,k+1)
              ptem      = ucdo(i,k) + ucdo(i,k+1)
              f1(i,k)   = f1(i,k) + (ptem - tem) *ptem1
              f1(i,k+1) = f1(i,k+1) - (ptem - tem) *ptem2
              tem       = v1_loc(i,k) + v1_loc(i,k+1)
              ptem      = vcdo(i,k) + vcdo(i,k+1)
              f2(i,k)   = f2(i,k) + (ptem - tem) * ptem1
              f2(i,k+1) = f2(i,k+1) - (ptem - tem) * ptem2
            endif
          endif
!
        enddo
      enddo
c
!> - Call tridi2() to solve tridiagonal problem for momentum
c
      call tridi2_satmedmf(im,km,al,ad,au,f1,f2,au_loc,f1_loc,f2_loc)
c
!> - Recover the tendencies of momentum
c
      do k = 1,km
         do i = 1,im
            ten_u(i,k) = real((f1(i,k)-u1_loc(i,k))*rdt, kind=kind_phys)
            ten_v(i,k) = real((f2(i,k)-v1_loc(i,k))*rdt, kind=kind_phys)
!           dusfc(i) = dusfc(i)+conw*del(i,k)*utend
!           dvsfc(i) = dvsfc(i)+conw*del(i,k)*vtend
         enddo
      enddo
      if (use_oceanuv) then
        do i = 1,im
          spd1_m = real(sqrt((u1(i,1)-usfco(i))**2+(v1(i,1)
     &           - vsfco(i))**2), kind=pbl_wp)
          dusfc(i) = -1.0_kind_phys * real(rho_a(i), kind=kind_phys)
     &             * stress(i) * (u1(i,1) - usfco(i)) / real(spd1_m,
     &                 kind=kind_phys)
          dvsfc(i) = -1.0_kind_phys * real(rho_a(i), kind=kind_phys)
     &             * stress(i) * (v1(i,1) - vsfco(i)) / real(spd1_m,
     &                 kind=kind_phys)
        enddo
      else
        do i = 1,im
          dusfc(i) = -1.0_kind_phys * real(rho_a(i), kind=kind_phys)
     &             * stress(i) *u1(i,1) / spd1(i)
          dvsfc(i) = -1.0_kind_phys * real(rho_a(i), kind=kind_phys)
     &             * stress(i) *v1(i,1) / max(spd1(i)
        enddo
      endif
!
      if(ldiag3d .and. .not. gen_tend) then
        idtend = dtidx(index_of_x_wind,index_of_process_pbl)
        if(idtend>=1) then
          do k = 1,km
            do i = 1,im
              dtend(i,k,idtend) = dtend(i,k,idtend) + ten_u(i,k)*delt
            enddo
          enddo
        endif

        idtend = dtidx(index_of_y_wind,index_of_process_pbl)
        if(idtend>=1) then
          do k = 1,km
            do i = 1,im
              dtend(i,k,idtend) = dtend(i,k,idtend) + ten_v(i,k)*delt
            enddo
          enddo
        endif
      endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> ## Save PBL height for diagnostic purpose
!
      do i = 1, im
         hpbl(i) = real(hpblx(i), kind=kind_phys)
         kpbl(i) = kpblx(i)
      enddo
      if(sa3dtke) then
        do k = 1, km
          do i = 1, im
            dku3d_h(i,k)  = real(dku_h(i,k), kind=kind_phys)   ! pass dku3d_h to dyn_core
            dku3d_e(i,k)  = real(dkq_h(i,k), kind=kind_phys)   ! pass dku3d_e to dyn_core
          enddo
        enddo
      endif  !sa3dtke
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (do_canopy) then
        if(allocated(EDDYVESTX)) deallocate(EDDYVESTX)
        if(allocated(ZCANX)) deallocate(ZCANX)
        if(allocated(ZOOOX)) deallocate(ZOOOX)
      endif
      
      return
      end subroutine satmedmfvdifq_run
!> @}
      end module satmedmfvdifq
