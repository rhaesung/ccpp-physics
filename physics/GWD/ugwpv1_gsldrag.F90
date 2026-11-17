!>  \file ugwpv1_gsldrag.F90

!> This module introduces two gravity wave drag schemes: UGWPv1 and orographic drag scheme
!!
!!      1) The "V1 CIRES UGWP" scheme as tested in the FV3GFSv16-127L atmosphere model and workflow, which includes:
!!            a) the orograhic gravity wave drag, flow blocking scheme and TOFD (Beljaars et al, 2004).
!!            b) the v1 CIRE ugwp non-stationary GW scheme, new revision that generate realistic climate of FV3GFS-127L
!!               in the strato-mesosphere in the multi-year simulations (Annual cycles, SAO and QBO in th tropical dynamics).
!!      2) The GSL orographic drag suite (drag_suite.F90), as implmeneted in the RAP/HRRR, which includes:
!!            a) large-scale gravity wave drag and low-level flow blocking -- active at horizontal scales
!!               down to ~5km (Kim and Arakawa, 1995 \cite kim_and_arakawa_1995; Kim and Doyle, 2005 \cite kim_and_doyle_2005)
!!            b) small-scale gravity wave drag scheme -- active typically in stable PBL at horizontal grid resolutions down to ~1km
!!               (Steeneveld et al, 2008 \cite steeneveld_et_al_2008; Tsiringakis et al, 2017 \cite tsiringakis_et_al_2017)
!!            c) turbulent orographic form drag -- active at horizontal grid ersolutions down to ~1km
!!               (Beljaars et al, 2004 \cite beljaars_et_al_2004)
!! See Valery Yudin's presentation at 2020 UFS User's meeting (Jul 2020):
!! Gravity waves (GWs): Mesoscale GWs transport momentum, energy (heat) , and create eddy mixing in the whole atmosphere domain; Breaking and dissipating GWs deposit: (a) momentum; (b) heat (energy); and create (c) turbulent mixing of momentum, heat, and tracers
!! To properly incorporate GW effects (a-c) unresolved by DYCOREs we need GW physics
!! "Unified": a) all GW effects due to both dissipation/breaking; b) identical GW solvers for all GW sources; c) ability to replace solvers.
!! Unified Formalism:
!! 1. GW Sources: Stochastic and physics based mechanisms for GW-excitations in the lower atmosphere, calibrated by the high-res analyses/forecasts, and observations (3 types of GW sources: orography, convection, fronts/jets).
!! 2. GW Propagation: Unified solver for "propagation, dissipation and breaking" excited from all type of GW sources.
!! 3. GW Effects: Unified representation of GW impacts on the "resolved" flow for all sources (energy-balanced schemes for momentum, heat and mixing).
!! https://www.weather.gov/media/sti/nggps/Presentations%202017/02%20NGGPS_VYUDIN_2017_.pdf
!!
!! The ugwpv1_gsldrag scheme is activated by gwd_opt = 2 in the namelist.
!! The choice of schemes is activated at runtime by the following namelist options (boolean):
!! NA    do_ugwp_v0           -- activates V0 CIRES UGWP scheme - both orographic and non-stationary GWD is not active (NA)
!! NA    do_ugwp_v0_orog_only -- activates V0 CIRES UGWP scheme - orographic GWD only
!!       do_gsl_drag_ls_bl    -- activates RAP/HRRR (GSL) large-scale OGWD and blocking
!!       do_gsl_drag_ss       -- activates RAP/HRRR (GSL) small-scale OGWD
!!       do_gsl_drag_tofd     -- activates RAP/HRRR (GSL) turbulent orographic drag
!!       do_ugwp_v1           -- activates V1 CIRES UGWP scheme - both orographic and non-stationary GWD
!!       do_ugwp_v1_orog_only -- activates V1 CIRES UGWP scheme - orographic GWD only
!!       do_ugwp_v1_w_gsldrag -- activates V1 CIRES UGWP scheme with orographic drag of GSL
!! Note that only one "large-scale" scheme can be activated at a time.
!!
module ugwpv1_gsldrag

    use machine,               only: kind_phys
    use mo_ugwpv1_kind,        only: wp

    use cires_ugwpv1_triggers, only:  slat_geos5_2020, slat_geos5_tamp_v1
    use cires_ugwpv1_module,   only:  cires_ugwpv1_init, ngwflux_update, calendar_ugwp
    use cires_ugwpv1_module,   only:  knob_ugwp_version, cires_ugwp_dealloc, tamp_mpa
    use cires_ugwpv1_solv2,    only:  cires_ugwpv1_ngw_solv2
!    use cires_ugwpv1_solv2,     only:  cires_ugwpv1_ngw_solv2, ecmwf_ngw
     use ecmwf_ngw,             only:  ecmwf_ngw_emc

    use cires_ugwpv1_oro,      only:  orogw_v1

    use drag_suite_ugwpv1,     only:  drag_suite_ugwpv1_run, drag_suite_ugwpv1_psl

    implicit none

    private

    public ugwpv1_gsldrag_init, ugwpv1_gsldrag_run, ugwpv1_gsldrag_final

    logical :: is_initialized = .False.

contains

! ------------------------------------------------------------------------
! CCPP entry points for CIRES Unified Gravity Wave Physics (UGWP) scheme v0
! ------------------------------------------------------------------------
!>@brief The subroutine initializes the unified UGWP
!> \section arg_table_ugwpv1_gsldrag_init Argument Table
!! \htmlinclude ugwpv1_gsldrag_init.html
!!
    subroutine ugwpv1_gsldrag_init  (                                          &
                me, master, nlunit, input_nml_file, logunit,                   &
                fn_nml2, jdat, lonr, levs, ak, bk, dtp,                        &
                con_pi, con_rerth, con_p0,                                     &
                con_g, con_omega,  con_cp, con_rd, con_rv,con_fvirt,           &
                do_ugwp,do_ugwp_v0, do_ugwp_v0_orog_only, do_gsl_drag_ls_bl,   &
                  do_gsl_drag_ss, do_gsl_drag_tofd, do_ngw_ec, do_ugwp_v1,     &
!!                do_gsl_drag_ss, do_gsl_drag_tofd, do_ugwp_v1,                  &
                do_ugwp_v1_orog_only, do_ugwp_v1_w_gsldrag, errmsg, errflg)

    use ugwp_common

!----  initialization of unified_ugwp
    implicit none

    integer,              intent (in) :: me
    integer,              intent (in) :: master
    integer,              intent (in) :: nlunit
    character(len=*),     intent (in) :: input_nml_file(:)
    integer,              intent (in) :: logunit
    integer,              intent (in) :: jdat(:)
    integer,              intent (in) :: lonr
    integer,              intent (in) :: levs
    real(kind=kind_phys), intent (in) :: ak(:), bk(:)
    real(kind=kind_phys), intent (in) :: dtp

    real(kind=kind_phys), intent (in) :: con_p0, con_pi, con_rerth
    real(kind=kind_phys), intent (in) :: con_g, con_cp, con_rd, con_rv, con_omega, con_fvirt
    logical,              intent (in) :: do_ugwp

    logical,              intent (in) :: do_ugwp_v0, do_ugwp_v0_orog_only,         &
                                         do_gsl_drag_ls_bl, do_gsl_drag_ss,        &
!!                                       do_gsl_drag_tofd, do_ugwp_v1,             &
                                         do_gsl_drag_tofd, do_ugwp_v1, do_ngw_ec,  &
                                         do_ugwp_v1_orog_only,do_ugwp_v1_w_gsldrag

    character(len=*), intent (in)  :: fn_nml2
    !character(len=*), parameter   :: fn_nml='input.nml'

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Local working-precision buffers for 32-bit internal initialization
    real(kind=wp) :: ak_loc(levs+1), bk_loc(levs+1), dtp_loc
    real(kind=wp) :: con_p0_loc, con_pi_loc, con_rerth_loc

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
!============================================================================
!
!     gwd_opt => "1 and 2, 3, 22, 33' see previous GSL-commits
!                 related to GSL-oro drag suite
!    for use of the new-GSL/old-GFS/EMC inputs for sub-grid orography
!  see details inside /ufs-weather-model/FV3/io/FV3GFS_io.F90
!  FV3GFS_io.F90:    if (Model%gwd_opt==3 .or. Model%gwd_opt==33 .or. &
!  FV3GFS_io.F90:        Model%gwd_opt==2 .or. Model%gwd_opt==22 ) then
!  FV3GFS_io.F90:          if ( (Model%gwd_opt==3 .or. Model%gwd_opt==33) .or.    &
!  FV3GFS_io.F90:               ( (Model%gwd_opt==2 .or. Model%gwd_opt==22) .and. &
!
! gwd_opt=1 -current 14-element GFS-EMC subgrid-oro input
! gwd_opt=2 and 3    24-element -current 14-element GFS-EMC subgrid-oro input
! GSL uses the gwd_opt flag to control "extra" diagnostics (22 and 33)
! CCPP may use gwd_opt to determine 14 or 24 variables for the input
!      but at present you work with "nmtvr"
! GFS_GWD_generic.F90:      integer, intent(in) :: im, levs, nmtvr
!GFS_GWD_generic.F90:      real(kind=kind_phys), intent(in) :: mntvar(im,nmtvr)
!GFS_GWD_generic.F90:      if (nmtvr == 14) then  ! gwd_opt=1 current operational - as of 2014
!GFS_GWD_generic.F90:      elseif (nmtvr == 10) then ????
!GFS_GWD_generic.F90:      elseif (nmtvr == 6) then  ????
!GFS_GWD_generic.F90:      elseif (nmtvr == 24) then   ! GSD_drag_suite and unified_ugwp gwd_opt=2,3
!
! 1) gsldrag:   do_gsl_drag_ls_bl, do_gsl_drag_ss, do_gsl_drag_tofd, do_ugwp_v1
! 2) CIRES-v1:  do_ugwp_v1,        do_ugwp_v1_orog_only,  do_tofd,   ldiag_ugwp
!==============================================================================
    ! Test to make sure that at most only one large-scale/blocking
    ! orographic drag scheme is chosen
    if ( (do_ugwp_v0.and.(do_ugwp_v0_orog_only.or.do_gsl_drag_ls_bl.or.    &
                          do_ugwp_v1.or.do_ugwp_v1_orog_only))        .or. &
         (do_ugwp_v0_orog_only.and.(do_gsl_drag_ls_bl.or.do_ugwp_v1.or.    &
                                    do_ugwp_v1_orog_only))            .or. &
         (do_gsl_drag_ls_bl.and.do_ugwp_v1_orog_only)  ) then

       write(errmsg,'(*(a))') "Logic error: Only one large-scale&
          &/blocking scheme (do_ugwp_v0,do_ugwp_v0_orog_only,&
          &do_gsl_drag_ls_bl,do_ugwp_v1 or &
          &do_ugwp_v1_orog_only) can be chosen"
       errflg = 1
       return

    end if

    if ( do_ugwp_v0_orog_only .or. do_ugwp_v0) then
       print *,  ' ccpp do_ugwp_v0 active ', do_ugwp_v0
       print *,  ' ccpp do_ugwp_v1_orog_only active ', do_ugwp_v0_orog_only
       write(errmsg,'(*(a))') " the CIRES <ugwpv1_gsldrag> CCPP-suite does not &
         support <ugwp_v0> schemes "
       errflg = 1
       return
    endif

    if (do_ugwp_v1_w_gsldrag .and. do_ugwp_v1_orog_only ) then

       print *,  '  do_ugwp_v1_w_gsldrag ', do_ugwp_v1_w_gsldrag
       print *,  '  do_ugwp_v1_orog_only ', do_ugwp_v1_orog_only
       print *,  '  do_gsl_drag_ls_bl ',do_gsl_drag_ls_bl
       write(errmsg,'(*(a))') " the CIRES <ugwpv1_gsldrag> CCPP-suite intend to &
         support <ugwp_v1> with <gsldrag>  but  has Logic error"
       errflg = 1
       return
    endif
!==========================
!
! initialize ugwp_common
!   con_pi, con_rerth, con_p0,  con_g, con_omega,  con_cp, con_rd, con_rv,con_fvirt
!
!==========================

    ak_loc        = real(ak(1:levs+1), kind=wp)
    bk_loc        = real(bk(1:levs+1), kind=wp)
    dtp_loc       = real(dtp,          kind=wp)
    con_p0_loc    = real(con_p0,       kind=wp)
    con_pi_loc    = real(con_pi,       kind=wp)
    con_rerth_loc = real(con_rerth,    kind=wp)

    pi    = con_pi_loc
    arad  = con_rerth_loc
    p0s   = con_p0_loc
    grav  = real(con_g,     kind=wp)
    omega1= real(con_omega, kind=wp)
    cpd   = real(con_cp,    kind=wp)
    rd    = real(con_rd,    kind=wp)
    rv    = real(con_rv,    kind=wp)
    fv    = real(con_fvirt, kind=wp)

    grav2  = grav + grav; rgrav  = 1.0_wp/grav ; rgrav2 = rgrav*rgrav
    rdi    = 1.0_wp / rd ; rcpd = 1.0_wp/cpd;  rcpd2 = 0.5_wp/cpd
    gor    = grav/rd
    gr2    = grav*gor
    grcp   = grav*rcpd
    gocp   = grcp
    rcpdl  = cpd*rgrav
    grav2cpd = grav*grcp

    pi2      = 2.0_wp*pi ;  pih = 0.5_wp*pi
    rad_to_deg=180.0_wp/pi
    deg_to_rad=pi/180.0_wp

    bnv2min = (pi2/1800.0_wp)*(pi2/1800.0_wp)
    bnv2max = (pi2/30.0_wp)*(pi2/30.0_wp)
    dw2min  = 1.0_wp
    velmin  = sqrt(dw2min)
    minvel  = 0.5_wp

    omega2  = 2.0_wp*omega1
    omega3  = 3.0_wp*omega1

    hpscale = 7000.0_wp ; hpskm = hpscale*1.0e-3_wp
    rhp     = 1.0_wp/hpscale
    rhp2 = 0.5_wp*rhp; rh4 = 0.25_wp*rhp
    rhp4 = rhp2 * rhp2
    khp  = rhp* rd/cpd
    mkzmin  = pi2/80.0e3_wp
    mkz2min = mkzmin*mkzmin
    mkzmax  = pi2/500.0_wp
    mkz2max = mkzmax*mkzmax
    cdmin   = 2.0e-2_wp/mkzmax

    rcpdt  = rcpd/dtp_loc

    if ( do_ugwp_v1 ) then
       call cires_ugwpv1_init (me, master, nlunit, logunit, jdat, con_pi_loc,  &
                               con_rerth_loc, fn_nml2, input_nml_file, lonr,   &
                               levs, ak_loc, bk_loc, con_p0_loc, dtp_loc,      &
                               errmsg, errflg)
       if (errflg/=0) return
    end if

    if (me == master) then
       print *,  ' ccpp: ugwpv1_gsldrag_init   '

       print *,  ' ccpp do_ugwp_v1  flag ', do_ugwp_v1
       print *,  ' ccpp do_gsl_drag_ls_bl  flag ',    do_gsl_drag_ls_bl
       print *,  ' ccpp do_gsl_drag_ss  flag ' ,      do_gsl_drag_ss
       print *,  ' ccpp do_gsl_drag_tofd  flag ',     do_gsl_drag_tofd

       print *, ' ccpp: ugwpv1_gsldrag_init  '
    endif

    is_initialized = .true.

    end subroutine ugwpv1_gsldrag_init

! -----------------------------------------------------------------------
! finalize of ugwpv1_gsldrag   (_final)
! -----------------------------------------------------------------------

!>@brief The subroutine finalizes the CIRES UGWP

!> \section arg_table_ugwpv1_gsldrag_final Argument Table
!! \htmlinclude ugwpv1_gsldrag_final.html
!!
    subroutine ugwpv1_gsldrag_final(errmsg, errflg)

    implicit none

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not.is_initialized) return

    call cires_ugwp_dealloc

    is_initialized = .false.

    end subroutine ugwpv1_gsldrag_final

! -----------------------------------------------------------------------
!    originally from ugwp_driver_v0.f
!    driver of cires_ugwp   (_driver)
! -----------------------------------------------------------------------
!   driver is called after pbl & before chem-parameterizations
! -----------------------------------------------------------------------
!  order = dry-adj=>conv=mp-aero=>radiation -sfc/land- chem -> vertdiff-> [rf-gws]=> ion-re
! -----------------------------------------------------------------------
!>\section gen_ugwpv1_gsldrag_run Unified Gravity Wave Physics General Algorithm
!! The physics of NGWs in the UGWP framework (Yudin et al. 2018 \cite yudin_et_al_2018) is represented by four GW-solvers, which is introduced in Lindzen (1981) \cite lindzen_1981, Hines (1997) \cite hines_1997, Alexander and Dunkerton (1999) \cite alexander_and_dunkerton_1999, and Scinocca (2003) \cite scinocca_2003. The major modification of these GW solvers is represented by the addition of the background dissipation of temperature and winds to the saturation criteria for wave breaking. This feature is important in the mesosphere and thermosphere for WAM applications and it considers appropriate scale-dependent dissipation of waves near the model top lid providing the momentum and energy conservation in the vertical column physics (Shaw and Shepherd 2009 \cite shaw_and_shepherd_2009). In the UGWP-v0, the modification of Scinocca (2003) \cite scinocca_2003 scheme for NGWs with non-hydrostatic and rotational effects for GW propagations and background dissipation is represented by the subroutine fv3_ugwp_solv2_v0. In the next release of UGWP, additional GW-solvers will be implemented along with physics-based triggering of waves and stochastic approaches for selection of GW modes characterized by horizontal phase velocities, azimuthal directions and magnitude of the vertical momentum flux (VMF).
!!
!! In UGWP-v0, the specification for the VMF function is adopted from the GEOS-5 global atmosphere model of GMAO NASA/GSFC, as described in Molod et al. (2015) \cite molod_et_al_2015 and employed in the MERRRA-2 reanalysis (Gelaro et al., 2017 \cite gelaro_et_al_2017). The Fortran subroutine \ref slat_geos5_tamp describes the latitudinal shape of VMF-function as displayed in Figure 3 of Molod et al. (2015) \cite molod_et_al_2015. It shows that the enhanced values of VMF in the equatorial region gives opportunity to simulate the QBO-like oscillations in the equatorial zonal winds and lead to more realistic simulations of the equatorial dynamics in GEOS-5 operational and MERRA-2 reanalysis products. For the first vertically extended version of FV3GFS in the stratosphere and mesosphere, this simplified function of VMF allows us to tune the model climate and to evaluate multi-year simulations of FV3GFS with the MERRA-2 and ERA-5 reanalysis products, along with temperature, ozone, and water vapor observations of current satellite missions. After delivery of the UGWP-code, the EMC group developed and tested approach to modulate the zonal mean NGW forcing by 3D-distributions of the total precipitation as a proxy for the excitation of NGWs by convection and the vertically-integrated  (surface - tropopause) Turbulent Kinetic Energy (TKE). The verification scores with updated NGW forcing, as reported elsewhere by EMC researchers, display noticeable improvements in the forecast scores produced by FV3GFS configuration extended into the mesosphere.
!!
!> \section arg_table_ugwpv1_gsldrag_run Argument Table
!! \htmlinclude ugwpv1_gsldrag_run.html
!!
     subroutine ugwpv1_gsldrag_run(me, master, im, levs, ak, bk, ntrac, lonr, dtp,      &
          kdt, ldiag3d, lssav, flag_for_gwd_generic_tend, do_gsl_drag_ls_bl,            &
          do_gsl_drag_ss, do_gsl_drag_tofd,                                             &
          do_gwd_opt_psl, psl_gwd_dx_factor,                                            &
          do_ngw_ec, do_ugwp_v1,  do_ugwp_v1_orog_only,                                 &
          do_ugwp_v1_w_gsldrag, gwd_opt, do_tofd, ldiag_ugwp, ugwp_seq_update,          &
          cdmbgwd, alpha_fd, jdat, nmtvr, hprime, oc, theta, sigma, gamma,              &
          elvmax, clx, oa4, varss,oc1ss,oa4ss,ol4ss, dx,  xlat, xlat_d, sinlat, coslat, &
          area, rain, br1, hpbl,vtype, kpbl, slmsk,                                     &
          ugrs, vgrs, tgrs, q1, prsi, prsl, prslk, phii, phil,  del, tau_amf,           &
          dudt_ogw, dvdt_ogw, du_ogwcol, dv_ogwcol,                                     &
          dudt_obl, dvdt_obl, du_oblcol, dv_oblcol,                                     &
          dudt_oss, dvdt_oss, du_osscol, dv_osscol,                                     &
          dudt_ofd, dvdt_ofd, du_ofdcol, dv_ofdcol,                                     &
          dudt_ngw, dvdt_ngw, dtdt_ngw, kdis_ngw, dudt_gw, dvdt_gw, dtdt_gw, dqdt_gw,   &
          kdis_gw, tau_ogw, tau_ngw,  tau_oss,                                          &
          zogw,  zlwb,  zobl,  zngw,   dusfcg, dvsfcg, rdxzb,                           &
          dtend, dtidx, index_of_x_wind, index_of_y_wind, index_of_temperature,         &
          index_of_process_orographic_gwd, index_of_process_nonorographic_gwd,          &
          lprnt, ipr, spp_wts_gwd, spp_gwd, errmsg, errflg)

!
!########################################################################
!  Attention New Arrays and Names must be ADDED inside
!
! a) /FV3/gfsphysics/GFS_layer/GFS_typedefs.meta
! b) /FV3/gfsphysics/GFS_layer/GFS_typedefs.F90
! c) /FV3/gfsphysics/GFS_layer/GFS_diagnostics.F90 "diag-cs is not tested"
!########################################################################

    use ugwp_common, only : con_pi => pi, con_g => grav,  con_rd   => rd,   &
                            con_rv => rv, con_cp => cpd,  con_fv   => fv,   &
                            con_rerth => arad, con_omega => omega1, rgrav

    implicit none

! Preference use    (im,levs) rather than (:,:) to avoid memory-leaks
!                    that found in Nov-Dec 2020
! order array-description control-logical
!                         other in-variables
!                         out-variables
!                         local-variables
!
!  unified GSL and CIRES diagnostics inside CCPP and GFS_typedefs.F90/GFS_diagnostics.F90
!
! interface variables
    logical,                 intent(in) :: ldiag3d, lssav
    logical,                 intent(in) :: flag_for_gwd_generic_tend
    logical,                 intent(in) :: lprnt

    integer,                 intent(in) :: ipr

! flags for choosing combination of GW drag schemes to run

    logical,  intent (in) :: do_gsl_drag_ls_bl, do_gsl_drag_ss, do_gsl_drag_tofd
    logical,  intent (in) :: do_ugwp_v1, do_ngw_ec, do_ugwp_v1_orog_only, do_tofd
    logical,  intent (in) :: ldiag_ugwp, ugwp_seq_update
    logical,  intent (in) :: do_ugwp_v1_w_gsldrag                              ! combination of ORO and NGW schemes

    integer,                 intent(in) :: me, master, im, levs, ntrac,lonr
    real(kind=kind_phys),    intent(in) :: dtp
    real(kind=kind_phys),    intent(in) :: ak(:), bk(:)
    integer,                 intent(in) :: kdt, jdat(:)
! option  for psl gwd
    logical, intent(in)              :: do_gwd_opt_psl      ! option for psl gravity wave drag
    real(kind=kind_phys), intent(in) :: psl_gwd_dx_factor   !
! SSO parameters and variables
    integer,                 intent(in) :: gwd_opt                         !gwd_opt  and nmtvr are "redundant" controls
    integer,                 intent(in) :: nmtvr
    real(kind=kind_phys),    intent(in) :: cdmbgwd(:), alpha_fd            ! for gsl_drag

    real(kind=kind_phys),    intent(in), dimension(:)       :: hprime, oc, theta, sigma, gamma

    real(kind=kind_phys),    intent(in), dimension(:)       :: elvmax
    real(kind=kind_phys),    intent(in), dimension(:,:)     :: clx, oa4

    real(kind=kind_phys),    intent(in), dimension(:)       :: dx
    real(kind=kind_phys),    intent(in), dimension(:)       :: varss,oc1ss
    real(kind=kind_phys),    intent(in), dimension(:,:)     :: oa4ss,ol4ss

!=====
!ccpp-style passing constants, I prefer to take them out from the "call-subr" list
!=====
!    real(kind=kind_phys),    intent(in) :: con_g, con_omega, con_pi, con_cp, con_rd, &
!                                           con_rv, con_rerth, con_fvirt
! grids

    real(kind=kind_phys),    intent(in), dimension(:)     :: xlat, xlat_d, sinlat, coslat, area

! State vars + PBL/slmsk +rain

    real(kind=kind_phys),    intent(in), dimension(:,:)   :: del, ugrs, vgrs, tgrs, prsl, prslk, phil
    real(kind=kind_phys),    intent(in), dimension(:,:)   :: prsi, phii
    real(kind=kind_phys),    intent(in), dimension(:,:)   :: q1
    integer,                 intent(in), dimension(:)     :: kpbl
    integer,                 intent(in), dimension(:)     :: vtype

    real(kind=kind_phys),    intent(in), dimension(:) :: rain
    real(kind=kind_phys),    intent(in), dimension(:) :: br1, hpbl,  slmsk
!
! moved to GFS_phys_time_vary
!    real(kind=kind_phys),    intent(in), dimension(:) :: ddy_j1tau, ddy_j2tau
!    integer,                 intent(in), dimension(:) :: jindx1_tau, jindx2_tau
     real(kind=kind_phys),    intent(in), dimension(:) :: tau_amf

!Output (optional):

    real(kind=kind_phys), intent(out), dimension(:)  ::         &
                            du_ogwcol,  dv_ogwcol,  du_oblcol, dv_oblcol, &
                            du_osscol,  dv_osscol,  du_ofdcol, dv_ofdcol
!
! we may add later but due to launch in the upper layes ~ mPa comparing to ORO Pa*(0.1)
!                            du_ngwcol, dv_ngwcol

    real(kind=kind_phys), intent(out), dimension(:)  :: dusfcg, dvsfcg
    real(kind=kind_phys), intent(out), dimension(:)  :: tau_ogw, tau_ngw, tau_oss

    real(kind=kind_phys), intent(out) , dimension(:,:) ::         &
                          dudt_ogw, dvdt_ogw, dudt_obl, dvdt_obl, &
                          dudt_oss, dvdt_oss, dudt_ofd, dvdt_ofd

    real(kind=kind_phys), intent(out) , dimension(:,:) :: dudt_ngw, dvdt_ngw, kdis_ngw, dtdt_ngw
    real(kind=kind_phys), intent(out) , dimension(:,:) :: dudt_gw,  dvdt_gw, dtdt_gw, kdis_gw
    real(kind=kind_phys), intent(out) , dimension(:,:,:) :: dqdt_gw

    real(kind=kind_phys), intent(out) , dimension(:)   :: zogw, zlwb, zobl, zngw
!
    real(kind=kind_phys), intent(inout), optional            :: dtend(:,:,:)
    integer, intent(in)                                      :: dtidx(:,:)
    integer, intent(in)                                 :: &
         index_of_x_wind, index_of_y_wind, index_of_temperature,              &
         index_of_process_orographic_gwd, index_of_process_nonorographic_gwd

    real(kind=kind_phys),    intent(out), dimension(:)      :: rdxzb     ! for stoch phys. mtb-level

    real(kind=kind_phys), intent(in), optional :: spp_wts_gwd(:,:)
    integer, intent(in) :: spp_gwd

    character(len=*),        intent(out) :: errmsg
    integer,                 intent(out) :: errflg

! Working-precision (_loc) buffers for CCPP host boundary isolation
    real(kind=wp) :: dtp_loc, psl_gwd_dx_factor_loc, alpha_fd_loc
    real(kind=wp), dimension(2)           :: cdmbgwd_loc
    real(kind=wp), dimension(im)          :: hprime_loc, oc_loc, theta_loc, sigma_loc, gamma_loc
    real(kind=wp), dimension(im)          :: elvmax_loc, dx_loc, varss_loc, oc1ss_loc
    real(kind=wp), dimension(im, 4)       :: clx_loc, oa4_loc, oa4ss_loc, ol4ss_loc
    real(kind=wp), dimension(im)          :: xlat_d_loc, sinlat_loc, coslat_loc, area_loc
    real(kind=wp), dimension(im)          :: rain_loc, br1_loc, hpbl_loc, slmsk_loc, tau_amf_loc
    real(kind=wp), dimension(im, levs)    :: del_loc, ugrs_loc, vgrs_loc, tgrs_loc, prsl_loc, prslk_loc, phil_loc
    real(kind=wp), dimension(im, levs+1)  :: prsi_loc, phii_loc
    real(kind=wp), dimension(im, levs)    :: q1_loc

    real(kind=wp), dimension(im)          :: du_ogwcol_loc, dv_ogwcol_loc, du_oblcol_loc, dv_oblcol_loc
    real(kind=wp), dimension(im)          :: du_osscol_loc, dv_osscol_loc, du_ofdcol_loc, dv_ofdcol_loc
    real(kind=wp), dimension(im)          :: dusfcg_loc, dvsfcg_loc
    real(kind=wp), dimension(im)          :: tau_ogw_loc, tau_ngw_loc, tau_oss_loc
    real(kind=wp), dimension(im, levs)    :: dudt_ogw_loc, dvdt_ogw_loc, dudt_obl_loc, dvdt_obl_loc
    real(kind=wp), dimension(im, levs)    :: dudt_oss_loc, dvdt_oss_loc, dudt_ofd_loc, dvdt_ofd_loc
    real(kind=wp), dimension(im, levs)    :: dudt_ngw_loc, dvdt_ngw_loc, kdis_ngw_loc, dtdt_ngw_loc
    real(kind=wp), dimension(im, levs)    :: dudt_gw_loc,  dvdt_gw_loc,  dtdt_gw_loc,  kdis_gw_loc
    real(kind=wp), dimension(im, levs, ntrac) :: dqdt_gw_loc
    real(kind=wp), dimension(im)          :: zogw_loc, zlwb_loc, zobl_loc, zngw_loc, rdxzb_loc

    real(kind=wp), allocatable            :: dtend_loc(:,:,:)
    real(kind=wp), allocatable            :: spp_wts_gwd_loc(:,:)

! local variables
    integer :: i, k
    real(kind=wp), dimension(im)          :: sgh30
    real(kind=wp), dimension(im, levs)    :: Pdvdt, Pdudt
    real(kind=wp), dimension(im, levs)    :: Pdtdt, Pkdis
!------------
!
! from ugwp_driver_v0.f -> cires_ugwp_initialize.F90 -> module ugwp_wmsdis_init
!  now in the namelist of cires_ugwp "knob_ugwp_tauamp" controls tamp_mpa
!
!        tamp_mpa =knob_ugwp_tauamp                         !amplitude for GEOS-5/MERRA-2
!------------
!    real(kind=kind_phys), parameter :: tamp_mpa_v0=30.e-3  ! large flux to help "GFS-ensembles" in July 2019

! switches that activate impact of OGWs and NGWs

!    integer :: nmtvr_temp

    real(kind=wp), dimension(im, levs)    :: zmet  ! geopotential height at model Layer centers
    real(kind=wp), dimension(im, levs+1)  :: zmeti ! geopotential height at model layer interfaces

! ugwp_v1 local variables

    integer :: y4, month, day,  ddd_ugwp, curdate, curday, idtend

!  ugwp_v1 temporary (local) diagnostic variables from cires_ugwp_solv2_v1
!  diagnostics for wind and temp rms to compare with space-borne data and metrics
!   in the Middle atmosphere: 20-110 km ( not active in CCPP-style, oct 2020)
!    real(kind=kind_phys) :: tauabs(im,levs), wrms(im,levs), trms(im,levs)

    ! Initialize CCPP error handling variables

    errmsg = ''
    errflg = 0

    ! Cast CCPP kind_phys inputs to internal working precision (wp) buffers
    dtp_loc               = real(dtp,               kind=wp)
    psl_gwd_dx_factor_loc = real(psl_gwd_dx_factor, kind=wp)
    alpha_fd_loc          = real(alpha_fd,          kind=wp)
    cdmbgwd_loc(1:2)      = real(cdmbgwd(1:2),      kind=wp)

    hprime_loc(1:im)      = real(hprime(1:im),      kind=wp)
    oc_loc(1:im)          = real(oc(1:im),          kind=wp)
    theta_loc(1:im)       = real(theta(1:im),       kind=wp)
    sigma_loc(1:im)       = real(sigma(1:im),       kind=wp)
    gamma_loc(1:im)       = real(gamma(1:im),       kind=wp)
    elvmax_loc(1:im)      = real(elvmax(1:im),      kind=wp)
    dx_loc(1:im)          = real(dx(1:im),          kind=wp)
    varss_loc(1:im)       = real(varss(1:im),       kind=wp)
    oc1ss_loc(1:im)       = real(oc1ss(1:im),       kind=wp)

    clx_loc(1:im, 1:4)    = real(clx(1:im, 1:4),    kind=wp)
    oa4_loc(1:im, 1:4)    = real(oa4(1:im, 1:4),    kind=wp)
    oa4ss_loc(1:im, 1:4)  = real(oa4ss(1:im, 1:4),  kind=wp)
    ol4ss_loc(1:im, 1:4)  = real(ol4ss(1:im, 1:4),  kind=wp)

    xlat_d_loc(1:im)      = real(xlat_d(1:im),      kind=wp)
    sinlat_loc(1:im)      = real(sinlat(1:im),      kind=wp)
    coslat_loc(1:im)      = real(coslat(1:im),      kind=wp)
    area_loc(1:im)        = real(area(1:im),        kind=wp)

    rain_loc(1:im)        = real(rain(1:im),        kind=wp)
    br1_loc(1:im)         = real(br1(1:im),         kind=wp)
    hpbl_loc(1:im)        = real(hpbl(1:im),        kind=wp)
    slmsk_loc(1:im)       = real(slmsk(1:im),       kind=wp)
    tau_amf_loc(1:im)     = real(tau_amf(1:im),     kind=wp)

    del_loc(1:im, 1:levs)   = real(del(1:im, 1:levs),   kind=wp)
    ugrs_loc(1:im, 1:levs)  = real(ugrs(1:im, 1:levs),  kind=wp)
    vgrs_loc(1:im, 1:levs)  = real(vgrs(1:im, 1:levs),  kind=wp)
    tgrs_loc(1:im, 1:levs)  = real(tgrs(1:im, 1:levs),  kind=wp)
    prsl_loc(1:im, 1:levs)  = real(prsl(1:im, 1:levs),  kind=wp)
    prslk_loc(1:im, 1:levs) = real(prslk(1:im, 1:levs), kind=wp)
    phil_loc(1:im, 1:levs)  = real(phil(1:im, 1:levs),  kind=wp)
    q1_loc(1:im, 1:levs)    = real(q1(1:im, 1:levs),    kind=wp)

    prsi_loc(1:im, 1:levs+1) = real(prsi(1:im, 1:levs+1), kind=wp)
    phii_loc(1:im, 1:levs+1) = real(phii(1:im, 1:levs+1), kind=wp)

    if (present(dtend)) then
       allocate(dtend_loc(size(dtend,1), size(dtend,2), size(dtend,3)))
       dtend_loc = real(dtend, kind=wp)
    else
       allocate(dtend_loc(1, 1, 1))
       dtend_loc = 0.0_wp
    endif

    if (present(spp_wts_gwd)) then
       allocate(spp_wts_gwd_loc(size(spp_wts_gwd,1), size(spp_wts_gwd,2)))
       spp_wts_gwd_loc = real(spp_wts_gwd, kind=wp)
    else
       allocate(spp_wts_gwd_loc(im, 1))
       spp_wts_gwd_loc = 0.0_wp
    endif

! 1) ORO stationary GWs
!    ------------------
!
! for all oro-suites can uze geo-meters having "hpbl"
!
! All GW-schemes operate with Zmet =phil*inv_g, passing Zmet/Zmeti can be more robust
! + rho*dz = =delp *  inv_g   can be also pre-comp for all "GW-schemes"
!
       zmeti  = phii_loc* rgrav
       zmet   = phil_loc* rgrav

!===============================================================
! ORO-diag

       if (do_ugwp_v1 .or. ldiag_ugwp) then
         dudt_ogw_loc(:,:)= 0.0_wp; dvdt_ogw_loc(:,:)=0.0_wp; dudt_obl_loc(:,:)=0.0_wp; dvdt_obl_loc(:,:)=0.0_wp
         dudt_oss_loc(:,:)= 0.0_wp; dvdt_oss_loc(:,:)=0.0_wp; dudt_ofd_loc(:,:)=0.0_wp; dvdt_ofd_loc(:,:)=0.0_wp
         du_ogwcol_loc(:)=0.0_wp  ; dv_ogwcol_loc(:)=0.0_wp ; du_oblcol_loc(:)=0.0_wp ; dv_oblcol_loc(:)=0.0_wp
         du_osscol_loc(:)=0.0_wp  ; dv_osscol_loc(:)=0.0_wp ; du_ofdcol_loc(:)=0.0_wp ; dv_ofdcol_loc(:)=0.0_wp
         dudt_ngw_loc(:,:)=0.0_wp; dvdt_ngw_loc(:,:)=0.0_wp; dtdt_ngw_loc(:,:)=0.0_wp; kdis_ngw_loc(:,:)=0.0_wp
       else
         dudt_ogw_loc(:,:)  = 0.0_wp
       end if

       dusfcg_loc (:)  = 0.0_wp  ;  dvsfcg_loc(:) = 0.0_wp

! ngw+ogw - diag

       dudt_gw_loc(:,:)=0.0_wp ;  dvdt_gw_loc(:,:)=0.0_wp  ; dtdt_gw_loc(:,:)=0.0_wp  ; dqdt_gw_loc(:,:,:)=0.0_wp ; kdis_gw_loc(:,:)=0.0_wp
! source fluxes

      tau_ogw_loc(:)=0.0_wp ; tau_ngw_loc(:)=0.0_wp ;  tau_oss_loc(:)=0.0_wp

! launch layers

      zlwb_loc(:)= 0.0_wp  ; zogw_loc(:)=0.0_wp ;  zobl_loc(:)=0.0_wp ;  zngw_loc(:)=0.0_wp
      rdxzb_loc(:)= 0.0_wp
!===============================================================
!  diag tendencies due to all-SSO schemes (ORO-physics)
!  ogw + obl + oss + ofd ..... no explicit "lee wave trapping"
!===============================================================
     do k=1,levs
        do i=1,im
          Pdvdt(i,k) = 0.0_wp
          Pdudt(i,k) = 0.0_wp
          Pdtdt(i,k) = 0.0_wp
          Pkdis(i,k) = 0.0_wp
        enddo
      enddo

    ! Run the appropriate large-scale (large-scale GWD + blocking) scheme
    ! Note:  In case of GSL drag_suite, this includes ss and tofd

    if ( do_gsl_drag_ls_bl.or.do_gsl_drag_ss.or.do_gsl_drag_tofd) then
!
! to do: the zero diag and tendency values assigned inside "drag_suite_run" can be skipped :
!
! dudt_ogw, dvdt_ogw, dudt_obl, dvdt_obl,dudt_oss, dvdt_oss, dudt_ofd, dvdt_ofd
! du_ogwcol, dv_ogwcol, du_oblcol, dv_oblcol, du_osscol, dv_osscol, du_ofdcol dv_ofdcol
! dusfcg,  dvsfcg
!
     if (do_gwd_opt_psl) then
       call drag_suite_ugwpv1_psl(im, levs, Pdvdt, Pdudt, Pdtdt,     &
                 ugrs_loc,vgrs_loc,tgrs_loc,q1_loc,                  &
                 kpbl,prsi_loc,del_loc,prsl_loc,prslk_loc,           &
                 phii_loc,phil_loc,dtp_loc,                          &
                 kdt,hprime_loc,oc_loc,oa4_loc,clx_loc,              &
                 varss_loc,oc1ss_loc,oa4ss_loc,                      &
                 ol4ss_loc,theta_loc,sigma_loc,gamma_loc,elvmax_loc, &
                 dudt_ogw_loc, dvdt_ogw_loc, dudt_obl_loc, dvdt_obl_loc, &
                 dudt_oss_loc, dvdt_oss_loc, dudt_ofd_loc, dvdt_ofd_loc, &
                 dusfcg_loc,  dvsfcg_loc,                            &
                 du_ogwcol_loc, dv_ogwcol_loc, du_oblcol_loc, dv_oblcol_loc, &
                 du_osscol_loc, dv_osscol_loc, du_ofdcol_loc, dv_ofdcol_loc, &
                 slmsk_loc,br1_loc,hpbl_loc,vtype,                   &
                 con_g,con_cp,con_rd,con_rv,                         &
                 con_fv, con_pi, lonr,                               &
                 cdmbgwd_loc(1:2),alpha_fd_loc,me,master,            &
                 lprnt,ipr,rdxzb_loc,dx_loc,gwd_opt,                 &
                 do_gsl_drag_ls_bl,do_gsl_drag_ss,do_gsl_drag_tofd,  &
                 psl_gwd_dx_factor_loc, flag_for_gwd_generic_tend,   &
                 dtend_loc, dtidx, index_of_process_orographic_gwd,  &
                 index_of_temperature, index_of_x_wind,              &
                 index_of_y_wind, ldiag3d, ldiag_ugwp,               &
                 ugwp_seq_update, spp_wts_gwd_loc, spp_gwd, errmsg, errflg)
     else
       call drag_suite_ugwpv1_run(im, levs, Pdvdt, Pdudt, Pdtdt, dqdt_gw_loc, &
                 ugrs_loc,vgrs_loc,tgrs_loc,q1_loc,                  &
                 kpbl,prsi_loc,del_loc,prsl_loc,prslk_loc,           &
                 phii_loc,phil_loc,dtp_loc,                          &
                 kdt,hprime_loc,oc_loc,oa4_loc,clx_loc,              &
                 varss_loc,oc1ss_loc,oa4ss_loc,                      &
                 ol4ss_loc,theta_loc,sigma_loc,gamma_loc,elvmax_loc, &
                 dudt_ogw_loc, dvdt_ogw_loc, dudt_obl_loc, dvdt_obl_loc, &
                 dudt_oss_loc, dvdt_oss_loc, dudt_ofd_loc, dvdt_ofd_loc, &
                 dusfcg_loc,  dvsfcg_loc,                            &
                 du_ogwcol_loc, dv_ogwcol_loc, du_oblcol_loc, dv_oblcol_loc, &
                 du_osscol_loc, dv_osscol_loc, du_ofdcol_loc, dv_ofdcol_loc, &
                 slmsk_loc,br1_loc,hpbl_loc,con_g,con_cp,con_rd,con_rv, &
                 con_fv, con_pi, lonr,                               &
                 cdmbgwd_loc(1:2),alpha_fd_loc,me,master,            &
                 lprnt,ipr,rdxzb_loc,dx_loc,gwd_opt,                 &
                 do_gsl_drag_ls_bl,do_gsl_drag_ss,do_gsl_drag_tofd,  &
                 flag_for_gwd_generic_tend,                          &
                 dtend_loc, dtidx, index_of_process_orographic_gwd,  &
                 index_of_temperature, index_of_x_wind,              &
                 index_of_y_wind, ldiag3d, ldiag_ugwp,               &
                 ugwp_seq_update, spp_wts_gwd_loc, spp_gwd, errmsg, errflg)
     endif
     if(errflg/=0) then
        if (allocated(dtend_loc))       deallocate(dtend_loc)
        if (allocated(spp_wts_gwd_loc)) deallocate(spp_wts_gwd_loc)
        return
     endif

!
! dusfcg = du_ogwcol + du_oblcol + du_osscol + du_ofdcol
!
!         if (kdt <= 2 .and. me == master) then
!      print *, ' unified drag_suite_run ', kdt
!      print *, ' GSL drag du/dt ', maxval(Pdudt)*86400, minval(Pdudt)*86400
!      print *, ' GSL drag dv/dt ', maxval(Pdvdt)*86400, minval(Pdvdt)*86400
!
! zero      print *, ' unified drag_GSL dT/dt ', maxval(Pdtdt)*86400, minval(Pdtdt)*86400
!
!      if (gwd_opt == 22 .or. gwd_opt == 33) then
!      print *, ' unified drag_GSL dUBL/dt ',  maxval(dudt_obl)*86400, minval(dudt_obl)*86400
!      print *, ' unified drag_GSL dVBL/dt ',  maxval(dvdt_obl)*86400, minval(dvdt_obl)*86400
!      print *, ' unified drag_GSL dUOGW/dt ', maxval(dudt_ogw)*86400, minval(dudt_ogw)*86400
!      print *, ' unified drag_GSL dVOGW/dt ', maxval(dvdt_ogw)*86400, minval(dvdt_ogw)*86400
!      print *, ' unified drag_GSL dUOss/dt ', maxval(dudt_oss)*86400, minval(dudt_oss)*86400
!      print *, ' unified drag_GSL dVOSS/dt ', maxval(dvdt_oss)*86400, minval(dvdt_oss)*86400
!      print *, ' unified drag_GSL dUOfd/dt ', maxval(dudt_ofd)*86400, minval(dudt_ofd)*86400
!      print *, ' unified drag_GSL dVOfd/dt ', maxval(dvdt_ofd)*86400, minval(dvdt_ofd)*86400
!      endif
!     endif

    endif
!
! not gsldrag large-scale oro-scheme for example "do_ugwp_v1_orog_only"
!

    if ( do_ugwp_v1_orog_only ) then
!
! for TOFD we use now "varss" of GSL-drag  /not sgh30=abs(oro-oro_f)/
! only sum of integrated ORO+GW effects (dusfcg and dvsfcg) = sum(ogw + obl + oss*0 + ofd + ngw)
!
! OROGW_V1 introduce "orchestration" between OGW-effects and Mountain Blocking
!      it starts to examines options for the Scale-Aware (SA)formulation of SSO-effects
!      if ( me == master .and. kdt == 1) print *, ' bf orogw_v1 nmtvr=', nmtvr, ' do_tofd=', do_tofd

         if (gwd_opt ==1 )sgh30 = 0.15_wp*hprime_loc    ! portion of the mesoscale SSO (~[oro_unfilt -oro_filt)
         if (gwd_opt >1 ) sgh30 = varss_loc             ! as in gsldrag: see drag_suite_run

       call orogw_v1 (im, levs,  lonr,  me, master,dtp_loc, kdt, do_tofd, &
                      xlat_d_loc, sinlat_loc, coslat_loc, area_loc,       &
                      cdmbgwd_loc(1:2), hprime_loc, oc_loc, oa4_loc, clx_loc, theta_loc, &
                      sigma_loc, gamma_loc, elvmax_loc,  sgh30,  kpbl, ugrs_loc, &
                      vgrs_loc, tgrs_loc, q1_loc, prsi_loc,del_loc,prsl_loc,prslk_loc, zmeti, zmet, &
                      Pdvdt, Pdudt, Pdtdt, Pkdis, dusfcg_loc, dvsfcg_loc, rdxzb_loc, &
                      zobl_loc, zlwb_loc, zogw_loc, tau_ogw_loc, dudt_ogw_loc, dvdt_ogw_loc, &
                      dudt_obl_loc, dvdt_obl_loc, dudt_ofd_loc, dvdt_ofd_loc, &
                      du_ogwcol_loc, dv_ogwcol_loc, du_oblcol_loc, dv_oblcol_loc, &
                      du_ofdcol_loc, dv_ofdcol_loc, errmsg, errflg )
       if(errflg/=0) then
          if (allocated(dtend_loc))       deallocate(dtend_loc)
          if (allocated(spp_wts_gwd_loc)) deallocate(spp_wts_gwd_loc)
          return
       endif

!
! orogw_v1: dusfcg = du_ogwcol + du_oblcol  + du_ofdcol                           only 3 terms
!
!          if (kdt <= 2 .and. me == master) then
!
!       print *, ' unified_ugwp orogw_v1 ', kdt, me,  nmtvr
!       print *, ' unified_ugwp orogw_v1 du/dt ', maxval(Pdudt)*86400, minval(Pdudt)*86400
!       print *, ' unified_ugwp orogw_v1 dv/dt ', maxval(Pdvdt)*86400, minval(Pdvdt)*86400
!       print *, ' unified_ugwp orogw_v1 dT/dt ', maxval(Pdtdt)*86400, minval(Pdtdt)*86400
!       print *, ' unified_ugwp orogw_v1 dUBL/dt ', maxval(dudt_obl)*86400, minval(dudt_obl)*86400
!       print *, ' unified_ugwp orogw_v1 dVBL/dt ', maxval(dvdt_obl)*86400, minval(dvdt_obl)*86400
!      endif

    end if
!
!  for  old-fashioned GFS-style diag-cs like dt3dt(:.:, 1:14) collections
!
    if(ldiag3d .and. lssav .and. .not. flag_for_gwd_generic_tend) then
      idtend = dtidx(index_of_x_wind,index_of_process_orographic_gwd)
      if(idtend>=1) then
         dtend_loc(:,:,idtend) = dtend_loc(:,:,idtend) + Pdudt*dtp_loc
      endif
      idtend = dtidx(index_of_y_wind,index_of_process_orographic_gwd)
      if(idtend>=1) then
         dtend_loc(:,:,idtend) = dtend_loc(:,:,idtend) + Pdvdt*dtp_loc
      endif
      idtend = dtidx(index_of_temperature,index_of_process_orographic_gwd)
      if(idtend>=1) then
         dtend_loc(:,:,idtend) = dtend_loc(:,:,idtend) + Pdtdt*dtp_loc
      endif
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Begin non-stationary GW schemes
! ugwp_v1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (do_ugwp_v1) then

!==================================================================
!       call slat_geos5_tamp_v1(im, tamp_mpa, xlat_d, tau_ngw)
!
! 2020 updates of MERRA/GEOS tau_ngw for the C96-QBO FV3GFS-127L runs
!==================================================================

       call  slat_geos5_2020(im, tamp_mpa, xlat_d_loc, tau_ngw_loc)

       y4 = jdat(1); month = jdat(2); day = jdat(3)
!
! hour = jdat(5)
! fhour = float(hour)+float(jdat(6))/60. + float(jdat(7))/3600.
!       fhour = (kdt-1)*dtp/3600.
!       fhrday  = fhour/24.  - nint(fhour/24.)

       call calendar_ugwp(y4, month, day, ddd_ugwp)
       curdate = y4*1000 + ddd_ugwp

       call ngwflux_update(me, master, im, levs, kdt, ddd_ugwp,curdate, &
         tau_amf_loc, xlat_d_loc, sinlat_loc, coslat_loc, rain_loc, tau_ngw_loc)

       if (do_ngw_ec) then

       call ecmwf_ngw_emc(me, master, im,   levs,  kdt, dtp_loc, dx_loc, &
                      tau_ngw_loc, tgrs_loc, ugrs_loc,  vgrs_loc,   q1_loc, prsl_loc, prsi_loc, &
                      zmet, zmeti, prslk_loc, xlat_d_loc, sinlat_loc, coslat_loc, &
                      dudt_ngw_loc, dvdt_ngw_loc, dtdt_ngw_loc, kdis_ngw_loc, zngw_loc)
       else

       call cires_ugwpv1_ngw_solv2(me, master, im, levs, kdt, dtp_loc,   &
                      tau_ngw_loc, tgrs_loc, ugrs_loc,  vgrs_loc,   q1_loc, prsl_loc, prsi_loc, &
                      zmet, zmeti, prslk_loc, xlat_d_loc, sinlat_loc, coslat_loc, &
                      dudt_ngw_loc, dvdt_ngw_loc, dtdt_ngw_loc, kdis_ngw_loc, zngw_loc)

      endif

!
! =>  con_g, con_cp, con_rd, con_rv, con_omega,  con_pi, con_fvirt
!
!       if (me == master .and. kdt <= 2) then
!         print *
!         write(6,*)'FV3GFS finished fv3_ugwp_solv2_v1   '
!         write(6,*) ' non-stationary GWs with GMAO/MERRA GW-forcing '
!         print *
!
!      print *, ' ugwp_v1 ', kdt
!      print *, ' ugwp_v1 du/dt ', maxval(dudt_ngw)*86400, minval(dudt_ngw)*86400
!      print *, ' ugwp_v1 dv/dt ', maxval(dvdt_ngw)*86400, minval(dvdt_ngw)*86400
!      print *, ' ugwp_v1 dT/dt ', maxval(dtdt_ngw)*86400, minval(dtdt_ngw)*86400
!       endif

    end if   ! do_ugwp_v1

    if(ldiag3d .and. lssav .and. .not. flag_for_gwd_generic_tend) then
      idtend = dtidx(index_of_x_wind,index_of_process_nonorographic_gwd)
      if(idtend>=1) then
         dtend_loc(:,:,idtend) = dtend_loc(:,:,idtend) + dudt_ngw_loc(:,:)*dtp_loc
      endif
      idtend = dtidx(index_of_y_wind,index_of_process_nonorographic_gwd)
      if(idtend>=1) then
         dtend_loc(:,:,idtend) = dtend_loc(:,:,idtend) + dvdt_ngw_loc(:,:)*dtp_loc
      endif
      idtend = dtidx(index_of_temperature,index_of_process_nonorographic_gwd)
      if(idtend>=1) then
         dtend_loc(:,:,idtend) = dtend_loc(:,:,idtend) + dtdt_ngw_loc(:,:)*dtp_loc
      endif
    endif

!
! get total sso-OGW + NGW
!
     if (do_ugwp_v1) then
        dudt_gw_loc =  Pdudt + dudt_ngw_loc
        dvdt_gw_loc =  Pdvdt + dvdt_ngw_loc
        dtdt_gw_loc =  Pdtdt + dtdt_ngw_loc
        kdis_gw_loc =  Pkdis + kdis_ngw_loc
     else
        dudt_gw_loc =  Pdudt
        dvdt_gw_loc =  Pdvdt
        dtdt_gw_loc =  Pdtdt
        kdis_gw_loc =  Pkdis
     end if

    ! Restore working-precision (_loc) outputs back to CCPP host kind_phys
    du_ogwcol(1:im) = real(du_ogwcol_loc(1:im), kind=kind_phys)
    dv_ogwcol(1:im) = real(dv_ogwcol_loc(1:im), kind=kind_phys)
    du_oblcol(1:im) = real(du_oblcol_loc(1:im), kind=kind_phys)
    dv_oblcol(1:im) = real(dv_oblcol_loc(1:im), kind=kind_phys)
    du_osscol(1:im) = real(du_osscol_loc(1:im), kind=kind_phys)
    dv_osscol(1:im) = real(dv_osscol_loc(1:im), kind=kind_phys)
    du_ofdcol(1:im) = real(du_ofdcol_loc(1:im), kind=kind_phys)
    dv_ofdcol(1:im) = real(dv_ofdcol_loc(1:im), kind=kind_phys)

    dusfcg(1:im)    = real(dusfcg_loc(1:im),    kind=kind_phys)
    dvsfcg(1:im)    = real(dvsfcg_loc(1:im),    kind=kind_phys)
    tau_ogw(1:im)   = real(tau_ogw_loc(1:im),   kind=kind_phys)
    tau_ngw(1:im)   = real(tau_ngw_loc(1:im),   kind=kind_phys)
    tau_oss(1:im)   = real(tau_oss_loc(1:im),   kind=kind_phys)

    dudt_ogw(1:im, 1:levs) = real(dudt_ogw_loc(1:im, 1:levs), kind=kind_phys)
    dvdt_ogw(1:im, 1:levs) = real(dvdt_ogw_loc(1:im, 1:levs), kind=kind_phys)
    dudt_obl(1:im, 1:levs) = real(dudt_obl_loc(1:im, 1:levs), kind=kind_phys)
    dvdt_obl(1:im, 1:levs) = real(dvdt_obl_loc(1:im, 1:levs), kind=kind_phys)
    dudt_oss(1:im, 1:levs) = real(dudt_oss_loc(1:im, 1:levs), kind=kind_phys)
    dvdt_oss(1:im, 1:levs) = real(dvdt_oss_loc(1:im, 1:levs), kind=kind_phys)
    dudt_ofd(1:im, 1:levs) = real(dudt_ofd_loc(1:im, 1:levs), kind=kind_phys)
    dvdt_ofd(1:im, 1:levs) = real(dvdt_ofd_loc(1:im, 1:levs), kind=kind_phys)

    dudt_ngw(1:im, 1:levs) = real(dudt_ngw_loc(1:im, 1:levs), kind=kind_phys)
    dvdt_ngw(1:im, 1:levs) = real(dvdt_ngw_loc(1:im, 1:levs), kind=kind_phys)
    dtdt_ngw(1:im, 1:levs) = real(dtdt_ngw_loc(1:im, 1:levs), kind=kind_phys)
    kdis_ngw(1:im, 1:levs) = real(kdis_ngw_loc(1:im, 1:levs), kind=kind_phys)

    dudt_gw(1:im, 1:levs)  = real(dudt_gw_loc(1:im, 1:levs),  kind=kind_phys)
    dvdt_gw(1:im, 1:levs)  = real(dvdt_gw_loc(1:im, 1:levs),  kind=kind_phys)
    dtdt_gw(1:im, 1:levs)  = real(dtdt_gw_loc(1:im, 1:levs),  kind=kind_phys)
    kdis_gw(1:im, 1:levs)  = real(kdis_gw_loc(1:im, 1:levs),  kind=kind_phys)
    dqdt_gw(1:im, 1:levs, 1:ntrac) = real(dqdt_gw_loc(1:im, 1:levs, 1:ntrac), kind=kind_phys)

    zogw(1:im)      = real(zogw_loc(1:im),      kind=kind_phys)
    zlwb(1:im)      = real(zlwb_loc(1:im),      kind=kind_phys)
    zobl(1:im)      = real(zobl_loc(1:im),      kind=kind_phys)
    zngw(1:im)      = real(zngw_loc(1:im),      kind=kind_phys)
    rdxzb(1:im)     = real(rdxzb_loc(1:im),     kind=kind_phys)

    if (present(dtend)) then
       dtend = real(dtend_loc, kind=kind_phys)
    endif

    if (allocated(dtend_loc))       deallocate(dtend_loc)
    if (allocated(spp_wts_gwd_loc)) deallocate(spp_wts_gwd_loc)

    end subroutine ugwpv1_gsldrag_run
end module ugwpv1_gsldrag
