!>  \file ugwpv1_gsldrag_post.F90

!> This module contains code to be executed after the UGWP v1 scheme
module ugwpv1_gsldrag_post

contains

!>\defgroup ugwpv1_gsldrag_post ugwpv1_gsldrag Scheme Post
!> \section arg_table_ugwpv1_gsldrag_post_run Argument Table
!! \htmlinclude ugwpv1_gsldrag_post_run.html
!!
     subroutine ugwpv1_gsldrag_post_run ( im, levs, ldiag_ugwp,       &
         dtf, dudt_gw, dvdt_gw,                                       &
         tau_ogw, tau_ngw, zobl, zlwb, zogw, dudt_obl, dvdt_obl,      &
         dudt_ofd, dvdt_ofd, dudt_ogw, dvdt_ogw,                      &
         dudt_oss, dvdt_oss, tot_zmtb, tot_zlwb, tot_zogw,            &
         tot_tofd, tot_mtb, tot_ogw, tot_ngw,                         &
         du3dt_mtb,du3dt_ogw, du3dt_tms, du3dt_ngw, dv3dt_ngw,        &
         dudt_ngw, dvdt_ngw, dtdt_ngw,                                &
         ldu3dt_ngw, ldv3dt_ngw, ldt3dt_ngw,                          &
         dws3dt_ogw, dws3dt_obl, dws3dt_oss, dws3dt_ofd,              &
         ldu3dt_ogw, ldu3dt_obl, ldu3dt_oss, ldu3dt_ofd,              &
         du_ogwcol, dv_ogwcol, du_oblcol, dv_oblcol, du_osscol,       &
         dv_osscol, du_ofdcol, dv_ofdcol, du3_ogwcol, dv3_ogwcol,     &
         du3_oblcol, dv3_oblcol, du3_osscol, dv3_osscol, du3_ofdcol,  &
         dv3_ofdcol, errmsg, errflg)

        use machine,                only: kind_phys
        use mo_ugwpv1_kind,         only: wp

        implicit none

        ! Interface variables
        integer,              intent(in) :: im, levs
        real(kind=kind_phys), intent(in) :: dtf
        logical,              intent(in) :: ldiag_ugwp      !< flag for CIRES UGWP Diagnostics

        real(kind=kind_phys), intent(in),    dimension(:)   :: zobl, zlwb, zogw
        real(kind=kind_phys), intent(in),    dimension(:)   :: tau_ogw, tau_ngw
        real(kind=kind_phys), intent(in),    dimension(:),optional   :: du_ofdcol, du_oblcol
        real(kind=kind_phys), intent(inout), dimension(:)   :: tot_mtb, tot_ogw, tot_tofd, tot_ngw
        real(kind=kind_phys), intent(inout), dimension(:)   :: tot_zmtb, tot_zlwb, tot_zogw

        real(kind=kind_phys), intent(in),    dimension(:,:) :: dudt_gw, dvdt_gw
        real(kind=kind_phys), intent(in),    dimension(:,:), optional :: dudt_obl, dvdt_obl, dudt_ogw
        real(kind=kind_phys), intent(in),    dimension(:,:), optional :: dvdt_ogw, dudt_ofd, dvdt_ofd
        real(kind=kind_phys), intent(in),    dimension(:,:), optional :: dudt_oss, dvdt_oss
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: du3dt_mtb, du3dt_ogw, du3dt_tms
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: du3dt_ngw, dv3dt_ngw
        real(kind=kind_phys), intent(in),    dimension(:,:), optional :: dudt_ngw, dvdt_ngw, dtdt_ngw
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: ldu3dt_ngw, ldv3dt_ngw, ldt3dt_ngw
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: dws3dt_ogw, dws3dt_obl
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: dws3dt_oss, dws3dt_ofd
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: ldu3dt_ogw, ldu3dt_obl
        real(kind=kind_phys), intent(inout), dimension(:,:), optional :: ldu3dt_oss, ldu3dt_ofd
        real(kind=kind_phys), intent(in),    dimension(:), optional   :: du_ogwcol, dv_ogwcol
        real(kind=kind_phys), intent(in),    dimension(:), optional   :: dv_oblcol
        real(kind=kind_phys), intent(in),    dimension(:), optional   :: du_osscol, dv_osscol
        real(kind=kind_phys), intent(in),    dimension(:), optional   :: dv_ofdcol
        real(kind=kind_phys), intent(inout), dimension(:), optional   :: du3_ogwcol, dv3_ogwcol
        real(kind=kind_phys), intent(inout), dimension(:), optional   :: du3_oblcol, dv3_oblcol
        real(kind=kind_phys), intent(inout), dimension(:), optional   :: du3_osscol, dv3_osscol
        real(kind=kind_phys), intent(inout), dimension(:), optional   :: du3_ofdcol, dv3_ofdcol

        character(len=*),        intent(out) :: errmsg
        integer,                 intent(out) :: errflg

        ! Working-precision (_loc) buffers for CCPP host boundary isolation
        real(kind=wp) :: dtf_loc
        real(kind=wp), dimension(im)       :: zobl_loc, zlwb_loc, zogw_loc
        real(kind=wp), dimension(im)       :: tau_ogw_loc, tau_ngw_loc
        real(kind=wp), dimension(im)       :: du_ofdcol_loc, du_oblcol_loc, dv_oblcol_loc
        real(kind=wp), dimension(im)       :: du_ogwcol_loc, dv_ogwcol_loc
        real(kind=wp), dimension(im)       :: du_osscol_loc, dv_osscol_loc, dv_ofdcol_loc
        real(kind=wp), dimension(im)       :: tot_mtb_loc, tot_ogw_loc, tot_tofd_loc, tot_ngw_loc
        real(kind=wp), dimension(im)       :: tot_zmtb_loc, tot_zlwb_loc, tot_zogw_loc
        real(kind=wp), dimension(im)       :: du3_ogwcol_loc, dv3_ogwcol_loc
        real(kind=wp), dimension(im)       :: du3_oblcol_loc, dv3_oblcol_loc
        real(kind=wp), dimension(im)       :: du3_osscol_loc, dv3_osscol_loc
        real(kind=wp), dimension(im)       :: du3_ofdcol_loc, dv3_ofdcol_loc

        real(kind=wp), dimension(im, levs) :: dudt_gw_loc, dvdt_gw_loc
        real(kind=wp), dimension(im, levs) :: dudt_obl_loc, dvdt_obl_loc, dudt_ogw_loc
        real(kind=wp), dimension(im, levs) :: dvdt_ogw_loc, dudt_ofd_loc, dvdt_ofd_loc
        real(kind=wp), dimension(im, levs) :: dudt_oss_loc, dvdt_oss_loc
        real(kind=wp), dimension(im, levs) :: du3dt_mtb_loc, du3dt_ogw_loc, du3dt_tms_loc
        real(kind=wp), dimension(im, levs) :: du3dt_ngw_loc, dv3dt_ngw_loc
        real(kind=wp), dimension(im, levs) :: dudt_ngw_loc, dvdt_ngw_loc, dtdt_ngw_loc
        real(kind=wp), dimension(im, levs) :: ldu3dt_ngw_loc, ldv3dt_ngw_loc, ldt3dt_ngw_loc
        real(kind=wp), dimension(im, levs) :: dws3dt_ogw_loc, dws3dt_obl_loc
        real(kind=wp), dimension(im, levs) :: dws3dt_oss_loc, dws3dt_ofd_loc
        real(kind=wp), dimension(im, levs) :: ldu3dt_ogw_loc, ldu3dt_obl_loc
        real(kind=wp), dimension(im, levs) :: ldu3dt_oss_loc, ldu3dt_ofd_loc

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0
!
! post creates the "time-averaged" diagnostics"
!

        if (ldiag_ugwp) then
          ! Cast CCPP kind_phys inputs/inouts to working-precision (_loc) buffers
          dtf_loc            = real(dtf,               kind=wp)
          zobl_loc(1:im)     = real(zobl(1:im),        kind=wp)
          zlwb_loc(1:im)     = real(zlwb(1:im),        kind=wp)
          zogw_loc(1:im)     = real(zogw(1:im),        kind=wp)
          tau_ogw_loc(1:im)  = real(tau_ogw(1:im),     kind=wp)
          tau_ngw_loc(1:im)  = real(tau_ngw(1:im),     kind=wp)

          du_ofdcol_loc(1:im)= real(du_ofdcol(1:im),   kind=wp)
          du_oblcol_loc(1:im)= real(du_oblcol(1:im),   kind=wp)
          dv_oblcol_loc(1:im)= real(dv_oblcol(1:im),   kind=wp)
          du_ogwcol_loc(1:im)= real(du_ogwcol(1:im),   kind=wp)
          dv_ogwcol_loc(1:im)= real(dv_ogwcol(1:im),   kind=wp)
          du_osscol_loc(1:im)= real(du_osscol(1:im),   kind=wp)
          dv_osscol_loc(1:im)= real(dv_osscol(1:im),   kind=wp)
          dv_ofdcol_loc(1:im)= real(dv_ofdcol(1:im),   kind=wp)

          tot_mtb_loc(1:im)  = real(tot_mtb(1:im),     kind=wp)
          tot_ogw_loc(1:im)  = real(tot_ogw(1:im),     kind=wp)
          tot_tofd_loc(1:im) = real(tot_tofd(1:im),    kind=wp)
          tot_ngw_loc(1:im)  = real(tot_ngw(1:im),     kind=wp)
          tot_zmtb_loc(1:im) = real(tot_zmtb(1:im),    kind=wp)
          tot_zlwb_loc(1:im) = real(tot_zlwb(1:im),    kind=wp)
          tot_zogw_loc(1:im) = real(tot_zogw(1:im),    kind=wp)

          du3_ogwcol_loc(1:im)= real(du3_ogwcol(1:im), kind=wp)
          dv3_ogwcol_loc(1:im)= real(dv3_ogwcol(1:im), kind=wp)
          du3_oblcol_loc(1:im)= real(du3_oblcol(1:im), kind=wp)
          dv3_oblcol_loc(1:im)= real(dv3_oblcol(1:im), kind=wp)
          du3_osscol_loc(1:im)= real(du3_osscol(1:im), kind=wp)
          dv3_osscol_loc(1:im)= real(dv3_osscol(1:im), kind=wp)
          du3_ofdcol_loc(1:im)= real(du3_ofdcol(1:im), kind=wp)
          dv3_ofdcol_loc(1:im)= real(dv3_ofdcol(1:im), kind=wp)

          dudt_gw_loc(1:im,1:levs)  = real(dudt_gw(1:im,1:levs),  kind=wp)
          dvdt_gw_loc(1:im,1:levs)  = real(dvdt_gw(1:im,1:levs),  kind=wp)
          dudt_obl_loc(1:im,1:levs) = real(dudt_obl(1:im,1:levs), kind=wp)
          dvdt_obl_loc(1:im,1:levs) = real(dvdt_obl(1:im,1:levs), kind=wp)
          dudt_ogw_loc(1:im,1:levs) = real(dudt_ogw(1:im,1:levs), kind=wp)
          dvdt_ogw_loc(1:im,1:levs) = real(dvdt_ogw(1:im,1:levs), kind=wp)
          dudt_ofd_loc(1:im,1:levs) = real(dudt_ofd(1:im,1:levs), kind=wp)
          dvdt_ofd_loc(1:im,1:levs) = real(dvdt_ofd(1:im,1:levs), kind=wp)
          dudt_oss_loc(1:im,1:levs) = real(dudt_oss(1:im,1:levs), kind=wp)
          dvdt_oss_loc(1:im,1:levs) = real(dvdt_oss(1:im,1:levs), kind=wp)
          dudt_ngw_loc(1:im,1:levs) = real(dudt_ngw(1:im,1:levs), kind=wp)
          dvdt_ngw_loc(1:im,1:levs) = real(dvdt_ngw(1:im,1:levs), kind=wp)
          dtdt_ngw_loc(1:im,1:levs) = real(dtdt_ngw(1:im,1:levs), kind=wp)

          du3dt_mtb_loc(1:im,1:levs) = real(du3dt_mtb(1:im,1:levs), kind=wp)
          du3dt_ogw_loc(1:im,1:levs) = real(du3dt_ogw(1:im,1:levs), kind=wp)
          du3dt_tms_loc(1:im,1:levs) = real(du3dt_tms(1:im,1:levs), kind=wp)
          du3dt_ngw_loc(1:im,1:levs) = real(du3dt_ngw(1:im,1:levs), kind=wp)
          dv3dt_ngw_loc(1:im,1:levs) = real(dv3dt_ngw(1:im,1:levs), kind=wp)
          ldu3dt_ngw_loc(1:im,1:levs)= real(ldu3dt_ngw(1:im,1:levs),kind=wp)
          ldv3dt_ngw_loc(1:im,1:levs)= real(ldv3dt_ngw(1:im,1:levs),kind=wp)
          ldt3dt_ngw_loc(1:im,1:levs)= real(ldt3dt_ngw(1:im,1:levs),kind=wp)
          dws3dt_ogw_loc(1:im,1:levs)= real(dws3dt_ogw(1:im,1:levs),kind=wp)
          dws3dt_obl_loc(1:im,1:levs)= real(dws3dt_obl(1:im,1:levs),kind=wp)
          dws3dt_oss_loc(1:im,1:levs)= real(dws3dt_oss(1:im,1:levs),kind=wp)
          dws3dt_ofd_loc(1:im,1:levs)= real(dws3dt_ofd(1:im,1:levs),kind=wp)
          ldu3dt_ogw_loc(1:im,1:levs)= real(ldu3dt_ogw(1:im,1:levs),kind=wp)
          ldu3dt_obl_loc(1:im,1:levs)= real(ldu3dt_obl(1:im,1:levs),kind=wp)
          ldu3dt_oss_loc(1:im,1:levs)= real(ldu3dt_oss(1:im,1:levs),kind=wp)
          ldu3dt_ofd_loc(1:im,1:levs)= real(ldu3dt_ofd(1:im,1:levs),kind=wp)

          ! Perform diagnostic accumulation in working precision (wp)
          tot_zmtb_loc =  tot_zmtb_loc + dtf_loc *zobl_loc
          tot_zlwb_loc =  tot_zlwb_loc + dtf_loc *zlwb_loc
          tot_zogw_loc =  tot_zogw_loc + dtf_loc *zogw_loc
    
          tot_tofd_loc  = tot_tofd_loc + dtf_loc *du_ofdcol_loc
          tot_mtb_loc   = tot_mtb_loc +  dtf_loc *du_oblcol_loc
          tot_ogw_loc   = tot_ogw_loc +  dtf_loc *tau_ogw_loc
          tot_ngw_loc   = tot_ngw_loc +  dtf_loc *tau_ngw_loc
    
          du3dt_mtb_loc = du3dt_mtb_loc + dtf_loc *dudt_obl_loc
          du3dt_tms_loc = du3dt_tms_loc + dtf_loc *dudt_ofd_loc
          du3dt_ogw_loc = du3dt_ogw_loc + dtf_loc *dudt_ogw_loc
          du3dt_ngw_loc = du3dt_ngw_loc + dtf_loc *dudt_gw_loc
          dv3dt_ngw_loc = dv3dt_ngw_loc + dtf_loc *dvdt_gw_loc

          dws3dt_ogw_loc = dws3dt_ogw_loc + dtf_loc *sqrt(dudt_ogw_loc**2+dvdt_ogw_loc**2)
          dws3dt_obl_loc = dws3dt_obl_loc + dtf_loc *sqrt(dudt_obl_loc**2+dvdt_obl_loc**2)

          du3_ogwcol_loc = du3_ogwcol_loc + dtf_loc *du_ogwcol_loc
          dv3_ogwcol_loc = dv3_ogwcol_loc + dtf_loc *dv_ogwcol_loc
          du3_oblcol_loc = du3_oblcol_loc + dtf_loc *du_oblcol_loc
          dv3_oblcol_loc = dv3_oblcol_loc + dtf_loc *dv_oblcol_loc

          dws3dt_oss_loc = dws3dt_oss_loc + dtf_loc *sqrt(dudt_oss_loc**2+dvdt_oss_loc**2)
          dws3dt_ofd_loc = dws3dt_ofd_loc + dtf_loc *sqrt(dudt_ofd_loc**2+dvdt_ofd_loc**2)

          ldu3dt_ogw_loc = ldu3dt_ogw_loc + dtf_loc *dudt_ogw_loc
          ldu3dt_obl_loc = ldu3dt_obl_loc + dtf_loc *dudt_obl_loc
          ldu3dt_oss_loc = ldu3dt_oss_loc + dtf_loc *dudt_oss_loc
          ldu3dt_ofd_loc = ldu3dt_ofd_loc + dtf_loc *dudt_ofd_loc

          du3_osscol_loc = du3_osscol_loc + dtf_loc*du_osscol_loc
          dv3_osscol_loc = dv3_osscol_loc + dtf_loc*dv_osscol_loc
          du3_ofdcol_loc = du3_ofdcol_loc + dtf_loc*du_ofdcol_loc
          dv3_ofdcol_loc = dv3_ofdcol_loc + dtf_loc*dv_ofdcol_loc

          ! Special treatment for non-stationary GWD diagnostics
          ldu3dt_ngw_loc = ldu3dt_ngw_loc + dtf_loc *dudt_ngw_loc
          ldv3dt_ngw_loc = ldv3dt_ngw_loc + dtf_loc *dvdt_ngw_loc
          ldt3dt_ngw_loc = ldt3dt_ngw_loc + dtf_loc *dtdt_ngw_loc

          ! Restore working-precision (_loc) outputs back to CCPP host kind_phys
          tot_mtb(1:im)     = real(tot_mtb_loc(1:im),     kind=kind_phys)
          tot_ogw(1:im)     = real(tot_ogw_loc(1:im),     kind=kind_phys)
          tot_tofd(1:im)    = real(tot_tofd_loc(1:im),    kind=kind_phys)
          tot_ngw(1:im)     = real(tot_ngw_loc(1:im),     kind=kind_phys)
          tot_zmtb(1:im)    = real(tot_zmtb_loc(1:im),    kind=kind_phys)
          tot_zlwb(1:im)    = real(tot_zlwb_loc(1:im),    kind=kind_phys)
          tot_zogw(1:im)    = real(tot_zogw_loc(1:im),    kind=kind_phys)

          du3_ogwcol(1:im)  = real(du3_ogwcol_loc(1:im),  kind=kind_phys)
          dv3_ogwcol(1:im)  = real(dv3_ogwcol_loc(1:im),  kind=kind_phys)
          du3_oblcol(1:im)  = real(du3_oblcol_loc(1:im),  kind=kind_phys)
          dv3_oblcol(1:im)  = real(dv3_oblcol_loc(1:im),  kind=kind_phys)
          du3_osscol(1:im)  = real(du3_osscol_loc(1:im),  kind=kind_phys)
          dv3_osscol(1:im)  = real(dv3_osscol_loc(1:im),  kind=kind_phys)
          du3_ofdcol(1:im)  = real(du3_ofdcol_loc(1:im),  kind=kind_phys)
          dv3_ofdcol(1:im)  = real(dv3_ofdcol_loc(1:im),  kind=kind_phys)

          du3dt_mtb(1:im,1:levs)  = real(du3dt_mtb_loc(1:im,1:levs),  kind=kind_phys)
          du3dt_ogw(1:im,1:levs)  = real(du3dt_ogw_loc(1:im,1:levs),  kind=kind_phys)
          du3dt_tms(1:im,1:levs)  = real(du3dt_tms_loc(1:im,1:levs),  kind=kind_phys)
          du3dt_ngw(1:im,1:levs)  = real(du3dt_ngw_loc(1:im,1:levs),  kind=kind_phys)
          dv3dt_ngw(1:im,1:levs)  = real(dv3dt_ngw_loc(1:im,1:levs),  kind=kind_phys)
          ldu3dt_ngw(1:im,1:levs) = real(ldu3dt_ngw_loc(1:im,1:levs), kind=kind_phys)
          ldv3dt_ngw(1:im,1:levs) = real(ldv3dt_ngw_loc(1:im,1:levs), kind=kind_phys)
          ldt3dt_ngw(1:im,1:levs) = real(ldt3dt_ngw_loc(1:im,1:levs), kind=kind_phys)
          dws3dt_ogw(1:im,1:levs) = real(dws3dt_ogw_loc(1:im,1:levs), kind=kind_phys)
          dws3dt_obl(1:im,1:levs) = real(dws3dt_obl_loc(1:im,1:levs), kind=kind_phys)
          dws3dt_oss(1:im,1:levs) = real(dws3dt_oss_loc(1:im,1:levs), kind=kind_phys)
          dws3dt_ofd(1:im,1:levs) = real(dws3dt_ofd_loc(1:im,1:levs), kind=kind_phys)
          ldu3dt_ogw(1:im,1:levs) = real(ldu3dt_ogw_loc(1:im,1:levs), kind=kind_phys)
          ldu3dt_obl(1:im,1:levs) = real(ldu3dt_obl_loc(1:im,1:levs), kind=kind_phys)
          ldu3dt_oss(1:im,1:levs) = real(ldu3dt_oss_loc(1:im,1:levs), kind=kind_phys)
          ldu3dt_ofd(1:im,1:levs) = real(ldu3dt_ofd_loc(1:im,1:levs), kind=kind_phys)
        endif

!=====================================================================
! Updates inside the ugwpv1_gsldrag.F90
!
!
!       "post" may  also create the "time-averaged" diagnostics"
!            
!     if(ldiag3d .and. lssav .and. .not. flag_for_gwd_generic_tend) then
!        do k=1,levs
!          do i=1,im
!             ldu3dt_ngw(i,k) = ldu3dt_ngw(i,k) + dudt_ngw(i,k)*dtf
!             ldv3dt_ngw(i,k) = ldv3dt_ngw(i,k) + dvdt_ngw(i,k)*dtf
!             ldt3dt_ngw(i,k) = ldt3dt_ngw(i,k) + dtdt_ngw(i,k)*dtf
!         
!             ldu3dt_ogw(i,k) = ldu3dt_ogw(i,k) + dudt_ogw(i,k)*dtf
!             ldv3dt_ogw(i,k) = ldv3dt_ogw(i,k) + dvdt_ogw(i,k)*dtf
!             ldt3dt_ogw(i,k) = ldt3dt_ogw(i,k) + dtdt_ogw(i,k)*dtf
!          enddo
!        enddo
!      endif
! 
!=====================================================================
      end subroutine ugwpv1_gsldrag_post_run      

end module ugwpv1_gsldrag_post
