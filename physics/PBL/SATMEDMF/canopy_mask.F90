   module canopy_mask_mod

   use machine , only : kind_phys
   use mo_pbl_kind, only : pbl_wp

   implicit none

   public :: canopy_mask_init, canopy_mask_run

   contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   subroutine canopy_mask_init(im, km, nkc, nkt, &
           claie, cfch, cfrt, cclu, cpopu, &  !in:
           FRT_mask,                       &  ! out
           errmsg,errflg)

   implicit none

! Horizontal arrays
   integer, intent(in)  :: im, km  ! horizontal & vertical domain specifications
   integer, intent(in)  :: nkc, nkt

   real(kind=pbl_wp), intent(in) :: claie(im), cfch(im), cfrt(im), &
                                    cclu(im), cpopu(im)
   real(kind=pbl_wp), intent(out) :: FRT_mask(im)

   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

!...local variables

! Initialize CCPP error handling variables
   errmsg = ''
   errflg = 0

!...Allocate and initialize new canopy arrays

! Initializations

   FRT_mask(:)=0.0_pbl_wp

   return
   end subroutine canopy_mask_init

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   subroutine canopy_mask_run (im, km, nkc, nkt, &  !in:
           claie, cfch, cfrt, cclu, cpopu, &  !in:
           FRT_mask,                       &  !out:
           errmsg,errflg)

   implicit none

!...Arguments:

! Horizontal arrays
   integer, intent(in)  :: im, km  ! horizontal & vertical domain specifications
   integer, intent(in) :: nkc, nkt

   real(kind=pbl_wp), intent(in) :: claie(im), cfch(im),  cfrt(im), &
                                    cclu(im), cpopu(im)
   real(kind=pbl_wp), intent(out) :: FRT_mask(im)

   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

!...local variables

   integer i

! Initialize CCPP error handling variables
   errmsg = ''
   errflg = 0

   do i=1,im

      !NOT a Continuous forest canopy
      if (   claie(i) < 0.1_pbl_wp                                       &
        .OR. cfch(i)  < 0.5_pbl_wp                                       &
!IVAI: modified contiguous canopy condition
!        .OR. MAX(0.0, 1.0 - cfrt(i)) .GT. 0.5 
        .OR. MAX(0.0_pbl_wp, 1.0_pbl_wp - cfrt(i)) > 0.75_pbl_wp         &
        .OR. cpopu(i) > 10000.0_pbl_wp                                   &
        .OR. (EXP(-0.5_pbl_wp * claie(i) * cclu(i)) > 0.45_pbl_wp        &
        .AND. cfch(i) < 18.0_pbl_wp) ) THEN

         FRT_mask(i) = -1.0_pbl_wp

      ! Continuous forest canopy
      ELSE

         FRT_mask(i) = 1.0_pbl_wp

      END IF ! Forest Canopy Mask

   end do

   return
   end subroutine canopy_mask_run

   end module canopy_mask_mod
