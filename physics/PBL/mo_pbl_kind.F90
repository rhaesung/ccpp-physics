!> This module provides the Fortran KIND parameters for REAL variables 
!! in the PBL schemes, synchronized with the model's machine definitions.
module mo_pbl_kind
  use machine, only : kind_phys, kind_sngl_prec, kind_dbl_prec
  implicit none
  public

  ! Define standard single and double precision kinds from machine module
  integer, parameter :: sp = kind_sngl_prec  ! 4
  integer, parameter :: dp = kind_dbl_prec   ! 8

  ! Floating point working precision (pbl_wp)
  ! 1. If PBL_USE_SP is defined, force 4-byte precision for PBL.
  ! 2. Otherwise, match kind_phys to ensure B4B with the legacy model.
#ifdef PBL_USE_SP
  integer, parameter :: pbl_wp = sp
#else
  integer, parameter :: pbl_wp = kind_phys
#endif

end module mo_pbl_kind
