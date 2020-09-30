!========================================================
!--------------------------------------------------------
! Module containing relevant SM and EW parameters 
! Values from PDG 
!--------------------------------------------------------
!========================================================

module ew_parameters

  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF

  implicit none
  
  ! Constant default parameters 
  real(dp), parameter, public :: mh       = 125.0_dp 
  real(dp), parameter, public :: mz       = 91.1876_dp 
  real(dp), parameter, public :: mw       = 80.385_dp 
  real(dp), parameter, public :: GF_GeVm2 = 0.116638e-04_dp 

  ! Unit conversion factors 
  real(dp), parameter, public :: invGev2_to_nb = 389379.323_dp
  real(dp), parameter, public :: invGev2_to_fb = 3.89379323e11_dp 

  ! Parameters which can be changed by user input 
  real(dp), public :: mz_in, mw_in, sinwsq_in, higgs_vev_in, GF_GeVm2_in

  real(dp), public :: gv2_ga2_u, gv2_ga2_d 
  real(dp), public :: gv2_ga2(6)

contains 

!=======================================================================================
! Function to set the vector axial-vector coupling ratios for quarks 

  subroutine set_gv2

    gv2_ga2_u = (half - four / three * sinwsq_in)**2 + one / four
    gv2_ga2_d = (-half + two / three * sinwsq_in)**2 + one / four
    gv2_ga2   = (/gv2_ga2_d, gv2_ga2_u, &
                & gv2_ga2_d, gv2_ga2_u, &
                & gv2_ga2_d, gv2_ga2_u /)

  end subroutine set_gv2

!=======================================================================================

end module ew_parameters


