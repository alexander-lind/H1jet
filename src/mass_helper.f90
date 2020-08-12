!========================================================
!--------------------------------------------------------
! Module containing various definitions of quark masses, 
! and relevant functions to handle them 
!--------------------------------------------------------
!========================================================

module mass_helper

  use hoppet_v1 

  implicit none

  private

  ! Constant default heavy quark masses 
  real(dp), parameter, public :: mb = 4.65_dp  ! from PDG 2012
  real(dp), parameter, public :: mt = 173.5_dp ! from PDG 2012

  ! Heavy quark masses which can be changed by user input 
  real(dp), public :: mb_in, mt_in 

  ! Anomalous Yukawa couplings
  real(dp), pointer, public :: yukawa(:)
  real(dp), pointer, public :: mass_array(:)
  logical, public :: cpodd = .false.
  
  real(dp), parameter :: gamma0 = one
  real(dp), parameter :: gamma1 = one / 16._dp * (202._dp / three - 20._dp / 9._dp * nf_def)

  public :: RunningMass

contains
  
!=======================================================================================
! Calculate the running quark mass 

  function RunningMass(muR, m0, as, as0) result(res)
  
    real(dp), intent(in) :: muR, m0, as, as0
    real(dp) :: A1, res

    A1 = -beta1 * gamma0 / beta0 + gamma1 / (pi * beta0) 

    res = m0 * (as / as0)**(gamma0 / (pi * beta0)) * (one + A1 * as / pi) / (one + A1 * as0 / pi)

    ! Hack to check with SusHi
    !res = m0 * ((one - two * 0.120_dp * beta0 * log(91.2_dp / m0)) / &
    !      & (one - two * 0.120_dp * beta0 * log(91.2_dp / muR)))**(gamma0 / (pi * beta0))

    ! New hack to check with SusHi 
    !res = m0 * (as / as0)**(gamma0 / (pi * beta0)) * ((one + (gamma1 / (pi * beta0) &
    !& - beta1 * gamma0 / (pi**2 * beta0**2)) * as / pi) / (one + (gamma1 / (pi * beta0) &
    !& - beta1 * gamma0 / (pi**2 * beta0**2)) * as0 / pi)) 

  end function RunningMass

!=======================================================================================

end module mass_helper


