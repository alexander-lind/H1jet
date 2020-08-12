!========================================================
!--------------------------------------------------------
! Module for calculating the cross-section for vector
! boson production with a jet 
!--------------------------------------------------------
!========================================================

module vboson

  use types
  use consts_dp
  use ew_parameters
  use qcd

  implicit none

contains

!=======================================================================================
! Cross-section 

  function vboson_cross_section(lumiqq) result(res)

    real(dp), intent(in) :: lumiqq  
    real(dp) :: res

    res = pi / three * sqrt(two) * GF_GeVm2_in * lumiqq * invGev2_to_nb

  end function vboson_cross_section

!=======================================================================================
! Matrix element squared 

  subroutine vboson_Msquared(s, t, u, wtgg, wtqg, wtgq, wtqq)

    real(dp), intent(in) :: s, t, u    
    real(dp), intent(out):: wtqq, wtqg, wtgq, wtgg

    wtqq = cf * ((s + t)**2 + (s + u)**2) / (u * t)
    wtqg = tr * ((s + t)**2 + (u + t)**2) / (-s * u)
    wtgq = tr * ((s + u)**2 + (t + u)**2) / (-s * t)
    wtgg = zero

  end subroutine vboson_Msquared

!=======================================================================================

end module vboson


