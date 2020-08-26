!========================================================
!--------------------------------------------------------
! Default module for the squared matrix elements from
! included user code for a custom process 
!--------------------------------------------------------
!========================================================

module user_interface 

  use ew_parameters
  use hoppet_v1 
  use common_vars 
  
  implicit none
  
  private 

  integer, parameter :: idef = selected_int_kind(8)
  
  public :: user_included, user_cross_section, user_Msquared, user_luminosities, &
            & user_help_message 

contains
  
!=======================================================================================
! Test function to see if the user interface is set up and included 
  
  integer function user_included() 

    ! This file is the default file where the user interface has not been set up 
    ! Please see the manual for the implementation 

    ! id_noImplementation   = user interface has not been set up 
    ! id_missingTotXsec     = total born-level cross-section not included 
    ! id_fullImplementation = full implementation 
    user_included = id_noImplementation 

  end function user_included 

!=======================================================================================
! Born-level cross-section from the generated user code 

  function user_cross_section(lumigg, lumiqg, lumigq, lumiqqbar) result(res) 

    real(dp) :: res 
    real(dp), intent(in) :: lumigg, lumiqg, lumigq, lumiqqbar 

    res = 0.0_dp 

  end function user_cross_section 

!=======================================================================================
! Squared matrix element from the generated user code
  
  subroutine user_Msquared(s, t, u, wtqq, wtqg, wtgq, wtgg)
  
    real(dp), intent(in) :: s, t, u ! Input
    real(dp), intent(out) :: wtqq, wtqg, wtgq, wtgg ! Output
  
    ! Default values
    wtqq = 0.0_dp
    wtqg = 0.0_dp
    wtgq = 0.0_dp
    wtgg = 0.0_dp
  
  end subroutine user_Msquared
  
!=======================================================================================
! Luminosities for the process from the generated user code
  
  subroutine user_luminosities(grid, pdf1, pdf2, lumi_gg, lumi_qg, lumi_gq, lumi_qqbar)
  
    type(grid_def), intent(in) :: grid
    real(dp), intent(in) :: pdf1(0:grid%ny,-6:7), pdf2(0:grid%ny,-6:7)
    real(dp), intent(inout) :: lumi_gg(0:), lumi_qg(0:), lumi_gq(0:), lumi_qqbar(0:)
  
    ! Default values
    lumi_gg = 0.0_dp 
    lumi_gq = 0.0_dp 
    lumi_qg = 0.0_dp 
    lumi_qqbar = 0.0_dp 
  
  end subroutine user_luminosities
  
!=======================================================================================
! Print any parameter additions to the help message
  
  subroutine user_help_message 

    ! Deliberately empty 
  
  end subroutine user_help_message
  
!=======================================================================================
  
end module user_interface 
  
  
