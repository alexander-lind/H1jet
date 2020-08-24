!========================================================
!--------------------------------------------------------
! Module containing the wrapper function for the Hoppet
! integrator function 
!--------------------------------------------------------
!========================================================

module gauss_integrator 

  use hoppet_v1
  use integrator ! From Hoppet 

  implicit none 

  private 

  public :: gauss_integrate 

contains 

!======================================================================================= 
! Adaptive Gaussian quadrature integration from Hoppet 

  function gauss_integrate(f, a, b, eps) result(res) 

    real(dp), intent(in) :: a, b, eps 
    real(dp) :: res 

    interface
      function f(x)
        use hoppet_v1 
        real(dp) :: f
        real(dp), intent(in) :: x
      end function f
    end interface

    res = ig_LinWeight(f, a, b, one, one, eps) 

  end function gauss_integrate

!======================================================================================= 

end module gauss_integrator


