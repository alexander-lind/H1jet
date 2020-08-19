!========================================================
!--------------------------------------------------------
! Module containing commonly used public variables  
!--------------------------------------------------------
!========================================================

module common_vars 

  use hoppet_v1

  implicit none

  ! Output level 
  integer, public, parameter :: stdout = 6
  integer, public, parameter :: stderr = 0

  ! Process id 
  integer, public, parameter :: id_H = 1
  integer, public, parameter :: id_Z = 2
  integer, public, parameter :: id_bbH = 3
  integer, public, parameter :: id_user = 4
  integer, public :: iproc

  ! User interface related id's 
  integer, public, parameter :: id_noImplementation = 0
  integer, public, parameter :: id_missingTotXsec = 1 
  integer, public, parameter :: id_fullImplementation = 2

  ! Collider type (pp or ppbar)
  ! This is an input of get_pdfs, so maybe avoid making it public
  character(len=5), public :: collider 

  ! Loop approximation
  ! AB we need to decide how to treat quark loops
!  character(len=8), public :: approx 

  ! Monte Carlo integration accuracy
  ! This is an input to other subroutines, so it should  not be made public
  real(dp), public :: accuracy

  ! Kinematics related variables 
  real(dp), public :: M
  real(dp), public :: roots, pt 

!  ! Model parameters 
  real(dp), pointer, public :: mb0
  real(dp), public :: mass

  ! Histogram related variables 
  integer, public :: nbins 
  real(dp), allocatable, public :: binmin(:), binmed(:), binmax(:)
  real(dp), public :: bin_width

  ! alpha_ew and alpha_s couplings 
  real(dp), public :: alphaw 
  real(dp), public :: alphas 
  integer, public :: as_pow = 0
  real(dp), public :: as0

  ! Scales related variables 
  character(len=2), public :: scale_strategy
  real(dp), public :: muF, muR
  real(dp), public :: xmur, xmuf

  ! tau = M^2 / s 
  type(gdval), public :: tau

  ! PDF related variables 
  character(len=30), public :: pdf_name
  integer, public :: pdf_mem

  ! Cross-sections 
  real(dp), allocatable, public :: dsigma_dpt(:), sigma(:)
  real(dp), public :: sigma0, sigma_bbH
  real(dp), public :: jakob, dsigmadpt, ew_prefactor

end module common_vars 


