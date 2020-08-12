!========================================================
!--------------------------------------------------------
! Module containing all public variables  
!--------------------------------------------------------
!========================================================

module common_vars 

  use hoppet_v1

  implicit none

  ! Output level 
  integer, public, parameter :: stdout = 6
  integer, public, parameter :: stderr = 0
  integer, public :: idev

  ! Warning parameter for wae_warn from HOPPET 
  integer, public, save :: warn_param = 3

  ! Process id 
  integer, public, parameter :: id_H = 1
  integer, public, parameter :: id_Z = 2
  integer, public, parameter :: id_bbH = 3
  integer, public, parameter :: id_user = 4
  integer, public :: iproc

  ! Array specifying the particles in the loops 
  integer, pointer, public :: iloop_array(:) 

  ! Collider type (pp or ppbar) 
  character(len=5), public :: collider 

  ! Loop approximation 
  character(len=8), public :: approx 

  ! Monte Carlo integration accuracy 
  real(dp), public :: accuracy

  ! Kinematics related variables 
  real(dp), public :: M
  real(dp), public :: roots, rootshat, pt 
  real(dp), public :: p(4,4), lnpt, y, yfrac
  real(dp), public :: xmom, ymin, ymax, ptmin, ptmax

  ! Model parameters 
  real(dp), public :: yt, yb, ytp
  real(dp), public :: mtp, sth2, cth2, s2th2
  real(dp), public :: yst1, yst2
  real(dp), public :: delta, alpha1, alpha2, tbeta, c2beta
  real(dp), public :: mst1, mst2 
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


