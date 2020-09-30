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

  ! Kinematics related variables 
  real(dp), public :: M
  real(dp), public :: roots, pt 

  ! Strong coupling alpha_s 
  real(dp), public :: alphas 

end module common_vars 



