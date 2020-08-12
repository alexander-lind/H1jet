!========================================================
!--------------------------------------------------------
! Module for calculating the cross-section for Higgs 
! production at large transverse momentum 
!--------------------------------------------------------
!========================================================

module hboson

  use hoppet_v1
  use ew_parameters
  use scalar_integrals
  use mass_helper
  
  implicit none

  private

  ! Some useful identifiers
  integer, public, parameter :: iloop_sm_fermion=0
  integer, public, parameter :: iloop_fm_fermion=1
  integer, public, parameter :: iloop_lm_fermion=2
  integer, public, parameter :: iloop_fm_scalar =3
  integer, public, parameter :: iloop_lm_scalar =4
  
  real(dp), public :: mh2 ! Higgs mass squared
  public :: hboson_cross_section, bbH_cross_section
  public :: hboson_msquared, hboson_cpodd_Msquared, bbH_msquared
  
contains

!=======================================================================================
! The exact +++ helicity amplitude for gg -> gH 
! See eq. (A.15) in Baur and Glover, Nucl.Phys. B339 (1990) 38-66 
  function ggxgH_pppamp_fermions(s, t, u, m2, mq2, i, j, k, i1, j1, k1, BI, CI, DI) result(res)
   
    real(dp) :: s, t, u, m2, mq2
    integer :: i, j, k, i1, j1, k1
    complex(dp) :: res

    real(dp) :: s1, t1, u1
    complex(dp) :: BI(4), CI(7), DI(3)

    s1 = s - m2
    t1 = t - m2
    u1 = u - m2
   
    res = mq2 * four * sqrt(two * s * t * u) * &
          & (-four * (one / (u * t) + one / (u * u1) + one / (t * t1)) &
          & - four * ((two * s + t) * BI(k) / u1**2 + (two * s + u) * BI(j) / t1**2) / s &
          & - (s - four * mq2) * (s1 * CI(i1) + (u - s) * CI(j1) + (t - s) * CI(k1)) / (s * t * u) &
          & - 8._dp * mq2 * (CI(j1) / (t * t1) + CI(k1) / (u * u1))&
          & + half * (s - four * mq2) * (s * t * DI(k) + u * s * DI(j) - u * t * DI(i)) / (s * t * u) &
          & + four * mq2 * DI(i) / s &
          & - two * (u * CI(k) + t * CI(j) + u1 * CI(k1) + t1 * CI(j1) - u * t * DI(i)) / s**2)

  end function ggxgH_pppamp_fermions

!=======================================================================================
! The exact ++- helicity amplitude for gg -> gH 
! See eq. (A.16) in Baur and Glover, Nucl.Phys. B339 (1990) 38-66 
  function ggxgH_ppmamp_fermions(s, t, u, m2, mq2, i, j, k, i1, j1, k1, BI, CI, DI) result(res)

    real(dp) :: s, t, u, m2, mq2
    integer :: i, j, k, i1, j1, k1
    complex(dp) :: res

    real(dp) :: s1, t1, u1
    complex(dp) :: BI(4), CI(7), DI(3)

    s1 = s - m2
    t1 = t - m2
    u1 = u - m2

    res = mq2 * four * sqrt(two * s * t * u) * (four * m2 &
          & + (m2 - four * mq2) * (s1 * CI(4) + t1 * CI(5) + u1 * CI(6)) &
          & - half * (m2 - four * mq2) * (s * t * DI(3) + u * s * DI(2) + u * t * DI(1))) / (s * t * u)

  end function ggxgH_ppmamp_fermions

!=======================================================================================
! The exact amplitude for qqbar -> gH 
! See eq. (A.19) in Baur and Glover, Nucl.Phys. B339 (1990) 38-66 
  function qqbarxgH_amp_fermions(s, t, u, m2, mq2, i, j, k, i1, j1, k1, BI, CI, DI) result(res)

    real(dp) :: s, t, u, m2, mq2
    integer :: i, j, k, i1, j1, k1
    complex(dp) :: res

    complex(dp) :: BI(4), CI(7), DI(3)

    res = mq2 * (dcmplx(two) + dcmplx(two * s / (s - m2)) * BI(i) &
          &     + dcmplx(four * mq2 - u - t) * CI(k))

  end function qqbarxgH_amp_fermions

!=======================================================================================
! The +++ helicity amplitude for gg -> gH in the small mass limit for the loops 
! See eq. (A.15) in Baur and Glover, Nucl.Phys. B339 (1990) 38-66 
  function ggxgH_pppamp_fermions_sml(s, t, u, m2, mq2, i, j, k, i1, j1, k1, BI, CI, DI) result(res)

    real(dp) :: s, t, u, m2, mq2
    integer :: i, j, k, i1, j1, k1
    complex(dp) :: res

    real(dp) :: s1, t1, u1
    complex(dp) :: BI(4), CI(7), DI(3)

    s1 = s - m2
    t1 = t - m2
    u1 = u - m2

    res = mq2 * four * sqrt(two * s * t * u) * &
          & (- four * (one / (u * t) + one / (u * u1) + one / (t * t1)) &
          & - four * ((two * s + t) * BI(k) / u1**2 + (two * s + u) * BI(j) / t1**2) / s &
          & - (s - four * mq2) * (s1 * CI(i1) + (u - s) * CI(j1) + (t - s) * CI(k1)) / (s * t * u) &
          & - 8._dp * mq2 * (CI(j1) / (t * t1) + CI(k1) / (u * u1)) &
          & + half * (s - four * mq2) * (s * t * DI(k) + u * s * DI(j) - u * t * DI(i)) / (s * t * u) &
          & + four * mq2 * DI(i) / s &
          & - two * (u * CI(k) + t * CI(j) + u1 * CI(k1) + t1 * CI(j1) - u * t * DI(i)) / s**2)

  end function ggxgH_pppamp_fermions_sml

!=======================================================================================
! The ++- helicity amplitude for gg -> gH in the small mass limit for the loops 
! See eq. (A.16) in Baur and Glover, Nucl.Phys. B339 (1990) 38-66 
  function ggxgH_ppmamp_fermions_sml(s, t, u, m2, mq2, i, j, k, i1, j1, k1, BI, CI, DI) result(res)

    real(dp) :: s, t, u, m2, mq2
    integer :: i, j, k, i1, j1, k1
    complex(dp) :: res

    real(dp) :: s1, t1, u1
    complex(dp) :: BI(4), CI(7), DI(3)

    s1 = s - m2
    t1 = t - m2
    u1 = u - m2

    res = mq2 * four * sqrt(two * s * t * u) * (four * m2 &
          & + (m2 - four * mq2) * (s1 * CI(4) + t1 * CI(5) + u1 * CI(6)) &
          & - half * (m2 - four * mq2) &
          & * (s * t * DI(3) + u * s * DI(2) + u * t * DI(1))) / (s * t * u)

  end function ggxgH_ppmamp_fermions_sml

!=======================================================================================
! The amplitude for qqbar -> gH in the small mass limit for the loops 
! See eq. (A.19) in Baur and Glover, Nucl.Phys. B339 (1990) 38-66 
  function qqbarxgH_amp_fermions_sml(s, t, u, m2, mq2, i, j, k, i1, j1, k1, BI, CI, DI) result(res)

    real(dp) :: s, t, u, m2, mq2
    integer :: i, j, k, i1, j1, k1
    complex(dp) :: res

    complex(dp) :: BI(4), CI(7), DI(3)

    res = dcmplx(two) + dcmplx(two * s / (s - m2)) * BI(i) &
          & + dcmplx(four * mq2 - u - t) * CI(k) 

  end function qqbarxgH_amp_fermions_sml

!=======================================================================================
! The exact +++ helicity amplitude for gg -> gH 
! with scalars instead of fermions in the loop 
  function ggxgH_pppamp_scalars(s, t, u, m2, mq2, i, j, k, i1, j1, k1, BI, CI, DI) result(res)

    real(dp) :: s, t, u, m2, mq2
    integer :: i, j, k, i1, j1, k1
    complex(dp) :: res

    real(dp) :: s1, t1, u1
    complex(dp) :: BI(4), CI(7), DI(3)

    s1 = s - m2
    t1 = t - m2
    u1 = u - m2

    res = mq2 * sqrt(two * s * t * u) * &
          & (four * (one / (u * t) + one / (u * u1) + one / (t * t1)) &
          & + four * ((two * s + t) * BI(k) / u1**2 + (two * s + u) * BI(j) / t1**2) / s &
          & - four * mq2 * (s1 * CI(i1) + (u - s) * CI(j1) + (t - s) * CI(k1)) / (s * t * u) &
          & + 8._dp * mq2 * (CI(j1) / (t * t1) + CI(k1) / (u * u1)) &
          & + two * (u * CI(k) + t * CI(j) + u1 * CI(k1) + t1 * CI(j1) - u * t * DI(i)) / s**2 & 
          & - four * mq2 * DI(i) / s &
          & + two * mq2 * (s * t * DI(k) + u * s * DI(j) - u * t * DI(i)) / (s * t * u) )

  end function ggxgH_pppamp_scalars

!=======================================================================================
! The exact ++- helicity amplitude for gg -> gH 
! with scalars instead of fermions in the loop 

  function ggxgH_ppmamp_scalars(s, t, u, m2, mq2, i, j, k, i1, j1, k1, BI, CI, DI) result(res)

    real(dp):: s, t, u, m2, mq2
    integer :: i, j, k, i1, j1, k1
    complex(dp) :: res

    real(dp) :: s1, t1, u1
    complex(dp) :: BI(4), CI(7), DI(3)

    s1 = s - m2
    t1 = t - m2
    u1 = u - m2

    res = mq2 * sqrt(two * s * t * u) * (-four * m2 & 
          & + four * mq2 * (s1 * CI(4) + t1 * CI(5) + u1 * CI(6)) &
          & - two * mq2 * (s * t * DI(3) + u * s * DI(2) + u * t * DI(1))) / (s * t * u)

  end function ggxgH_ppmamp_scalars

!=======================================================================================
! The exact amplitude for qqbar -> gH 
! with scalars instead of fermions in the loop 

  function qqbarxgH_amp_scalars(s, t, u, m2, mq2, i, j, k, i1, j1, k1, BI, CI, DI) result(res)

    real(dp):: s, t, u, m2, mq2
    integer :: i, j, k, i1, j1, k1
    complex(dp) :: res

    real(dp) :: s1, t1, u1
    complex(dp) :: BI(4), CI(7), DI(3)
    
    res = -mq2 * half * (dcmplx(one) + dcmplx(s / (s - m2)) * BI(i) &
          &             + dcmplx(two * mq2) * CI(k))

  end function qqbarxgH_amp_scalars

!====================================================================================================
! The exact +++ CP-odd helicity amplitude for gg -> gH 
  function ggxgH_cpodd_pppamp_fermions(s, t, u, m2, mq2, i, j, k, i1, j1, k1, BI, CI, DI) result(res)
   
    real(dp) :: s, t, u, m2, mq2, Gst, Gus, Gut
    integer :: i, j, k, i1, j1, k1
    complex(dp) :: res

    complex(dp) :: BI(4), CI(7), DI(3)

    Gst = s*t*DI(k)+two*s*CI(j1)+two*t*CI(i1)
    Gus = u*s*DI(j)+two*u*CI(i1)+two*s*CI(k1)
    Gut = u*t*DI(i)+two*u*CI(j1)+two*t*CI(k1)
   
    res = mq2 * sqrt(s * t * u)/(s*u)*(Gst-Gus+Gut)
    
  end function ggxgH_cpodd_pppamp_fermions


  !====================================================================================================
! The exact ++- CP-odd helicity amplitude for gg -> gH 
  function ggxgH_cpodd_ppmamp_fermions(s, t, u, m2, mq2, i, j, k, i1, j1, k1, BI, CI, DI) result(res)
   
    real(dp) :: s, t, u, m2, mq2, Gst, Gus, Gut
    integer :: i, j, k, i1, j1, k1
    complex(dp) :: res

    complex(dp) :: BI(4), CI(7), DI(3)

    Gst = s*t*DI(k)+two*s*CI(j1)+two*t*CI(i1)
    Gus = u*s*DI(j)+two*u*CI(i1)+two*s*CI(k1)
    Gut = u*t*DI(i)+two*u*CI(j1)+two*t*CI(k1)

    res = -mq2**2 * sqrt(s * t * u)/(s*t*u)*(Gst+Gus+Gut)
    
  end function ggxgH_cpodd_ppmamp_fermions

    !====================================================================================================
! The exact -++ CP-odd helicity amplitude for gg -> gH 
  function ggxgH_cpodd_mppamp_fermions(s, t, u, m2, mq2, i, j, k, i1, j1, k1, BI, CI, DI) result(res)
   
    real(dp) :: s, t, u, m2, mq2, Gst, Gus, Gut
    integer :: i, j, k, i1, j1, k1
    complex(dp) :: res

    complex(dp) :: BI(4), CI(7), DI(3)

    Gst = s*t*DI(k)+two*s*CI(j1)+two*t*CI(i1)
    Gus = u*s*DI(j)+two*u*CI(i1)+two*s*CI(k1)
    Gut = u*t*DI(i)+two*u*CI(j1)+two*t*CI(k1)

    res = mq2* sqrt(s * t * u)/(t*u)*(Gst+Gus-Gut)
    
  end function ggxgH_cpodd_mppamp_fermions

  !=======================================================================================
! The exact amplitude for CP-odd qqbar -> gH 
  function qqbarxgH_cpodd_amp_fermions(s, t, u, m2, mq2, i, j, k, i1, j1, k1, BI, CI, DI) result(res)

    real(dp) :: s, t, u, m2, mq2
    integer :: i, j, k, i1, j1, k1
    complex(dp) :: res

    complex(dp) :: BI(4), CI(7), DI(3)

    res = mq2 * CI(k)

  end function qqbarxgH_cpodd_amp_fermions
  
!=======================================================================================
! Cross-section for gg -> H + jet 

  function hboson_cross_section(lumigg, iloop_array, mass_array, yukawa) result(res)

    real(dp) :: nc = three  
    real(dp) :: res, lumigg
    integer, intent(in) :: iloop_array(:)
    real(dp), intent(in) :: mass_array(:), yukawa(:)

    res = one / pi / (nc**2 - 1) / 72._dp / higgs_vev_in**2
    res = res * born_ampsquare(iloop_array, mass_array, yukawa)
    res = res * lumigg
    res = res * invGev2_to_nb

  end function hboson_cross_section

!=======================================================================================
! Matrix element squared for pp -> H + jet 
! Subprocesses: 
!    gg -> gH 
!    qqbar -> gH 
!    qg -> qH (obtained by permutation of the qqbar amplitude) 
!    gqbar -> qbar H (obtained by permutation of the qqbar amplitude) 
! Modified code from Herwig 6.521 

  subroutine hboson_Msquared(s, t, u, iloop_in, mass_in, yukawa_in,  wtqq, wtqg, wtgq, wtgg)
    real(dp), intent(in) :: s, t, u
    integer, intent(in) :: iloop_in(:)
    real(dp), intent(in) :: mass_in(:), yukawa_in(:)
    real(dp), intent(out) :: wtqq, wtqg, wtgq, wtgg

    complex(dp) :: BI(4), CI(7), DI(3), amp(7) 

    real(dp) :: mw2, rngluon, rnquark, fluxgg, fluxgq, fluxqq, mq2, yukawa 
    real(dp) :: amp_imag(7), amp_real(7)
    integer :: i, imass, nc, iloop
    
    mw2 = mw_in**2  

    ! Spin and colour flux factors plus enhancement factor
    nc = three ! Number of colors 
    rngluon = one / float(nc**2 - 1) 
    rnquark = one / float(nc) 
    fluxgg = 0.25_dp * rngluon**2
    fluxgq = 0.25_dp * rngluon * rnquark
    fluxqq = 0.25_dp * rnquark**2
      
    amp_imag = zero
    amp_real = zero
    
    do imass = 1, size(iloop_in)
    
      mq2 = mass_in(imass)**2
      if (mq2 == zero) cycle

      yukawa = yukawa_in(imass)
      
      ! Selection for loop quarks 
      iloop = iloop_in(imass)

      ! Evaluate scalar loop integrals 
         
      ! Two-point function 
      BI(1) = B0(s, mq2) 
      BI(2) = B0(t, mq2) 
      BI(3) = B0(u, mq2) 
      BI(4) = B0(mh2, mq2) 
      BI(1) = BI(1) - BI(4) ! B1(s) 
      BI(2) = BI(2) - BI(4) ! B1(t) 
      BI(3) = BI(3) - BI(4) ! B1(u) 

      ! Three-point function 
      CI(1) = C0(s, mq2) 
      CI(2) = C0(t, mq2) 
      CI(3) = C0(u, mq2) 
      CI(7) = C0(mh2, mq2) 
      CI(4) = (s * CI(1) - mh2 * CI(7)) / (s - mh2) ! C1(s) 
      CI(5) = (t * CI(2) - mh2 * CI(7)) / (t - mh2) ! C1(t) 
      CI(6) = (u * CI(3) - mh2 * CI(7)) / (u - mh2) ! C1(u) 
       
      ! Four-point function 
      DI(1) = D0(u, t, mh2, mq2) 
      DI(2) = D0(s, u, mh2, mq2) 
      DI(3) = D0(s, t, mh2, mq2) 

      ! Compute the complex amplitudes
      ! PFM hack to obtain the right small mass limit

      select case(iloop)

      case(iloop_sm_fermion)

         ! Small mass limit in loops 
        amp(1) = ggxgH_pppamp_fermions_sml(s, t, u, mh2, mq2, 1, 2, 3, 4, 5, 6, BI, CI, DI)
        amp(2) = ggxgH_ppmamp_fermions_sml(s, t, u, mh2, mq2, 1, 2, 3, 0, 0, 0, BI, CI, DI)
        amp(3) = ggxgH_pppamp_fermions_sml(t, s, u, mh2, mq2, 2, 1, 3, 5, 4, 6, BI, CI, DI)
        amp(4) = ggxgH_pppamp_fermions_sml(u, t, s, mh2, mq2, 3, 2, 1, 6, 5, 4, BI, CI, DI)
        amp(5) = qqbarxgH_amp_fermions_sml(s, t, u, mh2, mq2, 1, 0, 4, 0, 0, 0, BI, CI, DI)
        amp(6) = qqbarxgH_amp_fermions_sml(t, s, u, mh2, mq2, 2, 0, 5, 0, 0, 0, BI, CI, DI)
        amp(7) = qqbarxgH_amp_fermions_sml(u, t, s, mh2, mq2, 3, 0, 6, 0, 0, 0, BI, CI, DI)

     case(iloop_fm_fermion)

        ! Exact result for loops 
        amp(1) = ggxgH_pppamp_fermions(s, t, u, mh2, mq2, 1, 2, 3, 4, 5, 6, BI, CI, DI)
        amp(2) = ggxgH_ppmamp_fermions(s, t, u, mh2, mq2, 1, 2, 3, 0, 0, 0, BI, CI, DI)
        amp(3) = ggxgH_pppamp_fermions(t, s, u, mh2, mq2, 2, 1, 3, 5, 4, 6, BI, CI, DI)
        amp(4) = ggxgH_pppamp_fermions(u, t, s, mh2, mq2, 3, 2, 1, 6, 5, 4, BI, CI, DI)
        amp(5) = qqbarxgH_amp_fermions(s, t, u, mh2, mq2, 1, 0, 4, 0, 0, 0, BI, CI, DI)
        amp(6) = qqbarxgH_amp_fermions(t, s, u, mh2, mq2, 2, 0, 5, 0, 0, 0, BI, CI, DI)
        amp(7) = qqbarxgH_amp_fermions(u, t, s, mh2, mq2, 3, 0, 6, 0, 0, 0, BI, CI, DI)

     case(iloop_lm_fermion)

        ! Infinite mass limit in loops 
        amp(1) = -two * four / three * sqrt(two * s * t * u) * s**2 / (s * t * u)
        amp(2) = two * four / three * sqrt(two * s * t * u) * mh2**2 / (s * t * u)
        amp(3) = -two * four / three * sqrt(two * s * t * u) * t**2 / (s * t * u)
        amp(4) = -two * four / three * sqrt(two * s * t * u) * u**2 / (s * t * u)
        amp(5) = -(s - mh2) / three
        amp(6) = -(t - mh2) / three
        amp(7) = -(u - mh2) / three

        case(iloop_fm_scalar)

        ! Exact result for scalar particles in loops 
        amp(1) = ggxgH_pppamp_scalars(s, t, u, mh2, mq2, 1, 2, 3, 4, 5, 6, BI, CI, DI)
        amp(2) = ggxgH_ppmamp_scalars(s, t, u, mh2, mq2, 1, 2, 3, 0, 0, 0, BI, CI, DI)
        amp(3) = ggxgH_pppamp_scalars(t, s, u, mh2, mq2, 2, 1, 3, 5, 4, 6, BI, CI, DI)
        amp(4) = ggxgH_pppamp_scalars(u, t, s, mh2, mq2, 3, 2, 1, 6, 5, 4, BI, CI, DI)
        amp(5) = qqbarxgH_amp_scalars(s, t, u, mh2, mq2, 1, 0, 4, 0, 0, 0, BI, CI, DI)
        amp(6) = qqbarxgH_amp_scalars(t, s, u, mh2, mq2, 2, 0, 5, 0, 0, 0, BI, CI, DI)
        amp(7) = qqbarxgH_amp_scalars(u, t, s, mh2, mq2, 3, 0, 6, 0, 0, 0, BI, CI, DI)

     case(iloop_lm_scalar)

        ! Infinite mass limit for scalar particles in loops 
        amp(1) = -one / three * sqrt(two * s * t * u) * s**2 / (s * t * u)
        amp(2) = one / three * sqrt(two * s * t * u) * mh2**2 / (s * t * u)
        amp(3) = -one / three * sqrt(two * s * t * u) * t**2 / (s * t * u)
        amp(4) = -one / three * sqrt(two * s * t * u) * u**2 / (s * t * u)
        amp(5) = -(s - mh2) / 24._dp
        amp(6) = -(t - mh2) / 24._dp
        amp(7) = -(u - mh2) / 24._dp
        
     case default
        call wae_error('hboson_Msquared',':unrecognised loop type',intval=iloop)
     end select

      do i = 1, 7
        amp_imag(i) = amp_imag(i) + dreal(amp(i)) * yukawa 
        amp_real(i) = amp_real(i) - dimag(amp(i)) * yukawa 
      end do

    end do
    
    ! Square and add prefactors 

    ! Process: gg -> gH 
    wtgg = 0.03125_dp * float(nc * (nc**2 - 1)) / mw2 &
           & * (amp_real(1)**2 + amp_imag(1)**2 + amp_real(2)**2 + amp_imag(2)**2 &
           & + amp_real(3)**2 + amp_imag(3)**2 + amp_real(4)**2 + amp_imag(4)**2) * fluxgg

    ! Process: q qbar -> gH 
    wtqq = 16._dp * (u**2 + t**2) / (u + t)**2 / (s * mw2) &
           & * (amp_real(5)**2 + amp_imag(5)**2) * fluxqq

    ! Crossed process: qg -> qH, eq. (A.21) 
    wtqg = -16._dp * (u**2 + s**2) / (u + s)**2 / (t * mw2) &
           & * (amp_real(6)**2 + amp_imag(6)**2) * fluxgq

    ! Crossed process: g qbar -> qbar H, eq. (A.22) 
    wtgq = -16._dp * (s**2 + t**2) / (s + t)**2 / (u * mw2) &
           & * (amp_real(7)**2 + amp_imag(7)**2) * fluxgq

  end subroutine hboson_Msquared

  !=======================================================================================
! Matrix element squared for CP-odd pp -> H + jet 
! Subprocesses: 
!    gg -> gH 
!    qqbar -> gH 
!    qg -> qH (obtained by permutation of the qqbar amplitude) 
!    gqbar -> qbar H (obtained by permutation of the qqbar amplitude) 
! Modified code from Herwig 6.521 
  subroutine hboson_cpodd_Msquared(s, t, u, iloop_in, mass_in, yukawa_in,  wtqq, wtqg, wtgq, wtgg)
    real(dp), intent(in) :: s, t, u
    integer, intent(in) :: iloop_in(:)
    real(dp), intent(in) :: mass_in(:), yukawa_in(:)
    real(dp), intent(out) :: wtqq, wtqg, wtgq, wtgg

    complex(dp) :: BI(4), CI(7), DI(3), amp(7) 

    real(dp) :: mw2, rngluon, rnquark, fluxgg, fluxgq, fluxqq, mq2, yukawa 
    real(dp) :: amp_imag(7), amp_real(7)
    integer :: i, imass, nc, iloop
    
    mw2 = mw_in**2  

    ! Spin and colour flux factors plus enhancement factor
    nc = three ! Number of colors 
    rngluon = one / float(nc**2 - 1) 
    rnquark = one / float(nc) 
    fluxgg = 0.25_dp * rngluon**2
    fluxgq = 0.25_dp * rngluon * rnquark
    fluxqq = 0.25_dp * rnquark**2
      
    amp_imag = zero
    amp_real = zero
    
    do imass = 1, size(iloop_in)
    
      mq2 = mass_in(imass)**2
      if (mq2 == zero) cycle

      yukawa = yukawa_in(imass)
      
      ! Selection for loop quarks 
      iloop = iloop_in(imass)

      ! Evaluate scalar loop integrals 
         
      ! Two-point function 
      BI(1) = B0(s, mq2) 
      BI(2) = B0(t, mq2) 
      BI(3) = B0(u, mq2) 
      BI(4) = B0(mh2, mq2) 
      BI(1) = BI(1) - BI(4) ! B1(s) 
      BI(2) = BI(2) - BI(4) ! B1(t) 
      BI(3) = BI(3) - BI(4) ! B1(u) 

      ! Three-point function 
      CI(1) = C0(s, mq2) 
      CI(2) = C0(t, mq2) 
      CI(3) = C0(u, mq2) 
      CI(7) = C0(mh2, mq2) 
      CI(4) = (s * CI(1) - mh2 * CI(7)) / (s - mh2) ! C1(s) 
      CI(5) = (t * CI(2) - mh2 * CI(7)) / (t - mh2) ! C1(t) 
      CI(6) = (u * CI(3) - mh2 * CI(7)) / (u - mh2) ! C1(u) 
       
      ! Four-point function 
      DI(1) = D0(u, t, mh2, mq2) 
      DI(2) = D0(s, u, mh2, mq2) 
      DI(3) = D0(s, t, mh2, mq2) 

      ! Compute the complex amplitudes
      select case(iloop)
         
      case(iloop_fm_fermion)

         ! Exact result for loops 
         amp(1) = ggxgH_cpodd_pppamp_fermions(s, t, u, mh2, mq2, 1, 2, 3, 4, 5, 6, BI, CI, DI)
         amp(2) = ggxgH_cpodd_ppmamp_fermions(s, t, u, mh2, mq2, 1, 2, 3, 0, 0, 0, BI, CI, DI)
         amp(3) = ggxgH_cpodd_pppamp_fermions(t, s, u, mh2, mq2, 2, 1, 3, 5, 4, 6, BI, CI, DI)
         amp(4) = ggxgH_cpodd_mppamp_fermions(u, t, s, mh2, mq2, 3, 2, 1, 6, 5, 4, BI, CI, DI)
         amp(5) = qqbarxgH_cpodd_amp_fermions(s, t, u, mh2, mq2, 1, 0, 4, 0, 0, 0, BI, CI, DI)
         amp(6) = qqbarxgH_cpodd_amp_fermions(t, s, u, mh2, mq2, 2, 0, 5, 0, 0, 0, BI, CI, DI)
         amp(7) = qqbarxgH_cpodd_amp_fermions(u, t, s, mh2, mq2, 3, 0, 6, 0, 0, 0, BI, CI, DI)
         
      case default
         call wae_error('hboson_Msquared',':unrecognised loop type',intval=iloop)
      end select

      do i = 1, 7
         amp_imag(i) = amp_imag(i) + dreal(amp(i)) * yukawa 
         amp_real(i) = amp_real(i) - dimag(amp(i)) * yukawa 
      end do

    end do
    
    ! Square and add prefactors 

    ! Process: gg -> gH 
    wtgg = 0.25_dp * float(nc * (nc**2 - 1)) / mw2 &
           & * (amp_real(1)**2 + amp_imag(1)**2 + amp_real(2)**2 + amp_imag(2)**2 &
           & + amp_real(3)**2 + amp_imag(3)**2 + amp_real(4)**2 + amp_imag(4)**2) * fluxgg

    ! Process: q qbar -> gH 
    wtqq = two*float(nc**2-1) * (u**2 + t**2) / (s * mw2) &
           & * (amp_real(5)**2 + amp_imag(5)**2) * fluxqq

    ! Crossed process: qg -> qH, eq. (A.21) 
    wtqg = -two*float(nc**2-1) * (u**2 + s**2)  / (t * mw2) &
           & * (amp_real(6)**2 + amp_imag(6)**2) * fluxgq

    ! Crossed process: g qbar -> qbar H, eq. (A.22) 
    wtgq = -two*float(nc**2-1) * (s**2 + t**2) / (u * mw2) &
           & * (amp_real(7)**2 + amp_imag(7)**2) * fluxgq

  end subroutine hboson_cpodd_Msquared

!=======================================================================================
! Eq. (2.4) 

  function f(tau) result(res)

    real(dp), intent(in) :: tau
    complex(dp) :: res
  
    if(tau < one) then
      res = asin(sqrt(tau))**2
    else
      res = -one / four * (log((one + sqrt(one - one / tau)) / (one - sqrt(one - one / tau))) & 
            & - dcmplx(zero, one) * pi)**2
    endif
  
  end function f

!=======================================================================================

  function F12(t) result(res)
    real(dp), intent(in) :: t
    complex(dp) :: res

    res = two * (t + (t - one) * f(t)) / t**2

  end function F12

  function F12_cpodd(t) result(res)
    real(dp), intent(in) :: t
    complex(dp) :: res

    res = two*f(t)/t

  end function F12_cpodd

!=======================================================================================

  function F0(t) result(res)

    real(dp), intent(in) :: t
    complex(dp) :: res

    res = -(t - f(t)) / t**2

  end function F0

!=======================================================================================

  function born_ampsquare(iloop_in, mass_in, yukawa_in) result(res)

    integer, intent(in) :: iloop_in(:)
    real(dp), intent(in) :: mass_in(:), yukawa_in(:)
    real(dp) :: res  

    real(dp) :: tau
    integer :: i
    complex(dp) :: amp

    amp = zero

    do i = 1, size(iloop_in)

      if (mass_in(i) == zero) cycle

      tau = mh2 / (four * mass_in(i)**2)
   
      select case(iloop_in(i))
        case(iloop_fm_fermion) ! fermion
          if (cpodd) then
             amp = amp + F12_cpodd(tau) * yukawa_in(i)
          else
             amp = amp + F12(tau) * yukawa_in(i) 
          end if
       case(iloop_lm_fermion)
          if (cpodd) then
             amp = amp + two * yukawa_in(i)
          else
             amp = amp + four / three * yukawa_in(i) 
          end if
       case(iloop_fm_scalar) ! scalar 
          if (cpodd) call wae_error('born_ampsquare',&
               &'CP-odd Higgs with scalar loops not yet implemented')
          amp = amp + half * F0(tau) * yukawa_in(i) 
       case(iloop_lm_scalar)
          if (cpodd) call wae_error('born_ampsquare',&
               &'CP-odd Higgs with scalar loops not yet implemented')
          amp = amp + one / three * yukawa_in(i) 
       case default
          call wae_error('born_ampsquare','unknown iloop',intval=iloop_in(i))
      end select

    end do

    ! Divide off the large-mt limit
    amp = amp / (four / three)

    res = real(amp)**2 + aimag(amp)**2

  end function born_ampsquare

!=======================================================================================
! Cross-section for bbH 

  function bbH_cross_section(lumiqqbar) result(res)

    real(dp) :: lumiqqbar
    real(dp) :: res

    res = pi / 6._dp / higgs_vev_in**2 / mh2 
    res = res * lumiqqbar 
    res = res * invGev2_to_nb 

  end function bbH_cross_section

!=======================================================================================
! Matrix element squared for bbH  

  subroutine bbH_Msquared(s, t, u, wtqq, wtqg, wtgq, wtgg)

    real(dp), intent(in) :: s, t, u

    real(dp), intent(out) :: wtqq, wtqg, wtgq, wtgg

    wtgg = zero
    wtqq = cf * (s**2 + mh2**2) / (t * u)
    wtqg = tr * (t**2 + mh2**2) / (-s * u)
    wtgq = tr * (u**2 + mh2**2) / (-t * s)

  end subroutine bbH_Msquared

!=======================================================================================

end module hboson


