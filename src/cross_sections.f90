!========================================================
!--------------------------------------------------------
! Module containing general functions to handle the 
! calculation of the cross-section  
!--------------------------------------------------------
!========================================================

module cross_sections

  use hoppet_v1 
  use pdfs_tools
  use ew_parameters 
  use mass_helper
  use vboson
  use hboson
  use input 
  use user_interface 
  use common_vars 
  
  implicit none

  real(dp), pointer :: lumi_gg(:), lumi_qg(:), lumi_gq(:), lumi_qqbar(:)

contains

!=======================================================================================
! Dot product between two four-vectors 

  function dot(p, q) result(res)

    real(dp), intent(in) :: p(4), q(4)
    real(dp) :: res

    res = p(4) * q(4) - dot_product(p(:3), q(:3))

  end function dot

!=======================================================================================
! Process-dependent Electroweak Prefactor 

  subroutine h1jet_prefactor(factor)

    real(dp), intent(out) :: factor

    select case(iproc)
      case (id_H)
        ! alpha_W = g_W^2 / (4 * pi) 
        ! and M_W = (1/2) * g_W * vev 
        factor = mw_in**2 / (pi * higgs_vev_in**2)
      case (id_Z)
        factor = 8.0_dp * pi / three * sqrt(two) * GF_GeVm2_in * mz_in**2  
      case (id_bbH) 
        factor = four * pi / higgs_vev_in**2 / three
      case (id_user)
        ! EW Prefactor already included in user code 
        factor = 1.0_dp 
      case default
        call wae_error('h1jet_prefactor', 'Unrecognised process')
    end select
      
    !factor = factor * invGev2_to_nb 
    factor = factor * invGev2_to_fb 

  end subroutine h1jet_prefactor

!=======================================================================================
! Calculate the born-level total cross-section 

  function cross_section(lumi_gg, lumi_qg, lumi_gq, lumi_qqbar, tau) result(res)

    real(dp), intent(in) :: lumi_gg(0:), lumi_qg(0:), lumi_gq(0:), lumi_qqbar(0:)
    type(gdval), intent(in) :: tau

    real(dp) :: res

    select case(iproc)
      case (id_H)
        mh2 = M**2
        res = hboson_cross_section(lumi_gg .atx. tau, iloop_array, mass_array, yukawa)
      case (id_Z)
        res = vboson_cross_section(lumi_qqbar .atx. tau)
      case (id_bbH)
        mh2 = M**2
        res = bbH_cross_section(lumi_qqbar .atx. tau)
      case (id_user) 
        res = user_cross_section(lumi_gg .atx. tau, lumi_qg .atx. tau, lumi_gq .atx. tau, lumi_qqbar .atx. tau) 
      case default
        call wae_error('cross_section', 'Unrecognised process')
    end select

  end function cross_section

!=======================================================================================
! Calculate the differential cross-section dsigma / dpT dy 
! as a function of rapidity y for the specified process 

  function dsigma_dptdy(y) result(res)
    
    real(dp), intent(in) :: y

    real(dp) :: res
    real(dp) :: p(4,4), rootshat, Eb
    type(gdval) :: tauhat
    real(dp) :: s, t, u
    real(dp) :: lumigg, lumiqg, lumigq, lumiqqbar
    real(dp) :: wtqq, wtqg, wtgq, wtgg
    real(dp) :: jakob

    call gen_momenta(y, p)

    rootshat = two * p(1,4)
    Eb = p(4,4)
    jakob = pt / (8.0_dp * pi * Eb * rootshat**3)
    
    tauhat    = (rootshat**2 / roots**2) .with. grid
    lumigg    = lumi_gg .atx. tauhat
    lumiqg    = lumi_qg .atx. tauhat
    lumigq    = lumi_gq .atx. tauhat
    lumiqqbar = lumi_qqbar .atx. tauhat

    s =  two * dot(p(1,:), p(2,:))
    t = -two * dot(p(1,:), p(3,:))
    u = -two * dot(p(2,:), p(3,:))

    select case(iproc)
      case (id_H)
        if (cpodd) then
          call hboson_cpodd_Msquared(s, t, u, iloop_array, & 
                 & mass_array, yukawa, wtqq, wtqg, wtgq, wtgg)
        else
          call hboson_Msquared(s, t, u, iloop_array, & 
                 & mass_array, yukawa, wtqq, wtqg, wtgq, wtgg)
        end if
      case (id_Z)
        call vboson_Msquared(s, t, u, wtgg, wtqg, wtgq, wtqq)
      case (id_bbH)
        call bbH_Msquared(s, t, u, wtqq, wtqg, wtgq, wtgg)
      case (id_user)
        call user_Msquared(s, t, u, wtqq, wtqg, wtgq, wtgg) 
      case default
        call wae_error('dsigma_dptdy', 'Unrecognised process')
    end select

    res = (wtgg * lumigg + lumiqg * wtqg + lumigq * wtgq + wtqq * lumiqqbar) * jakob 

  end function dsigma_dptdy

!=======================================================================================

end module cross_sections 


