!========================================================
!--------------------------------------------------------
! Module handling various user inputs  
!--------------------------------------------------------
!========================================================

module input

  use hoppet_v1 
  use ew_parameters 
  use mass_helper
  use sub_defs_io
  use hboson
  use common_vars 

  implicit none

  private

  character(len=8) :: proc 

  public :: handle_input, gen_momenta, set_scale, print_settings 

contains

!======================================================================================= 
! Handle all user input 

  subroutine handle_input 

    ! Print commandline options
    write(idev,'(a)') 'Provided command options: '
    write(idev,'(a)') '# '//trim(command_line()) 
    write(idev,*) ! Blank line for nicer output 

    ! Process 
    proc = trim(string_val_opt('--proc', 'H'))
    iproc = proc_to_id(proc)

    ! Set the alpha_s power depending on the process 
    if (iproc == id_H) then 
      ! Result multiplied by factor alpha_s^(as_pow + 1) = alpha_s^3 
      as_pow = 2 
    else if (iproc == id_user) then 
      ! Result multiplied by factor alpha_s^0 = 1, since the 
      ! appropriate alpha_s factor is already included 
      as_pow = -1 
    end if 

    ! Centre-of-mass collision energy 
    roots = dble_val_opt('--roots', 13000._dp) 
    
    ! Set up collider type 
    collider = trim(string_val_opt('--collider', 'pp'))
    if (collider /= 'pp' .and. collider /= 'ppbar') then
      call wae_error('handle_input', 'Unrecognised collider type &
                &-- Please specify collider as pp or ppbar')
    end if 

    ! Set SM and EW parameters if changed (from ew_parameters.f90) 
    ! G_mu scheme 
    mz_in = dble_val_opt('--mZ', mz)
    mw_in = dble_val_opt('--mW', mw)
    GF_GeVm2_in = dble_val_opt('--gf', GF_GeVm2)
    sinwsq_in = one - (mw_in**2) / (mz_in**2) 
    higgs_vev_in = one / dsqrt(dsqrt(two) * GF_GeVm2_in) 
    call set_gv2 

    ! Set heavy quark masses if changed (from mass_helper.f90) 
    mt_in = dble_val_opt('--mt', mt)
    mb_in = dble_val_opt('--mb', mb)
    
    ! Set PDF set 
    pdf_name = trim(string_val_opt('--pdf_name', 'MSTW2008nlo68cl')) !//'.LHgrid'
    pdf_mem  = int_val_opt('--pdf_mem', 0)

    ! Approximation for quark loops 
    approx = trim(string_val_opt('--approx', 'none'))
    if (approx /= 'none' .and. approx /= 'sml' .and. approx /= 'iml') then
      call wae_error('handle_input', 'Unrecognised approximation type: ', approx)
    end if 

    call set_masses_and_couplings

    ! Set up factorisation and renormalisation scales 
    scale_strategy = trim(string_val_opt('--scale_strategy', 'HT'))

    muF = set_scale(scale_strategy, zero)
    muR = muF

    xmur = dble_val_opt('--xmur', 0.5_dp)
    xmuf = dble_val_opt('--xmuf', 0.5_dp)

    muF = muF * xmuf
    muR = muR * xmur

    ! Histogram options 
    nbins = int_val_opt('--nbins', 400)
    ptmin = dble_val_opt('--ptmin', zero)
    ptmax = dble_val_opt('--ptmax', 4e3_dp)

    if (log_val_opt('--log') .and. (ptmin == zero)) then 
      call wae_error('handle_input', 'Set minimum pT value (with --ptmin) &
                    &when using --log')
    end if

    ! Set the desired Monte Carlo accuracy 
    accuracy = dble_val_opt('--accuracy', 1e-3_dp)

  end subroutine handle_input
  
!======================================================================================= 
! Convert --proc command option argument to internal id 

  function proc_to_id(proc) result(res)

    character(len=*) :: proc
    integer :: res
    
    if (proc == 'H') then
      res = id_H
    else if (proc == 'Z') then
      res = id_Z
    else if (proc == 'bbH') then
      res = id_bbH
    else  if (proc == 'user') then
      res = id_user
    else
      call wae_error('proc_to_id', 'Unrecognised process ', proc)
    end if

  end function proc_to_id
  
!=======================================================================================
! Generate initial momenta 

  subroutine gen_momenta(y, p)

    real(dp), intent(in) :: y
    real(dp), intent(out) :: p(4,4)

    real(dp) :: E, pz, Eb, Ebeam

    E = pt * cosh(y) 
    pz = pt * sinh(y)
    p(3,:) = (/zero, pt, pz, E/)
    Eb = sqrt(M**2 + E**2)
    p(4,:) = (/zero, -pt, -pz, Eb/)

    Ebeam = E + Eb

    p(1,:) = Ebeam / two * (/zero, zero, one, one/)
    p(2,:) = Ebeam / two * (/zero, zero, -one, one/)
  
    !s = mh**2 + two * pt * cosh(y) * sqrt(pt**2 * cosh(y)**2 + mh**2) &
    !    & + two * pt**2 * cosh(y)**2 

  end subroutine gen_momenta

!=======================================================================================
! Function to set the factorisation and renormalisation scales according to the 
! specified scale strategy and given pT 

  function set_scale(scale_strategy, ptval) result(res)

    character(len=*), intent(in) :: scale_strategy
    real(dp), intent(in) :: ptval
    real(dp) :: res

    if (scale_strategy == 'M') then
      res = M
    else if (scale_strategy == 'HT') then
      res = ptval + sqrt(ptval**2 + M**2)
    else if (scale_strategy == 'MT') then
      res = sqrt(ptval**2 + M**2)
    else
      call wae_error('set_scale', 'Unrecognised scale strategy: ', scale_strategy)
    end if

  end function set_scale

!=======================================================================================
! Function to set the relevant masses and parameters depending on the selected process 

  subroutine set_masses_and_couplings
    
    select case(iproc)

      ! pp -> H + jet 
      case(id_H)

        ! Check for CP odd Higgs 
        cpodd = log_val_opt('--cpodd') 
        
        ! Standard parameters 
        M = dble_val_opt('--mH', mh) ! Higgs mass 
        yt = dble_val_opt('--yt', one) ! Top Yukawa factor 
        yb = dble_val_opt('--yb', one) ! Bottom Yukawa factor 
        sth2 = dble_val_opt('--sth2', zero) ! Stop/Top-partner mixing angle 

        ! Top partner mass and Yukawa coupling 
        mtp = dble_val_opt('--mtp', zero)
        ytp = dble_val_opt('--ytp', one)
        
        ! Stop squark mass and SUSY parameters 
        mst1 = dble_val_opt('--mst', zero) ! Stop mass 
        delta = dble_val_opt('--delta', zero) ! SUSY fine-tuning parameter 
        mst2 = sqrt(mst1**2 + delta**2)
        tbeta = dble_val_opt('--tbeta', zero) ! Ratio of VEV's in MHDMs, tan(beta) 
        c2beta = (one-tbeta**2) / (one + tbeta**2) ! cos(2 * beta) 

        if ((mst1 * mtp) /= zero) then 
          call wae_error('set_masses_and_couplings', 'Simoultaneous top partners &
                    &and stops not allowed') 
        end if 

        if (mtp /= zero) then
          ! Include top partner in quark loops 
          allocate(mass_array(3), yukawa(3), iloop_array(3))
          yt = yt * (one - sth2)
          ytp = ytp * sth2
          mass_array = (/mt_in, mb_in, mtp/)
          yukawa = (/yt, yb, ytp/)
          select case(approx)
            case('sml')
              ! Small mass limit 
              iloop_array = (/0, 0, 0/)
            case('iml')
              ! Infinite mass limit 
              iloop_array = (/2, 2, 2/)
            case default 
              ! Exact result -- equivalent to approx=none 
              iloop_array = (/1, 1, 1/)
          end select
        else if (mst1 /= zero) then
          ! Include SUSY stop squark in loops 
          cth2 = one - sth2
          s2th2 = four * sth2 * cth2
          alpha1 = mz_in**2 / mt_in**2 * c2beta * (one - four / three * sinwsq_in)
          alpha2 = mz_in**2 / mt_in**2 * c2beta * four / three * sinwsq_in
          yst1 = (mt_in / mst1)**2 * (alpha1 * cth2 + alpha2 * sth2 + two - delta**2 / two / mt_in**2 * s2th2)
          yst2 = (mt_in / mst2)**2 * (alpha1 * sth2 + alpha2 * cth2 + two + delta**2 / two / mt_in**2 * s2th2)
          allocate(mass_array(4), yukawa(4), iloop_array(4))
          mass_array = (/mt_in, mb_in, mst1, mst2/)
          yukawa = (/yt, yb, yst1, yst2/)
          select case(approx)
            case('sml')
              ! Small mass limit 
              iloop_array = (/0, 0, 3, 3/)
            case('iml')
              ! Infinite mass limit 
              iloop_array = (/2, 2, 4, 4/)
            case default 
              ! Exact result -- equivalent to approx=none 
              iloop_array = (/1, 1, 3, 3/)
          end select 
        else
          ! SM, only bottom and top in quark loops 
          allocate(mass_array(2), yukawa(2), iloop_array(2))
          mass_array = (/mt_in, mb_in/)
          yukawa = (/yt, yb/)
          select case(approx)
            case('sml')
              ! Small mass limit 
              iloop_array = (/0, 0/)
            case('iml')
              ! Infinite mass limit 
              iloop_array = (/2, 2/)
            case default 
              ! Exact result -- equivalent to approx=none 
              iloop_array = (/1, 1/)
          end select 
        end if

      ! pp -> Z + jet   
      case(id_Z)

        M = mz_in
      
      ! bbbar -> H + jet 
      case(id_bbH)

        M = dble_val_opt('--mH', mh) ! Higgs mass 
        allocate(mass_array(1), yukawa(1))
        mass_array(1) = dble_val_opt('--mbmb', 4.18_dp)

      ! User specified process 
      case (id_user)

        ! Relevant mass for the user specified amplitude 
        ! Important to set for correct scales muR and muF 
        M = dble_val_opt('--mass', mh) 

      case default

        call wae_error('proc_to_id', 'Unrecognised process ', proc)

    end select

  end subroutine set_masses_and_couplings

!=======================================================================================
! Print the settings from user input 

  subroutine print_settings 

    write(idev,*) ! Blank line for nicer output 
    write(idev,*) 'collider       = ', collider
    write(idev,*) 'roots(GeV)     =', roots
    write(idev,*) 'process        = ', proc
    if (iproc == id_user) then 
      write(idev,*) 'model          = User-defined from included custom code' 
    else if (mtp /= zero) then 
      write(idev,*) 'model          = SM with top partners' 
    else if (mst1 /= zero) then 
      write(idev,*) 'model          = SUSY' 
    else if (cpodd) then 
      write(idev,*) 'model          = Pseudoscalar Higgs'
    else 
      write(idev,*) 'model          = SM' 
    end if
    write(idev,*) 'mass           =', M
    write(idev,*) 'muR            =', muR
    write(idev,*) 'muF            =', muF
    write(idev,*) 'alphas(muR)    =', alphas
    write(idev,*) 'scale strategy = ', scale_strategy 
    if (iproc /= id_user) then 
      write(idev,*) 'masses         =', mass_array
    end if 
    write(idev,*) 'tan(beta)      =', tbeta
    write(idev,*) 'sin^2(theta)   =', sth2
    write(idev,*) 'sigma0 [nb]    =', sigma0 
    write(idev,*) 'log            =', log_val_opt('--log')
    write(idev,*) ! Blank line for nicer output 

  end subroutine print_settings 

!=======================================================================================

end module input 


