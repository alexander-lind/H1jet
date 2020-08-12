!========================================================
!--------------------------------------------------------
! Module implementing the scalar loop integrals  
! for two-, three-, and four-point functions 
! appearing in the amplitudes for Higgs + jet production 
!--------------------------------------------------------
!========================================================

module scalar_integrals

  use types 
  use consts_dp
  use hoppet_v1 
  
  implicit none

  complex(dp), parameter :: ii = dcmplx(zero, one)
  complex(dp), parameter :: ipi = ii * pi

contains

!=======================================================================================
! Scalar two-point function 
! See eq. (B.3) and (B.4) in Baur and Glover, Nucl.Phys. B339 (1990) 38-66 

  function B0(s, mq2) result(res)
    
    real(dp), intent(in) :: s, mq2

    complex(dp):: res, z
    real(dp) :: rat
    complex(dp) :: HPL1
    
    rat = (four * mq2) / s

    if (rat > one) then
       z = half * (one + ii * sqrt(rat - one))
       res = two - ii * sqrt(rat - one) * HPL1(0, z / (z - one))
    else
       z = half * (one + sqrt(one - rat))
       res = two + sqrt(one - rat) * HPL1(0, (z - one) / z)
    end if

  end function B0

!=======================================================================================
! Scalar three-point function with two massless external lines 
! See eq. (B.6) in Baur and Glover, Nucl.Phys. B339 (1990) 38-66 

 function C0(s, mq2) result(res)
    
    real(dp), intent(in) :: s, mq2
    
    complex(dp) :: res, z
    real(dp) :: rat
    complex(dp) :: HPL2
   
    rat = (four * mq2) / s 

    if (rat > one) then
       z = half * (one + ii * sqrt(rat - one))
    else
       z = half * (one + sqrt(one - rat))
    end if

    res = HPL2(0, 0, (z - one) / z) / s
    
 end function C0

!=======================================================================================
! Scalar four-point function with three massless and one massive external line 
! See eq. (B.8) and (B.9) in Baur and Glover, Nucl.Phys. B339 (1990) 38-66 

  function D0(s, t, m2, mq2) result(res)
    
    real(dp), intent(in) :: s, t, m2, mq2

    complex(dp) :: res 
    real(dp) :: u 

    u = m2 - s - t 

    res = CI2(s, t, u, s, mq2) + CI2(s, t, u, t, mq2) - CI2(s, t, u, m2, mq2) 

    res = two * res / (t * s) / sqrt(one + four * mq2 * u / (t * s)) 

  end function D0

  ! Ensure a +i*epsilon prescription for CHAPLIN 
  ! by using either the reflection properties or simply H(0,-1,z) 
  function hpl_prescription_helper(z) result(res)

    complex(dp), intent(in) :: z

    complex(dp) :: res
    complex(dp) :: z_in 
    complex(dp) :: HPL1, HPL2
    complex(dp) :: tmp1, tmp2, tmp3 

    ! Set imaginary part to zero if it's much smaller than the real part 
    z_in = z 
    if (abs(aimag(z)/real(z)) <= 1.0e-60_dp) then 
        z_in = dcmplx(real(z), zero)  
    endif 

    ! Each of the three tmp's are equal to Li2(z) 
    ! However, they use different regions of the argument so that CHAPLIN 
    ! will use different methods for evaluation  
    tmp1 = -HPL2(0,1,one/z_in) - pi**two/6._dp - HPL2(0,0,-z_in)
    tmp2 = -HPL2(0,-1,-z_in)
    tmp3 = -HPL2(0,1,one-z_in) + pi**two/6._dp + HPL1(1,-z_in+one) * HPL1(-1,-z_in) 

    ! We will now check which of the methods give the smallest imaginary part 
    ! This does not matter if the imaginary part is non-zero
    ! However, in the case that we have a infinitesimal imaginary part this is necessary 
    if (abs(aimag(tmp1)) < abs(aimag(tmp2))) then
      if (abs(aimag(tmp1)) < abs(aimag(tmp3))) then
        res = tmp1
      else 
        res = tmp3 
      end if
    else 
      if (abs(aimag(tmp2)) < abs(aimag(tmp3))) then
        res = tmp2
      else 
        res = tmp3 
      end if
    end if 

  end function hpl_prescription_helper

  ! Eq. (B.9) in Baur and Glover 
  function CI2(s, t, u, v, mq2) result(res)

    real(dp), intent(in) :: s, t, u, v, mq2 

    complex(dp) :: res
    real(dp) :: xp, xm, rat
    complex(dp) :: y, HPL1, HPL2
    complex(dp) :: z1, z2, z3, z4
    complex(dp) :: sp1, sp2, sp3, sp4

    ! Initialise all variables
    y = (zero, zero)
    z1 = (zero, zero); z2 = (zero, zero); z3 = (zero, zero); z4 = (zero, zero)

    rat = four * mq2 / v

    if (v > 0 .and. v < (four * mq2)) then
      y = half * (one + ii * dsqrt(rat - one))
    else
      y = half * (one + dsqrt(one - rat))
    end if

    xp = half * (one + dsqrt(one + four * mq2 * u / (t * s)))
    xm = one - xp

    z1 = xm / (xm - y)
    z2 = xp / (xp - y)
    z3 = xm / (y - xp)
    z4 = xp / (y - xm)

    ! Spence's Functions (dilogs) 
    sp1 = HPL2(0, 1, z1)
    sp2 = HPL2(0, 1, z2)
    sp3 = HPL2(0, 1, z3)
    sp4 = HPL2(0, 1, z4)

    ! Ensure a +i*epsilon prescription for CHAPLIN 
    if (v < zero) then 
        sp2 = hpl_prescription_helper(z2) 
        sp3 = hpl_prescription_helper(z3) 
    else if (v > (four * mq2)) then 
        sp1 = hpl_prescription_helper(z1) 
        sp4 = hpl_prescription_helper(z4) 
    endif 

    res = sp1 - sp2 + sp3 - sp4 

    res = res + log(dcmplx(-xm / xp, zero)) * log(dcmplx(one - v * xm * xp / mq2, zero)) 

  end function CI2

!=======================================================================================

end module scalar_integrals


