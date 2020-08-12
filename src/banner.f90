!========================================================
!--------------------------------------------------------
! Module containing relevant subroutines for the welcome 
! banner, help message, and version number 
!--------------------------------------------------------
!========================================================

module banner

  implicit none

  private

  ! Version number 
  character(len=5), parameter :: version_number = '1.0.0' 

  ! Publically available subroutines  
  public :: print_help_message, print_welcome_banner, print_version_number 

contains
  
!=======================================================================================
! Print the help message (invoked with the --help or -h command options) 

  subroutine print_help_message(idev)

    use user_interface

    integer, intent(in) :: idev

    write(idev,'(a)') 'Usage: h1jet [options]'
    write(idev,'(a)') 'Calculate the transverse momentum distribution for Higgs + jet production.'
    write(idev,'(a)') 
    write(idev,'(a)') 'Options:'
    write(idev,'(a)') '  -o, --out <file>       Direct output to <file>'
    write(idev,'(a)') '  --collider <arg>       Specify collider type (pp, ppbar), default = pp'
    write(idev,'(a)') '  --roots <value>        Center-of-mass energy in [GeV], default = 13 TeV'
    write(idev,'(a)') '  --proc <arg>           Specify the process (H, bbH, Z, user)'
    write(idev,'(a)') '                         H    : pp -> Higgs + jet (default)' 
    write(idev,'(a)') '                         bbH  : b bbar -> Higgs + jet' 
    write(idev,'(a)') '                         Z    : pp -> Z + jet' 
    write(idev,'(a)') '                         user : Amplitude from custom code, see manual for implementation' 
    write(idev,'(a)') '  --approx <arg>         Specify loop approximation (none, sml, iml)'
    write(idev,'(a)') '                         none : Exact result for loops (default)' 
    write(idev,'(a)') '                         sml  : Small mass limit for quarks in the loops' 
    write(idev,'(a)') '                         iml  : Infinite mass limit for quarks in the loops' 
    write(idev,'(a)') '  --pdf_name <arg>       Specify the PDF set name from LHAPDF, default = MSTW2008nlo68cl'
    write(idev,'(a)') '  --pdf_mem <value>      Integer value specifying the PDF member, default = 0'
    write(idev,'(a)') '  --scale_strategy <arg> Set the scale strategy, i.e. set the dynamic muR = muF value (M, HT, MT)'
    write(idev,'(a)') '                         M  : muR = muF = M' 
    write(idev,'(a)') '                         HT : muR = muF = pT + sqrt(pT^2 + M^2) (default)' 
    write(idev,'(a)') '                         MT : muR = muF = sqrt(pT**2 + M^2)' 
    write(idev,'(a)') '                         M is --mH for proc = H / bbH, --mZ for proc = Z, --mass for proc = user' 
    write(idev,'(a)') '  --xmur <value>         Factor for the renormalization scale, muR = muR * xmur, default = 0.5'
    write(idev,'(a)') '  --xmuf <value>         Factor for the factorization scale, muF = muF * xmuf, default = 0.5'
    write(idev,'(a)') '  --nbins <value>        Number of histogram bins in the output, default = 400'
    write(idev,'(a)') '  --log                  Enables logarithmic x-axis of histogram, i.e. logarithmic bins and pT'
    write(idev,'(a)') '  --ptmin <value>        Minimum pT value, default = 0 (set to nonzero if --log is enabled)'
    write(idev,'(a)') '  --ptmax <value>        Maximum pT value, default = 4000 GeV'
    write(idev,'(a)') '  --accuracy <value>     The desired Monte Carlo integration accuracy, default = 0.001'
    write(idev,'(a)') '  --cppodd               Toggle for pseudoscalar Higgs' 
    write(idev,'(a)') '  --mass <value>         Relevant mass in user process for muR/muF [GeV], default = 125 GeV' 
    write(idev,'(a)') '  --mH <value>           Higgs mass [GeV], default = 125 GeV'
    write(idev,'(a)') '  --mZ <value>           Z boson mass [GeV], default = 91.1876 GeV'
    write(idev,'(a)') '  --mW <value>           W boson mass [GeV], default = 80.385 GeV'
    write(idev,'(a)') '  --gf <value>           Fermi coupling constant [GeV^(-2)], default = 0.116638e-4 GeV^(-2)'
    write(idev,'(a)') '  --mt <value>           Top quark mass [GeV], default = 173.5 GeV'
    write(idev,'(a)') '  --mb <value>           On-shell bottom quark mass [GeV], default = 4.65 GeV'
    write(idev,'(a)') '  --yt <value>           Top Yukawa factor [GeV], default = 1'
    write(idev,'(a)') '  --yb <value>           Bottom Yukawa factor [GeV], default = 1'
    write(idev,'(a)') '  --mtp <value>          Top partner mass [GeV], default = 0 GeV'
    write(idev,'(a)') '  --ytp <value>          Top partner Yukawa factor, default = 1'
    write(idev,'(a)') '  --mst <value>          MSSM stop mass [GeV], default = 0 GeV'
    write(idev,'(a)') '  --delta <value>        MSSM stop mass separation [GeV], default = 0 GeV'
    write(idev,'(a)') '  --sth2 <value>         Stop (or top-partner) mixing angle, sin^2(theta), default = 0'
    write(idev,'(a)') '  --tbeta <value>        Ratio of VEVs in MHDMs, tan(beta), default = 0'
    write(idev,'(a)') '  --mbmb <value>         MSbar bottom mass mb(mb) [GeV], default = 4.18 GeV'
    write(idev,'(a)') 

    call user_help_message(idev) 

    write(idev,'(a)') '  -h, --help             Print this help message and exit'
    write(idev,'(a)') '  -v, --version          Print the version number and exit'
    write(idev,'(a)') ! Blank space for cleaner output 

    stop 

  end subroutine print_help_message

!=======================================================================================
! Print the welcome banner at the start of the program 

  subroutine print_welcome_banner(idev)

    integer, intent(in) :: idev
    
    write(idev,'(a)') '============================================================'
    write(idev,'(a)') '                    H1jet version '//trim(version_number) 
    write(idev,'(a)') '    Tool for the fast and accurate calculation of the pT    '
    write(idev,'(a)') ' distribution in Higgs + jet production at hadron colliders '
    write(idev,'(a)') '============================================================'

    write(idev,*) ! Blank space for cleaner output 

  end subroutine print_welcome_banner

!=======================================================================================

  subroutine print_version_number(idev) 

    integer, intent(in) :: idev

    write(idev,'(a)') 'H1jet version '//trim(version_number) 

    stop 

  end subroutine print_version_number 

!=======================================================================================

end module banner


