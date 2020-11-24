!========================================================
!--------------------------------------------------------
! Module containing relevant subroutines for the welcome 
! banner, help message, and version number 
!--------------------------------------------------------
!========================================================

module banner

  use common_vars 

  implicit none

  private

  ! Version number 
  character(len=5), parameter :: version_number = '1.0.0' 

  ! Publically available subroutines  
  public :: print_help_message, print_welcome_banner, print_version_number 

contains
  
!=======================================================================================
! Print the help message (invoked with the --help or -h command options) 

  subroutine print_help_message

    use user_interface

    write(stderr,'(a)') 'Usage: h1jet [options]'
    write(stderr,'(a)') 'Calculate the transverse momentum distribution for Higgs + jet production.'
    write(stderr,'(a)') 
    write(stderr,'(a)') 'Options:'
    write(stderr,'(a)') '  -o, --out <file>       Direct output to <file>'
    write(stderr,'(a)') '  --collider <arg>       Specify the collider type (pp, ppbar), default = pp'
    write(stderr,'(a)') '  --roots <value>        Center-of-mass energy in [GeV], default = 13 TeV'
    write(stderr,'(a)') '  --proc <arg>           Specify the process (H, bbH, Z, user)'
    write(stderr,'(a)') '                         H    : pp -> Higgs + jet (default)' 
    write(stderr,'(a)') '                         bbH  : b bbar -> Higgs + jet' 
    write(stderr,'(a)') '                         Z    : pp -> Z + jet' 
    write(stderr,'(a)') '                         user : Amplitude from custom code, see manual for implementation' 
    write(stderr,'(a)') '  --pdf_name <arg>       Specify the PDF set name from LHAPDF, default = MSTW2008nlo68cl'
    write(stderr,'(a)') '  --pdf_mem <value>      Integer value specifying the PDF member, default = 0'
    write(stderr,'(a)') '  --scale_strategy <arg> Set the scale strategy, i.e. set the dynamic muR = muF value (M, HT, MT)'
    write(stderr,'(a)') '                         M  : muR = muF = M' 
    write(stderr,'(a)') '                         HT : muR = muF = pT + sqrt(pT^2 + M^2) (default)' 
    write(stderr,'(a)') '                         MT : muR = muF = sqrt(pT**2 + M^2)' 
    write(stderr,'(a)') '                         M is --mH for proc = H / bbH, --mZ for proc = Z, --mass for proc = user' 
    write(stderr,'(a)') '  --xmur <value>         Factor for the renormalization scale, muR = muR * xmur, default = 0.5'
    write(stderr,'(a)') '  --xmuf <value>         Factor for the factorization scale, muF = muF * xmuf, default = 0.5'
    write(stderr,'(a)') '  --nbins <value>        Number of histogram bins in the output, default = 400'
    write(stderr,'(a)') '  --log                  Enables logarithmic x-axis of histogram, i.e. logarithmic bins and pT'
    write(stderr,'(a)') '  --ptmin <value>        Minimum pT value, default = 0 (set to nonzero if --log is enabled)'
    write(stderr,'(a)') '  --ptmax <value>        Maximum pT value, default = 4000 GeV'
    write(stderr,'(a)') '  --accuracy <value>     The desired Monte Carlo integration accuracy, default = 0.001'
    write(stderr,'(a)') '  --cpodd                Toggle for CP-odd Higgs' 
    write(stderr,'(a)') '  -M, --mass <value>     Relevant mass in user process [GeV], default = 0 GeV' 
    write(stderr,'(a)') '  --mH <value>           Higgs mass [GeV], default = 125 GeV'
    write(stderr,'(a)') '  --mZ <value>           Z boson mass [GeV], default = 91.1876 GeV'
    write(stderr,'(a)') '  --mW <value>           W boson mass [GeV], default = 80.385 GeV'
    write(stderr,'(a)') '  --GF <value>           Fermi coupling constant [GeV^(-2)], default = 0.116638e-4 GeV^(-2)'
    write(stderr,'(a)') '  --mt <value>           Top quark mass [GeV], default = 173.5 GeV'
    write(stderr,'(a)') '  --mb <value>           On-shell bottom quark mass [GeV], default = 4.65 GeV'
    write(stderr,'(a)') '  --yt <value>           Top Yukawa factor [GeV], default = 1'
    write(stderr,'(a)') '  --yb <value>           Bottom Yukawa factor [GeV], default = 1 (0 for CP-odd Higgs)'
    write(stderr,'(a)') '  --mtp <value>          Top partner mass [GeV], default = 0 GeV'
    write(stderr,'(a)') '  --ytp <value>          Top partner Yukawa factor, default = 1'
    write(stderr,'(a)') '  --mst <value>          SUSY stop mass [GeV], default = 0 GeV'
    write(stderr,'(a)') '  --delta <value>        SUSY stop mass separation [GeV], default = 0 GeV'
    write(stderr,'(a)') '  --sth2 <value>         Stop (or top-partner) mixing angle, sin^2(theta), default = 0'
    write(stderr,'(a)') '  --tbeta <value>        Ratio of VEVs of the two SUSY Higgs fields, tan(beta), default = 0'
    write(stderr,'(a)') '  --mbmb <value>         MSbar bottom mass mb(mb) [GeV], default = 4.18 GeV'
    write(stderr,'(a)') '  -i, --in <file>        Include input file with top-partner masses and Yukawas' 
    write(stderr,'(a)') '                         See SM.dat for the SM case, and the manual for more information'
    write(stderr,'(a)') '  --model <arg>          Specify top-partner model' 
    write(stderr,'(a)') '                         M1_5  : Light top-partner transforming as a 1_2/3 of SO(4) and the'
    write(stderr,'(a)') '                                 SM top-bottom doublet embedded in a 5 of SO(5) (default)'
    write(stderr,'(a)') '                         M1_14 : Light top-partner transforming as a 1_2/3 of SO(4) and the'
    write(stderr,'(a)') '                                 SM top-bottom doublet embedded in a 14 of SO(5)'
    write(stderr,'(a)') '                         M4_5  : Light top-partner transforming as a 4_2/3 of SO(4) and the'
    write(stderr,'(a)') '                                 SM top-bottom doublet embedded in a 5 of SO(5)'
    write(stderr,'(a)') '                         M4_14 : Light top-partner transforming as a 4_2/3 of SO(4) and the'
    write(stderr,'(a)') '                                 SM top-bottom doublet embedded in a 14 of SO(5)'
    write(stderr,'(a)') '                         See arXiv:1905.12747 for more information on the models'
    write(stderr,'(a)') '  --imc1 <value>         Imaginary part of the c_1 coefficients, default = 0'
    write(stderr,'(a)') '  -f, --fscale <value>   Top-partner confinement scale [GeV], default = 0 GeV'
    write(stderr,'(a)') 

    call user_help_message 

    write(stderr,'(a)') '  -h, --help             Print this help message and exit'
    write(stderr,'(a)') '  -v, --version          Print the version number and exit'
    write(stderr,'(a)') ! Blank space for cleaner output 

    stop 

  end subroutine print_help_message

!=======================================================================================
! Print the welcome banner at the start of the program 

  subroutine print_welcome_banner(idev)
    integer, intent(in) :: idev
    
    write(idev,'(a)') '============================================================'
    write(idev,'(a)') '                     H1jet version '//trim(version_number) 
    write(idev,'(a)') '    Tool for the fast and accurate calculation of the pT    '
    write(idev,'(a)') ' distribution in Higgs + jet production at hadron colliders '
    write(idev,'(a)') 
    write(idev,'(a)') '         Written by Alexander Lind and Andrea Banfi'
    write(idev,'(a)') '                  arXiv:2011.04694 [hep-ph]'
    write(idev,'(a)') '============================================================'

    write(idev,*) ! Blank space for cleaner output 

  end subroutine print_welcome_banner

!=======================================================================================

  subroutine print_version_number 

    write(stderr,'(a)') 'H1jet version '//trim(version_number) 

    stop 

  end subroutine print_version_number 

!=======================================================================================

end module banner


