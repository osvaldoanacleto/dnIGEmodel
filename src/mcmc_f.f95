module dnIGEmcmc
    use, intrinsic :: iso_c_binding
    implicit none
    private
    public :: runMCMC_f

contains


  subroutine runMCMC_f(xvector, xvector_shifted) bind(C, name = "runMCMC_f_")
      !real(kind = c_double), intent(in)               :: mean, sd
      real(kind = c_double), intent(in), dimension(100) :: xvector
      real(kind = c_double), intent(out), dimension(100) :: xvector_shifted

      !llc = 0.0_c_double

      xvector_shifted = xvector + 70.0_c_double

  end subroutine runMCMC_f

end module dnIGEmcmc
