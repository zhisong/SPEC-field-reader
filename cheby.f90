! the Chebyshev polynomials
module cheby

contains

  subroutine get_cheby(n, lrad, s, T, dT)
    implicit none
    integer, intent(in) :: n, lrad
    real, dimension(n) :: s
    real, dimension(n, 0:lrad), intent(out) :: T
    real, dimension(n, 0:lrad), intent(out), optional :: dT

    integer :: ll

    T(:,0) = 1.0
    T(:,1) = s(:)
    do ll = 2, lrad
      T(:,ll) = 2.0 * s(:) * T(:,ll-1) - T(:,ll-2)
    end do

    if (PRESENT(dT)) then
      dT(:,0) = 0.0
      dT(:,1) = 1.0
      do ll = 2, lrad
        dT(:,ll) = 2.0 * T(:,ll-1) + 2.0 * s(:) * dT(:,ll-1) - dT(:, ll-2)
      end do
    end if
    
  end subroutine get_cheby


end module cheby