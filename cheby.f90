! the Chebyshev polynomials
module cheby

  interface get_cheby
    ! We have implemented two version of get_cheby, with single input and multiple inputs 
    ! We have only used the subroutine for single input, the multiple inputs version
    ! may be used in future version of the code
    module procedure get_cheby_single !, get_cheby_multi 
  end interface get_cheby

contains

  subroutine get_cheby_single(lrad, s, T, dT, ddT)
  ! Get the value and the derivatives of Chebyshev polynomials at position s
  ! with one input at one time
  ! INPUTS:
  ! lrad - INTEGER, the order of polynomials
  ! s    - REAL, the input list of position s
  !
  ! RETURNS:
  ! T    - REAL(0:lrad), the value of the polynomials
  ! dT   - REAL(0:lrad), optional, the value of the first derivative over s
  ! ddT  - REAL(0:lrad), optional, the value of the second derivative over s
  
    implicit none
    integer, intent(in) :: lrad
    real, intent(in) :: s
    real, dimension(0:lrad), intent(out) :: T
    real, dimension(0:lrad), optional, intent(out) :: dT, ddT

    integer :: ll

    T(0) = 1.0
    T(1) = s
    do ll = 2, lrad
      T(ll) = 2.0 * s * T(ll-1) - T(ll-2)
    end do
    
    if (PRESENT(dT)) then
      dT(0) = 0.0
      dT(1) = 1.0
      do ll = 2, lrad
        dT(ll) = 2.0 * T(ll-1) + 2.0 * s * dT(ll-1) - dT(ll-2)
      end do
    end if

    if (PRESENT(ddT)) then
      ddT(0) = 0
      ddT(1) = 0
      do ll = 2, lrad
        ddT(ll) = 4.0 * dT(ll-1) + 2.0 * s * ddT(ll-1) - ddT(ll-2)
      end do
    end if

  end subroutine get_cheby_single

  ! subroutine get_cheby_multi(n, lrad, s, T, dT)
  ! ! Get the value and the derivatives of Chebyshev polynomials at position s
  ! ! with n inputs at one time
  ! ! INPUTS:
  ! ! n    - INTEGER, the length of the list of input
  ! ! lrad - INTEGER, the order of polynomials
  ! ! s    - REAL(n), the input list of position s
  ! !
  ! ! RETURNS:
  ! ! T    - REAL(n, 0:lrad), the value of the polynomials
  ! ! dT   - REAL(n, 0:lrad), the value of the first derivative over s
  ! ! ddT  - REAL(n, 0:lrad), the value of the second derivative over s

  !   implicit none
  !   integer, intent(in) :: n, lrad
  !   real, dimension(n) :: s
  !   real, dimension(n, 0:lrad), intent(out) :: T
  !   real, dimension(n, 0:lrad), intent(out) :: dT

  !   integer :: ll

  !   T(:,0) = 1.0
  !   T(:,1) = s(:)
  !   do ll = 2, lrad
  !     T(:,ll) = 2.0 * s(:) * T(:,ll-1) - T(:,ll-2)
  !   end do

  !   dT(:,0) = 0.0
  !   dT(:,1) = 1.0
  !   do ll = 2, lrad
  !     dT(:,ll) = 2.0 * T(:,ll-1) + 2.0 * s(:) * dT(:,ll-1) - dT(:, ll-2)
  !   end do
    
  ! end subroutine get_cheby_multi

end module cheby