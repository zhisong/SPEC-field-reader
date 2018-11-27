! this module defines the data structure of A and the subroutine that obtains the field
module spec_field

  type, public :: Afield

    real, dimension(:,:), allocatable :: Ate, Ato, Aze, Azo
    real, dimension(:), allocatable :: im, in
    integer :: lrad, mn, Nfp
    logical :: isym, isingular

  end type Afield

  !type(Afield), dimension(:), allocatable :: A

contains
  subroutine get_spec_field(af, lrad, s, theta, xi, a, gb)
  ! Get the magnetic vector potential and/or magnetic field
  ! INPUTS:
  ! af     - TYPE(Afield), the object contains the field information
  ! lrad   - interger, the order of Chebyshev polynomials
  ! s      - REAL, the s coordinate
  ! theta  - REAL, the theta coordinate
  ! xi     - REAL, the xi coordinate
  !
  ! RETURNS:
  ! a      - REAL(3), the vector potential
  ! gb     - REAL(3), the magnetic field times the Jacobian g

    use cheby
    implicit none
    type(Afield), intent(in) :: af
    real, intent(in) :: s, theta, xi
    real, dimension(3), intent(out) :: a
    real, dimension(3), intent(out) :: gb

    real, dimension(0:lrad) :: T, dT, Tbar, dTbar
    real, dimension(0:lrad) :: work
    integer :: ii, lrad, mn
    real :: alphai, sinai, cosai, sbar, sbarmi

    mn = af%mn     ! shorthand

    gb(:) = 0
    a(:) = 0

    sbar = (s + 1.0) / 2.0
    call get_cheby(lrad, s, T, dT)

    do ii = 1, mn
      ! alphai = m * theta - n * xi
      alphai = af%im(ii) * theta - af%in(ii) * xi
      cosai = COS(alphai)
      sinai = SIN(alphai)

      if (af%isingular) then
        sbarmi = sbar**(af%im(ii) / 2.0)
        Tbar(:) = T(:) * sbarmi
        dTbar(:) = dT(:) * sbarmi + T(:) * sbarmi / sbar * af%im(ii) / 4.0
      else
        Tbar(:) = T(:)
        dTbar(:) = dT(:)
      end if ! if Lcoordinatesingularity == .true.

      ! now we first calculate vector potential
      a(2) = a(2) + SUM(Tbar(:) * af%Ate(:,ii)) * cosai
      a(3) = a(3) + SUM(Tbar(:) * af%Aze(:,ii)) * cosai

      if (.not. af%isym) then
        a(2) = a(2) + SUM(Tbar(:) * af%Ato(:,ii)) * sinai
        a(3) = a(3) + SUM(Tbar(:) * af%Azo(:,ii)) * sinai
      end if
   
      ! then we need to calculate the field
      gb(1) = gb(1) - (af%im(ii) * SUM(Tbar(:) * af%Aze(:,ii)) &
                    +  af%in(ii) * SUM(Tbar(:) * af%Ate(:,ii))) * sinai
      gb(2) = gb(2) - SUM(dTbar(:) * af%Aze(:,ii)) * cosai
      gb(3) = gb(3) + SUM(dTbar(:) * af%Ate(:,ii)) * cosai

      if (.not. af%isym) then
        gb(1) = gb(1) + (af%im(ii) * SUM(Tbar(:) * af%Azo(:,ii)) &
                      +  af%in(ii) * SUM(Tbar(:) * af%Ato(:,ii))) * cosai
        gb(2) = gb(2) - SUM(dTbar(:) * af%Azo(:,ii)) * sinai
        gb(3) = gb(3) + SUM(dTbar(:) * af%Ato(:,ii)) * sinai
      end if ! if Lcoordinatesingularity == .true.
    end do

  end subroutine get_spec_field

  subroutine destroy_field(f)
    implicit none
    type(Afield) :: f
    ! free the space for field f
    if (ALLOCATED(f%Ate)) deallocate(f%Ate)
    if (ALLOCATED(f%Ato)) deallocate(f%Ato)
    if (ALLOCATED(f%Aze)) deallocate(f%Aze)
    if (ALLOCATED(f%Azo)) deallocate(f%Azo)
    if (ALLOCATED(f%im)) deallocate(f%im)
    if (ALLOCATED(f%im)) deallocate(f%in)
  end subroutine destroy_field

end module spec_field