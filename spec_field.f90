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
  subroutine get_spec_field(af, s, theta, xi, a, gb, dgb)
  ! Get the magnetic vector potential and/or magnetic field
  ! INPUTS:
  ! af     - TYPE(Afield), the object contains the field information
  ! s      - REAL, the s coordinate
  ! theta  - REAL, the theta coordinate
  ! xi     - REAL, the xi coordinate
  !
  ! RETURNS:
  ! a      - REAL(3), the vector potential
  ! gb     - REAL(3), the magnetic field times the Jacobian g
  ! dgb    - REAL(3,3), the derivative of gb with respect to (s, theta, xi)
  !                     first index: component; second index: derivative

    use cheby
    implicit none
    type(Afield), intent(in) :: af
    real, intent(in) :: s, theta, xi
    real, dimension(3), intent(out) :: a
    real, dimension(3), intent(out) :: gb
    real, dimension(3,3), intent(out) :: dgb

    real, dimension(0:af%lrad) :: T, dT, ddT, Tbar, dTbar, ddTbar
    real, dimension(0:af%lrad) :: work

    integer :: ii, lrad, mn
    real :: alphai, sinai, cosai, sbar, sbarmi
    real :: sumTbarAte, sumTbarAze, sumTbarAto, sumTbarAzo
    real :: sumdTbarAte, sumdTbarAze, sumdTbarAto, sumdTbarAzo    

    lrad = af%lrad
    mn = af%mn     ! shorthand

    gb(:) = 0
    a(:) = 0

    sbar = (s + 1.0) / 2.0
    call get_cheby(lrad, s, T, dT, ddT)

    do ii = 1, mn
      ! alphai = m * theta - n * xi
      alphai = af%im(ii) * theta - af%in(ii) * xi
      cosai = COS(alphai)
      sinai = SIN(alphai)

      if (af%isingular) then
        sbarmi = sbar**(af%im(ii) / 2.0)
        Tbar(:) = T(:) * sbarmi
        dTbar(:) = dT(:) * sbarmi + T(:) * sbarmi / sbar * af%im(ii) / 4.0
        ddTbar(:) = ddT(:) * sbarmi + dT(:) * sbarmi / sbar * af%im(ii) / 2.0 &
                  + T(:) * sbarmi/sbar**2 * af%im(ii) * (af%im(ii)/2.0 - 1) / 8.0
      else
        Tbar(:) = T(:)
        dTbar(:) = dT(:)
        ddTbar(:) = ddT(:)
      end if ! if Lcoordinatesingularity == .true.

      ! to speed up calculations, first we calculate some quantities that can be reused
      sumTbarAte = SUM(Tbar(:) * af%Ate(:,ii))
      sumTbarAze = SUM(Tbar(:) * af%Aze(:,ii))
      sumdTbarAte = SUM(dTbar(:) * af%Ate(:,ii))
      sumdTbarAze = SUM(dTbar(:) * af%Aze(:,ii))

      if (.not. af%isym) then
        sumTbarAto = SUM(Tbar(:) * af%Ato(:,ii))
        sumTbarAzo = SUM(Tbar(:) * af%Azo(:,ii))
        sumdTbarAto = SUM(dTbar(:) * af%Ato(:,ii))
        sumdTbarAzo = SUM(dTbar(:) * af%Azo(:,ii))
      end if ! if not stellarator symmetric

      ! now we first calculate vector potential
      a(2) = a(2) + sumTbarAte * cosai
      a(3) = a(3) + sumTbarAze * cosai

      if (.not. af%isym) then
        a(2) = a(2) + sumTbarAto * sinai
        a(3) = a(3) + sumTbarAzo * sinai
      end if ! if not stellarator symmetric
   
      ! then we need to calculate the field
      gb(1) = gb(1) - (af%im(ii) * sumTbarAze &
                    +  af%in(ii) * sumTbarAte) * sinai
      gb(2) = gb(2) - sumdTbarAze * cosai
      gb(3) = gb(3) + sumdTbarAte * cosai

      if (.not. af%isym) then
        gb(1) = gb(1) + (af%im(ii) * sumTbarAzo &
                      +  af%in(ii) * sumTbarAto) * cosai
        gb(2) = gb(2) - sumdTbarAzo * sinai
        gb(3) = gb(3) + sumdTbarAto * sinai
      end if ! if not stellarator symmetric
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