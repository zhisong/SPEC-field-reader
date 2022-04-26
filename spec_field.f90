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
  subroutine get_spec_field(af, s, theta, xi, mpol, a, gb, dgb)
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

    use poly
    implicit none
    type(Afield), intent(in) :: af
    REAL, INTENT(IN) :: s, theta, xi
    INTEGER, INTENT(IN) :: mpol

    real, dimension(3), intent(out) :: a
    real, dimension(3), intent(out) :: gb
    real, dimension(3,3), intent(out) :: dgb

    real, dimension(0:af%lrad) :: T, dT, ddT, Tbar, dTbar, ddTbar
    real, dimension(0:af%lrad) :: work

    real, dimension(0:2) :: TT

    REAL, DIMENSION(0:af%lrad,0:mpol,0:2) :: zernike
    REAL, DIMENSION(0:af%lrad,0:2) :: cheby

    integer :: ii, ll, lrad, mn
    real :: im, in
    real :: alphai, sinai, cosai, sbar, sbarmi
    real :: sumTbarAte, sumTbarAze, sumTbarAto, sumTbarAzo
    real :: sumdTbarAte, sumdTbarAze, sumdTbarAto, sumdTbarAzo
    real :: sumddTbarAte, sumddTbarAze, sumddTbarAto, sumddTbarAzo  

    lrad = af%lrad
    mn = af%mn     ! shorthand

    dgb(:,:) = 0
    gb(:) = 0
    a(:) = 0


    IF (af%isingular) THEN
      sbar = MAX((s + 1.0)/2.0, 0.0)
      CALL get_zernike(sbar,lrad,mpol,zernike)
    ELSE
      CALL get_cheby(s,lrad,cheby)
    END IF


    DO ii = 1, mn
      im = af%im(ii)
      in = af%in(ii)

      alphai = im*theta - in*xi
      cosai = COS(alphai)
      sinai = SIN(alphai)


      ! if (af%isingular) then
      !   sbar = (s + 1.0) / 2.0
      !   sbarmi = sbar**(im / 2.0)
      !   Tbar(:) = T(:) * sbarmi
      !   dTbar(:) = dT(:) * sbarmi + T(:) * sbarmi / sbar * im / 4.0
      !   ddTbar(:) = ddT(:) * sbarmi + dT(:) * sbarmi / sbar * im / 2.0 &
      !             + T(:) * sbarmi/sbar**2 * im * (im/2.0 - 1) / 8.0
      ! else
      !   Tbar(:) = T(:)
      !   dTbar(:) = dT(:)
      !   ddTbar(:) = ddT(:)
      ! end if ! if Lcoordinatesingularity == .true.


      DO ll = 0, lrad

        IF (af%isingular) THEN
          TT(0) = zernike(ll,Int(im),0)
          TT(1) = zernike(ll,Int(im),1)/2.0
          TT(2) = zernike(ll,Int(im),2)/4.0
        ELSE
          TT(0) = cheby(ll,0)
          TT(1) = cheby(ll,1)
          TT(2) = cheby(ll,2)
        ENDIF

        ! Vector potential
        a(2) = a(2) + af%Ate(ll,ii) * TT(0) * cosai
        a(3) = a(3) + af%Aze(ll,ii) * TT(0) * cosai
        ! Magnetic field
        gb(1) = gb(1) - (im * af%Aze(ll,ii) + in * af%Ate(ll,ii)) * TT(0) * sinai
        gb(2) = gb(2) - af%Aze(ll,ii) * TT(1) * cosai
        gb(3) = gb(3) + af%Ate(ll,ii) * TT(1) * cosai

        !--- DERIVATIVES
        ! Derivatives (ds)
        dgb(1,1) = dgb(1,1) - (im * af%Aze(ll,ii) + in * af%Ate(ll,ii)) * TT(1) * sinai
        dgb(2,1) = dgb(2,1) - af%Aze(ll,ii) * TT(2) * cosai
        dgb(3,1) = dgb(3,1) + af%Ate(ll,ii) * TT(2) * cosai
        ! Derivatives (dtheta)
        dgb(1,2) = dgb(1,2) - im * (im*af%Aze(ll,ii) + in*af%Ate(ll,ii)) * TT(0) * cosai
        dgb(2,2) = dgb(2,2) + im * af%Aze(ll,ii) * TT(1) * sinai
        dgb(3,2) = dgb(3,2) - im * af%Ate(ll,ii) * TT(1) * sinai
        ! Derivatives (dzeta)
        dgb(1,3) = dgb(1,3) + in * (im*af%Aze(ll,ii) + in*af%Ate(ll,ii)) * TT(0) * cosai
        ! dgb(1,3) = dgb(1,3) - in * af%Aze(ll,ii)  * TT(1) * sinai
        ! dgb(1,3) = dgb(1,3) + in * af%Ate(ll,ii)  * TT(1) * sinai
        dgb(2,3) = dgb(2,3) - in * af%Aze(ll,ii)  * TT(1) * sinai
        dgb(3,3) = dgb(3,3) + in * af%Ate(ll,ii)  * TT(1) * sinai

        IF (.NOT.af%isym) THEN

          a(2) = a(2) + af%Ato(ll,ii) * TT(0) * sinai
          a(3) = a(3) + af%Azo(ll,ii) * TT(0) * sinai

          gb(1) = gb(1) + (im * af%Azo(ll,ii) + in * af%Ato(ll,ii)) * TT(0) * cosai
          gb(2) = gb(2) - af%Azo(ll,ii) * TT(1) * sinai
          gb(3) = gb(3) + af%Ato(ll,ii) * TT(1) * sinai

          dgb(1,1) = dgb(1,1) + (im * af%Azo(ll,ii) + in * af%Ato(ll,ii)) * TT(1) * cosai
          dgb(2,1) = dgb(2,1) - af%Azo(ll,ii) * TT(2) * sinai
          dgb(3,1) = dgb(3,1) + af%Ato(ll,ii) * TT(2) * sinai

          dgb(1,2) = dgb(1,2) - im * (im*af%Azo(ll,ii) + in*af%Ato(ll,ii)) * TT(0) * sinai
          dgb(2,2) = dgb(2,2) - im * af%Azo(ll,ii)  * TT(1) * cosai
          dgb(3,2) = dgb(3,2) + im * af%Ato(ll,ii)  * TT(1) * cosai

          dgb(1,3) = dgb(1,3) + in * (im*af%Azo(ll,ii) + in*af%Ato(ll,ii)) * TT(0) * sinai
          ! dgb(1,3) = dgb(1,3) + in * af%Azo(ll,ii)  * TT(1) * cosai
          ! dgb(1,3) = dgb(1,3) - in * af%Ato(ll,ii)  * TT(1) * cosai
          dgb(2,3) = dgb(2,3) + in * af%Azo(ll,ii)  * TT(1) * cosai
          dgb(3,3) = dgb(3,3) - in * af%Ato(ll,ii)  * TT(1) * cosai
        ENDIF

      END DO

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