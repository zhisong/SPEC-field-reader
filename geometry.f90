! this module defines the data structure of geometry and the calculation of gij/dgij
module geometry

  type, public :: volume

    real, dimension(:,:), allocatable :: Rbc, Rbs, Zbc, Zbs
    logical :: icoordinatesingularity
    logical :: isym
    integer, dimension(:), allocatable :: im, in
    integer :: mn, nfp, mvol, igeometry

  end type volume

contains

  subroutine get_coordinate(v, lvol, mn, s, theta, xi, x, gij, dgij)
  ! Obtain the coordinate quantities
  ! INPUTS:
  ! v     - TYPE(volume), the volume object read from file
  ! lvol  - INTEGER, which volume we are looking at,
  ! mn    - number of harmonics
  ! s     - s coordinate
  ! theta - theta coordinate
  ! xi    - xi coordinate
  ! RETURNS:
  ! x     - REAL(3), the coordinates (x, y, z)
  ! gij   - REAL(3,3), metric tensor with lower indices
  ! dgij  - REAL(3,3,3), the derivative of gij with respect to (s, theta, xi)
    implicit none
    type(volume) :: v
    integer, intent(in) :: lvol, mn
    real, intent(in) :: s, theta, xi

    real, dimension(3), intent(out) :: x
    real, dimension(3,3), intent(out) :: gij
    real, dimension(3,3,3), intent(out) :: dgij

    real :: Rij(0:3,0:3), Zij(0:3,0:3)
    real :: sbar, alss, blss
    real, dimension(mn) :: alphai, cosai, sinai
    real, dimension(mn) :: t1, t2, t3, t4, ddt1, ddt3, fj, dfj, ddfj
    integer :: ii, jj, kk

    Rij(:,:) = 0
    Zij(:,:) = 0
    x(:) = 0
    gij(:,:) = 0
    dgij(:,:,:) = 0

    if (v%icoordinatesingularity .and. lvol == 1) then
      ! compute R
      sbar = (1.0 + s) / 2.0
      fj(1) = sbar
      dfj(1) = 0.5
      ddfj(1) = 0.0

      ddfj(2:mn) = sbar**(v%im(2:mn)/2.0 - 2.0)
      dfj(2:mn) = ddfj(2:mn) * sbar * v%im(2:mn) / 4.0
      fj(2:mn) = ddfj(2:mn) * sbar**2
      ddfj(2:mn) = ddfj(2:mn) * v%im(2:mn) * (v%im(2:mn) - 2.0) / 16.0

      t1(:) = v%Rbc(ii,0) + (v%Rbc(ii,1) - v%Rbc(ii,0)) * fj(:)
      t2(:) = (v%Rbc(ii,1) - v%Rbc(ii,0)) * dfj(:)
      ddt1(:) = (v%Rbc(ii,1) - v%Rbc(ii,0)) * ddfj(:)
      if (.not. v%isym) then
        t3(:) = v%Rbs(ii,0) + (v%Rbs(ii,1) - v%Rbs(ii,0)) * fj(:)
        t4(:) = (v%Rbs(ii,1) - v%Rbs(ii,0)) * dfj(:)
        ddt3(:) = (v%Rbs(ii,1) - v%Rbs(ii,0)) * ddfj(:)
      end if

    else
      alss = 0.5 * ( 1.0 - s )
      blss = 0.5 * ( 1.0 + s )
      t1(:) = (alss * v%Rbc(ii,lvol-1) + blss * v%Rbc(ii,lvol))
      t2(:) = (-0.5 * v%Rbc(ii,lvol-1) + 0.5 * v%Rbc(ii,lvol))
      ddt1(:) = 0.0
      if (.not. v%isym) then
        t3(:) = (alss * v%Rbs(ii,lvol-1) + blss * v%Rbs(ii,lvol))
        t4(:) = (-0.5 * v%Rbs(ii,lvol-1) + 0.5 * v%Rbs(ii,lvol))
        ddt3(:) = 0.0
      end if
    end if


    alphai = v%im * theta - v%in * xi
    cosai = COS(alphai)
    sinai = SIN(alphai)

    Rij(0,0) = SUM(t1 * cosai)
    Rij(0,1) = SUM(t2 * cosai)
    Rij(0,2) = SUM(t1 * (-v%im * sinai))
    Rij(0,3) = SUM(t1 * (v%in * sinai))
    Rij(1,1) = SUM(ddt1 * cosai)
    Rij(1,2) = SUM(t2 * (-v%im * sinai))
    Rij(1,3) = SUM(t2 * (v%in * sinai))
    Rij(2,2) = SUM(t1 * (-v%im**2 * cosai))
    Rij(2,3) = SUM(t1 * (v%im * v%in * cosai))
    Rij(3,3) = SUM(t1 * (-v%in**2 * cosai))

    if (.not. v%isym) then
      Rij(0,0) = Rij(0,0) + SUM(t3 * sinai)
      Rij(0,1) = Rij(0,1) + SUM(t4 * sinai)
      Rij(0,2) = Rij(0,2) + SUM(t3 * (v%im * cosai))
      Rij(0,3) = Rij(0,3) + SUM(t3 * (-v%in(ii) * cosai))
      Rij(1,1) = Rij(1,1) + SUM(ddt3 * sinai)
      Rij(1,2) = Rij(1,2) + SUM(t4 * (v%im(ii) * cosai))
      Rij(1,3) = Rij(1,3) + SUM(t4 * (-v%in(ii) * cosai))
      Rij(2,2) = Rij(2,2) + SUM(t3 * (-v%im(ii)**2 * sinai))
      Rij(2,3) = Rij(2,3) + SUM(t3 * (v%im(ii) * v%in(ii) * sinai))
      Rij(3,3) = Rij(3,3) + SUM(t3 * (-v%in(ii)**2 * sinai))
    end if

    Rij(2,1) = Rij(1,2)
    Rij(3,1) = Rij(1,3)
    Rij(3,2) = Rij(2,3)

    if (v%igeometry == 1) then ! Cartesian
      x(1) = theta
      x(2) = xi
      x(3) = Rij(0,0)

      do ii = 1, 3
        do jj = 1, 3
          gij(jj,ii) = gij(jj,ii) + Rij(0,jj) * Rij(0,ii)
          do kk = 1, 3
            dgij(jj,ii,kk) = dgij(jj,ii,kk) + Rij(kk,jj) * Rij(0,ii) + Rij(0,jj) * Rij(kk,ii)
          end do ! kk
        end do ! jj
      end do ! ii

      gij(2,2) = gij(2,2) + 1.0
      gij(3,3) = gij(3,3) + 1.0

    elseif (v%igeometry == 2) then ! Cylindrical
    elseif (v%igeometry == 3) then ! Toroidal
    else                           ! Unknown, error
      STOP 'Igeometry must be 1 to 3'
    end if ! igeometry

  end subroutine get_coordinate

  subroutine destroy_volume(v)
    implicit none
    type(volume) :: v
    if (ALLOCATED(v%im)) deallocate(v%im)
    if (ALLOCATED(v%im)) deallocate(v%im)
  end subroutine destroy_volume

end module geometry