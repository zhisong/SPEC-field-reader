! this module defines the data structure of geometry and the calculation of gij/dgij
module spec_geometry

  type, public :: volume

    real, dimension(:,:), allocatable :: Rbc, Rbs, Zbc, Zbs
    logical :: icoordinatesingularity
    logical :: isym
    integer, dimension(:), allocatable :: im, in
    integer :: mn, nfp, mvol, igeometry, mpol, ntor

  end type volume

contains

  subroutine get_spec_coord(v, lvol, mn, s, theta, xi, jac, x, gij, dgij)
  ! Obtain the coordinate quantities
  ! INPUTS:
  ! v     - TYPE(volume), the volume object read from file
  ! lvol  - INTEGER, which volume we are looking at,
  ! mn    - number of harmonics
  ! s     - s coordinate
  ! theta - theta coordinate
  ! xi    - xi coordinate
  ! RETURNS:
  ! jac   - REAL, Jacobian
  ! x     - REAL(3), the coordinates (x, y, z)
  ! gij   - REAL(3,3), metric tensor with lower indices
  ! dgij  - REAL(3,3,3), the derivative of gij with respect to (s, theta, xi)
    implicit none
    type(volume) :: v
    integer, intent(in) :: lvol, mn
    real, intent(in) :: s, theta, xi

    real, intent(out) :: jac
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

    alphai = v%im * theta - v%in * xi
    cosai = COS(alphai)
    sinai = SIN(alphai)

    if (v%icoordinatesingularity .and. lvol == 1) then
      sbar = (1.0 + s) / 2.0
      fj(1:v%mpol+1) = sbar
      dfj(1:v%mpol+1) = 0.5
      ddfj(1:v%mpol+1) = 0.0

      if (v%ntor > 0) then
        ddfj(v%mpol+2:mn) = sbar**(v%im(v%mpol+2:mn)/2.0 - 2.0)
        dfj(v%mpol+2:mn) = ddfj(v%mpol+2:mn) * sbar * v%im(v%mpol+2:mn) / 4.0
        fj(v%mpol+2:mn) = ddfj(v%mpol+2:mn) * sbar**2
        ddfj(v%mpol+2:mn) = ddfj(v%mpol+2:mn) * v%im(v%mpol+2:mn) * (v%im(v%mpol+2:mn) - 2.0) / 16.0
      end if

      t1(:) = v%Rbc(:,0) + (v%Rbc(:,1) - v%Rbc(:,0)) * fj(:)
      t2(:) = (v%Rbc(:,1) - v%Rbc(:,0)) * dfj(:)
      ddt1(:) = (v%Rbc(:,1) - v%Rbc(:,0)) * ddfj(:)
      if (.not. v%isym) then
        t3(:) = v%Rbs(:,0) + (v%Rbs(:,1) - v%Rbs(:,0)) * fj(:)
        t4(:) = (v%Rbs(:,1) - v%Rbs(:,0)) * dfj(:)
        ddt3(:) = (v%Rbs(:,1) - v%Rbs(:,0)) * ddfj(:)
      end if

    else
      alss = 0.5 * ( 1.0 - s )
      blss = 0.5 * ( 1.0 + s )
      t1(:) = (alss * v%Rbc(:,lvol-1) + blss * v%Rbc(:,lvol))
      t2(:) = (-0.5 * v%Rbc(:,lvol-1) + 0.5 * v%Rbc(:,lvol))
      ddt1(:) = 0.0
      if (.not. v%isym) then
        t3(:) = (alss * v%Rbs(:,lvol-1) + blss * v%Rbs(:,lvol))
        t4(:) = (-0.5 * v%Rbs(:,lvol-1) + 0.5 * v%Rbs(:,lvol))
        ddt3(:) = 0.0
      end if
    end if

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
      Rij(0,3) = Rij(0,3) + SUM(t3 * (-v%in * cosai))
      Rij(1,1) = Rij(1,1) + SUM(ddt3 * sinai)
      Rij(1,2) = Rij(1,2) + SUM(t4 * (v%im * cosai))
      Rij(1,3) = Rij(1,3) + SUM(t4 * (-v%in * cosai))
      Rij(2,2) = Rij(2,2) + SUM(t3 * (-v%im**2 * sinai))
      Rij(2,3) = Rij(2,3) + SUM(t3 * (v%im * v%in * sinai))
      Rij(3,3) = Rij(3,3) + SUM(t3 * (-v%in**2 * sinai))
    end if

    Rij(2,1) = Rij(1,2)
    Rij(3,1) = Rij(1,3)
    Rij(3,2) = Rij(2,3)

    if (v%igeometry == 1) then ! Cartesian
      jac = Rij(0,1)

      x(1) = theta
      x(2) = xi
      x(3) = Rij(0,0)

      do ii = 1, 3
        do jj = ii, 3
          gij(jj,ii) = gij(jj,ii) + Rij(0,jj) * Rij(0,ii)
          do kk = 1, 3
            dgij(jj,ii,kk) = dgij(jj,ii,kk) + Rij(kk,jj) * Rij(0,ii) + Rij(0,jj) * Rij(kk,ii)
          end do ! kk
        end do ! jj
      end do ! ii

      gij(2,2) = gij(2,2) + 1.0
      gij(3,3) = gij(3,3) + 1.0

    elseif (v%igeometry == 2) then ! Cylindrical

      jac = Rij(0,0) * Rij(0,1)
      x(1) = Rij(0,0) * COS(theta)
      x(2) = Rij(0,0) * SIN(theta)
      x(3) = xi

      do ii = 1, 3
        do jj = ii, 3
          gij(jj,ii) = gij(jj,ii) + Rij(0,jj) * Rij(0,ii)
          do kk = 1, 3
            dgij(jj,ii,kk) = dgij(jj,ii,kk) + Rij(kk,jj) * Rij(0,ii) + Rij(0,jj) * Rij(kk,ii)
          end do ! kk
        end do ! jj
      end do ! ii

      gij(2,2) = gij(2,2) + Rij(0,0)**2
      gij(3,3) = gij(3,3) + 1.0
      do ii = 1, 3
        dgij(2,2,ii) = dgij(2,2,ii) + 2.0 * Rij(0,0) * Rij(0,ii)
      end do ! ii

    elseif (v%igeometry == 3) then ! Toroidal
    ! compute Zij
      if (v%icoordinatesingularity .and. lvol == 1) then
  
        t1(:) = v%Zbs(:,0) + (v%Zbs(:,1) - v%Zbs(:,0)) * fj(:)
        t2(:) = (v%Zbs(:,1) - v%Zbs(:,0)) * dfj(:)
        ddt1(:) = (v%Zbs(:,1) - v%Zbs(:,0)) * ddfj(:)
        if (.not. v%isym) then
          t3(:) = v%Zbc(:,0) + (v%Zbc(:,1) - v%Zbc(:,0)) * fj(:)
          t4(:) = (v%Zbc(:,1) - v%Zbc(:,0)) * dfj(:)
          ddt3(:) = (v%Zbc(:,1) - v%Zbc(:,0)) * ddfj(:)
        end if

      else
        alss = 0.5 * ( 1.0 - s )
        blss = 0.5 * ( 1.0 + s )
        t1(:) = (alss * v%Zbs(:,lvol-1) + blss * v%Zbs(:,lvol))
        t2(:) = (-0.5 * v%Zbs(:,lvol-1) + 0.5 * v%Zbs(:,lvol))
        ddt1(:) = 0.0
        if (.not. v%isym) then
          t3(:) = (alss * v%Zbc(:,lvol-1) + blss * v%Zbc(:,lvol))
          t4(:) = (-0.5 * v%Zbc(:,lvol-1) + 0.5 * v%Zbc(:,lvol))
          ddt3(:) = 0.0
        end if
      end if

      Zij(0,0) = SUM(t1 * sinai)
      Zij(0,1) = SUM(t2 * sinai)
      Zij(0,2) = SUM(t1 * (v%im * cosai))
      Zij(0,3) = SUM(t1 * (-v%in * cosai))
      Zij(1,1) = SUM(ddt1 * sinai)
      Zij(1,2) = SUM(t2 * (v%im * cosai))
      Zij(1,3) = SUM(t2 * (-v%in * cosai))
      Zij(2,2) = SUM(t1 * (-v%im**2 * sinai))
      Zij(2,3) = SUM(t1 * (v%im * v%in * sinai))
      Zij(3,3) = SUM(t1 * (-v%in**2 * sinai))

      if (.not. v%isym) then
        Zij(0,0) = Zij(0,0) + SUM(t3 * cosai)
        Zij(0,1) = Zij(0,1) + SUM(t4 * cosai)
        Zij(0,2) = Zij(0,2) + SUM(t3 * (-v%im * sinai))
        Zij(0,3) = Zij(0,3) + SUM(t3 * (v%in * sinai))
        Zij(1,1) = Zij(1,1) + SUM(ddt3 * cosai)
        Zij(1,2) = Zij(1,2) + SUM(t4 * (-v%im * sinai))
        Zij(1,3) = Zij(1,3) + SUM(t4 * (v%in * sinai))
        Zij(2,2) = Zij(2,2) + SUM(t3 * (-v%im**2 * cosai))
        Zij(2,3) = Zij(2,3) + SUM(t3 * (v%im * v%in * cosai))
        Zij(3,3) = Zij(3,3) + SUM(t3 * (-v%in**2 * cosai))
      end if

      Zij(2,1) = Zij(1,2)
      Zij(3,1) = Zij(1,3)
      Zij(3,2) = Zij(2,3)

      jac = Rij(0,0) * (Zij(0,1)*Rij(0,2) - Rij(0,1)*Zij(0,2))
      x(1) = Rij(0,0) * COS(xi)
      x(2) = Rij(0,0) * SIN(xi)
      x(3) = Zij(0,0)

      do ii = 1, 3
        do jj = ii, 3
          gij(jj,ii) = gij(jj,ii) + Rij(0,jj) * Rij(0,ii) + Zij(0,jj) * Zij(0,ii)
          do kk = 1, 3
            dgij(jj,ii,kk) = dgij(jj,ii,kk) + Rij(kk,jj) * Rij(0,ii) + Rij(0,jj) * Rij(kk,ii) &
                                            + Zij(kk,jj) * Zij(0,ii) + Zij(0,jj) * Zij(kk,ii)
          end do ! kk
        end do ! jj
      end do ! ii

      gij(3,3) = gij(3,3) + Rij(0,0)**2
      do ii = 1, 3
        dgij(3,3,ii) = dgij(3,3,ii) + 2.0 * Rij(0,0) * Rij(0,ii)
      end do ! ii
    else                           ! Unknown, error
      STOP 'Igeometry must be 1 to 3'
    end if ! igeometry

    gij(1,2) = gij(2,1)
    gij(1,3) = gij(3,1)
    gij(2,3) = gij(3,2)
    do ii = 1, 3
      dgij(1,2,ii) = dgij(2,1,ii)
      dgij(1,3,ii) = dgij(3,1,ii)
      dgij(2,3,ii) = dgij(3,2,ii)
    end do 

  end subroutine get_spec_coord

  subroutine destroy_volume(v)
    implicit none
    type(volume) :: v
    if (ALLOCATED(v%im)) deallocate(v%im)
    if (ALLOCATED(v%im)) deallocate(v%im)
  end subroutine destroy_volume

end module spec_geometry