MODULE poly

  CONTAINS

  SUBROUTINE get_cheby(s, lrad, cheby)
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
  
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: lrad
    REAL, INTENT(IN) :: s

    REAL, INTENT(INOUT), DIMENSION(0:lrad,0:2) :: cheby

    INTEGER :: ll

    cheby(0,0:2) = (/ 1.0, 0.0, 0.0 /)
    cheby(1,0:2) = (/ s, 1.0, 0.0 /)

    DO ll = 2, lrad
      cheby(ll,0) = 2.0 * s * cheby(ll-1,0) - cheby(ll-2,0) - (-1.0)**ll  ! chebychev
      cheby(ll,1) = 2.0 * cheby(ll-1,0) + 2.0 * s * cheby(ll-1,1) - cheby(ll-2,1) ! derivative
      cheby(ll,2) = 4.0 * cheby(ll-1,1) + 2.0 * s * cheby(ll-1,2) - cheby(ll-2,2) ! second derivative
    END DO

    DO ll = 1, lrad
      cheby(ll, 0) = cheby(ll,0) - (-1.0)**ll
    END DO

    DO ll = 0, lrad
      cheby(ll, 0:2) = cheby(ll, 0:2) / FLOAT(ll+1)
    END DO
      ! cheby(ll,0:2) = cheby(ll,0:2) / FLOAT(ll+1)
    ! cheby(1,0) = ( cheby(1,0) + 1.0 ) / 2.0
    ! cheby(0,0:2) = cheby(0,0:2) / 1.0

  END SUBROUTINE get_cheby




!----------  ZERNIKE POLYNOMIALS ----------!
  !> \brief Get the Zernike polynomials  \f$\hat{R}^{m}_{l}\f$ with zeroth, first, second derivatives
  !>
  !> See get_zernike for more detail.
  !>
  !> @param[in] r coordinate input, note that this is normalized to \f$[0, 1]\f$
  !> @param[in] lrad radial resolution
  !> @param[in] mpol poloidal resolution
  !> @param[out] zernike the value, first/second derivative of Zernike polynomial
  SUBROUTINE get_zernike(r, lrad, mpol, zernike)

    IMPLICIT NONE

    REAL,INTENT(IN) :: r
    INTEGER, intent(in) :: lrad, mpol
    REAL, intent(inout) :: zernike(0:lrad,0:mpol,0:2)

    REAL ::    rm, rm1, rm2  ! r to the power of m'th, m-1'th and m-2'th
    REAL ::    factor1, factor2, factor3, factor4
    INTEGER :: m, n  ! Zernike R^m_n
    
    rm = 1.0  ! r to the power of m'th
    rm1 = 0. ! r to the power of m-1'th
    rm2 = 0. ! r to the power of m-2'th
    zernike(:,:,:) = 0.
    do m = 0, mpol
      if (lrad >= m) then
        zernike(m,m,0:2) = (/ rm, real(m)*rm1, real(m*(m-1))*rm2 /)
        !write(0, *) m, m, r, zernike(m,m,:)
      endif

      if (lrad >= m+2) then
        zernike(m+2,m,0) = real(m+2)*rm*r**2 - real(m+1)*rm
        zernike(m+2,m,1) = real((m+2)**2)*rm*r - real((m+1)*m)*rm1
        zernike(m+2,m,2) = real((m+2)**2*(m+1))*rm - real((m+1)*m*(m-1))*rm2
        !write(0, *) m+2, m, r, zernike(m+2,m,:)
      endif

      do n = m+4, lrad, 2
        factor1 = real(n)/real(n**2 - m**2)
        factor2 = real(4 * (n-1))
        factor3 = real((n-2+m)**2)/real(n-2) + real((n-m)**2)/real(n)
        factor4 = real((n-2)**2-m**2) / real(n-2)

        zernike(n, m, 0) = factor1 * ((factor2*r**2 - factor3)*zernike(n-2,m,0) - factor4*zernike(n-4,m,0))
        zernike(n, m, 1) = factor1 * (2.0*factor2*r*zernike(n-2,m,0) + &
          (factor2*r**2 - factor3)*zernike(n-2,m,1) - factor4*zernike(n-4,m,1))
        zernike(n, m, 2) = factor1 * (2.0*factor2*(2.0*r*zernike(n-2,m,1) + zernike(n-2,m,0)) &
          +(factor2*r**2 - factor3)*zernike(n-2,m,2) - factor4*zernike(n-4,m,2))
        !write(0, *) n, m, r, zernike(n,m,:)
      enddo

      rm2 = rm1
      rm1 = rm
      rm = rm * r

    enddo
    do n = 2, lrad, 2
      zernike(n,0,0) = zernike(n,0,0) - (-1)**(n/2)
    enddo
    if (mpol >= 1) then
      do n = 3, lrad, 2
        zernike(n,1,0) = zernike(n,1,0) - (-1)**((n-1)/2) * real((n+1)/2) * r
        zernike(n,1,1) = zernike(n,1,1) - (-1)**((n-1)/2) * real((n+1)/2)
      enddo
    end if

    do m = 0, mpol
      do n = m, lrad, 2
        zernike(n,m,:) = zernike(n,m,:) / real(n+1)
      end do
    end do
  END SUBROUTINE get_zernike



END MODULE poly