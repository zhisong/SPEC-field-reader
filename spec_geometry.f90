! this module defines the data structure of geometry and the calculation of gij/dgij
MODULE spec_geometry
  IMPLICIT NONE

  TYPE, PUBLIC :: spec_metric
     REAL :: jac                       !< Jacobian
     REAL, DIMENSION(3) :: x           !<  the coordinates, for Igeometry==1 (Cartesian): (x, y, z)
     !<                                                         Igeometry==2 (Cylinder):  (r, theta, z),
     !<                                                         Igeometry==3 (Toroidal):  (R, phi, Z)
     REAL, DIMENSION(3) :: grad_jac    !< gradient of the Jacobian with respect (s, theta, xi)

     REAL, DIMENSION(3,3) :: gij       !< metric tensor (lower indices)
     REAL, DIMENSION(3,3) :: jacmat    !< Jacobian matrix, for Igeometry==1 (Cartesian): d (x,y,z)/d(s,theta,xi)
     !<                                                        Igeometry==2 (Cylinder):  d(r,theta,z)/d(s,theta,xi)
     !<                                                        Igeometry==3 (Toroidal):  d(R,phi,Z)/d(s,theta,xi)
     REAL, DIMENSION(3,3,3) :: dgij    !< dgij  the derivative of gij with respect to (s, theta, xi)
     !< first and second indices: i,j; third index: derivative

     REAL :: Rij(0:3,0:3), Zij(0:3,0:3)
  END TYPE spec_metric

  TYPE, PUBLIC :: volume
     REAL, DIMENSION(:,:), ALLOCATABLE :: Rbc, Rbs, Zbc, Zbs
     LOGICAL :: icoordinatesingularity
     LOGICAL :: isym
     INTEGER, DIMENSION(:), ALLOCATABLE :: im, in
     INTEGER :: mn, nfp, mvol, igeometry, mpol, ntor
  END TYPE volume

CONTAINS

  !> Obtain the coordinate quantities
  !> INPUTS:
  !> v     - TYPE(volume), the volume object read from file
  !> lvol  - INTEGER, which volume we are looking at,
  !> s     - s coordinate
  !> theta - theta coordinate
  !> xi    - xi coordinate
  !> RETURNS: a SPEC metric type
  ! 
  FUNCTION get_spec_metric(v, lvol, s, theta, xi)
    TYPE(spec_metric), TARGET, SAVE :: cache
    TYPE(spec_metric), POINTER :: get_spec_metric

    TYPE(volume),INTENT(IN) :: v
    INTEGER, INTENT(in) :: lvol
    REAL, INTENT(in) :: s, theta, xi
    REAL, SAVE :: s_save=-300.,theta_save=-300.,xi_save=-300.
    
    INTEGER :: ii, jj, kk

    IF (s .NE. s_save .OR. theta .NE. theta_save .OR. xi .NE. xi_save) THEN
       
       cache%gij=0.
       cache%jacmat=0.
       cache%dgij=0.
       
       SELECT CASE(v%igeometry)
       CASE (1) ! Cartesian
          CALL get_spec_derivatives(v,lvol,s,theta,xi,cache%rij)
          cache%x(1) = theta
          cache%x(2) = xi
          cache%x(3) = cache%Rij(0,0)

          cache%jac = cache%Rij(0,1)
          cache%grad_jac = cache%Rij(1:3,1)

          ! only non-vanishing terms
          cache%jacmat(1,2) = 1.
          cache%jacmat(2,3) = 1.
          cache%jacmat(3,1:3) = cache%Rij(0,1:3)

          DO ii = 1, 3
             DO jj = ii, 3
                cache%gij(jj,ii) = cache%gij(jj,ii) + cache%Rij(0,jj)*cache%Rij(0,ii)
                DO kk = 1, 3
                   cache%dgij(jj,ii,kk) = cache%dgij(jj,ii,kk) + cache%Rij(kk,jj)*cache%Rij(0,ii) + cache%Rij(0,jj)*cache%Rij(kk,ii)
                END DO ! kk
             END DO ! jj
          END DO ! ii

          cache%gij(2,2) = cache%gij(2,2) + 1.
          cache%gij(3,3) = cache%gij(3,3) + 1.

       CASE (2)  ! Cylindrical
          CALL get_spec_derivatives(v,lvol,s,theta,xi,cache%Rij)

          cache%x(1) = cache%Rij(0,0)
          cache%x(2) = theta
          cache%x(3) = xi

          cache%jac = cache%Rij(0,0)*cache%Rij(0,1)
          cache%grad_jac = cache%Rij(0,1:3)*cache%Rij(0,1) + cache%Rij(0,0)*cache%Rij(1:3,1)

          ! only non-vanishing terms
          cache%jacmat(1,1:3) = cache%Rij(0,1:3)
          cache%jacmat(2,2) = 1.
          cache%jacmat(3,3) = 1.

          DO ii = 1, 3
             DO jj = ii, 3
                cache%gij(jj,ii) = cache%gij(jj,ii) + cache%Rij(0,jj)*cache%Rij(0,ii)
                DO kk = 1, 3
                   cache%dgij(jj,ii,kk) = cache%dgij(jj,ii,kk) + cache%Rij(kk,jj)*cache%Rij(0,ii) + cache%Rij(0,jj)*cache%Rij(kk,ii)
                END DO ! kk
             END DO ! jj
          END DO ! ii

          cache%gij(2,2) = cache%gij(2,2) + cache%Rij(0,0)**2
          cache%gij(3,3) = cache%gij(3,3) + 1.
          DO ii = 1, 3
             cache%dgij(2,2,ii) = cache%dgij(2,2,ii) + 2.*cache%Rij(0,0)*cache%Rij(0,ii)
          END DO ! ii


       CASE (3) ! Toroidal
          CALL get_spec_derivatives(v,lvol,s,theta,xi,cache%Rij,cache%Zij)

          cache%x(1) = cache%Rij(0,0)
          cache%x(2) = xi
          cache%x(3) = cache%Zij(0,0)

          cache%jac = cache%Rij(0,0)*(cache%Zij(0,1)*cache%Rij(0,2) - cache%Rij(0,1)*cache%Zij(0,2))
          cache%grad_jac = cache%Rij(0,1:)*(cache%Zij(0,1)*cache%Rij(0,2) - cache%Rij(0,1)*cache%Zij(0,2))&
               + cache%Rij(0,0)*(cache%Zij(1:,1)*cache%Rij(0,2) - cache%Rij(1:,1)*cache%Zij(0,2))&
               + cache%Rij(0,0)*(cache%Zij(0,1)*cache%Rij(1:,2) - cache%Rij(0,1)*cache%Zij(1:,2))

          ! only non-vanishing terms
          cache%jacmat(1,1:3) = cache%Rij(0,1:3)
          cache%jacmat(2,3) = 1.
          cache%jacmat(3,1:3) = cache%Zij(0,1:3)

          DO ii = 1, 3
             DO jj = ii, 3
                cache%gij(jj,ii) = cache%gij(jj,ii) + cache%Rij(0,jj)*cache%Rij(0,ii) + cache%Zij(0,jj)*cache%Zij(0,ii)
                DO kk = 1, 3
                   cache%dgij(jj,ii,kk) = cache%dgij(jj,ii,kk)&
                        + cache%Rij(kk,jj)*cache%Rij(0,ii)&
                        + cache%Rij(0,jj)*cache%Rij(kk,ii)&
                        + cache%Zij(kk,jj)*cache%Zij(0,ii)&
                        + cache%Zij(0,jj)*cache%Zij(kk,ii)
                END DO ! kk
             END DO ! jj
          END DO ! ii

            cache%gij(3,3) = cache%gij(3,3) + cache%Rij(0,0)**2
            DO ii = 1, 3
               cache%dgij(3,3,ii) = cache%dgij(3,3,ii) + 2.*cache%Rij(0,0)*cache%Rij(0,ii)
            END DO ! ii


         CASE DEFAULT
            STOP 'Igeometry must be 1 to 3'
         END SELECT

         cache%gij(1,2) = cache%gij(2,1)
         cache%gij(1,3) = cache%gij(3,1)
         cache%gij(2,3) = cache%gij(3,2)

         !Mirror about the axis
         cache%dgij(1,2,1) = cache%dgij(2,1,1)
         cache%dgij(2,3,1) = cache%dgij(3,2,1)
         cache%dgij(1,3,1) = cache%dgij(3,1,1)

         cache%dgij(1,2,2) = cache%dgij(2,1,2)
         cache%dgij(2,3,2) = cache%dgij(3,2,2)
         cache%dgij(1,3,2) = cache%dgij(3,1,2)
         
         cache%dgij(1,2,3) = cache%dgij(2,1,3)
         cache%dgij(2,3,3) = cache%dgij(3,2,3)
         cache%dgij(1,3,3) = cache%dgij(3,1,3)
         
         s_save = s
         theta_save = theta
         xi_save = xi
      END IF
      get_spec_metric => cache

   ! PRINT*,"metric derivatives s",cache%dgij(:,:,1)
   ! PRINT*,"metric derivatives t",cache%dgij(:,:,2)
   ! PRINT*,"metric derivatives z",cache%dgij(:,:,3)

   END FUNCTION get_spec_metric

  SUBROUTINE get_spec_derivatives(v,lvol,s,theta,xi,Rij,Zij)
    TYPE(volume),INTENT(IN) :: v
    INTEGER, INTENT(in) :: lvol
    REAL, INTENT(in) :: s, theta, xi
    REAL,INTENT(OUT) :: Rij(0:3,0:3)
    REAL,INTENT(OUT),OPTIONAL :: zij(0:3,0:3)

    REAL :: sbar, alss, blss
    REAL, DIMENSION(v%mn) :: alphai, cosai, sinai
    REAL, DIMENSION(v%mn) :: t1, t2, t3, t4, ddt1, ddt3, fj, dfj, ddfj

    INTEGER :: ii

   alphai = v%im*theta - v%in*xi
   cosai = COS(alphai)
   sinai = SIN(alphai)

   IF (v%icoordinatesingularity .AND. lvol == 1) THEN
      ! See Z S Qu et al 2020 Plasma Phys. Control Fusion for the Zernike coordinate parameterisation
      ! PRINT*,"s=",s
         sbar = (1. + s) / 2.
      ! PRINT*,"s=",sbar

      !!!!
      fj(1:v%ntor+1) = sbar
      dfj(1:v%ntor+1) = 0.5
      ddfj(1:v%ntor+1) = 0.0

      fj(v%ntor+2:v%mn) = sbar**(v%im(v%ntor+2:v%mn)/2.0)
      dfj(v%ntor+2:v%mn) = 0.5 * (v%im(v%ntor+2:v%mn)/2.0) * fj(v%ntor+2:v%mn) / sbar
      ddfj(v%ntor+2:v%mn) = 0.5 * (v%im(v%ntor+2:v%mn)/2.0 - 1.0) * dfj(v%ntor+2:v%mn) / sbar

      !!!!
      !    fj(1)    = sbar
      !    dfj(1)   = 0.5
      !    ddfj(1)  = 0.

      !    fj(2:v%mpol+1) = sbar**(v%im(2:v%mpol+1)/2.)
      !    dfj(2:v%mpol+1) = (v%im(2:v%mpol+1)/4.) * sbar**(v%im(2:v%mpol+1)/2. - 1.)
      !    ddfj(2:v%mpol+1) = 0.

      ! IF (v%ntor.GT.0) THEN
      !    ddfj(v%mpol+2:v%mn) = sbar**(v%im(v%mpol+2:v%mn)/2. - 2.)
      !    dfj(v%mpol+2:v%mn) = ddfj(v%mpol+2:v%mn)*sbar*v%im(v%mpol+2:v%mn)/4.
      !    fj(v%mpol+2:v%mn) = ddfj(v%mpol+2:v%mn)*sbar**2
      !    ddfj(v%mpol+2:v%mn) = ddfj(v%mpol+2:v%mn)*v%im(v%mpol+2:v%mn)*(v%im(v%mpol+2:v%mn) - 2.)/16.
      ! END IF


      t1(:) = v%Rbc(:,0) + (v%Rbc(:,1) - v%Rbc(:,0))*fj(:)
      t2(:) = (v%Rbc(:,1) - v%Rbc(:,0))*dfj(:)
      ddt1(:) = (v%Rbc(:,1) - v%Rbc(:,0))*ddfj(:)
      IF (.NOT. v%isym) THEN
         t3(:) = v%Rbs(:,0) + (v%Rbs(:,1) - v%Rbs(:,0))*fj(:)
         t4(:) = (v%Rbs(:,1) - v%Rbs(:,0))*dfj(:)
         ddt3(:) = (v%Rbs(:,1) - v%Rbs(:,0))*ddfj(:)
      END IF

    ELSE !Use Chebychev basis for lvol > 1
       alss = 0.5*( 1. - s )
       blss = 0.5*( 1. + s )
       t1(:) = (alss*v%Rbc(:,lvol-1) + blss*v%Rbc(:,lvol))
       t2(:) = (-0.5*v%Rbc(:,lvol-1) + 0.5*v%Rbc(:,lvol))
       ddt1(:) = 0.
       IF (.NOT. v%isym) THEN
          t3(:) = (alss*v%Rbs(:,lvol-1) + blss*v%Rbs(:,lvol))
          t4(:) = (-0.5*v%Rbs(:,lvol-1) + 0.5*v%Rbs(:,lvol))
          ddt3(:) = 0.
       END IF
   END IF

    Rij(0,0) = SUM(t1*cosai)                  ! R
    Rij(0,1) = SUM(t2*cosai)                  ! dRds
    Rij(0,2) = SUM(t1*(-v%im*sinai))        ! dRdu
    Rij(0,3) = SUM(t1*(v%in*sinai))         ! dRdv
    Rij(1,1) = SUM(ddt1*cosai)                ! d2Rds2
    Rij(1,2) = SUM(t2*(-v%im*sinai))        ! d2Rdsdu
    Rij(1,3) = SUM(t2*(v%in*sinai))         ! d2Rdsdv
    Rij(2,2) = SUM(t1*(-v%im**2*cosai))     ! d2Rdu2
    Rij(2,3) = SUM(t1*(v%im*v%in*cosai))  ! d2Rdudv
    Rij(3,3) = SUM(t1*(-v%in**2*cosai))     ! d2Rdv2

    IF (.NOT. v%isym) THEN
       Rij(0,0) = Rij(0,0) + SUM(t3*sinai)
       Rij(0,1) = Rij(0,1) + SUM(t4*sinai)
       Rij(0,2) = Rij(0,2) + SUM(t3*(v%im*cosai))
       Rij(0,3) = Rij(0,3) + SUM(t3*(-v%in*cosai))
       Rij(1,1) = Rij(1,1) + SUM(ddt3*sinai)
       Rij(1,2) = Rij(1,2) + SUM(t4*(v%im*cosai))
       Rij(1,3) = Rij(1,3) + SUM(t4*(-v%in*cosai))
       Rij(2,2) = Rij(2,2) + SUM(t3*(-v%im**2*sinai))
       Rij(2,3) = Rij(2,3) + SUM(t3*(v%im*v%in*sinai))
       Rij(3,3) = Rij(3,3) + SUM(t3*(-v%in**2*sinai))
    END IF

    Rij(2,1) = Rij(1,2)
    Rij(3,1) = Rij(1,3)
    Rij(3,2) = Rij(2,3)

    IF(PRESENT(zij)) THEN
       IF (v%icoordinatesingularity .AND. lvol == 1) THEN
          t1(:) = v%Zbs(:,0) + (v%Zbs(:,1) - v%Zbs(:,0))*fj(:)
          t2(:) = (v%Zbs(:,1) - v%Zbs(:,0))*dfj(:)
          ddt1(:) = (v%Zbs(:,1) - v%Zbs(:,0))*ddfj(:)
          IF (.NOT. v%isym) THEN
             t3(:) = v%Zbc(:,0) + (v%Zbc(:,1) - v%Zbc(:,0))*fj(:)
             t4(:) = (v%Zbc(:,1) - v%Zbc(:,0))*dfj(:)
             ddt3(:) = (v%Zbc(:,1) - v%Zbc(:,0))*ddfj(:)
          END IF

       ELSE
          alss = 0.5*( 1. - s )
          blss = 0.5*( 1. + s )
          t1(:) = (alss*v%Zbs(:,lvol-1) + blss*v%Zbs(:,lvol))
          t2(:) = (-0.5*v%Zbs(:,lvol-1) + 0.5*v%Zbs(:,lvol))
          ddt1(:) = 0.
          IF (.NOT. v%isym) THEN
             t3(:) = (alss*v%Zbc(:,lvol-1) + blss*v%Zbc(:,lvol))
             t4(:) = (-0.5*v%Zbc(:,lvol-1) + 0.5*v%Zbc(:,lvol))
             ddt3(:) = 0.
          END IF
       END IF

       Zij(0,0) = SUM(t1*sinai)
       Zij(0,1) = SUM(t2*sinai)
       Zij(0,2) = SUM(t1*(v%im*cosai))
       Zij(0,3) = SUM(t1*(-v%in*cosai))
       Zij(1,1) = SUM(ddt1*sinai)
       Zij(1,2) = SUM(t2*(v%im*cosai))
       Zij(1,3) = SUM(t2*(-v%in*cosai))
       Zij(2,2) = SUM(t1*(-v%im**2*sinai))
       Zij(2,3) = SUM(t1*(v%im*v%in*sinai))
       Zij(3,3) = SUM(t1*(-v%in**2*sinai))

       IF (.NOT. v%isym) THEN
          Zij(0,0) = Zij(0,0) + SUM(t3*cosai)
          Zij(0,1) = Zij(0,1) + SUM(t4*cosai)
          Zij(0,2) = Zij(0,2) + SUM(t3*(-v%im*sinai))
          Zij(0,3) = Zij(0,3) + SUM(t3*(v%in*sinai))
          Zij(1,1) = Zij(1,1) + SUM(ddt3*cosai)
          Zij(1,2) = Zij(1,2) + SUM(t4*(-v%im*sinai))
          Zij(1,3) = Zij(1,3) + SUM(t4*(v%in*sinai))
          Zij(2,2) = Zij(2,2) + SUM(t3*(-v%im**2*cosai))
          Zij(2,3) = Zij(2,3) + SUM(t3*(v%im*v%in*cosai))
          Zij(3,3) = Zij(3,3) + SUM(t3*(-v%in**2*cosai))
       END IF

       Zij(2,1) = Zij(1,2)
       Zij(3,1) = Zij(1,3)
       Zij(3,2) = Zij(2,3)
    END IF
  END SUBROUTINE get_spec_derivatives


  SUBROUTINE destroy_volume(v)
    IMPLICIT NONE
    TYPE(volume) :: v
    IF (ALLOCATED(v%im)) DEALLOCATE(v%im)
    IF (ALLOCATED(v%im)) DEALLOCATE(v%im)
  END SUBROUTINE destroy_volume

END MODULE spec_geometry
