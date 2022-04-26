!! Test program
PROGRAM test_specio
    USE spec_state
    USE spec_io
    USE spec_field
    USE spec_geometry
    IMPLICIT NONE

    TYPE(state) :: s_spec
    TYPE(spec_metric) :: met
    CHARACTER(*), PARAMETER :: filename = "eq.spec.h5"

    ! INTEGER :: i
    REAL, DIMENSION(3) :: x = (/0.0, 0.0, 0.0/)
    INTEGER :: lvol = 1

    REAL, DIMENSION(3) :: A
    REAL, DIMENSION(3) :: B, H
    REAL, DIMENSION(3,3) :: dB

    REAL :: modB, modH

    ! OPEN(1,file='out.txt',status='scratch')
    ! WRITE(1,*) x(1), x(2), x(3), lvol
    ! CLOSE(1)


    ! PRINT*,"Input s, theta, zeta, lvol:"
    ! READ*, x(1), x(2) , x(3), lvol
    x = (/0.3,0.4,0.5/)
    lvol = 1


    CALL read_spec_h5(filename,s_spec)

    ! Get the metric info
    met = get_spec_metric(s_spec%Ri,lvol,x(1),x(2),x(3))

    ! Get the magnetic field and vector potential
    CALL get_spec_field(s_spec%A(lvol), x(1), x(2), x(3), s_spec%Mpol, A, B, dB)

    H = B/met%jac
    modB = SQRT(DOT_PRODUCT(B,B)/met%jac)
    modH = SQRT(DOT_PRODUCT(H,H)/met%jac)

    
    PRINT*,"R, Z : ",met%Rij(0,0),met%Zij(0,0)
    PRINT*,"jacobian: ",met%jac
    PRINT*,"metric: ",met%gij(1,:)
    PRINT*,met%gij(2,:)
    PRINT*,met%gij(3,:)

    PRINT*,"jac der: ",met%grad_jac

    PRINT*,"dgij (s): ",met%dgij(1,:,1)
    PRINT*,met%dgij(2,:,1)
    PRINT*,met%dgij(3,:,1)

    PRINT*,"dgij (t): ",met%dgij(1,:,2)
    PRINT*,met%dgij(2,:,2)
    PRINT*,met%dgij(3,:,2)

    PRINT*,"dgij (z): ",met%dgij(1,:,3)
    PRINT*,met%dgij(2,:,3)
    PRINT*,met%dgij(3,:,3)

    PRINT*,"modB: ",modB
    PRINT*,"modH: ",modH
    PRINT*,"B: ",B
    PRINT*,"H: ",H
    PRINT*,"dB: ",dB
    PRINT*,"A: ",A



END PROGRAM