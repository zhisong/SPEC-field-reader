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

    REAL :: modB

    OPEN(1,file='out.txt',status='scratch')
    WRITE(1,*) x(1), x(2), x(3), lvol
    CLOSE(1)


    PRINT*,"Input s, theta, zeta, lvol:"
    READ*, x(1), x(2) , x(3), lvol


    CALL read_spec_h5(filename,s_spec)

    ! Get the metric info
    met = get_spec_metric(s_spec%Ri,lvol,x(1),x(2),x(3))

    ! Get the magnetic field and vector potential
    CALL get_spec_field(s_spec%A(lvol), x(1), x(2), x(3), s_spec%Mpol, A, B, dB)

    ! H = MATMUL(met%gij,B)/met%jac
    ! modB = SQRT(DOT_PRODUCT(B,H)/met%jac)
    modB = SQRT(DOT_PRODUCT(B,B))
    




END PROGRAM