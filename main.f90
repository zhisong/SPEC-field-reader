!! Test program
PROGRAM test_specio
    USE spec_state
    USE spec_io
    IMPLICIT NONE

    TYPE(state) :: s_spec
    CHARACTER(*), PARAMETER :: filename = "eq.spec.h5"

    CALL read_spec_h5(filename,s_spec)

END PROGRAM