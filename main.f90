program SPEC_field_reader
  use spec_state
  use io

  type(state) :: ss

  call read_spec_h5('G3V02L1Fi.001.sp.h5',ss)
  call read_spec_field('.G3V02L1Fi.001.sp.A',ss)
  write(*,*) ss%Ri%isym, ss%Ri%icoordinatesingularity

  call destroy_state(ss)

end program SPEC_field_reader