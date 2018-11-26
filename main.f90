program SPEC_field_reader
  use spec_state
  use field
  use io

  type(state) :: ss
  integer :: lrad
  real :: a(3), gb(3)

  call read_spec_h5('G3V02L1Fi.001.sp.h5',ss)
  call read_spec_field('.G3V02L1Fi.001.sp.A',ss)
  
  lrad = ss%lrad(2)
  call get_field(ss%A(2), lrad, -0.3, 0.7, 2.3, a, gb)
  write(*,*) ss%A(2)%isingular
  write(*,*) a
  write(*,*) gb

  call destroy_state(ss)

end program SPEC_field_reader