program SPEC_field_reader
  use spec_state
  use spec_field
  use spec_geometry
  use spec_io

  type(state) :: ss
  integer :: lrad
  real :: a(3), gb(3)
  real :: jac, x(3), gij(3,3), dgij(3,3,3)

  call read_spec_h5('G3V02L1Fi.001.sp.h5',ss)
  call read_spec_field('.G3V02L1Fi.001.sp.A',ss)
  
  lrad = ss%lrad(2)
  call get_spec_field(ss%A(2), lrad, -0.3, 0.7, 2.3, a, gb)
  call get_spec_coord(ss%Ri, 2, ss%mn, -0.3, 0.7, 2.3, jac, x, gij, dgij)
  
  write(*,*) jac
  write(*,*) x
  write(*,*)
  do ii = 1, 3
    write(*,*) gij(:,ii)
  enddo
  write(*,*)
  do ii = 1, 3
    do jj = 1,3
      write(*,*) dgij(:,jj,ii)
    end do
    write(*,*)
  end do

  call destroy_spec_state(ss)

end program SPEC_field_reader