program SPEC_field_reader
  use spec_state
  use spec_field
  use spec_geometry
  use spec_io

  type(state) :: ss
  integer :: lrad
  integer :: lvol
  real :: s, theta, xi
  real :: s2, theta2, xi2, delta
  real :: a(3), gb(3), dgb(3,3)
  real :: a2(3), gb2(3), dgb2(3,3)
  real :: jac, x(3), gij(3,3), dgij(3,3,3)

  call read_spec_h5('G3V02L1Fi.001.sp.h5',ss)
  call read_spec_field('.G3V02L1Fi.001.sp.A',ss)
  
  lrad = ss%lrad(2)

  lvol = 2 ! in volume 2
  s = -0.3
  theta = 0.7
  xi = 2.3
  call get_spec_field(ss%A(lvol), s, theta, xi, a, gb, dgb)
  call get_spec_coord(ss%Ri, lvol, s, theta, xi, jac, x, gij, dgij)
  ! need jac for each
  ! J_down or deriv of B
  
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
  write(*,*) a
  write(*,*) gb
  write(*,*)

  ! now we verify the derivatives
  ! the s derivatives
  delta = 0.00001
  s2 = s + delta
  theta2 = theta
  xi2 = xi
  call get_spec_field(ss%A(lvol), s2, theta2, xi2, a2, gb2, dgb2)
  dgb2(:,1) = (gb2 - gb) / delta
  write(*,*) 's derivatives of field, direct'
  write(*,*) dgb(:,1)
  write(*,*) 's derivatives of field, finite difference'
  write(*,*) dgb2(:,1)
  write(*,*)

  ! the theta derivatives
  delta = 0.00001
  s2 = s 
  theta2 = theta + delta
  xi2 = xi
  call get_spec_field(ss%A(lvol), s2, theta2, xi2, a2, gb2, dgb2)
  dgb2(:,2) = (gb2 - gb) / delta
  write(*,*) 'theta derivatives of field, direct'
  write(*,*) dgb(:,2)
  write(*,*) 'theta derivatives of field, finite difference'
  write(*,*) dgb2(:,2)
  write(*,*)

  ! the xi derivatives
  delta = 0.00001
  s2 = s 
  theta2 = theta 
  xi2 = xi + delta
  call get_spec_field(ss%A(lvol), s2, theta2, xi2, a2, gb2, dgb2)
  dgb2(:,3) = (gb2 - gb) / delta
  write(*,*) 'theta derivatives of field, direct'
  write(*,*) dgb(:,3)
  write(*,*) 'theta derivatives of field, finite difference'
  write(*,*) dgb2(:,3)

  call destroy_spec_state(ss)

end program SPEC_field_reader