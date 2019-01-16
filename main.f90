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
  real :: jac, djac(3), x(3), gij(3,3), dgij(3,3,3)
  real :: jac2, djac2(3), x2(3), gij2(3,3), dgij2(3,3,3)

  call read_spec_h5('G3V02L1Fi.001.sp.h5',ss)
  call read_spec_field('.G3V02L1Fi.001.sp.A',ss)
  
  lrad = ss%lrad(2)

  lvol = 2 ! in volume 2
  s = -0.3
  theta = 0.7
  xi = 2.3
  delta = 0.00001
  call get_spec_field(ss%A(lvol), s, theta, xi, a, gb, dgb)
  call get_spec_coord(ss%Ri, lvol, s, theta, xi, jac, djac, x, gij, dgij)
  ! need jac for each
  ! J_down or deriv of B
  
  write(*,*) 'Jacobian'
  write(*,*) jac
  write(*,*) 'derivatives direct'
  write(*,*) djac
  write(*,*) 'derivatives finite difference'
  s2 = s + delta
  theta2 = theta
  xi2 = xi
  call get_spec_coord(ss%Ri, lvol, s2, theta2, xi2, jac2, djac2, x2, gij2, dgij2)
  djac2(1) = (jac2 - jac) / delta
  dgij2(:,:,1) = (gij2 - gij2) / delta

  s2 = s 
  theta2 = theta + delta
  xi2 = xi
  call get_spec_coord(ss%Ri, lvol, s2, theta2, xi2, jac2, djac2, x2, gij2, dgij2)
  djac2(2) = (jac2 - jac) / delta
  dgij2(:,:,2) = (gij2 - gij2) / delta
  s2 = s 
  theta2 = theta
  xi2 = xi + delta
  call get_spec_coord(ss%Ri, lvol, s2, theta2, xi2, jac2, djac2, x2, gij2, dgij2)
  djac2(3) = (jac2 - jac) / delta
  dgij2(:,:,3) = (gij2 - gij) / delta

  write(*,*) djac2
  write(*,*)

  write(*,*) 'x'
  write(*,*) x
  write(*,*) 'gij'
  do ii = 1, 3
    write(*,*) gij(:,ii)
  enddo
  write(*,*) 'dgij direct'
  do ii = 1, 3
    do jj = 1,3
      write(*,*) dgij(:,jj,ii)
    end do
    write(*,*)
  end do

  write(*,*) 'dgij finite difference'
  do ii = 1, 3
    do jj = 1,3
      write(*,*) dgij2(:,jj,ii)
    end do
    write(*,*)
  end do

  write(*,*) 'A'
  write(*,*) a
  write(*,*) 'Jacobian * B'
  write(*,*) gb
  write(*,*)

  ! now we verify the derivatives
  ! the s derivatives
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