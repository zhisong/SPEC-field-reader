program SPEC_field_reader
  use spec_state
  use spec_field
  use spec_geometry
  use spec_io

  type(state) :: ss
  integer :: lrad
  integer :: lvol
  real :: s, theta, xi
  real :: delta
  real :: a(3), gb(3), dgb(3,3)
  real :: a2(3), gb2(3), dgb2(3,3)
  real :: jac, djac(3), x(3), gij(3,3),jacmat(3,3)
  
  TYPE(spec_metric), POINTER :: met
  
  integer :: ii,jj

  call read_spec_h5('G3V02L1Fi.001.sp.h5',ss)
  call read_spec_field('.G3V02L1Fi.001.sp.A',ss)
  
  lrad = ss%lrad(2)

  lvol = 2 ! in volume 2
  s = -0.3
  theta = 0.7
  xi = 2.3
  delta = 0.00001
  call get_spec_field(ss%A(lvol), s, theta, xi, a, gb, dgb)
  met=>get_spec_metric(ss%Ri, lvol, s, theta, xi)

  ! save for comparison
  jac=met%jac
  x=met%x
  gij=met%gij
  
  write(*,*) 'Jacobian'
  write(*,*) met%jac
  write(*,*) 'derivatives direct'
  write(*,*) met%grad_jac
  write(*,*) 'jacobian matrix direct'
  write(*,*) met%jacmat(1,:)
  write(*,*) met%jacmat(2,:)
  write(*,*) met%jacmat(3,:)
  write(*,*) 'derivatives finite difference'

  ! s deriv
  met=>get_spec_metric(ss%Ri, lvol, s+delta, theta, xi)
  djac(1) = (met%jac - jac) / delta
  jacmat(:,1) = (met%x - x)/delta

  ! theta deriv
  met=>get_spec_metric(ss%Ri, lvol, s, theta+delta, xi)
  djac(2) = (met%jac - jac) / delta
  jacmat(:,2) = (met%x - x)/delta

  ! xi deriv
  met=>get_spec_metric(ss%Ri, lvol, s, theta, xi+delta)
  djac(3) = (met%jac - jac) / delta
  jacmat(:,3) = (met%x- x)/delta

  write(*,*) djac
  write(*,*)
  write(*,*) jacmat(1,:)
  write(*,*) jacmat(2,:)
  write(*,*) jacmat(3,:)
  write(*,*)

  write(*,*) 'x'
  write(*,*) x
  write(*,*) 'gij'
  do ii = 1, 3
    write(*,*) gij(:,ii)
  enddo

  write(*,*) 'A'
  write(*,*) a
  write(*,*) 'Jacobian * B'
  write(*,*) gb
  write(*,*)

  ! now we verify the derivatives
  ! the s derivatives
  call get_spec_field(ss%A(lvol), s+delta, theta, xi, a2, gb2, dgb2)
  dgb2(:,1) = (gb2 - gb) / delta
  write(*,*) 's derivatives of field, direct'
  write(*,*) dgb(:,1)
  write(*,*) 's derivatives of field, finite difference'
  write(*,*) dgb2(:,1)
  write(*,*)

  ! the theta derivatives
  call get_spec_field(ss%A(lvol), s, theta+delta, xi, a2, gb2, dgb2)
  dgb2(:,2) = (gb2 - gb) / delta
  write(*,*) 'theta derivatives of field, direct'
  write(*,*) dgb(:,2)
  write(*,*) 'theta derivatives of field, finite difference'
  write(*,*) dgb2(:,2)
  write(*,*)

  ! the xi derivatives
  call get_spec_field(ss%A(lvol), s, theta, xi+delta, a2, gb2, dgb2)
  dgb2(:,3) = (gb2 - gb) / delta
  write(*,*) 'theta derivatives of field, direct'
  write(*,*) dgb(:,3)
  write(*,*) 'theta derivatives of field, finite difference'
  write(*,*) dgb2(:,3)

  call destroy_spec_state(ss)

end program SPEC_field_reader
