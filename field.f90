! this module defines the data structure of A
module field

  type, public :: Afield

    real, dimension(:,:), allocatable :: Ate, Ato, Aze, Azo
    real, dimension(:), allocatable :: im, in
    integer :: lrad, mn, Nfp
    logical :: isym

  end type Afield

  type(Afield), dimension(:), allocatable :: A

contains
  ! subroutine getfield(s, theta, phi, a, gb)


  ! end subroutine getfield

  subroutine destroy_field(f)
    implicit none
    type(Afield) :: f
    ! free the space for field f
    if (ALLOCATED(f%Ate)) deallocate(f%Ate)
    if (ALLOCATED(f%Ato)) deallocate(f%Ato)
    if (ALLOCATED(f%Aze)) deallocate(f%Aze)
    if (ALLOCATED(f%Azo)) deallocate(f%Azo)
    if (ALLOCATED(f%im)) deallocate(f%im)
    if (ALLOCATED(f%im)) deallocate(f%in)
  end subroutine destroy_field

end module field