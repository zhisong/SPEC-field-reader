! this module defines the data structure of geometry
module geometry

  type, public :: volume

    real, dimension(:,:), allocatable :: Rbc, Rbs, Zbc, Zbs
    logical :: icoordinatesingularity
    logical :: isym
    integer, dimension(:), allocatable :: im, in
    integer :: mn, nfp, mvol

  end type volume

contains

  subroutine destroy_volume(v)
    implicit none
    type(volume) :: v
    if (ALLOCATED(v%im)) deallocate(v%im)
    if (ALLOCATED(v%im)) deallocate(v%im)
  end subroutine destroy_volume
  
end module geometry