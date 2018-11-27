! this module defines the state of a SPEC run
module spec_state
  use spec_field, only : Afield
  use spec_geometry, only : volume

  type, public :: state

    integer :: Igeometry ! the degree of freedom (1=slab, 2=cylinder, 3=toroidal)
    integer :: Istellsym ! stellarator geometry (1=yes, others=no)
    integer :: Lfreebound ! free boundary? (0=no, 1=yes)
    integer :: Nvol  ! number of volume (without freebound)
    integer :: Mvol  ! number of volume (with freebound)
    integer :: Nfp   ! the base toroidal harmonic (toroidal harmonic number can only be 0, Nfp, 2*Nfp, ...)
    integer :: Mpol  ! number of poloidal harmonics
    integer :: Ntor  ! number of toroidal harmonics
    integer :: mn    ! total number of poloidal + toroidal harmonics
    !integer, dimension(:), allocatable :: im, in ! m and n harmonics
    integer, dimension(:), allocatable :: Lrad  ! the number of radial bases functions in each volume

    type(Afield), dimension(:), allocatable :: A
    type(volume) :: Ri

  end type state

contains

  subroutine destroy_spec_state(ss)
    use spec_field, only : destroy_field
    use spec_geometry, only : destroy_volume
    implicit none
    type(state) :: ss
    integer :: i1

    if (ALLOCATED(ss%Lrad)) deallocate(ss%Lrad)

    do i1 = 1, ss%Mvol
      call destroy_field(ss%A(i1))
    end do
    if (ALLOCATED(ss%A)) deallocate(ss%A)

    call destroy_volume(ss%Ri)

  end subroutine destroy_spec_state

end module