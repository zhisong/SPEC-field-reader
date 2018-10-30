! this module defines the state of a SPEC run
module spec_state
  use field, only : Afield
  use geometry, only : interfaces

  type, public :: state

    integer :: Igeometry ! the degree of freedom (1=slab, 2=cylinder, 3=toroidal)
    integer :: Istellsym ! stellarator geometry (1=yes, others=no)
    integer :: Lfreebound ! free boundary? (0=no, 1=yes)
    integer :: Nvol  ! number of volume
    integer :: Nfp   ! the base toroidal harmonic (toroidal harmonic number can only be 0, Nfp, 2*Nfp, ...)
    integer :: Mpol  ! number of poloidal harmonics
    integer :: Ntor  ! number of toroidal harmonics
    integer, dimension(:), allocatable :: Lrad  ! the number of radial bases functions in each volume

    type(Afield), dimension(:), allocatable :: A


  end type state

end module