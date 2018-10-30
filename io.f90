! the module reads the hdf5 file
module io

contains
  subroutine read_spec_h5(filename, sstate)
  ! read the SPEC hdf5 output
  ! INPUTS:
  ! filename - character(:), the hdf5 filename
  ! OUTPUTS:
  ! sstate
    use hdf5
    use spec_state
    implicit none
    character(LEN=*), intent(in) :: filename
    type(state) :: sstate

    integer(HID_T) :: file_id,  dset_id
    integer :: hdferr

    call h5open_f(hdferr)
    call check(hdferr)

    call h5fopen_f(TRIM(filename), H5F_ACC_RDONLY_F, file_id, hdferr)
    call check(hdferr)

    call read_h5_int(file_id, 'Igeometry', sstate%Igeometry)
    call read_h5_int(file_id, 'Istellsym', sstate%Istellsym)
    call read_h5_int(file_id, 'Lfreebound', sstate%Lfreebound)
    call read_h5_int(file_id, 'Nvol', sstate%Nvol)
    call read_h5_int(file_id, 'Nfp', sstate%Nfp)
    call read_h5_int(file_id, 'Mpol', sstate%Mpol)
    call read_h5_int(file_id, 'Ntor', sstate%Ntor)

    allocate(sstate%Lrad(sstate%Nvol))
    call read_h5_int_array(file_id, 'Lrad', sstate%Lrad, sstate%Nvol)

    call h5fclose_f(file_id, hdferr) 

  end subroutine read_spec_h5

  subroutine read_spec_field(filename, sstate)
  ! read the SPEC field output
  ! INPUTS:
  ! filename - character(:), the .A filename
  ! OUTPUTS:
  ! sstate
    use spec_state
    implicit none
    character(LEN=*), intent(in) :: filename
    type(state) :: sstate

    integer :: aunit
    integer :: Mvol, Mpol, Ntor, mn, Nfp, Lrad
    integer :: vvol, ii
    integer, dimension(:), allocatable :: im, in
    !!! be careful here, one should really care about SMALL/BIG ENDIAN here
    open(aunit, file=filename, action="read", form="unformatted" )

    read(aunit) Mvol, Mpol, Ntor, mn, Nfp

    if (sstate%Nvol .ne. Mvol) stop "please check consistency between .h5 and .A files"
    if (sstate%Mpol .ne. Mpol) stop "please check consistency between .h5 and .A files"
    if (sstate%Ntor .ne. Ntor) stop "please check consistency between .h5 and .A files"
    if (sstate%Nfp .ne. Nfp) stop "please check consistency between .h5 and .A files"

    allocate(im(mn))
    allocate(in(mn))
    allocate(sstate%A(Mvol))

    read(aunit) im(1:mn)
    read(aunit) in(1:mn)

    do vvol = 1, Mvol
      read(aunit) Lrad
      if (sstate%Lrad(vvol) .ne. Lrad) stop "please check consistency between .h5 and .A files"
      sstate%A(vvol)%lrad = Lrad

      sstate%A(vvol)%mn = mn
      sstate%A(vvol)%Nfp = Nfp
      allocate(sstate%A(vvol)%im(mn))
      allocate(sstate%A(vvol)%in(mn))
      sstate%A(vvol)%im(:) = im(:)
      sstate%A(vvol)%in(:) = in(:)

      if (sstate%Istellsym == 1) then
        sstate%A(vvol)%isym = .true.
      else
        sstate%A(vvol)%isym = .false.
      end if

      allocate(sstate%A(vvol)%Ate(0:Lrad, mn))
      allocate(sstate%A(vvol)%Aze(0:Lrad, mn))
      allocate(sstate%A(vvol)%Ato(0:Lrad, mn))
      allocate(sstate%A(vvol)%Azo(0:Lrad, mn))

      do ii = 1, mn ! loop over Fourier harmonics;

        read(aunit) sstate%A(vvol)%Ate(0:Lrad, ii)
        read(aunit) sstate%A(vvol)%Aze(0:Lrad, ii)
        read(aunit) sstate%A(vvol)%Ato(0:Lrad, ii)
        read(aunit) sstate%A(vvol)%Azo(0:Lrad, ii)
      
      enddo

    end do

    deallocate(im)
    deallocate(in)

    close(aunit)
  end subroutine read_spec_field

  subroutine read_h5_int(file_id, name, data)
    use hdf5
    implicit none
    integer(HID_T), intent(in) :: file_id
    character(LEN=*), intent(in) :: name
    integer, intent(out) ::  data

    integer(HID_T) :: dset_id
    integer :: hdferr
    integer, dimension(1) :: data_out
    integer(HSIZE_T), dimension(1) :: len

    len = (/1/)
    call h5dopen_f(file_id, name, dset_id, hdferr)
    call check(hdferr)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data_out, len, hdferr)
    call check(hdferr)
    data = data_out(1)
  end subroutine read_h5_int

  subroutine read_h5_int_array(file_id, name, data, nsize)
    use hdf5
    implicit none
    integer(HID_T), intent(in) :: file_id
    character(LEN=*), intent(in) :: name
    integer, dimension(nsize), intent(out) ::  data
    integer, intent(in) :: nsize

    integer(HID_T) :: dset_id
    integer :: hdferr
    integer, dimension(nsize) :: data_out
    integer(HSIZE_T), dimension(1) :: len

    if (nsize .ne. SIZE(data)) stop 'internal consistency when reading data'
    len = (/nsize/)
    call h5dopen_f(file_id, name, dset_id, hdferr)
    call check(hdferr)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data_out, len, hdferr)
    call check(hdferr)
    data(:) = data_out(:)
  end subroutine read_h5_int_array

  subroutine check(ierr)
    implicit none
    integer, intent(in) :: ierr
    if (ierr.ne.0) stop "hdf5 error!"
  end subroutine

end module io