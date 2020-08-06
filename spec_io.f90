! Module for generating SPEC fields for LEVIS
! 
MODULE spec_io

PRIVATE

PUBLIC :: read_spec_h5

CONTAINS
    
    SUBROUTINE read_spec_h5(filename,sstate)
        USE hdf5
        USE spec_state ! The data structure for SPEC

        IMPLICIT NONE
        
        CHARACTER(len=*), INTENT(IN) :: filename
        TYPE(state) :: sstate

        INTEGER(HID_T) :: file_id, dset_id
        INTEGER :: hdferr
        INTEGER :: ivol, LradSize, Ll, Lu
        INTEGER, DIMENSION(:), ALLOCATABLE :: Lrad
        REAL, DIMENSION(:,:), ALLOCATABLE :: Ate, Ato, Aze, Azo

        ! CHECK HDF5 is working
        CALL h5open_f(hdferr)
        CALL checkh5(hdferr)

        ! Open the hdf5 file
        CALL h5fopen_f(TRIM(filename), H5F_ACC_RDONLY_F, file_id, hdferr)
        CALL checkh5(hdferr)
        
        !! Reading in base equilbrium data needed for reading/constructing all other data
        ! Input data
        CALL read_h5_int(file_id, 'input/physics/Igeometry', sstate%Igeometry)
        CALL read_h5_int(file_id, 'input/physics/Istellsym', sstate%Istellsym)
        CALL read_h5_int(file_id, 'input/physics/Lfreebound', sstate%Lfreebound)
        CALL read_h5_int(file_id, 'input/physics/Nvol', sstate%Nvol)
        CALL read_h5_int(file_id, 'input/physics/Nfp', sstate%Nfp)
        CALL read_h5_int(file_id, 'input/physics/Mpol', sstate%Mpol)
        CALL read_h5_int(file_id, 'input/physics/Ntor', sstate%Ntor)
        ! Output data
        CALL read_h5_int(file_id, 'output/Mvol', sstate%Mvol)
        CALL read_h5_int(file_id, 'output/mn', sstate%mn)
        
        ALLOCATE(sstate%Lrad(sstate%Mvol))
        CALL read_h5_int_array(file_id, 'input/physics/Lrad', sstate%Lrad, sstate%Mvol)
        
        !!!!!!!!!!!! VOLUME
        ! Harmonics
        sstate%Ri%Nfp = sstate%Nfp
        sstate%Ri%Mvol = sstate%Mvol
        sstate%Ri%icoordinatesingularity = (sstate%Igeometry .EQ. 2) .OR. (sstate%Igeometry .EQ. 3)
        sstate%Ri%Igeometry = sstate%Igeometry
        sstate%Ri%isym = (sstate%Istellsym .EQ. 1)
        sstate%Ri%mn = sstate%mn
        sstate%Ri%Mpol = sstate%Mpol
        sstate%Ri%Ntor = sstate%Ntor

        ALLOCATE(sstate%Ri%im(sstate%mn))
        ALLOCATE(sstate%Ri%in(sstate%mn))
        CALL read_h5_int_array(file_id, 'output/im', sstate%Ri%im, sstate%mn)
        CALL read_h5_int_array(file_id, 'output/in', sstate%Ri%in, sstate%mn)

        ! Amplitudes
        ALLOCATE(sstate%Ri%Rbc(sstate%mn, 0:sstate%Mvol))
        ALLOCATE(sstate%Ri%Rbs(sstate%mn, 0:sstate%Mvol))
        ALLOCATE(sstate%Ri%Zbc(sstate%mn, 0:sstate%Mvol))
        ALLOCATE(sstate%Ri%Zbs(sstate%mn, 0:sstate%Mvol))
        CALL read_h5_real_array(file_id, 'output/Rbc', sstate%Ri%Rbc, (/sstate%mn, sstate%Mvol+1/))
        CALL read_h5_real_array(file_id, 'output/Rbs', sstate%Ri%Rbs, (/sstate%mn, sstate%Mvol+1/))
        CALL read_h5_real_array(file_id, 'output/Zbc', sstate%Ri%Zbc, (/sstate%mn, sstate%Mvol+1/))
        CALL read_h5_real_array(file_id, 'output/Zbs', sstate%Ri%Zbs, (/sstate%mn, sstate%Mvol+1/))
        

        !!!!!!!!!!!! VECTOR POTENTIAL
        ! Allocate space
        LradSize = SIZE(sstate%Lrad)+SUM(sstate%Lrad)
        ALLOCATE(sstate%A(sstate%Mvol))

        ! Read in the vector potentials
        ALLOCATE(Ate(LradSize,sstate%mn))
        ALLOCATE(Aze(LradSize,sstate%mn))
        ALLOCATE(Ato(LradSize,sstate%mn))
        ALLOCATE(Azo(LradSize,sstate%mn))
        CALL read_h5_real_array(file_id, 'vector_potential/Ate', Ate, (/LradSize,sstate%mn/))
        CALL read_h5_real_array(file_id, 'vector_potential/Aze', Aze, (/LradSize,sstate%mn/))
        CALL read_h5_real_array(file_id, 'vector_potential/Ato', Ato, (/LradSize,sstate%mn/))
        CALL read_h5_real_array(file_id, 'vector_potential/Azo', Azo, (/LradSize,sstate%mn/))

        DO ivol = 1, sstate%Mvol
            ! Index
            Ll = SUM(sstate%Lrad(1:ivol-1))+SIZE(sstate%Lrad(1:ivol))
            Lu = SUM(sstate%Lrad(1:ivol))+SIZE(sstate%Lrad(1:ivol))
            ! Vector potential
            ALLOCATE(sstate%A(ivol)%Ate(0:sstate%Lrad(ivol), sstate%mn))
            ALLOCATE(sstate%A(ivol)%Aze(0:sstate%Lrad(ivol), sstate%mn))
            ALLOCATE(sstate%A(ivol)%Ato(0:sstate%Lrad(ivol), sstate%mn))
            ALLOCATE(sstate%A(ivol)%Azo(0:sstate%Lrad(ivol), sstate%mn))
            sstate%A(ivol)%Ate(:,:) = Ate(Ll:Lu,1:sstate%mn)
            sstate%A(ivol)%Aze(:,:) = Aze(Ll:Lu,1:sstate%mn)
            sstate%A(ivol)%Ato(:,:) = Ato(Ll:Lu,1:sstate%mn)
            sstate%A(ivol)%Azo(:,:) = Azo(Ll:Lu,1:sstate%mn)
            ! Lrad
            sstate%A(ivol)%Lrad = sstate%Lrad(ivol)

            ! im and in
            ALLOCATE(sstate%A(ivol)%im(sstate%mn))
            ALLOCATE(sstate%A(ivol)%in(sstate%mn))
            sstate%A(ivol)%im = (sstate%Ri%im)
            sstate%A(ivol)%in = (sstate%Ri%in)
            
            sstate%A(ivol)%Isingular = .FALSE. !Correct for this after loop
            sstate%A(ivol)%Isym = (sstate%Istellsym .EQ. 1)
            sstate%A(ivol)%Nfp = sstate%Nfp
            sstate%A(ivol)%mn = sstate%mn

        END DO
        ! Coodinate singularity
        sstate%A(1)%isingular = (sstate%Ri%icoordinatesingularity .EQV. .TRUE.)
        
        ! Deallocate temp variables
        DEALLOCATE(Ate)
        DEALLOCATE(Aze)
        DEALLOCATE(Ato)
        DEALLOCATE(Azo)

    END SUBROUTINE

    !!!! READING SUBROUTINES !!!!
    SUBROUTINE read_h5_int(file_id, name, data_out)
        ! USED FOR READING IN INTEGER VALUES
        USE HDF5
        IMPLICIT NONE
        INTEGER(HID_T), INTENT(IN) :: file_id
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER, DIMENSION(1) :: tmp_data
        INTEGER, INTENT(OUT) ::  data_out
    
        INTEGER(HID_T) :: dset_id
        INTEGER :: hdferr
        INTEGER(HSIZE_T), DIMENSION(1) :: len = (/1/)

        ! Open dataset
        CALL h5dopen_f(file_id, name, dset_id, hdferr)
        CALL checkh5(hdferr)
        ! Read dataset
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tmp_data, len, hdferr)
        CALL checkh5(hdferr)
        ! Close dataset
        CALL h5dclose_f(dset_id,hdferr)
        CALL checkh5(hdferr)
        ! Output data
        data_out = tmp_data(1)

    END SUBROUTINE read_h5_int

    SUBROUTINE read_h5_int_array(file_id, name, data_out, nsize)
        use hdf5
        IMPLICIT NONE
        INTEGER(HID_T), INTENT(IN) :: file_id
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER, DIMENSION(nsize), INTENT(OUT) ::  data_out
        INTEGER, INTENT(IN) :: nsize
    
        INTEGER(HID_T) :: dset_id
        INTEGER :: hdferr
        INTEGER(HSIZE_T), DIMENSION(1) :: len
    
        len = (/nsize/)
        ! Open dataset
        CALL h5dopen_f(file_id, name, dset_id, hdferr)
        CALL checkh5(hdferr)
        ! Read dataset
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, data_out, len, hdferr)
        CALL checkh5(hdferr)
        ! Close dataset
        CALL h5dclose_f(dset_id,hdferr)
        CALL checkh5(hdferr)

    END SUBROUTINE read_h5_int_array
    
    SUBROUTINE read_h5_real_array(file_id, name, data_out, nsize)
        USE hdf5
        IMPLICIT NONE
        INTEGER(HID_T), INTENT(IN) :: file_id
        CHARACTER(LEN=*), INTENT(IN) :: name
        REAL, DIMENSION(:,:), INTENT(OUT) ::  data_out
        INTEGER, DIMENSION(2),INTENT(IN) :: nsize
    
        INTEGER(HID_T) :: dset_id
        INTEGER :: hdferr
        INTEGER(HSIZE_T), DIMENSION(2) :: len
        
        len(:) = nsize(:)
        ! Open dataset
        CALL h5dopen_f(file_id, name, dset_id, hdferr)
        CALL checkh5(hdferr)
        ! Read dataset
        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, len, hdferr)
        CALL checkh5(hdferr)
        ! Close dataset
        CALL h5dclose_f(dset_id,hdferr)
        CALL checkh5(hdferr)

    END SUBROUTINE read_h5_real_array

    SUBROUTINE checkh5(ierr)
        ! Check to make sure no errors occured with hdf5 interface
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ierr
        IF (ierr.NE.0) THEN
            STOP "hdf5 error"
        END IF
    END SUBROUTINE

END MODULE


