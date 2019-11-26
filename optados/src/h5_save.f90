  !SUBROUTINE for saving array with HDF5

MODULE h5_save
  use HDF5
  use H5LT
  use od_io, only: filename_len, seedname, stdout
  use od_constants, only: dp
  public SAVE_OPTMAT

  CHARACTER(LEN=filename_len) :: filename ! File name
  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER        :: error
  LOGICAL        :: h5_initialized=.false.
  LOGICAL        :: h5_finalized=.false.
  CONTAINS

SUBROUTINE H5_INIT
  !! Initialized the HDF5 file interface and create the new file
  LOGICAL :: is_hdf5

  write(stdout,'(1x,a78)') 'initializing hdf5 interface'
  filename = trim(seedname)//".h5"
  CALL h5open_f(error)
  !
  ! Create a new file using default properties.
  !
  CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

  !INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
  !
  ! Create the dataspace.
  !CALL h5screate_simple_f(rank, dims, dspace_id, error)
  ! These are templates of the standard API
  !CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dspace_id, &
  !    dset_id, error)
  !CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, array_in, dims, error)

END SUBROUTINE H5_INIT

SUBROUTINE H5_FINALIZE
  ! This subroutine closes the file
  write(stdout,'(1x,a78)') 'Closing hdf5 interface'
  CALL h5fclose_f(file_id, error)
  CALL h5close_f(error)

END SUBROUTINE H5_FINALIZE

SUBROUTINE SAVE_MATRIX_WEIGHTS(weights)
  IMPLICIT NONE
  REAL(kind=dp), dimension(:,:,:,:,:), intent(in) :: weights
  CHARACTER(LEN=*), PARAMETER :: weight_array_name="opt_weight"
  INTEGER :: rank = 5
  INTEGER(HSIZE_T), dimension(5)       :: dims

  dims = int(shape(weights), kind=HSIZE_T)
  write(stdout,'(1x,a78)') 'Saving matrix weights'
  CALL h5ltmake_dataset_double_f(file_id, weight_array_name, rank, dims, weights, error)
  call h5ltset_attribute_string_f(file_id, weight_array_name, "unit", "eV", error)

END SUBROUTINE SAVE_MATRIX_WEIGHTS

SUBROUTINE SAVE_PDOS_ARRAY(weights)
  !! Subroutine for saving the PDOS array
  !! This saves the PDOS of orbitals 
  use od_dos_utils, only: E
  use od_electronic, only: efermi
  IMPLICIT NONE
  REAL(kind=dp), dimension(:,:,:), intent(in) :: weights
  CHARACTER(LEN=*), PARAMETER :: weight_array_name="dos_partial"
  CHARACTER(LEN=*), PARAMETER :: eng_array_name="engs"
  INTEGER :: rank = 3
  INTEGER(HSIZE_T), dimension(3)       :: dims
  REAL(kind=dp), dimension(1)       :: ef

  dims = int(shape(weights), kind=HSIZE_T)
  write(stdout,'(1x,a78)') 'Saving pdos weights'
  CALL h5ltmake_dataset_double_f(file_id, weight_array_name, rank, dims, weights, error)
  CALL SAVE_1D_DOUBLE_ARRAY(E, eng_array_name)
  ef=efermi
  call h5ltset_attribute_double_f(file_id, eng_array_name, "efermi", ef, int(1, kind=HSIZE_T), error)
END SUBROUTINE SAVE_PDOS_ARRAY

SUBROUTINE SAVE_RAW_WEIGHTS(weights)
  !! Subroutine for saving the PDOS weights
  !! This save the weight of each band for each orbital
  use od_dos_utils, only: E
  IMPLICIT NONE
  REAL(kind=dp), dimension(:,:,:,:), intent(in) :: weights
  CHARACTER(LEN=*), PARAMETER :: weight_array_name="raw_weights"
  INTEGER :: rank = 4
  INTEGER(HSIZE_T), dimension(4)       :: dims

  dims = int(shape(weights), kind=HSIZE_T)
  write(stdout,'(1x,a78)') 'Saving raw pdos weights'
  CALL h5ltmake_dataset_double_f(file_id, weight_array_name, rank, dims, weights, error)
END SUBROUTINE SAVE_RAW_WEIGHTS


SUBROUTINE SAVE_ORBITAL_INFO(pdos_symbol)
  !! Save definition of the project orbitals
  use od_electronic, only: pdos_orbital, orbitals
  use od_constants, only: periodic_table_name
  IMPLICIT NONE
  CHARACTER(len=3), allocatable, intent(in) :: pdos_symbol(:)
  INTEGER, allocatable :: orbital_species_atomic_number(:)   ! Array of the symbol for each orbital
  INTEGER, allocatable :: orbital_am_m(:)   ! Array of the angular momentum channel number
  INTEGER, allocatable :: orbital_am_count(:,:,:)   ! Aux array to construct the am_index
  INTEGER   :: num_orb
  INTEGER   :: loop, loop2, ns, nr, na

  num_orb = size(pdos_orbital%species_no)


  ! Construct an array for the Z of the specie for each orbital
  allocate(orbital_species_atomic_number(num_orb))
  do loop=1, size(pdos_orbital%species_no)
     do loop2=1, 109
        ! Find the atomic number of the symbol
        if (periodic_table_name(loop2)==pdos_symbol(pdos_orbital%species_no(loop))) then
           orbital_species_atomic_number(loop) = loop2
        end if
     end do
  end do

  ! Construct an array for the AM count for each orbital
  allocate(orbital_am_m(num_orb))
  allocate(orbital_am_count(maxval(pdos_orbital%species_no), maxval(pdos_orbital%rank_in_species), &
       maxval(pdos_orbital%am_channel) + 1))  ! am_channel starts from 0 

  ! Assign the m number for each am channel
  orbital_am_count = 0
  do loop=1, num_orb
     ns = pdos_orbital%species_no(loop)
     nr = pdos_orbital%rank_in_species(loop)
     na = pdos_orbital%am_channel(loop) + 1 ! am_channel starts from 0 but array index is from 1
     orbital_am_m(loop) = orbital_am_count(ns, nr, na)
     orbital_am_count(ns, nr, na) = orbital_am_count(ns, nr, na) + 1
  end do

  deallocate(orbital_am_count)

  CALL SAVE_1D_INTEGER_ARRAY(orbital_species_atomic_number, 'orbital_species_Z')
  CALL SAVE_1D_INTEGER_ARRAY(orbital_am_m, 'orbital_am_m')
  CALL SAVE_1D_INTEGER_ARRAY(pdos_orbital%rank_in_species, 'orbital_species_rank')
  CALL SAVE_1D_INTEGER_ARRAY(pdos_orbital%am_channel, 'orbital_am')

  deallocate(orbital_am_m)
  deallocate(orbital_species_atomic_number)

END SUBROUTINE SAVE_ORBITAL_INFO

SUBROUTINE SAVE_1D_INTEGER_ARRAY(array, name)
! SAVE an 1D array with givne name
  IMPLICIT NONE
  CHARACTER(LEN=*), intent(in) :: name
  INTEGER, dimension(:), intent(in):: array

  INTEGER(HSIZE_T), dimension(1) :: dims
  INTEGER :: rank=1

  dims = int(shape(array), kind=HSIZE_T)
  CALL h5ltmake_dataset_int_f(file_id, name, rank, dims, array, error)
END SUBROUTINE SAVE_1D_INTEGER_ARRAY

SUBROUTINE SAVE_1D_DOUBLE_ARRAY(array, name)
  ! SAVE an 1D array with givne name
  IMPLICIT NONE
  CHARACTER(LEN=*), intent(in) :: name
  REAL(kind=dp), dimension(:), intent(in):: array

  INTEGER(HSIZE_T), dimension(1) :: dims
  INTEGER :: rank=1

  dims = int(shape(array), kind=HSIZE_T)
  CALL h5ltmake_dataset_double_f(file_id, name, rank, dims, array, error)

END SUBROUTINE SAVE_1D_DOUBLE_ARRAY


SUBROUTINE SAVE_BAND_ENERGY
  ! Save energies
  use od_electronic, only: band_energy
  IMPLICIT NONE
  !REAL, dimension(jdos_nbins) :: engs
  CHARACTER(LEN=*), PARAMETER :: e_array_name="band_engs"
  INTEGER     ::   rank = 3                     ! Dataset rank
  INTEGER(HSIZE_T), dimension(3)     :: dims

  !engs = E
  write(stdout,'(1x,a78)') 'Saving band energy array'
  dims = int(shape(band_energy), kind=HSIZE_T)
  CALL h5ltmake_dataset_double_f(file_id, e_array_name, rank, dims, band_energy, error)
  call h5ltset_attribute_string_f(file_id, e_array_name, "unit", "eV", error)

END SUBROUTINE SAVE_BAND_ENERGY

SUBROUTINE SAVE_WEIGHTED_JDOS(weighted_jdos)
  ! Save energies
  use od_jdos_utils, only: E, jdos_nbins
  IMPLICIT NONE
  !REAL, dimension(jdos_nbins) :: engs
  CHARACTER(LEN=*), PARAMETER :: wjdos_name="weighted_jdos"
  CHARACTER(LEN=*), PARAMETER :: E_name="jdos_E"
  REAL(kind=dp), dimension(:,:,:), allocatable, intent(in) :: weighted_jdos

  INTEGER     ::   wjrank = 3                     ! Dataset rank
  INTEGER(HSIZE_T), dimension(3)     :: wjdims

  INTEGER     ::   erank = 1                   ! Dataset rank
  INTEGER(HSIZE_T), dimension(1)     :: Edims
  !engs = E
  write(stdout,'(1x,a78)') 'Saving weighted jdos'
  wjdims = int(shape(weighted_jdos), kind=HSIZE_T)
  Edims = int(jdos_nbins, kind=HSIZE_T)

  CALL h5ltmake_dataset_double_f(file_id, wjdos_name, wjrank, wjdims, weighted_jdos, error)
  CALL h5ltmake_dataset_double_f(file_id, E_name, erank, Edims, E, error)
  call h5ltset_attribute_string_f(file_id, E_name, "unit", "eV", error)

END SUBROUTINE SAVE_WEIGHTED_JDOS


SUBROUTINE SAVE_OPTMAT
  !! Example of saving rank 5 array
  !! Use either the LITE API or full API

  use od_electronic, only : optical_mat
  IMPLICIT NONE
  CHARACTER(LEN=*), PARAMETER :: rdsetname="opt_real" ! Dataset name
  CHARACTER(LEN=*), PARAMETER :: idsetname="opt_complex" ! Dataset name
  INTEGER(HSIZE_T), dimension(5)       :: dims

  !COMPLEX(8), dimension(dims(1), dims(2), dims(3), dims(4), dims(5))    :: array_in
  REAL(kind=dp), dimension(:,:,:,:,:), allocatable   :: array_real
  REAL(kind=dp), dimension(:,:,:,:,:), allocatable   :: array_imag

  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  !INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER     ::   error ! Error flag
  INTEGER     ::   i, j
  INTEGER     ::   rank = 5                        ! Dataset rank

  dims = int(shape(optical_mat), kind=HSIZE_T)

  allocate(array_real(dims(1), dims(2), dims(3), dims(4), dims(5)))
  allocate(array_imag(dims(1), dims(2), dims(3), dims(4), dims(5)))

  array_real = real(optical_mat)
  array_imag = aimag(optical_mat)
    !
  ! Create the dataspace.
  !CALL h5screate_simple_f(rank, dims, dspace_id, error)


  !CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dspace_id, &
  !    dset_id, error)

  !CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, array_in, dims, error)

  write(stdout,'(1x,a78)') 'Saving full optical matrix'
  CALL h5ltmake_dataset_double_f(file_id, rdsetname, rank, dims, array_real, error)
  CALL h5ltmake_dataset_double_f(file_id, idsetname, rank, dims, array_imag, error)

  deallocate(array_real)
  deallocate(array_imag)

  !CALL h5dclose_f(dset_id, error)

  END SUBROUTINE SAVE_OPTMAT

  END MODULE h5_save
