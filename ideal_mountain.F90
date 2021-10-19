!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

module ideal_mtn_mod
  !-------------------------------------------------------------------------------
  !>@brief The module 'multi_vtx_mod' contains members related to
  !! the idealized simulation of warm-core vortices (tropical cyclones, TCs).
  !>@details This module includes public routines that read the test-case
  !! namelist, and establish the initial state of the model. Options regarding
  !! the structure(s) of the initial vortices are availble, including those
  !! regulating initial location, intensity, size, and thermal structure.
  !! This module is based on test_case == -55 in test_cases_mod (DCMIP16_TC).
  ! <table>
  !   <tr>
  !     <th>Module Name</th>
  !     <th>Functions Included</th>
  !   </tr>
  !   <tr>
  !     <td>test_control_mod</td>
  !   </tr>
  !   <tr>
  !     <td>constants_mod</td>
  !     <td>radius, pi_8, grav</td>
  !   </tr>
  !   <tr>
  !     <td>fms_mod</td>
  !     <td>check_nml_error</td>
  !   </tr>
  !   <tr>
  !     <td>fv_arrays_mod</td>
  !     <td>fv_grid_type, R_GRID</td>
  !   </tr>
  !   <tr>
  !     <td>mpp_domains_mod</td>
  !     <td>domain2d</td>
  !   </tr>
  !   <tr>
  !     <td>init_hydro_mod</td>
  !     <td>hydro_eq</td>
  !   </tr>
  !   <tr>
  !     <td>fv_grid_utils_mod</td>
  !     <td>great_circle_dist</td>
  !   </tr>
  !   <tr>
  !     <td>fv_mp_mod</td>
  !     <td>is_master</td>
  !   </tr>
  !   <tr>
  !     <td>mpp_mod</td>
  !     <td>input_nml_file, mpp_error, FATAL, WARNING, NOTE</td>
  !   </tr>
  ! </table>

  use test_control_mod
  use constants_mod,      only: const_A=>radius, pi=>pi_8, grav
  use fms_mod,            only: check_nml_error
  use fv_arrays_mod,      only: fv_grid_type, R_GRID
  use mpp_domains_mod,    only: domain2d
  use init_hydro_mod,     only: hydro_eq
  use fv_grid_utils_mod,  only: great_circle_dist
  use fv_mp_mod,          only: is_master
  use mpp_mod,            only: input_nml_file, mpp_error, FATAL, WARNING, NOTE

  implicit none
  private

  !---------------!
  ! ---Options--- !
  !---------------!
  !> Required variables to set in the namelist are indicated with a (V).
  !> Recommended variables to set are indicated with an (R).
  !> Optional variables to set are indicated with an (O).
  !> Discouraged variables to set are indicated with an (N).

  real :: peak_hgt   !< (R) Height of the peak of the idealized mountain (m).
  real :: base_size  !< (R) Size of the mountain's base (km).
  real :: lon        !< (R) Central longitude of the mountain's peak
                     !!     (degrees positive east of the Prime Meridian).
  real :: lat        !< (R) Central latitude of the mountain's peak
                     !!     (degrees positive north of the Equator).

  !----------------!
  ! ---Routines--- !
  !----------------!

  public  :: read_ideal_mtn
  public  :: intern_init_ideal_mtn
  public  :: init_ideal_mtn
  public  :: ic_ideal_mtn

contains

  subroutine read_ideal_mtn(fname)
    implicit none
    character(*), intent(in) :: fname
    integer :: ierr, f_unit, ios !, unit

    namelist /test_case_nml/ &
      peak_hgt, base_size, lon, lat

    ! unit = stdlog()

    ! Default settings
    peak_hgt   = 5960.
    base_size  = 2000.
    lon        = 0.
    lat        = 0.

    ! Read user-defined values in namelist (specified values replace default)
#ifdef INTERNAL_FILE_NML
    read(input_nml_file,test_case_nml,iostat=ios)
    ierr = check_nml_error(ios,'test_case_nml')
#else
    f_unit = open_namelist_file(fname)
    rewind(f_unit)
    read(f_unit,test_case_nml,iostat=ios)
    ierr = check_nml_error(ios,'test_case_nml')
    call close_file(f_unit)
#endif
    ! write(unit, nml=test_case_nml)

    if(is_master()) then
      write(*,*) '(ideal_mtn,read) peak_hgt=',peak_hgt
      write(*,*) '(ideal_mtn,read) base_size=',base_size
      write(*,*) '(ideal_mtn,read) lon=',lon
      write(*,*) '(ideal_mtn,read) lat=',lat
    endif

  end subroutine read_ideal_mtn

  !-------------------------------------------------------------------------------
  !>@brief intern_nit_multi_vtx reads an internal namelist to set case settings.
  !>@details The case's namelist is declared and its members are initialized with
  !!         default values. The FV3 input namelist file is read, which overwrites
  !!         default settings. This allows user-friendly modifications to the
  !!         simulation parameters without necessitating recompilation.
  integer function intern_init_ideal_mtn()
    implicit none
    integer :: ios

    namelist /test_case_nml/ &
      peak_hgt, base_size, lon, lat

    ! Default settings
    peak_hgt   = 5960.
    base_size  = 2000.
    lon        = 0.
    lat        = 0.

    ! Read user-defined values in namelist
    read(input_nml_file,test_case_nml,iostat=ios)
    intern_init_ideal_mtn = check_nml_error(ios,'test_case_nml')

    if(is_master()) then
      write(*,*) '(ideal_mtn,init) peak_hgt=',peak_hgt
      write(*,*) '(ideal_mtn,init) base_size=',base_size
      write(*,*) '(ideal_mtn,init) lon=',lon
      write(*,*) '(ideal_mtn,init) lat=',lat
      write(*,*) '(ideal_mtn) exit intern_init_ideal_mtn=',intern_init_ideal_mtn
    endif

  end function intern_init_ideal_mtn

  !-------------------------------------------------------------------------------
  !>@brief init_multi_vtx reads a namelist to initialize case settings.
  !>@details The case's namelist is declared and its members are initialized with
  !!         default values. The FV3 input namelist file is read, which overwrites
  !!         default settings. This allows user-friendly modifications to the
  !!         simulation parameters without necessitating recompilation.
  integer function init_ideal_mtn(f_unit)
    implicit none
    integer, intent(in) :: f_unit
    integer :: ios

    namelist /test_case_nml/ &
      peak_hgt, base_size, lon, lat

    ! Default settings
    peak_hgt   = 5960.
    base_size  = 2000.
    lon        = 0.
    lat        = 0.

    ! Read user-defined values in namelist (specified values replace default)
    read(f_unit,test_case_nml,iostat=ios)
    init_ideal_mtn = check_nml_error(ios,'test_case_nml')

    if(is_master()) then
      write(*,*) '(ideal_mtn,init) peak_hgt=',peak_hgt
      write(*,*) '(ideal_mtn,init) base_size=',base_size
      write(*,*) '(ideal_mtn,init) lon=',lon
      write(*,*) '(ideal_mtn,init) lat=',lat
      write(*,*) '(ideal_mtn) exit init_ideal_mtn=',init_ideal_mtn
    endif

  end function init_ideal_mtn

  !-------------------------------------------------------------------------------
  !>@brief ic_multi_vtx initializes the model state with idealized warm-core cyclones.
  !>@details The model's initial state is set in this routine, using the parameters
  !!         specified in the user's namelist (default settings if none specified).
  subroutine ic_ideal_mtn(u, v, w, pt, delp, delz, q, &
                          phis, ps, ak, bk, gridstruct, domain)

    implicit none
    type(fv_grid_type), target, intent(IN) :: gridstruct
    type(domain2d), intent(IN) :: domain
    real, intent(IN), dimension(npz+1) :: ak, bk

    real, intent(OUT), dimension(isd:ied,jsd:jed+1,npz) :: u  !< abscissa-sided
    real, intent(OUT), dimension(isd:ied+1,jsd:jed,npz) :: v  !< ordinate-sided
    real, intent(OUT), dimension(isd:ied,jsd:jed,npz)   :: w, pt, delp !, delz  !< cell-centered
    real, intent(OUT), dimension(is:,js:,1:) :: delz
    real, intent(INOUT), dimension(isd:ied,jsd:jed,npz,ncnst) :: q !< cell-centered, dim-4 for species
    real, intent(OUT), dimension(isd:ied,jsd:jed) :: phis, ps

    real(kind=R_GRID), pointer, dimension(:,:,:) :: agrid
    real(kind=R_GRID), pointer, dimension(:,:)   :: area

    real(kind=R_GRID), dimension(2) :: mtn_center        !< mountain center coordinates
    real(kind=R_GRID), parameter    :: earth_r = const_A !< Earth's radius

    real :: peak_gz   !< geopotential (gpm)
    real :: r         !< working distance from mountain center

    integer :: i,j

    if (peak_hgt < 0.0) then
      call mpp_error(NOTE, '[ERROR] (ideal_mountain.F90, ic_ideal_mtn):')
      call mpp_error(NOTE, '    Negative mountain height specified.')
      ! return/stop
      call mpp_error(FATAL, 'FATAL, ic_ideal_mtn: Invalid setting of peak_hgt.')
    endif

    if (base_size <= 0.0) then
      call mpp_error(NOTE, '[ERROR] (ideal_mountain.F90, ic_ideal_mtn):')
      call mpp_error(NOTE, '    Non-positive mountain size specified.')
      ! return/stop
      call mpp_error(FATAL, 'FATAL, ic_ideal_mtn: Invalid setting of base_size.')
    endif

    agrid => gridstruct%agrid_64
    area  => gridstruct%area_64

    q(:,:,:,:) = 3.e-6
    u(:,:,:)   = 0.0
    v(:,:,:)   = 0.0
    phis       = 0.0
    if (.not.hydrostatic) w(:,:,:) = 0.0

    peak_gz    = grav * peak_hgt
    mtn_center = (/ lon * pi / 180.0, lat * pi / 180.0 /)

    if(is_master()) then
      write(*,*) '(ideal_mtn,ic) is=',is
      write(*,*) '(ideal_mtn,ic) ie=',ie
      write(*,*) '(ideal_mtn,ic) js=',js
      write(*,*) '(ideal_mtn,ic) je=',je
      write(*,*) '(ideal_mtn,ic) is2=',is2
      write(*,*) '(ideal_mtn,ic) ie2=',ie2
      write(*,*) '(ideal_mtn,ic) js2=',js2
      write(*,*) '(ideal_mtn,ic) je2=',je2
    endif

    do j=js2,je2
      do i=is2,ie2
        r = great_circle_dist(agrid(i,j,:), mtn_center, earth_r) / 1000.0
        if (base_size > r) then
          phis(i,j) = peak_gz * (1.0 - (r/base_size))
        else
          phis(i,j) = 0.0
        endif
      enddo
    enddo
    if (is_master()) then
      write(*,*) '(ideal_mtn,ic) phis=',phis
    endif
    call hydro_eq(npz, is, ie, js, je, ps, phis, dry_mass, &
              delp, ak, bk, pt, delz, area, ng, mountain, hydrostatic, hybrid_z, domain)
!    call hydro_eq(npz, is, ie, js, je, ps, phis, 1.E5, &
!              delp, ak, bk, pt, delz, area, ng, .false., hydrostatic, hybrid_z, domain)

    nullify(agrid)
    nullify(area)

  end subroutine ic_ideal_mtn

end module ideal_mtn_mod
