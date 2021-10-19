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

module multi_vtx_mod
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
  !     <td>radius, pi_8, omega, grav, rdgas, kappa</td>
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
  !     <td>great_circle_dist, mid_pt_sphere, get_unit_vect2, get_latlon_vector,
  !         inner_prod</td>
  !   </tr>
  !   <tr>
  !     <td>tracer_manager_mod</td>
  !     <td>get_tracer_index</td>
  !   </tr>
  !   <tr>
  !     <td>field_manager_mod</td>
  !     <td>MODEL_ATMOS</td>
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
  use constants_mod,      only: const_A=>radius, pi=>pi_8, omega, grav, rdgas, kappa
  use fms_mod,            only: check_nml_error
  use fv_arrays_mod,      only: fv_grid_type, R_GRID
  use mpp_domains_mod,    only: domain2d
  use init_hydro_mod,     only: hydro_eq
  use fv_grid_utils_mod,  only: great_circle_dist, mid_pt_sphere, get_unit_vect2, &
                                get_latlon_vector, inner_prod
  use tracer_manager_mod, only: get_tracer_index
  use field_manager_mod,  only: MODEL_ATMOS
  use fv_mp_mod,          only: is_master
  use mpp_mod,            only: input_nml_file, mpp_error, FATAL, WARNING, NOTE

  implicit none
  private
  integer, parameter :: MAX_VORTEX = 16 !< The maximum number of initial vortices
                                        !! allowed (arbitrary limit).

  !---------------!
  ! ---Options--- !
  !---------------!
  !> Required variables to set in the namelist are indicated with a (V).
  !> Recommended variables to set are indicated with an (R).
  !> Optional variables to set are indicated with an (O).
  !> Discouraged variables to set are indicated with an (N).

  real, dimension(MAX_VORTEX) :: lon    !< (V) Central longitude(s) of the vortex(vortices)
                                        !!     in degrees (positive east from Prime Meridian).
  real, dimension(MAX_VORTEX) :: lat    !< (V) Central latitude(s) of the vortex(vortices)
                                        !!     in degrees (positive north from Equator).
  real, dimension(MAX_VORTEX) :: dp     !< (R) Central surface pressure depression relative to
                                        !!     the background state (in Pa).
  real, dimension(MAX_VORTEX) :: rsize  !< (R) A size parameter (km) governing the
                                        !!     radial structure of the vortex(vortices).
                                        !!     Acts as an e-folding scale, where r=rsize is
                                        !!     where the the surface pressure depression is
                                        !!     dp/e .

  integer :: num_vortex !< (V) The number of vortices to initialize.
  logical :: aqua       !< (R) Whether to initialize the model as an aquaplanet.
  real    :: q_0        !< (O) Specific humidity (q) at the surface (kg/kg).
  real    :: z_q1       !< (O) A scaling parameter (m) for the vertical profile of q.
                        !!     Reduces q from q_0 with exp(-z/z_q1).
  real    :: z_q2       !< (O) A scaling parameter (m) for the vertical profile of q.
                        !!     Reduces q from q_0 with exp[-(z/z_q2)^2].
                        !!     By default, z_q1 < z_q2 .
  real    :: q_trop     !< (O) Specific humidity at the tropopause and above.
  real    :: T_0        !< (O) Surface temperature (K).
  real    :: lapse_rate !< (O) Atmospheric lapse rate, -dT/dz (K/m).
  real    :: p_0        !< (O) Background/base surface pressure (Pa).
  real    :: z_trop     !< (O) Tropopause height (m).
  real    :: z_p        !< (O) Shape parameter (m) for the vertical profile of mass
                        !!     and internal energy.
  real    :: z_conv     !< (N) A numerical convergence threshold for iterative solving.

  !----------------!
  ! ---Routines--- !
  !----------------!

  public  :: read_multi_vtx
  public  :: intern_init_multi_vtx
  public  :: init_multi_vtx
  public  :: ic_multi_vtx

contains

  !-------------------------------------------------------------------------------
  !>@brief intern_nit_multi_vtx reads an internal namelist to set case settings.
  !>@details The case's namelist is declared and its members are initialized with
  !!         default values. The FV3 input namelist file is read, which overwrites
  !!         default settings. This allows user-friendly modifications to the
  !!         simulation parameters without necessitating recompilation.
  subroutine read_multi_vtx(fname)
    implicit none
    character(*), intent(in) :: fname
    integer :: ierr, f_unit, ios

    namelist /test_case_nml/ &
      num_vortex, aqua, &
      lon, lat, &
      p_0, dp, rsize, &
      q_0, z_q1, z_q2, q_trop, &
      T_0, &
      lapse_rate, &
      z_trop, z_p, z_conv

    ! Default settings
    num_vortex = 0
    aqua       = .true.
    lon(:)     = 0.
    lat(:)     = 0.
    dp(:)      = 1115.
    rsize(:)   = 282000.
    z_q1       = 3000.
    z_q2       = 8000.
    q_0        = 0.021
    q_trop     = 1.e-11
    T_0        = 302.15
    lapse_rate = 7.e-3
    p_0        = 101500.
    z_trop     = 15000.
    z_p        = 7000.
    z_conv     = 1.e-6

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

    if(is_master()) then
      write(*,*) '(multi_vortex,read) num_vortex=',num_vortex
      write(*,*) '(multi_vortex,read) aqua=',aqua
      write(*,*) '(multi_vortex,read) lon=',lon
      write(*,*) '(multi_vortex,read) lat=',lat
      write(*,*) '(multi_vortex,read) dp=',dp
      write(*,*) '(multi_vortex,read) rsize=',rsize
      write(*,*) '(multi_vortex,read) z_q1=',z_q1
      write(*,*) '(multi_vortex,read) z_q2=',z_q2
      write(*,*) '(multi_vortex,read) q_0=',q_0
      write(*,*) '(multi_vortex,read) q_trop=',q_trop
      write(*,*) '(multi_vortex,read) T_0=',T_0
      write(*,*) '(multi_vortex,read) lapse_rate=',lapse_rate
      write(*,*) '(multi_vortex,read) p_0=',p_0
      write(*,*) '(multi_vortex,read) z_trop=',z_trop
      write(*,*) '(multi_vortex,read) z_p=',z_p
      write(*,*) '(multi_vortex,read) z_conv=',z_conv
    endif

  end subroutine read_multi_vtx

  !-------------------------------------------------------------------------------
  !>@brief intern_nit_multi_vtx reads an internal namelist to set case settings.
  !>@details The case's namelist is declared and its members are initialized with
  !!         default values. The FV3 input namelist file is read, which overwrites
  !!         default settings. This allows user-friendly modifications to the
  !!         simulation parameters without necessitating recompilation.
  integer function intern_init_multi_vtx()
    implicit none
    integer :: ios

    namelist /test_case_nml/ &
      num_vortex, aqua, &
      lon, lat, &
      p_0, dp, rsize, &
      q_0, z_q1, z_q2, q_trop, &
      T_0, &
      lapse_rate, &
      z_trop, z_p, z_conv

    ! Default settings
    num_vortex = 0
    aqua       = .true.
    lon(:)     = 0.
    lat(:)     = 0.
    dp(:)      = 1115.
    rsize(:)   = 282000.
    z_q1       = 3000.
    z_q2       = 8000.
    q_0        = 0.021
    q_trop     = 1.e-11
    T_0        = 302.15
    lapse_rate = 7.e-3
    p_0        = 101500.
    z_trop     = 15000.
    z_p        = 7000.
    z_conv     = 1.e-6

    ! Read user-defined values in namelist
    read(input_nml_file,test_case_nml,iostat=ios)
    intern_init_multi_vtx = check_nml_error(ios,'test_case_nml')

    if(is_master()) then
      write(*,*) '(multi_vortex,init) num_vortex=',num_vortex
      write(*,*) '(multi_vortex,init) aqua=',aqua
      write(*,*) '(multi_vortex,init) lon=',lon
      write(*,*) '(multi_vortex,init) lat=',lat
      write(*,*) '(multi_vortex,init) dp=',dp
      write(*,*) '(multi_vortex,init) rsize=',rsize
      write(*,*) '(multi_vortex,init) z_q1=',z_q1
      write(*,*) '(multi_vortex,init) z_q2=',z_q2
      write(*,*) '(multi_vortex,init) q_0=',q_0
      write(*,*) '(multi_vortex,init) q_trop=',q_trop
      write(*,*) '(multi_vortex,init) T_0=',T_0
      write(*,*) '(multi_vortex,init) lapse_rate=',lapse_rate
      write(*,*) '(multi_vortex,init) p_0=',p_0
      write(*,*) '(multi_vortex,init) z_trop=',z_trop
      write(*,*) '(multi_vortex,init) z_p=',z_p
      write(*,*) '(multi_vortex,init) z_conv=',z_conv
      write(*,*) '(multi_vortex) exit intern_init_multi_vtx=',intern_init_multi_vtx
    endif

  end function intern_init_multi_vtx

  !-------------------------------------------------------------------------------
  !>@brief init_multi_vtx reads a namelist to initialize case settings.
  !>@details The case's namelist is declared and its members are initialized with
  !!         default values. The FV3 input namelist file is read, which overwrites
  !!         default settings. This allows user-friendly modifications to the
  !!         simulation parameters without necessitating recompilation.
  integer function init_multi_vtx(f_unit)
    implicit none
    integer, intent(in) :: f_unit
    integer :: ios

    namelist /test_case_nml/ &
      num_vortex, aqua, &
      lon, lat, &
      p_0, dp, rsize, &
      q_0, z_q1, z_q2, q_trop, &
      T_0, &
      lapse_rate, &
      z_trop, z_p, z_conv

    ! Default settings
    num_vortex = 0
    aqua       = .true.
    lon(:)     = 0.
    lat(:)     = 0.
    dp(:)      = 1115.
    rsize(:)   = 282000.
    z_q1       = 3000.
    z_q2       = 8000.
    q_0        = 0.021
    q_trop     = 1.e-11
    T_0        = 302.15
    lapse_rate = 7.e-3
    p_0        = 101500.
    z_trop     = 15000.
    z_p        = 7000.
    z_conv     = 1.e-6

    ! Read user-defined values in namelist (specified values replace default)
    read(f_unit,test_case_nml,iostat=ios)
    init_multi_vtx = check_nml_error(ios,'test_case_nml')

    if(is_master()) then
      write(*,*) '(multi_vortex,init) num_vortex=',num_vortex
      write(*,*) '(multi_vortex,init) aqua=',aqua
      write(*,*) '(multi_vortex,init) lon=',lon
      write(*,*) '(multi_vortex,init) lat=',lat
      write(*,*) '(multi_vortex,init) dp=',dp
      write(*,*) '(multi_vortex,init) rsize=',rsize
      write(*,*) '(multi_vortex,init) z_q1=',z_q1
      write(*,*) '(multi_vortex,init) z_q2=',z_q2
      write(*,*) '(multi_vortex,init) q_0=',q_0
      write(*,*) '(multi_vortex,init) q_trop=',q_trop
      write(*,*) '(multi_vortex,init) T_0=',T_0
      write(*,*) '(multi_vortex,init) lapse_rate=',lapse_rate
      write(*,*) '(multi_vortex,init) p_0=',p_0
      write(*,*) '(multi_vortex,init) z_trop=',z_trop
      write(*,*) '(multi_vortex,init) z_p=',z_p
      write(*,*) '(multi_vortex,init) z_conv=',z_conv
      write(*,*) '(multi_vortex) exit init_multi_vtx=',init_multi_vtx
    endif

  end function init_multi_vtx

  !-------------------------------------------------------------------------------
  !>@brief ic_multi_vtx initializes the model state with idealized warm-core cyclones.
  !>@details The model's initial state is set in this routine, using the parameters
  !!         specified in the user's namelist (default settings if none specified).
  subroutine ic_multi_vtx(u, v, w, pt, delp, delz, q, &
                          pk, peln, pe, pkz, phis, ps, &
                          ak, bk, gridstruct, domain)

    implicit none
    real, intent(IN), dimension(npz+1) :: ak, bk
    type(fv_grid_type), target, intent(IN) :: gridstruct
    type(domain2d), intent(IN) :: domain

    real, intent(OUT), dimension(isd:ied,jsd:jed+1,npz) :: u  !< abscissa-sided
    real, intent(OUT), dimension(isd:ied+1,jsd:jed,npz) :: v  !< ordinate-sided
    real, intent(OUT), dimension(isd:ied,jsd:jed,npz)   :: w, pt, delp !, delz  !< cell-centered
    real, intent(INOUT), dimension(is:,js:,1:) :: delz
    !real, intent(INOUT), dimension(is:ie,js:je,1:npz) :: delz
    real, intent(INOUT), dimension(isd:ied,jsd:jed,npz,ncnst) :: q !< cell-centered, dim-4 for species

    real, intent(OUT), dimension(is:ie,js:je,npz+1) :: pk
    real, intent(OUT), dimension(is:ie,npz+1,js:je) :: peln
    real, intent(OUT), dimension(is-1:ie+1,npz+1,js-1:je+1) :: pe
    real, intent(OUT), dimension(is:ie,js:je,npz) :: pkz
    real, intent(OUT), dimension(isd:ied,jsd:jed) :: phis, ps

    real(kind=R_GRID), pointer, dimension(:,:,:) :: grid, agrid
    real(kind=R_GRID), pointer, dimension(:,:)   :: area

    real(kind=R_GRID), parameter :: earth_r = const_A !< Earth's radius
    real, parameter :: rdgrav = rdgas/grav  !< R_d / g
    real, parameter :: rrdgrav = grav/rdgas !< g / R_d (reciprocal rdgrav)

    real :: Tv0 !< Surface virtual temperature (K)
    real :: Tvt !< Tropopause virtual temperature (K)
    real :: ptt !< Tropopause pressure (Pa)

    real, dimension(isd:ied,jsd:jed,npz+1) :: gz      !< geopotential height (m)

    real, dimension(num_vortex,isd:ied,jsd:jed) :: rc     !< distance from the vortex center
    real :: rc_total                                  !< sum of all rc from all vortices
    real(kind=R_GRID), dimension(2) :: vortcenter     !< working vortex coordinates
    real, dimension(num_vortex,isd:ied,jsd:jed) :: num_weight !< 1 / (r_n)^p
    real, dimension(isd:ied,jsd:jed)            :: den_weight !< SUM[1 / (r_n)^p]
    real, dimension(num_vortex,isd:ied,jsd:jed) :: weight !< num_weight / den_weight

    integer :: i,j,k,n
    integer :: iter
    integer :: sphum !< specific humidity tracer index

    real :: p, pl, z, z0        !< working variables for pressure, log-p, gz
    real :: p_lower, p_upper
    real :: ziter, piter, titer !< iterative terms
    real :: uu, vv

    real(kind=R_GRID), dimension(2) :: pa
    real(kind=R_GRID), dimension(3) :: e1, e2, ex, ey

    real, dimension(num_vortex,isd:ied,jsd:jed+1) :: num_weight_u
    real, dimension(isd:ied,jsd:jed+1)            :: den_weight_u
    real, dimension(num_vortex,isd:ied,jsd:jed+1) :: rc_u, weight_u
    real, dimension(isd:ied,jsd:jed+1) :: gz_u, p_u, peln_u, ps_u, u1, u2
    real(kind=R_GRID), dimension(isd:ied,jsd:jed+1) :: lat_u, lon_u

    real, dimension(num_vortex,isd:ied+1,jsd:jed) :: num_weight_v
    real, dimension(isd:ied+1,jsd:jed)            :: den_weight_v
    real, dimension(num_vortex,isd:ied+1,jsd:jed) :: rc_v, weight_v
    real, dimension(isd:ied+1,jsd:jed) :: gz_v, p_v, peln_v, ps_v, v1, v2
    real(kind=R_GRID), dimension(isd:ied+1,jsd:jed) :: lat_v, lon_v

    ! These members are for debugging
    logical, parameter :: calc_rc = .true.
    logical, parameter :: calc_ps = .true.
    logical, parameter :: calc_delp = .true.
    logical, parameter :: calc_misc_p = .true.
    logical, parameter :: calc_height = .true.
    logical, parameter :: calc_t = .true.
    logical, parameter :: calc_uv = .true.
    logical, parameter :: calc_q = .true.
    logical, parameter :: calc_nonhydro = .true.
    real :: max_real, min_real
    logical :: has_warned

    has_warned = .false.

    max_real = SQRT(HUGE(0.0) * 0.1)
    min_real = SQRT(TINY(0.0) * 10.)

    if (num_vortex < 0) then
      call mpp_error(NOTE, '[ERROR] (multi_vortex.F90, test_multi_vtx):')
      call mpp_error(NOTE, '    Non-physical number of vortices specified (num_vortex < 0).')
      ! return/stop
      call mpp_error(FATAL, 'FATAL, test_multi_vtx: Invalid setting of num_vortex.')
    else if (num_vortex .EQ. 0) then
      call mpp_error(NOTE, '[WARN] (multi_vortex.F90, test_multi_vtx):')
      call mpp_error(NOTE, '    Number of vortices specified is ZERO (num_vortex == 0).')
      ! continue if allowed
      call mpp_error(WARNING, 'WARNING, test_multi_vtx: num_vortex is set to ZERO.')
    endif

    grid  => gridstruct%grid_64
    agrid => gridstruct%agrid_64
    area  => gridstruct%area_64

    ! Compute ps, phis, delp, aux pressure variables, Temperature, winds
    ! (with or without perturbation), moisture, w, delz

    ! Compute p, z, T on both the staggered and unstaggered grids. Then compute the zonal
    !  and meridional winds on both grids, and rotate as needed

    ! Use hydrostatic assumption, ideal gas law:
    Tv0 = T_0*(1.+0.608*q_0)
    Tvt = Tv0 - lapse_rate*z_trop
    ptt = p_0*(Tvt/Tv0)**(grav/rdgas/lapse_rate)

    den_weight   = 0.0
    num_weight   = 0.0
    den_weight_u = 0.0
    num_weight_u = 0.0
    den_weight_v = 0.0
    num_weight_v = 0.0

    if (aqua) then
      ! Initialize dry hydrostatic atmosphere at rest
      q(:,:,:,:) = 3.e-6
      u(:,:,:)   = 0.0
      v(:,:,:)   = 0.0
      if (.not. hydrostatic) w(:,:,:) = 0.0
      phis = 0.0
      call hydro_eq(npz, is, ie, js, je, ps, phis, p_0, &
                delp, ak, bk, pt, delz, area, ng, .false., hydrostatic, hybrid_z, domain)
    endif

    if (calc_rc) then
      ! Calculate r and distance-weighting to use if num_vortex > 1
      do j=js2,je2
        do i=is2,ie2
          do n=1,num_vortex
            vortcenter = (/ lon(n) * pi / 180.0, lat(n) * pi / 180.0 /)
            rc(n,i,j) = great_circle_dist(agrid(i,j,:), vortcenter, earth_r)
          enddo
          if (num_vortex > 1) then
            do n=1,num_vortex
              num_weight(n,i,j) = 1. / (rc(n,i,j)**2)
              den_weight(i,j)   = den_weight(i,j) + 1. / (rc(n,i,j)**2)
            enddo
            do n=1,num_vortex
              weight(n,i,j)     = num_weight(n,i,j) / den_weight(i,j)
            enddo
          endif
        enddo
      enddo
    endif

    if (calc_ps) then
      ! Calculate surface pressure
      do j=js2,je2
        do i=is2,ie2
          do n=1,num_vortex
            if (num_vortex > 1) then
              if (n .EQ. 1) then
                ps(i,j) = weight(n,i,j) * (p_0 - dp(n)*EXP( -SQRT((rc(n,i,j)/rsize(n))**3) ))
              else
                ps(i,j) = ps(i,j) + (weight(n,i,j) * (p_0 - dp(n)*EXP( -SQRT((rc(n,i,j)/rsize(n))**3) )))
              endif
            else
              ps(i,j) = p_0 - dp(n)*EXP( -SQRT((rc(n,i,j)/rsize(n))**3) )
            endif
          enddo
        enddo
      enddo
    endif

    if (calc_delp) then
      ! Calculate layer thicknesses in pressure
      do k=1,npz
        do j=js2,je2
          do i=is2,ie2
            delp(i,j,k) = ak(k+1)-ak(k) + ps(i,j)*(bk(k+1)-bk(k))
          enddo
        enddo
      enddo
    endif

    if (calc_misc_p) then
      ! Calculate misc/aux pressure variables
      ! Calcs for model top first
      do j=js,je
        do i=is,ie
          pe(i,1,j) = ptop
        enddo
        do i=is,ie
          peln(i,1,j) = LOG(ptop)
          pk(i,j,1) = ptop**kappa
        enddo
        ! Calcs for everywhere else in z
        do k=2,npz+1
          do i=is,ie
            pe(i,k,j) = ak(k) + ps(i,j)*bk(k)
          enddo
          do i=is,ie
            ! pk(i,j,k) = EXP(kappa*LOG(pe(i,k,j)))
            pk(i,j,k)   = pe(i,k,j)**kappa
            peln(i,k,j) = LOG(pe(i,k,j))
          enddo
        enddo
      enddo
      do k=1,npz
        do j=js,je
          do i=is,ie
            pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(kappa*(peln(i,k+1,j)-peln(i,k,j)))
          enddo
        enddo
      enddo
    endif

    if (calc_height) then
      ! Height: Use Newton's method
      ! Cell centered
      do j=js2,je2
        do i=is2,ie2
          phis(i,j) = 0.
          gz(i,j,npz+1) = 0.
        enddo
      enddo
      do k=npz,1,-1
        do j=js2,je2
          do i=is2,ie2
            has_warned = .false.
            do n=1,num_vortex
              ! p = pe(i,k,j)
              if (k .EQ. 1) then
                p = ptop
              else
                p = ak(k) + ps(i,j)*bk(k)
              endif
              z = gz(i,j,k+1)
              do iter=1,30
                ziter = z
                piter = TC_pressure(ziter, rc(n,i,j), dp(n), rsize(n))
                titer = TC_temperature(ziter, rc(n,i,j), dp(n), rsize(n))
                z = ziter + (piter - p)*rdgrav*titer/piter

                if (ABS(z - ziter) < z_conv) exit
              enddo
              if (num_vortex > 1) then
                if (n .EQ. 1) then
                  gz(i,j,k) = weight(n,i,j) * z
                else
                  gz(i,j,k) = gz(i,j,k) + (weight(n,i,j) * z)
                endif
              else
                gz(i,j,k) = z
              endif
            enddo
          enddo
        enddo
      enddo
    endif

    if (calc_t) then
      ! Potential Temperature: Compute from hydro balance
      do k=1,npz
        do j=js2,je2
          do i=is2,ie2
!            pt(i,j,k) = rrdgrav * ( gz(i,j,k) - gz(i,j,k+1) ) / ( peln(i,k+1,j) - peln(i,k,j))

            if (k .EQ. 1) then
              p_upper = LOG(ptop)
            else
              p_upper = LOG(ak(k) + ps(i,j)*bk(k))
            endif
            p_lower = LOG(ak(k+1) + ps(i,j)*bk(k+1))

            pt(i,j,k) = rrdgrav * ( gz(i,j,k) - gz(i,j,k+1) ) / (p_lower - p_upper)
            if ((i.EQ.1) .AND. (j.EQ.1)) then
              write(*,*) '[KA] (calc_t) gz=',gz(i,j,k),'k=',k,'T=',pt(i,j,k)
            endif
          enddo
        enddo
      enddo
    endif

    if (calc_uv) then
      ! Compute height and temperature for grids with u and v points also,
      ! to be able to compute the local winds...
      ! Use temporary 2d arrays for this purpose
      do j=js2,je2+1
        do i=is2,ie2
          call mid_pt_sphere(grid(i,j,:),grid(i+1,j,:),pa)
          lat_u(i,j) = pa(2)
          lon_u(i,j) = pa(1)
          call get_unit_vect2(grid(i,j,:),grid(i+1,j,:),e1)
          call get_latlon_vector(pa,ex,ey)
          u1(i,j) = inner_prod(e1,ex) !u components
          u2(i,j) = inner_prod(e1,ey)
          do n=1,num_vortex
            vortcenter = (/ lon(n) * pi / 180.0, lat(n) * pi / 180.0 /)
            rc_u(n,i,j) = great_circle_dist(pa, vortcenter, earth_r)
          enddo
          if (num_vortex > 1) then
            do n=1,num_vortex
              num_weight_u(n,i,j) = 1. / (rc_u(n,i,j)**2)
              den_weight_u(i,j)   = den_weight_u(i,j) + 1. / (rc_u(n,i,j)**2)
            enddo
            do n=1,num_vortex
              weight_u(n,i,j)     = num_weight_u(n,i,j) / den_weight_u(i,j)
            enddo
          endif
        enddo
      enddo

      do n=1,num_vortex
        vortcenter = (/ lon(n) * pi / 180.0, lat(n) * pi / 180.0 /)
        do j=js2,je2+1
          do i=is2,ie2
            p_u(i,j) = p_0 - dp(n) * EXP( -SQRT((rc_u(n,i,j)/rsize(n))**3) )
            peln_u(i,j) = LOG(p_u(i,j))
            ps_u(i,j) = p_u(i,j)
            gz_u(i,j) = 0.
            do k=npz,1,-1
              ! Pressure (Top of interface)
              p = ak(k) + ps_u(i,j)*bk(k)
              pl = LOG(p)
              ! Height (top of interface); use newton's method
              z = gz_u(i,j) !first guess, height of lower level
              z0 = z
              do iter=1,30
                ziter = z
                piter = TC_pressure(ziter,rc_u(n,i,j), dp(n), rsize(n))
                titer = TC_temperature(ziter,rc_u(n,i,j), dp(n), rsize(n))
                z = ziter + (piter - p)*rdgrav*titer/piter
                if (ABS(z - ziter) < z_conv) exit
              enddo
              ! Now compute winds
              call TC_uwind_pert(0.5*(z+z0),rc_u(n,i,j),dp(n),rsize(n),vortcenter(1),vortcenter(2),lon_u(i,j),lat_u(i,j), uu, vv)

              if (num_vortex > 1) then
                u(i,j,k) = u(i,j,k) + (weight_u(n,i,j) * (u1(i,j)*uu + u2(i,j)*vv))
              else
                u(i,j,k) = u1(i,j)*uu + u2(i,j)*vv
              endif

              gz_u(i,j) = z
              p_u(i,j) = p
              peln_u(i,j) = pl
            enddo
          enddo
        enddo
      enddo

      do j=js2,je2
        do i=is2,ie2+1
          call mid_pt_sphere(grid(i,j,:),grid(i,j+1,:),pa)
          lat_v(i,j) = pa(2)
          lon_v(i,j) = pa(1)
          call get_unit_vect2(grid(i,j,:),grid(i,j+1,:),e2)
          call get_latlon_vector(pa,ex,ey)
          v1(i,j) = inner_prod(e2,ex) !v components
          v2(i,j) = inner_prod(e2,ey)
          do n=1,num_vortex
            vortcenter = (/ lon(n) * pi / 180.0, lat(n) * pi / 180.0 /)
            rc_v(n,i,j) = great_circle_dist(pa, vortcenter, earth_r)
          enddo
          if (num_vortex > 1) then
            do n=1,num_vortex
              num_weight_v(n,i,j) = 1. / (rc_v(n,i,j)**2)
              den_weight_v(i,j)   = den_weight_v(i,j) + 1. / (rc_v(n,i,j)**2)
            enddo
            do n=1,num_vortex
              weight_v(n,i,j)     = num_weight_v(n,i,j) / den_weight_v(i,j)
            enddo
          endif
        enddo
      enddo

      do n=1,num_vortex
        vortcenter = (/ lon(n) * pi / 180.0, lat(n) * pi / 180.0 /)
        do j=js2,je2
          do i=is2,ie2+1
            p_v(i,j) = p_0 - dp(n)*EXP( -SQRT((rc_v(n,i,j)/rsize(n))**3) )
            peln_v(i,j) = LOG(p_v(i,j))
            ps_v(i,j) = p_v(i,j)
            gz_v(i,j) = 0.
            do k=npz,1,-1
              ! Pressure (Top of interface)
              p = ak(k) + ps_v(i,j)*bk(k)
              pl = LOG(p)
              ! Height (top of interface); use newton's method
              z = gz_v(i,j) !first guess, height of lower level
              z0 = z
              do iter=1,30
                ziter = z
                piter = TC_pressure(ziter,rc_v(n,i,j), dp(n), rsize(n))
                titer = TC_temperature(ziter,rc_v(n,i,j), dp(n), rsize(n))
                z = ziter + (piter - p)*rdgrav*titer/piter
                if (ABS(z - ziter) < z_conv) exit
              enddo
              ! Now compute winds
              call TC_uwind_pert(0.5*(z+z0),rc_v(n,i,j),dp(n),rsize(n),vortcenter(1),vortcenter(2),lon_v(i,j),lat_v(i,j), uu, vv)

              if (num_vortex > 1) then
                v(i,j,k) = v(i,j,k) + (weight_v(n,i,j) * (v1(i,j)*uu + v2(i,j)*vv))
              else
                v(i,j,k) = v1(i,j)*uu + v2(i,j)*vv
              endif
              gz_v(i,j) = z
              p_v(i,j) = p
              peln_v(i,j) = pl
            enddo
          enddo
        enddo
      enddo
    endif !< calc_uv

    if (calc_q) then
      ! Compute moisture and other tracer fields, as desired
      do n=1,ncnst
        do k=1,npz
          do j=js2,je2
            do i=is2,ie2
              q(i,j,k,n) = 0.
            enddo
          enddo
        enddo
      enddo
      if (.not. adiabatic) then
        sphum = get_tracer_index(MODEL_ATMOS, 'sphum')
        do k=1,npz
          do j=js2,je2
            do i=is2,ie2
              z = 0.5*(gz(i,j,k) + gz(i,j,k+1))
              q(i,j,k,sphum) = TC_sphum(z, z_q1, z_q2)
            enddo
          enddo
        enddo
      endif
    endif

    if (is_master()) then
      write(*,*) '[KA] (delz) is=',is
      write(*,*) '[KA] (delz) ie=',ie
      write(*,*) '[KA] (delz) js=',js
      write(*,*) '[KA] (delz) je=',je
      write(*,*) '[KA] (delz) is2=',is2
      write(*,*) '[KA] (delz) ie2=',ie2
      write(*,*) '[KA] (delz) js2=',js2
      write(*,*) '[KA] (delz) je2=',je2
      write(*,*) '[KA] (delz) shape=',SHAPE(delz)
      write(*,*) '[KA] (delz) size=',SIZE(delz)
!      write(*,*) '[KA] (delz) delz=',delz
    endif

    if (calc_nonhydro) then
      ! Compute nonhydrostatic variables, if needed
      if (.not. hydrostatic) then
        do k=1,npz
          do j=js2,je2
            do i=is2,ie2
              w(i,j,k) = 0.
            enddo
          enddo
          do j=js,je
            do i=is,ie
              if (is_master()) then
                p_upper = gz(i,j,k) - gz(i,j,k+1)
                write(*,*) '[KA] i=',i,'j=',j,'k=',k,'dz=',delz(i,j,k),'clc=',p_upper
              endif
!              delz(i,j,k) = gz(i,j,k) - gz(i,j,k+1)
              delz(i,j,k) = gz(i,j,k+1) - gz(i,j,k)
            enddo
          enddo
        enddo
      endif
    endif

    if (is_master()) then
      write(*,*) '[KA] (delz) is=',is
      write(*,*) '[KA] (delz) ie=',ie
      write(*,*) '[KA] (delz) js=',js
      write(*,*) '[KA] (delz) je=',je
      write(*,*) '[KA] (delz) is2=',is2
      write(*,*) '[KA] (delz) ie2=',ie2
      write(*,*) '[KA] (delz) js2=',js2
      write(*,*) '[KA] (delz) je2=',je2
      write(*,*) '[KA] (delz) shape=',SHAPE(delz)
      write(*,*) '[KA] (delz) size=',SIZE(delz)
!      write(*,*) '[KA] (delz) delz=',delz
    endif
!    delz = 0.0

    nullify(grid)
    nullify(agrid)
    nullify(area)

  contains

    !-------------------------------------------------------------------------------
    !>@brief TC_temperature solves for virtual temperature at radius r and height z.
    real function TC_temperature(z, r, dp, rsize)

      real, intent(IN) :: z, r, dp, rsize
      real :: Tv, term1, term2
      real :: exp_term
      real :: max_real, ln_max_real, max_exp

      if (z > z_trop) then
        TC_temperature = Tvt
        return
      endif

      max_real    = HUGE(0.0)           !< maximum value of real types
      ln_max_real = LOG(max_real)       !< natural-log of the max value
      max_exp     = 0.5 * ln_max_real   !< max value allowable for exp_term

      Tv = Tv0 - lapse_rate*z
      exp_term = SQRT(r/rsize)**3 + (z/z_p)**2

      if (exp_term < max_exp) then !< this is to prevent overflow
        term1 = (grav*z_p*z_p) * ( 1. - (p_0/dp) * EXP( exp_term ) )
        term2 = 2*rdgas*Tv*z
        TC_temperature = Tv + (Tv * ( 1./(1 + term2/term1) - 1.))
      else
        TC_temperature = Tv
      endif

    end function TC_temperature

    !-------------------------------------------------------------------------------
    !>@brief TC_pressure solves for pressure at radius r and height z.
    real function TC_pressure(z, r, dp, rsize)

      real, intent(IN) :: z, r, dp, rsize

      if (z <= z_trop) then
        TC_pressure = p_0*EXP(grav/(rdgas*lapse_rate) * LOG( (Tv0-lapse_rate*z)/Tv0) ) &
        - dp*EXP(-SQRT((r/rsize)**3) - (z/z_p)**2) &
        * EXP( grav/(rdgas*lapse_rate) * LOG( (Tv0-lapse_rate*z)/Tv0) )
      else
        TC_pressure = ptt*EXP(grav*(z_trop-z)/(rdgas*Tvt))
      endif

    end function TC_pressure

    !-------------------------------------------------------------------------------
    !>@brief TC_uwind_pert calculates horizontal wind at radius r and height z.
    subroutine TC_uwind_pert(z,r,dp,rsize,c_lon,c_lat,lon,lat,uu,vv)

      real, intent(IN) :: z, r, dp, rsize
      real(kind=R_GRID), intent(IN) :: c_lon, c_lat, lon, lat
      real, intent(OUT) :: uu, vv
      real :: rfac, Tvrd, vt, fr5, d1, d2, d
      real :: fc
      real :: exp_term
      real :: max_real, ln_max_real, max_exp

      if (z > z_trop) then
        uu = 0.
        vv = 0.
        return
      endif

      max_real    = HUGE(0.0)           !< maximum value of real types
      ln_max_real = LOG(max_real)       !< natural-log of the max value
      max_exp     = 0.5 * ln_max_real   !< max value allowable for exp_term

      fc = 2.*omega*sin(c_lat)
      rfac = SQRT(r/rsize)**3

      fr5 = 0.5*fc*r
      Tvrd = (Tv0 - lapse_rate*z)*rdgas

      exp_term = rfac + (z/z_p)**2

      if (exp_term < max_exp) then
        vt = -fr5 + SIGN(1.0,c_lat) * SQRT( fr5**2 - (1.5 * rfac * Tvrd) / &
             ( 1. + 2.*Tvrd*z/(grav*z_p**2) - (p_0/dp)*EXP(exp_term) ) )
      else
        vt = 0.
      endif

      d1 = SIN(c_lat)*COS(lat) - COS(c_lat)*SIN(lat)*COS(lon - c_lon)
      d2 = COS(c_lat)*SIN(lon - c_lon)
      d = MAX(1.e-25,SQRT(d1*d1 + d2*d2))

      uu = vt * d1/d
      vv = vt * d2/d

    end subroutine TC_uwind_pert

    !-------------------------------------------------------------------------------
    !>@brief TC_sphum solves for specific humidity at height z (no structure in r).
    real function TC_sphum(z, z_q1, z_q2)
      real, intent(IN) :: z, z_q1, z_q2

      TC_sphum = q_trop
      if (z < z_trop) then
        TC_sphum = q_0 * EXP(-z/z_q1) * EXP(-(z/z_q2)**2)
      endif

    end function TC_sphum

  end subroutine ic_multi_vtx

end module multi_vtx_mod
