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

module test_control_mod
  !-------------------------------------------------------------------------------
  !>@brief The module 'test_control_mod' contains information
  !!       essential or shared amongst most/all idealized cases.
  !>@details This module holds several variables for use in
  !!         most, if not all, idealized cases. These variables include
  !!         flags indicating certain assumptions in the initial state
  !!         (such as whether to assume hydrostatic balance), and
  !!         values for the basic structure of the model grid.
  ! <table>
  !   <tr>
  !     <th>Module Name</th>
  !     <th>Functions Included</th>
  !   </tr>
  !   <tr>
  !     <td>fv_arrays_mod</td>
  !     <td>fv_grid_type, fv_flags_type, fv_grid_bounds_type</td>
  !   </tr>
  ! </table>

  use fv_arrays_mod, only: fv_grid_type, fv_flags_type, fv_grid_bounds_type

  implicit none
  integer :: selected_case = 0
  real :: alpha = 0.0

  !< Enumerated list of test/idealized cases
  integer, parameter :: CASE_DIV_CONS  = 1   !< Divergence conservation test
  integer, parameter :: CASE_NL_DEFORM = 2   !< Ideal non-linear deformation flow
  integer, parameter :: CASE_COS_BELL  = 3   !< Cosine bell advection
  integer, parameter :: CASE_ZONAL     = 4   !< Zonal geostrophic-balance flow
  integer, parameter :: CASE_ZONAL_MTN = 5   !< Zonal geostrophic-balance flow over a mountain
  integer, parameter :: CASE_NR_POT    = 6   !< Non-rotating potential flow
  integer, parameter :: CASE_ROSSBY_4  = 7   !< Rossby wavenumber-4
  integer, parameter :: CASE_BT_INSTAB = 8   !< Barotropic instability
  integer, parameter :: CASE_SOLITON   = 9   !< Soliton propagation
  integer, parameter :: CASE_POLAR_VTX = 10  !< Polar vortex
  integer, parameter :: CASE_IDEAL_MTN = 11  !< Hydrostatic balanced state with a mountain
  ! etc. ...
  integer, parameter :: CASE_MULTI_VTX = 99


  !< Note: no idea what "ng" is
  public  :: set_test_hgrid
  integer :: is, ie, js, je, isd, ied, jsd, jed, npx, npy, ng
  logical :: bounded
  integer :: is2, ie2, js2, je2

  public  :: set_test_vgrid
  integer :: npz
  real    :: ptop
  logical :: hybrid_z

  public  :: set_test_tracers
  integer :: ncnst, nwat

  public  :: set_test_phys
  logical :: hydrostatic, adiabatic, mountain, moist_phys
  real    :: dry_mass

contains

  !> The following set-and-get system of functions is made to prevent circular dependencies.
  !> The variables themselves are stored to cut down on excessively large argument lists
  !> and thus improve clarity.
  !-------------------------------------------------------------------------------
  !>@brief set_test_hgrid stores information for the working horizontal grid.
  !>@details This routine stores information from this PE's working
  !!         horizontal grid, which is then used by ideal case subroutines.
  subroutine set_test_hgrid(gridstruct, bounds, num_x, num_y, n_g)
    implicit none
    type(fv_grid_type), intent(in) :: gridstruct
    type(fv_grid_bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_x, num_y
    integer, intent(in) :: n_g
    bounded = gridstruct%bounded_domain
    is  = bounds%is
    ie  = bounds%ie
    js  = bounds%js
    je  = bounds%je
    isd = bounds%isd
    ied = bounds%ied
    jsd = bounds%jsd
    jed = bounds%jed
    npx = num_x
    npy = num_y
    ng  = n_g !< what is this variable?

    if (bounded) then
      is2 = isd
      ie2 = ied
      js2 = jsd
      je2 = jed
    else
      is2 = is
      ie2 = ie
      js2 = js
      je2 = je
    endif

  end subroutine set_test_hgrid

  !-------------------------------------------------------------------------------
  !>@brief set_test_vgrid stores information for the working vertical grid.
  !>@details This routine stores information from this PE's working
  !!         vertical grid, which is then used by ideal case subroutines.
  subroutine set_test_vgrid(vert_levs, press_top, is_hybrid)
    implicit none
    integer, intent(in) :: vert_levs
    real, intent(in)    :: press_top
    logical, intent(in) :: is_hybrid
    npz       = vert_levs
    ptop      = press_top
    hybrid_z  = is_hybrid
  end subroutine set_test_vgrid

  !-------------------------------------------------------------------------------
  !>@brief set_test_tracers stores information about atmospheric tracers.
  !>@details This routine stores information about atmospheric tracers,
  !!         which is then used by ideal case subroutines.
  subroutine set_test_tracers(n_tracers, n_water_species)
    implicit none
    integer, intent(in) :: n_tracers, n_water_species
    ncnst = n_tracers
    nwat  = n_water_species
  end subroutine set_test_tracers

  !-------------------------------------------------------------------------------
  !>@brief set_test_phys stores information about basic model physics.
  !>@details This routine stores information about basic model physics,
  !!         such as whether the model is hydrostatic or uses terrain,
  !!         which is then used by ideal case subroutines.
  subroutine set_test_phys(flagstruct)
    implicit none
    type(fv_flags_type), intent(in) :: flagstruct
    hydrostatic = flagstruct%hydrostatic
    adiabatic   = flagstruct%adiabatic
    mountain    = flagstruct%mountain
    moist_phys  = flagstruct%moist_phys
    dry_mass    = flagstruct%dry_mass
  end subroutine set_test_phys

end module test_control_mod
