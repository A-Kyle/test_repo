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

module test_init_mod
  !-------------------------------------------------------------------------------
  !>@brief The module 'test_init_mod' is used to set-up and handle the
  !!       initialization/forcings of idealized tests.
  !>@details This module calls the appropriate routines for a given test case,
  !!         which should often come from an external module specific to the given case.
  !!         The routines include those for loading namelists/options for a case,
  !!         initializing the model state, and adding forcing over time.
  !!         Initialization methods shared amongst ALL cases should also occur here.
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
  !     <td>omega, grav, kappa</td>
  !   </tr>
  !   <tr>
  !     <td>fv_arrays_mod</td>
  !     <td>fv_grid_type, R_GRID</td>
  !   </tr>
  !   <tr>
  !     <td>mpp_domains_mod</td>
  !     <td>mpp_update_domains, domain2d</td>
  !   </tr>
  !   <tr>
  !     <td>mpp_mod</td>
  !     <td>input_nml_file</td>
  !   </tr>
  !   <tr>
  !     <td>fv_mp_mod</td>
  !     <td>is_master, fill_corners, XDir, YDir</td>
  !   </tr>
  !   <tr>
  !     <td>fv_grid_utils_mod</td>
  !     <td>g_sum, ptop_min</td>
  !   </tr>
  !   <tr>
  !     <td>init_hydro_mod</td>
  !     <td>p_var</td>
  !   </tr>
  !   <tr>
  !     <td>mpp_parameter_mod</td>
  !     <td>DGRID_NE</td>
  !   </tr>
  !   <tr>
  !     <td>multi_vtx_mod</td>
  !     <td>intern_init_multi_vtx, init_multi_vtx, ic_multi_vtx</td>
  !   </tr>
  ! </table>

  use test_control_mod
  use constants_mod,      only: omega, grav, kappa
  use fv_arrays_mod,      only: fv_grid_type, R_GRID
  use mpp_domains_mod,    only: mpp_update_domains, domain2d
  use mpp_mod,            only: input_nml_file
  use fv_mp_mod,          only: is_master, fill_corners, XDir, YDir
  use fv_grid_utils_mod,  only: g_sum, ptop_min
  use init_hydro_mod,     only: p_var
  use mpp_parameter_mod,  only: DGRID_NE

  use ideal_mtn_mod,      only: read_ideal_mtn, intern_init_ideal_mtn, init_ideal_mtn, ic_ideal_mtn
  use multi_vtx_mod,      only: read_multi_vtx, intern_init_multi_vtx, init_multi_vtx, ic_multi_vtx

  implicit none
  private
  public :: read_test_nml, intern_init_test, init_test, ic_test

contains

  !-------------------------------------------------------------------------------
  !>@brief intern_init_test directs a case's module to read an internal namelist.
  subroutine read_test_nml(fname)
    implicit none
    character(*), intent(in) :: fname
    ! integer :: ierr
    ! rewind(f_unit)

    select case(selected_case)
      ! case(case_list%DIV_CONS)
      ! case(case_list%NL_DEFORM)
      ! case(case_list%COS_BELL)
      ! case(case_list%ZONAL)
      ! case(case_list%ZONAL_MTN)
      ! case(case_list%NR_POT)
      ! case(case_list%ROSSBY_4)
      ! case(case_list%BT_INSTAB)
      ! case(case_list%SOLITON)
      ! case(case_list%POLAR_VTX)
      case(CASE_IDEAL_MTN)
        call read_ideal_mtn(fname)
      case(CASE_MULTI_VTX)
        call read_multi_vtx(fname)
    end select

  end subroutine read_test_nml

  !-------------------------------------------------------------------------------
  !>@brief intern_init_test directs a case's module to read an internal namelist.
  subroutine intern_init_test()
    implicit none
    integer :: ierr
    ! rewind(f_unit)

    select case(selected_case)
      ! case(case_list%DIV_CONS)
      ! case(case_list%NL_DEFORM)
      ! case(case_list%COS_BELL)
      ! case(case_list%ZONAL)
      ! case(case_list%ZONAL_MTN)
      ! case(case_list%NR_POT)
      ! case(case_list%ROSSBY_4)
      ! case(case_list%BT_INSTAB)
      ! case(case_list%SOLITON)
      ! case(case_list%POLAR_VTX)
      case(CASE_IDEAL_MTN)
        ierr = intern_init_ideal_mtn()
      case(CASE_MULTI_VTX)
        ierr = intern_init_multi_vtx()
    end select

  end subroutine intern_init_test

  !-------------------------------------------------------------------------------
  !>@brief init_test directs a case's module to read an external namelist.
  subroutine init_test(f_unit)
    implicit none
    integer, intent(in) :: f_unit
    integer :: ierr

    rewind(f_unit)

    select case(selected_case)
      ! case(case_list%DIV_CONS)
      ! case(case_list%NL_DEFORM)
      ! case(case_list%COS_BELL)
      ! case(case_list%ZONAL)
      ! case(case_list%ZONAL_MTN)
      ! case(case_list%NR_POT)
      ! case(case_list%ROSSBY_4)
      ! case(case_list%BT_INSTAB)
      ! case(case_list%SOLITON)
      ! case(case_list%POLAR_VTX)
      case(CASE_IDEAL_MTN)
        ierr = intern_init_ideal_mtn()
      case(CASE_MULTI_VTX)
        ierr = init_multi_vtx(f_unit)
    end select

  end subroutine init_test

  !-------------------------------------------------------------------------------
  !>@brief ic_test uses a case's module to initialize the model state.
  !>@details The ic_test routine is used to set the initial state of the model.
  !!         This includes initializing some basic fields, and then calling
  !!         routines specific to the chosen test case. Some updates and
  !!         data logging are done after the case-specific initialization occurs.
  subroutine ic_test(u, v, w, pt, delp, delz, q, &
                          pk, peln, pe, pkz, phis, ps, &
                          uc, vc, ua, va, ze0, ak, bk, gridstruct, domain)
    implicit none
    real, intent(INOUT) ::    u(isd:ied  ,jsd:jed+1,npz)
    real, intent(INOUT) ::    v(isd:ied+1,jsd:jed  ,npz)
    real, intent(INOUT) ::    w(isd:  ,jsd:  ,1:)
    real, intent(INOUT) ::   pt(isd:ied  ,jsd:jed  ,npz)
    real, intent(INOUT) :: delp(isd:ied  ,jsd:jed  ,npz)
!    real, intent(INOUT) :: delz(isd:,jsd:,1:)
    real, intent(INOUT) :: delz(is:,js:,1:)
    real, intent(INOUT) ::    q(isd:ied  ,jsd:jed  ,npz, ncnst)

    real, intent(INOUT) ::   pk(is:ie    ,js:je    ,npz+1)
    real, intent(INOUT) :: peln(is :ie   ,npz+1    ,js:je)
    real, intent(INOUT) ::   pe(is-1:ie+1,npz+1,js-1:je+1)
    real, intent(INOUT) ::  pkz(is:ie    ,js:je    ,npz)
    real, intent(INOUT) :: phis(isd:ied  ,jsd:jed)
    real, intent(INOUT) ::   ps(isd:ied  ,jsd:jed)

    real, intent(INOUT) ::   uc(isd:ied+1,jsd:jed  ,npz)
    real, intent(INOUT) ::   vc(isd:ied  ,jsd:jed+1,npz)
    real, intent(INOUT) ::   ua(isd:ied  ,jsd:jed  ,npz)
    real, intent(INOUT) ::   va(isd:ied  ,jsd:jed  ,npz)
    real, intent(INOUT) ::  ze0(is:,js:,1:)

    real, intent(INOUT) ::   ak(npz+1)
    real, intent(INOUT) ::   bk(npz+1)

    type(fv_grid_type), target, intent(INOUT) :: gridstruct
    type(domain2d), intent(INOUT) :: domain

    real(kind=R_GRID), pointer, dimension(:,:,:) :: agrid, grid
    real(kind=R_GRID), pointer, dimension(:,:)   :: area
    real, pointer, dimension(:,:) :: fC, f0
    logical, pointer :: cubed_sphere

    integer :: i, j
    real :: ftop

    grid          => gridstruct%grid_64
    agrid         => gridstruct%agrid_64
    area          => gridstruct%area_64
    fC            => gridstruct%fC
    f0            => gridstruct%f0
    cubed_sphere  => gridstruct%cubed_sphere

    pe(:,:,:) = 0.0
    pt(:,:,:) = 1.0
    f0(:,:)   = HUGE(0.0)
    fC(:,:)   = HUGE(0.0)

    ! set Coriolis terms
    do j=jsd,jed+1
      do i=isd,ied+1
        fC(i,j) = 2.*omega*( -1.*COS(grid(i,j,1))*COS(grid(i,j,2))*SIN(alpha) + &
                                 SIN(grid(i,j,2))*COS(alpha) )
      enddo
    enddo
    do j=jsd,jed
      do i=isd,ied
        f0(i,j) = 2.*omega*( -1.*COS(agrid(i,j,1))*COS(agrid(i,j,2))*SIN(alpha) + &
                                 SIN(agrid(i,j,2))*COS(alpha) )
      enddo
    enddo

    call mpp_update_domains( f0, domain )
    if (cubed_sphere) call fill_corners(f0, npx, npy, YDir)

    delp(isd:is-1,jsd:js-1,1:npz)=0.0
    delp(isd:is-1,je+1:jed,1:npz)=0.0
    delp(ie+1:ied,jsd:js-1,1:npz)=0.0
    delp(ie+1:ied,je+1:jed,1:npz)=0.0

    ! Initialize the model state according to the specified
    ! test case.
    select case(selected_case)
      ! case(CASE_DIV)
      ! case(CASE_NL_DEFORM)
      ! case(CASE_COS_BELL)
      ! case(CASE_ZONAL)
      ! case(CASE_ZONAL_MTN)
      ! case(CASE_NR_POT)
      ! case(CASE_ROSSBY_4)
      ! case(CASE_BT_INSTAB)
      ! case(CASE_SOLITON)
      ! case(CASE_POLAR_VTX)
      case(CASE_IDEAL_MTN)
        call ic_ideal_mtn(u, v, w, pt, delp, delz, q, &
                                phis, ps, ak, bk, gridstruct, domain)
      case(CASE_MULTI_VTX)
        call ic_multi_vtx(u, v, w, pt, delp, delz, q, &
                                pk, peln, pe, pkz, phis, ps, &
                                ak, bk, gridstruct, domain)
    end select

    call mpp_update_domains( phis, domain )

    ftop = g_sum(domain, phis(is:ie,js:je), is, ie, js, je, ng, area, 1)
    if(is_master()) write(*,*) 'mean terrain height (m)=', ftop/grav

    ! The flow is initially hydrostatic
#ifndef SUPER_K
    call p_var(npz, is, ie, js, je, ptop, ptop_min, delp, delz, pt, ps,   &
               pe, peln, pk, pkz, kappa, q, ng, ncnst, area, dry_mass, .false., mountain, &
               moist_phys, hydrostatic, nwat, domain, .not.hydrostatic)
#endif

#ifdef COLUMN_TRACER

#endif

    call mpp_update_domains(u, v, domain, gridtype=DGRID_NE, complete=.true.)

    nullify(agrid)
    nullify(grid)
    nullify(area)
    nullify(fC)
    nullify(f0)
    nullify(cubed_sphere)

  end subroutine ic_test

end module test_init_mod
