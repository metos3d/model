!
! Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D
! Copyright (C) 2014  Jaroslaw Piwonski, CAU, jpi@informatik.uni-kiel.de
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!

!
!   metos3dbgcinit
!
subroutine metos3dbgcinit(n, nz, m, nbc, ndc, dt, q, t, y, u, bc, dc)
    implicit none
    ! input variables
    integer :: n, nz, m, nbc, ndc
    real*8  :: dt, q(nz, n), t, y(nz, n), u(m), bc(nbc), dc(nz, ndc)

!    ! translate metos3d names to mops names
!    Nrloc                           = nz                ! # of layers
!    DeltaT                          = dt*360.0*86400    ! convert fraction of year (dt) to fraction of day in seconds (DeltaT), [dt] = yr, [DeltaT] = s
!    ! ~/.metos3d/model/model/MOPS-2.0/option/test.MOPS-2.0.option.txt
!    ! ...
!    ! -Metos3DDomainConditionName LayerDepth,LayerHeight,Temperature,Salinity
!    ! ...
!    drFloc                          = dc(:, 1)          ! layer bottom depth
!    dzloc                           = dc(:, 2)          ! layer thickness
!    thetaloc                        = dc(:, 3)          ! temperature
!    saltloc                         = dc(:, 4)          ! salinity
!    nzmax                           = 15                ! maximum # of layers
!    nzeuph                          = 2                 ! # of euphotic layers
!    numBiogeochemStepsPerOceanStep  = 8                 ! # of time
!    setDefaults                     = .false.

    INTEGER bgc_ktotal
    PARAMETER(bgc_ktotal=100)
    INTEGER bgc_ntracer
    PARAMETER(bgc_ntracer=7)
    REAL*8 bgc_tracer
    COMMON/BGC/bgc_tracer(bgc_ktotal,bgc_ntracer)
    real*8 bgc_dt
    COMMON/BGCCONTROL/bgc_dt

    ! translate metos3d names to mops names, if possible
    integer :: nzmax = 15
    integer :: nzeuph = 2
    integer :: numBiogeochemStepsPerOceanStep = 1
    logical :: setDefaults = .true.

    bgc_tracer(:,:)       = 0.d0
    bgc_tracer(1:nz, 1:n) = y(1:nz, 1:n)

    ! call original init routine
    call MOPS_BIOGEOCHEM_INI( &
        nz,                             &   ! Nrloc,      &
        dt*360.0*86400,                 &   ! DeltaT,     &
        dc(:, 3),                       &   ! thetaloc,   &
        dc(:, 4),                       &   ! saltloc,    &
        dc(:, 2),                       &   ! dzloc,      &
        dc(:, 1),                       &   ! drFloc,     &
        nzmax,                          &
        nzeuph,                         &
        numBiogeochemStepsPerOceanStep, &
        setDefaults                     &
    )

!    print *, 'n, nz:', n, nz
!    print *, 'bgc_dt:', bgc_dt

end subroutine

!
!   metos3dbgcfinal
!
subroutine metos3dbgcfinal(n, nz, m, nbc, ndc, dt, q, t, y, u, bc, dc)
    implicit none
    ! input variables
    integer :: n, nz, m, nbc, ndc
    real*8  :: dt, q(nz, n), t, y(nz, n), u(m), bc(nbc), dc(nz, ndc)
end subroutine

!
!   metos3dbgc
!
subroutine metos3dbgc(n, nz, m, nbc, ndc, dt, q, t, y, u, bc, dc)
    implicit none
    ! input variables
    integer :: n, nz, m, nbc, ndc
    real*8  :: dt, q(nz, n), t, y(nz, n), u(m), bc(nbc), dc(nz, ndc)

    INTEGER bgc_ktotal
    PARAMETER(bgc_ktotal=100)
    INTEGER bgc_ntracer
    PARAMETER(bgc_ntracer=7)
    REAL*8 bgc_tracer
    COMMON/BGC/bgc_tracer(bgc_ktotal,bgc_ntracer)

    real*8 bgc_dt
    COMMON/BGCCONTROL/bgc_dt

!    ! translate metos3d names to mops names
!    Nrloc                           = nz                ! # of layers
!    DeltaT                          = dt*360.0*86400    ! convert fraction of year (dt) to fraction of day in seconds (DeltaT)
!    ! ~/.metos3d/model/model/MOPS-2.0/option/test.MOPS-2.0.option.txt
!    ! ...
!    ! -Metos3DDomainConditionName LayerDepth,LayerHeight,Temperature,Salinity
!    ! ...
!    drFloc                          = dc(:, 1)          ! layer bottom depth
!    dzloc                           = dc(:, 2)          ! layer thickness
!    thetaloc                        = dc(:, 3)          ! temperature
!    saltloc                         = dc(:, 4)          ! salinity
!    nzmax                           = 15                ! maximum # of layers
!    nzeuph                          = 2                 ! # of euphotic layers
!    numBiogeochemStepsPerOceanStep  = 8                 ! # of time
!    setDefaults                     = .false.

    ! original insolation routine,
    ! computes all profiles at once,
    ! is embedded in C code,
    ! here, we compute only one profile and use the fortran call

    ! translate metos3d names to mops names, if possible, otherwise define them
    ! ~/.metos3d/model/model/MOPS-2.0/option/test.MOPS-2.0.option.txt
    ! ...
    ! -Metos3DBoundaryConditionName Latitude,IceCover,WindSpeed,AtmospherePressure
    ! ...
    real*8 :: SWRADloc
    real*8 :: TAUloc
    real*8 :: localburial = 0.d0
    real*8 :: globalrunoff = 0.0d0
    call insolation( &
        1,                  &   ! &lNumProfiles,
        dt*360.0*86400,     &   ! &myTime,
        bc(1),              &   ! &locallatitude[0],
        SWRADloc,           &   ! &localswrad[0],
        TAUloc              &   ! &localtau[0]);
    )

!    print *, SWRADloc, TAUloc
!

!    print *, 'bgc_dt:', bgc_dt

!    print *, 'in:'
!    print *, 'bgc_tracer(1:nz,1:n)'
!    print *, bgc_tracer(1:nz,1:n)
!!    call exit(1)

    bgc_tracer(1:nz, 1:n) = y(1:nz, 1:n)

    ! original model routine
    call MOPS_BIOGEOCHEM_MODEL( &
        nz,                 &   ! Nrloc,          &
        dt*360.0*86400,     &   ! DeltaT,         &

        dc(:, 3),           &   ! thetaloc,       &
        dc(:, 4),           &   ! saltloc,        &

        bc(2),              &   ! FIceloc,        &

        SWRADloc,           &
        TAUloc,             &

        bc(3),              &   ! WINDloc,        &
        bc(4),              &   ! ATMOSPloc,      &

        dc(:, 2),           &   ! dzloc,          &

        localburial,        &
        globalrunoff,       &
        dc(:, 2),           &   ! localrunoffloc, &

        .true.              &   ! useSeparateBiogeochemTS &, not used
    )

    q(1:nz, 1:n) = bgc_tracer(1:nz, 1:n) - y(1:nz, 1:n)

!    print *, 'out:'
!    print *, 'bgc_tracer(1:nz,1:n)'
!    print *, bgc_tracer(1:nz,1:n)
!    call exit(1)

end subroutine


