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

#include "insolation.F"
#include "BGC_INI.F"
#include "mops_biogeochem_ini.F"
#include "BGC_MODEL.F"
#include "mops_biogeochem_model.F"

      !
      !   metos3dbgcinit
      !
      subroutine metos3dbgcinit(ny, nx, nu, nb, nd, dt,
     &  q, t, y, u, b, d, ndiag, diag)
        implicit none
        ! input variables
        integer :: ny, nx, nu, nb, nd, ndiag
        real(8) :: dt, q(nx, ny), t, y(nx, ny), u(nu)
        real(8) :: b(nb), d(nx, nd), diag(nx, ndiag)

        ! mops variables
        INTEGER bgc_ktotal
        PARAMETER(bgc_ktotal=100)
        INTEGER bgc_ntracer
        PARAMETER(bgc_ntracer=7)
        REAL*8 bgc_tracer
        COMMON/BGC/bgc_tracer(bgc_ktotal,bgc_ntracer)
        real*8 bgc_dt
        COMMON/BGCCONTROL/bgc_dt

        ! translate metos3d names to mops names, if possible
        ! use option file ...
        ! ~/.metos3d/model/model/MOPS-2.0/option/test.MOPS-2.0.option.txt
        ! ... and therein look for the domain data ...
        ! -Metos3DDomainConditionName LayerDepth,LayerHeight,Temperature,Salinity
        integer :: nzmax = 15
        integer :: nzeuph = 2
        integer :: numBiogeochemStepsPerOceanStep = 1
        logical :: setDefaults = .true.

        ! zero all tracers
        ! copy metos3d tracers to mops tracers
        bgc_tracer(:,:)       = 0.d0
        bgc_tracer(1:nx, 1:ny) = y(1:nx, 1:ny)

        ! call original init routine
        call MOPS_BIOGEOCHEM_INI(
     &       nx,                                ! Nrloc,      &
     &       dt*360.0*86400,                    ! DeltaT,     &
     &       d(:, 3),                           ! thetaloc,   &
     &       d(:, 4),                           ! saltloc,    &
     &       d(:, 2),                           ! dzloc,      &
     &       d(:, 1),                           ! drFloc,     &
     &       nzmax,
     &       nzeuph,
     &       numBiogeochemStepsPerOceanStep,
     &       setDefaults
     &   )
      end subroutine

!
!   metos3dbgcfinal
!
      subroutine metos3dbgcfinal(
     &    n, nz, m, nbc, ndc, dt, q, t, y, u, bc, dc)
        implicit none
        ! input variables
        integer :: n, nz, m, nbc, ndc
        real*8  :: dt, q(nz, n), t, y(nz, n), u(m), bc(nbc), dc(nz, ndc)
      end subroutine

!
!   metos3dbgc
!
      subroutine metos3dbgc(
     &    n, nz, m, nbc, ndc, dt, q, t, y, u, bc, dc)
        implicit none
        ! input variables
        integer :: n, nz, m, nbc, ndc
        real*8  :: dt, q(nz, n), t, y(nz, n), u(m), bc(nbc), dc(nz, ndc)

        ! mops variables
        INTEGER bgc_ktotal
        PARAMETER(bgc_ktotal=100)
        INTEGER bgc_ntracer
        PARAMETER(bgc_ntracer=7)
        REAL*8 bgc_tracer
        COMMON/BGC/bgc_tracer(bgc_ktotal,bgc_ntracer)

        ! original insolation routine,
        ! computes all profiles at once,
        ! is embedded in C code,
        ! here, we compute only one profile and use the fortran call

        ! translate metos3d names to mops names, if possible, otherwise define them
        ! use option file ...
        ! ~/.metos3d/model/model/MOPS-2.0/option/test.MOPS-2.0.option.txt
        ! ... and therein look for the boundary data ...
        ! -Metos3DBoundaryConditionName Latitude,IceCover,WindSpeed,AtmospherePressure
        real*8 :: SWRADloc
        real*8 :: TAUloc
        real*8 :: localburial = 0.d0
        real*8 :: globalrunoff = 0.0d0
        call insolation(
     &       1,                     ! &lNumProfiles,
     &       dt*360.0*86400,        ! &myTime,
     &       bc(1),                 ! &locallatitude[0],
     &       SWRADloc,              ! &localswrad[0],
     &       TAUloc                 ! &localtau[0]);
     &   )

        ! copy metos3d tracers to mops tracers
        bgc_tracer(1:nz, 1:n) = y(1:nz, 1:n)

        ! original model routine
        call MOPS_BIOGEOCHEM_MODEL(
     &       nz,                    ! Nrloc,          &
     &       dt*360.0*86400,        ! DeltaT,         &
     &       dc(:, 3),              ! thetaloc,       &
     &       dc(:, 4),              ! saltloc,        &
     &       bc(2),                 ! FIceloc,        &
     &       SWRADloc,
     &       TAUloc,
     &       bc(3),                 ! WINDloc,        &
     &       bc(4),                 ! ATMOSPloc,      &
     &       dc(:, 2),              ! dzloc,          &
     &       localburial,
     &       globalrunoff,
     &       dc(:, 2),              ! localrunoffloc, &
     &       .true.                 ! useSeparateBiogeochemTS &, not used
     &   )

        ! copy mops tracers back to metos3d tracers
        ! subtract input from mops result,
        ! since we only want the source minus sink term
        q(1:nz, 1:n) = bgc_tracer(1:nz, 1:n) - y(1:nz, 1:n)

      end subroutine


