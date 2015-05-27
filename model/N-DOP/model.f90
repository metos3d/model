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

    ! N-DOP model
    call bgc_po4_dop(n, nz, m, dt, q, t, y(1,1), y(1,2), u, bc(1), bc(2), dc(1,1), dc(1,2))

end subroutine

#include "insolation.f90"

!
!   N-DOP
!
subroutine bgc_po4_dop(n, nz, m, dt, q, t, po4, dop, u, latitude, icecover, z, dz)
    implicit none
    ! input variables
    integer :: n, nz, m
    real*8  :: dt, q(nz, n), t, po4(nz), dop(nz), u(m), latitude, icecover, z(nz), dz(nz)

    ! constants
    integer, parameter :: jeuphotic = 2

    ! work vars
    integer :: iz
    real*8  :: po4_i, dop_i, iswr

    ! compute insolation
    call insolation(t, latitude, icecover, iswr)

    print *, t, latitude, icecover, iswr

    ! euphotic zone
    do iz = 1, min(jeuphotic, nz)
        ! ensure positive values (or zero)
        po4_i = max(po4(iz), 0.d0)
        dop_i = max(dop(iz), 0.d0)
        !


    end do




end subroutine



