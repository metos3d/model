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
subroutine metos3dbgcinit(ny, nx, nu, nb, nd, dt, q, t, y, u, b, d, ndiag, diag)
    implicit none
    ! input variables
    integer :: ny, nx, nu, nb, nd, ndiag
    real(8) :: dt, q(nx, ny), t, y(nx, ny), u(nu), b(nb), d(nx, nd), diag(nx, ndiag)
end subroutine

!
!   metos3dbgcfinal
!
subroutine metos3dbgcfinal(ny, nx, nu, nb, nd, dt, q, t, y, u, b, d, ndiag, diag)
    implicit none
    ! input variables
    integer :: ny, nx, nu, nb, nd, ndiag
    real(8) :: dt, q(nx, ny), t, y(nx, ny), u(nu), b(nb), d(nx, nd), diag(nx, ndiag)
end subroutine

!
!   metos3dbgc
!
subroutine metos3dbgc(ny, nx, nu, nb, nd, dt, q, t, y, u, b, d, ndiag, diag)
    implicit none
    ! input variables
    integer :: ny, nx, nu, nb, nd, ndiag
    real(8) :: dt, q(nx, ny), t, y(nx, ny), u(nu), b(nb), d(nx, nd), diag(nx, ndiag)
end subroutine

