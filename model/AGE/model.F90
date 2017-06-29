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
!   metos3dbgc
!
subroutine metos3dbgc(ny, nx, nu, nb, nd, dt, q, t, y, u, b, d, ndiag, diag)
    implicit none
    ! input variables
    integer :: ny, nx, nu, nb, nd, ndiag
    real(8) :: dt, q(nx, ny), t, y(nx, ny), u(nu), b(nb), d(nx, nd), diag(nx, ndiag)

    !
    !  AGE tracer, d/dt y = 1
    !
    !  means:
    !  we just add the time step to the prognostic variable,
    !  i.e. adding up to plus one (+1.0) at the end of the year,
    !  and thus aging while model years are computed / simulated
    !  and,
    !  we set the tracer to zero at the surface,
    !  to mark the water as 'young'
    !

    ! set the source minus sink term to the time step
    ! this will be added to the tracer
    q = dt

    ! we assume, we only have one tracer / prognostic variable
    ! first layer, first tracer
    ! set it exact to the opposite of the input,
    ! thus enforce zero
    q(1,1) = - y(1,1)

end subroutine


