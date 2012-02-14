!
! Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D
! Copyright (C) 2012  Jaroslaw Piwonski, CAU, jpi@informatik.uni-kiel.de
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
!   model.F
!
!   Model for Iodine and Caesium
!

      subroutine metos3dbgcinit(n, nz, m, nbc, ndc, dt,
     &  q, t, y, u, bc, dc)
        implicit none
        ! input variables
        integer :: n, nz, m, nbc, ndc
        real*8  :: dt, q(nz, n), t, y(nz, n), u(m)
        real*8  :: bc(nbc), dc(nz, ndc)
      end

      subroutine metos3dbgcfinal(n, nz, m, nbc, ndc, dt,
     &  q, t, y, u, bc, dc)
        implicit none
        ! input variables
        integer :: n, nz, m, nbc, ndc
        real*8  :: dt, q(nz, n), t, y(nz, n), u(m)
        real*8  :: bc(nbc), dc(nz, ndc)
      end

      subroutine metos3dbgc(n, nz, m, nbc, ndc, dt,
     &  q, t, y, u, bc, dc)
        implicit none
        ! input variables
        integer :: n, nz, m, nbc, ndc
        real*8  :: dt, q(nz, n), t, y(nz, n), u(m)
        real*8  :: bc(nbc), dc(nz, ndc)

        call bgc(
     &      t, dt, n, nz, nbc, bc, ndc, dc, y, q, m, u)

      end

      subroutine bgc(t, dt, ntracer, nlayer,
     &              nbc, bc, ndc, dc, y, ybgc, nparam, u)
        implicit none
        ! input variables
        real*8  :: t, dt
        integer :: ntracer, nlayer, nbc, ndc, nparam
        real*8  :: bc(nbc), dc(nlayer, ndc)
        real*8  :: y(nlayer, ntracer)
        real*8  :: ybgc(nlayer, ntracer)
        real*8  :: u(nparam)

        real*8  :: indicator, IOD, CAE
        real*8  :: daystart, daylength, entry
        integer :: ilayer

        do ilayer = 1, nlayer
            IOD = max(y(ilayer,1), 0.d0)
            CAE = max(y(ilayer,2), 0.d0)
            ybgc(ilayer,1) = dt * log(0.5d0)*360.d0/8.02070d0 * IOD
            ybgc(ilayer,2) = dt * log(0.5d0)*1.d0/30.17d0 * CAE        
        end do
        
        if (nbc > 0) then
            indicator = bc(1)
            if (indicator == 1.d0) then
                daystart = u(1)
                daylength = u(2)
                entry = u(3)
                
        if ((t > dt*8.d0*daystart)
     &      .and.(t < dt*8.d0*(daystart+daylength))) then
                    ybgc(1,1) = ybgc(1,1) + dt*entry
                    ybgc(1,2) = ybgc(1,2) + dt*entry
                end if
            end if
        end if
        
      end
