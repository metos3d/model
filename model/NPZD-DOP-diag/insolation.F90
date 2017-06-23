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
!   insolation
!
subroutine insolation(t, latitude, iswr)
    implicit none
    ! arguments
    real*8  :: t, latitude, iswr
    ! constants
    real*8, parameter :: pi = 3.141592653589793d0
    real*8, parameter :: solar  = 1360.d0
    real*8, parameter :: albedo = 0.6d0
    ! work vars
    real*8  :: dayrad, delta, latrad, sun, dayhrs, cosz, dayfrac
    ! day of year in radians
    dayrad = 2.d0 * pi * t
    ! declination
    delta = (0.006918d0                             &
                - 0.399912d0 * cos(dayrad)          &
                + 0.070257d0 * sin(dayrad)          &
                - 0.006758d0 * cos(2.d0 * dayrad)   &
                + 0.000907d0 * sin(2.d0 * dayrad)   &
                - 0.002697d0 * cos(3.d0 * dayrad)   &
                + 0.001480d0 * sin(3.d0 * dayrad))
    ! latitude in radians
    latrad = latitude / 180.d0 * pi
    ! sun angle
    sun = - sin(delta) / cos(delta) * sin(latrad) / cos(latrad)
    ! bound values
    sun = max(sun, - 0.999d0)
    sun = min(sun, + 0.999d0)
    ! day hours
    dayhrs = abs(acos(sun))
    ! average zenith angle
    cosz = sin(delta) * sin(latrad) + cos(delta) * cos(latrad) * sin(dayhrs) / dayhrs
    ! bound value
    cosz = max(cosz, 0.005d0)
    ! fraction of daylight
    dayfrac = dayhrs / pi
    ! scale with constants
    iswr = solar * (1.0d0 - albedo) * dayfrac * cosz
    ! bound value
    iswr = max(iswr, 0.00001d0)
end subroutine









