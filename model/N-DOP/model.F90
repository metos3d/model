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
subroutine metos3dbgc(n, nz, m, nbc, ndc, dt, q, t, y, u, bc, dc, ndiag, diag)
    implicit none
    ! input variables
    integer :: n, nz, m, nbc, ndc, ndiag
    real*8  :: dt, q(nz, n), t, y(nz, n), u(m), bc(nbc), dc(nz, ndc), diag(nz, ndiag)

    ! N-DOP model
    call NDOPmodel(n, nz, m, dt, q, t, y(1,1), y(1,2), u, bc(1), bc(2), dc(1,1), dc(1,2))

end subroutine

!
!   N-DOP model
!
subroutine NDOPmodel(n, nz, m, dt, q, t, yN, yDOP, u, phi, sigmaice, z, dz)
    implicit none
    ! input variables
    integer :: n, nz, m
    real*8  :: dt, q(nz, n), t, yN(nz), yDOP(nz), u(m), phi, sigmaice, z(nz), dz(nz)

    ! constants
    integer, parameter :: jeuphotic = 2                         ! number of euphotic layers
    real*8, parameter  :: sigmaPAR  = 0.4d0                     ! photosynthetically available radiation (PAR)
    real*8, parameter  :: yPstar    = 0.002777777777777778d0    ! [mmolP/m^3]

    ! work vars
    integer :: j, k
    real*8  :: yNj, yDOPj, ISWR, Ij, Ijprime, Ikm1, Ik
    real*8  :: kw, muP, KN, KI, sigmaDOP, lambdaDOPprime, b
    real*8  :: fP, E
    real*8  :: sigmaDOPbar, dtbio

    ! retrieve and scale parameters
    kw              = u(1)          ! attenuation of water      [1/m]
    muP             = u(2)          ! maximum groth rate P      [1/d]
    KN              = u(3)          ! N half saturation         [mmolP/m^3]
    KI              = u(4)          ! I half satuartion         [W/m^2]
    sigmaDOP        = u(5)          ! fraction of DOP           [1]
    lambdaDOPprime  = u(6)/360.d0   ! DOP reminalization rate   [1/y]
    b               = u(7)          ! power law coefficient     [1]

    ! compute insolation
    ! take PAR and ice cover into account
    ! initialize product
    call insolation(t, phi, ISWR)
    ISWR = sigmaPAR * (1.d0 - sigmaice) * ISWR
    Ikm1 = 1.d0
    Ik   = 1.d0

    ! euphotic zone
    sigmaDOPbar = (1.d0 - sigmaDOP)
    do j = 1, min(jeuphotic, nz)

        ! ensure positive values (or zero)
        yNj   = max(yN(j), 0.d0)
        yDOPj = max(yDOP(j), 0.d0)

        ! attenaution of water
        ! build product for layers above current layer
        ! store value of current *full* layer
        ! set value of current *half* layer
        Ik      = Ik * Ikm1
        Ikm1    = exp(-kw * dz(j))
        Ijprime = exp(-kw * 0.5d0 * dz(j))
        ! combine factors
        Ij = ISWR * Ijprime * Ik

        ! production
        fP = muP * yPstar * yNj / (KN + yNj) * Ij / (KI + Ij)

        ! uptake
        q(j, 1) = q(j, 1) - fP
        q(j, 2) = q(j, 2) + sigmaDOP * fP

        ! remineralization of last *euphotic* layer
        E = sigmaDOPbar * fP
        if (j == nz) then
            q(j, 1) = q(j, 1) + E
        else
            ! export to layers below
            do k = j+1, nz
                ! approximate derivative d/dz
                if (k == nz) then
                    ! last layer
                    q(k, 1) = q(k, 1) + E * dz(j) * (z(k-1)/z(j))**(-b) / dz(k)
                else
                    ! layers in between
                    q(k, 1) = q(k, 1) + E * dz(j) * ((z(k-1)/z(j))**(-b) - (z(k)/z(j))**(-b)) / dz(k)
                end if
            end do
        end if
    end do

    ! all layers
    do j = 1, nz
        ! ensure positive values (or zero)
        yNj   = max(yN(j), 0.d0)
        yDOPj = max(yDOP(j), 0.d0)
        ! reminalization
        q(j, 1) = q(j, 1) + lambdaDOPprime * yDOPj
        q(j, 2) = q(j, 2) - lambdaDOPprime * yDOPj
    end do

    ! scale with *bio* time step
    dtbio = dt * 360.d0
    do j = 1, nz
        q(j, 1) = q(j, 1) * dtbio
        q(j, 2) = q(j, 2) * dtbio
    end do

end subroutine


