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
!   metos3dbgcbegin
!
subroutine metos3dbgcbegin(ny, nx, nu, nb, nd, dt, q, t, y, u, b, d, ndiag, diag)
    implicit none
    ! input variables
    integer :: ny, nx, nu, nb, nd, ndiag
    real(8) :: dt, q(nx, ny), t, y(nx, ny), u(nu), b(nb), d(nx, nd), diag(nx, ndiag)

    ! diagnostics
    ! init variables to zero at the beginning
    diag(1:nx, 1:ndiag) = 0.0d0

end subroutine

!
!   metos3dbgc
!
subroutine metos3dbgc(ny, nx, nu, nb, nd, dt, q, t, y, u, b, d, ndiag, diag)
    implicit none
    ! input variables
    integer :: ny, nx, nu, nb, nd, ndiag
    real*8  :: dt, q(nx, ny), t, y(nx, ny), u(nu), b(nb), d(nx, nd), diag(nx, ndiag)

    ! N model
    call Nmodel(ny, nx, nu, dt, q, t, y(1,1), u, b(1), b(2), d(1,1), d(1,2), ndiag, diag)

end subroutine

!
!   N model
!
subroutine Nmodel(n, nz, m, dt, q, t, yN, u, phi, sigmaice, z, dz, ndiag, diag)
    implicit none
    ! input variables
    integer :: n, nz, m
    real*8  :: dt, q(nz, n), t, yN(nz), u(m), phi, sigmaice, z(nz), dz(nz)

    ! diagnostics
    integer  :: ndiag
    integer, parameter :: igome = 1
    integer, parameter :: igpp  = 2
    real(8) :: diag(nz, ndiag)

    ! constants
    integer, parameter :: jeuphotic = 2                         ! number of euphotic layers
    real*8, parameter  :: sigmaPAR  = 0.4d0                     ! photosynthetically available radiation (PAR)
    real*8, parameter  :: yPstar    = 0.002777777777777778d0    ! [mmolP/m^3]

    ! work vars
    integer :: j, k
    real*8  :: yNj, ISWR, Ij, Ijprime, Ikm1, Ik
    real*8  :: kw, muP, KN, KI, b
    real*8  :: fP, E
    real*8  :: dtbio

    ! retrieve and scale parameters
    kw  = u(1)          ! attenuation of water      [1/m]
    muP = u(2)          ! maximum groth rate P      [1/d]
    KN  = u(3)          ! N half saturation         [mmolP/m^3]
    KI  = u(4)          ! I half satuartion         [W/m^2]
    b   = u(5)          ! power law coefficient     [1]

    ! compute insolation
    ! take PAR and ice cover into account
    ! initialize product
    call insolation(t, phi, ISWR)
    ISWR = sigmaPAR * (1.d0 - sigmaice) * ISWR
    Ikm1 = 1.d0
    Ik   = 1.d0

    ! euphotic zone
    do j = 1, min(jeuphotic, nz)

        ! attenaution of water
        ! build product for layers above current layer
        ! store value of current *full* layer
        ! set value of current *half* layer
        Ik      = Ik * Ikm1
        Ikm1    = exp(-kw * dz(j))
        Ijprime = exp(-kw * 0.5d0 * dz(j))
        ! combine factors
        Ij = ISWR * Ijprime * Ik

        ! ensure positive values (or zero)
        yNj = max(yN(j), 0.d0)

        ! production
        fP = muP * yPstar * yNj / (KN + yNj) * Ij / (KI + Ij)

        ! diagnostics
        diag(j, igpp) = diag(j, igpp) + fP * dt * 360.d0

        ! uptake
        q(j, 1) = q(j, 1) - fP

        ! remineralization
        E = fP
        if (j == nz) then
            ! last *euphotic* layer
            q(j, 1) = q(j, 1) + E
        else
            ! export to layers below
            do k = j+1, nz
                ! approximate derivative d/dz
                if (k == nz) then
                    ! last layer
                    q(k, 1) = q(k, 1) + E * dz(j) * (z(k-1)/z(j))**(-b) / dz(k)
if (k > jeuphotic) then
diag(j, igome) = diag(j, igome) + E * dz(j) * (z(k-1)/z(j))**(-b) / dz(k) * dt * 360.d0
endif
                else
                    ! layers in between
                    q(k, 1) = q(k, 1) + E * dz(j) * ((z(k-1)/z(j))**(-b) - (z(k)/z(j))**(-b)) / dz(k)
if (k > jeuphotic) then
diag(j, igome) = diag(j, igome) + E * dz(j) * ((z(k-1)/z(j))**(-b) - (z(k)/z(j))**(-b)) / dz(k) * dt * 360.d0
endif
                end if
            end do
        end if
    end do

    ! scale with *bio* time step
    dtbio = dt * 360.d0
    do j = 1, nz
        q(j, 1) = q(j, 1) * dtbio
    end do

end subroutine


