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

    ! NP-DOP model
    call NPDOPmodel(n, nz, m, dt, q, t, y(1,1), y(1,2), y(1,3), u, bc(1), bc(2), dc(1,1), dc(1,2))

end subroutine

!
!   NP-DOP model
!
subroutine NPDOPmodel(n, nz, m, dt, q, t, yN, yP, yDOP, u, phi, sigmaice, z, dz)
    implicit none
    ! input variables
    integer :: n, nz, m
    real*8  :: dt, q(nz, n), t, yN(nz), yP(nz), yDOP(nz), u(m), phi, sigmaice, z(nz), dz(nz)

    ! constants
    integer, parameter :: jeuphotic = 2                         ! number of euphotic layers
    real*8, parameter  :: sigmaPAR  = 0.4d0                     ! photosynthetically available radiation (PAR)
    real*8, parameter  :: yZstar    = 0.01d0                    ! [mmolP/m^3]

    ! work vars
    integer :: j, k
    real*8  :: yNj, yPj, yDOPj, ISWR, IPj, IPjprime, IPkm1, IPk
    real*8  :: kw, kc, muP, muZ, KN, KP, KI, sigmaDOP, lambdaP, kappaP, lambdaPprime, lambdaDOPprime, b
    real*8  :: fP, fZ, E
    real*8  :: sigmaDOPbar, dtbio

    ! retrieve and scale parameters
    kw              = u(1)          ! attenuation of water                  [1/m]
    kc              = u(2)          ! attenuation of chlorophyll (P)        [1/m (m^3/mmolP)]
    muP             = u(3)          ! maximum groth rate P                  [1/d]
    muZ             = u(4)          ! maximum groth rate Z                  [1/d]
    KN              = u(5)          ! N half saturation                     [mmolP/m^3]
    KP              = u(6)          ! P half saturation                     [mmolP/m^3]
    KI              = u(7)          ! I half satuartion                     [W/m^2]
    sigmaDOP        = u(8)          ! fraction of DOP                       [1]
    lambdaP         = u(9)          ! linear loss rate P (euphotic)         [1/d]
    kappaP          = u(10)         ! quadratic loss rate P (euphotic)      [1/d (m^3/mmolP)]
    lambdaPprime    = u(11)         ! linear loss rate Z (all layers)       [1/d]
    lambdaDOPprime  = u(12)/360.d0  ! DOP reminalization rate (all layers)  [1/y]
    b               = u(13)         ! power law coefficient                 [1]

    ! compute insolation
    ! take PAR and ice cover into account
    ! initialize product
    call insolation(t, phi, ISWR)
    ISWR  = sigmaPAR * (1.d0 - sigmaice) * ISWR
    IPkm1 = 1.d0
    IPk   = 1.d0

    ! euphotic zone
    sigmaDOPbar = (1.d0 - sigmaDOP)
    do j = 1, min(jeuphotic, nz)

        ! ensure positive values (or zero)
        yNj   = max(yN(j), 0.d0)
        yPj   = max(yP(j), 0.d0)
        yDOPj = max(yDOP(j), 0.d0)

        ! attenaution of water
        ! build product for layers above current layer
        ! store value of current *full* layer
        ! set value of current *half* layer
        IPk      = IPk * IPkm1
        IPkm1    = exp(-(kw + kc * yPj) * dz(j))
        IPjprime = exp(-(kw + kc * yPj) * 0.5d0 * dz(j))
        ! combine factors
        IPj = ISWR * IPjprime * IPk

        ! production
        fP = muP * yPj * yNj / (KN + yNj) * IPj / (KI + IPj)
        fZ = muZ * yZstar * yPj*yPj / (KP*KP + yPj*yPj)

        ! uptake
        q(j, 1) = q(j, 1) - fP
        q(j, 2) = q(j, 2) + fP - fZ            - lambdaP * yPj - kappaP * yPj*yPj
        q(j, 3) = q(j, 3)      + sigmaDOP * fZ + lambdaP * yPj + kappaP * yPj*yPj

        ! remineralization of last *euphotic* layer
        E = sigmaDOPbar * fZ
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
        yPj   = max(yP(j), 0.d0)
        yDOPj = max(yDOP(j), 0.d0)
        ! reminalization
        q(j, 1) = q(j, 1)                      + lambdaDOPprime * yDOPj
        q(j, 2) = q(j, 2) - lambdaPprime * yPj
        q(j, 3) = q(j, 3) + lambdaPprime * yPj - lambdaDOPprime * yDOPj
    end do

    ! scale with *bio* time step
    dtbio = dt * 360.d0
    do j = 1, nz
        q(j, 1) = q(j, 1) * dtbio
        q(j, 2) = q(j, 2) * dtbio
        q(j, 3) = q(j, 3) * dtbio
    end do

end subroutine


