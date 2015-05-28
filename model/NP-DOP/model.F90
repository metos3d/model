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

    ! NP-DOP model
    call bgc_po4_phy_dop(n, nz, m, dt, q, t, y(1,1), y(1,2), y(1,3), u, bc(1), bc(2), dc(1,1), dc(1,2))

end subroutine

#include "insolation.F90"

!
!   NP-DOP model
!
subroutine bgc_po4_phy_dop(n, nz, m, dt, q, t, po4, phy, dop, u, latitude, icecover, z, dz)
    implicit none
    ! input variables
    integer :: n, nz, m
    real*8  :: dt, q(nz, n), t, po4(nz), dop(nz), phy(nz), u(m), latitude, icecover, z(nz), dz(nz)

    ! constants
    integer, parameter :: jeuphotic = 2
    real*8, parameter  :: sig_par   = 0.4d0                     ! photosynthetically available radiation (PAR)
    real*8, parameter  :: y_z_star  = 0.01d0                    ! [mmolP/m^3]

    ! work vars
    integer :: j, k
    real*8  :: po4_j, dop_j, phy_j, iswr, swr_j
    real*8  :: k_w, k_p, mu_p, mu_z, k_po4, k_z, k_swr, sig_dop, lam_p, kap_p, lamp_p, lamp_dop, b
    real*8  :: f_p, f_z

    ! retrieve and scale parameters
    k_w         = u(1)          ! attenuation of water          [1/m]
    k_p         = u(2)          ! attenuation of phytoplankton  [1/m]
    mu_p        = u(3)          ! maximum groth rate PHY        [1/d]
    mu_z        = u(4)          ! maximum groth rate ZOO        [1/d]
    k_po4       = u(5)          ! po4 half saturation           [mmolP/m^3]
    k_z         = u(6)          ! zoo half saturation           [mmolP/m^3]
    k_swr       = u(7)          ! light half satuartion         [W/m^2]
    sig_dop     = u(8)          ! fraction of dop               [1]
    lam_p       = u(9)          ! linear loss rate PHY
    kap_p       = u(10)         ! quadratic loss rate PHY
    lamp_p      = u(11)         ! linear loss rate (all layers)
    lamp_dop    = u(12)/360.d0  ! dop reminalization rate       [1/y]
    b           = u(13)         ! power law coefficient         [1]

    ! compute insolation, take PAR and ice cover into account
    call insolation(t, latitude, iswr)
    iswr = sig_par * (1.d0 - icecover) * iswr

    ! euphotic zone
    do j = 1, min(jeuphotic, nz)
        ! ensure positive values (or zero)
        po4_j = max(po4(j), 0.d0)
        phy_j = max(phy(j), 0.d0)
        dop_j = max(dop(j), 0.d0)
        ! attenaution of water
        if (j == 1) then
            ! first layer
            swr_j = exp(-(k_w + k_p * phy_j)* 0.5d0 * dz(j)) * iswr
        else
            ! other layers
            swr_j = exp(-(k_w + k_p * phy_j) * dz(j)) * swr_j
        end if
        ! production
        f_p = mu_p * phy_j * po4_j / (po4_j + k_po4) * swr_j / (swr_j + k_swr)
        f_z = mu_z * y_z_star * phy_j**2 / (phy_j**2 + k_z**2)

        ! uptake
        q(j, 1) = q(j, 1) - f_p
        q(j, 2) = q(j, 2) + f_p - f_z           - lam_p * phy_j - kap_p * phy_j**2
        q(j, 3) = q(j, 3)       + sig_dop * f_z + lam_p * phy_j + kap_p * phy_j**2

        ! remineralization of last *euphotic* layer
        if (j == nz) then
            q(j, 1) = q(j, 1) + (1.d0 - sig_dop) * f_z
        else
            ! export to layers below
            do k = j+1, nz
                ! approximate derivative d/dz
                if (k == nz) then
                    ! last layer
                    q(k, 1) = q(k, 1) + (1.d0 - sig_dop) * f_z * dz(j) * (z(k-1)/z(j))**(-b) / dz(k)
                else
                    ! layers in between
                    q(k, 1) = q(k, 1) + (1.d0 - sig_dop) * f_z * dz(j) * ((z(k-1)/z(j))**(-b) - (z(k)/z(j))**(-b)) / dz(k)
                end if
            end do
        end if
    end do

    ! all layers
    do j = 1, nz
        ! ensure positive values (or zero)
        po4_j = max(po4(j), 0.d0)
        phy_j = max(phy(j), 0.d0)
        dop_j = max(dop(j), 0.d0)
        ! reminalization
        q(j, 1) = q(j, 1) + lamp_dop * dop_j
        q(j, 2) = q(j, 2)                    - lamp_p * phy_j
        q(j, 3) = q(j, 3) - lamp_dop * dop_j + lamp_p * phy_j
    end do

    ! scale with *bio* time step
    do j = 1, nz
        q(j, 1) = q(j, 1) * dt * 360.d0
        q(j, 2) = q(j, 2) * dt * 360.d0
        q(j, 3) = q(j, 3) * dt * 360.d0
    end do

end subroutine


