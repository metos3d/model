C
C Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D
C Copyright (C) 2012  Jaroslaw Piwonski, CAU, jpi@informatik.uni-kiel.de
C
C This program is free software: you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published by
C the Free Software Foundation, either version 3 of the License, or
C (at your option) any later version.
C
C This program is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C GNU General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this program.  If not, see <http://www.gnu.org/licenses/>.
C
C This file incorporates work covered by the following copyright and
C permission notice:
C
C     Copyright (C) 2012, MITgcm, http://mitgcm.org/
C     Copyright (C) 2012, Samar Khatiwala, spk@ldeo.columbia.edu
C     Copyright (C) 2012, Iris Kriest, ikriest@ifm-geomar.de
C
C     Permission is hereby granted, free of charge, to any person
C     obtaining a copy of this software and associated documentation
C     files (the "Software"), to deal in the Software without
C     restriction, including without limitation the rights to use,
C     copy, modify, merge, publish, distribute, sublicense, and/or
C     sell copies of the Software, and to permit persons to whom
C     the Software is furnished to do so, subject to the following
C     conditions:
C
C     The above copyright notice and this permission notice shall
C     be included in all copies or substantial portions of the Software.
C
C     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
C     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
C     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
C     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
C     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
C     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
C     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
C     OTHER DEALINGS IN THE SOFTWARE.
C
C   model.F
C
C   Model for nutrients (phosphate), N and dissolved organic phosphorous, DOP
C

#include "insolation.F"

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
        ! work variables
        integer :: j, k
        
        ! constants
        integer :: jeuphotic = 2
        real*8  :: biofrac   = 0.4d0            ! fraction of light spectrum that is used
        
        ! to compute, stau not used
        ! bc(1), latitude
        ! bc(2), ice cover
        real*8  :: insol, stau
        call insolation(1, t*360.d0*86400.d0, bc(1), insol, stau)
        insol = (1.0d0-bc(2))*biofrac*insol

        ! real value parameters
        ! to optimize
        real*8  :: lambda                       ! remineralization rate            : 1 / yr
        real*8  :: alpha                        ! maximum export rate              : uM / yr
        real*8  :: sigma                        ! fraction of DOP
        real*8  :: b                            ! Martin exponent
        ! not to optimize
        real*8  :: K_N                          ! half saturation constant: N      : uM
        real*8  :: K_insol                      ! half saturation constant: light  : W / m^2
        real*8  :: K_H2O                        ! attenuationOfWater

        ! set parameters
        lambda  = u(1)
        alpha   = u(2)
        nu      = u(3)
        b       = u(4)
        K_N     = 0.5d0
        K_insol = 30.d0
        K_H2O   = 0.02d0
        
        ! prepare model
        ! dc(:,1), layer heights
        ! dc(:,2), layer depths
        sigma_bar   = (1.0d0-sigma)
        minlayer    = min(nz, jeuphotic);
        
        ! fraction of parameters
        dtbio       = dt*360.0d0
        lambdabio   = lambda/360.0d0
        alphabio    = alpha/360.0d0

        ! euphotic layers
        do j = 1, minlayer
        
            ! attenaution of water
            if (j == 1) then
                ! first layer
                insol_j = insol
            else
                ! other layers
                insol_j = exp(-dc(j-1,2)*K_H2O)*insol
            end if

            ! take only positive values into account
            N_j = max(y(j,1), 0.d0)
            
            ! community production
            gamma_j =
     &  alphabio*N_j*insol_j/(N_j+K_N)/(insol_j+K_insol)
            
            ! consumption
            q(j,1) = q(j,1) - gamma_j
            q(j,2) = q(j,2) + sigma*gamma_j

            ! two cases
            if ((nz <= jeuphotic) .AND. (j == nz))then
                ! remineralization of last layer
                q(j,1) = q(j,1) + nu_sigma*gamma_j
            else
                ! export to layers below
                do k = j+1, nz
                    ! last layer
                    if (k == nz) then
                        q(k,1) = q(k,1)
      & + gamma_j*sigma_bar*dc(j,1)*(dc(k-1,2)/dc(j,2))**(-b)/dc(k,1)
                    ! layers in between
                    else
                        q(k,1) = q(k, 1)
      & + gamma_j*sigma_bar*dc(j,1)
      & * ((dc(k-1,2)/dc(j,2))**(-b)-(dc(k,2)/dc(j,2))**(-b))/dc(k,1)
                    end if
                end do
            end if
        end do
        
        ! remineralization
        do j = 1, nz
            ! take only positive values into account
            DOP_j = max(y(j,2),0.d0)
            ! update
            q(j,1) = q(j,1) + lambdabio*DOP_j
            q(j,2) = q(j,2) - lambdabio*DOP_j
        end do
        
        ! scale
        do j = 1, nz
            q(j,1) = q(j,1)*dtbio
            q(j,2) = q(j,2)*dtbio
        end do

      end
