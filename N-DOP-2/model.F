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
C   Model for Phosphate, PO4 and dissolved organic phosphorous, DOP
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

        call bgc(
     &      t, dt, n, nz, nbc, bc, ndc, dc, y, q, m, u)

      end

      subroutine bgc(t, dt, ntracer, nlayer, nbc, bc,
     &      ndc, dc, y, ybgc, nparam, u)
            implicit none
            ! input variables
            real*8  :: t, dt
            integer :: ntracer, nlayer, nbc, ndc, nparam
            real*8  :: bc(nbc), dc(nlayer, ndc)
            real*8  :: y(nlayer, ntracer);
            real*8  :: ybgc(nlayer, ntracer), u(nparam)
            ! dummy
            real*8  :: dayfrac, wind
            real*8  :: salt(nlayer), temp(nlayer)
            ! to compute
            real*8  :: insol, stau
            real*8  :: d(nlayer)
            ! work
            integer :: ilayer
            d(1) = dc(1,1)
            do ilayer = 2, nlayer
                d(ilayer) = d(ilayer-1) + dc(ilayer,1)
            end do
            call insolation(1, t*360.d0*86400.d0, bc(1), insol, stau)
            call bgcmodel(t, dt, bc(1), insol, dayfrac, bc(2), wind,
     &          ntracer, nlayer, dc(:,1), d, salt, temp, y, ybgc,
     &          nparam, u)
      end

      subroutine bgcmodel(t, dt, lat, insol, dayfrac, ice, wind,
     &  ntracer, nlayer, h, d, salt, temp, y, ybgc, nparam, u)
            implicit none
            ! input variables
            real*8  :: t, dt, lat, insol, dayfrac, ice, wind
            integer :: ntracer, nlayer, nparam
            real*8  :: h(0:nlayer-1), d(0:nlayer-1)
            real*8  :: salt(0:nlayer-1), temp(0:nlayer-1)
            real*8  :: y(0:nlayer-1, 0:ntracer-1)
            real*8  :: ybgc(0:nlayer-1, 0:ntracer-1)
            real*8  :: u(0:nparam-1)
            ! integer value parameters
            integer :: nbiostep    = 1
            integer :: jeuphotic   = 2
            ! real value parameters
            real*8  :: lambda                       ! remineralization rate            : 1 / yr
            real*8  :: alpha                        ! maximum export rate              : uM / yr
            real*8  :: nu                           ! fraction of DOP
            real*8  :: K_PO4                        ! half saturation constant: PO4    : uM
            real*8  :: K_insol                      ! half saturation constant: light  : W / m^2
            real*8  :: K_H2O                        ! attenuationOfWater
            real*8  :: b                            ! Martin exponent
            real*8  :: biofrac      =  0.4d0        ! fraction of light spectrum that is used
            ! work vars
            real*8  :: ywork(0:nlayer-1, 0:ntracer-1)
            real*8  :: dtbio, lambdabio, alphabio
            integer :: ibiostep, j, k, minlayer
            real*8  :: PO4_j, DOP_j, insol_j, gamma_j, nu_bar
            ! set parameters
            lambda  = u(0)
            alpha   = u(1)
            nu      = u(2)
            K_PO4   = u(3)
            K_insol = u(4)
            K_H2O   = u(5)
            b       = u(6)
            ! prepare model
            insol_j     = 1.0d0
            nu_bar      = (1.0d0-nu)
            minlayer    = nlayer
            !minlayer    = min(nlayer, jeuphotic);
            ! fraction of parameters
            dtbio       = (dt*360.0d0)/nbiostep
            lambdabio   = lambda/360.0d0
            alphabio    = alpha/360.0d0
            insol       = (1.0d0-ice)*biofrac*insol
            ! bio steps per ocean step
            do ibiostep = 0, nbiostep-1
                ! zero work memory
                do j = 0, nlayer-1
                    ywork(j, 0) = 0.0d0
                    ywork(j, 1) = 0.0d0
                end do
                ! euphotic layers
                do j = 0, minlayer-1
                    ! attenaution of water
                    if (j == 0) then
                        ! first layer
                        insol_j = exp(-d(j)*0.5d0*K_H2O)*insol
                    else
                        ! other layers
                        insol_j = exp(-(d(j-1)+0.5d0*h(j))*K_H2O)*insol
                    end if
                    ! take previous steps into account
                    PO4_j   = y(j, 0)+ybgc(j, 0)
                    ! bio uptake amount
                    gamma_j =
     &  alphabio*PO4_j*insol_j/(PO4_j+K_PO4)/(insol_j+K_insol)
                    ! store in work variable
                    ! consumption
                    ywork(j, 0) = ywork(j, 0)-gamma_j
                    ywork(j, 1) = ywork(j, 1)+nu*gamma_j
                    ! two cases
                    if (j == nlayer-1) then
                        ! remineralization of last layer
                        ywork(j, 0) = ywork(j, 0)+nu_bar*gamma_j
                    else
                        ! export to layers below
                        do k = j+1, nlayer-1
                            ! last layer
                            if (k == nlayer-1) then
                        ywork(k, 0) =
     &  ywork(k, 0)+gamma_j*nu_bar*h(j)*(d(k-1)/d(j))**(-b)/h(k)
                            ! layers in between
                            else
                        ywork(k, 0) =
     &  ywork(k, 0)+gamma_j*nu_bar*h(j)
     &  *((d(k-1)/d(j))**(-b)-(d(k)/d(j))**(-b))/h(k)
                            end if
                        end do
                    end if
                end do
                ! remineralization
                do j = 0, nlayer-1
                    ! take previous steps into account
                    DOP_j = y(j, 1)+ybgc(j, 1)
                    ! update
                    ywork(j, 0) =
     &  ywork(j, 0)+lambdabio*max(DOP_j, 0.0d0)
                    ywork(j, 1) =
     &  ywork(j, 1)-lambdabio*max(DOP_j, 0.0d0)
                end do
                ! scale with bio time step
                ! and add to BGC part
                do j = 0, nlayer-1
                    ! scale
                    ywork(j, 0) = ywork(j, 0)*dtbio
                    ywork(j, 1) = ywork(j, 1)*dtbio
                    ! add
                    ybgc (j, 0) = ybgc(j, 0)+ywork(j, 0)
                    ybgc (j, 1) = ybgc(j, 1)+ywork(j, 1)
                end do
            end do
            ! bio steps per ocean step
      end
