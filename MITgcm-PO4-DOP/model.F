!
! Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D
! Copyright (C) 2012  Jaroslaw Piwonski, CAU, jpi@informatik.uni-kiel.de
! Copyright (C) 2012  Thomas Slawig, CAU, ts@informatik.uni-kiel.de
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
!   Model for PO4, DOP from MITgcm
!

#include "insol.F"
#include "phos_flux.F"
#include "bio_export.F"
#include "dic_biotic_forcing.F"

#include "DIC_OPTIONS.h"

      subroutine metos3dbgcinit(n, nz, m, nbc, ndc, dt,
     &  q, t, y, u, bc, dc)
        implicit none
        ! input variables
        integer :: n, nz, m, nbc, ndc
        real*8  :: dt, q(nz, n), t, y(nz, n), u(m)
        real*8  :: bc(nbc), dc(nz, ndc)

#include "SIZE.h"
#include "GRID.h"
#include "DIC_VARS.h"
#include "PARAMS.h"
#include "PTRACERS_PARAMS.h"

        ! unused
        R_CP            = 117. _d 0                  ! bgcinit
        R_NP            = 16. _d 0                   ! bgcinit
        ! used
        parfrac         = 0.4 _d 0                   ! bgcinit
        PTRACERS_dTLev  = dt*360.d0                  ! bgcinit
        ! depths and heights
        ! make sure depths are: dc(:,1)
        !          heights are: dc(:,2)
        rF(1)           = 0.d0                       ! bgcinit
        rF(2:nz+1)      = - dc(1:nz,1)
        drF(1:nz)       = dc(1:nz,2)
        recip_drF(1:nz) = 1.d0/drF(1:nz)             ! bgcinit

      end

      subroutine metos3dbgc(n, nz, m, nbc, ndc, dt,
     &  q, t, y, u, bc, dc)
        implicit none
        ! input variables
        integer :: n, nz, m, nbc, ndc
        real*8  :: dt, q(nz, n), t, y(nz, n), u(m)
        real*8  :: bc(nbc), dc(nz, ndc)

#include "SIZE.h"
#include "GRID.h"
#include "DIC_VARS.h"
#include "PARAMS.h"
#include "PTRACERS_PARAMS.h"

        _RL  PTR_DIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
        _RL  PTR_ALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
        _RL  PTR_PO4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
        _RL  PTR_DOP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)

        INTEGER myIter
        real*8 myTime
        INTEGER myThid
        integer bi, bj, imin, imax, jmin, jmax

        integer icp

        ! set variables
        ! unused
        bi = 1
        bj = 1
        imin = 1
        imax = 1
        jmin = 1
        jmax = 1
        myIter = 0
        myThid = 0

        ! point in time
        myTime = t*360.d0*86400.d0

        ! set parameters
        KDOPRemin                                   = u(1)/360.d0
        alpha(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)  = u(2)/360.d0
        DOPfraction                                 = u(3)
        KPO4                                        = u(4)
        lit0                                        = u(5)
        k0                                          = u(6)
        KRemin                                      = u(7)

        ! unused
        PTR_DIC = 0.d0
        PTR_ALK = 0.d0

        ! set nlayer/nlevel
        nlev = nz

        ! set boundary data
        ! make sure latitude is: bc(1)
        !          ice cover is: bc(2)
        yC   = bc(1)
        FIce = bc(2)

        ! set MITgcm masks
        maskC = 0.d0
        maskC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,1:nlev,nSx,nSy) = 1.d0

        hFacC = 0.d0
        hFacC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,1:nlev,nSx,nSy) = 1.d0
        recip_hFacC = 0.d0
        recip_hFacC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,1:nlev,nSx,nSy) = 1.d0

        ! set PO4, DOP
        PTR_PO4 = 0.d0
        PTR_DOP = 0.d0
        do icp = 1, nz
            PTR_PO4(1,1,icp) = y(icp,1)
            PTR_DOP(1,1,icp) = y(icp,2)
        enddo

        ! call model
        call DIC_BIOTIC_FORCING( PTR_DIC, PTR_ALK, PTR_PO4,
     &                            PTR_DOP,
     &                            bi,bj,imin,imax,jmin,jmax,
     &                             myIter,myTime,myThid)

        ! set ybgc
        do icp = 1, nz
            q(icp, 1) = PTR_PO4(1,1,icp) - y(icp, 1)
            q(icp, 2) = PTR_DOP(1,1,icp) - y(icp, 2)
        enddo

      end

        subroutine metos3dbgcfinal(n, nz, m, nbc, ndc, dt,
     &  q, t, y, u, bc, dc)
            implicit none
            ! input variables
            integer :: n, nz, m, nbc, ndc
            real*8  :: dt, q(nz, n), t, y(nz, n), u(m)
            real*8  :: bc(nbc), dc(nz, ndc)
        end

        SUBROUTINE ALK_SURFFORCING( PTR_ALK , GALK,
     I           bi,bj,imin,imax,jmin,jmax,
     I           myIter,myTime,myThid)

#include "SIZE.h"

        INTEGER myIter
        real*8 myTime
        INTEGER myThid

        real*8 PTR_ALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
        real*8 GALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

        integer bi, bj, imin, imax, jmin, jmax


        end subroutine

        SUBROUTINE CAR_FLUX( CAR_S, cflux,
     I           bi,bj,imin,imax,jmin,jmax,
     I           myIter,myTime,myThid)

#include "SIZE.h"

        INTEGER myIter
        real*8 myTime
        INTEGER myThid

        real*8 PTR_ALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)

        integer bi, bj, imin, imax, jmin, jmax

        end subroutine

        SUBROUTINE DIC_SURFFORCING( PTR_CO2 , PTR_ALK, PTR_PO4, GDC,
     I           bi,bj,imin,imax,jmin,jmax,
     I           myIter,myTime,myThid)

#include "SIZE.h"

        INTEGER myIter
        real*8 myTime
        INTEGER myThid

        real*8 PTR_CO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
        real*8 PTR_ALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
        real*8 PTR_PO4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)

        real*8 GDC(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

        integer bi, bj, imin, imax, jmin, jmax

        end subroutine

