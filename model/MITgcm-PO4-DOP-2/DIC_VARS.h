
        COMMON /CARBON_NEEDS/
     &      FIce

        _RL  FIce(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

        COMMON /BIOTIC_NEEDS/
     &      DOPfraction, KDOPremin, R_cp, R_NP, rain_ratio, k0, Kpo4,
     &      lit0, nlev, QSW_underice, alpha, parfrac, KRemin

        _RL DOPfraction
        _RL KDOPremin
        _RL R_cp
        _RL R_NP
        _RL rain_ratio(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
        _RL k0
        _RL Kpo4
        _RL lit0
        INTEGER nlev
        LOGICAL QSW_underice
        _RL alpha(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
        _RL parfrac
        _RL KRemin

c        DOPfraction = 0.67 _d 0
c        KDOPRemin   = 1. _d 0/(6. _d 0*30. _d 0*86400. _d 0)
c        R_CP        = 117. _d 0
c        R_NP        = 16. _d 0

c        rain_ratio
c        rainRatioUniform = 7. _d -2

c        k0          = 0.02 _d 0
c        KPO4        = 5. _d -4
c        lit0        = 30. _d 0

c        alphaUniform     = 2. _d -3/(360. _d 0 * 86400. _d 0)
c        parfrac     = 0.4 _d 0
c        KRemin      = 0.9 _d 0
