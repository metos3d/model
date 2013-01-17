
        Real*8 PI
        PARAMETER ( PI    = 3.14159265358979323844D0   )
        Real*8 deg2rad
        PARAMETER ( deg2rad = 2.D0*PI/360.D0           )

        COMMON /PARM_R/
     &      omega

        _RL omega

        COMMON /PARM_L/
     &      rotateGrid, usingCartesianGrid, usingCurvilinearGrid,
     &      usingCylindricalGrid

        LOGICAL rotateGrid
        LOGICAL usingCartesianGrid
        LOGICAL usingCurvilinearGrid
        LOGICAL usingCylindricalGrid
