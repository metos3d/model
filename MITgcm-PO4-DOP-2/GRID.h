
        COMMON /GRID_RS/
     &      maskC, hFacC, drF, yC, fCori, rF, recip_hFacC, recip_drF

        _RS maskC          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
        _RS hFacC          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
        _RS drF            (Nr)
        _RS yC             (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
        _RS fCori          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
        _RS rF             (Nr+1)
        _RS recip_hFacC    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
        _RS recip_drF      (Nr)
