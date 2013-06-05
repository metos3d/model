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
C   insolation.F
C

CBOP
C !ROUTINE: INSOL

C !INTERFACE: ==========================================================
      SUBROUTINE insolation(N,Time,Yc,sfac,stau)

C !DESCRIPTION:
C find shortwave radiation as function of date and latitude
c based on paltridge and parson
c modified by SPK from MITGCM

C !USES: ===============================================================
      IMPLICIT NONE
C     === Global variables ===

C !INPUT PARAMETERS: ===================================================
C time                 :: current time
       integer N
       real*8 time,Yc(N)

C !OUPUT PARAMETERS: ===================================================
       real*8 sfac(N),stau(N)

C !LOCAL VARIABLES: ====================================================
       real*8  solar, albedo, par
       real*8  dayfrac, yday, delta
       real*8  lat, sun1, dayhrs
       real*8  cosz, frac, fluxi, fracmin
       integer j
CEOP

c
      solar = 1360.   !solar constant
      albedo = 0.6    !planetary albedo
C      par = 0.4       !photosynthetically reactive frac

c
c find day (****NOTE for year starting in winter*****)
        dayfrac=mod(Time,360.d0 *86400.d0)
     &                    /(360.d0 *86400.d 0)          !fraction of year
        yday = 2.0*3.1416*dayfrac                         !convert to radians
        delta = (0.006918 - (0.399912*cos(yday))      !cosine zenith angle
     &          +(0.070257*sin(yday))                 !(paltridge+platt)
     &          -(0.006758*cos(2.0*yday))
     &          +(0.000907*sin(2.0*yday))
     &          -(0.002697*cos(3.0*yday))
     &          +(0.001480*sin(3.0*yday)) )
       do j=1,N
c latitude in radians
          lat=YC(j)/180.d0*3.1416
          sun1 = -sin(delta)/cos(delta) * sin(lat)/cos(lat)
          if (sun1.le.-0.999) sun1=-0.999
          if (sun1.ge. 0.999) sun1= 0.999
          dayhrs = abs(acos(sun1))
          cosz = ( sin(delta)*sin(lat)+              !average zenith angle
     &            (cos(delta)*cos(lat)*sin(dayhrs)/dayhrs) )
          if (cosz.le.0.005) cosz=0.005
          frac = dayhrs/3.1416               !fraction of daylight in day
c daily average photosynthetically active solar radiation just below surface
          fluxi = solar*(1.0-albedo)*cosz*frac
c
c convert to sfac
          if (fluxi.gt.0.0) sfac(j)=fluxi
c very large for polar night
          if (fluxi.lt.0.00001) sfac(j)=0.00001
c daylength; ensure that it lies between 0 and 1 (may be slightly
c out of this range in high latitudes)
          fracmin = MIN(frac,1.0)
          stau(j) = MAX(fracmin,0.0)

       enddo !j
      return
      end
c
C==========================================================================
