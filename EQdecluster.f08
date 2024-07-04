!__________________________________________________________________________________________________________________________________
!
      PROGRAM EQDECLUSTER
!__________________________________________________________________________________________________________________________________
!
!     Copyright (c) 2024 by Álvaro González
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!     Contact info: alvaro@geonaut.eu
!__________________________________________________________________________________________________________________________________
!
      IMPLICIT NONE
!
!     ----------
!     PARAMETRES
!     ----------
!
      INTEGER, PARAMETER    :: Kind   = SELECTED_REAL_KIND(15) ! Kind value of a real data type with decimal precision greater of at least 15 digits
      REAL(Kind), PARAMETER :: Radius = 6371.007 ! Authalic Earth radius (km)
!
!     ---------
!     FUNCTIONS 
!     ---------
!
      REAL(Kind) DISTANCE  ! Great-circle distance between two points on the sphere
      REAL(Kind) JD        ! Julian date (days)
!
!     -----------------
!     GENERIC VARIABLES
!     -----------------
! 
      INTEGER(16)   i         ! Counting index
      INTEGER(16)   j         ! Counting index
!
!     ----------------------------------------------------
!     VARIABLES REGARDING ALL EARTHQUAKES IN CATALOGUE.DAT
!     ----------------------------------------------------
!
      INTEGER(16)             eTotal         ! Total number of earthquakes in catalogue.dat
      INTEGER, ALLOCATABLE    :: y(:)        ! Year
      INTEGER, ALLOCATABLE    :: mon(:)      ! Month
      INTEGER, ALLOCATABLE    :: d(:)        ! Day
      INTEGER, ALLOCATABLE    :: h(:)        ! Hour
      INTEGER, ALLOCATABLE    :: minut(:)    ! Minute
      REAL(Kind), ALLOCATABLE :: s(:)        ! Second
      REAL(Kind), ALLOCATABLE :: mag(:)      ! Magnitude
      REAL(Kind), ALLOCATABLE :: julian(:)   ! Julian date (days) of origin time
      REAL(Kind), ALLOCATABLE :: lon(:)      ! Epicentral longitude
      REAL(Kind), ALLOCATABLE :: lat(:)      ! Epicentral latitude
      INTEGER, ALLOCATABLE    :: GK(:)       ! Gardner-Knopoff (1=mainshock, 0=aftershock or foreshock)
      INTEGER, ALLOCATABLE    :: GR(:)       ! Grünthal (1=mainshock, 0=aftershock or foreshock)
      INTEGER, ALLOCATABLE    :: PE(:)       ! Peláez et al. (1=mainshock, 0=aftershock or foreshock)
      INTEGER, ALLOCATABLE    :: IGN(:)      ! IGN (1=mainshock, 0=aftershock or foreshock)
!
!     --------------------------------------------------------
!     VARIABLES FOR TIME & SPACE DISTANCES BETWEEN EARTHQUAKES
!     --------------------------------------------------------
!
      REAL(Kind) interval ! Interval (days) between origin times
      REAL(Kind) dist     ! Great-circle distance from each epicentre to the next one
!
!     ----------------------------------------------------------------------------
!     FIRST READING OF CATALOGUE.DAT - COUNT THE TOTAL NUMBER OF EARTHQUAKES
!     ----------------------------------------------------------------------------
!     In order to allocate the minimum possible memory space, catalogue.dat is read two times in total.
!
      OPEN (22, FILE="catalogue.dat", STATUS="OLD")
      eTotal=0
      DO
        READ(22,*,END=1000)
        eTotal=eTotal+1
      END DO
1000  CONTINUE
      PRINT *
      IF(eTotal == 0) THEN
        PRINT *, "Error: The file catalogue.dat contains no data."
        STOP
      ELSE
        PRINT *, "Earthquakes in catalogue.dat: ", eTotal
      END IF  
!
!     --------------------------------------------------------------------
!     ALLOCATE MEMORY FOR RECORDING EARTHQUAKE DATA, AND INITIALIZE VALUES
!     --------------------------------------------------------------------
!
      ALLOCATE (y(eTotal))
      ALLOCATE (mon(eTotal))
      ALLOCATE (d(eTotal))
      ALLOCATE (h(eTotal))
      ALLOCATE (minut(eTotal))
      ALLOCATE (s(eTotal))
      ALLOCATE (julian(eTotal))
      ALLOCATE (lon(eTotal))
      ALLOCATE (lat(eTotal))
      ALLOCATE (mag(eTotal))
      ALLOCATE (GK(eTotal))
      ALLOCATE (GR(eTotal))
      ALLOCATE (PE(eTotal))
      ALLOCATE (IGN(eTotal))
      GK=1  ! (All are mainshocks by default)
      GR=1  ! (All are mainshocks by default)
      PE=1  ! (All are mainshocks by default)
      IGN=1 ! (All are mainshocks by default) 
!
!     --------------------------------------------------------
!     SECOND READING OF CATALOGUE.DAT - RECORD EARTHQUAKE DATA
!     --------------------------------------------------------
!
      REWIND(22)   ! Start again from the beginning of catalogue.dat
      DO i=1, eTotal
        READ(22,*)  y(i), mon(i), d(i), h(i), minut(i), s(i), lon(i), lat(i), mag(i)
        julian(i)  = JD(y(i), mon(i), d(i), h(i), minut(i), s(i))
      END DO
!
!     ---------------------------------------------------------------
!     MAIN PROGRAM LOOP (EARTHQUAKE BY EARTHQUAKE, IN TEMPORAL ORDER)
!     ---------------------------------------------------------------
!     Considers an earthquake in turn (i-th) and each of the next ones (j-th).
!
      PRINT *
      PRINT *, "Please wait while processing earthquake data..."
      PRINT *
      DO i=1, (eTotal-1)
        DO j=(i+1), eTotal
          interval=julian(j)-julian(i)
          dist=DISTANCE(lat(i),lon(i),lat(j),lon(j),Radius)
!
!         Gardner & Knopoff (1974)
! 
          IF(GK(i).EQ.1) THEN ! Do not consider aftershocks of aftershocks
            IF (mag(i).GE.6.5) THEN
              IF (interval.LE.(10.0**(0.032*mag(i)+2.7389))) THEN
                IF (dist.LE.(10.0**(0.1238*mag(i)+0.983))) THEN
                  IF (mag(i).GE.mag(j)) THEN ! Second earthquake is an aftershock
                    GK(j)=0
                  ELSE ! First earthquake is a foreshock
                    GK(i)=0
                  END IF
                END IF
              END IF
            ELSE
              IF (interval.LE.(10.0**(0.5409*mag(i)-0.547))) THEN
                IF (dist.LE.(10.0**(0.1238*mag(i)+0.983))) THEN
                  IF (mag(i).GE.mag(j)) THEN ! Second earthquake is an aftershock
                    GK(j)=0
                  ELSE ! First earthquake is a foreshock
                    GK(i)=0
                  END IF
                END IF
              END IF
            END IF
          END IF
!
!         Grünthal (1985)
!
          IF(GR(i).EQ.1) THEN ! Do not consider aftershocks of aftershocks
            IF (mag(i).LT.6.65) THEN
              IF (interval.LE.ABS(EXP(-3.95+SQRT(0.62+17.32*mag(i))))) THEN
                IF (dist.LE.(EXP(1.77+SQRT(0.037+1.02*mag(i))))) THEN
                  IF (mag(i).GE.mag(j)) THEN ! Second earthquake is an aftershock
                    GR(j)=0
                  ELSE ! First earthquake is a foreshock
                    GR(i)=0
                  END IF
                END IF
              END IF   
            ELSE
              IF (interval.LE.(10.0**(2.8+0.024*mag(i)))) THEN
                IF (dist.LE.(EXP(1.77+SQRT(0.037+1.02*mag(i))))) THEN
                  IF (mag(i).GE.mag(j)) THEN  ! Second earthquake is an aftershock
                    GR(j)=0
                  ELSE ! First earthquake is a foreshock
                    GR(i)=0
                  END IF
                END IF
              END IF
            END IF
          END IF
!
!         Peláez et al. (2007)
!
          IF (PE(i).EQ.1) THEN  ! Do not consider aftershocks of aftershocks
            IF (interval.LE.10.**(0.3908485*mag(i)-0.17255)) THEN
              IF (dist.LE.10.**(0.139794*mag(i)+0.881648)) THEN
                IF (mag(i).GE.mag(j)) THEN ! Second earthquake is an aftershock
                  PE(j)=0
                ELSE ! First earthquake is a foreshock
                  PE(i)=0
                END IF
              END IF
            END IF
          END IF
!
!         IGN (2013)
!
          IF (IGN(i).EQ.1) THEN  ! Do not consider aftershocks of aftershocks
            IF (mag(i).GE.5.9) THEN
              IF (interval.LE.(10.0**(0.04*mag(i)+2.6))) THEN
                IF (dist.LE.(10.0**(0.2*mag(i)+0.4))) THEN
                  IF (mag(i).GE.mag(j)) THEN ! Second earthquake is an aftershock
                    IGN(j)=0
                  ELSE ! First earthquake is a foreshock
                    IGN(i)=0
                  END IF
                END IF
              END IF   
            ELSE
              IF (interval.LE.(10.0**(0.6*mag(i)-0.7))) THEN
                IF (dist.LE.(10.0**(0.2*mag(i)+0.4))) THEN
                  IF (mag(i).GE.mag(j)) THEN ! Second earthquake is an aftershock
                    IGN(j)=0
                  ELSE ! First earthquake is a foreshock
                    IGN(i)=0
                  END IF
                END IF
              END IF
            END IF
          END IF
        END DO
      END DO
      OPEN (24, FILE = "results.dat",   STATUS="UNKNOWN")
      WRITE (24,*) "Gardner-Knopoff", "Gruenthal", "Pelaez", "IGN"
      DO i=1, eTotal
        WRITE (24,*) GK(i), GR(i), PE(i), IGN(i)
      END DO
      CLOSE (24)
      PRINT *
      PRINT *, "Calculation finished correctly."
      END PROGRAM
!
!__________________________________________________________________________________________________________________________________
!
      FUNCTION DISTANCE(lat1,lon1,lat2,lon2,Radius)
!
!     Calculates the shortest distance between two points on the surface of
!     a sphere, measured along a great circle.
!     Uses the spherical law of cosines.
!     ______________________________________________________________________
!
!     Copyright (c) 2010-2024 Alvaro Gonzalez
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!     Contact info: alvaro@geonaut.eu
!     ______________________________________________________________________
!   
      IMPLICIT NONE
      INTEGER, PARAMETER :: Kind=SELECTED_REAL_KIND(15) ! Precision control
      REAL(Kind), PARAMETER :: Pi = 3.1415926535897932_Kind ! pi mathematical constant
      REAL(Kind) :: radius    ! Radius of the sphere
      REAL(Kind) :: lat1      ! Latitude of first point (degrees)
      REAL(Kind) :: lon1      ! Longitude of first point (degrees)
      REAL(Kind) :: lat2      ! Latitude of second point (degrees)
      REAL(Kind) :: lon2      ! Longitude of second point (degrees)
      REAL(Kind) :: DISTANCE  ! Output value of the function
      REAL(Kind) :: cosine    ! Of the angle between points and the sphere centre

      IF((lat1 == lat2).AND.(lon1 == lon2)) THEN ! Points are identical
        DISTANCE=0
      ELSE                                       ! Points are different
        cosine = COS(lat1*Pi/180.)*COS(lat2*Pi/180.) * COS((lon1-lon2)*Pi/180.) + SIN(lat1*Pi/180.)*SIN(lat2*Pi/180.)
        IF (cosine > 1D+0) cosine=1D+0   ! Prevent finite-precission error
        IF (cosine < -1D+0) cosine=-1D+0 ! Prevent finite-precission error
        DISTANCE=radius*ACOS(cosine)
      END IF
      RETURN
      END FUNCTION DISTANCE
!
!__________________________________________________________________________________________________________________________________
!
      FUNCTION JD(year,month,day,hour,minute,second)
!
!     Computes the Julian Date (JD),
!      given a Gregorian calendar date (year, month, day)
!      and Universal Time -UT- (hour, minute, second).
!
!     Uses the algorithm by:
!      Fliegel, H. F. & van Flandern, T. C. (1968):
!      A machine algorithm for processing calendar dates.
!      Communications of the ACM, Vol. 11, No. 10, p. 657.
!
!     Minor correction by Alvaro Gonzalez: 0.5 is substracted to the date,
!     because 0h UT corresponds to a Julian date fraction of 0.5.
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: Kind=SELECTED_REAL_KIND(15) ! Precision control
      INTEGER year, month, day, hour, minute
      REAL(Kind) second
      REAL(Kind) JD
      JD = day-32075+1461*(year+4800+(month-14)/12)/4 + 367*(month-2-(month-14)/12*12)/12
      JD = JD - 3*((year+4900+(month-14)/12)/100)/4 + hour/24D+0 + minute/1440D+0 + second/86400D+0 - 0.5D+0
      RETURN
      END FUNCTION JD
