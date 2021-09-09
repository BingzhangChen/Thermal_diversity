PROGRAM MAIN
USE TOPT_MOD
IMPLICIT NONE

!Seasonal Amplitude of temperature: (Tmax - Tmin)/2
INTEGER, PARAMETER :: NAm     = 2
REAL,    PARAMETER :: Am(NAm) = [3., 6.]


!A vector of mean temperature (three)
integer, parameter :: NavgT = 2
REAL    :: avgtemp(NavgT)   = [6., 24.]
REAL    :: cur_day ! Current simulation time (day)
INTEGER :: j,k, NI,v
REAL    :: t1, t2, Ew
character(len=10), parameter :: format_string = "(A3,I0)"

!Filename
character(LEN=20) :: fname = 'Mod.out'

!Start 
CALL CPU_TIME(T1)

! Number of integrations:
NI  = NDAYS*INT(1.d0/dtdays)
 
!Initialise phytoplankton label
DO i = 1, NPHY
   WRITE(PHYLAB(I), FORMAT_STRING) 'PHY',I
   WRITE(muMLAB(I), FORMAT_STRING) 'muM',I
ENDDO

DO v = 1, NAm
   DO i = 1, NavgT
     do k = 1, NaG
   
      !Initialise variables:
      NO3  = 1.; PHY(:)=1D-3; 
      PHYt = sum(PHY(:))
      ZOO  = PHYt*0.1
    
      !Create file name:
      WRITE(fname,100) 'Mod', int(avgtemp(i)),'_',k,'_',v,'.out'
      
      !Create file for saving results
      OPEN (UNIT=10, file = fname, status = 'replace')
      WRITE(10, 101) 'Day', 'Temp', 'NO3', 'PHY', 'ZOO', 'avgZ', &
                     'sigmaZ', 'mumax', PHYLAB, muMLAB, "Ew" 
      
      !Calculate temperature
      cur_day   = 0.
      temp      = avgtemp(i) - Am(v)*cosine(cur_day)
      mumax_avg = 0.
      mumax(:)  = 0.
   
      call moments !Calculate trait variance
   
      !Record initial conditions:
      WRITE(10, 102) 0., temp, NO3, PHYt, ZOO, avgZ, sigmaZ, 0., PHY, mumax, 0.65
      
      DO j = 1, NI
      
        !Current time (in days)
        cur_day = dtdays*dble(j)
        temp    = avgtemp(i) - Am(v)*cosine(cur_day)
   
        !Save data per day:
        IF (mod(j,int(1./dtdays)) .eq. 0) THEN
          call moments
          call SENSITIVITY(Ew)
          WRITE(10, 102) cur_day, temp, NO3, PHYt, ZOO, avgZ, sigmaZ,mumax_avg, PHY,mumax, Ew 
        ENDIF
      
        call biology(alphaG(k)) 
      ENDDO
      
      !Close file
      CLOSE(10)
     enddo
   ENDDO
ENDDO

CALL CPU_TIME(T2)
WRITE(6, '("TIME = ",F8.3," MINS.")') (T2-T1)/60.0

100 format(A3, I0,2(A1,I0),A4) 
101 format(6x,440(A7, 6x))
102 format(F12.1, 1x, 439(1pe10.2, 1x))   
END PROGRAM
