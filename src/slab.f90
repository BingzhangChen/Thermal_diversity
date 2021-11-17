MODULE TOPT_MOD

INTEGER, PARAMETER :: NPHY = 41
REAL               :: ZOO=0., NO3=0., PHY(NPHY) = 0.
REAL               :: avgZ = 0., sigmaZ = 0.

!Total phytoplankton biomass
REAL               :: PHYt = 0.

!Maximal growth rate of each species
REAL               :: mumax(NPHY) = 0.
INTEGER        :: i

!Optimal temperature in celcius
REAL,    PARAMETER :: Topt(NPHY) = [(-5. + dble(i), i = 1, NPHY)]

!Labels for phytoplankton
CHARACTER(LEN=5)   :: PHYLAB(NPHY), muMLAB(NPHY)

!Number of alphaG values
REAL,        PARAMETER :: pi = 3.141592653589798 

! Temperature
REAL :: temp = 0., mumax_avg=0.d0
REAL,    PARAMETER :: par_= 10000.  !Saturating light

! Time step
integer, parameter :: d_per_y = 365   !How many days a year
integer, parameter :: s_per_d = 86400 !How many seconds a day

!Number of days running:
integer, parameter :: NDAYS   = d_per_y*5

!Timestep (unit: day)
REAL,    PARAMETER :: dtdays  = 300.d0/DBLE(S_PER_D)
REAL,    PARAMETER :: KN0     = 0.5       !Half-saturation constant for phytoplankton growth (mM N) normalized to 15 C      
REAL,    PARAMETER :: aI0       = 0.1
CONTAINS

REAL FUNCTION cosine(day) result(y) 
IMPLICIT NONE
real, intent(in) :: day
real, parameter  :: offset = 60.d0
  y   = cos(2.d0 * pi * (day-offset) / dble(d_per_y) )
END FUNCTION cosine

SUBROUTINE BIOLOGY(ALPHAG_, E_KN_)
IMPLICIT NONE
REAL, INTENT(IN):: ALPHAG_, E_KN_
REAL, PARAMETER :: RDN  = 0.1    !Regeneration rate from DET to NO3
REAL, PARAMETER :: gmax = 1.d0  !Maximal zooplankton ingestion rate
REAL, PARAMETER :: GGE  = 0.3    !Gross growth efficiency

INTEGER         :: j
REAL            :: total_uptake,INGES,RES, Zmort,tf_z,gbar
REAL            :: pp_ZP, pp_NZ
REAL            :: mu,  PHYtot2
REAL            :: Graz         = 0. 
REAL            :: KN           = 0. 
REAL            :: mu_avg    = 0. 
REAL, PARAMETER :: Kp  = 0.5d0
REAL, PARAMETER :: mz  = 0.15d0
REAL, PARAMETER :: Ez  = 0.65d0

PHYt    = sum(PHY(:))
PHYtot2 = 0d0   !For Kill-the-winner grazing
do j = 1, NPHY
  PHYtot2 = PHYtot2 + PHY(j)**alphaG_
enddo

total_uptake = 0.d0

! Mean growth rate of the community
mu_avg       = 0.

! Mean µmax of the community
mumax_avg    = 0.

!temperature coefficient of zooplankton grazing
tf_z  = exp(TK(temp)*Ez)

! The grazing dependence on total prey
gbar  = PHYt**2/(PHYt**2 + Kp**2)

!Zooplankton per capita total ingestion rate
INGES = tf_z*gmax*gbar

!Zooplankton excretion rate (-> DIN)
RES = INGES*(1.0d0 - GGE)

Zmort = ZOO**2*mz*tf_z  !Mortality term for ZOO

do j = 1, NPHY
               
   call GROWTH(temp, NO3, par_,  E_KN_, Topt(j),  mu,  mumax(j))

   mu_avg   = mu_avg    + mu*PHY(j)
   mumax_avg= mumax_avg + mumax(j)*PHY(j)

   !Maintain diversity using kill-the-winner
   !Calculate the mortality of each phyto species by ZOO
   Graz = (INGES*ZOO/PHYtot2)*PHY(j)**alphaG_

   !Compute total phyto. uptake 
   total_uptake = total_uptake + mu*PHY(j)

   PHY(j) = PHY(j) + dtdays*(mu*PHY(j) - Graz)
enddo
mu_avg        =    mu_avg/PHYt           !Phytoplankton community average growth rate
mumax_avg = mumax_avg/PHYt      !Phytoplankton community averaged maximal growth rate

! For production/destruction matrix:
pp_NZ = ZOO*RES + Zmort        
pp_ZP = ZOO*INGES      
NO3   = NO3 + (pp_NZ-total_uptake)*dtdays
ZOO   = ZOO + dtdays*(pp_ZP - pp_NZ)
END SUBROUTINE BIOLOGY

SUBROUTINE SENSITIVITY(E_KN_, Einst)
IMPLICIT NONE
real, intent(in)    :: E_KN_
real, intent(out)  :: Einst
!Calculate the temperature sensitivity of the community for each day
!without considering the changes of community composition

!Suppose we already have mumax_avg
!Then we calculate mumax_avg at (temp -1)
real    :: t_, mumax_, mu_, mumax_avg_
integer :: j

t_ = temp - 1.d0
mumax_avg_ = 0.d0

do j = 1, NPHY
   call GROWTH(t_, NO3, par_,  E_KN_, Topt(j),  mu_,  mumax_)
   mumax_avg_ = mumax_avg_ + mumax_ * PHY(j)
enddo
mumax_avg_ = mumax_avg_/PHYt
Einst  = log(mumax_avg/mumax_avg_)/(TK(temp)-TK(t_))
END SUBROUTINE SENSITIVITY

!Eq. (4) of the main manuscript
REAL FUNCTION JOHNSON(tC, mumax, Ea, Eh, Topt_) RESULT(y)
IMPLICIT NONE
!Both tC and Topt_ are in ºC
real,   intent(in)  :: tC, mumax, Ea, Eh, Topt_
real,   parameter   :: kb   = 8.62D-5
real,   parameter   :: T0   = 273.15D0
real,   parameter   :: Tref = 15D0
real                       :: ED, x, theta

   ED = Eh-Ea
   x  = TK(TC)
   theta = TK(Topt_)
   y = mumax*Eh/ED*exp(Ea*(x-theta))/(1.D0+Ea/ED*exp(Eh*(x-theta)))   
return
END FUNCTION JOHNSON

PURE REAL FUNCTION TK(TC)
IMPLICIT NONE
!DESCRIPTION:
!The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
! tC: in situ temperature
! Tr: reference temperature
!
!INPUT PARAMETERS:
REAL, INTENT (IN) :: TC
! boltzman constant constant [ eV /K ]
REAL, PARAMETER   :: kb = 8.62d-5, Tr = 15.0

TK = -(1./kb)*(1./(273.15 + tC) - 1./(273.15 + Tr))
return 
END FUNCTION TK

PURE REAL FUNCTION alloscale(Topt_, mu0p, alpha)
IMPLICIT NONE
real, intent(in) :: Topt_     !Topt in ºC
real, intent(in) :: mu0p  !Normalized growth rate
real, intent(in) :: alpha    !Exponent of thermal traits normalized to z
alloscale =  mu0p * exp(TK(Topt_) * alpha) 
END FUNCTION alloscale

!phytoplankton growth rate, Chl:C and N:C function
SUBROUTINE GROWTH(t, N, par,  E_KN_, Topt_,  gr,  um)
IMPLICIT NONE
real, intent(in)  :: t     !Environmental Temperature (ºC)
real, intent(in)  :: N     !Nitrate (mmol m-3)
real, intent(in)  :: par   !light (mol E m-2 d-1)
real, intent(in)  :: E_KN_   !E of KN
real, intent(in)  :: Topt_  !Optimal temperature (ºC)
real, intent(out) :: gr    !Specific growth rate (d-1)
real, intent(out) :: um    !Phyto. maximal growth rate (excluding N and light limitation)
real, parameter   :: Ea0     = 0.98  
real, parameter   :: Eh0     = 4.43  
real, parameter   :: Ei      = 0.22  
real, parameter   :: beta    =-0.2  !Exponent for Ea0
real, parameter   :: phi     = 0.21  !Exponent for Eh
real, PARAMETER :: mumax0 = 0.59  !Normalized growth rate for mumax (d-1)

real :: par_, SI, fN, KN
real :: mumax, Ea, Eh

if (par <= 0.) stop "PAR is not positive!"

par_  = par/.4  !Change unit to W m-2
mumax = alloscale(Topt_, mumax0,  Ei) 
Ea    = alloscale(Topt_, Ea0,  beta) 
Eh    = alloscale(Topt_, Eh0,  phi) 
um    = JOHNSON(t, mumax, Ea, Eh, Topt_)

!Calculate KN
KN  = exp(TK(t)*E_KN_) * KN0
fN    = N / (N + KN)
SI    = 1.d0 - exp(-aI0*par_/um)
gr    = um*SI*fN  !Growth rate

RETURN
END SUBROUTINE GROWTH

SUBROUTINE MOMENTS
IMPLICIT NONE
INTEGER i

avgZ   =0.d0
sigmaZ=0.d0
do i = 1, NPHY
   avgZ = avgZ + TK(Topt(i)) * PHY(i)
enddo
avgZ = avgZ/PHYt

do i = 1, NPHY
   sigmaZ = sigmaZ + (TK(Topt(i)) - avgZ)**2  * PHY(i)
enddo
sigmaZ = sigmaZ/PHYt

END SUBROUTINE MOMENTS
END MODULE TOPT_MOD
