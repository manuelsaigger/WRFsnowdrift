MODULE SNOWSUB

!-----------------------------------------------------------------------
! Module to calculate the sublimation of snowdrift particles
!-----------------------------------------------------------------------

CONTAINS
SUBROUTINE snowsublimation(dt, config_flags, tracer, rho_phy,  &
                             dz8w, t_phy, p_phy, rmean,         &
                             alt, snowsubl, tempch, moistch,    &
                             moist, swdown, psis,sigma, rhum,   &
                             ustar, zmid, zlow, setvel, snowsublacc, &
                             ids,ide,jds,jde,kds,kde,           &
                             ims,ime,jms,jme,kms,kme,           &
                             its,ite,jts,jte,kts,kte)
 USE module_configure         !/frame
 USE module_state_description !/frame

 IMPLICIT NONE
 TYPE(grid_config_rec_type),  INTENT(IN   ) :: config_flags
 INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde,             &
                        ims,ime, jms,jme, kms,kme,             &
                        its,ite, jts,jte, kts,kte

 REAL,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer), INTENT(INOUT) :: tracer
 REAL,DIMENSION(ims:ime,kms:kme,jms:jme,num_moist), INTENT(INOUT) ::  moist
 REAL,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t_phy, snowsublacc
 REAL,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(OUT) :: tempch, moistch, psis, sigma, rhum, snowsubl
 REAL,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: rho_phy, alt, p_phy,zmid, zlow, dz8w, setvel, rmean
 REAL,DIMENSION(ims:ime,jms:jme), INTENT(IN) :: ustar, swdown
 REAL,INTENT(IN) :: dt
 INTEGER         :: i, k, j

!Local variables
 REAL :: rhoi, nue, Kt, Ls, Rv, D, cp, pi
 REAL :: mmean, Re, Qr, N, Z, Nu, dm, qvs
 REAL :: esat, ew, wvsden, wvden, rho_dry, ttemp, moisttemp, snowsubl_limit, snowsubl_limit_tracer, snowsubl_limit_moist

!---------------CONSTANTS------------------------------------------------
 rhoi = 900 !ice density [kg/m³]
 nue = 1.2E-5 !kinematic viscosity of air [m²/s]
 Kt = 0.024 !thermal conductivity of air [W/mK]
 Ls = 2.838E6 !latent heat of sublimation [J/kg]
 Rv = 461.5 !gas constant for water vapor [J/kg K]
 cp = 1005 !specific heat capacity of dry air [J/kg K]
 D = 2.25E-5 !molecular diffusivity of water vapor in air [m²/s]
 pi = 3.14159265359
!-----------------------------------------------------------------------

DO j=jts,jte
 DO i=its,ite
  DO k=kts,kte !Height loop
  
  !if the snow concentration is very low, no sublimation of particles
  IF (tracer(i,k,j,P_snowconc) .le. 1E-8) THEN
   snowsubl(i,k,j)=0.0
   GO TO 105 !termination of snowdrift sublimation
  ELSE
   !calculate dry air mass density
   rho_dry = 1.0/alt(i,k,j)
   !saturation water vapor pressure (Clausius-Clapeyron)
   esat = 611.2 * exp((Ls/Rv)*((1/273.15)-(1/t_phy(i,k,j))))
   !saturation mixing ratio
   qvs = (0.622 * esat) / (p_phy(i,k,j) - esat)
   !cut off negative values of water vapor mixing ratio
   moist(i,k,j,p_qv) = MAX(0.0, moist(i,k,j,p_qv))
   rhum(i,k,j) = MIN(100., (moist(i,k,j,p_qv) / qvs) * 100)
   !water vapor pressure
   ew = (rhum(i,k,j) * esat)/100
   !water vapor deficit with respect to ice
   sigma(i,k,j) = (ew / esat) -1

  IF (sigma(i,k,j) .gt. -0.001) THEN !if relative humidity > 99.9 % no sublimation
   snowsubl(i,k,j)=0.0
   GO TO 105 !termination of snowdrift sublimation
  ELSE
   !calculation of mean particle mass, Schmid 2021, equation (29)
   mmean = (4/3) * pi * rhoi * rmean(i,k,j)**3
   !reynolds particle number, Schmid 2021, equation (26)
   Re = 2 * rmean(i,k,j) * setvel(i,k,j) / nue
   IF (Re .le. 10.) THEN
   Nu = 1.79 + 0.606 * Re**(0.5) !Schmid 2021, equation (24)
   ELSE
   Nu = 1.88 + 0.580 * Re**(0.5) !Schmid 2021, equation (25)
   ENDIF

   !Qr = pi * rmean**2 *(1-0.1) * swdown(i,j) !radiation term, neglected in this thesis
   !calculate numerator Z and denominator N, Schmid 2021, equation (27)
   Z = 2*3.1415 * rmean(i,k,j) * sigma(i,k,j) !- (Qr/(Kt*Nu*t_phy(i,k,j)))*(Ls/(Rv*t_phy(i,k,j))-1)
   N = (Ls / (Kt*Nu*t_phy(i,k,j))) * (Ls/(Rv*t_phy(i,k,j))-1) + (Rv*t_phy(i,k,j)) / (Nu*D*esat)
   dm = (Z/N) * dt !mass loss [kg]
   !sublimation loss rate, Schmid 2021, equation (28)
   psis(i,k,j) = dm / mmean
   !loss of sublimation per time step, Schmid 2021, equation (30)
   snowsubl(i,k,j) = psis(i,k,j) * tracer(i,k,j,P_snowconc)

   ! bugfix flux limiter
   ! limit snowsubl with tracer and qvs at upper and lower boundary to snowsubl_limit, use that to update snowsubl, moist, tracer,
   ! temperature.., limit snowsubl to <0, so that abs() in update moist and t_phy can be left out (ambigous..)
   snowsubl(i,k,j) = MIN(0.0, snowsubl(i,k,j)) ! limit snowsubl to negative values, msaigger
   snowsubl_limit_tracer  = MAX(-tracer(i,k,j,P_snowconc), snowsubl(i,k,j)) ! limit snowsubl to available snow particles, msaigger
   IF (config_flags%subl_opt==2) then
   snowsubl_limit_moist = MAX(snowsubl(i,k,j), -(moist(i,k,j,p_qv)+qvs)*rho_dry)
   snowsubl_limit = MAX(snowsubl_limit_tracer, snowsubl_limit_moist)
   ELSE
   snowsubl_limit = snowsubl_limit_tracer
   ENDIF

   ! accumulate snowsublimation
   snowsublacc(i,k,j) = snowsublacc(i,k,j) + snowsubl_limit
   
   !update snow particle concentration
   tracer(i,k,j,P_snowconc) = MAX(0.0, tracer(i,k,j,P_snowconc) + snowsubl_limit)

   !following lines calculate the feedback of snowdrift sublimation
   !if subl_opt = 1, no feedback is calculated
   IF (config_flags%subl_opt==2) then
   !save the current water vapor mixing ratio
   moisttemp = moist(i,k,j,p_qv)
   !water vapor mixing ratio increases due to snowdrift sublimation, Schmid 2021, equation (34)
   !can not go higher than the saturation value qvs
   moist(i,k,j,p_qv) = MIN(qvs, moist(i,k,j,p_qv) - snowsubl_limit/rho_dry)
   !change of moisture due to snowdrift
  ! moistch(i,k,j) = moist(i,k,j,p_qv) – moisttemp !ms 
   moistch(i,k,j)=moist(i,k,j,p_qv)-moisttemp 
  !suggestion: accumulate it over time, the above is the rate of change per time step
   !moistch(i,k,j) = moistch(i,k,j) + (moist(i,k,j,p_qv) – moisttemp)
   !temperature decreases due to snowdrift sublimation
   ttemp = t_phy(i,k,j)
   !Schmid 2021, equation (35) & (36)
   t_phy(i,k,j) = t_phy(i,k,j) - (Ls * (-snowsubl_limit) / (cp * rho_dry))
  ! tempch(i,k,j) = t_phy(i,k,j) – ttemp !ms
   tempch(i,k,j)=t_phy(i,k,j)-ttemp
   !suggestion: accumulate it over time, the above is the rate of change per time step
   !tempch(i,k,j) = tempch(i,k,j) + (t_phy(i,k,j) – ttemp)

   ENDIF !end subl_opt=2, feedback
  ENDIF !end sigma > 0
 ENDIF !end snowconc < d-8

 105 continue

   ENDDO !end k
 ENDDO !end i
ENDDO !end j

END SUBROUTINE snowsublimation
END MODULE SNOWSUB

