MODULE SNOWSUB

!-----------------------------------------------------------------------
! Module to calculate the sublimation of snowdrift particles
!-----------------------------------------------------------------------

CONTAINS
SUBROUTINE snowsublimation(dt, config_flags, tracer, rho_phy,                          &
                             dz8w, t_phy, p_phy, rmean,                                &
                             alt, snowsubl, tempch, moistch, tempchacc, moistchacc,    &
                             moist, swdown, psis,sigma, rhum,                          &
                             ustar, zmid, zlow, setvel, snowsublacc,                   &
                             ids,ide,jds,jde,kds,kde,                                  &
                             ims,ime,jms,jme,kms,kme,                                  &
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
 REAL,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t_phy, snowsublacc, tempchacc, moistchacc
 REAL,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(OUT) :: tempch, moistch, psis, sigma, rhum, snowsubl
 REAL,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: rho_phy, alt, p_phy,zmid, zlow, dz8w, setvel, rmean
 REAL,DIMENSION(ims:ime,jms:jme), INTENT(IN) :: ustar, swdown
 REAL,INTENT(IN) :: dt
 INTEGER         :: i, k, j

!Local variables
 REAL :: rhoi, nue, Kt, Ls, Rv, D, cp, pi
 REAL :: mmean, Re, Qr, N, Z, Nu, dm, qvs
 REAL :: esat, ew, rho_dry, snowsubl_limit, snowsubl_limit_tracer, snowsubl_tmp,        &
     moist_tmp, T_tmp, delta_qv_limit_0, snowsubl_limit_moist_tmp, delta_T_limit_tmp,                         & 
     T_limit_tmp, esat_tmp, qvs_tmp, delta_qv_limit
   
   


!---------------CONSTANTS------------------------------------------------
 rhoi = 900 !ice density [kg m-3]
 nue = 1.2E-5 !kinematic viscosity of air [m2 s-2]
 Kt = 0.024 !thermal conductivity of air [J s-1 m-1 K-1]
 Ls = 2.838E6 !latent heat of sublimation [J kg-1]
 Rv = 461.5 !gas constant for water vapor [J kg-1 K-1]
 cp = 1005 !specific heat capacity of dry air [J kg-1 K-1]
 D = 2.25E-5 !molecular diffusivity of water vapor in air [m2 s-1]
 pi = 3.14159265359
!-----------------------------------------------------------------------

DO j=jts,jte
 DO i=its,ite
  DO k=kts,kte !Height loop
  
  ! ----------------------------------------
  ! PREPARE AIR MOISTURE
  ! ----------------------------------------
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
   
   !limit water vapor mixing ratio to posivive values, relative humidity to 100%
   moist(i,k,j,p_qv) = MAX(0.0, moist(i,k,j,p_qv))
   rhum(i,k,j) = MIN(100., (moist(i,k,j,p_qv) / qvs) * 100)

   !calculate water vapor pressure
   ew = (rhum(i,k,j) * esat) / 100
   
   !water vapor deficit with respect to ice
   sigma(i,k,j) = (ew / esat) - 1
   
   IF (sigma(i,k,j) .gt. -0.001) THEN !if relative humidity > 99.9 % no sublimation
     snowsubl(i,k,j) = 0.0
     GO TO 105 !termination of snowdrift sublimation
   ELSE
   
   ! ---------------------------------------
   ! SNOW PARTICLES
   ! ---------------------------------------
   
     !calculation of mean particle mass, Schmid 2021, equation (29)
     mmean = (4/3) * pi * rhoi * rmean(i,k,j)**3
   
     !reynolds particle number, Schmid 2021, equation (26)
     Re = 2 * rmean(i,k,j) * setvel(i,k,j) / nue
     IF (Re .le. 10.) THEN
         Nu = 1.79 + 0.606 * Re**(0.5) !Schmid 2021, equation (24)
     ELSE
         Nu = 1.88 + 0.580 * Re**(0.5) !Schmid 2021, equation (25)
     ENDIF
     
     !calculate numerator Z and denominator N, Schmid 2021, equation (27)
     Z = 2 * pi * rmean(i,k,j) * sigma(i,k,j) !- (Qr/(Kt*Nu*t_phy(i,k,j)))*(Ls/(Rv*t_phy(i,k,j))-1)
     N = (Ls / (Kt*Nu*t_phy(i,k,j))) * (Ls/(Rv*t_phy(i,k,j))-1) + (Rv*t_phy(i,k,j)) / (Nu*D*esat)
     dm = (Z/N) * dt ! mass loss of single particle (kg)
     
     !sublimation loss rate, Schmid 2021, equation (28)
     psis(i,k,j) = dm / mmean
    
     !loss of sublimation per time step, Schmid 2021, equation (30)
     snowsubl_tmp = psis(i,k,j) * tracer(i,k,j,P_snowconc)
     
     ! limit snowsubl with tracer and qvs as upper and lower boundary to snowsubl_limit, use that to update snowsubl, moist, tracer,
     ! temperature.., limit snowsubl to <0
     snowsubl_tmp = MIN(0.0, snowsubl_tmp) ! limit snowsubl to negative values, msaigger
     
     snowsubl_limit_tracer  = MAX(-tracer(i,k,j,P_snowconc), snowsubl_tmp) ! limit snowsubl to available snow particles, msaigger
   
     ! ---------------------------------------
     ! TEMPERATURE AND MOISTURE FEEDBACK
     ! ---------------------------------------
     !following lines calculate the feedback of snowdrift sublimation
     IF (config_flags%subl_opt==2) then 
     
     ! conservative approach: calcualte first moisture change and temperature change out of sublimation
     ! -> restrict sublimation to saturation mixing ratio of new temperature
     ! -> use that restricted sublimation for moisture change and temperature change
         
         ! save current water vapor mixing ratio and temperature
         moist_tmp = moist(i,k,j,p_qv)
         T_tmp = t_phy(i,k,j)
         
         ! limit sublimation in first step to current saturation mixing ratio
         delta_qv_limit_0 = qvs - moist_tmp
         snowsubl_limit_moist_tmp = MAX(snowsubl_limit_tracer, -delta_qv_limit_0*rho_dry)
         
         ! calculate maximum temperature change
         delta_T_limit_tmp = (Ls*(-snowsubl_limit_moist_tmp)) / (cp*rho_dry) 
         
         ! calculate saturation mixing ratio for sublimation-cooled air
         T_limit_tmp = T_tmp - delta_T_limit_tmp
         esat_tmp = 611.2 * exp((Ls/Rv)*((1/273.15)-(1/T_limit_tmp)))
         qvs_tmp = (0.622 * esat_tmp) / (p_phy(i,k,j) - esat_tmp)
         
         ! limit sublimation to updated saturation mixing ratio
         delta_qv_limit = qvs_tmp - moist_tmp
         snowsubl_limit = MAX(snowsubl_limit_tracer, -delta_qv_limit*rho_dry)
         
         !! use limited sublimation to update mixing ratio and temperature
         ! update moisture
         moist(i,k,j,p_qv) = MIN(qvs_tmp, moist_tmp - snowsubl_limit/rho_dry)
   
         !write out rate of change of moisture due to snowdrift
         moistch(i,k,j) = (moist(i,k,j,p_qv) - moist_tmp) / dt ! moisture change rate (kg kg-1 s-1)
       
         ! accumulated moisture change
         moistchacc(i,k,j) = moistchacc(i,k,j) + moistch(i,k,j)*dt ! acc moisture change (kg kg-1)
         
         ! Update temperaure (Schmid 2021, equation (35) & (36))
         t_phy(i,k,j) = T_tmp - (Ls * (-snowsubl_limit) / (cp * rho_dry))

         ! write out rate of temperature change 
         tempch(i,k,j) = (t_phy(i,k,j) - T_tmp) / dt ! temperature change rate (K s-1)
      
         ! accumulate temperature change over time
         tempchacc(i,k,j) = tempchacc(i,k,j) + tempch(i,k,j)*dt ! acc temperature change (K)
          
     ELSE ! no sublimation feedback
         snowsubl_limit = snowsubl_limit_tracer
     ENDIF

     ! write out sublimation rate 
     snowsubl(i,k,j) = snowsubl_limit / dt  ! sublimation rate (kg m-3 s-1)
     
     ! accumulate snowsublimation
     snowsublacc(i,k,j) = snowsublacc(i,k,j) + snowsubl_limit ! accumulated sublimation (kg m-3)
   
     !update snow particle concentration
     tracer(i,k,j,P_snowconc) = MAX(0.0, tracer(i,k,j,P_snowconc) + snowsubl_limit)

   ENDIF !end sigma > 0
 ENDIF !end snowconc < d-8

 105 continue

   ENDDO !end k
 ENDDO !end i
ENDDO !end j

END SUBROUTINE snowsublimation
END MODULE SNOWSUB

   
   
   
