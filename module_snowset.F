MODULE SNOWSET

!---------------------------------------------------------------------------
! Module to calculate the settling of snowdrift particles
!---------------------------------------------------------------------------

CONTAINS
SUBROUTINE snowsettling(dt, tracer, dz8w, fdep,              & 
           snoweracc, snowdepacc, snowhdrift, snowload,      &
           uth, ustar, ustarcor, snowh, snow,rhodrift,       &
           setvel, snowdepacckg,                             & 
           ids, ide, jds, jde, kds, kde,                     &
           ims, ime, jms, jme, kms, kme,                     &
           its, ite, jts, jte, kts, kte)
 USE module_configure         !/frame
 USE module_state_description !/frame 

 IMPLICIT NONE
 INTEGER,INTENT(IN) :: ids,ide, jds,jde, kds,kde,                  &
                       ims,ime, jms,jme, kms,kme,                  &
                       its,ite, jts,jte, kts,kte
 REAL,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer), INTENT(INOUT) :: tracer
! REAL,DIMENSION(ims:ime,-2:0,jms:jme),INTENT(INOUT) :: snicexy
 REAL,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) ::dz8w, setvel
 REAL,DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: snowdepacc, snow, snowh, snowload, &
                                                 snowdepacckg
 REAL,DIMENSION(ims:ime,jms:jme), INTENT(OUT) :: fdep, snowhdrift
 REAL,DIMENSION(ims:ime,jms:jme), INTENT(IN) :: uth, ustar, ustarcor, snoweracc, rhodrift
 REAL,INTENT(IN) :: dt

!Local variables
 INTEGER         :: i, k, j
 REAL            :: fgrad(kte), downward, tracer_tmp, snow_old, &
         mass_in_kplus, mass_out_max, mass_out_k, mass_out_k_limit, mass_change_k

!---------------VARIABLES-------------------------------------------------
 !snow:  snow water equivalent [kg/m²]
 !snowh: snowheight [m]
 !phis:  particle concentration in saltation layer [kg/m³]
 !dz8w:  layer thickness [m]
 !dt:    time step [s]
 !uth, ustar: threshold friction velocity and friction velocity [m/s]
 !fdep: snow deposition flux (lowest layer) [kg/m²s]
 !snowdepacc: accumulated snow change due to deposition [m]
 !setvel = settling velocity [m/s]
 !snowload: snow concentratin integrated over height [kg/m2] 
!------------------------------------------------------------------------

 snowload(its:ite,jts:jte)   = 0.0 !initial value

 DO j=jts,jte
   DO i=its,ite

!set lower boundary condition
    ! downward=0.0
     mass_in_kplus=0.0 ! mass into layer k out of layer k+1 (=0 at upper boundary), msaigger

!loop starts from the highest model layer
     DO k = kte,kts,-1         

!save the current snow concentration
       tracer_tmp = tracer(i,k,j,P_snowconc)  !temporal storage of snow concentration

!calculate the concentration that falls into the layer below
 ! fgrad(k) = MIN(tracer(i,k,j,P_snowconc), tracer(i,k,j,P_snowconc) * dt * setvel(i,k,j) * (1/dz8w(i,k,j)))   ! in orig.

  ! new settling msaigger: calculate mass in/out not concentration in/out
       mass_out_max = tracer(i,k,j,P_snowconc)*dz8w(i,k,j) ! mass in layer k
       mass_out_k = tracer(i,k,j,P_snowconc)*setvel(i,k,j)*dt ! flux out of layer k
       mass_out_k_limit = MIN(mass_out_max, mass_out_k) ! limited flux out of layer k
  
       IF (k==kts) THEN  ! deposition flux at lower boundary is calculated later
         mass_out_k_limit = 0.0 
       ENDIF
       mass_change_k = mass_in_kplus - mass_out_k_limit
       tracer(i,k,j,P_snowconc) = tracer(i,k,j,P_snowconc) + (mass_change_k/dz8w(i,k,j))
       mass_in_kplus = mass_out_k_limit
       
  
  !for the lowest layer this will 0, for now, because its calculated later
 ! fgrad(kts) = 0.0      ! in orig.

!update the snow concentration, downward is the concentration from above
!for the highest layer: downward is 0
 ! tracer(i,k,j,P_snowconc) = MAX(tracer(i,k,j,P_snowconc) - fgrad(k)  + downward, 0.0)  ! in orig.
!for the layer below
 ! downward = fgrad(k)    ! in orig.

!old code!, used for the thesis, don’t use it: can produce artificial mass
! IF (k .eq. kts) THEN
! fdep(i,j) = tracer_tmp * setvel(i,kts,j) * MAX((uth(i,j)**2-ustarcor(i,j)**2)/(uth(i,j)**2),0.0)
! if ((fdep(i,j) * dt)/dz8w(i,kts,j) .gt. tracer(i,kts,j,P_snowconc)) then
! endif
! tracer(i,kts,j,P_snowconc) = MAX(0.0, tracer(i,kts,j,P_snowconc) - (fdep(i,j) * dt)/dz8w(i,kts,j))
! ENDIF !endif k=kts
!use this update instead!
!deposition flux from the lowest model layer to the ground

       IF (k .eq. kts) THEN
   !friction velocity must be smaller than its threshold
         IF (uth(i,j) .gt. ustarcor(i,j)) THEN
   !deposition flux is now restricted to the maximal concentration in the layer
           fdep(i,j) = MIN(tracer_tmp*dz8w(i,kts,j)/dt, tracer_tmp*setvel(i,kts,j)*((uth(i,j)**2-ustarcor(i,j)**2)/(uth(i,j)**2)))
         ELSE
           fdep(i,j) = 0.0
         ENDIF

  !update the snow concentration
         tracer(i,kts,j,P_snowconc) = MAX(0.0, tracer(i,kts,j,P_snowconc)- (fdep(i,j) * dt)/dz8w(i,kts,j))
       ENDIF !endif k=kts

  !calculate vertical integrated snow concentration, updated with layer thickness, msaigger
       snowload(i,j) = snowload(i,j) + (tracer(i,k,j,P_snowconc)*dz8w(i,k,j)) 

     ENDDO !end k

 !update snow height
   snow_old = snow(i,j)
   snow(i,j) = snow(i,j) + fdep(i,j) * dt
 ! snicexy(i,0,j) = snicexy(i,0,j) + fdep(i,j) * dt

  !accumulated snowdrift deposition, delete /rhodrift to get it in mm snow water equivalent
   snowdepacc(i,j) = snowdepacc(i,j) + (snow(i,j) - snow_old)/rhodrift(i,j) ![m]
   snowhdrift(i,j) = snowdepacc(i,j) - snoweracc(i,j)
   snowdepacckg(i,j) = snowdepacckg(i,j) + (snow(i,j) - snow_old) ! (kg m-2)

   ENDDO !end i
 ENDDO !end j

END SUBROUTINE snowsettling
END MODULE SNOWSET

