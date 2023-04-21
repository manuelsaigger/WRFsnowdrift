MODULE SNOWSET

!---------------------------------------------------------------------------
! Module to calculate the settling of snowdrift particles
!---------------------------------------------------------------------------

CONTAINS
SUBROUTINE snowsettling(dt, tracer, dz8w, fdep,              & 
           snoweracc, snowdepacc, snowhdrift, snowload,      &
           uth, ustar, snowh, snow, rhodrift,    &
           setvel, snowdepacckg,                             & !snicexy, &
           ids, ide, jds, jde, kds, kde,                     &
           ims, ime, jms, jme, kms, kme,                     &
           its, ite, jts, jte, kts, kte)
 USE module_configure         !/frame
 USE module_state_description !/frame 

 IMPLICIT NONE
 INTEGER,INTENT(IN) :: ids,ide, jds,jde, kds,kde,                                          &
                       ims,ime, jms,jme, kms,kme,                                          &
                       its,ite, jts,jte, kts,kte
 REAL,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer), INTENT(INOUT) :: tracer
! REAL,DIMENSION(ims:ime,-2:0,jms:jme),INTENT(INOUT) :: snicexy ! uncomment this line for real
 REAL,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) ::dz8w, setvel
 REAL,DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: snowdepacc, snow, snowh, snowload,      &
                                                 snowdepacckg, fdep
 REAL,DIMENSION(ims:ime,jms:jme), INTENT(OUT) :: snowhdrift
 REAL,DIMENSION(ims:ime,jms:jme), INTENT(IN) :: uth, ustar, snoweracc, rhodrift
 REAL,INTENT(IN) :: dt

!Local variables
 INTEGER         :: i, k, j
 REAL            :: fgrad(kte), downward, tracer_tmp, fdep_tmp ,                           &
         mass_in_kplus, mass_out_max, mass_out_k, mass_out_k_limit, mass_change_k, ustar_use

!---------------VARIABLES-------------------------------------------------
 !snow:  snow water equivalent [kg m-2]
 !snowh: snowheight [m]
 !phis:  particle concentration in saltation layer [kg m-3]
 !dz8w:  layer thickness [m]
 !dt:    time step [s]
 !uth, ustar: threshold friction velocity and friction velocity [m s-1]
 !fdep: snow deposition flux (lowest layer) [kg m-2 s-1]
 !snowdepacc: accumulated snow change due to deposition [m]
 !setvel = settling velocity [m s-1]
 !snowload: snow concentratin integrated over height [kg m-2] 
!------------------------------------------------------------------------

 snowload(its:ite,jts:jte)   = 0.0 !initial value

 DO j=jts,jte
   DO i=its,ite
     
     ! -------------------------------------
     ! SETTLING
     ! -------------------------------------
     !set upper boundary condition
     mass_in_kplus=0.0 ! mass into layer k out of layer k+1 (=0 at upper boundary), msaigger

     !loop starts from the highest model layer
     DO k = kte,kts,-1         

         ! save the current snow concentration
         tracer_tmp = tracer(i,k,j,P_snowconc)  !temporal storage of snow concentration

         ! calculate mass in/out of layer
         mass_out_max = tracer(i,k,j,P_snowconc)*dz8w(i,k,j) ! mass in layer k is maximum mass to go out
         mass_out_k = tracer(i,k,j,P_snowconc)*setvel(i,k,j)*dt ! flux out of layer k
         mass_out_k_limit = MIN(mass_out_max, mass_out_k) ! limited flux out of layer k to available mass
  
         IF (k==kts) THEN  ! deposition flux at lower boundary is calculated later
             mass_out_k_limit = 0.0 
         ENDIF
         mass_change_k = mass_in_kplus - mass_out_k_limit ! total mass change in layer k
         
         ! update concentration
         tracer(i,k,j,P_snowconc) = tracer(i,k,j,P_snowconc) + (mass_change_k/dz8w(i,k,j))
       
         ! prepare for next loop: mass out in k is mass in k+1
         mass_in_kplus = mass_out_k_limit
       
  
     ! -------------------------------------
     ! DEPOSITION
     ! -------------------------------------
         fdep_tmp = fdep(i,j) 
         IF (k .eq. kts) THEN
             !friction velocity must be smaller than its threshold
             IF (uth(i,j) .gt. ustar(i,j)) THEN
                 !deposition flux is now restricted to the maximal concentration in the layer
                 fdep(i,j) = MIN(tracer_tmp*dz8w(i,kts,j)/dt, tracer_tmp*setvel(i,kts,j)*((uth(i,j)**2-ustar(i,j)**2)/(uth(i,j)**2)))
             ELSE
                 fdep(i,j) = 0.0
             ENDIF
             fdep(i,j) = fdep(i,j) + fdep_tmp

             !update the snow concentration
             tracer(i,kts,j,P_snowconc) = MAX(0.0, tracer(i,kts,j,P_snowconc)- (fdep(i,j) * dt)/dz8w(i,kts,j))
         ENDIF !endif k=kts

     ! -------------------------------------
     ! OUTPUT FIELDS
     ! -------------------------------------
     
     ! calculate vertical integrated snow concentration (kg m-2)
     snowload(i,j) = snowload(i,j) + (tracer(i,k,j,P_snowconc)*dz8w(i,k,j)) 

     ENDDO !end k

     !update snow height
     snow(i,j) = snow(i,j) + fdep(i,j) * dt
     snowh(i,j) = snowh(i,j) + (fdep(i,j)*dt) / rhodrift(i,j) ! update snowheight as well, comment this line for real (???)
     ! snicexy(i,0,j) = snicexy(i,0,j) + fdep(i,j) * dt ! uncomment this line for real

     snowdepacc(i,j) = snowdepacc(i,j) + (fdep(i,j)*dt)/rhodrift(i,j) ! accumulated depositional snow height change (m)
     snowhdrift(i,j) = snowdepacc(i,j) - snoweracc(i,j) ! accumulated total snow height change (m)
     snowdepacckg(i,j) = snowdepacckg(i,j) + (fdep(i,j)*dt) ! accumulated depositional snow mass change (kg m-2)

   ENDDO !end i
 ENDDO !end j

END SUBROUTINE snowsettling
END MODULE SNOWSET

