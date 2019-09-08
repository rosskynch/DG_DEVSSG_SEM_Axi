!    DG_DEVSSG_SEM
!    Copyright (C) 2008-2019 Ross M Kynch
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE fene_p_mp_module
  USE shared_data
  USE functions_module
  USE IO_module
  USE viscoelastic_module

  IMPLICIT NONE
  CONTAINS

  SUBROUTINE applyElasticStress_FENE_PMP
    IMPLICIT NONE
    stress_cont_to_stokes_xNm1=stress_cont_to_stokes_x
    stress_cont_to_stokes_yNm1=stress_cont_to_stokes_y

    CALL calcStress_weakform_FENE_PMP
    CALL integrate_divergence_of_local_stress

! EXJ scheme:
    f_x = f_x + time_beta_0*stress_cont_to_stokes_x + time_beta_1*stress_cont_to_stokes_xNm1
    f_y = f_y + time_beta_0*stress_cont_to_stokes_y + time_beta_1*stress_cont_to_stokes_yNm1

  END SUBROUTINE applyElasticStress_FENE_PMP

  DOUBLE PRECISION FUNCTION calculatePsi_FENE_PMP(lambdaD, localGradUxx_in, localGradUxy_in, &
    localGradUyx_in, localGradUyy_in, localGradUzz_in)
    IMPLICIT NONE
    DOUBLE PRECISION :: Dxy, I2, I3, extensionRate
    DOUBLE PRECISION, INTENT(IN) :: lambdaD, localGradUxx_in, localGradUxy_in, &
      localGradUyx_in, localGradUyy_in, localGradUzz_in
! Psi(e) = (cosh(lamda_d*e) - 1) / 2
!
! e = 3I_3(D) / I_2(D)
!
! where D = (1/2)*(gradU + gradU') (' denotes transpose).
! Note Dxx = gradUxx, Dyy = gradUyy, Dzz = gradUzz
! but Dxy = Dyx = (1/2)*(gradUxy+gradUyx)
!
! with principal invariants:
! I1(D) = tr(D) = Dxx + Dyy + Dzz (not required)
! I2(D) = (1/2)*(tr(D^2) - tr(D)^2) = Dxy*Dxy - Dxx*Dyy - Dxx*Dzz - Dyy*Dzz
! I3(D) = det(D) = Dxx*Dyy*Dzz - Dxy*Dyx*Dzz
!
! tr denotes trace, det denotes determinant, ' denotes transpose.

    Dxy = 0.5*(localGradUxy_in + localGradUyx_in)
    I2 = Dxy*Dxy - (localGradUxx_in*localGradUyy_in + localGradUxx_in*localGradUzz_in + localGradUyy_in*localGradUzz_in)
    I3 = localGradUxx_in*localGradUyy_in*localGradUzz_in - Dxy*Dxy*localGradUzz_in

    IF (abs(I2).lt.1d-15) THEN
      calculatePsi_FENE_PMP = 0d0;
    ELSE
      extensionRate = 3d0*I3 / I2
      calculatePsi_FENE_PMP = 0.5*(cosh(lambdaD*extensionRate) - 1d0)
    ENDIF
    
    RETURN
  END FUNCTION calculatePsi_FENE_PMP

  DOUBLE PRECISION FUNCTION calculateF_FENE_PMP(bValue, Cxx_in, Cyy_in, Czz_in)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: bValue, Cxx_in, Cyy_in, Czz_in
    DOUBLE PRECISION :: trace, denominator
! F(tr(C)) = 1 / (1 - (tr(C)/b^2))
!
! where b is the square of the maximum extension of the dumbbell (considered in the model).
    trace = Cxx_in + Cyy_in + Czz_in
    IF (abs(trace).lt.1d-15) THEN
      calculateF_FENE_PMP = 1d0
      RETURN
    ENDIF
    
    denominator = 1d0 - (trace / (bValue*bValue))
    IF (abs(denominator).lt.1d-15) THEN
      calculateF_FENE_PMP = 0d0
      RETURN
    ENDIF   

    calculateF_FENE_PMP = 1d0 / denominator
    RETURN

  END FUNCTION calculateF_FENE_PMP
  

  SUBROUTINE calcStress_weakform_FENE_PMP
    IMPLICIT NONE
    INTEGER :: i,ij,el,sum_counter,&
! LAPACK bits:
    info,ipiv(4)

    DOUBLE PRECISION :: Cxx_approx, Cxy_approx, Cyy_approx, Czz_approx,&
      tauWeConstant, psiValue, psiTimesDxy, f_of_trC, tempWeCoefficient,&
      sum_temp, sumxx, sumxy, sumyy, sumzz,&
! LAPACK bits:
      temp_matrix(4,4), temp_rhs(4)

    DOUBLE PRECISION :: tempVx(1:nptot),&
      convective_contrib_xx(0:NP1SQM1,numelm),&
      convective_contrib_xy(0:NP1SQM1,numelm),&
      convective_contrib_yy(0:NP1SQM1,numelm),&
      convective_contrib_zz(0:NP1SQM1,numelm),&
      temp_contrib_xx(0:NP1SQM1,numelm),&
      temp_contrib_xy(0:NP1SQM1,numelm),&
      temp_contrib_yy(0:NP1SQM1,numelm),&
      temp_contrib_zz(0:NP1SQM1,numelm),&
      temp_contrib_xxNm1(0:NP1SQM1,numelm),&
      temp_contrib_xyNm1(0:NP1SQM1,numelm),&
      temp_contrib_yyNm1(0:NP1SQM1,numelm),&
      temp_contrib_zzNm1(0:NP1SQM1,numelm),&
      tempCxx(0:NP1SQM1,numelm),&
      tempCxy(0:NP1SQM1,numelm),&
      tempCyy(0:NP1SQM1,numelm),&
      tempCzz(0:NP1SQM1,numelm),&
      tempCxxPrevious(0:NP1SQM1,numelm),&
      tempCxyPrevious(0:NP1SQM1,numelm),&
      tempCyyPrevious(0:NP1SQM1,numelm),&
      tempCzzPrevious(0:NP1SQM1,numelm)

    IF (coordflag.eq.0) THEN
! CARTESIAN CASE !
! Not implemented. Will need to base it on the OldroydB/Giesekus version within the main viscoelastic_module.
      print*, 'Error: The FENE-P-MP model is not yet implemented for Cartesian co-ordinates...'
      print*, 'Stopping'
      STOP
    ELSEIF (coordflag.eq.1) THEN
! CYLINDERICAL POLAR (AXISYMMETRIC) CASE !
!
! NOTE:
! The x co-ordinate in the cartesian case is equivalent to the z component
! and y is equivalent to r component.
! in the cylindrical polar coordinate case.
!
! Semi-Implicit iteration.
      IF (param_iterative_convection) THEN
        tempCxx = localCxx
        tempCxy = localCxy
        tempCyy = localCyy
        tempCzz = localCzz
        tempVx = V_x - mesh_velocity

        sum_counter = 0
        sum_temp = 1d0;        
        DO WHILE (sum_temp.gt.1d-9) 
          sum_counter = sum_counter + 1

          CALL calc_axisymm_convective_term(tempVx, V_y, tempCxx, tempCxy, tempCyy, tempCzz, &
            convective_contrib_xx, convective_contrib_xy, &
            convective_contrib_yy, convective_contrib_zz)
          DO el=1,numelm
! Now loop over all stress nodes and compute solutions or apply boundary conditions, etc.
            DO ij=0,NP1SQM1
              i = mapg(ij,el)
              IF (inflowflag(i)) THEN 
!  We will apply inflow boundary conditions on stress, so nothing to do.
                CYCLE
              ENDIF
! Calculate entries of the matrix for FENE-P-MP:
! DEVSS-G:
! Where G tensor is used for deformation terms as well as Strain Rate.
! We copy this into a matrix and solve with LAPACK.
              f_of_trC = calculateF_FENE_PMP(param_fene_b, tempCxx(ij,el), tempCyy(ij,el), tempCzz(ij,el))
              psiValue = calculatePsi_FENE_PMP(param_fene_lambdaD, localGradUxx(ij,el), localGradUxy(ij,el), &
                localGradUyx(ij,el), localGradUyy(ij,el), localGradUzz(ij,el))

! We make use of both gradU and D = (1/2)*(gradU + gradU') (' = transpose), but they're identical
! in all components except the Dxy and Dyx components (and Dxy==Dyx), so we only need to compute 
! Dxy for use in the off-diagonal terms. All other uses of D will be replaced by the gradU component
! which is already stored.

! First compute the diagonal terms of the linear system for the conformation tensor.
              temp_matrix = 0d0
              tempWeCoefficient = We*(psiValue - 1d0);

              temp_matrix(1,1) = f_of_trC + Wetime_constant1 + 2d0*tempWeCoefficient*localGradUxx(ij,el)

              temp_matrix(2,2) = f_of_trC + Wetime_constant1
              ! Should be:
              ! + tempWeCoefficient*(localGradUxx(ij,el) + localGradUyy(ij,el))
              ! but this is zero by mass conservation.

              temp_matrix(3,3) = f_of_trC + Wetime_constant1 + 2d0*tempWeCoefficient*localGradUyy(ij,el)

              temp_matrix(4,4) = f_of_trC + Wetime_constant1 + 2d0*tempWeCoefficient*localGradUzz(ij,el)

! Now the off-diagonal terms.
              psiTimesDxy = 0.5*psiValue*(localGradUxy(ij,el) + localGradUyx(ij,el))

              temp_matrix(1,2) = 2d0*We*(psiTimesDxy - localGradUyx(ij,el))

              temp_matrix(2,1) = We*(psiTimesDxy - localGradUxy(ij,el)) 

              temp_matrix(2,3) = We*(psiTimesDxy - localGradUyx(ij,el)) 

              temp_matrix(3,2) = 2d0*We*(psiTimesDxy - localGradUxy(ij,el))

! Calculate RHS entries
! BDFJ:
              temp_rhs(1) = 1d0 + Wetime_constant2*( time_alpha_0*localCxx(ij,el) + time_alpha_1*localCxxNm1(ij,el) ) &
                - We*convective_contrib_xx(ij,el)

              temp_rhs(2) = Wetime_constant2*( time_alpha_0*localCxy(ij,el) + time_alpha_1*localCxyNm1(ij,el) ) &
                - We*convective_contrib_xy(ij,el)

              temp_rhs(3) = 1d0 + Wetime_constant2*( time_alpha_0*localCyy(ij,el) + time_alpha_1*localCyyNm1(ij,el) ) &
                - We*convective_contrib_yy(ij,el)

               temp_rhs(4) = 1d0 + Wetime_constant2*( time_alpha_0*localCzz(ij,el) + time_alpha_1*localCzzNm1(ij,el) ) &
                - We*convective_contrib_zz(ij,el)

! Using LAPACK to solve the 4x4 matrix:
              call dgetrf( 4, 4, temp_matrix, 4, ipiv, info )
              IF (info.ne.0) THEN
                write(*,*) 'Error in calcStress_weakform_FENE_PMP:',el,info
                STOP
              ENDIF

              call dgetrs( 'N', 4, 1, temp_matrix, 4, ipiv, temp_rhs, 4, info )  
              IF (info.ne.0) THEN
                write(*,*) 'Error in calcStress_weakform_FENE_PMP:',el,info
                STOP
              ENDIF

! Copy solution from lapack:
              tempCxxPrevious = tempCxx
              tempCxyPrevious = tempCxy
              tempCyyPrevious = tempCyy
              tempCzzPrevious = tempCzz
              tempCxx(ij,el) = temp_rhs(1)
              tempCxy(ij,el) = temp_rhs(2)
              tempCyy(ij,el) = temp_rhs(3)
              tempCzz(ij,el) = temp_rhs(4)
            ENDDO
          ENDDO
          sumxx = 0d0
          sumxy = 0d0
          sumyy = 0d0
          sumzz = 0d0
          DO el=1,numelm
            DO ij=0,NP1SQM1
              sumxx = sumxx + abs(tempCxx(ij,el) - tempCxxPrevious(ij,el))
              sumxy = sumxy + abs(tempCxy(ij,el) - tempCxyPrevious(ij,el))
              sumyy = sumyy + abs(tempCyy(ij,el) - tempCyyPrevious(ij,el))
              sumzz = sumzz + abs(tempCzz(ij,el) - tempCzzPrevious(ij,el))
            ENDDO
          ENDDO
          sumxx = sumxx / (numelm*NP1SQ)
          sumxy = sumxy / (numelm*NP1SQ)
          sumyy = sumyy / (numelm*NP1SQ)
          sumzz = sumzz / (numelm*NP1SQ)

          sum_temp = (sumxx + sumxy + sumyy + sumzz)/4d0
          IF(sum_temp.gt.1d10.or.sum_counter.gt.1000) THEN
            print*,'Iterative scheme in constitutive equation failed to converge! Sum counter = ',sum_counter
            STOP
          ENDIF
        ENDDO
      ELSE
! EXJ VERSION
        tempVx = V_x-mesh_velocity
        CALL calc_axisymm_convective_term(tempVx, V_y, &
          localCxx, localCxy, localCyy, localCzz, &
          temp_contrib_xx, temp_contrib_xy, temp_contrib_yy, temp_contrib_zz)
        IF (param_time_order.eq.2) THEN
          tempVx=V_xNm1-mesh_velocityNm1
          CALL calc_axisymm_convective_term(tempVx, V_yNm1, &
            localCxxNm1, localCxyNm1, localCyyNm1, localCzzNm1, &
            temp_contrib_xxNm1, temp_contrib_xyNm1, temp_contrib_yyNm1, temp_contrib_zzNm1)
        ENDIF

        convective_contrib_xx = time_beta_0*temp_contrib_xx + time_beta_1*temp_contrib_xxNm1
        convective_contrib_xy = time_beta_0*temp_contrib_xy + time_beta_1*temp_contrib_xyNm1
        convective_contrib_yy = time_beta_0*temp_contrib_yy + time_beta_1*temp_contrib_yyNm1
        convective_contrib_zz = time_beta_0*temp_contrib_zz + time_beta_1*temp_contrib_zzNm1

        DO el=1,numelm

! Now loop over all stress nodes and compute solutions or apply boundary conditions, etc.
          DO ij=0,NP1SQM1
            i=mapg(ij,el)
            IF (inflowflag(i)) THEN 
!  We will apply inflow boundary conditions on stress, so nothing to do.
              CYCLE
            ENDIF
! Calculate entries of the matrix for FENE-P-MP:

! DEVSS-G:
! Where G tensor is used for deformation terms as well as Strain Rate.
! We copy this into a matrix and solve with LAPACK.
            temp_matrix=0d0

            Cxx_approx = time_beta_0*localCxx(ij,el) + time_beta_1*localCxxNm1(ij,el)
            Cyy_approx = time_beta_0*localCyy(ij,el) + time_beta_1*localCyyNm1(ij,el)
            Czz_approx = time_beta_0*localCzz(ij,el) + time_beta_1*localCzzNm1(ij,el)
            Cxy_approx = time_beta_0*localCxy(ij,el) + time_beta_1*localCxyNm1(ij,el)

            f_of_trC = calculateF_FENE_PMP(param_fene_b, Cxx_approx, Cyy_approx, Czz_approx)
            psiValue = calculatePsi_FENE_PMP(param_fene_lambdaD, localGradUxx(ij,el), localGradUxy(ij,el), &
              localGradUyx(ij,el), localGradUyy(ij,el), localGradUzz(ij,el))


! We make use of both gradU and D = (1/2)*(gradU + gradU') (' = transpose), but they're identical
! in all components except the Dxy and Dyx components (and Dxy==Dyx), so we only need to compute 
! Dxy for use in the off-diagonal terms. All other uses of D will be replaced by the gradU component
! which is already stored.

! First compute the diagonal terms of the linear system for the conformation tensor.
            tempWeCoefficient = We*(psiValue - 1d0);

            temp_matrix(1,1) = f_of_trC + Wetime_constant1 + 2d0*tempWeCoefficient*localGradUxx(ij,el)

            temp_matrix(2,2) = f_of_trC + Wetime_constant1 + tempWeCoefficient*(localGradUxx(ij,el) + localGradUyy(ij,el))

            temp_matrix(3,3) = f_of_trC + Wetime_constant1 + 2d0*tempWeCoefficient*localGradUyy(ij,el)

            temp_matrix(4,4) = f_of_trC + Wetime_constant1 + 2d0*tempWeCoefficient*localGradUzz(ij,el)

! Now the off-diagonal terms.
            psiTimesDxy = 0.5*psiValue*(localGradUxy(ij,el) + localGradUyx(ij,el))

            temp_matrix(1,2) = 2d0*We*(psiTimesDxy - localGradUyx(ij,el))

            temp_matrix(2,1) = We*(psiTimesDxy - localGradUxy(ij,el)) 

            temp_matrix(2,3) = We*(psiTimesDxy - localGradUyx(ij,el)) 

            temp_matrix(3,2) = 2d0*We*(psiTimesDxy - localGradUxy(ij,el))

! Calculate RHS entries
! BDFJ:
            temp_rhs(1) = 1d0 + Wetime_constant2*( time_alpha_0*localCxx(ij,el) + time_alpha_1*localCxxNm1(ij,el) ) &
              - We*convective_contrib_xx(ij,el)

            temp_rhs(2) = Wetime_constant2*( time_alpha_0*localCxy(ij,el) + time_alpha_1*localCxyNm1(ij,el) ) &
              - We*convective_contrib_xy(ij,el)

            temp_rhs(3) = 1d0 + Wetime_constant2*( time_alpha_0*localCyy(ij,el) + time_alpha_1*localCyyNm1(ij,el) ) &
              - We*convective_contrib_yy(ij,el)

            temp_rhs(4) = 1d0 + Wetime_constant2*( time_alpha_0*localCzz(ij,el) + time_alpha_1*localCzzNm1(ij,el) ) &
              - We*convective_contrib_zz(ij,el)

! Using LAPACK to solve the 4x4 matrix:
            call dgetrf( 4, 4, temp_matrix, 4, ipiv, info )
            IF (info.ne.0) THEN
              write(*,*) 'Error in calcStress_weakform_FENE_PMP:',el,info
              STOP
            ENDIF

            call dgetrs( 'N', 4, 1, temp_matrix, 4, ipiv, temp_rhs, 4, info )  
            IF (info.ne.0) THEN
              write(*,*) 'Error in calcStress_weakform_FENE_PMP:',el,info
              STOP
            ENDIF

! Copy solution from lapack:
            tempCxx(ij,el) = temp_rhs(1)
            tempCxy(ij,el) = temp_rhs(2)
            tempCyy(ij,el) = temp_rhs(3)
            tempCzz(ij,el) = temp_rhs(4)
          ENDDO
        ENDDO
      ENDIF

! Update stored valeus of Cpq and Tpq 
! ie move time along by one for these, so that localTpq and localCpq now hold the current stress/conformation values
      localCxxNm2 = localCxxNm1
      localCxyNm2 = localCxyNm1
      localCyyNm2 = localCyyNm1
      localCzzNm2 = localCzzNm1

      localCxxNm1 = localCxx
      localCxyNm1 = localCxy
      localCyyNm1 = localCyy
      localCzzNm1 = localCzz

      localCxx = tempCxx
      localCxy = tempCxy
      localCyy = tempCyy
      localCzz = tempCzz
      
      localTxxNm2 = localTxxNm1
      localTxyNm2 = localTxyNm1
      localTyyNm2 = localTyyNm1
      localTzzNm2 = localTzzNm1

      localTxxNm1 = localTxx
      localTxyNm1 = localTxy
      localTyyNm1 = localTyy
      localTzzNm1 = localTzz
! Calculate the new value of tau from the conformation tensor that's just been computed (currently stored as temp).
      tauWeConstant = (1d0 - param_beta)/We
      DO el=1,numelm
        DO ij=0,NP1SQM1
          IF (inflowflag(mapg(ij,el))) THEN 
            localTxx(ij,el) = boundary_stress_xx(i)
            localTxy(ij,el) = boundary_stress_xy(i)
            localTyy(ij,el) = boundary_stress_yy(i)
            localTzz(ij,el) = boundary_stress_zz(i)
          ELSE
            f_of_trC = calculateF_FENE_PMP(param_fene_b, localCxx(ij,el), localCyy(ij,el), localCzz(ij,el))
            localTxx(ij,el) = tauWeConstant*(f_of_trC*localCxx(ij,el) - 1d0)
            localTxy(ij,el) = tauWeConstant*f_of_trC*localCxy(ij,el)
            localTyy(ij,el) = tauWeConstant*(f_of_trC*localCyy(ij,el) - 1d0)
            localTzz(ij,el) = tauWeConstant*(f_of_trC*localCzz(ij,el) - 1d0)
          ENDIF
        ENDDO
      ENDDO
    ELSE
      print*, 'Error: No co-ordinate system specified...'
      print*, 'Stopping'
      STOP
    ENDIF
  END SUBROUTINE calcStress_weakform_FENE_PMP

END MODULE fene_p_mp_module