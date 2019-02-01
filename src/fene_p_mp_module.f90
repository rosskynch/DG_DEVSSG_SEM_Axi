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
! I2(D) = (1/2)*(tr(D)^2 - tr(D^2)) = Dxx*Dyy + Dxx*Dzz + Dyy*Dzz - Dxy*Dxy 
! I3(D) = det(D) = Dxx*Dyy*Dzz - Dxy*Dyx*Dzz
!
! tr denotes trace, det denotes determinant, ' denotes transpose.

    Dxy = 0.5*(localGradUxy_in + localGradUyx_in)
    I2 = localGradUxx_in*localGradUyy_in + localGradUxx_in*localGradUzz_in + localGradUyy_in*localGradUzz_in - Dxy*Dxy
    I3 = localGradUxx_in*localGradUyy_in*localGradUzz_in - Dxy*Dxy*localGradUzz_in

    extensionRate = 3d0*I3 / I2

    calculatePsi_FENE_PMP = 0.5*(cosh(lambdaD*extensionRate) - 1d0)  
    RETURN

  END FUNCTION calculatePsi_FENE_PMP

  DOUBLE PRECISION FUNCTION calculateF_FENE_PMP(bValue, Cxx_in, Cyy_in, Czz_in)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: bValue, Cxx_in, Cyy_in, Czz_in
! F(tr(C)) = 1 / (1 - (tr(C)/b^2))
!
! where b is the square of the maximum extension of the dumbbell (considered in the model).

    calculateF_FENE_PMP = 1d0 / (1d0 - (Cxx_in + Cyy_in + Czz_in) / (bValue*bValue))
    RETURN

  END FUNCTION calculateF_FENE_PMP
  

  SUBROUTINE calcStress_weakform_FENE_PMP
    IMPLICIT NONE
    INTEGER :: i,ij,el,&
! LAPACK bits:
    info,ipiv(4)

    DOUBLE PRECISION :: Cxx_approx, Cxy_approx, Cyy_approx, Czz_approx,&
      TauWeConstant, PsiValue, WePsiValue, Dxy, f_of_trC, bValue, lambdaD,&
! LAPACK bits:
      temp_matrix(4,4),temp_rhs(4)

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
      tempCxxNext(0:NP1SQM1,numelm),&
      tempCxyNext(0:NP1SQM1,numelm),&
      tempCyyNext(0:NP1SQM1,numelm),&
      tempCzzNext(0:NP1SQM1,numelm),&
      tempTxx(0:NP1SQM1,numelm),&
      tempTxy(0:NP1SQM1,numelm),&
      tempTyy(0:NP1SQM1,numelm),&
      tempTzz(0:NP1SQM1,numelm)!,&
      !tempTxxNext(0:NP1SQM1,numelm),&
      !tempTxyNext(0:NP1SQM1,numelm),&
      !tempTyyNext(0:NP1SQM1,numelm),&
      !tempTzzNext(0:NP1SQM1,numelm)

    temp_contrib_xx=0d0
    temp_contrib_xy=0d0
    temp_contrib_yy=0d0
    temp_contrib_zz=0d0
    temp_contrib_xxNm1=0d0
    temp_contrib_xyNm1=0d0
    temp_contrib_yyNm1=0d0
    temp_contrib_zzNm1=0d0


! If using the semi-implicit iterative method, then we will iterate with tempTpq as the solution at time N+1,
! and localTpq as Tpq at time N, and localTpqNm1 as Txx at time N-1. We then update Txx with the converged value.
! 
! For now this is disabled.
!    IF (param_iterative_convection) THEN
!      tempTxx = localTxx
!      tempTxy = localTxy
!      tempTyy = localTyy
!      tempTzz = localTzz
!
!      tempTxxNext = localTxx
!      tempTxyNext = localTxy
!      tempTyyNext = localTyy
!      tempTzzNext = localTzz
!
!      tempCxx = localCxx
!      tempCxy = localCxy
!      tempCyy = localCyy
!      tempCzz = localCzz
!
!      tempCxxNext = localCxx
!      tempCxyNext = localCxy
!      tempCyyNext = localCyy
!      tempCzzNext = localCzz
!    ENDIF


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
! Not implemented. Will need to base it on the OldroydB/Giesekus version within the main viscoelastic_module.
        print*, 'Error: The FENE-P-MP model is not yet implemented for semi-implicit iterative scheme...'
        print*, 'Stopping'
        STOP
      ELSE
! EXJ VERSION
        tempVx=V_x-mesh_velocity
        CALL calc_axisymm_convective_term(tempVx,V_y,localCxx,localCxy,localCyy, localCzz, &
          temp_contrib_xx,temp_contrib_xy, &
          temp_contrib_yy,temp_contrib_zz)
        IF (param_time_order.eq.2) THEN
          tempVx=V_xNm1-mesh_velocityNm1
          CALL calc_axisymm_convective_term(tempVx,V_yNm1,localCxxNm1,localCxyNm1,localCyyNm1,localCzzNm1, &
            temp_contrib_xxNm1,temp_contrib_xyNm1, &
            temp_contrib_yyNm1,temp_contrib_zzNm1)
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
!  Apply inflow boundary conditions on stress.
              tempCxx(ij,el) = boundary_stress_xx(i)
              tempCxy(ij,el) = boundary_stress_xy(i)
              tempCyy(ij,el) = boundary_stress_yy(i)
              tempCzz(ij,el) = boundary_stress_yy(i)
              CYCLE
            ENDIF
! Calculate entries of the matrix for OldroydB:

! DEVSS-G:
! Where G tensor is used for deformation terms as well as Strain Rate.
! We copy this into a matrix and solve with LAPACK.
            temp_matrix=0d0

            Cxx_approx = time_beta_0*localCxx(ij,el) + time_beta_1*localCxxNm1(ij,el)
            Cyy_approx = time_beta_0*localCyy(ij,el) + time_beta_1*localCyyNm1(ij,el)
            Czz_approx = time_beta_0*localCzz(ij,el) + time_beta_1*localCzzNm1(ij,el)
            Cxy_approx = time_beta_0*localCxy(ij,el) + time_beta_1*localCxyNm1(ij,el)

! TODO: Ask Tim where these values should be defined?
            bValue = 1d0; 
            lambdaD = 1d0;

            f_of_trC = calculateF_FENE_PMP(bValue, Cxx_approx, Cyy_approx, Czz_approx)
            psiValue = calculatePsi_FENE_PMP(lambdaD, localGradUxx(ij,el), localGradUxy(ij,el), &
              localGradUyx(ij,el), localGradUyy(ij,el), localGradUzz(ij,el))

            WePsiValue = We*psiValue

! TODO The extra terms for FENE P-MP can be simplified and combined in the terms used for OldB/Giesekus.
! Have not done for now so readability isn't harmed.
            temp_matrix(1,1) = f_of_trC + Wetime_constant1 + 2d0*We*localGradUxx(ij,el) &
              + WePsiValue*2d0*localGradUxx(ij,el) !Dxx

            temp_matrix(2,2) = f_of_trC + Wetime_constant1 -  We*(localGradUxx(ij,el) + localGradUyy(ij,el)) &
              + WePsiValue*(localGradUxx(ij,el) + localGradUyy(ij,el)) !Dxx + Dyy

            temp_matrix(3,3) = f_of_trC + Wetime_constant1 - 2d0*We*localGradUyy(ij,el) &
              + WePsiValue*2d0*localGradUyy(ij,el) !Dyy

            temp_matrix(4,4) = f_of_trC + Wetime_constant1 - 2d0*We*localGradUzz(ij,el) &
              + WePsiValue*2d0*localGradUzz(ij,el) !Dzz

! We make use of both gradU and D = (1/2)*(gradU + gradU') (' = transpose), but they're identical
! in all components except the Dxy and Dyx components (and Dxy==Dyx), so we only need to compute 
! Dxy for use below, all other uses of D will be replaced by the gradU component which is already stored.

            Dxy = 0.5*(localGradUxy(ij,el) + localGradUyx(ij,el))

            temp_matrix(1,2) = -2d0*We*localGradUyx(ij,el) &
              + 2d0*WePsiValue*Dxy
              
            temp_matrix(2,1) = -We*localGradUxy(ij,el) &
             + WePsiValue*Dxy

            temp_matrix(2,3) = -We*localGradUyx(ij,el) &
              + WePsiValue*Dxy

            temp_matrix(3,2) = -2d0*We*localGradUxy(ij,el) &
              + 2d0*WePsiValue*Dxy

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
              write(*,*) 'Error in calcStress_weakform:',el,info
              STOP
            ENDIF

            call dgetrs( 'N', 4, 1, temp_matrix, 4, ipiv, temp_rhs, 4, info )  
            IF (info.ne.0) THEN
              write(*,*) 'Error in calcStress_weakform:',el,info
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
      
! Calculate the new value of tau from the conformation tensor that's just been computed (currently stored as temp).
      TauWeConstant = (1d0 - param_beta)/We
      DO el=1,numelm
        DO ij=0,NP1SQM1
          f_of_trC = calculateF_FENE_PMP(bValue, tempCxx(ij,el), tempCzz(ij,el), tempCzz(ij,el))
          tempTxx(ij,el) = TauWeConstant*(f_of_trC*tempCxx(ij,el) - 1d0)
          tempTxy(ij,el) = TauWeConstant*f_of_trC*tempCxy(ij,el)
          tempTyy(ij,el) = TauWeConstant*(f_of_trC*tempCyy(ij,el) - 1d0)
          tempTzz(ij,el) = TauWeConstant*(f_of_trC*tempCzz(ij,el) - 1d0)
        ENDDO        
      ENDDO

! Update stored valeus of Cpq and Tpq 
! ie move time along by one for these, so that localTpq and localCpq now hold the current stress/conformation values
      localCxxNm2=localCxxNm1
      localCxyNm2=localCxyNm1
      localCyyNm2=localCyyNm1
      localCzzNm2=localCzzNm1

      localCxxNm1=localCxx
      localCxyNm1=localCxy
      localCyyNm1=localCyy
      localCzzNm1=localCzz

      localCxx=tempCxx
      localCxy=tempCxy
      localCyy=tempCyy
      localCzz=tempCzz
      
      localTxxNm2=localTxxNm1
      localTxyNm2=localTxyNm1
      localTyyNm2=localTyyNm1
      localTzzNm2=localTzzNm1

      localTxxNm1=localTxx
      localTxyNm1=localTxy
      localTyyNm1=localTyy
      localTzzNm1=localTzz

      localTxx=tempTxx
      localTxy=tempTxy
      localTyy=tempTyy
      localTzz=tempTzz
    ELSE
      print*, 'Error: No co-ordinate system specified...'
      print*, 'Stopping'
      STOP
    ENDIF
  END SUBROUTINE calcStress_weakform_FENE_PMP

END MODULE fene_p_mp_module