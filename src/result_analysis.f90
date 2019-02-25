!    DG_DEVSSG_SEM
!    Copyright (C) 2008-2017 Ross M Kynch
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

MODULE result_analysis
! OLDER VERSION IS IN SNAPSHOT DIRECTORY - PRE 03/11/12
! THIS USED GLOBAL GRADIENTS AND OLDER ROUTINES TO CALC ANALYTICAL SOLNS.
  USE constants
  USE shared_data
  USE functions_module
  IMPLICIT NONE

  CONTAINS
! OLD VERSION (Re/We):
  SUBROUTINE calc_drag
    IMPLICIT NONE
    
! REMOVED ADAPTIVE DEVSS(G) PARTS FROM THIS - they can be found in result_analysis.f90 in snapshots pre 03/12/2012
  
    INTEGER :: i,j,k,l,el,m,ij,jj
    DOUBLE PRECISION :: thetam,temp_beta,&
          p_m,dudx,dudy,dvdx,temp1,temp2,temp3,drag_m
    
    dragNm2=dragNm1
    dragNm1=drag
    drag=0d0
    
    IF (coordflag.eq.0) THEN
! Cartesian co-ordinates, so we are computing drag force on a cylinder.
      DO el=1,numelm
        IF (circflag(el)) THEN
          drag_m=0d0

          temp_beta = param_beta_s - param_beta
          DO m=0,N
! Use blending functions to compute theta_m (G-L point between Theta1 and Theta2)
! This is our mapping from theta into a 1-d G-L point so we can use our usual numerical quadrature.
! Will need to multiply by the differential when completing the integration.

            thetam = PI - ((1d0-gl(m))*theta(1,el)+(1d0+gl(m))*theta(2,el))*0.5d0

            IF (param_beta_s.gt.0d0) THEN
              drag_m = drag_m + ( ( localpressure(m,el) &
                - 2d0*param_beta_s*localGradUxx(m,el) - localTxx(m,el) + 2d0*temp_beta*devss_soln_xx(m,el) &
                )*dcos(thetam) + &
                ( param_beta_s*( localGradUxy(m,el) + localGradUyx(m,el) ) + localTxy(m,el) &
                - temp_beta*(devss_soln_xy(m,el) + devss_soln_yx(m,el)) &
                )*dsin(thetam) )*w(m)
            ELSE
              drag_m = drag_m + ( ( localpressure(m,el) - 2d0*(param_beta)*localGradUxx(m,el) - localTxx(m,el) )*dcos(thetam) + &
                ( (param_beta)*( localGradUxy(m,el) + localGradUyx(m,el) ) + localTxy(m,el) )*dsin(thetam) )*w(m)
            ENDIF
          ENDDO
          drag = drag - (theta(2,el)-theta(1,el))*drag_m*0.5d0 ! Multiply by the differential of the mapping from theta space into G-L point space.
        ENDIF
      ENDDO
      drag=drag*2d0
      drag_star = -1d0
    ELSEIF (coordflag.eq.1) THEN
! Axisymmetric Cylindrical co-ordinates so we are computing drag force on a sphere.
      DO el=1,numelm
        IF (circflag(el)) THEN
          drag_m=0d0
          temp_beta=param_beta_s - param_beta
          DO m=0,N
! Use blending functions to compute theta_m (G-L point between Theta1 and Theta2)
! This is our mapping from theta into a 1-d G-L point so we can use our usual numerical quadrature.
! Will need to multiply by the differential when completing the integration.

            thetam = PI - ((1d0-gl(m))*theta(1,el)+(1d0+gl(m))*theta(2,el))*0.5d0
            
! Now calculate the contribution to the drag for this node
            IF (param_beta_s.gt.0d0) THEN
!             
              drag_m = drag_m + ( ( &
                localpressure(m,el) - 2d0*(param_beta_s)*localGradUxx(m,el) - localTxx(m,el) + 2d0*temp_beta*devss_soln_xx(m,el) &
                )*dcos(thetam) &
                + ( (param_beta_s)*( localGradUxy(m,el) + localGradUyx(m,el) ) + localTxy(m,el) &
                - temp_beta*( devss_soln_xy(m,el) + devss_soln_yx(m,el) ) &
                )*dsin(thetam) )*dsin(thetam)*w(m)
            ELSE
              drag_m = drag_m + ((localpressure(m,el) - 2d0*param_beta*localGradUxx(m,el) - localTxx(m,el) )*dcos(thetam) + &
                (param_beta*(localGradUxy(m,el) + localGradUyx(m,el)) + localTxy(m,el) )*dsin(thetam))*dsin(thetam)*w(m)
            ENDIF
          ENDDO
          
          drag = drag - (theta(2,el)-theta(1,el))*drag_m*0.5d0 ! Multiply by the differential of the mapping from theta space into G-L point space.
        ENDIF
      ENDDO
      drag=drag*(2d0)*PI*rad_sphere**2
      drag_star = drag/(6d0*pi*V_sphere*rad_sphere)
    ELSE
      print*, 'Error: No co-ordinate system specified...'
      print*, 'Stopping'
      STOP
    ENDIF
  END SUBROUTINE calc_drag
  
  SUBROUTINE extrapolate_pressure_local(out_pressure)
    
    IMPLICIT NONE
    INTEGER :: i,j,jj,ij,intij,k,l,ll,intkl,el
    DOUBLE PRECISION, INTENT(OUT) :: out_pressure(0:NP1SQM1,numelm)
    DOUBLE PRECISION :: temp

    DO el=1,numelm
! Extrapolate edge of element values.
! Edge 1
      j=0
      DO i=0,N
        ij=i
        temp = 0d0
        DO l=1,NM1
          ll=(l-1)*NM1
          DO k=1,NM1
            intkl= k+ll
            temp = temp + pressure(mapg_pressure(intkl,el))*evalh(k,i)*evalh(l,j)
          ENDDO
        ENDDO
        out_pressure(ij,el)=temp
      ENDDO

      i=N
      DO j=1,N
        ij=i+j*NP1
        temp = 0d0
        DO l=1,NM1
          ll=(l-1)*NM1
          DO k=1,NM1
            intkl= k+ll
            temp = temp + pressure(mapg_pressure(intkl,el))*evalh(k,i)*evalh(l,j)
          ENDDO
        ENDDO
        out_pressure(ij,el)=temp
      ENDDO
    
      j=N
      jj=j*NP1
      DO i=0,N-1
        ij=i+jj
        temp = 0d0
        DO l=1,NM1
          ll=(l-1)*NM1
          DO k=1,NM1
            intkl= k+ll
            temp = temp + pressure(mapg_pressure(intkl,el))*evalh(k,i)*evalh(l,j)
          ENDDO
        ENDDO
        out_pressure(ij,el)=temp
      ENDDO
    
      i=0
      DO j=1,NM1
        ij=j*NP1
        temp = 0d0
        DO l=1,NM1
          ll=(l-1)*NM1
          DO k=1,NM1
            intkl= k+ll
            temp = temp + pressure(mapg_pressure(intkl,el))*evalh(k,i)*evalh(l,j)
          ENDDO
        ENDDO
        out_pressure(ij,el)=temp
      ENDDO
! Fill in internal nodes:
      DO intij=1,NM1SQ
        ij=interior_to_local_node(intij)
        out_pressure(ij,el) = pressure(mapg_pressure(intij,el))
      ENDDO   
    ENDDO
  END SUBROUTINE extrapolate_pressure_local
  
  SUBROUTINE extrapolate_velocity_local(out_x,out_y)
    IMPLICIT NONE
    INTEGER :: ij,globi,el
    DOUBLE PRECISION, INTENT(OUT) :: out_x(0:NP1SQM1,numelm),out_y(0:NP1SQM1,numelm)
    
    DO el=1,numelm
      DO ij=0,NP1SQM1
        out_x(ij,el)=V_x(mapg(ij,el))
        out_y(ij,el)=V_y(mapg(ij,el))
      ENDDO
    ENDDO
  END SUBROUTINE extrapolate_velocity_local
  
  
  SUBROUTINE calc_analytical_solutions
    IMPLICIT NONE
    INTEGER :: ij,globi,el
    
    DO el=1,numelm
      DO ij=0,NP1SQM1
        globi=mapg(ij,el)
        V_x_analytic(ij,el) = func_analyticalV_x(globi)
        V_y_analytic(ij,el) = func_analyticalV_y(globi)
        pressure_analytic(ij,el) = func_analytical_pressure(globi)
        gradUxx_analytic(ij,el) = func_analyticalGradVxx(globi)
        gradUyx_analytic(ij,el) = func_analyticalGradVyx(globi)
        gradUxy_analytic(ij,el) = func_analyticalGradVxy(globi)
        gradUyy_analytic(ij,el) = func_analyticalGradVyy(globi)    
      ENDDO
    ENDDO
  END SUBROUTINE calc_analytical_solutions
  
  SUBROUTINE calc_analytical_errors!(outV_x,outV_y,out_gradxx,out_gradyx,out_gradxy,out_gradyy,out_pressure)
! OLD VERSION IS IN SNAPSHOT DIRECTORY - PRE 03/11/12
! Calculates the error between computed and known solutions.
! All output is in local form as we will wish to use these for local integration when computing norms.
    IMPLICIT NONE
    INTEGER :: el,ij,globi
    
! Calc locally stored errors from the known solution:
    DO el=1,numelm
      DO ij=0,NP1SQM1
        globi=mapg(ij,el)
        V_x_error(ij,el) = V_x(globi) - V_x_analytic(ij,el)
        V_y_error(ij,el) = V_y(globi) - V_y_analytic(ij,el)
        
        pressure_error(ij,el) = localpressure(ij,el) - pressure_analytic(ij,el)
      
        gradUxx_error(ij,el)=localGradUxx(ij,el) - gradUxx_analytic(ij,el)
        gradUyx_error(ij,el)=localGradUyx(ij,el) - gradUyx_analytic(ij,el)
        gradUxy_error(ij,el)=localGradUxy(ij,el) - gradUxy_analytic(ij,el)
        gradUyy_error(ij,el)=localGradUyy(ij,el) - gradUyy_analytic(ij,el) ! Don't bother with GradUzz as we don't have analytical solns for this.
      ENDDO
    ENDDO
    
    IF ( param_function_choice.eq.6) THEN ! 3-d newt waters soln, uzz should be 0 (as u_r = 0)
      gradUzz_error=localGradUzz
    ELSEIF (param_function_choice.eq.7) THEN
      DO el=1,numelm
        DO ij=0,NP1SQM1
          globi=mapg(ij,el)
          Txx_error(ij,el)=localTxx(ij,el) - transient_txx(globi)
          Txy_error(ij,el)=localTxy(ij,el) - transient_txy(globi)
        ENDDO
      ENDDO
      Tyy_error=localTyy !should be 0d0
      Tzz_error=localTzz ! should be 0d0
    ELSEIF (param_function_choice.eq.8) THEN
      gradUzz_error=localGradUzz
      DO el=1,numelm
        DO ij=0,NP1SQM1
          globi=mapg(ij,el)
          Txx_error(ij,el)=localTxx(ij,el) - transient_txx(globi)
          Txy_error(ij,el)=localTxy(ij,el) - transient_txy(globi)
        ENDDO
      ENDDO
      Tyy_error=localTyy !should be 0d0
      Tzz_error=localTzz ! should be 0d0
    ENDIF
  END SUBROUTINE calc_analytical_errors


! "Known" solutions - used for validation in conjunction
! with the above routine.
  DOUBLE PRECISION FUNCTION func_analyticalV_x(i) RESULT(a)
! Known V_x
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i  
    
    IF (param_function_choice.eq.3) THEN ! known model
      a = model_soln_vel_x(i)
    ELSEIF (param_function_choice.eq.4) THEN ! known cylinder
      a = cylinder_soln_vel_x(i)
    ELSEIF (param_function_choice.eq.10) THEN ! known model 2 (Tim/Tom paper)
      a = model_soln2_vel_x(i)
    ELSEIF (param_waters) THEN ! known transient poiseuille flow
      a = transient_u(i)
    ELSE
      a = -999d0
! ADD MORE AS APPROPRIATE.
    ENDIF

  END FUNCTION func_analyticalV_x
  
  DOUBLE PRECISION FUNCTION func_analyticalGradVxx(i) RESULT(a)
! differential of V_x function above wrt x
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    
    IF     (param_function_choice.eq.3) THEN
      a = model_soln_Grad_vel_xx(i)
    ELSEIF (param_function_choice.eq.4) THEN
      a = cylinder_soln_Grad_vel_xx(i)
  ELSEIF (param_function_choice.eq.10) THEN ! known model 2 (Tim/Tom paper)
      a = model_soln2_Grad_vel_xx(i)
    ELSEIF (param_waters) THEN ! known transient poiseuille flow
      a = 0d0
    ELSE
      a = -999d0
! ADD MORE AS APPROPRIATE.
    ENDIF

  END FUNCTION func_analyticalGradVxx

  DOUBLE PRECISION FUNCTION func_analyticalGradVyx(i) RESULT(a)
! differential of V_x function above wrt y
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i

    IF     (param_function_choice.eq.3) THEN
      a = model_soln_Grad_vel_yx(i)
    ELSEIF (param_function_choice.eq.4) THEN
      a = cylinder_soln_Grad_vel_yx(i)
    ELSEIF (param_function_choice.eq.10) THEN ! known model 2 (Tim/Tom paper)
      a = model_soln2_Grad_vel_yx(i)
    ELSEIF (param_waters) THEN ! known transient poiseuille flow
      a = transient_gradUyx(i)
    ELSE
      a = -999d0
! ADD MORE AS APPROPRIATE.
    ENDIF

  END FUNCTION func_analyticalGradVyx
  
  DOUBLE PRECISION FUNCTION func_analyticalV_y(i) RESULT(a)
! Known V_y
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i    
    
    IF     (param_function_choice.eq.3) THEN
      a = model_soln_vel_y(i)
    ELSEIF (param_function_choice.eq.4) THEN
      a = cylinder_soln_vel_y(i)
    ELSEIF (param_function_choice.eq.10) THEN ! known model 2 (Tim/Tom paper)
      a = model_soln2_vel_y(i)
    ELSEIF (param_waters) THEN ! known transient poiseuille flow
      a = 0d0
    ELSE
      a = -999d0
! ADD MORE AS APPROPRIATE.
    ENDIF
    
  END FUNCTION func_analyticalV_y

  DOUBLE PRECISION FUNCTION func_analyticalGradVxy(i) RESULT(a)
! differential of V_y function above wrt x
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    
    IF     (param_function_choice.eq.3) THEN
      a = model_soln_Grad_vel_xy(i)
    ELSEIF (param_function_choice.eq.4) THEN
      a = cylinder_soln_Grad_vel_xy(i)
    ELSEIF (param_function_choice.eq.10) THEN ! known model 2 (Tim/Tom paper)
      a = model_soln2_Grad_vel_xy(i)
    ELSEIF (param_waters) THEN ! known transient poiseuille flow
      a = 0d0
    ELSE
      a = -999d0
! ADD MORE AS APPROPRIATE.
    ENDIF
  
  END FUNCTION func_analyticalGradVxy

  DOUBLE PRECISION FUNCTION func_analyticalGradVyy(i) RESULT(a)
! differential of V_y function above wrt y
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    
    IF     (param_function_choice.eq.3) THEN
      a = model_soln_Grad_vel_yy(i)
    ELSEIF (param_function_choice.eq.4) THEN
      a = cylinder_soln_Grad_vel_yy(i)
    ELSEIF (param_function_choice.eq.10) THEN ! known model 2 (Tim/Tom paper)
      a = model_soln2_Grad_vel_yy(i)
    ELSEIF (param_waters) THEN ! known transient poiseuille flow
      a = 0d0
    ELSE
      a = -999d0
! ADD MORE AS APPROPRIATE.
    ENDIF

  END FUNCTION func_analyticalGradVyy

  DOUBLE PRECISION FUNCTION func_analytical_pressure(i) RESULT(a)
! Known pressure
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i    
    
    IF     (param_function_choice.eq.3) THEN
      a = model_soln_pressure(i)
    ELSEIF (param_function_choice.eq.4) THEN
      a = cylinder_soln_pressure(i)
    ELSEIF (param_function_choice.eq.5) THEN ! known transient poiseuille flow
      a = -8d0*nodeCoord(i,1) ! MAYBE WRONG ??? !
    ELSEIF (param_function_choice.eq.6) THEN ! known transient poiseuille flow
      a = -4d0*nodeCoord(i,1) ! MAYBE WRONG ??? !
    ELSEIF (param_function_choice.eq.7) THEN ! known transient poiseuille flow
      a = -8d0*nodeCoord(i,1) ! MAYBE WRONG ??? !
    ELSEIF (param_function_choice.eq.8) THEN ! known transient poiseuille flow
      a = -4d0*nodeCoord(i,1) ! MAYBE WRONG ??? !
    ELSEIF (param_function_choice.eq.10) THEN ! known model 2 (Tim/Tom paper)
      a = model_soln2_pressure(i)
    ELSE
      a = -999d0
! ADD MORE AS APPROPRIATE.
    ENDIF

  END FUNCTION func_analytical_pressure
  
  DOUBLE PRECISION FUNCTION L2Norm_scalar_field_globalinput(x_in) RESULT(answer)
! Calculates the L2 norm of a vector x
! that is SQRT((integral over domain omega) of x^2)
    IMPLICIT NONE
    INTEGER :: i,j,jj,ij,el
    DOUBLE PRECISION, DIMENSION(nptot), INTENT(IN) :: x_in
    DOUBLE PRECISION, DIMENSION(0:NP1SQM1) :: tempvec
    DOUBLE PRECISION :: temp_sum, temp_ans

    temp_ans=0d0
    IF (coordflag.eq.0) THEN
      DO el=1,numelm
        temp_sum=0d0
        CALL veclocalrestrict(x_in,el,tempvec)
        DO j=0,N
          jj=j*NP1
          DO i=0,N
            ij=i+jj
            temp_sum = temp_sum + jac(i,j,el)*w(i)*w(j)*tempvec(ij)**2
          ENDDO
        ENDDO
        temp_ans = temp_ans + temp_sum
      ENDDO
    ELSEIF (coordflag.eq.1) THEN
      DO el=1,numelm
        temp_sum=0d0
        CALL veclocalrestrict(x_in,el,tempvec)
        DO j=0,N
          jj=j*NP1
          DO i=0,N
            ij=i+jj
            temp_sum = temp_sum + jac(i,j,el)*w(i)*w(j)*nodeCoord(mapg(ij,el),2)*tempvec(ij)**2
          ENDDO
        ENDDO
        temp_ans = temp_ans + temp_sum
      ENDDO
    ENDIF
    answer=dsqrt(temp_ans)
    RETURN
  END FUNCTION L2Norm_scalar_field_globalinput
  
  DOUBLE PRECISION FUNCTION L2Norm_scalar_field_localinput(x_in) RESULT(answer)
! Calculates the L2 norm of a vector x
! that is SQRT((integral over domain omega) of x^2)
    IMPLICIT NONE
    INTEGER :: i,j,jj,ij,el
    DOUBLE PRECISION, DIMENSION(0:NP1SQM1,numelm), INTENT(IN) :: x_in
    DOUBLE PRECISION :: temp,temp_ans

    temp_ans=0d0
    IF (coordflag.eq.0) THEN
      DO el=1,numelm
        temp=0d0
        DO j=0,N
          jj=j*NP1
          DO i=0,N
            ij=i+jj
            temp = temp + jac(i,j,el)*w(i)*w(j)*x_in(ij,el)**2
          ENDDO
        ENDDO
        temp_ans = temp_ans + temp
      ENDDO
    ELSEIF (coordflag.eq.1) THEN
      DO el=1,numelm
        temp=0d0
        DO j=0,N
          jj=j*NP1
          DO i=0,N
            ij=i+jj
            temp = temp + jac(i,j,el)*w(i)*w(j)*nodeCoord(mapg(ij,el),2)*x_in(ij,el)**2
          ENDDO
        ENDDO
        temp_ans = temp_ans + temp
      ENDDO
    ENDIF
    answer = dsqrt(temp_ans)
    RETURN
  END FUNCTION L2Norm_scalar_field_localinput
  
  DOUBLE PRECISION FUNCTION L2Norm_vector_field_globalinput(x_in,y_in) RESULT(answer)
! Calculates the L2 norm of a vector x
! that is SQRT((integral over domain omega) of x^2)
    IMPLICIT NONE
    INTEGER :: i,j,jj,ij,el
    DOUBLE PRECISION, DIMENSION(nptot), INTENT(IN) :: x_in,y_in
    DOUBLE PRECISION, DIMENSION(0:NP1SQM1) :: tempvec_x,tempvec_y
    DOUBLE PRECISION :: temp_sum, temp_ans

    temp_ans=0d0
    IF (coordflag.eq.0) THEN
      DO el=1,numelm
        temp_sum=0d0
        CALL veclocalrestrict(x_in,el,tempvec_x)
        CALL veclocalrestrict(y_in,el,tempvec_y)
        DO j=0,N
          jj=j*NP1
          DO i=0,N
            ij=i+jj
            temp_sum = temp_sum + jac(i,j,el)*w(i)*w(j)*( tempvec_x(ij)**2 + tempvec_y(ij)**2 )
          ENDDO
        ENDDO
        temp_ans = temp_ans + temp_sum
      ENDDO
    ELSEIF (coordflag.eq.1) THEN
      DO el=1,numelm
        temp_sum=0d0
        CALL veclocalrestrict(x_in,el,tempvec_x)
        CALL veclocalrestrict(y_in,el,tempvec_y)
        DO j=0,N
          jj=j*NP1
          DO i=0,N
            ij=i+jj
            temp_sum = temp_sum + jac(i,j,el)*w(i)*w(j)*( tempvec_x(ij)**2 + tempvec_y(ij)**2 )*nodeCoord(mapg(ij,el),2)
          ENDDO
        ENDDO
        temp_ans = temp_ans + temp_sum
      ENDDO
    ENDIF
    answer=dsqrt(temp_ans)
    RETURN
  END FUNCTION L2Norm_vector_field_globalinput
  
  DOUBLE PRECISION FUNCTION L2Norm_vector_field_localinput(x_in,y_in) RESULT(answer)
! Calculates the L2 norm of a vector x
! that is SQRT((integral over domain omega) of x^2)
    IMPLICIT NONE
    INTEGER :: i,j,jj,ij,el
    DOUBLE PRECISION, DIMENSION(0:NP1SQM1,numelm), INTENT(IN) :: x_in,y_in
    DOUBLE PRECISION :: temp,temp_ans

    temp_ans=0d0
    IF (coordflag.eq.0) THEN
      DO el=1,numelm
        temp=0d0
        DO j=0,N
          jj=j*NP1
          DO i=0,N
            ij=i+jj
            temp = temp + jac(i,j,el)*w(i)*w(j)*( x_in(ij,el)**2 + y_in(ij,el)**2 )
          ENDDO
        ENDDO
        temp_ans = temp_ans + temp
      ENDDO
    ELSEIF (coordflag.eq.1) THEN
      DO el=1,numelm
        temp=0d0
        DO j=0,N
          jj=j*NP1
          DO i=0,N
            ij=i+jj
            temp = temp + jac(i,j,el)*w(i)*w(j)*( x_in(ij,el)**2 + y_in(ij,el)**2 )*nodeCoord(mapg(ij,el),2)
          ENDDO
        ENDDO
        temp_ans = temp_ans + temp
      ENDDO
    ENDIF
    answer = dsqrt(temp_ans)
    RETURN
  END FUNCTION L2Norm_vector_field_localinput
  
  
  DOUBLE PRECISION FUNCTION L2Norm_tensor_field_localinput(xx_in, yx_in, xy_in, yy_in, zz_in) RESULT(answer)
! Calculates the L2 norm of a tensor (x)
! that is SQRT((integral over domain omega) of x^2)
    IMPLICIT NONE
    INTEGER :: i,j,jj,ij,el
    DOUBLE PRECISION, DIMENSION(0:NP1SQM1,numelm), INTENT(IN) :: xx_in, yx_in, xy_in, yy_in
    DOUBLE PRECISION, DIMENSION(0:NP1SQM1,numelm), INTENT(IN), OPTIONAL :: zz_in
    DOUBLE PRECISION :: temp,temp_ans

    temp_ans=0d0
    IF (coordflag.eq.0) THEN
      DO el=1,numelm
        temp=0d0
        DO j=0,N
          jj=j*NP1
          DO i=0,N
            ij=i+jj
            temp = temp + jac(i,j,el)*w(i)*w(j)*( xx_in(ij,el)**2 + yx_in(ij,el)**2 + xy_in(ij,el)**2 + yy_in(ij,el)**2 )
          ENDDO
        ENDDO
        temp_ans = temp_ans + temp
      ENDDO
    ELSEIF (coordflag.eq.1) THEN
      DO el=1,numelm
        temp=0d0
        DO j=0,N
          jj=j*NP1
          DO i=0,N
            ij=i+jj
            temp = temp + jac(i,j,el)*w(i)*w(j)*( xx_in(ij,el)**2 + yx_in(ij,el)**2 + xy_in(ij,el)**2 + yy_in(ij,el)**2 + zz_in(ij,el)**2 )*nodeCoord(mapg(ij,el),2)
          ENDDO
        ENDDO
        temp_ans = temp_ans + temp
      ENDDO
    ENDIF
!     print*,temp_ans
!     print*,sqrt(temp_ans)
!     print*,dsqrt(temp_ans)
    answer = dsqrt(temp_ans)
    
    RETURN
  END FUNCTION L2Norm_tensor_field_localinput
   
  DOUBLE PRECISION FUNCTION H1Norm_vector_field_localinput(x_in,y_in,xx_in,yx_in,xy_in,yy_in,zz_in) RESULT(answer)
! Calculates the H1 norm of a vector field with x and y components
! x_in and y_in are the values at the GLL points
! xx_in, .., yy_in are the gradient values at those GLL points
    IMPLICIT NONE
    INTEGER :: i,j,jj,ij,el
    DOUBLE PRECISION, DIMENSION(0:NP1SQM1,1:numelm), INTENT(IN) :: x_in,y_in,xx_in,yx_in,xy_in,yy_in
    DOUBLE PRECISION, DIMENSION(0:NP1SQM1,1:numelm), INTENT(IN), OPTIONAL :: zz_in
    DOUBLE PRECISION :: temp_ans, temp
    
    
    
    temp_ans=0d0
    IF (coordflag.eq.0) THEN
      DO el=1,numelm
        temp=0d0
        DO j=0,N
          jj=j*NP1
          DO i=0,N
            ij=i+jj
            temp = temp + jac(i,j,el)*w(i)*w(j)*(x_in(ij,el)**2 + y_in(ij,el)**2 + &
              xx_in(ij,el)**2 + yx_in(ij,el)**2 + &
              xy_in(ij,el)**2 + yy_in(ij,el)**2 )
          ENDDO
        ENDDO
        temp_ans = temp_ans + temp
      ENDDO
    ELSEIF (coordflag.eq.1) THEN
! NOTE WE HAVE MISSED OUT THE ZZ COMPONENT OF GRADIENT HERE. ARE NOT REALLY USING THIS ROUTINE FOR AXISYMMETRY AT THE MOMENT.
      DO el=1,numelm
        temp=0d0
        DO j=0,N
          jj=j*NP1
          DO i=0,N
            ij=i+jj
            temp = temp + jac(i,j,el)*w(i)*w(j)*(x_in(ij,el)**2 + y_in(ij,el)**2 + &
              xx_in(ij,el)**2 + yx_in(ij,el)**2 + &
              xy_in(ij,el)**2 + yy_in(ij,el)**2 +&
              zz_in(ij,el)**2 )*nodeCoord(mapg(ij,el),2)      
          ENDDO
        ENDDO
        temp_ans = temp_ans + temp
      ENDDO    
    ENDIF
    
    answer = dsqrt(temp_ans)
    RETURN
  END FUNCTION H1Norm_vector_field_localinput
  
  SUBROUTINE calc_stopping_criteria(out,drag_out)
! Returns the L2Norm of the difference between field variables computed at the current and previous timesteps.
! The answer is normalised by the L2norm of the current time solution
    IMPLICIT NONE
    DOUBLE PRECISION :: top,bottom,&
      l2_current_vel,l2_current_pressure,l2_current_stress, &
      l2_diff_vel, l2_diff_pressure, l2_diff_stress, &
      tempvel_x(nptot), tempvel_y(nptot), &
      temp_pressure(0:NP1SQM1,numelm), &
      tempstress_xx(0:NP1SQM1,numelm), tempstress_xy(0:NP1SQM1,numelm), &
      tempstress_yy(0:NP1SQM1,numelm), tempstress_zz(0:NP1SQM1,numelm)
    DOUBLE PRECISION, INTENT(OUT) :: out,drag_out
 
    IF (param_problem_choice.lt.20) THEN
! Stokes - no timestep so stop it now !
      out=0d0
      drag_out=0d0
      RETURN
    ELSEIF (mod(param_problem_choice,10).eq.3.or.mod(param_problem_choice,10).eq.4) THEN
! Newtonian/Viscoelastic start-up Poiseuille Flow - Want to stop it after time limit only, so always return 1 for each.
! Time limit is adjust in main.
      out=1d0
      drag_out=1d0
      RETURN
    ELSEIF (param_problem_choice.eq.41.or.param_problem_choice.eq.42 &
      .or.param_problem_choice.eq.51.or.param_problem_choice.eq.52) THEN ! moving mesh, not interested in change of fields, only velocity and drag for sphere.
      out=abs((V_sphere-V_sphereNm1)/V_sphere)
      drag_out=abs((drag-dragNm1)/drag)
    ELSE
    
      l2_current_vel = L2Norm_vector_field_globalinput(V_x,V_y)
      tempvel_x = V_x - V_xNm1
      tempvel_y = V_y - V_yNm1
      l2_diff_vel = L2Norm_vector_field_globalinput(tempvel_x,tempvel_y)
    
      l2_current_pressure = L2Norm_scalar_field_localinput(localpressure)
      temp_pressure = localpressure - localpressureNm1
      l2_diff_pressure = L2Norm_scalar_field_localinput(temp_pressure)

      top = l2_diff_vel**2 + l2_diff_pressure**2
      bottom = l2_current_vel**2 + l2_current_pressure**2
      IF (param_beta.ne.1d0) THEN
        l2_current_stress = L2Norm_tensor_field_localinput(localTxx,localTxy,localTxy,localTyy,localTzz)
        tempstress_xx = localTxx - localTxxNm1
        tempstress_xy = localTxy - localTxyNm1
        tempstress_yy = localTyy - localTyyNm1
        tempstress_zz = localTzz - localTzzNm1
        l2_diff_stress = L2Norm_tensor_field_localinput(tempstress_xx,tempstress_xy,tempstress_xy,tempstress_yy,tempstress_zz)

        top = top + l2_diff_stress**2
        bottom = bottom + l2_current_stress**2
      ENDIF

      out = dsqrt(top/bottom)

      drag_out=abs((drag-dragNm1)/drag)

    ENDIF 

  END SUBROUTINE calc_stopping_criteria
  

  SUBROUTINE run_error_analysis
    IMPLICIT NONE
    

    CALL calc_analytical_solutions
    CALL calc_analytical_errors

! Always want velocity error norm.
    
! Now calculate appropriate error measures depending on problem.
    IF (param_function_choice.eq.3.or.param_function_choice.eq.4.or.param_function_choice.eq.10) THEN
! Stokes model soln. h1-norm of velocity error. L2-norm of pressure error.
! AND
! Stokes 2d cylinder. h1-norm of velocity error. L2-norm of pressure error.
! AND
! Stokes model soln2 For Tim/Tom paper. l2norm and H1-norm of velocity error. L2-norm of pressure error.
      L2norm_vel_err = L2Norm_vector_field_globalinput(V_x_error,V_y_error)
      H1norm_vel_err = H1Norm_vector_field_localinput(V_x_error,V_y_error,gradUxx_error,gradUyx_error,gradUxy_error,gradUyy_error)
      L2norm_press_err = L2Norm_scalar_field_localinput(pressure_error)
    
    ELSEIF (param_function_choice.eq.5) THEN
! 2D Newtonian. H1-norm of velocity error.
! Point error of V_x and GradVyx.
      L2norm_vel_err = L2Norm_vector_field_globalinput(V_x_error,V_y_error)
      H1norm_vel_err = H1Norm_vector_field_localinput(V_x_error,V_y_error,gradUxx_error,gradUyx_error,gradUxy_error,gradUyy_error)
      transient_U_error = V_x(transient_velocity_testnode) - transient_u(transient_velocity_testnode)
      transient_gradUyx_error = gradUyx_error(transient_stress_testnode,transient_stress_testelement)
      L2norm_press_err = L2Norm_scalar_field_localinput(pressure_error)
    
    ELSEIF (param_function_choice.eq.6) THEN
! 3D Axisymm Newtonian. H1-norm of velocity error.
! Point error of V_x and GradVyx.
      L2norm_vel_err = L2Norm_vector_field_globalinput(V_x_error,V_y_error)
      H1norm_vel_err = H1Norm_vector_field_localinput(V_x_error,V_y_error,gradUxx_error,gradUyx_error,gradUxy_error,gradUyy_error,gradUzz_error)
      transient_U_error = V_x(transient_velocity_testnode) - transient_u(transient_velocity_testnode)
      transient_gradUyx_error = gradUyx_error(transient_stress_testnode,transient_stress_testelement)
      L2norm_press_err = L2Norm_scalar_field_localinput(pressure_error)
    
    ELSEIF (param_function_choice.eq.7) THEN
! 2D Viscoelastic. H1-norm of velocity error.  L2norm of elastic stress.
! Point error of V_x, GradVyx, Tao_xx and Tao_xy
      L2norm_vel_err = L2Norm_vector_field_globalinput(V_x_error,V_y_error)
      H1norm_vel_err = H1Norm_vector_field_localinput(V_x_error,V_y_error,gradUxx_error,gradUyx_error,gradUxy_error,gradUyy_error)
      L2norm_stress_err = L2Norm_tensor_field_localinput(Txx_error,Txy_error,Txy_error,Tyy_error)
      L2norm_press_err = L2Norm_scalar_field_localinput(pressure_error)
      
      transient_U_error = V_x(transient_velocity_testnode) - transient_u(transient_velocity_testnode)
      transient_gradUyx_error = gradUyx_error(transient_stress_testnode,transient_stress_testelement)    
      transient_Txx_error = Txx_error(transient_stress_testnode,transient_stress_testelement)
      transient_Txy_error = Txy_error(transient_stress_testnode,transient_stress_testelement)
    
    
    ELSEIF (param_function_choice.eq.8) THEN
! 3D Axisymm Viscoelastic. H1-norm of velocity error. L2norm of elastic stress.
! Point error of V_x, GradVyx, Tao_xx and Tao_xy
      L2norm_vel_err = L2Norm_vector_field_globalinput(V_x_error,V_y_error)
      H1norm_vel_err = H1Norm_vector_field_localinput(V_x_error,V_y_error,gradUxx_error,gradUyx_error,gradUxy_error,gradUyy_error,gradUzz_error)
      L2norm_stress_err = L2Norm_tensor_field_localinput(Txx_error,Txy_error,Txy_error,Tyy_error,Tzz_error)
      L2norm_press_err = L2Norm_scalar_field_localinput(pressure_error)
      
      transient_U_error = V_x(transient_velocity_testnode) - transient_u(transient_velocity_testnode)
      transient_gradUyx_error = gradUyx_error(transient_stress_testnode,transient_stress_testelement)
      transient_Txx_error = Txx_error(transient_stress_testnode,transient_stress_testelement)
      transient_Txy_error = Txy_error(transient_stress_testnode,transient_stress_testelement)

    ENDIF
  END SUBROUTINE run_error_analysis

END MODULE result_analysis


