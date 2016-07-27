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

MODULE initial_conditions_module
  USE shared_data
  USE functions_module
  USE SEM_module
  USE waters_solution_module
  IMPLICIT NONE
  CONTAINS
  
  SUBROUTINE apply_initial_conditions
! Need to have called setup_problem_style first
    IMPLICIT NONE
    INTEGER :: i,el
    
    IF (param_waters) THEN ! V/E Transient solution - need the constants Cxx and Cxy
      CALL initialise_waters_solution
! Assign the test points for the Transient solution error.

      IF (param_function_choice.eq.6.or.param_function_choice.eq.8) THEN
! This is the penultimate node in the x-direction for both, and the node at the top of the channel for stress, and the bottom of the channel for velocity.
	CALL find_transient_error_points_axisymm
      ELSE
! This is the penultimate node in the x-direction for both, and the node nearest y-coord 0.5 for velocity and nearest 0.0 for stress.
	CALL find_transient_error_points
      ENDIF
    ENDIF

    IF (param_problem_choice.lt.20) THEN ! Stokes
! time not included so not required.    
      V_sphere=1d0
      
      timeNm2 = 0d0
      timeNm1 = 0d0
      timeN = 0d0
      mesh_velocity = 0d0
      mesh_velocityNm1 = 0d0
      mesh_velocityNm2 = 0d0
      
      V_x=0d0
      V_xNm1=0d0

      V_y=0d0
      V_yNm1=0d0

      pressure=0d0
      localpressure=0d0
      localpressureNm1=0d0
      
      localGradUxx=0d0
      localGradUxy=0d0
      localGradUyx=0d0
      localGradUyy=0d0
      localGradUzz=0d0
      
      boundary_x=0d0
      boundary_y=0d0
    ELSEIF (param_problem_choice.lt.30) THEN ! Newtonian (fixed mesh)
    
      V_sphere = 1d0
      
      timeNm2 = 0d0
      timeNm1 = 0d0
      timeN = 0d0!deltat
      mesh_velocity = 0d0
      mesh_velocityNm1 = 0d0
      mesh_velocityNm2 = 0d0
      
      V_x=0d0
      V_xNm1=0d0 
      V_xNm2=0d0
      
      V_y=0d0
      V_yNm1=0d0
      V_yNm2=0d0
      
      pressure=0d0
      localpressure=0d0
      localpressureNm1=0d0
      
      localGradUxx=0d0
      localGradUxy=0d0
      localGradUyx=0d0
      localGradUyy=0d0
      localGradUzz=0d0
      
      boundary_x=0d0
      boundary_y=0d0
      
    ELSEIF (param_problem_choice.lt.40) THEN ! Viscoelastic (fixed mesh)
    
      V_sphere = 1d0
      timeNm2 = 0d0
      timeNm1 = 0d0
      timeN = 0d0!deltat
      
      mesh_velocity = 0d0
      mesh_velocityNm1 = 0d0
      mesh_velocityNm2 = 0d0
      
      V_x=0d0
      V_xNm1=0d0
      V_xNm2=0d0
      
      V_y=0d0
      V_yNm1=0d0
      V_yNm2=0d0
      
      pressure=0d0
      localpressure=0d0
      localpressureNm1=0d0
      
      localGradUxx=0d0
      localGradUxy=0d0
      localGradUyx=0d0
      localGradUyy=0d0
      localGradUzz=0d0

      boundary_x=0d0
      boundary_y=0d0
      
      stress_cont_to_stokes_xNm1=0d0
      stress_cont_to_stokes_yNm1=0d0
      stress_cont_to_stokes_x=0d0
      stress_cont_to_stokes_y=0d0

      devss_soln_xx=0d0
      devss_soln_xy=0d0
      devss_soln_yx=0d0
      devss_soln_yy=0d0
      devss_soln_zz=0d0
      
      localTxx=0d0
      localTxy=0d0
      localTyy=0d0
      localTzz=0d0
      localTxxNm1=0d0
      localTxyNm1=0d0
      localTyyNm1=0d0
      localTzzNm1=0d0
      localTxxNm2=0d0
      localTxyNm2=0d0
      localTyyNm2=0d0
      localTzzNm2=0d0
      

      
    ELSEIF (param_problem_choice.lt.50) THEN ! Newtonian (moving mesh)
    
      V_sphere=1d-5
      V_sphereNm1=1d-5
      V_sphereNm2=1d-5
      
      timeNm2 = 0d0
      timeNm1 = 0d0
      timeN = 0d0!deltat
      
      mesh_velocity = 0d0
      mesh_velocityNm1 = 0d0
      mesh_velocityNm2 = 0d0
      
      V_x=0d0
      V_xNm1=0d0
      V_xNm2=0d0
      
      V_y=0d0
      V_yNm1=0d0
      V_yNm2=0d0
      
      pressure=0d0
      localpressure=0d0
      localpressureNm1=0d0
      
      localGradUxx=0d0
      localGradUxy=0d0
      localGradUyx=0d0
      localGradUyy=0d0
      localGradUzz=0d0
      
      boundary_x=0d0
      boundary_y=0d0
      
    ELSEIF (param_problem_choice.lt.60) THEN ! Viscoelastic (moving mesh)
    
      V_sphere=1d-5
      V_sphereNm1=1d-5
      V_sphereNm2=1d-5
      
      mesh_velocity = 0d0
      mesh_velocityNm1 = 0d0
      mesh_velocityNm2 = 0d0
      
      timeNm2 = 0d0
      timeNm1 = 0d0
      timeN = 0d0!deltat
      
      V_x=0d0
      V_xNm1=0d0
      V_xNm2=0d0
      
      V_y=0d0
      V_yNm1=0d0
      V_yNm2=0d0
      
      pressure=0d0
      localpressure=0d0
      localpressureNm1=0d0
      
      localGradUxx=0d0
      localGradUxy=0d0
      localGradUyx=0d0
      localGradUyy=0d0
      localGradUzz=0d0
      
      boundary_x=0d0
      boundary_y=0d0
      
      stress_cont_to_stokes_xNm1=0d0
      stress_cont_to_stokes_yNm1=0d0
      stress_cont_to_stokes_x=0d0
      stress_cont_to_stokes_y=0d0
      
      devss_soln_xx=0d0
      devss_soln_xy=0d0
      devss_soln_yx=0d0
      devss_soln_yy=0d0
      devss_soln_zz=0d0
      
      localTxx=0d0
      localTxy=0d0
      localTyy=0d0
      localTzz=0d0
      localTxxNm1=0d0
      localTxyNm1=0d0
      localTyyNm1=0d0
      localTzzNm1=0d0
      localTxxNm2=0d0
      localTxyNm2=0d0
      localTyyNm2=0d0
      localTzzNm2=0d0
      

      
    ENDIF
    


    
    
! Steady Poiseuille Flow for pipe from [-1,1]:
!     DO i=1,nptot
! 	V_xNm1(i) = 1d0-nodeCoord(i,2)**2!0.25d0*(4d0-nodeCoord(i,2)**2)0.25d0*(4d0-nodeCoord(i,2)**2)
! 	V_x(i) = 1d0-nodeCoord(i,2)**2!0.25d0*(4d0-nodeCoord(i,2)**2)
!     ENDDO


! Transient ICs (Waters & King solution):
! Newtonian:
!     DO i=1,nptot
! 	V_xNm1(i) = calcWaters_newtonian_U(nodeCoord(i,1),nodeCoord(i,2),timeNm1)
! 	V_x(i) = calcWaters_newtonian_U(nodeCoord(i,1),nodeCoord(i,2),timeN)
!     ENDDO

! Oldroyd B:
!     DO i=1,nptot
!     
!       V_xNm1(i) = calcWaters_analytical_U(nodeCoord(i,1),nodeCoord(i,2),timeNm1)
!       V_x(i) = calcWaters_analytical_U(nodeCoord(i,1),nodeCoord(i,2),timeN)
!     
!       TxxNm1(i) = calcWaters_analytical_Txx(nodeCoord(i,1),nodeCoord(i,2),timeNm1)
!       TxyNm1(i) = calcWaters_analytical_Txy(nodeCoord(i,1),nodeCoord(i,2),timeNm1)
!     
!       Txx(i) = calcWaters_analytical_Txx(nodeCoord(i,1),nodeCoord(i,2),timeN)
!       Txy(i) = calcWaters_analytical_Txy(nodeCoord(i,1),nodeCoord(i,2),timeN)
!    
!       uxy(i) = calcWaters_analytical_dUdy(nodeCoord(i,1),nodeCoord(i,2),timeN)
! 
!     ENDDO

!     CALL calc_waters_solution(nodeCoord(1:nptot,2),timeNm1,V_xNm1,TxxNm1,TxyNm1)
!     CALL calc_waters_solution(nodeCoord(1:nptot,2),timeN,V_x,Txx,Txy)

!     CALL calc_gradient_of_U
!     CALL calc_local_gradient_of_U
! !     DO i=1,nptot
! !       uxy(i) = calcWaters_analytical_dUdy(nodeCoord(i,1),nodeCoord(i,2),timeN)
! !     ENDDO
! !     DO el=1,numelm
! !     DO i=0,NP1SQ-1
! !       localTxxNm1(i,el)=TxxNm1(mapg(i,el))
! !       localTxyNm1(i,el)=TxyNm1(mapg(i,el))
! !       localTyyNm1(i,el)=TyyNm1(mapg(i,el))
! !       localTzzNm1(i,el)=TzzNm1(mapg(i,el))
! !       localTxx(i,el)=Txx(mapg(i,el))
! !       localTxy(i,el)=Txy(mapg(i,el))
! !       localTyy(i,el)=Tyy(mapg(i,el))
! !       localTzz(i,el)=Tzz(mapg(i,el))
! !       localGradUxx(i,el)=0d0
! !       localGradUyx(i,el)=uxy(mapg(i,el))
! !       localGradUxy(i,el)=0d0
! !       localGradUyy(i,el)=0d0
! !       localGradUzz(i,el)=0d0
! !     ENDDO
! !     ENDDO
  
  END SUBROUTINE apply_initial_conditions
  
  SUBROUTINE setup_problem_style
    IMPLICIT NONE
    

    
    IF (param_problem_choice.eq.11.or.param_problem_choice.eq.21.or.param_problem_choice.eq.31) THEN ! Poiseuille flow past fixed cylinder.
      IF (coordflag.eq.1) THEN
	write(*,*) 'ERROR: coordflag inconsistent with problem style'
	STOP
      ENDIF
      param_function_choice = 1
    ELSEIF (param_problem_choice.eq.12.or.param_problem_choice.eq.22.or.param_problem_choice.eq.32) THEN ! Uniform flow past fixed sphere.
      IF (coordflag.eq.0) THEN
	write(*,*) 'ERROR: coordflag inconsistent with problem style'
	STOP
      ENDIF
      param_function_choice = 2
    ELSEIF (param_problem_choice.eq.13) THEN ! Stokes Model solution
      IF (coordflag.eq.1) THEN
	write(*,*) 'ERROR: coordflag inconsistent with problem style'
	STOP
      ENDIF
      param_function_choice = 3
      param_error=.true.
    ELSEIF (param_problem_choice.eq.14) THEN ! Stokes cylinder solution
      IF (coordflag.eq.1) THEN
	write(*,*) 'ERROR: coordflag inconsistent with problem style'
	STOP
      ENDIF
      param_function_choice = 4
      param_error=.true.
! THESE ARE YET TO BE MADE FINAL:
    ELSEIF (param_problem_choice.eq.23) THEN ! Newtonian 2-D Transient Waters solution
      IF (coordflag.eq.1) THEN
	write(*,*) 'ERROR: coordflag inconsistent with problem style'
	STOP
      ENDIF
      param_function_choice = 5
      param_error=.true.
      param_waters=.true.
    ELSEIF (param_problem_choice.eq.24) THEN ! Newtonian 3-D Axisymmetric Transient Waters solution
      IF (coordflag.eq.0) THEN
	write(*,*) 'ERROR: coordflag inconsistent with problem style'
	STOP
      ENDIF
      param_function_choice = 6
      param_error=.true.
      param_waters=.true.
    ELSEIF (param_problem_choice.eq.33) THEN ! Viscoelastic 2-D Transient Waters solution
      IF (coordflag.eq.1) THEN
	write(*,*) 'ERROR: coordflag inconsistent with problem style'
	STOP
      ENDIF
      param_function_choice = 7
      param_error=.true.
      param_waters=.true.
    ELSEIF (param_problem_choice.eq.34) THEN ! Viscoelastic 3-D Axisymmetric Transient Waters solution
      IF (coordflag.eq.0) THEN
	write(*,*) 'ERROR: coordflag inconsistent with problem style'
	STOP
      ENDIF
      param_function_choice = 8
      param_error=.true.
      param_waters=.true.
    ELSEIF (param_problem_choice.eq.41.or.param_problem_choice.eq.42.or.&  ! Moving mesh Newtonian/Viscoelastic 2-D
	    param_problem_choice.eq.51.or.param_problem_choice.eq.52) THEN ! Moving mesh Newtonian/Viscoelastic 3-D Axisymmetric.

      param_function_choice = 9

    ELSE
      write(*,*) 'ERROR: No problem choice made'
      write(*,*) 'param_function_choice is ', param_function_choice
      STOP
    ENDIF
    
! Check movingmeshflag is correct:
    IF (param_problem_choice.lt.40) THEN
      IF (movingmeshflag.eq.1) THEN
	write(*,*) 'ERROR: movingmeshflag inconsistent with problem style'
	STOP
      ENDIF
    ELSEIF (param_problem_choice.lt.60) THEN
      IF (movingmeshflag.eq.0) THEN
	write(*,*) 'ERROR: movingmeshflag inconsistent with problem style'
	STOP
      ENDIF
    ENDIF
    
! Check fluid parameters match problem choices:
    IF (param_problem_choice.lt.20) THEN ! Stokes
      IF (Re.gt.0d0) THEN
	write(*,*) 'ERROR: Re inconsistent with problem style'
	STOP
      ENDIF
      IF (We.gt.0d0) THEN
	write(*,*) 'ERROR: We inconsistent with problem style'
	STOP
      ENDIF
      IF (param_beta.ne.1d0) THEN
	write(*,*) 'ERROR: param_beta inconsistent with problem style'
	STOP
      ENDIF    

    ELSEIF (param_problem_choice.lt.30) THEN ! Newtonian
      IF (Re.lt.1d-6) THEN
	write(*,*) 'ERROR: Re inconsistent with problem style'
	STOP
      ENDIF
      IF (We.gt.0d0) THEN
	write(*,*) 'ERROR: We inconsistent with problem style'
	STOP
      ENDIF
      IF (param_beta.ne.1d0) THEN
	write(*,*) 'ERROR: param_beta inconsistent with problem style'
	STOP
      ENDIF 
   
    ELSEIF (param_problem_choice.lt.40) THEN ! Viscoelastic
      IF (We.lt.1d-6) THEN
	write(*,*) 'ERROR: We inconsistent with problem style'
	STOP
      ENDIF
      IF (param_beta.ge.1d0.or.param_beta.le.0d0) THEN
	write(*,*) 'ERROR: param_beta inconsistent with problem style'
	STOP
      ENDIF 
   
    ELSEIF (param_problem_choice.lt.50) THEN ! Moving Newtonian
      IF (Re.lt.1d-6) THEN
	write(*,*) 'ERROR: Re inconsistent with problem style'
	STOP
      ENDIF
      IF (We.gt.0d0) THEN
	write(*,*) 'ERROR: We inconsistent with problem style'
	STOP
      ENDIF
      IF (param_beta.ne.1d0) THEN
	write(*,*) 'ERROR: param_beta inconsistent with problem style'
	STOP
      ENDIF 
   
    ELSEIF (param_problem_choice.lt.60) THEN ! Moving Viscoelastic
      IF (Re.lt.1d-6) THEN
	write(*,*) 'ERROR: Re inconsistent with problem style'
	STOP
      ENDIF
      IF (We.lt.1d-6) THEN
	write(*,*) 'ERROR: We inconsistent with problem style'
	STOP
      ENDIF
      IF (param_beta.ge.1d0.or.param_beta.le.0d0) THEN
	write(*,*) 'ERROR: param_beta inconsistent with problem style'
	STOP
      ENDIF 
    
    ENDIF
  END SUBROUTINE setup_problem_style
  
  
  SUBROUTINE find_transient_error_points
    IMPLICIT NONE
    INTEGER :: el,i,j,xi_i,eta_j,tempij,&
		chan_end_glob_node,chan_end_loc_vert
    DOUBLE PRECISION :: channel_end,xi_out,eta_out,dist_xi(0:N),dist_eta(0:N),&
			temp_y,temp_xi,temp_eta,test_node_xcoord
    
! find x-coord of the end of channel: which is the largest x co-ord of vertex nodes.
    channel_end=-999d0
    temp_y=999d0
    DO i=1,numnp
      IF (nodeCoord(i,1).gt.channel_end.and.nodeCoord(i,2).lt.temp_y) THEN
	channel_end=nodeCoord(i,1)
	temp_y=nodeCoord(i,2)
	chan_end_glob_node=i
	ENDIF
    ENDDO
    IF (channel_end.lt.-998d0) THEN
      print*, 'ERROR: computing end of channel failed'
      STOP
    ENDIF
! find element and local vertex of the node which lies on the end of the channel (and is at the "bottom")
    DO el=1,numelm
      DO i=1,4
	IF (node(el,i).eq.chan_end_glob_node) THEN
	  transient_stress_testelement=el
	  chan_end_loc_vert=i
	ENDIF
      ENDDO
    ENDDO
! Use the vertex to work out the local numbering and assign the transient stress test point.
! The error checking point is the penultimate GLL point..
    IF (chan_end_loc_vert.eq.1) THEN
      transient_stress_testnode=NM1 ! GLL: (N-1,0)
    ELSEIF (chan_end_loc_vert.eq.2) THEN
      transient_stress_testnode=NP1 ! GLL: (0,1)
    ELSEIF (chan_end_loc_vert.eq.3) THEN
      transient_stress_testnode=1+N*NP1  ! GLL: (1,N)
    ELSEIF (chan_end_loc_vert.eq.4) THEN
      transient_stress_testnode=N+NM1*NP1  ! GLL: (N,N-1)
    ELSE
      print*,'ERROR: computing chan_end_loc_vert failed'
      STOP
    ENDIF

! Now need the x-coord of this point so we can find the corresponding velocity test node (@y=0.5 above)
    test_node_xcoord=nodeCoord(mapg(transient_stress_testnode,transient_stress_testelement),1)

! Find the element containing this point and local co-ordinates within the element.
    CALL map_x_y_to_xi_eta(test_node_xcoord,0.5d0, el, xi_out, eta_out)
     
! Now find closest GLL point to this local co-ordinate.
    DO i=0,N
      dist_xi(i)=abs((gl(i)-xi_out))
      dist_eta(i)=abs((gl(i)-eta_out))
    ENDDO
    
    temp_xi=dist_xi(0)
    temp_eta=dist_eta(0)
    xi_i=0
    eta_j=0
    DO i=1,N
      
      IF (dist_xi(i).lt.temp_xi) THEN
	temp_xi=dist_xi(i)
	xi_i=i
      ENDIF
      IF (dist_eta(i).lt.temp_eta) THEN
	temp_eta=dist_eta(i)
	eta_j=i
      ENDIF
    ENDDO
    
    tempij=xi_i+eta_j*NP1
    
! Now translate to the global point
    transient_velocity_testnode = mapg(tempij,el)

      
      
  END SUBROUTINE find_transient_error_points
  
  SUBROUTINE find_transient_error_points_axisymm
    IMPLICIT NONE
    INTEGER :: el,i,temp_vel_node,temp_vel_element,&
		chan_end_glob_node,chan_end_loc_vert
    DOUBLE PRECISION :: channel_end,temp_y

! Stress node:
! Find point which lies at end (largest positive x) and top (largest positive y) of the channel:
    channel_end=-999d0
    temp_y=-999d0
    DO i=1,numnp
      IF (nodeCoord(i,1).gt.channel_end.and.nodeCoord(i,2).gt.temp_y) THEN
	channel_end=nodeCoord(i,1)
	temp_y=nodeCoord(i,2)
	chan_end_glob_node=i
	ENDIF
    ENDDO
    IF (channel_end.lt.-998d0.or.temp_y.lt.-998d0) THEN
      print*, 'ERROR: computing end of channel failed (stress)'
      STOP
    ENDIF
! Find the local element and local vertex within the element which corresponds to the point found.
! will only be a single element which this lies in as it is in a global corner.
    DO el=1,numelm
      DO i=1,4
	IF (node(el,i).eq.chan_end_glob_node) THEN
	  transient_stress_testelement=el
	  chan_end_loc_vert=i
	ENDIF
      ENDDO
    ENDDO
! Use the vertex to work out the local numbering and assign the transient stress test point.
! The error checking point is the penultimate GLL point on the top edge.
    IF (chan_end_loc_vert.eq.1) THEN
      transient_stress_testnode=1 ! GLL: (1,0)
    ELSEIF (chan_end_loc_vert.eq.2) THEN
      transient_stress_testnode=N+1*NP1 ! GLL: (N,1)
    ELSEIF (chan_end_loc_vert.eq.3) THEN
      transient_stress_testnode=NM1 + N*NP1 ! GLL: (N-1,N)
    ELSEIF (chan_end_loc_vert.eq.4) THEN
      transient_stress_testnode=NM1*NP1 ! GLL: (0,N-1)
    ELSE
      print*,'ERROR: computing chan_end_loc_vert failed (stress)'
      STOP
    ENDIF
! ===============================================================================================
! Velocity Node:
! Find point which lies at end (largest positive x) and top (smallest positive y) of the channel:
    channel_end=-999d0
    temp_y=999d0
    DO i=1,numnp
      IF (nodeCoord(i,1).gt.channel_end.and.nodeCoord(i,2).lt.temp_y) THEN
	channel_end=nodeCoord(i,1)
	temp_y=nodeCoord(i,2)
	chan_end_glob_node=i
	ENDIF
    ENDDO
    IF (channel_end.lt.-998d0.or.temp_y.gt.998d0) THEN
      print*, 'ERROR: computing end of channel failed (velocity)'
      STOP
    ENDIF
! Find the local element and local vertex within the element which corresponds to the point found.
! will only be a single element which this lies in as it is in a global corner.
    DO el=1,numelm
      DO i=1,4
	IF (node(el,i).eq.chan_end_glob_node) THEN
	  temp_vel_element=el
	  chan_end_loc_vert=i
	ENDIF
      ENDDO
    ENDDO
! Use the vertex to work out the local numbering and assign the transient velocity test node.
! The error checking point is the penultimate GLL point on the bottom edge.
    IF (chan_end_loc_vert.eq.1) THEN
      temp_vel_node=NM1 ! GLL: (N-1,0)
    ELSEIF (chan_end_loc_vert.eq.2) THEN
      temp_vel_node=NP1 ! GLL: (0,1)
    ELSEIF (chan_end_loc_vert.eq.3) THEN
      temp_vel_node=1+N*NP1  ! GLL: (1,N)
    ELSEIF (chan_end_loc_vert.eq.4) THEN
      temp_vel_node=N+NM1*NP1  ! GLL: (N,N-1)
    ELSE
      print*,'ERROR: computing chan_end_loc_vert failed (velocity)'
      STOP
    ENDIF
    
! Now translate to the global point
    transient_velocity_testnode = mapg(temp_vel_node,temp_vel_element)
      
  END SUBROUTINE find_transient_error_points_axisymm


  SUBROUTINE initialise_2nd_timestep
! Use a 2nd order scheme (modified version of pg 156 in Tim's book) to
! predict the 2nd timestep required for the first iteration of our main scheme.
    IMPLICIT NONE
    INTEGER :: i,j,el
    DOUBLE PRECISION :: temp
    DOUBLE PRECISION, DIMENSION(nptot) :: u_half,v_half,tau_xx_half,tau_xy_half,tau_yy_half,tau_zz_half
    DOUBLE PRECISION, DIMENSION(npint) :: old_pressure

! Step 1: calculate u and v at the half timestep
  
    

  
  END SUBROUTINE initialise_2nd_timestep

END MODULE initial_conditions_module