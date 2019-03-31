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

PROGRAM main
  USE constants
  USE shared_data
  USE geometry_module
  USE SEM_module
  USE boundary_module
  USE initial_conditions_module
  USE pardiso_solver
  USE OIFS_module
  USE viscoelastic_module
  USE fene_p_mp_module
  USE devss_module
  USE IO_module
  USE result_analysis
  USE movingSphere_module
  
  
  IMPLICIT NONE
  INTEGER :: edge,i,j,ij,k,l,kl,ijkl,el,tempint,minresconv=0,internalij,kk,ll,&
       timestep,numtimesteps,rowcount,printoutcount,stressnode,velnode,localstressnode,printel,meh,&
       print_threshold,ierror
  DOUBLE PRECISION :: time_limit,test,msqe,temp1,temp2,cputime_temp1,cputime_temp2,cputime_section1,cputime_section2,&
          V_xH1norm,V_yH1norm,pressure_L2norm,stopping_tol,stopcheck,xi_out,eta_out,drag_criteria
  LOGICAL :: printed_to_tecplot

! Read in values (This will assign required memory)
  CALL read_input

! Create local to global map
  CALL create_global_map
! assign remaining memory
  CALL assignMem2

! Create dimension and map to deal with dirichlet nodes
  CALL map_non_dirichlet_edge_nodes
  
! parameter to set the boundary conditions and result analysis:
! GUIDE:
! 1x = Stokes (Fixed mesh)
! 2x = Newtonian (Fixed mesh)
! 3x = Viscoelastic (Fixed mesh)
! 4x = Newtonian (Moving mesh)
! 5x = Viscoelastic (Moving mesh)
!
! x1 = Poiseuille flow (past cylinder)
! x2 = Uniform flow (past sphere)
! x3 = Stokes: model solution | Newt/Visco: transient Poiseuille (2-D)
! x4 = Stokes: known cylinder solution | Newt/Visco: transient Poiseuille (Axisymm 3-D)
! x9 = Viscoelastic only (39): FENE-P-MP Model uniform flow past sphere.

  CALL setup_problem_style  
  
  
! parameter to set order of temporal scheme (default 2)
! Haven't streamlined to allow for other orders in all routines.. wouldn't be that difficult using sums up to param_time_order-1
  CALL setup_time_constants

  stopping_tol = 1d-6*deltat ! should be timestep-dependent ?? 
  
! Set the threshold iterations between printing to tecplot.
! The aim is to print 10 per time unit for a fixed mesh and
! 20 per time unit for a moving mesh simulation.
  IF (movingmeshflag.eq.1) THEN
    print_threshold = INT(5d-2/deltat)
  ELSE
    print_threshold = INT(1d-1/deltat)
  ENDIF
 
  RK4_timesteps=8
  h_RK4=deltat/dfloat(RK4_timesteps)


! SEM SETUP:
! calc time-indepedent/geometry-indepedent parts of SEM:
  CALL initialiseSEM
! Geometry Setup:
  CALL initialise_geometry
 
! Find co-ords of all points for the initial mesh (requires SEM routines first):
  CALL calcNodeCoords
! In case of moving mesh we are starting from rest:
  nodeCoordNm1=nodeCoord
! Calculate all geometry related parts
  CALL update_geometry

! Calculate time-dependent/geometry parts of SEM
  CALL initialiseSEMarrays
  CALL calc_domain_area
 
  
! ALL INITIAL CONDITIONS ARE DONE HERE:
!!! CHANGES FOR PROBLEM CHOICE HERE !!!
  CALL apply_initial_conditions
  
  IF (param_waters) THEN
    IF (param_problem_choice.lt.30) THEN
! Newtonian:
      time_limit = 10d0
    ELSE
! Viscoelastic:
      time_limit = 40d0 
    ENDIF
  ELSEIF (movingmeshflag.eq.1) THEN
    time_limit=20d0
  ELSE
! All others:
    time_limit = 40d0
  ENDIF
  

! INITIALISE THE DEVSS SCHEME:
  IF (param_beta_s.gt.0d0) THEN
    CALL initialise_devss 
  ENDIF
 

! Open tecplot output file & write the ICs set.
! MODIFIED TO PRINT OUT FOR BOTH MOVING AND NON-MOVING
  IF (enable_interim_output_to_file) THEN
    OPEN(tecplot_output_fileid,FILE=tecplot_output_filename,IOSTAT=ierror)
    CALL output_to_tecplot
  ENDIF

  printoutcount=print_threshold ! Forces code to print first timestep.

  
  numtimesteps=0
  
! INITIALISE THE SOLVER:
  CALL initialise_pardiso_stokes

  CALL cpu_time(cputime_temp1)
  cputime_initialise=cputime_temp1
  
  CALL cpu_time(cputime_temp2)
  DO ! Main timestepping loop
      

! Newtonian part, calculates time(N+1) contribution:
    IF (Re.gt.1d-5) THEN
      CALL cpu_time(cputime_section1)
      CALL applyOIFS_to_stokes
      CALL cpu_time(cputime_section2)
      cputime_newtonian = cputime_section2-cputime_section1
    ENDIF

 ! Time-increment:
    numtimesteps = numtimesteps + 1
    timeNm2 = timeNm1
    timeNm1 = timeN
    timeN = (dfloat(numtimesteps))*deltat
    
    IF (param_waters) THEN
      CALL calc_waters_solution
    ENDIF
! Apply Boundary conditions & Remove from the linear system.
    CALL updateBoundary
    
! Viscoelastic part, calculates time(N+1) contribution:
    IF (param_beta.ne.1d0) THEN
      CALL cpu_time(cputime_section1)
      CALL calc_Upwinding_local_edge_nodes(V_x,V_y)
! DEVSS part:
      IF (param_beta_s.gt.0d0) THEN
        CALL applyDEVSS
      ENDIF

! Calculate Tau and add integration of Div(Tau) to RHS:
      IF (param_ve_model_choice.eq.1) THEN ! OldB / Giesekus
        CALL applyElasticStress
      ELSEIF (param_ve_model_choice.eq.2) THEN ! FENE-P-MP
        CALL applyElasticStress_FENE_PMP
      ELSE
        print*, 'Error: Unknown viscoelastic model selected: ', param_ve_model_choice, ' ', param_problem_choice
        print*, 'Stopping'
        STOP
      ENDIF
      CALL cpu_time(cputime_section2)
      cputime_viscoelastic = cputime_section2-cputime_section1
    ENDIF



   CALL removeBC_from_arrays ! Remove all traces of BCs from RHS and matrices
! Move velocity back one step before calculating latest.
    V_xNm2 = V_xNm1
    V_yNm2 = V_yNm1
    V_xNm1 = V_x
    V_yNm1 = V_y

    
    CALL cpu_time(cputime_temp1)
    cputime_setup=cputime_temp1-cputime_temp2
    
!!!! NEW PARDISO SOLVER !!!
    CALL pardiso_solve_stokes
    
    CALL cpu_time(cputime_temp2)
    cputime_solve=cputime_temp2-cputime_temp1

! Extrapolate local versions of all calculated variables.
    localpressureNm1=localpressure
    CALL extrapolate_pressure_local(localpressure)
    localV_xNm1=localV_x
    localV_yNm1=localV_y
    CALL extrapolate_velocity_local(localV_x,localV_y)
    
    CALL calc_local_gradient_of_U
    
! Calculate drag and sphere velocity (from drag) if required.    
    IF (circelm.ne.0) THEN
      CALL calc_drag
      IF (movingmeshflag.eq.1) THEN
        CALL calcSphereVelocity
      ENDIF
    ENDIF

! Calculate stopping check
!     IF (movingmeshflag.eq.0) THEN
      CALL calc_stopping_criteria(stopping_criteria,drag_criteria)
!     ENDIF

! Calculate any error analysis required
    IF (param_error) THEN
      CALL run_error_analysis
    ENDIF

!!! THINK ABOUT ORDER HERE.. should result analysis/output/calc of stopcheck come before updating SEM/Geom??? !!!


! OUTPUT SECTION:
    CALL cpu_time(cputime_total)
    CALL output_to_screen ! Always wanted.
! Added output to tecplot (non-fine) for all runs
! TODO: Add an input parameter (via CLI or input file) which turns this output on/off.
    printoutcount=printoutcount+1
    printed_to_tecplot = .FALSE.
    IF (enable_interim_output_to_file.AND.printoutcount.ge.print_threshold) THEN 
      CALL output_to_tecplot
      printed_to_tecplot = .TRUE.
      printoutcount=0
    ENDIF

! Check if solution has converged:
    IF (movingmeshflag.eq.0) THEN 
      IF ((stopping_criteria.lt.stopping_tol.and.drag_criteria.lt.stopping_tol) &
        .or.timeN.gt.time_limit) THEN
      IF (enable_interim_output_to_file.AND..NOT.printed_to_tecplot) THEN
        CALL output_to_tecplot
        ENDIF
        EXIT
      ENDIF
    ELSEIF (movingmeshflag.eq.1) THEN
      IF ((stopping_criteria.lt.stopping_tol.and.drag_criteria.lt.stopping_tol) &
        .or.timeN.gt.time_limit) THEN
      IF (enable_interim_output_to_file.AND..NOT.printed_to_tecplot) THEN
        CALL output_to_tecplot
      ENDIF
      EXIT
      ENDIF
    ENDIF

! Update time & geometry depedent parts for SEM and geometry.
    IF (movingmeshflag.eq.1) THEN ! Need to separate f_x, f_y and g and use them as references when updating RHS vectors.
      CALL applyMovingSphere
      CALL update_geometry
    ENDIF
    CALL updateSEM   


!!! NEEDS REVISING !!!
! Check for NaN problems   
    DO i=1,nptot
      IF (IsNaN(V_x(i))) THEN
        print*, 'ERROR! V_x(',i,') is NaN'
        STOP
      ELSEIF (IsNaN(V_y(i))) THEN
        print*, 'ERROR! V_y(',i,') is NaN'
        STOP
      ENDIF
    ENDDO
  ENDDO
  
  
  
! If required, output final solution in both fine and node form:
  IF(enable_final_output_to_file) THEN
! If the non-fine tecplot file isn't open yet, then open it.
    IF(.NOT.enable_interim_output_to_file) THEN
      OPEN(tecplot_output_fileid,FILE=tecplot_output_filename,IOSTAT=ierror)
    ENDIF
    CALL initialise_fine_grid ! setup fine node points
    CALL create_fine_solution ! generate solution(s) on these points

    CALL final_stage_output ! output to tecplot and along central axis for matlab.

    CLOSE(tecplot_output_fileid)
  ENDIF


! Memory release
  CALL release_pardiso_memory  
  CALL deassign_devss_memory
  
  CALL deassignMem

END PROGRAM main
