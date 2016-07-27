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

MODULE boundary_module
  USE shared_data
  USE constants
  USE functions_module
  USE waters_solution_module
  IMPLICIT NONE
  CONTAINS
  
  SUBROUTINE updateBoundary
    IMPLICIT NONE
    INTEGER :: i
    IF (param_beta.lt.1d0) THEN
      CALL generate_boundary_stress
    ENDIF
    IF (movingmeshflag.eq.0) THEN
      V_sphere=1d0
      CALL calcBoundaryConditionsFIXED
    ELSE
      CALL calcBoundaryConditions
    ENDIF
!    CALL applyBCs_to_solution ! NOT NEEDED UNTIL AFTER SOLVER ?????? - check and update as appropriate!
    CALL calcBoundaryRHS
    CALL calcBoundaryRHS2
    CALL applyBC_to_RHS
    
  END SUBROUTINE updateBoundary
  
  SUBROUTINE applyBCs_to_solution
! Simply subs-in the boundary values on dirichlet nodes onto the global solution
    IMPLICIT NONE
    INTEGER :: i
    
    DO i=1,npedg

      IF (bdflag(1,i)) THEN
	V_x(i) = boundary_x(i)
      ENDIF
      IF (bdflag(2,i)) THEN
	V_y(i) = boundary_y(i)
      ENDIF
    ENDDO

  END SUBROUTINE applyBCs_to_solution
  

  SUBROUTINE calcBoundaryConditions
! Calculates the vector of known boundary terms for a moving mesh.
    IMPLICIT NONE
    INTEGER :: i
    boundary_x=0d0
    boundary_y=0d0
    DO i=1,npedg
      IF (wallflag(i)) THEN
	boundary_x(i) = 0d0
	boundary_y(i) = 0d0
      ELSEIF (circnodeflag(i)) THEN
	boundary_x(i) = V_sphere
	boundary_y(i) = 0d0
      ELSEIF (inflowflag(i)) THEN
	boundary_x(i) = 0d0
	boundary_y(i) = 0d0
      ELSEIF (outflowflag(i)) THEN
	boundary_x(i) = 0d0
	boundary_y(i) = 0d0
      ELSEIF (wallsymmflag(i)) THEN
	boundary_x(i) = 0d0
	boundary_y(i) = 0d0
      ELSE !In case we use a non-boundary node for anything elsewhere in the code, this will make it easier to track.
	boundary_x(i) = -1d99
	boundary_y(i) = -1d99
      ENDIF
    ENDDO

  END SUBROUTINE calcBoundaryConditions

  DOUBLE PRECISION FUNCTION wallterm_x(x,y)
! x component of the wall boundary condition 
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x,y
    
    IF (movingmeshflag.eq.0) THEN
      wallterm_x = 1d0
!       wallterm_x = 1d0-y**2
!       wallterm_x = -4d0*y*(y-1d0)
    ELSE
      wallterm_x = 0d0
    ENDIF

  END FUNCTION wallterm_x
  
  DOUBLE PRECISION FUNCTION wallterm_y(x,y)
! y component of the wall boundary condition
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x,y
		
    wallterm_y = 0d0
		
  END FUNCTION wallterm_y
  
  DOUBLE PRECISION FUNCTION circterm_x(x,y)
! x component of the circular boundary condition
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x,y
    
    circterm_x = V_sphere!0d0
    
  END FUNCTION circterm_x
  
  DOUBLE PRECISION FUNCTION circterm_y(x,y)
! y component of the circular boundary condition
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x,y
		
    circterm_y = 0d0
		
  END FUNCTION circterm_y
  
  DOUBLE PRECISION FUNCTION inflowterm_x(x,y)
! x component of the inflow condition
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x,y
! Poiseuille flow:		
!    inflowterm_x = -4d0*y*(y-1d0)!1d0-y**2
!     inflowterm_x = 1d0-y**2
! Coette flow:
!    inflowterm_x = y
! BCs for bounded cylinder
    IF (movingmeshflag.eq.0) THEN
      inflowterm_x = 1d0
!       inflowterm_x = 1d0-y**2
!       inflowterm_x = -4d0*y*(y-1d0)
    ELSE
      inflowterm_x = 0d0
    ENDIF
  
  END FUNCTION inflowterm_x
  
  DOUBLE PRECISION FUNCTION inflowterm_y(x,y)
! y component of the inflow condition
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x,y
    
    inflowterm_y = 0d0

  END FUNCTION inflowterm_y
  
  DOUBLE PRECISION FUNCTION outflowterm_x(x,y)
! x component of the outflow condition
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x,y
!Poiseuille flow:		
!     outflowterm_x = -4d0*y*(y-1d0)!1d0-y**2
!     outflowterm_x = 1d0-y**2
! Coette flow:
!    outflowterm_x = y
! BCs for bounded cylinder
    IF (movingmeshflag.eq.0) THEN
      outflowterm_x = 1d0
!      outflowterm_x = -4d0*y*(y-1d0)
!       outflowterm_x = 1d0-y**2
    ELSE
      outflowterm_x = 0d0
    ENDIF
    
  END FUNCTION outflowterm_x
  

  
  DOUBLE PRECISION FUNCTION outflowterm_y(x,y)
! y component of the outflow condition
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x,y
    
    outflowterm_y = 0d0

  END FUNCTION outflowterm_y
  
  
  DOUBLE PRECISION FUNCTION wallsymmterm_x(x,y)
! x component of the symmetric wall condition
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x,y
		
    wallsymmterm_x = 0d0!-1d6
		
  END FUNCTION wallsymmterm_x
  
  DOUBLE PRECISION FUNCTION wallsymmterm_y(x,y)
! y component of the symmetric wall condition
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x,y
		
    wallsymmterm_y = 0d0
		
  END FUNCTION wallsymmterm_y
  
  SUBROUTINE calcBoundaryConditionsFIXED
! Calculates the vector of known boundary terms for a fixed mesh.
    IMPLICIT NONE
    INTEGER :: i
    boundary_x=0d0
    boundary_y=0d0
    
   
! GENERAL WAY:   
    DO i=1,npedg
      IF (wallflag(i)) THEN
	boundary_x(i) = FIXEDwallterm_x(i)
	boundary_y(i) = FIXEDwallterm_y(i)
      ELSEIF (circnodeflag(i)) THEN
	boundary_x(i) = FIXEDcircterm_x(i)
	boundary_y(i) = FIXEDcircterm_y(i)
      ELSEIF (inflowflag(i)) THEN
	boundary_x(i) = FIXEDinflowterm_x(i)
	boundary_y(i) = FIXEDinflowterm_y(i)
      ELSEIF (outflowflag(i)) THEN
	boundary_x(i) = FIXEDoutflowterm_x(i)
	boundary_y(i) = FIXEDoutflowterm_y(i)
      ELSEIF (wallsymmflag(i)) THEN
	boundary_x(i) = FIXEDwallsymmterm_x(i)
	boundary_y(i) = FIXEDwallsymmterm_y(i)
      ELSE !In case we use a non-boundary node for anything elsewhere in the code, this will make it easier to track.
	boundary_x(i) = -1d99
	boundary_y(i) = -1d99
      ENDIF
    ENDDO
  END SUBROUTINE calcBoundaryConditionsFIXED
  
  DOUBLE PRECISION FUNCTION FIXEDwallterm_x(i) RESULT(a)
! x component of the wall boundary condition 
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i

    IF (param_function_choice.eq.2) THEN
      a = V_sphere
    ELSEIF (param_function_choice.eq.3) THEN
      a = model_soln_vel_x(i)
    ELSEIF (param_function_choice.eq.4) THEN
      a = cylinder_soln_vel_x(i)
    ELSE
      a = 0d0 ! Includes 1 5 6 7 8 9.
! ADD MORE AS APPROPRIATE.
    ENDIF

  END FUNCTION FIXEDwallterm_x
  
  DOUBLE PRECISION FUNCTION FIXEDwallterm_y(i) RESULT(a)
! y component of the wall boundary condition
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    
    IF (param_function_choice.eq.3) THEN
      a = model_soln_vel_y(i)
    ELSEIF (param_function_choice.eq.4) THEN
      a = cylinder_soln_vel_y(i)
    ELSE
      a = 0d0 ! Includes 1 2 5 6 7 8.
! ADD MORE AS APPROPRIATE.
    ENDIF
    


  END FUNCTION FIXEDwallterm_y
  
  DOUBLE PRECISION FUNCTION FIXEDcircterm_x(i) RESULT(a)
! x component of the circular boundary condition
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i

    IF (param_function_choice.eq.3) THEN
      a = model_soln_vel_x(i)
    ELSEIF (param_function_choice.eq.4) THEN
      a = cylinder_soln_vel_x(i) ! SHOULD BE ZERO?!?!?!
    ELSE
      a = 0d0 ! Includes 1 2 5 6 7 8.

! ADD MORE AS APPROPRIATE.
    ENDIF
    
  END FUNCTION FIXEDcircterm_x
  
  DOUBLE PRECISION FUNCTION FIXEDcircterm_y(i) RESULT(a)
! y component of the circular boundary condition
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    
    IF (param_function_choice.eq.3) THEN
      a = model_soln_vel_y(i)
    ELSEIF (param_function_choice.eq.4) THEN
      a = cylinder_soln_vel_y(i)
    ELSE
      a = 0d0 ! Includes 1 2 5 6 7 8.
! ADD MORE AS APPROPRIATE.
    ENDIF
		
  END FUNCTION FIXEDcircterm_y
  
  DOUBLE PRECISION FUNCTION FIXEDinflowterm_x(i) RESULT(a)
! x component of the inflow condition
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    
    IF (param_function_choice.eq.1) THEN
      a = 3d0/8d0*(4d0-nodeCoord(i,2)**2)
    ELSEIF (param_function_choice.eq.2) THEN
      a = V_sphere
    ELSEIF (param_function_choice.eq.3) THEN
      a = model_soln_vel_x(i)
    ELSEIF (param_function_choice.eq.4) THEN
      a = cylinder_soln_vel_x(i)
    ELSEIF (param_waters) THEN ! Includes 5 6 7 8.
      a = transient_u(i)
    ELSE
      a = 0d0 ! catch all.
! ADD MORE AS APPROPRIATE.
    ENDIF

    
  
  END FUNCTION FIXEDinflowterm_x
  
  DOUBLE PRECISION FUNCTION FIXEDinflowterm_y(i) RESULT(a)
! y component of the inflow condition
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    

    IF (param_function_choice.eq.3) THEN
      a = model_soln_vel_y(i)
    ELSEIF (param_function_choice.eq.4) THEN
      a = cylinder_soln_vel_y(i)
    ELSE
      a = 0d0 ! Includes 1 2 5 6 7 8.   
! ADD MORE AS APPROPRIATE.
    ENDIF

  END FUNCTION FIXEDinflowterm_y
  
  DOUBLE PRECISION FUNCTION FIXEDoutflowterm_x(i) RESULT(a)
! x component of the outflow condition, FIXED MESH CASE
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    
    IF (param_function_choice.eq.1) THEN
      a = 3d0/8d0*(4d0-nodeCoord(i,2)**2)
    ELSEIF (param_function_choice.eq.2) THEN
      a = V_sphere
    ELSEIF (param_function_choice.eq.3) THEN
      a = model_soln_vel_x(i)
    ELSEIF (param_function_choice.eq.4) THEN
      a = cylinder_soln_vel_x(i)
    ELSEIF (param_waters) THEN ! Includes 5 6 7 8.
      a = transient_u(i)
    ELSE
      a = 0d0 ! catch all.
! ADD MORE AS APPROPRIATE.
    ENDIF
    
  END FUNCTION FIXEDoutflowterm_x
  
  DOUBLE PRECISION FUNCTION FIXEDoutflowterm_y(i) RESULT(a)
! y component of the outflow condition
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    
    IF (param_function_choice.eq.3) THEN
      a = model_soln_vel_y(i)
    ELSEIF (param_function_choice.eq.4) THEN
      a = cylinder_soln_vel_y(i)
    ELSE
      a = 0d0 ! Includes 1 2 5 6 7 8.
! ADD MORE AS APPROPRIATE.
    ENDIF
    
END FUNCTION FIXEDoutflowterm_y
  
  
  DOUBLE PRECISION FUNCTION FIXEDwallsymmterm_x(i) RESULT(a)
! x component of the symmetric wall condition
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    
    IF (param_function_choice.eq.3) THEN
      a = model_soln_vel_x(i)
    ELSEIF (param_function_choice.eq.4) THEN
      a = cylinder_soln_vel_x(i)
    ELSE
      a = 0d0 ! Includes 1 2 5 6 7 8.
! ADD MORE AS APPROPRIATE.
    ENDIF
		
  END FUNCTION FIXEDwallsymmterm_x
  
  DOUBLE PRECISION FUNCTION FIXEDwallsymmterm_y(i) RESULT(a)
! y component of the symmetric wall condition
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i


    IF (param_function_choice.eq.3) THEN
      a = model_soln_vel_y(i)
    ELSEIF (param_function_choice.eq.4) THEN
      a = cylinder_soln_vel_y(i)
    ELSE
      a = 0d0 ! Includes 1 2 5 6 7 8.
! ADD MORE AS APPROPRIATE.
    ENDIF
		
  END FUNCTION FIXEDwallsymmterm_y
  
  
  SUBROUTINE calcBoundaryRHS
! Calculates
  ! gl are the G-L points
  ! N - as usual
  ! bonundary = output
  ! uses the function bdterms below, which is the known function on the boundary
  ! will need to update this to include non-dirichlet conditions!
  
    IMPLICIT NONE
    INTEGER :: el,kl,ij,i
    DOUBLE PRECISION :: temp
    DOUBLE PRECISION, DIMENSION(0:NP1SQM1) :: bdloc
    DOUBLE PRECISION, DIMENSION(nptot) :: tempglob

    boundaryContribution_x=0d0
    DO el=1,numelm
      tempglob=0d0
      bdloc=0d0
      DO kl=0,NP1SQM1
! 	IF (bdflag(1,mapg(kl,el))) CYCLE ! Not needed as we zero the boundary node rows/columns anyway !
	temp = 0d0
	DO ij=0,NP1SQM1
	  IF (bdflag(1,mapg(ij,el))) THEN
	    temp = temp + boundary_x(mapg(ij,el))*A_x(kl,ij,el)
	  ENDIF
	ENDDO
	bdloc(kl) = temp
      ENDDO
      CALL vecglobalprolongation(bdloc,el,tempglob)
      boundaryContribution_x = boundaryContribution_x + tempglob
    ENDDO
    
    boundaryContribution_y=0d0
    DO el=1,numelm
      tempglob=0d0
      bdloc=0d0
      DO kl=0,NP1SQM1
! 	IF (bdflag(2,mapg(kl,el))) CYCLE ! Not needed as we zero the boundary node rows/columns anyway !
	temp = 0d0
	DO ij=0,NP1SQM1
	  IF (bdflag(2,mapg(ij,el))) THEN
	    temp = temp + boundary_y(mapg(ij,el))*A_y(kl,ij,el)
	  ENDIF
	ENDDO
	bdloc(kl) = temp
      ENDDO
      CALL vecglobalprolongation(bdloc,el,tempglob)
      boundaryContribution_y = boundaryContribution_y + tempglob
    ENDDO
    
  END SUBROUTINE calcBoundaryRHS
  
  SUBROUTINE calcBoundaryRHS2
! Takes the known terms which multiply with the B matrix & puts them onto the right hand side, into a vector of size npint, called g.
    IMPLICIT NONE
    
    INTEGER :: i,j,k,l,ij,kl,el
    DOUBLE PRECISION :: temp1,temp2
    DOUBLE PRECISION, DIMENSION(NM1SQ) :: bdloc
    DOUBLE PRECISION, DIMENSION(npint) :: tempglob
! Contribution from the Divergence part of the stiffness matrix.
! possible g should be replaced by something else and then updated similar to f_x and f_y ??


    boundaryContribution_p=0d0
    DO el=1,numelm
      tempglob=0d0
      bdloc=0d0
      DO kl=1,NM1SQ
	temp1=0d0
	temp2=0d0
	DO ij=0,NP1SQM1
	  IF (bdflag(1,mapg(ij,el))) THEN
	    temp1 = temp1 + boundary_x(mapg(ij,el))*B_x(kl,ij,el)
	  ENDIF
	  IF (bdflag(2,mapg(ij,el))) THEN
	    temp2 = temp2 + boundary_y(mapg(ij,el))*B_y(kl,ij,el)
	  ENDIF
	ENDDO
	bdloc(kl)=temp1+temp2
      ENDDO
      CALL vecglobalprolongation_internal_nodes(bdloc,el,tempglob)
      boundaryContribution_p = boundaryContribution_p + tempglob
    ENDDO

  END SUBROUTINE calcBoundaryRHS2
    
  
  SUBROUTINE applyBC_to_RHS
    IMPLICIT NONE
    
    f_x = f_x - boundaryContribution_x
    f_y = f_y - boundaryContribution_y
    g = g - boundaryContribution_p
    
  END SUBROUTINE applyBC_to_RHS
  
  SUBROUTINE removeBC_from_arrays
    IMPLICIT NONE
    INTEGER :: el,ij,kl
!
! After transferring known terms to RHS, we put zeros in the boundary entries of our matrices
! 
    DO el=1,numelm
      DO ij=0,NP1SQM1
	IF(bdflag(1,mapg(ij,el))) THEN
	  f_x(mapg(ij,el))=0d0
	  DO kl=0,NP1SQM1	  	  
	    A_x(ij,kl,el) = 0d0
	    A_x(kl,ij,el) = 0d0
	  ENDDO
	  A_x(ij,ij,el) = 1d0
	  DO kl=1,NM1SQ
	    B_x(kl,ij,el) = 0d0
	  ENDDO
	ENDIF
	IF(bdflag(2,mapg(ij,el))) THEN
	  f_y(mapg(ij,el))=0d0
	  DO kl=0,NP1SQM1  
	    A_y(ij,kl,el) = 0d0
	    A_y(kl,ij,el) = 0d0
	  ENDDO
	  A_y(ij,ij,el) = 1d0
	  DO kl=1,NM1SQ
	    B_y(kl,ij,el) = 0d0
	  ENDDO
	ENDIF
      ENDDO
    ENDDO
    
  END SUBROUTINE removeBC_from_arrays
  

 
!
! Pretty much all other functions are redundant from here on:
!

  DOUBLE PRECISION FUNCTION bdterm_x(x,y)
    IMPLICIT NONE
    DOUBLE PRECISION :: x,y
    
    

! Planar flow for unit sphere
!     IF ((x.eq.1d0.OR.x.eq.-1d0).AND.y.eq.0d0) THEN
!       bdterm_x = 0d0
!     ELSEIF (y.eq.0d0.OR.y.eq.2d0.OR.x.eq.-8d0.OR.x.eq.8d0) THEN
!       bdterm_x = 1d0
!     ELSE
!       bdterm_x = 0d0
!     ENDIF
      
! Planar flow for half unit sphere @(0,0)
!     IF ((x.eq.5d-1.OR.x.eq.-5d-1).AND.y.eq.0d0) THEN
!       bdterm_x = 0d0
!     ELSEIF (y.eq.0d0.OR.y.eq.1d0.OR.x.eq.-4d0.OR.x.eq.4d0) THEN
!       bdterm_x = 1d0
!     ELSE
!       bdterm_x = 0d0
!     ENDIF

! Planar flow for half unit sphere @(1,0)
    IF ((x.eq.5d-1.OR.x.eq.15d-1).AND.y.eq.0d0) THEN
      bdterm_x = 0d0
    ELSEIF (y.eq.0d0.OR.y.eq.1d0.OR.x.eq.-3d0.OR.x.eq.5d0) THEN
      bdterm_x = 1d0
    ELSE
      bdterm_x = 0d0
    ENDIF
!     IF ((x.eq.sphereXleft.OR.x.eq.sphereXright).AND.y.eq.centre_of_sphere(2)) THEN
!       bdterm_x = 0d0
!     ELSEIF (y.eq.centre_of_sphere(2).OR.y.eq.1d0.OR.x.eq.inflowX.OR.x.eq.outflowX) THEN
!       bdterm_x = 1d0
!     ELSE
!       bdterm_x = 0d0
!     ENDIF
    
!Poiseuille flow:
!     IF (x.eq.-4d0.or.x.eq.4d0) THEN
!       bdterm_x = -4d0*y*(y-1d0)
!       bdterm_x = 1d0-y**2
!     ELSE
!       bdterm_x = 0d0
!     ENDIF

! Poiseuille Axisymm case:
!     bdterm_x = 1d0-y**2
    
! Non-trivial soln:
!     bdterm_x = sin(PI*x)*cos(PI*y)

!     bdterm_x = exp(x+y)
!     bdterm_x = sin(PI*x)*sin(PI*y)
!     bdterm_x = cos(PI*x)*cos(PI*y)
!     bdterm_x = cos(0.5*PI*x)*cos(0.5*PI*y)
!     bdterm_x = sin(2*PI*x)*sin(2*PI*y)
!     bdterm_x = (1-x**2)*(1-y**2)*exp(x+y-2)
!     bdterm_x = x*(x-1)*(x+1-y)*(x+y) ! for geometry x1 = 0,0 x2= 1,-1 x3=1,2 x4=0,1
    RETURN
  END FUNCTION bdterm_x
  
  DOUBLE PRECISION FUNCTION bdterm_y(x,y)
    IMPLICIT NONE
    DOUBLE PRECISION :: x,y
!Poiseuille flow
    bdterm_y = 0d0
    
! Non-trivial soln:
!     bdterm_y = -cos(PI*x)*sin(PI*y)
!     bdterm_y=-exp(x+y)
!     bdterm_y = sin(PI*x)*sin(PI*y)
!     bdterm_y = cos(PI*x)*cos(PI*y)
!     bdterm_y = cos(0.5*PI*x)*cos(0.5*PI*y)
!     bdterm_y = sin(2*PI*x)*sin(2*PI*y)
!     bdterm_y = (1-x**2)*(1-y**2)*exp(x+y-2)
!     bdterm_y = x*(x-1)*(x+1-y)*(x+y) ! for geometry x1 = 0,0 x2= 1,-1 x3=1,2 x4=0,1
    RETURN
  END FUNCTION bdterm_y
  

  
! Function to define the actual solution, u(x,y) - UPDATE TO REFLECT THE PROBLEM
! USED FOR CODE VALIDATION!!!
  DOUBLE PRECISION FUNCTION uterm_x(x,y)
    IMPLICIT NONE
    DOUBLE PRECISION :: x,y

! Planar flow for unit sphere
!     IF ((x.eq.1d0.OR.x.eq.-1d0).AND.y.eq.0d0) THEN
!       uterm_x = 0d0
!     ELSEIF (y.eq.0d0.OR.y.eq.2d0.OR.x.eq.-8d0.OR.x.eq.8d0) THEN
!       uterm_x = 1d0
!     ELSE
!       uterm_x = 0d0
!     ENDIF

! Planar flow for half unit sphere @ (0,0)
!     IF ((x.eq.5d-1.OR.x.eq.-5d-1).AND.y.eq.0d0) THEN
!       uterm_x = 0d0
!     ELSEIF (y.eq.0d0.OR.y.eq.1d0.OR.x.eq.-4d0.OR.x.eq.4d0) THEN
!       uterm_x = 1d0
!     ELSE
!       uterm_x = 0d0
!     ENDIF
    
! Planar flow for half unit sphere @(1,0)
    IF ((x.eq.5d-1.OR.x.eq.15d-1).AND.y.eq.0d0) THEN
      uterm_x = 0d0
    ELSEIF (y.eq.0d0.OR.y.eq.1d0.OR.x.eq.-3d0.OR.x.eq.5d0) THEN
      uterm_x = 1d0
    ELSE
      uterm_x = 0d0
    ENDIF
      
    
!Poiseuille flow
!     IF (x.eq.-4d0.OR.x.eq.4d0) THEN
!       uterm_x=-4d0*y*(y-1d0)
! ! 	uterm_x=1d0-y**2
!     ELSE
!       uterm_x=0d0
!     ENDIF
!      uterm_x=-4d0*y*(y-1d0)


! Non-trivial soln:
!     uterm_x = sin(PI*x)*cos(PI*y)

!     uterm_x = exp(x+y)
!     uterm_x = sin(PI*x)*sin(PI*y)
!     uterm_x = cos(PI*x)*cos(PI*y)
!     uterm_x = cos(0.5*PI*x)*cos(0.5*PI*y)
!     uterm_x = sin(2*PI*x)*sin(2*PI*y)
!     uterm_x = (1-x**2)*(1-y**2)*exp(x+y-2)
!     uterm_x = x*(x-1)*(x+1-y)*(x+y) ! for geometry x1 = 0,0 x2= 1,-1 x3=1,2 x4=0,1
    RETURN
  END FUNCTION uterm_x
  
  DOUBLE PRECISION FUNCTION uterm_y(x,y)
    IMPLICIT NONE
    DOUBLE PRECISION :: x,y
!Poiseuille flow
    uterm_y = 0d0
    
! Non-trivial soln:
!     uterm_y = -cos(PI*x)*sin(PI*y)


!     uterm_y = -exp(x+y)
!     uterm_y = sin(PI*x)*sin(PI*y)
!     uterm_y = cos(PI*x)*cos(PI*y)
!     uterm_y = cos(0.5*PI*x)*cos(0.5*PI*y)
!     uterm_y = sin(2*PI*x)*sin(2*PI*y)
!     uterm_y = (1-x**2)*(1-y**2)*exp(x+y-2)
!     uterm_y = x*(x-1)*(x+1-y)*(x+y) ! for geometry x1 = 0,0 x2= 1,-1 x3=1,2 x4=0,1
    RETURN
  END FUNCTION uterm_y
  
  DOUBLE PRECISION FUNCTION pressureterm(x,y)
    IMPLICIT NONE
    DOUBLE PRECISION :: x,y
!Poiseuille flow
!     pressureterm = -2d0*x 
    pressureterm = 0d0!-8d0*x
    
!Non-Trivial soln
!     pressureterm = sin(PI*x)*sin(PI*y)

!For non-trig pressure term (in order to avoid problems?)
!     pressureterm = x**2d0!x!x**3 + y**3
    RETURN
  END FUNCTION pressureterm
  
  SUBROUTINE generate_boundary_stress
    IMPLICIT NONE
    INTEGER :: i
	  
    IF (param_function_choice.eq.1) THEN ! Steady Poiseuille flow (2-D)
      DO i=1,npedg
	boundary_stress_xx(i) = 2d0*We*(1d0-param_beta)*((-3d0/4d0*nodeCoord(i,2))**2)
	boundary_stress_xy(i) = (1d0-param_beta)*(-3d0/4d0*nodeCoord(i,2))
	boundary_stress_yy(i) = 0d0
      ENDDO
    ELSEIF (param_function_choice.eq.2) THEN ! Steady Uniform flow (3-D)
      DO i=1,npedg
	boundary_stress_xx(i) = 0d0
	boundary_stress_xy(i) = 0d0
	boundary_stress_yy(i) = 0d0
	boundary_stress_zz(i) = 0d0
      ENDDO  
    ELSEIF (param_function_choice.eq.7) THEN ! 2-D Waters solution
      DO i=1,npedg
	boundary_stress_xx(i)=transient_txx(i)
	boundary_stress_xy(i)=transient_txy(i)
	boundary_stress_yy(i)=0d0
      ENDDO
    ELSEIF (param_function_choice.eq.8) THEN ! 3-D Waters solution.. extrapolated from known velocity and gradV
      DO i=1,npedg
	boundary_stress_xx(i)=transient_txx(i)
	boundary_stress_xy(i)=transient_txy(i)
	boundary_stress_yy(i)=0d0
	boundary_stress_zz(i)=0d0
      ENDDO
    ELSEIF (param_function_choice.eq.9) THEN ! Moving mesh - zero inflow velocity.
      DO i=1,npedg
	boundary_stress_xx(i) = 0d0
	boundary_stress_xy(i) = 0d0
	boundary_stress_yy(i) = 0d0
	boundary_stress_zz(i) = 0d0
      ENDDO  
    ENDIF	
	
  END SUBROUTINE generate_boundary_stress
  
END MODULE boundary_module

