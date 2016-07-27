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
  IMPLICIT NONE
  CONTAINS



!     IF (param_problem_choice.eq.1) THEN
!       out_x = model_soln_vel_x(i)
!       out_y = model_soln_vel_y(i)
!     ELSEIF (param_problem_choice.eq.2) THEN
!       out_x = cylinder_soln_vel_x(i)
!       out_y = cylinder_soln_vel_y(i)
!     ELSEIF (param_problem_choice.eq.3) THEN
!       out_x = static_sphere_bcs_x(i)
!       out_y
!     ENDIF

! BEGIN FUNCTIONS FOR STOKES MODEL SOLUTION
  DOUBLE PRECISION FUNCTION model_soln_vel_x(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    
    out = sin(PI*nodeCoord(glob_i_in,1))*cos(PI*nodeCoord(glob_i_in,2))
    
  END FUNCTION model_soln_vel_x
  
  DOUBLE PRECISION FUNCTION model_soln_vel_y(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    
    out = -cos(PI*nodeCoord(glob_i_in,1))*sin(PI*nodeCoord(glob_i_in,2))
    
  END FUNCTION model_soln_vel_y
  
  DOUBLE PRECISION FUNCTION model_soln_f_x(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    DOUBLE PRECISION :: x,y
    
    x=PI*nodeCoord(glob_i_in,1)
    y=PI*nodeCoord(glob_i_in,2)
    
    out = 2d0*sin(x)*cos(y)*PI*PI + PI*cos(x)*sin(y)
    
  END FUNCTION model_soln_f_x
  
  DOUBLE PRECISION FUNCTION model_soln_f_y(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    DOUBLE PRECISION :: x,y
    
    x=PI*nodeCoord(glob_i_in,1)
    y=PI*nodeCoord(glob_i_in,2)
    
    out = -2d0*cos(x)*sin(y)*PI*PI + PI*sin(x)*cos(y)
    
  END FUNCTION model_soln_f_y
  
  DOUBLE PRECISION FUNCTION model_soln_pressure(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    
    x=PI*nodeCoord(glob_i_in,1)
    y=PI*nodeCoord(glob_i_in,2)
    
    out = sin(x)*sin(y)
    
  END FUNCTION model_soln_pressure
  
!	GRADIENTS

  DOUBLE PRECISION FUNCTION model_soln_Grad_vel_xx(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    
    out = PI*cos(PI*nodeCoord(glob_i_in,1))*cos(PI*nodeCoord(glob_i_in,2))
    
  END FUNCTION model_soln_Grad_vel_xx
  
  DOUBLE PRECISION FUNCTION model_soln_Grad_vel_yx(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    
    out = -PI*sin(PI*nodeCoord(glob_i_in,1))*sin(PI*nodeCoord(glob_i_in,2))
    
  END FUNCTION model_soln_Grad_vel_yx
 
  DOUBLE PRECISION FUNCTION model_soln_Grad_vel_xy(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    
    out = PI*sin(PI*nodeCoord(glob_i_in,1))*sin(PI*nodeCoord(glob_i_in,2))
    
  END FUNCTION model_soln_Grad_vel_xy
  
  DOUBLE PRECISION FUNCTION model_soln_Grad_vel_yy(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    
    out = -PI*cos(PI*nodeCoord(glob_i_in,1))*cos(PI*nodeCoord(glob_i_in,2))
    
  END FUNCTION model_soln_Grad_vel_yy
  
! END FUNCTIONS FOR STOKES MODEL SOLUTION

! BEGIN FUNCTIONS FOR CARTESIAN FLOW PAST A CYLINDER
  DOUBLE PRECISION FUNCTION cylinder_soln_vel_x(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    DOUBLE PRECISION :: x,y,r,z
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
    
! Flow past cylinder
! convert to spherical co-ordinates,
! NOTE: r & z cannot be 0 for the cylinder case - (0,0) is not part of domain.
    r = SQRT(x**2+y**2)
    th = acos(x/r)
    out = (rad_sphere**2/r**2 - 1d0)*(cos(th)**2 - 5d-1) + log(r/rad_sphere)
    
  END FUNCTION cylinder_soln_vel_x
  
  DOUBLE PRECISION FUNCTION cylinder_soln_vel_y(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    DOUBLE PRECISION :: x,y,r,z
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
    
! Flow past cylinder
! convert to spherical co-ordinates,
! NOTE: r & z cannot be 0 for the cylinder case - (0,0) is not part of domain.
    r = SQRT(x**2+y**2)
    th=acos(x/r)
    a = (rad_sphere**2/r**2 - 1d0)*cos(th)*sin(th)
    
  END FUNCTION cylinder_soln_vel_y
  
  DOUBLE PRECISION FUNCTION cylinder_soln_f_x(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in

    out=0d0
    
  END FUNCTION cylinder_soln_f_x
  
    DOUBLE PRECISION FUNCTION cylinder_soln_f_y(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    
    out = 0d0
    
  END FUNCTION cylinder_soln_f_y
  
  DOUBLE PRECISION FUNCTION cylinder_soln_pressure(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    DOUBLE PRECISION :: x,y,r,z
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
    
! Flow past cylinder
! convert to spherical co-ordinates,
! NOTE: r & z cannot be 0 for the cylinder case - (0,0) is not part of domain.
    r = dsqrt(x**2+y**2)
    th = dacos(x/r)
    out = -2d0*dcos(th)/r
    
  END FUNCTION cylinder_soln_pressure
  
  !	GRADIENTS

  DOUBLE PRECISION FUNCTION cylinder_soln_Grad_vel_xx(glob_i_in) RESULT(out)
! Diff of vel_x wrt x
    IMPLICIT NONE
    INTEGER :: glob_i_in
    DOUBLE PRECISION :: x,y,r,z
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
    
! Flow past cylinder
! convert to spherical co-ordinates,
! NOTE: r & z cannot be 0 for the cylinder case - (0,0) is not part of domain.
    r = dsqrt(x**2+y**2)
    th = dacos(x/r)
    
    a = (x/r)*((-2d0*rad_sphere**2/r**3)*(dcos(th)**2 - 5d-1) + 1d0/r) + &
	  (2d0*y/r**2)*((rad_sphere/r)**2 - 1d0)*dcos(th)*dsin(th)
    
  END FUNCTION cylinder_soln_Grad_vel_xx
  
  DOUBLE PRECISION FUNCTION cylinder_soln_Grad_vel_yx(glob_i_in) RESULT(out)
! Diff of vel_x wrt y
    IMPLICIT NONE
    INTEGER :: glob_i_in
    DOUBLE PRECISION :: x,y,r,z
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)

! Flow past cylinder
! convert to spherical co-ordinates,
! NOTE: r & z cannot be 0 for the cylinder case - (0,0) is not part of domain.
    r = dsqrt(x**2+y**2)
    th = dacos(x/r)

    out =  (y/r)*((-2d0*rad_sphere**2/r**3)*(dcos(th)**2 - 5d-1) + 1d0/r) - &
	  (2d0*x/r**2)*((rad_sphere/r)**2 - 1d0)*dcos(th)*dsin(th)
    
    
  END FUNCTION cylinder_soln_Grad_vel_yx
 
  DOUBLE PRECISION FUNCTION cylinder_soln_Grad_vel_xy(glob_i_in) RESULT(out)
! Diff of vel_y wrt x
    IMPLICIT NONE
    INTEGER :: glob_i_in
    DOUBLE PRECISION :: x,y,r,z
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
    
! Flow past cylinder
! convert to spherical co-ordinates,
! NOTE: r & z cannot be 0 for the cylinder case - (0,0) is not part of domain.
    r = dsqrt(x**2+y**2)
    th=dacos(x/r)

    out = (x/r)*(-2d0*rad_sphere**2/r**3)*dcos(th)*dsin(th) -&
	  (y/r**2)*((rad_sphere/r)**2 - 1d0)*(dcos(th)**2 - dsin(th)**2)
    
  END FUNCTION cylinder_soln_Grad_vel_xy
  
  DOUBLE PRECISION FUNCTION cylinder_soln_Grad_vel_yy(glob_i_in) RESULT(out)
! Diff of vel_y wrt y
    IMPLICIT NONE
    INTEGER :: glob_i_in
    DOUBLE PRECISION :: x,y,r,z
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
    
! Flow past cylinder
! convert to spherical co-ordinates,
! NOTE: r & z cannot be 0 for the cylinder case - (0,0) is not part of domain.
    r = dsqrt(x**2+y**2)
    th = dacos(x/r)

    out = (y/r)*(-2d0*rad_sphere**2/r**3)*dcos(th)*dsin(th) + &
	  (x/r**2)*((rad_sphere/r)**2 -1d0)*(dcos(th)**2 - dsin(th)**2)
	  
  END FUNCTION cylinder_soln_Grad_vel_yy
  


! BEGIN FUNCTIONS FOR STOKES BENCHMARK: UNIFORM FLOW PAST A FIXED SPHERE
  DOUBLE PRECISION FUNCTION uniform_fixed_sphere_vel_x(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    
    
    out = 1d0
    
  END FUNCTION uniform_fixed_sphere_vel_x
  
  DOUBLE PRECISION FUNCTION uniform_fixed_sphere_vel_y(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    
    out = 0d0
    
  END FUNCTION uniform_fixed_sphere_vel_y
  
  DOUBLE PRECISION FUNCTION uniform_fixed_sphere_f_x(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    
    out = 0d0
    
  END FUNCTION uniform_fixed_sphere_f_x
  
  DOUBLE PRECISION FUNCTION uniform_fixed_sphere_f_y(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in

    
    out = 0d0
    
  END FUNCTION uniform_fixed_sphere_f_y
  
  DOUBLE PRECISION FUNCTION uniform_fixed_sphere_pressure(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    
    out = -999d0
    
  END FUNCTION uniform_fixed_sphere_pressure
  
!	GRADIENTS

  DOUBLE PRECISION FUNCTION uniform_fixed_sphere_Grad_vel_xx(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    
    out = -999d0
    
  END FUNCTION uniform_fixed_sphere_Grad_vel_xx
  
  DOUBLE PRECISION FUNCTION uniform_fixed_sphere_Grad_vel_yx(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    
    out = -999d0
    
  END FUNCTION uniform_fixed_sphere_Grad_vel_yx
 
  DOUBLE PRECISION FUNCTION uniform_fixed_sphere_Grad_vel_xy(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    
    out = -999d0
    
  END FUNCTION uniform_fixed_sphere_Grad_vel_xy
  
  DOUBLE PRECISION FUNCTION uniform_fixed_sphere_Grad_vel_yy(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER :: glob_i_in
    
    out = -999d0
    
  END FUNCTION uniform_fixed_sphere_Grad_vel_yy
  
! END FUNCTIONS FOR STOKES MODEL SOLUTION

  
  
  

END MODULE boundary_module