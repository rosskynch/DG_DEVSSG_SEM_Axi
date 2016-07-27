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

MODULE movingSphere_module
  USE constants
  USE shared_data
  
  IMPLICIT NONE
  CONTAINS


  SUBROUTINE applyMovingSphere
    IMPLICIT NONE
    INTEGER :: i
    DOUBLE PRECISION :: temp
! MUST ALREADY HAVE CALCULATED THE SPHERE VELOCITY
! 
!  First order:
!     temp=deltat*V_sphere
! Update the mesh co-ordinates using 2nd order AB scheme.
! Must make sure not to move nodes at each end of the pipe.
    temp = (0.5d0*deltat)*(3d0*V_sphere - V_sphereNm1)

    DO i=1,numnp
      IF (.not.fixed_node(i)) THEN  
	vertx(i,1)=vertx(i,1) + temp 
      ENDIF
    ENDDO

    centre_of_sphere(1)=centre_of_sphere(1) + temp
  
  END SUBROUTINE applyMovingSphere
  
  
  
  
  
  SUBROUTINE calcSphereVelocity
    IMPLICIT NONE
    DOUBLE PRECISION :: temp1,temp2,temp3,densratio
    
    V_sphereNm2=V_sphereNm1
    V_sphereNm1=V_sphere
    densratio=rho_f/rho_s !usually would be rho_s/rho_f
  
!     rho_s = 5d0 ! Sphere density
!     rho_f = 1d0 ! Fluid density
!     gravity_const = 1d0!9.81d0 ! Gravity constant (non-dimensional)

!     temp = (4d0/3d0)*pi*((rho_s/rho_f) - 1d0)*g - (1d0/Re)*drag
!     temp = temp*(3d0*rho_f)/(4d0*pi*rho_s)
! 3 step adams bashforth:
! OLD VERSION (Re/We):
!     temp1 = ((4d0/3d0)*PI*((rho_s/rho_f)-1d0)*gravity_const + (1d0/Re)*drag)*(3d0*rho_f)/(4d0*PI*rho_s)
!     temp2 = ((4d0/3d0)*PI*((rho_s/rho_f)-1d0)*gravity_const + (1d0/Re)*dragNm1)*(3d0*rho_f)/(4d0*PI*rho_s)
!     temp3 = ((4d0/3d0)*PI*((rho_s/rho_f)-1d0)*gravity_const + (1d0/Re)*dragNm2)*(3d0*rho_f)/(4d0*PI*rho_s)


! Change of non-dim to use radius of cylinder as the length.
! not sure if this works
! For now the radius of the cylinder is 1. If it was not we would also need to divide other terms by it, have put this in now
! but it may need checking. (rad_cyl/rad_sphere)**2

!     temp1 = ((4d0/3d0)*PI*(densratio-1d0)*gravity_const + (rad_cyl/rad_sphere**3)*(1d0/Re)*drag)*3d0/(4d0*PI*densratio)
!     temp2 = ((4d0/3d0)*PI*(densratio-1d0)*gravity_const + (rad_cyl/rad_sphere**3)*(1d0/Re)*dragNm1)*3d0/(4d0*PI*densratio)
!     temp3 = ((4d0/3d0)*PI*(densratio-1d0)*gravity_const + (rad_cyl/rad_sphere**3)*(1d0/Re)*dragNm2)*3d0/(4d0*PI*densratio)

!     temp1 = (1d0-1d0/densratio)*gravity_const + 3d0*drag/(densratio*4d0*PI*Re)
!     temp2 = (1d0-1d0/densratio)*gravity_const + 3d0*dragNm1/(densratio*4d0*PI*Re)
!     temp3 = (1d0-1d0/densratio)*gravity_const + 3d0*dragNm2/(densratio*4d0*PI*Re)
    
    temp1 = 0.75d0*densratio*drag/(PI*Re)
!     if (timeN.lt.2.5d0*deltat) then ! fix problem with very high initial drag jump
!       temp2=temp1
!       temp3=temp1
!     else
      temp2 = 0.75d0*densratio*dragNm1/(PI*Re)
      temp3 = 0.75d0*densratio*dragNm2/(PI*Re)
!     endif

!First order:
!     V_sphere = V_sphere + ((1d0-densratio)*gravity_const + temp1)*deltat
! AB3:
!     V_sphere = V_sphere + (1d0-densratio)*gravity_const*deltat + (deltat/12d0)*(23d0*temp1 - 16d0*temp2 + 5d0*temp3)
! AB2:
    V_sphere = V_sphere + (1d0-densratio)*gravity_const*deltat + (0.5d0*deltat)*(3d0*temp1 - temp2)

!           print*,drag,dragNm1,dragNm2, V_sphere,V_sphereNm1,V_sphereNm2
! BDF2:
!     V_sphere = (deltat/time_gamma_0)*(time_beta_0*temp1 + time_beta_1*temp2) + time_alpha_0*V_sphereNm1 + time_alpha_0*V_sphereNm2
    
  END SUBROUTINE calcSphereVelocity
  
  

END MODULE movingSphere_module