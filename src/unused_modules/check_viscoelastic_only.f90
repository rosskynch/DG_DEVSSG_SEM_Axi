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

PROGRAM check_viscoelastic_only
  USE constants
  USE shared_data
  USE geometry_module
  USE SEM_module
  USE boundary_module
  USE initial_conditions_module
!   USE mkl_solver ! Uncomment to use MKL library - see others later in prog.
  USE solver
  USE OIFS_module
  USE viscoelastic_module
  USE utils
  USE IO_module
  USE result_analysis
  USE movingSphere_module
  IMPLICIT NONE
  INTEGER :: i,j,ij,k,l,kl,ijkl,el,tempint,minresconv,internalij,kk,ll,&
	     timestep,numtimesteps,rowcount,printoutcount,stressnode,velnode
  DOUBLE PRECISION :: test,msqe,temp,cputime1,cputime2,&
		      V_xH1norm,V_yH1norm,pressure_L2norm

! Read in values (This will assign required memory)
  CALL read_input
! Create local to global map
  CALL create_global_map
! assign remaining memory
  CALL assignMem2
  
! avoid error caused by dealloc later
  ALLOCATE(global_matrix(1,1),globalRHS(1),globalSOL(1))
! If using MKL, comment out above and uncomment below:
!  CALL assignGlobalDimension

! GLOBAL PARAMETERS - MUST ADD TO INPUT LATER
  param_beta = 1d0/9d0
  We=1d0
  Re=1d0

  IF (Re.lt.1d-6.AND.movingmeshflag.eq.1) THEN
    print*,'ERROR Moving mesh with Re = 0?'
    STOP
  ENDIF

! Set up parameters and initial values:
  param_alphaZ = 5d-3
  rho_s = 5d0
  rho_f = 1d0
  gravity_const = 9.8d0
  
  deltat = 1d-3

  
  RK4_timesteps=8
  h_RK4=deltat/dfloat(RK4_timesteps)
  
  
!   timeN=2d0*deltat
!   timeNm1=time-deltat
!   timeNm2=timeNm1-deltat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TESTING SECTION
! DONT FORGET TO DELETE THIS!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  timeN=deltat
  timeNm1=0d0

  CALL initialiseSEM
  CALL calcNodeCoords
  CALL update_geometry
  CALL updateSEM
  
  DO i=1,nptot
    V_xNm1(i) = calcWaters_analytical_U(nodeCoord(i,1),nodeCoord(i,2),timeNm1)
    V_yNm1(i) = 0d0
    
    V_x(i) = calcWaters_analytical_U(nodeCoord(i,1),nodeCoord(i,2),timeN)
    V_y(i) = 0d0
    
    TxxNm1(i) = calcWaters_analytical_Txx(nodeCoord(i,1),nodeCoord(i,2),timeNm1)
    TxyNm1(i) = calcWaters_analytical_Txy(nodeCoord(i,1),nodeCoord(i,2),timeNm1)
    TyyNm1(i) = 0d0
    
    Txx(i) = calcWaters_analytical_Txx(nodeCoord(i,1),nodeCoord(i,2),timeN)
    Txy(i) = calcWaters_analytical_Txy(nodeCoord(i,1),nodeCoord(i,2),timeN)
    Tyy(i) = 0d0
  ENDDO
  
velnode = mapg((N+1)*N/2 + N/2,numelm)
stressnode=2
  DO j=1,50000
!     CALL calcOIFS_stress
    DO i=1,nptot
      OIFSstress_xx(i)=We/(2d0*deltat)*(-4d0*Txx(i) + TxxNm1(i))
      OIFSstress_xy(i)=We/(2d0*deltat)*(-4d0*Txy(i) + TxyNm1(i))
      OIFSstress_yy(i)=We/(2d0*deltat)*(-4d0*Tyy(i) + TyyNm1(i))
    ENDDO
    
    
    timeNm1=timeN
    timeN=timeN+deltat
    V_xNm1=V_x
    V_yNm1=V_y
    
    DO i=1,nptot
      V_x(i) = calcWaters_analytical_U(nodeCoord(i,1),nodeCoord(i,2),timeN)
      V_y(i) = 0d0
      uxy(i) = calcWaters_analytical_dUdy(nodeCoord(i,1),nodeCoord(i,2),timeN)
      uxx(i) = 0d0
      uyx(i) = 0d0
      uyy(i) = 0d0
    ENDDO
!     CALL calc_Upwinding_nodes(V_x,V_y)
!     CALL calc_upwinded_GradVec(V_x,uxx,uxy)
!     CALL calc_upwinded_GradVec(V_y,uyx,uyy)
    CALL applyElasticStress
    print*,timeN,&
	    V_x(velnode),calcWaters_analytical_U(nodeCoord(velnode,1),nodeCoord(velnode,2),timeN),&
	    Txx(2),calcWaters_analytical_Txx(0d0,0d0,timeN),&
	    Txy(2),calcWaters_analytical_Txy(0d0,0d0,timeN),&
	    Tyy(2),0d0,&
	    uxy(2),calcWaters_analytical_dUdy(nodeCoord(stressnode,1),nodeCoord(stressnode,2),timeN)
    
      
      
  ENDDO

END PROGRAM check_viscoelastic_only