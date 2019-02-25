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

MODULE OIFS_module
  USE shared_data
  USE functions_module
  IMPLICIT NONE
  CONTAINS

! Not used:
!   SUBROUTINE calcOIFS_stress
!     IMPLICIT NONE
!     integer:: i
!     DOUBLE PRECISION, DIMENSION(nptot) :: temp1,temp2
! 
! !xx part:    
!     CALL RK4_OIFS(Txx,h_RK4,RK4_timesteps,timeN,temp1)
!     CALL RK4_OIFS(TxxNm1,h_RK4,2*RK4_timesteps,timeNm1,temp2)
!     OIFSstress_xx = Wetime_constant2*(time_alpha_0*temp1 + time_alpha_1*temp2)
! 
! !xy part:
!     CALL RK4_OIFS(Txy,h_RK4,RK4_timesteps,timeN,temp1)
!     CALL RK4_OIFS(TxyNm1,h_RK4,2*RK4_timesteps,timeNm1,temp2)
!     OIFSstress_xy = Wetime_constant2*(time_alpha_0*temp1 + time_alpha_1*temp2)
!   
! !yy part:
!     CALL RK4_OIFS(Tyy,h_RK4,RK4_timesteps,timeN,temp1)
!     CALL RK4_OIFS(TyyNm1,h_RK4,2*RK4_timesteps,timeNm1,temp2)
!     OIFSstress_yy = Wetime_constant2*(time_alpha_0*temp1 + time_alpha_1*temp2)
! 
! !zz part:
!     CALL RK4_OIFS(Tzz,h_RK4,RK4_timesteps,timeN,temp1)
!     CALL RK4_OIFS(TzzNm1,h_RK4,2*RK4_timesteps,timeNm1,temp2)
!     OIFSstress_zz = Wetime_constant2*(time_alpha_0*temp1 + time_alpha_1*temp2)
! 
!   END SUBROUTINE calcOIFS_stress
  
  SUBROUTINE applyOIFS_to_stokes
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(nptot) :: OIFS_soln_x,OIFS_soln_y,&
            temp1,temp2,temp3
    
! Calculate the RHS contribution from the u^tilde and u^tilde,tilde terms
! x part:
    CALL RK4_OIFS(V_x,h_RK4,RK4_timesteps,timeN,temp1)
    CALL RK4_OIFS(V_xNm1,h_RK4,2*RK4_timesteps,timeNm1,temp2)
!     CALL RK4_OIFS(V_xNm2,h_RK4,3*RK4_timesteps,timeNm2,temp3) ! doesn't work
    OIFS_soln_x= time_alpha_0*temp1 + time_alpha_1*temp2 ! + time_alpha_2*temp3
! y part:
    CALL RK4_OIFS(V_y,h_RK4,RK4_timesteps,timeN,temp1)
    CALL RK4_OIFS(V_yNm1,h_RK4,2*RK4_timesteps,timeNm1,temp2)
!     CALL RK4_OIFS(V_xNm2,h_RK4,3*RK4_timesteps,timeNm2,temp3) ! doesn't work
    OIFS_soln_y= time_alpha_0*temp1 + time_alpha_1*temp2 ! + time_alpha_2*temp3

! Integrate so it fits into the weak form of our problem:
    CALL integrateOIFS_solution_stokes(Mv_x,OIFS_soln_x,OIFS_stokes_contrib_x)
    CALL integrateOIFS_solution_stokes(Mv_y,OIFS_soln_y,OIFS_stokes_contrib_y)    
! OLD VERSION (Re/We):

    f_x = f_x + Retime_constant2*OIFS_stokes_contrib_x
    f_y = f_y + Retime_constant2*OIFS_stokes_contrib_y

  END SUBROUTINE applyOIFS_to_stokes
  


!   SUBROUTINE integrateOIFS_solution_stokes(massmatrix,in,out)
!     IMPLICIT NONE
!     INTEGER :: el,ij
!     DOUBLE PRECISION, DIMENSION(nptot) :: in,out,tempglob
!     DOUBLE PRECISION, DIMENSION(0:NP1SQM1) :: templocal
!     DOUBLE PRECISION, DIMENSION(0:NP1SQM1,0:NP1SQM1,numelm) :: massmatrix
! ! Integrate the final solutions in each element
!       out=0d0
!       DO el=1,numelm
!   DO ij=0,NP1SQM1
!     out(mapg(ij,el)) = out(mapg(ij,el)) + in(mapg(ij,el))*massmatrix(ij,ij,el)
!   ENDDO
!       ENDDO
!       DO ij=1,npedg
!   IF (mult(ij).ne.1) THEN
!     out(ij)=out(ij)!/dfloat(mult(i)) ! needed? - check!!!
!   ENDIF
!       ENDDO
!   
!     END SUBROUTINE integrateOIFS_solution_stokes

! OLD VERSION USING vector restriction/prolongation
  SUBROUTINE integrateOIFS_solution_stokes(massmatrix,in,out)
    IMPLICIT NONE
    INTEGER :: el,kl
    DOUBLE PRECISION, DIMENSION(nptot) :: in,out,tempglob
    DOUBLE PRECISION, DIMENSION(0:NP1SQM1) :: templocal
    DOUBLE PRECISION, DIMENSION(0:NP1SQM1,0:NP1SQM1,numelm), INTENT(IN) :: massmatrix
! Integrate the final solutions in each element
    out=0d0
    DO el=1,numelm
      templocal=0d0  
      CALL veclocalrestrict(in,el,templocal)
      DO kl=0,NP1SQM1
        templocal(kl)=templocal(kl)*massmatrix(kl,kl,el)
      ENDDO
      tempglob=0d0
      CALL vecglobalprolongation(templocal,el,tempglob)
      out = out + tempglob
    ENDDO
  END SUBROUTINE integrateOIFS_solution_stokes

  
  SUBROUTINE RK4_OIFS(U_in,h,M,tstart,soln)
! Use the RK4 method to solve the time dependent problem
! represented by:
! dG/dt = -UStar * grad(G)
!
! where G has initial condition G0 = U_in
!
! and Ustar = ((t - t_n-1)/Dt)*U_n + (1-(t-t_n-1)/Dt)*U_n-1
!
! where U_i is the velocity solution at outer timestep i
!
! For RK4, the initial condition is the vector [U_in,V_in]
!
! The initial time value is tstart (typically t_n-1 or t_n depending
! on which problem is being solved.
!
!
! Note the - mesh_velocity term: this is due to the use of ALE - see Giancarlo Russo's Thesis. (pdf pg 83)
!
    IMPLICIT NONE
    INTEGER :: i,j,M,k
    DOUBLE PRECISION :: h,t_curr,tstart,t_n
      
      
    DOUBLE PRECISION, DIMENSION(nptot) :: U_in,&
            k1,k2,k3,k4,&
            dkdx,dkdy,&
            soln,tempvec
    
    DOUBLE PRECISION, DIMENSION(nptot,2) :: Ustar

! Setup initial parameters
!     Dt=dfloat(M)*h
    t_n=tstart ! set inner time to be the beginning outer time step
    soln=U_in

! Iterate over M inner timesteps of RK4
    DO i=1,M 
! Calculate k1
      CALL calcUstar(t_n,Ustar)
!       CALL calc_Upwinding_nodes(Ustar(1:nptot,1),UStar(1:nptot,2)) 
!       CALL calc_upwinded_GradVec(soln,dkdx,dkdy)
      CALL calcGradVec(soln,dkdx,dkdy)
      
      DO k=1,nptot
        k1(k)=-((Ustar(k,1))*dkdx(k) + Ustar(k,2)*dkdy(k))
      ENDDO

! Calculate k2
      CALL calcUstar(t_n+h/2d0,Ustar)
!       CALL calc_Upwinding_nodes(Ustar(1:nptot,1),UStar(1:nptot,2))
!       CALL calc_upwinded_GradVec(soln+(h*k1/2d0),dkdx,dkdy)
      tempvec=soln+(h*k1/2d0)
      CALL calcGradVec(tempvec,dkdx,dkdy)
      
      DO k=1,nptot
        k2(k)=-((Ustar(k,1))*dkdx(k) + Ustar(k,2)*dkdy(k))
      ENDDO

! Calculate k3
!       CALL calc_upwinded_GradVec(soln+(h*k2/2d0),dkdx,dkdy)
      tempvec=soln+(h*k2/2d0) 
      CALL calcGradVec(tempvec,dkdx,dkdy)
      
      DO k=1,nptot  
        k3(k)=-((Ustar(k,1))*dkdx(k) + Ustar(k,2)*dkdy(k))
      ENDDO

! Calculate k4
      CALL calcUstar(t_n+h,Ustar)
!       CALL calc_Upwinding_nodes(Ustar(1:nptot,1),UStar(1:nptot,2))
!       CALL calc_upwinded_GradVec(soln+h*k3,dkdx,dkdy)
      tempvec=soln+h*k3
      CALL calcGradVec(tempvec,dkdx,dkdy)

      DO k=1,nptot    
        k4(k)=-((Ustar(k,1))*dkdx(k) + Ustar(k,2)*dkdy(k))
      ENDDO

      soln=soln+h*(k1+2d0*k2+2d0*k3+k4)/6d0
      DO k=1,nptot
        IF (IsNaN(soln(k))) THEN
          print*,'OIFS ERROR!',timeN,k,t_n,k1(k),k2(k),k3(k),k4(k),U_in(k),Ustar(k,1),Ustar(k,2),dkdx(k),dkdy(k)
          STOP
        ENDIF
      ENDDO

      t_n=t_n+h
    ENDDO
  END SUBROUTINE RK4_OIFS

  


! OLD - NOT USING RESTRICT/PROLONG ROUTINES.
!   SUBROUTINE calcGradVec(in,outx,outy)
! ! Calculates the gradient of one component of a vector field
! ! This is done by restricting into each local element according to the upwinding scheme.
! !
! ! E.G. if we have G(u,v), then input would be u at all GL points
! !      and outputs would be:
! !      outx = the differential of u w.r.t x at all GL points
! !      outy = the differential of u w.r.t y at all GL points
! 
!     IMPLICIT NONE
!     INTEGER :: el,i,j,ij,p,pj,ip
!     DOUBLE PRECISION :: temp1,temp2
!     DOUBLE PRECISION, DIMENSION(nptot) :: in,outx,outy
! 
!   
!     outx=0d0
!     outy=0d0
!     DO el=1,numelm
!      DO j=0,N
!       DO i=0,N
!   ij=i+j*NP1
!   IF (upwinded_element(mapg(ij,el)).ne.el) CYCLE
!   temp1=0d0
!   temp2=0d0
!   DO p=0,N
!     pj=p+j*NP1
!     ip=i+p*NP1
!     temp1 = temp1 + in(mapg(pj,el))*d(i,p)
!     temp2 = temp2 + in(mapg(ip,el))*d(j,p) ! switched p and q to save a new loop
!   ENDDO
!   outx(mapg(ij,el)) = (dyde(i,j,el)*temp1 - dydp(i,j,el)*temp2)/jac(i,j,el)
!   outy(mapg(ij,el)) = (dxdp(i,j,el)*temp2 - dxde(i,j,el)*temp1)/jac(i,j,el)
!       ENDDO
!     ENDDO
!     ENDDO
! 
!     
!   END SUBROUTINE calcGradVec

 
  SUBROUTINE calcUstar(t_in,out)
! Calculates U^* as defined in Van Os' thesis pg 71 (pdf pg 93)
! Output is Ustar in both vector components at all GL points
! This is done in order to avoid problems in passing different old and current
! velocities - see difference between the two problems, tilde and 2tilde
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: t_in
    DOUBLE PRECISION, DIMENSION(nptot,2), INTENT(OUT) :: out

    out(1:nptot,1) = (t_in - timeNm1)*(V_x(1:nptot)- mesh_velocity(1:nptot))/deltat + (1d0 - (t_in-timeNm1)/deltat)*(V_xNm1(1:nptot)- mesh_velocity(1:nptot))
    out(1:nptot,2) = (t_in - timeNm1)*V_y(1:nptot)/deltat + (1d0 - (t_in-timeNm1)/deltat)*V_yNm1(1:nptot)

! trying to get it working for 3rd order.. not working atm:
!   IF ( t_in.ge.timeNm1) THEN
!     out(1:nptot,1) = (t_in - timeNm1)*V_x(1:nptot)/deltat + (1d0 - (t_in-timeNm1)/deltat)*V_xNm1(1:nptot)
!     out(1:nptot,2) = (t_in - timeNm1)*V_y(1:nptot)/deltat + (1d0 - (t_in-timeNm1)/deltat)*V_yNm1(1:nptot)
!   ELSE
!     out(1:nptot,1) = (t_in - timeNm2)*V_xNm1(1:nptot)/deltat + (1d0 - (t_in-timeNm2)/deltat)*V_xNm2(1:nptot)
!     out(1:nptot,2) = (t_in - timeNm2)*V_yNm1(1:nptot)/deltat + (1d0 - (t_in-timeNm2)/deltat)*V_yNm2(1:nptot)
!   ENDIF
!  OR?:
!   out(1:nptot,1) = sqrt((t_in - timeNm1)/deltat*((t_in-timeNm2)/deltat -1d0))*V_x(1:nptot)+ sqrt((1d0 - (t_in-timeNm1)/deltat)*(t_in-timeNm2)/deltat)*V_xNm1(1:nptot) &
!         + sqrt((t_in-timeNm1)/deltat*((t_in-timeNm2)/deltat -1))*V_xNm2(1:nptot)
!   out(1:nptot,2) = sqrt((t_in - timeNm1)/deltat*((t_in-timeNm2)/deltat -1d0))*V_y(1:nptot)+ sqrt((1d0 - (t_in-timeNm1)/deltat)*(t_in-timeNm2)/deltat)*V_yNm1(1:nptot) &
!         + sqrt((t_in-timeNm1)/deltat*((t_in-timeNm2)/deltat -1))*V_yNm2(1:nptot)
  END SUBROUTINE calcUstar

  SUBROUTINE setup_time_constants
    IMPLICIT NONE
    
    IF (param_time_order.eq.1) THEN
      time_gamma_0=1d0
      time_alpha_0=1d0

      time_beta_0=1d0
! All others zero:
      time_alpha_1=0d0
      time_alpha_2=0d0
      time_beta_1=0d0
      time_beta_2=0d0
    ELSEIF (param_time_order.eq.2) THEN
      time_gamma_0=3d0/2d0
      time_alpha_0=2d0
      time_alpha_1=-1d0/2d0
      time_beta_0=2d0
      time_beta_1=-1d0
! All others zero:
      time_alpha_2=0d0
      time_beta_2=0d0
    ELSEIF (param_time_order.eq.3) THEN
      time_gamma_0=11d0/6d0
      time_alpha_0=3d0
      time_alpha_1=-3d0/2d0
      time_alpha_2=1d0/3d0
      time_beta_0=3d0
      time_beta_1=-3d0
      time_beta_2=1d0
    ENDIF
  END SUBROUTINE setup_time_constants
  
END MODULE OIFS_module
