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

MODULE functions_module
  USE shared_data
  USE constants
  IMPLICIT NONE
  CONTAINS
  
!   SUBROUTINE update_sphere_values
! ! NOT CURRENTLY NEEDED
! ! Updates the maximum and minimum values of X & Y co-ordinates for the sphere
! !
!     IMPLICIT NONE
!     
!     sphereXleft = centre_of_sphere(1)-rad_sphere
!     sphereXright = centre_of_sphere(1)+rad_sphere
!     sphereYtop = centre_of_sphere(2)+rad_sphere
!     sphereYbottom = centre_of_sphere(2)-rad_sphere
!     
!   END SUBROUTINE update_sphere_values



  DOUBLE PRECISION FUNCTION interpolate_solution_hbasis(el_in,soln_in,xi_in,eta_in)
    IMPLICIT NONE
    INTEGER :: el_in,i,j,ij
    DOUBLE PRECISION :: soln_in(nptot),xi_in,eta_in,temp
    
    temp=0d0
    DO i=0,N
      DO j=0,N
        ij=i+j*(N+1)
        temp = temp + soln_in(mapg(ij,el_in))*hbasis(i,xi_in)*hbasis(j,eta_in)
      ENDDO
    ENDDO
    interpolate_solution_hbasis=temp
    RETURN
  END FUNCTION interpolate_solution_hbasis
  
!! BOTH OF THE FOLLOWING ROUTINES NEED PROPERLY CHECKING !!
  DOUBLE PRECISION FUNCTION hbasis(i,x_in)
    IMPLICIT NONE
    INTEGER :: i,j
    DOUBLE PRECISION :: S,x_in

!       if(xz.lt.-(1d0+eps))then!!.or.(xz.gt.(1d0+eps)))then
!        !!  write(*,*)'h out of range',xz
!        xz = -1d0
!       elseif(xz.gt.(1d0+eps))then
!        xz = 1d0
! 
!       else
    S=1d0
    DO j=0,N
      IF ( gl(i).ne.gl(j) ) S=S*(x_in-gl(j))/(gl(i)-gl(j))
    ENDDO
    hbasis=S
!       ENDIF
    RETURN
  END FUNCTION hbasis
  
  DOUBLE PRECISION FUNCTION hbasis_tilde(i,x_in)
    IMPLICIT NONE
    INTEGER          :: i,j
    DOUBLE PRECISION :: S,x_in

!     if(xz.lt.-(1d0+eps))then!!.or.(xz.gt.(1d0+eps)))then
!      !!  write(*,*)'h out of range',xz
!      xz = -1d0
!     elseif(xz.gt.(1d0+eps))then
!      xz = 1d0
!   
!     else
    S=1d0
    DO j=1,N-1
      IF ( gl(i).ne.gl(j) ) S=S*(x_in-gl(j))/(gl(i)-gl(j))
    ENDDO
    hbasis_tilde=S
!     ENDIF
    RETURN
  END FUNCTION hbasis_tilde  
  
  
  DOUBLE PRECISION FUNCTION calc_gl_htilde(i,j)
! NOTE ORDER OF i AND j MAY DIFFER FROM NOTATION IN THESIS.
!  calculates tilde{h_i(xj)} - used in constructB routine. (TILDE)
    IMPLICIT NONE
    INTEGER :: i,j

    IF (i.eq.0.OR.i.eq.N) THEN
      print*,'ERROR: Called non-existant pressure test function!'
      STOP
    ELSEIF (j.eq.0.OR.j.eq.N) THEN
      calc_gl_htilde = -(1d0 - gl(i)*gl(i))*Leg1(j)/(dfloat(N)*(dfloat(N) + 1d0)*Leg(i)*(gl(j) - gl(i)))
    ELSEIF (i.eq.j) THEN
      calc_gl_htilde = 1d0
    ELSE
      calc_gl_htilde = 0d0
    ENDIF
    RETURN
  END FUNCTION calc_gl_htilde
  
! Functions to define the RHS term, f:
  DOUBLE PRECISION FUNCTION fterm_x(i) RESULT(a)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
! ADD MORE AS APPROPRIATE.    
    IF     (param_function_choice.eq.3) THEN
      a = model_soln_f_x(i)
  ELSEIF (param_function_choice.eq.10) THEN
    a = model_soln2_f_x(i)
    ELSE ! ALL OTHERS HAVE NO BODY FORCE
      a = 0d0
    ENDIF
    RETURN
  END FUNCTION fterm_x
  
  DOUBLE PRECISION FUNCTION fterm_y(i) RESULT(a)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
! ADD MORE AS APPROPRIATE.
    IF     (param_function_choice.eq.3) THEN
      a = model_soln_f_y(i)
    ELSEIF (param_function_choice.eq.10) THEN
      a = model_soln2_f_y(i)
    ELSE ! ALL OTHERS HAVE NO BODY FORCE
      a = 0d0
    ENDIF
    RETURN
  END FUNCTION fterm_y

  SUBROUTINE veclocalrestrict(in,k,out)

! uses mapg to map from a global vector of size nptot
! to a local vector of size (N+1)^2
    IMPLICIT NONE
    INTEGER :: ij,k
    DOUBLE PRECISION, DIMENSION(nptot) :: in
    DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1) :: out

    DO ij=0,(N+1)**2-1
      out(ij)=in(mapg(ij,k))
    ENDDO
  END SUBROUTINE veclocalrestrict

  SUBROUTINE veclocalrestrict_internal_nodes(in,k,out)

! uses mapg to map from a global vector of size nptot
! to a local vector of size (N+1)^2
    IMPLICIT NONE
    INTEGER :: i,j,ij,k
    DOUBLE PRECISION, DIMENSION(npint) :: in
    DOUBLE PRECISION, DIMENSION(1:(N-1)**2) :: out

    DO j=1,N-1
      DO i=1,N-1
        ij = i+(j-1)*(N-1)
        out(ij)=in(mapg_pressure(ij,k))
      ENDDO
    ENDDO
  END SUBROUTINE veclocalrestrict_internal_nodes

  SUBROUTINE matlocalrestrict(in,k,out)

! uses mapg to map from a global matrix of size nptot by nptot
! to a local matrix of size (N+1)^2 by (N+1)^2
    IMPLICIT NONE
    INTEGER :: ij,kl,k
    DOUBLE PRECISION, DIMENSION(nptot,nptot) :: in
    DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1,0:(N+1)**2-1) :: out
  
    DO kl=0,(N+1)**2-1
      DO ij=0,(N+1)**2-1
        out(ij,kl)=in(mapg(ij,k),mapg(kl,k))
      ENDDO
    ENDDO
  END SUBROUTINE matlocalrestrict

  SUBROUTINE matglobalprolongation(in,k,out)

! uses mapg to map from a local vector of size (N+1)^2
! to a global vector of size nptot
    IMPLICIT NONE
    INTEGER :: ij,kl,k
    DOUBLE PRECISION, DIMENSION(nptot,nptot) :: out
    DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1,0:(N+1)**2-1) :: in
  
    DO ij=0,(N+1)**2-1
      DO kl=0,(N+1)**2-1
        out(mapg(ij,k),mapg(kl,k))=in(ij,kl)
      ENDDO
    ENDDO  
  END SUBROUTINE matglobalprolongation

  SUBROUTINE vecglobalprolongation(in,k,out)

! uses mapg to map from a local vector of size (N+1)^2
! to a global vector of size nptot
    IMPLICIT NONE
    INTEGER :: ij,k
    DOUBLE PRECISION, DIMENSION(nptot) :: out
    DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1) :: in

    DO ij=0,(N+1)**2-1
      out(mapg(ij,k))=in(ij)
    ENDDO
  END SUBROUTINE vecglobalprolongation

  SUBROUTINE vecglobalprolongation_internal_nodes(in,k,out)

! uses mapg to map from a local vector of size (N-1)^2
! to a global vector of size npint - ie pressure term nodes
    IMPLICIT NONE
    INTEGER :: i,j,ij,internalij,k,temp
    DOUBLE PRECISION, DIMENSION(npint) :: out
    DOUBLE PRECISION, DIMENSION(1:(N-1)**2) :: in

    temp=nptot-npint
    DO j=1,N-1
      DO i=1,N-1
        ij=i+j*(N+1)
        internalij=i+(j-1)*(N-1)
        out(mapg(ij,k)-temp)=in(internalij)
      ENDDO
    ENDDO
  END SUBROUTINE vecglobalprolongation_internal_nodes

! routine only for u actual in order to avoid overlap of actual solution between elements
  SUBROUTINE vecglobalprolongationUACT(in,k,out)

! uses mapg to map from a local vector of size (N+1)^2
! to a global vector of size nptot
    IMPLICIT NONE
    INTEGER :: ij,k,np
    DOUBLE PRECISION, DIMENSION(nptot) :: out
    DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1) :: in
  
    DO ij=0,(N+1)**2-1
      np=mapg(ij,k)
      out(np)=in(ij)
      IF (mult(np).ne.1) THEN
        out(np)=out(np)/dfloat(mult(np))
      ENDIF
    ENDDO
  END SUBROUTINE vecglobalprolongationUACT
  
  SUBROUTINE calc_Upwinding_nodes(in_x,in_y)
! Decides which element should be used within the upwinding scheme according to the input velocity field)
    IMPLICIT NONE
    INTEGER :: el,i,j,ij,local_vertex,k
    DOUBLE PRECISION :: in_angle,tempx,tempy,temp,ang1,ang2
    DOUBLE PRECISION, DIMENSION(nptot) :: in_x,in_y

    upwinded_element=0
    
    DO el=1,numelm
      ! All internal nodes have multiplicity 1:    
      DO i=1,N
        DO j=1,N
          ij=i+j*(N+1)
          k=mapg(ij,el)
          upwinded_element(k)=el
        ENDDO
      ENDDO

! First deal with the edge nodes (non-vertex)      
! Edge 1 
      DO ij=1,N-1 !(j=0)
        k=mapg(ij,el)
        IF (mult(k).eq.1) THEN
          upwinded_element(k)=el
        ELSE
          temp = norm_to_edge(1,1,el)*in_x(k) + norm_to_edge(2,1,el)*in_y(k)
          IF (temp.ge.0d0) THEN
! CHECK !
! We take less than or equal here because if the velocity runs along a shared edge, then it does not matter which element
! we use, they should have the same gradient (ie, the last element checked will be the one we use)
! CHECK !
            upwinded_element(k)=el
          ENDIF
        ENDIF
      ENDDO
! Edge 2
      i=N
      DO j=1,N-1
        ij=i+j*(N+1)
        k=mapg(ij,el)
        IF (mult(k).eq.1) THEN
          upwinded_element(k)=el
        ELSE
          temp = norm_to_edge(1,2,el)*in_x(k) + norm_to_edge(2,2,el)*in_y(k)
          IF (temp.ge.0d0) THEN
! see above
            upwinded_element(k)=el
          ENDIF
        ENDIF
      ENDDO
! Edge 3
      j=N
      DO i=1,N-1
        ij=i+j*(N+1)
        k=mapg(ij,el)
        IF (mult(k).eq.1) THEN
          upwinded_element(k)=el
        ELSE
          temp = norm_to_edge(1,3,el)*in_x(mapg(ij,el)) + norm_to_edge(2,3,el)*in_y(mapg(ij,el))
          IF (temp.ge.0d0) THEN
! see above
            upwinded_element(mapg(ij,el))=el
          ENDIF
        ENDIF
      ENDDO
! Edge 4
      DO j=1,N-1
        ij=j*(N+1)
        k=mapg(ij,el)
        IF (mult(k).eq.1) THEN
          upwinded_element(k)=el
        ELSE
          temp = norm_to_edge(1,4,el)*in_x(mapg(ij,el)) + norm_to_edge(2,4,el)*in_y(mapg(ij,el))
          IF (temp.ge.0d0) THEN
! see above
            upwinded_element(mapg(ij,el))=el
          ENDIF
        ENDIF
      ENDDO

! Now deal with vertices of elements
! These are dealt with by checking if the angle of the velocity lies within the bounds of the 2 edges
! which start at the vertex
      DO local_vertex=1,4
        k=node(el,local_vertex)
        IF (mult(k).eq.1) THEN
          upwinded_element(k)=el
          CYCLE
        ENDIF
! Set k to be the global node number of the vertex:
  
        IF (inflowflag(k)) THEN
          tempx=in_x(k)
          tempy=in_y(k)
        ELSEIF(wallsymmflag(k)) THEN
          tempx=-in_x(k)
          tempy=0d0!in_y(k) !should be 0 anyway
        ELSE
          tempx=-in_x(k)
          tempy=-in_y(k)
        ENDIF
    
! check if magnitude of the velocity is zero.. I.e the velocity is zero.
        IF (SQRT(tempx**2+tempy**2).lt.1d-15) THEN
          IF (upwinded_element(k).eq.0) THEN
            upwinded_element(k)=el
            CYCLE
          ELSE
! Check upwinded node on either side of the vertex node, should keep the upwinding grid less jagged.
            IF (local_vertex.eq.1) THEN
              IF( upwinded_element(mapg(NP1,el)).eq.el.and.upwinded_element(mapg(1,el)).eq.el ) THEN
          upwinded_element(k)=el
              ENDIF
            ELSEIF (local_vertex.eq.2) THEN
              IF( upwinded_element(mapg(N-1,el)).eq.el.and.upwinded_element(mapg(N+NP1,el)).eq.el ) THEN
          upwinded_element(k)=el
              ENDIF
            ELSEIF (local_vertex.eq.3) THEN
              IF( upwinded_element(mapg(N+NM1*NP1,el)).eq.el.and.upwinded_element(mapg(N-1+N*NP1,el)).eq.el ) THEN
          upwinded_element(k)=el
              ENDIF
            ELSEIF (local_vertex.eq.4) THEN
              IF( upwinded_element(mapg(1+N*NP1,el)).eq.el.and.upwinded_element(mapg(NM1*NP1,el)).eq.el ) THEN
          upwinded_element(k)=el
              ENDIF
            ENDIF
            CYCLE
          ENDIF
        ENDIF

        IF (tempx.gt.0d0) THEN
          in_angle = atan(tempy/tempx)
        ELSEIF (tempx.lt.0d0.AND.tempy.ge.0d0) THEN
          in_angle = atan(tempy/tempx) + PI
        ELSEIF (tempx.lt.0d0.AND.tempy.lt.0d0) THEN
          in_angle = atan(tempy/tempx) - PI
        ELSEIF (tempx.eq.0d0.AND.tempy.gt.0d0) THEN
          in_angle = PI/2d0
        ELSEIF (tempx.eq.0d0.AND.tempy.lt.0d0) THEN
          in_angle = -PI/2d0
        ELSE
          in_angle = 0d0
        ENDIF

        IF (k.eq.3) THEN
          print*,k,el,tempx,tempy,in_angle,angles_at_vertex(1,local_vertex,el),angles_at_vertex(2,local_vertex,el)
          print*,(angles_at_vertex(1,local_vertex,el)-in_angle).lt.-1d-16,(angles_at_vertex(2,local_vertex,el)-in_angle).gt.1d-16
          print*,ABS(in_angle).eq.PI,timeN
        ENDIF
!     stop

        ang1=angles_at_vertex(1,local_vertex,el)
        ang2=angles_at_vertex(2,local_vertex,el)
        IF (in_angle.ge.ang1.AND.in_angle.le.ang2) THEN
          upwinded_element(k)=el
        ELSEIF (ang1.ge.0d0.AND.ang2.lt.0d0) THEN ! ie. the largest angle is actually first (where the two angles are between the PI/-PI divide)
          IF (in_angle.le.ang2.OR.in_angle.gt.ang1) THEN
            upwinded_element(k)=el
          ENDIF
        ENDIF
!         ELSEIF (ABS(in_angle-PI).lt.1d-15) THEN
!           IF (-in_angle.ge.ang1-1d-15.AND. &
!               -in_angle.le.ang2+1d-15) THEN
!             upwinded_element(k)=el
!           ENDIF
      ENDDO
    ENDDO
    
! Check that all nodes have an associated element via upwinding:
    DO i=1,nptot
      IF (upwinded_element(i).eq.0) THEN
        print*,'ERROR: Upwinding Scheme has no element for node',i,'time is',timeN
        print*,nodeCoord(i,1),nodeCoord(i,2),V_x(i),V_y(i)
        STOP
      ENDIF
    ENDDO
  END SUBROUTINE calc_Upwinding_nodes
  
  SUBROUTINE calc_Upwinding_local_edge_nodes(in_x,in_y)
    IMPLICIT NONE
    INTEGER :: el,i,ij,k,edge!,edgeM1,edgeP1
    DOUBLE PRECISION :: in_x(nptot),in_y(nptot),&
      temp,tol
      
    tol = 1d-6
    upwind_local_edge_node = -1
    DO el=1,numelm
      DO ij=0,NP1SQM1
        upwinded_element(mapg(ij,el))=el
      ENDDO
    ENDDO
! Calculates upwinding for all edges according to the dot product of the velocity and outward normal.
    DO el=1,numelm
      DO edge=1,4
        DO i=0,N
          ij=local_edge_node(i,edge)
          upwind_local_edge_node(i,edge,el)=el
          k=mapg(ij,el)
          IF (.not.inflowflag(k).and.mult(k).ne.1) THEN
! calculate U.N where N is the outward normal to the edge
            temp = in_x(k)*norm_to_edge(1,edge,el) + in_y(k)*norm_to_edge(2,edge,el)
            IF (temp.lt.-tol) THEN
              IF (conelm(el,edge).ne.0) upwind_local_edge_node(i,edge,el) = conelm(el,edge)
            ENDIF
          ENDIF
          IF (upwinded_element(k).eq.el) upwinded_element(k)=upwind_local_edge_node(i,edge,el)
        ENDDO
      ENDDO
    ENDDO

! Additional array which stores the 1-D upwinding for an axisymmetric wall:   
! In 1-D the upwinding is only needed for 2 points: at the start and end "vertices"
! We assume that the anti-clockwise numbering is matched in the integration.
!
! this may mean that the boundary terms are backwards in some cases. Although not for the axisymmetric wall (where the edge is always at the "bottom").


! FAIRLY SURE THIS ROUTINE IS WRONG SOMEHOW..
! FOR NOW AM USING STRONG FORM OF CONS EQN SO ITS NOT NEEDED:
!     IF (coordflag.eq.1) THEN
!       DO el=1,numelm
!   DO edge=1,4
!     upwind_wallsymm_1d(1:2,edge,el) = el
!     IF (is_wallsymm_edge(edge,el)) THEN
! ! First point on the edge (i.e. 1st vertex node on the edge)
!       k = node(el,edge)
!       ij = global_to_local_map(k,el)
!       
!       edgeM1 = edge-1
!       IF (edge.eq.1) edgeM1=4
!       edgeP1=edge+1
!       IF (edge.eq.4) edgeP1=1
! 
!       i = local_ij_to_edge_node(ij,edgeM1)
!       upwind_wallsymm_1d(1,edge,el) = upwind_local_edge_node(i,edgeM1,el)
!       
!       i = local_ij_to_edge_node(ij,edgeP1)
!       upwind_wallsymm_1d(2,edge,el) = upwind_local_edge_node(i,edgeP1,el)
!       
!       upwinded_element(k)=upwind_wallsymm_1d(1,edge,el)
!       upwinded_element(node(el,edgeP1))=upwind_wallsymm_1d(2,edge,el)
!     ENDIF
!   ENDDO
!       ENDDO
!     ENDIF



  END SUBROUTINE calc_Upwinding_local_edge_nodes  
  
  SUBROUTINE calc_upwinded_GradVec(in,outx,outy)
! Calculates the gradient of one component of a vector field
! This is done by restricting into each local element according to the upwinding scheme.
!
! E.G. if we have G(u,v), then input would be u at all GL points
!      and outputs would be:
!      outx = the differential of u w.r.t x at all GL points
!      outy = the differential of u w.r.t y at all GL points
    IMPLICIT NONE
    INTEGER :: el,i,j,ij,p,pj,ip
    DOUBLE PRECISION :: temp1,temp2
    DOUBLE PRECISION, DIMENSION(nptot) :: in,outx,outy,globalx,globaly
    DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1) :: local_in,localx_out,localy_out

  
    DO el=1,numelm
      CALL veclocalrestrict(in,el,local_in)
      DO j=0,N
        DO i=0,N
          ij=i+j*(N+1)
          temp1=0d0
          temp2=0d0
          DO p=0,N
            pj=p+j*(N+1)
            ip=i+p*(N+1)
            temp1 = temp1 + local_in(pj)*d(i,p)
            temp2 = temp2 + local_in(ip)*d(j,p) ! switched p and q to save a new loop
          ENDDO
          localx_out(ij) = (dyde(i,j,el)*temp1 - dydp(i,j,el)*temp2)/jac(i,j,el)
          localy_out(ij) = (dxdp(i,j,el)*temp2 - dxde(i,j,el)*temp1)/jac(i,j,el)
        ENDDO
      ENDDO
      globalx=0d0
      globaly=0d0
      CALL vecglobalprolongation(localx_out,el,globalx)
      CALL vecglobalprolongation(localy_out,el,globaly)
      DO i=1,nptot
        IF (upwinded_element(i).eq.el) THEN
          outx(i) = globalx(i)
          outy(i) = globaly(i)
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE calc_upwinded_GradVec

  SUBROUTINE calcGradVec(in,outx,outy)
! Calculates the gradient of one component of a vector field
! This is done by restricting into each local element.
!
! In the case of a node shared by more than 1 element, we take equal
! contributions and divide by the multiplicity of the global node.
!
! E.G. if we have G(u,v), then input would be u at all GL points
!      and outputs would be:
!      outx = the differential of in w.r.t x at all GL points
!      outy = the differential of in w.r.t y at all GL points
!
! Equivalent in Axisymmetric would be:
!      outx = diff of in w.r.t. z at each GL point
!      outy = diff of in w.r.t. r at each GL point
!
! Note: that the differential of in w.r.t. theta is not given as this routine
!       is for use in OIFS where it is not needed since u_z (equiv u_theta) = 0

    IMPLICIT NONE
    INTEGER :: el,i,j,ij,p,pj,ip
    DOUBLE PRECISION :: temp1,temp2
    DOUBLE PRECISION, DIMENSION(nptot) :: in,outx,outy,globalx,globaly
    DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1) :: local_in,localx_out,localy_out

    outx=0d0
    outy=0d0
    DO el=1,numelm
      CALL veclocalrestrict(in,el,local_in)
      DO j=0,N
        DO i=0,N
          ij=i+j*(N+1)
          temp1=0d0
          temp2=0d0
          DO p=0,N
            pj=p+j*(N+1)
            ip=i+p*(N+1)
            temp1 = temp1 + local_in(pj)*d(i,p)
            temp2 = temp2 + local_in(ip)*d(j,p) ! switched p and q to save a new loop
          ENDDO
          localx_out(ij) = (dyde(i,j,el)*temp1 - dydp(i,j,el)*temp2)/jac(i,j,el)
          localy_out(ij) = (dxdp(i,j,el)*temp2 - dxde(i,j,el)*temp1)/jac(i,j,el)
        ENDDO
      ENDDO
      globalx=0d0
      globaly=0d0
      CALL vecglobalprolongation(localx_out,el,globalx)
      CALL vecglobalprolongation(localy_out,el,globaly)
      outx = outx + globalx
      outy = outx + globaly
    ENDDO

    DO i=1,npedg
      outx(i)=outx(i)/dfloat(mult(i))
      outy(i)=outy(i)/dfloat(mult(i))
    ENDDO
  END SUBROUTINE calcGradVec
  
 
  DOUBLE PRECISION FUNCTION calcWaters_analytical_U(y,t) RESULT(out_U)

    IMPLICIT NONE
    INTEGER :: i,n_crit
    DOUBLE PRECISION, INTENT(IN) :: y,t
    DOUBLE PRECISION :: Ay,El,&
      add_on_u,temp,crit,&
      bigN,alphaN,betaN,temp_betaN,gammaN,alphaNStar,betaNstar,aN,bN,pN,qN,bigGN
    
    crit = 1d-6
    n_crit=0
    Ay = -4d0*y*(y - 1d0)
    El = We/Re
    
    i=0
    add_on_u = 0d0

    DO WHILE (n_crit.le.3)
    
      i=i+1
      bigN = (2d0*dfloat(i) - 1d0)*pi
      alphaN = 1d0 + param_beta*El*bigN**2
      temp_betaN = alphaN**2 - 4d0*El*bigN**2
      IF (temp_betaN.lt.0d0) THEN 
        betaN=sqrt(-temp_betaN)
      ELSE ! betaN^2 >= 0 !
        betaN=sqrt(temp_betaN)
      ENDIF
      gammaN = 1d0 + bigN**2*El*(param_beta - 2d0)
      alphaNStar = alphaN/2d0/El
      betaNStar = betaN/2d0/El
      aN = 1d0 + gammaN/betaN
      bN = 1d0 - gammaN/betaN
      pN = -alphaNStar + betaNStar
      qN = -alphaNStar - betaNStar

      IF (temp_betaN.ge.0d0) THEN
        bigGN = 0.5d0*( aN*exp(pN*t) + bN*exp(qN*t) )
      ELSE ! betaN^2 < 0 !
        bigGN = exp(-alphaNStar*t)*(cos(betaNStar*t)+gammaN/betaN*sin(betaNStar*t))
      ENDIF
      
      temp = (sin(bigN*y)/bigN**3)*bigGN
      add_on_u = add_on_u + temp
      
      IF (abs(temp).lt.crit) THEN
        n_crit= n_crit + 1
      ELSE
        n_crit = 0
      ENDIF
    ENDDO
    
    out_U = Ay - 32d0*add_on_u
    
  END FUNCTION calcWaters_analytical_U
  
    DOUBLE PRECISION FUNCTION calcWaters_analytical_dUdy(y,t) RESULT(out_Uxy)
  
    IMPLICIT NONE
    INTEGER :: i,n_crit
    DOUBLE PRECISION, INTENT(IN) ::y,t
    DOUBLE PRECISION :: Ay1,El,&
      add_on_uxy,temp,crit,&
      bigN,alphaN,betaN,temp_betaN,gammaN,alphaNStar,betaNstar,aN,bN,pN,qN,bigGN
    
    crit=1d-6
    n_crit=0
    
    Ay1= -8d0*y + 4d0
    El = We/Re
    
    
    add_on_uxy = 0d0
    i=0
    temp=1d0
    DO WHILE (n_crit.le.3)
      i=i+1
      bigN = (2d0*dfloat(i) - 1d0)*pi
      alphaN = 1d0 + param_beta*El*bigN**2
      temp_betaN = alphaN**2 - 4d0*bigN**2*El
      IF (temp_betaN.lt.0d0) THEN 
        betaN=sqrt(-temp_betaN)
      ELSE ! betaN^2 >= 0 !
        betaN=sqrt(temp_betaN)
      ENDIF
      gammaN = 1d0 + bigN**2*El*(param_beta - 2d0)
      alphaNStar = alphaN/2d0/El
      betaNStar = betaN/2d0/El
      aN = 1d0 + gammaN/betaN
      bN = 1d0 - gammaN/betaN
      pN = -alphaNStar + betaNStar
      qN = -alphaNStar - betaNStar

      IF (temp_betaN.ge.0d0) THEN
        bigGN = 0.5d0*( aN*exp(pN*t) + bN*exp(qN*t) )
      ELSE ! betaN^2 < 0 !
        bigGN = exp(-alphaNStar*t)*(cos(betaNStar*t)+gammaN/betaN*sin(betaNStar*t))
      ENDIF
      
      temp = (cos(bigN*y)/bigN**2)*bigGN
      
      add_on_uxy = add_on_uxy + temp 
      
      IF (abs(temp).lt.crit) THEN
        n_crit= n_crit + 1
      ELSE
        n_crit = 0
      ENDIF      
    ENDDO

    out_Uxy = Ay1 - 32d0*add_on_uxy

  END FUNCTION calcWaters_analytical_dUdy
    
  DOUBLE PRECISION FUNCTION calcWaters_analytical_Txx(y,t) RESULT(out_Txx)
    IMPLICIT NONE
    INTEGER :: i,j,n_crit1,n_crit2,n_crit3,n_crit4
    DOUBLE PRECISION, INTENT(IN) :: y,t
    DOUBLE PRECISION :: Ay1,El,Cxx,Cxy,temp1,temp2,temp3,temp4,crit,&
      add_on_txx1,add_on_txx2,add_on_txx3,add_on_txx4,add_on_KNM,&
      bigN,alphaN,betaN,temp_betaN,gammaN,alphaNStar,betaNstar,aN,bN,pN,qN,hN,cNsq,bigHN,bigIN,bigJN,&
      bigM,alphaM,betaM,temp_betaM,gammaM,alphaMStar,betaMstar,aM,bM,pM,qM,hM,cMsq,&
      bigKNM,hNM,betaNMplus,betaNMminus,cNMsqplus,cNMsqminus,ANM,BNM,DNM,ENM
    
    crit = 1d-6
    n_crit1 = 0
    n_crit2 = 0
    n_crit3 = 0
    n_crit4 = 0
    Ay1= -8d0*y + 4d0
    El = We/Re
    
! Calculate the constant from integration, Cxx and Cxy. These are constants which satisfy Txx=Txy=0 at t=0.
! The subroutine uses the same loops as below, but with all time-dependent terms enforced to t=0, eg, exp(x*t) = 1, sin(x*t) = 0, etc.
    Cxy = calc_Cxy(y)
    Cxx = calc_Cxx(Cxy,y)

    add_on_txx1 = 0d0
    add_on_txx2 = 0d0
    add_on_txx3 = 0d0
    add_on_txx4 = 0d0

    i=0

    DO WHILE (n_crit1.le.3.AND.n_crit2.le.3.AND.n_crit3.le.3)
      i=i+1

      bigN = (2d0*dfloat(i) - 1d0)*pi
      alphaN = 1d0 + param_beta*El*bigN**2
      temp_betaN = alphaN**2 - 4d0*bigN**2*El
      IF (temp_betaN.lt.0d0) THEN 
        betaN=sqrt(-temp_betaN)
      ELSE ! betaN^2 >= 0 !
        betaN=sqrt(temp_betaN)
      ENDIF
      gammaN = 1d0 + bigN**2*El*(param_beta - 2d0)
      alphaNStar = alphaN/2d0/El
      betaNStar = betaN/2d0/El
      aN = 1d0 + gammaN/betaN
      bN = 1d0 - gammaN/betaN
      pN = -alphaNStar + betaNStar
      qN = -alphaNStar - betaNStar
    
      hN = -alphaNStar + 1d0/El
      cNsq = hN**2 + betaNStar**2
    
      IF (temp_betaN.ge.0d0) THEN

        bigHN = 0.5d0*( aN/(pN + 1d0/El)*exp(pN*t) + bN/(qN + 1d0/El)*exp(qN*t) )

        bigIN = 0.5d0*( aN/pN*exp((pN - 1d0/El)*t) + bN/qN*exp((qN - 1d0/El)*t) )

        bigJN = 0.5d0*( (aN/((pN + 1d0/El)**2))*exp(pN*t) + (bN/((qN + 1d0/El)**2))*exp(qN*t) )

      ELSE ! betaN^2 < 0 !

        bigHN = exp(-alphaNStar*t)/cNsq*( (betaNStar + hN*gammaN/betaN)*sin(betaNStar*t) + &
                                       (hN - betaNStar*gammaN/betaN)*cos(betaNStar*t))
        
        bigIN = exp((-alphaNStar - 1d0/El)*t)/(alphaNStar**2 + betaNStar**2)*( &
                (betaNStar - alphaNStar*gammaN/betaN)*sin(betaNStar*t) + &
                (-alphaNStar - betaNStar*gammaN/betaN)*cos(betaNStar*t) )
            
        bigJN = exp(-alphaNStar*t)/(cNsq**2)*(&
                (hN*(betaNStar + hN*gammaN/betaN) + betaNStar*(hN - betaNStar*gammaN/betaN))*sin(betaNStar*t) &
              + (hN*(hN - betaNStar*gammaN/betaN) - betaNStar*(betaNStar + hN*gammaN/betaN))*cos(betaNStar*t) )
          
      ENDIF
      add_on_KNM=0d0
      
      j=0
      n_crit4=0
      DO WHILE (n_crit4.le.3)
        j=j+1

        bigM = (2d0*dfloat(j) - 1d0)*pi
        alphaM = 1d0 + param_beta*El*bigM**2
        temp_betaM = alphaM**2 - 4d0*bigM**2*El
        IF (temp_betaM.lt.0d0) THEN
            betaM = sqrt(-temp_betaM)
        ELSE ! temp_betaM >= 0 !
            betaM = sqrt(temp_betaM)
        ENDIF
        
        gammaM = 1d0 + bigM**2*El*(param_beta - 2d0)
        alphaMStar = alphaM/2d0/El
        betaMStar = betaM/2d0/El
        aM = 1d0 + gammaM/betaM
        bM = 1d0 - gammaM/betaM
        pM = -alphaMStar + betaMStar
        qM = -alphaMStar - betaMStar
        
        hM = -alphaMStar + 1d0/El
        cMsq = hM**2 + betaMStar**2
        
        IF (temp_betaM.ge.0d0.AND.temp_betaN.ge.0d0) THEN
            
            bigKNM = 0.25d0*( &
              aN*aM/((pM + 1d0/El)*(pN + pM + 1d0/El))*exp((pN + pM)*t) + &
              aN*bM/((qM + 1d0/El)*(pN + qM + 1d0/El))*exp((pN + qM)*t) + &
              bN*aM/((pM + 1d0/El)*(qN + pM + 1d0/El))*exp((qN + pM)*t) + &
              bN*bM/((qM + 1d0/El)*(qN + qM + 1d0/El))*exp((qN + qM)*t) )

        ELSEIF (temp_betaM.lt.0d0.AND.temp_betaN.lt.0d0) THEN

          hNM = -alphaNStar - alphaMStar + 1d0/El
          betaNMplus = betaNStar + betaMStar
          betaNMminus = betaNStar - betaMStar
          cNMsqplus = hNM**2 + betaNMplus**2
          cNMsqminus = hNM**2 + betaNMminus**2
            
          ANM = 0.5d0*(betaMStar + hM*gammaM/betaM + hM*gammaN/betaN - betaMStar*gammaN*gammaM/betaN/betaM)
          BNM = 0.5d0*(-betaMStar - hM*gammaM/betaM + hM*gammaN/betaN - betaMStar*gammaN*gammaM/betaN/betaM)
          DNM = 0.5d0*(hM - betaMStar*gammaM/betaM - betaMStar*gammaN/betaN - hM*gammaN*gammaM/betaN/betaM)
          ENM = 0.5d0*(hM - betaMStar*gammaM/betaM + betaMStar*gammaN/betaN + hM*gammaN*gammaM/betaN/betaM)
            
          bigKNM = exp(-(alphaNStar + alphaMStar)*t)/cMsq*( &
            (ANM*hNM + DNM*betaNMplus)*sin(betaNMplus*t)/cNMsqplus + &
            (DNM*hNM - ANM*betaNMplus)*cos(betaNMplus*t)/cNMsqplus + &
            (BNM*hNM + ENM*betaNMminus)*sin(betaNMminus*t)/cNMsqminus + &
            (ENM*hNM - BNM*betaNMminus)*cos(betaNMminus*t)/cNMsqminus )
        ELSE
          bigKNM = 0d0
        ENDIF

        temp4 = cos(bigN*y)*cos(bigM*y)*bigKNM/bigM**2/bigN**2
        add_on_txx4 = add_on_txx4 + temp4
        IF (abs(temp4).lt.crit) THEN
          n_crit4 = n_crit4 + 1
        ELSE
          n_crit4 = 0
        ENDIF
      ENDDO

      temp1 = cos(bigN*y)*bigIN/bigN**2
      temp2 = cos(bigN*y)*bigHN/bigN**2
      temp3 = cos(bigN*y)*bigJN/bigN**2

      IF (n_crit1.le.3) THEN
        add_on_txx1 = add_on_txx1 + temp1
        IF (abs(temp1).lt.crit) THEN
          n_crit1 = n_crit1 + 1
        ELSE 
          n_crit1 = 0
        ENDIF
      ENDIF
      IF (n_crit2.le.3) THEN
        add_on_txx2 = add_on_txx2 + temp2
        IF (abs(temp2).lt.crit) THEN
          n_crit2 = n_crit2 + 1
        ELSE 
          n_crit2 = 0
        ENDIF
      ENDIF
      IF (n_crit3.le.3) THEN
        add_on_txx3 = add_on_txx3 + temp3
        IF (abs(temp3).lt.crit) THEN
          n_crit3 = n_crit3 + 1
        ELSE 
          n_crit3 = 0
        ENDIF
      ENDIF
     
      
    ENDDO
    
    out_Txx = 2d0*Re*Cxy*( Ay1*exp(-t/El)*t - 32d0*add_on_txx1) + &
      2d0*Re*(1d0 - param_beta)*Ay1*(Ay1*El - 32d0*add_on_txx2) - &
      64d0*Re*Ay1*(1d0 - param_beta)*add_on_txx3/El + &
      2048d0*Re*(1d0 - param_beta)*add_on_txx4/El + &
      Cxx*exp(-t/El)

  END FUNCTION calcWaters_analytical_Txx
    
  DOUBLE PRECISION FUNCTION calcWaters_analytical_Txy(y,t) RESULT(out_Txy)
    IMPLICIT NONE
    INTEGER :: i,n_crit
    DOUBLE PRECISION, INTENT(IN) :: y,t
    DOUBLE PRECISION :: Ay1,El,Cxy,&
      add_on_txy,temp,crit,&
      bigN,alphaN,betaN,temp_betaN,gammaN,alphaNStar,betaNstar,aN,bN,pN,qN,hN,cNsq,bigHN

    crit = 1d-6
    n_crit = 0

    Ay1= -8d0*y + 4d0

    El = We/Re

! Calculate the constant from integration, Cxx and Cxy. These are constants which satisfy Txx=Txy=0 at t=0.
! The subroutine uses the same loops as below, but with all time-dependent terms enforced to t=0, eg, exp(x*t) = 1, sin(x*t) = 0, etc.
    Cxy = calc_Cxy(y)
    
    add_on_txy = 0d0
    i=0
    temp=1d0
    DO WHILE (n_crit.le.3)
      i=i+1
      bigN = (2d0*dfloat(i) - 1d0)*pi
      alphaN = 1d0 + param_beta*El*bigN**2
      temp_betaN = alphaN**2 - 4d0*bigN**2*El
      
      IF (temp_betaN.lt.0d0) THEN 
        betaN=sqrt(-temp_betaN)
      ELSE ! betaN^2 >= 0 !
        betaN=sqrt(temp_betaN)
      ENDIF
      
      gammaN = 1d0 + bigN**2*El*(param_beta - 2d0)
      alphaNStar = alphaN/2d0/El
      betaNStar = betaN/2d0/El
      aN = 1d0 + gammaN/betaN
      bN = 1d0 - gammaN/betaN
      pN = -alphaNStar + betaNStar
      qN = -alphaNStar - betaNStar
      hN = -alphaNStar + 1d0/El
      cNsq = hN**2 + betaNStar**2
    
      IF (temp_betaN.ge.0d0) THEN

        bigHN = 0.5d0*( aN/(pN + 1d0/El)*exp(pN*t) + bN/(qN + 1d0/El)*exp(qN*t) )

      ELSE ! betaN^2 < 0 !

        bigHN = exp(-alphaNStar*t)/cNsq*( (betaNStar + hN*gammaN/betaN)*sin(betaNStar*t) + &
          (hN - betaNStar*gammaN/betaN)*cos(betaNStar*t) )

      ENDIF
      temp = cos(bigN*y)/bigN**2*bigHN
      add_on_txy = add_on_txy + temp
      
      IF (abs(temp).lt.crit) THEN
        n_crit = n_crit + 1
      ELSE
        n_crit = 0
      ENDIF
    ENDDO
    
    out_Txy = (1d0 - param_beta)/El*(El*Ay1 - 32d0*add_on_txy) + Cxy*exp(-t/El)

  END FUNCTION calcWaters_analytical_Txy


  DOUBLE PRECISION FUNCTION calc_Cxx(Cxy,y) RESULT(out_Cxx)
  
    IMPLICIT NONE
    INTEGER :: i,j,n_crit1,n_crit2,n_crit3,n_crit4
    DOUBLE PRECISION, INTENT(IN) :: y,Cxy
    DOUBLE PRECISION :: Ay1,El,t,temp1,temp2,temp3,temp4,crit,&
      add_on_txx1,add_on_txx2,add_on_txx3,add_on_txx4,add_on_KNM,&
      bigN,alphaN,betaN,temp_betaN,gammaN,alphaNStar,betaNstar,aN,bN,pN,qN,hN,cNsq,bigHN,bigIN,bigJN,&
      bigM,alphaM,betaM,temp_betaM,gammaM,alphaMStar,betaMstar,aM,bM,pM,qM,hM,cMsq,&
      bigKNM,hNM,betaNMplus,betaNMminus,cNMsqplus,cNMsqminus,ANM,BNM,DNM,ENM

    t=0d0
    crit = 1d-6
    n_crit1 = 0
    n_crit2 = 0
    n_crit3 = 0
    n_crit4 = 0
    Ay1= -8d0*y + 4d0
    El = We/Re

    add_on_txx1 = 0d0
    add_on_txx2 = 0d0
    add_on_txx3 = 0d0
    add_on_txx4 = 0d0

    i=0

    DO WHILE (n_crit1.le.3.AND.n_crit2.le.3.AND.n_crit3.le.3)
      i=i+1

      bigN = (2d0*dfloat(i) - 1d0)*pi
      alphaN = 1d0 + param_beta*El*bigN**2
      temp_betaN = alphaN**2 - 4d0*bigN**2*El
      IF (temp_betaN.lt.0d0) THEN 
        betaN=sqrt(-temp_betaN)
      ELSE ! betaN^2 >= 0 !
        betaN=sqrt(temp_betaN)
      ENDIF
      gammaN = 1d0 + bigN**2*El*(param_beta - 2d0)
      alphaNStar = alphaN/2d0/El
      betaNStar = betaN/2d0/El
      aN = 1d0 + gammaN/betaN
      bN = 1d0 - gammaN/betaN
      pN = -alphaNStar + betaNStar
      qN = -alphaNStar - betaNStar
    
      hN = -alphaNStar + 1d0/El
      cNsq = hN**2 + betaNStar**2
    
      IF (temp_betaN.ge.0d0) THEN

        bigHN = 0.5d0*( aN/(pN + 1d0/El)*exp(pN*t) + bN/(qN + 1d0/El)*exp(qN*t) )

        bigIN = 0.5d0*( aN/pN*exp((pN - 1d0/El)*t) + bN/qN*exp((qN - 1d0/El)*t) )

        bigJN = 0.5d0*( (aN/((pN + 1d0/El)**2))*exp(pN*t) + (bN/((qN + 1d0/El)**2))*exp(qN*t) )

      ELSE ! betaN^2 < 0 !
      
        bigHN = exp(-alphaNStar*t)/cNsq*( (betaNStar + hN*gammaN/betaN)*sin(betaNStar*t) + &
          (hN - betaNStar*gammaN/betaN)*cos(betaNStar*t))

        bigIN = exp((-alphaNStar - 1d0/El)*t)/(alphaNStar**2 + betaNStar**2)*( &
          (betaNStar - alphaNStar*gammaN/betaN)*sin(betaNStar*t) + &
          (-alphaNStar - betaNStar*gammaN/betaN)*cos(betaNStar*t) )

        bigJN = exp(-alphaNStar*t)/(cNsq**2)*(&
          (hN*(betaNStar + hN*gammaN/betaN) + betaNStar*(hN - betaNStar*gammaN/betaN))*sin(betaNStar*t) &
          + (hN*(hN - betaNStar*gammaN/betaN) - betaNStar*(betaNStar + hN*gammaN/betaN))*cos(betaNStar*t) )

      ENDIF
      add_on_KNM=0d0

      j=0
      n_crit4=0
      DO WHILE (n_crit4.le.3)
        j=j+1

        bigM = (2d0*dfloat(j) - 1d0)*pi
        alphaM = 1d0 + param_beta*El*bigM**2
        temp_betaM = alphaM**2 - 4d0*bigM**2*El
        IF (temp_betaM.lt.0d0) THEN
            betaM = sqrt(-temp_betaM)
        ELSE ! temp_betaM >= 0 !
            betaM = sqrt(temp_betaM)
        ENDIF
        
        gammaM = 1d0 + bigM**2*El*(param_beta - 2d0)
        alphaMStar = alphaM/2d0/El
        betaMStar = betaM/2d0/El
        aM = 1d0 + gammaM/betaM
        bM = 1d0 - gammaM/betaM
        pM = -alphaMStar + betaMStar
        qM = -alphaMStar - betaMStar

        hM = -alphaMStar + 1d0/El
        cMsq = hM**2 + betaMStar**2

        IF (temp_betaM.ge.0d0.AND.temp_betaN.ge.0d0) THEN

          bigKNM = 0.25d0*( &
            aN*aM/((pM + 1d0/El)*(pN + pM + 1d0/El))*exp((pN + pM)*t) + &
            aN*bM/((qM + 1d0/El)*(pN + qM + 1d0/El))*exp((pN + qM)*t) + &
            bN*aM/((pM + 1d0/El)*(qN + pM + 1d0/El))*exp((qN + pM)*t) + &
            bN*bM/((qM + 1d0/El)*(qN + qM + 1d0/El))*exp((qN + qM)*t) )

        ELSEIF (temp_betaM.lt.0d0.AND.temp_betaN.lt.0d0) THEN

          hNM = -alphaNStar - alphaMStar + 1d0/El
          betaNMplus = betaNStar + betaMStar
          betaNMminus = betaNStar - betaMStar
          cNMsqplus = hNM**2 + betaNMplus**2
          cNMsqminus = hNM**2 + betaNMminus**2

          ANM = 0.5d0*(betaMStar + hM*gammaM/betaM + hM*gammaN/betaN - betaMStar*gammaN*gammaM/betaN/betaM)
          BNM = 0.5d0*(-betaMStar - hM*gammaM/betaM + hM*gammaN/betaN - betaMStar*gammaN*gammaM/betaN/betaM)
          DNM = 0.5d0*(hM - betaMStar*gammaM/betaM - betaMStar*gammaN/betaN - hM*gammaN*gammaM/betaN/betaM)
          ENM = 0.5d0*(hM - betaMStar*gammaM/betaM + betaMStar*gammaN/betaN + hM*gammaN*gammaM/betaN/betaM)
            
          bigKNM = exp(-(alphaNStar + alphaMStar)*t)/cMsq*( &
            (ANM*hNM + DNM*betaNMplus)*sin(betaNMplus*t)/cNMsqplus + &
            (DNM*hNM - ANM*betaNMplus)*cos(betaNMplus*t)/cNMsqplus + &
            (BNM*hNM + ENM*betaNMminus)*sin(betaNMminus*t)/cNMsqminus + &
            (ENM*hNM - BNM*betaNMminus)*cos(betaNMminus*t)/cNMsqminus )
        ELSE
          bigKNM = 0d0
        ENDIF
  
        temp4 = cos(bigN*y)*cos(bigM*y)*bigKNM/bigM**2/bigN**2
        add_on_txx4 = add_on_txx4 + temp4
        IF (abs(temp4).lt.crit) THEN
          n_crit4 = n_crit4 + 1
        ELSE
          n_crit4 = 0
        ENDIF
      ENDDO

      temp1 = cos(bigN*y)*bigIN/bigN**2
      temp2 = cos(bigN*y)*bigHN/bigN**2
      temp3 = cos(bigN*y)*bigJN/bigN**2
      
      IF (n_crit1.le.3) THEN
        add_on_txx1 = add_on_txx1 + temp1
        IF (abs(temp1).lt.crit) THEN
          n_crit1 = n_crit1 + 1
        ELSE 
          n_crit1 = 0
        ENDIF
      ENDIF
      IF (n_crit2.le.3) THEN
        add_on_txx2 = add_on_txx2 + temp2
        IF (abs(temp2).lt.crit) THEN
          n_crit2 = n_crit2 + 1
        ELSE 
          n_crit2 = 0
        ENDIF
      ENDIF
      IF (n_crit3.le.3) THEN
        add_on_txx3 = add_on_txx3 + temp3
        IF (abs(temp3).lt.crit) THEN
          n_crit3 = n_crit3 + 1
        ELSE 
          n_crit3 = 0
        ENDIF
      ENDIF
    ENDDO

    out_Cxx =-( 2d0*Re*Cxy*( Ay1*exp(-t/El)*t - 32d0*add_on_txx1) + &
      2d0*Re*(1d0 - param_beta)*Ay1*(Ay1*El - 32d0*add_on_txx2) - &
      64d0*Re*Ay1*(1d0 - param_beta)/El*add_on_txx3 + &
      2048d0*Re*(1d0 - param_beta)/El*add_on_txx4 )

!     IMPLICIT NONE
!     INTEGER, INTENT(IN) :: nterms
!     INTEGER :: i,j
!     DOUBLE PRECISION, INTENT(IN) :: Cxy,y
!     DOUBLE PRECISION :: Ay1,El,&
!       add_on_cxx1,add_on_cxx2,add_on_cxx3,add_on_cxx4,add_on_KNM,&
!       bigN,alphaN,betaN,temp_betaN,gammaN,alphaNStar,betaNstar,aN,bN,pN,qN,hN,cNsq,bigHN,bigIN,bigJN,&
!       bigM,alphaM,betaM,temp_betaM,gammaM,alphaMStar,betaMstar,aM,bM,pM,qM,hM,cMsq,&
!       bigKNM,hNM,betaNMplus,betaNMminus,cNMsqplus,cNMsqminus,ANM,BNM,DNM,ENM
!     
!     Ay1 = -8d0*y + 4d0
!     El = We/Re
!     
!     add_on_cxx1 = 0d0
!     add_on_cxx2 = 0d0
!     add_on_cxx3 = 0d0
!     add_on_cxx4 = 0d0
! 
!     DO i=1,nterms
!     
!       bigN = (2d0*i - 1d0)*pi
!       alphaN = 1d0 + param_beta*El*bigN**2
!       temp_betaN = alphaN**2 - 4d0*bigN**2*El
!       
!       IF (temp_betaN.lt.0d0) THEN 
!       
!   betaN=sqrt(-temp_betaN)
!   
!       ELSE ! betaNsq >= 0 !
!       
!   betaN=sqrt(temp_betaN)
!   
!       ENDIF
!       
!       gammaN = 1d0 + bigN**2*El*(param_beta - 2d0)
!       alphaNStar = alphaN/2d0/El
!       betaNStar = betaN/2d0/El
!       aN = 1d0 + gammaN/betaN
!       bN = 1d0 - gammaN/betaN
!       pN = -alphaNStar + betaNStar
!       qN = -alphaNStar - betaNStar
!       hN = -alphaNStar + 1d0/El
!       cNsq = hN**2 + betaNStar**2
!     
!       IF (temp_betaN.ge.0d0) THEN
!       
!         bigHN = 0.5d0*( aN/(pN + 1d0/El) + bN/(qN + 1d0/El) )
!         
!         bigIN = 0.5d0*( aN/pN + bN/qN )
!         
!         bigJN = 0.5d0*( aN/(pN + 1d0/El)**2 + bN/(qN + 1d0/El)**2 )
!         
!        
!       ELSE ! betaN^2 < 0 !
!       
!         bigHN = 1d0/cNsq*( hN - betaNStar*gammaN/betaN )
!         
!         bigIN = 1d0/(alphaNStar**2 + betaNStar**2)*(-alphaNStar - betaNStar*gammaN/betaN)
!             
!         bigJN = 1d0/cNsq**2*( (hN*(hN - betaNStar*gammaN/betaN) - betaNStar*(betaNStar + hN*gammaN/betaN)) )
!           
!       ENDIF
! 
!       add_on_KNM = 0d0
!       DO j=1,nterms
! 
!   bigM = (2d0*j - 1d0)*pi
!         alphaM = 1d0 + param_beta*El*bigM**2
!         temp_betaM = alphaM**2 - 4d0*bigM**2*El
!         
!         IF (temp_betaM.lt.0d0) THEN
!         
!             betaM = sqrt(-temp_betaM)
!             
!         ELSE ! betaMsq >= 0 !
!         
!             betaM = sqrt(temp_betaM)
!             
!         ENDIF
!         
!         gammaM = 1d0 + bigM**2*El*(param_beta - 2d0)
!         alphaMStar = alphaM/2d0/El
!         betaMStar = betaM/2d0/El
!         aM = 1d0 + gammaM/betaM
!         bM = 1d0 - gammaM/betaM
!         pM = -alphaMStar + betaMStar
!         qM = -alphaMStar - betaMStar
!         hM = -alphaMStar + 1d0/El
!         cMsq = hM**2 + betaMStar**2
! 
!         IF (temp_betaM.ge.0d0.AND.temp_betaN.ge.0d0) THEN
!             
!             bigKNM = 0.25d0*( &
!         aN*aM/((pM + 1d0/El)*(pN + pM + 1d0/El)) + &
!         aN*bM/((qM + 1d0/El)*(pN + qM + 1d0/El)) + &
!         bN*aM/((pM + 1d0/El)*(qN + pM + 1d0/El)) + &
!         bN*bM/((qM + 1d0/El)*(qN + qM + 1d0/El))   &
!        )
!             
!         ELSEIF (temp_betaM.lt.0d0) THEN
!             
!     hNM = -alphaNStar - alphaMStar + 1d0/El
!           betaNMplus = betaNStar + betaMStar
!           betaNMminus = betaNStar - betaMStar
!           cNMsqplus = hNM**2 + betaNMplus**2
!           cNMsqminus = hNM**2 + betaNMminus**2
!             
!           ANM = 0.5d0*(betaMStar + hM*gammaM/betaM + hM*gammaN/betaN - betaMStar*gammaN*gammaM/betaN/betaM)
!           BNM = 0.5d0*(-betaMStar - hM*gammaM/betaM + hM*gammaN/betaN - betaMStar*gammaN*gammaM/betaN/betaM)
!           DNM = 0.5d0*(hM - betaMStar*gammaM/betaM - betaMStar*gammaN/betaN - hM*gammaN*gammaM/betaN/betaM)
!           ENM = 0.5d0*(hM - betaMStar*gammaM/betaM + betaMStar*gammaN/betaN + hM*gammaN*gammaM/betaN/betaM)
!             
!           bigKNM = 1d0/cMsq*( &
!         (DNM*hNM - ANM*betaNMplus)/cNMsqplus + &
!         (ENM*hNM - BNM*betaNMminus)/cNMsqminus &
!        )
!        
!         ELSE
!         
!           bigKNM = 0d0
!           
!         ENDIF
! 
!         add_on_KNM = add_on_KNM + cos(bigN*y)*cos(bigM*y)/bigN**2/bigM**2*bigKNM
!         
!       ENDDO
!     
!    
!       add_on_cxx1 = add_on_cxx1 + cos(bigN*y)/bigN**2*bigIN
!       add_on_cxx2 = add_on_cxx2 + cos(bigN*y)/bigN**2*bigHN
!       add_on_cxx3 = add_on_cxx3 + cos(bigN*y)/bigN**2*bigJN
!       add_on_cxx4 = add_on_cxx4 + add_on_KNM
!       
!     ENDDO
!        
!     out_Cxx = - (2d0*Re*Cxy*(-32d0*add_on_cxx1) + &
!      2d0*Re*(1d0 - param_beta)*Ay1*(Ay1*El - 32d0*add_on_cxx2) - &
!      64d0*Re*Ay1*(1d0 - param_beta)/El*add_on_cxx3 + &
!      2048d0*Re*(1d0 - param_beta)/El*add_on_cxx4 )
  END FUNCTION calc_Cxx

  DOUBLE PRECISION FUNCTION calc_Cxy(y) RESULT(out_Cxy)

    IMPLICIT NONE
    INTEGER :: i,n_crit
    DOUBLE PRECISION, INTENT(IN) :: y
    DOUBLE PRECISION :: Ay1,El,t,&
      add_on_txy,temp,crit,&
      bigN,alphaN,betaN,temp_betaN,gammaN,alphaNStar,betaNstar,aN,bN,pN,qN,hN,cNsq,bigHN

    crit = 1d-6
    n_crit = 0
    t=0d0

    Ay1= -8d0*y + 4d0

    El = We/Re

    add_on_txy = 0d0
    i=0
    temp=1d0
    DO WHILE (n_crit.le.3)
      i=i+1
      bigN = (2d0*dfloat(i) - 1d0)*pi
      alphaN = 1d0 + param_beta*El*bigN**2
      temp_betaN = alphaN**2 - 4d0*bigN**2*El
      
      IF (temp_betaN.lt.0d0) THEN
        betaN=sqrt(-temp_betaN)
      ELSE ! betaN^2 >= 0 !      
        betaN=sqrt(temp_betaN)
      ENDIF
      
      gammaN = 1d0 + bigN**2*El*(param_beta - 2d0)
      alphaNStar = alphaN/2d0/El
      betaNStar = betaN/2d0/El
      aN = 1d0 + gammaN/betaN
      bN = 1d0 - gammaN/betaN
      pN = -alphaNStar + betaNStar
      qN = -alphaNStar - betaNStar
      hN = -alphaNStar + 1d0/El
      cNsq = hN**2 + betaNStar**2
    
      IF (temp_betaN.ge.0d0) THEN
      
        bigHN = 0.5d0*( aN/(pN + 1d0/El)*exp(pN*t) + bN/(qN + 1d0/El)*exp(qN*t) )

      ELSE ! betaN^2 < 0 !

        bigHN = exp(-alphaNStar*t)/cNsq*( (betaNStar + hN*gammaN/betaN)*sin(betaNStar*t) + &
          (hN - betaNStar*gammaN/betaN)*cos(betaNStar*t))

      ENDIF
      temp = cos(bigN*y)/bigN**2*bigHN
      add_on_txy = add_on_txy + temp
      
      IF (abs(temp).lt.crit) THEN
        n_crit = n_crit + 1
      ELSE
        n_crit = 0
      ENDIF
    ENDDO

    out_Cxy = (1d0 - param_beta)/El*(El*Ay1 - 32d0*add_on_txy)

!       IMPLICIT NONE
!       INTEGER :: i,j
!       INTEGER, INTENT(IN) :: nterms
!       DOUBLE PRECISION, INTENT(IN) :: y
!       DOUBLE PRECISION :: Ay1,El,&
!         add_on_cxy,&
!         bigN,alphaN,betaN,temp_betaN,gammaN,alphaNStar,betaNstar,aN,bN,pN,qN,hN,cNsq,bigHN
!       
!       Ay1 = -8d0*y + 4d0
!       El = We/Re
! 
!       add_on_cxy = 0d0
! 
!       DO i=1,nterms
!     
!   bigN = (2d0*i - 1d0)*pi
!   alphaN = 1d0 + param_beta*El*bigN**2
!   temp_betaN = alphaN**2 - 4d0*bigN**2*El
!   
!   IF (temp_betaN.lt.0d0) THEN 
!   
!     betaN=sqrt(-temp_betaN)
!     
!   ELSE ! betaN^2 >= 0 !
!   
!     betaN=sqrt(temp_betaN)
!     
!   ENDIF
!   
!   gammaN = 1d0 + bigN**2*El*(param_beta - 2d0)
!   alphaNStar = alphaN/2d0/El
!   betaNStar = betaN/2d0/El
!   aN = 1d0 + gammaN/betaN
!   bN = 1d0 - gammaN/betaN
!   pN = -alphaNStar + betaNStar
!   qN = -alphaNStar - betaNStar
!   hN = -alphaNStar + 1d0/El
!   cNsq = hN**2 + betaNStar**2
!     
!   IF (temp_betaN.ge.0d0) THEN
!   
!     bigHN = 0.5d0*( aN/(pN + 1d0/El) + bN/(qN + 1d0/El) )
!     
!   ELSE ! betaN^2 < 0 !
!   
!     bigHN = 1d0/cNsq*( hN - betaNStar*gammaN/betaN )
!     
!   ENDIF
!   
!   add_on_cxy = add_on_cxy + cos(bigN*y)/bigN**2*bigHN
!   
!       ENDDO
!     
!       out_Cxy = -(1d0 - param_beta)/El*(El*Ay1 - 32d0*add_on_cxy)
  END FUNCTION calc_Cxy


  SUBROUTINE map_x_y_to_xi_eta(x_in, y_in, el_out, xi_out, eta_out)
     IMPLICIT NONE    
     INTEGER :: el_out,el,i    
     DOUBLE PRECISION  :: xi_out, eta_out, x_in, y_in, A, cornerx(4), cornery(4)


!! Finds the element inwhich a point (x,y) resides
    el_out=-1
    elloop: DO el=1,numelm
      DO i=1,4
        A = (vertx(node(el,i+1),1) - vertx(node(el,i),1))*(y_in - vertx(node(el,i),2)) &
          - (x_in - vertx(node(el,i),1))*(vertx(node(el,i+1),2) - vertx(node(el,i),2))
!! Determines if point is on left (A>0) of a line (see Computational Geometry in C. J. Roukre)

        IF (A.lt.-(1d-14)) EXIT
!! If point is left of all sides of element, then
!! point is within that element
        IF (i.eq.4) el_out=el

      ENDDO
    ENDDO elloop

    IF (el_out.lt.0) THEN
      write(*,*) 'could not find point within elements',x_in,y_in
      STOP
    ENDIF
    
!! Found the correct element now find xi & eta
    DO i=1,4
      cornerx(i) = vertx(node(el_out,i),1)
      cornery(i) = vertx(node(el_out,i),2)
    ENDDO

!! Inverts the transfinite map to find (xi, eta), given (x, y)
    CALL quad2squ(x_in, y_in, cornerx, cornery, xi_out, eta_out)

!! If cant find a (xi,eta) - error
    IF (xi_out.lt.-(1d0).or.xi_out.gt.(1d0)) THEN
      write(*,*)'error xi=',xi_out,'el=',el_out
      STOP
    ENDIF
    IF (eta_out.lt.-(1d0).or.eta_out.gt.(1d0)) THEN
      write(*,*)'error eta=',eta_out,'el=',el_out
      STOP
    ENDIF

    RETURN
  END SUBROUTINE map_x_y_to_xi_eta


  SUBROUTINE quad2squ(x_in, y_in, cornerx, cornery, xi_out, eta_out)
    IMPLICIT NONE

    DOUBLE PRECISION :: a,b,c,d,e,f,g,h,aa,bb,cc,dd,ee,ff,&
      x_pos, x_neg, xi_out, eta_out, cornerx(4), cornery(4),&
      x_in, y_in, y_pos, y_neg, eps

    eps=1d-8

    a = (cornerx(1) + cornerx(2) + cornerx(3) + cornerx(4))/4d0
    b = (-cornerx(1) + cornerx(2) + cornerx(3) -cornerx(4))/4d0
    c = (-cornerx(1) - cornerx(2) + cornerx(3) + cornerx(4))/4d0
    d = (cornerx(1) - cornerx(2) + cornerx(3) - cornerx(4))/4d0

    e = (cornery(1) + cornery(2) + cornery(3) + cornery(4))/4d0
    f = (-cornery(1) + cornery(2) + cornery(3) - cornery(4))/4d0
    g = (-cornery(1) - cornery(2) + cornery(3) + cornery(4))/4d0
    h = (cornery(1) - cornery(2) + cornery(3) - cornery(4))/4d0

    aa = -f*d + h*b
    bb = (y_in - e)*d + (a - x_in)*h - c*f + g*b 
    cc = (y_in - e)*c + (a - x_in)*g

    dd = h*c - g*d
    ee = (y_in - e)*d + (a - x_in)*h + c*f - g*b
    ff = (y_in - e)*b + (a - x_in)*f

!! solves quadratic equation seen in the inverse map using
!! simple formula
    CALL quadratic_eqn_solver(aa, bb, cc, x_pos, x_neg, eps)
    CALL quadratic_eqn_solver(dd, ee, ff, y_pos, y_neg, eps)

!! write(111,*) al, be,'x',x_in,'y',y_in, x_pos, x_neg, y_pos, y_neg
!! write(111,*) 'a',aa,'b',bb,'c',cc,'d',dd,'e',ee,'f',ff

!! Of the solutions, determine which "make sense"
!! (xi,eta) should both both be between [-1,1]

    CALL root_decider(x_pos, x_neg, xi_out, eps)
    CALL root_decider(y_pos, y_neg, eta_out, eps)

!! This does always work as the map IS bijective over the 
!! element domain.
    IF (abs(xi_out-1d0).lt.eps) THEN
      xi_out=1d0
    ELSEIF (abs(xi_out+1d0).lt.eps) THEN
      xi_out=-1d0
    ENDIF
    IF (abs(eta_out-1d0).lt.eps) THEN
      eta_out=1d0
    ELSEIF (abs(eta_out+1d0).lt.eps) THEN
      eta_out=-1d0
    ENDIF      

    RETURN
  END SUBROUTINE quad2squ

   SUBROUTINE quadratic_eqn_solver(a1, b1, c1, s_pos, s_neg, eps)
    IMPLICIT NONE

    DOUBLE PRECISION :: a1, b1, c1, s_pos, s_neg, determinant, eps, abs_a1

    determinant = b1**2 - 4d0*a1*c1
    abs_a1 = abs(a1)

    IF (abs_a1.le.eps) THEN 
      s_pos = -c1/b1
      s_neg = 2d0
    ELSEIF (determinant.lt.0d0) THEN
      write(*,*) 'xi, eta, complex roots'
      STOP
    ELSE
      s_pos = (-b1 + dsqrt(determinant))/(2d0*a1)
      s_neg = (-b1 - dsqrt(determinant))/(2d0*a1)

    ENDIF
    RETURN
  END SUBROUTINE quadratic_eqn_solver


  SUBROUTINE root_decider(s_pos, s_neg, soln, eps)
    IMPLICIT NONE

    DOUBLE PRECISION :: s_pos, s_neg, soln, eps

    IF ((((s_pos.gt.(1d0+eps))).or.(s_pos.lt.(-(1d0+eps)))).and.(((s_neg.gt.(1d0+eps))).or.(s_neg.lt.(-(1d0+eps))))) THEN
      write(*,*) 'xi, zeta out of limits',s_pos,s_neg
      STOP

    ELSEIF (((s_neg.gt.(1d0+eps))).or.(s_neg.lt.(-(1d0+eps)))) THEN
      soln = s_pos
    ELSEIF (((s_pos.gt.(1d0+eps))).or.(s_pos.lt.(-(1d0+eps)))) THEN
      soln = s_neg
    ELSE
      write(*,*) 'multiple xi, eta within parent element'
      write(*,*) s_pos,s_neg
      STOP
    ENDIF
  END SUBROUTINE root_decider
    
! BEGIN FUNCTIONS FOR STOKES MODEL SOLUTION
! 
! We only require non-zero functions. They are:
! vel_x/y
! RHS f_x/y
! pressure
! gradients
! 
  DOUBLE PRECISION FUNCTION model_soln_vel_x(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    
    out = sin(PI*nodeCoord(glob_i_in,1))*cos(PI*nodeCoord(glob_i_in,2))
    
  END FUNCTION model_soln_vel_x
  
  DOUBLE PRECISION FUNCTION model_soln_vel_y(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    
    out = -cos(PI*nodeCoord(glob_i_in,1))*sin(PI*nodeCoord(glob_i_in,2))
    
  END FUNCTION model_soln_vel_y
  
  DOUBLE PRECISION FUNCTION model_soln_f_x(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y
    
    x=PI*nodeCoord(glob_i_in,1)
    y=PI*nodeCoord(glob_i_in,2)
    
    out = 2d0*sin(x)*cos(y)*PI*PI + PI*cos(x)*sin(y)
    
  END FUNCTION model_soln_f_x
  
  DOUBLE PRECISION FUNCTION model_soln_f_y(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y
    
    x=PI*nodeCoord(glob_i_in,1)
    y=PI*nodeCoord(glob_i_in,2)
    
    out = -2d0*cos(x)*sin(y)*PI*PI + PI*sin(x)*cos(y)
    
  END FUNCTION model_soln_f_y
  
  DOUBLE PRECISION FUNCTION model_soln_pressure(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    
    out = sin(PI*nodeCoord(glob_i_in,1))*sin(PI*nodeCoord(glob_i_in,2))
    
  END FUNCTION model_soln_pressure
  
!  GRADIENTS

  DOUBLE PRECISION FUNCTION model_soln_Grad_vel_xx(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    
    out = PI*cos(PI*nodeCoord(glob_i_in,1))*cos(PI*nodeCoord(glob_i_in,2))
    
  END FUNCTION model_soln_Grad_vel_xx
  
  DOUBLE PRECISION FUNCTION model_soln_Grad_vel_yx(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    
    out = -PI*sin(PI*nodeCoord(glob_i_in,1))*sin(PI*nodeCoord(glob_i_in,2))
    
  END FUNCTION model_soln_Grad_vel_yx
 
  DOUBLE PRECISION FUNCTION model_soln_Grad_vel_xy(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    
    out = PI*sin(PI*nodeCoord(glob_i_in,1))*sin(PI*nodeCoord(glob_i_in,2))
    
  END FUNCTION model_soln_Grad_vel_xy
  
  DOUBLE PRECISION FUNCTION model_soln_Grad_vel_yy(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    
    out = -PI*cos(PI*nodeCoord(glob_i_in,1))*cos(PI*nodeCoord(glob_i_in,2))
    
  END FUNCTION model_soln_Grad_vel_yy
  
! Model Solution 2 for Tim:
! u = -sin(2*PI*y)*(cos(2*PI*x) - 1)
! v = sin(2*PI*x)*(cos(2*PI*y) - 1)
! p = sin(PI*x*y)
!
! f_x = PI*y*cos(PI*x*y) - 4*PI^2*sin(2*PI*y)*(2*cos(2*PI*x) - 1) ! NOTE -/+ on the second term compared to Tim's PDF.
! f_y = PI*x*cos(PI*x*y) + 4*PI^2*sin(2*PI*x)*(2*cos(2*PI*y) - 1)

  DOUBLE PRECISION FUNCTION model_soln2_vel_x(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
  
    out = -sin(2d0*PI*y)*(cos(2d0*PI*x) - 1d0)
    
  END FUNCTION model_soln2_vel_x
  
  DOUBLE PRECISION FUNCTION model_soln2_vel_y(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
    
    out = sin(2d0*PI*x)*(cos(2d0*PI*y) - 1d0)
    
  END FUNCTION model_soln2_vel_y
  
  DOUBLE PRECISION FUNCTION model_soln2_f_x(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
    
    out = PI*y*cos(PI*x*y) - 4d0*PI*PI*sin(2d0*PI*y)*(2d0*cos(2d0*PI*x) - 1d0)
    
  END FUNCTION model_soln2_f_x
  
  DOUBLE PRECISION FUNCTION model_soln2_f_y(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
    
    out = PI*x*cos(PI*x*y) + 4d0*PI*PI*sin(2d0*PI*x)*(2d0*cos(2d0*PI*y) - 1d0)
    
  END FUNCTION model_soln2_f_y
  
  DOUBLE PRECISION FUNCTION model_soln2_pressure(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)  
    
    out = sin(PI*x*y)
    
  END FUNCTION model_soln2_pressure
  
!  GRADIENTS

  DOUBLE PRECISION FUNCTION model_soln2_Grad_vel_xx(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
  
    out = 2d0*PI*sin(2d0*PI*x)*sin(2d0*PI*y)
    
  END FUNCTION model_soln2_Grad_vel_xx
  
  DOUBLE PRECISION FUNCTION model_soln2_Grad_vel_yx(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
  
    out = -2d0*PI*cos(2d0*PI*y)*(cos(2d0*PI*x) - 1d0)
    
  END FUNCTION model_soln2_Grad_vel_yx
 
  DOUBLE PRECISION FUNCTION model_soln2_Grad_vel_xy(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
    
    out = 2d0*PI*cos(2d0*PI*x)*(cos(2d0*PI*y) - 1d0)
    
  END FUNCTION model_soln2_Grad_vel_xy
  
  DOUBLE PRECISION FUNCTION model_soln2_Grad_vel_yy(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
    
    out = -2d0*PI*sin(2d0*PI*x)*sin(2d0*PI*y)
    
  END FUNCTION model_soln2_Grad_vel_yy
  
! END FUNCTIONS FOR STOKES MODEL SOLUTION

! BEGIN FUNCTIONS FOR CARTESIAN FLOW PAST A CYLINDER
! 
! We only require non-zero functions. They are:
! vel_x/y
! pressure
! gradients
! 
  DOUBLE PRECISION FUNCTION cylinder_soln_vel_x(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y,r,th
    
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
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y,r,th
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
    
! Flow past cylinder
! convert to spherical co-ordinates,
! NOTE: r & z cannot be 0 for the cylinder case - (0,0) is not part of domain.
    r = SQRT(x**2+y**2)
    th=acos(x/r)
    out = (rad_sphere**2/r**2 - 1d0)*cos(th)*sin(th)
    
  END FUNCTION cylinder_soln_vel_y
  
  
  DOUBLE PRECISION FUNCTION cylinder_soln_pressure(glob_i_in) RESULT(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y,r,th
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
    
! Flow past cylinder
! convert to spherical co-ordinates,
! NOTE: r & z cannot be 0 for the cylinder case - (0,0) is not part of domain.
    r = dsqrt(x**2+y**2)
    th = dacos(x/r)
    out = -2d0*dcos(th)/r
    
  END FUNCTION cylinder_soln_pressure
  
  !  GRADIENTS

  DOUBLE PRECISION FUNCTION cylinder_soln_Grad_vel_xx(glob_i_in) RESULT(out)
! Diff of vel_x wrt x
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y,r,th
    
    x=nodeCoord(glob_i_in,1)
    y=nodeCoord(glob_i_in,2)
    
! Flow past cylinder
! convert to spherical co-ordinates,
! NOTE: r & z cannot be 0 for the cylinder case - (0,0) is not part of domain.
    r = dsqrt(x**2+y**2)
    th = dacos(x/r)
    
    out = (x/r)*((-2d0*rad_sphere**2/r**3)*(dcos(th)**2 - 5d-1) + 1d0/r) + &
      (2d0*y/r**2)*((rad_sphere/r)**2 - 1d0)*dcos(th)*dsin(th)
    
  END FUNCTION cylinder_soln_Grad_vel_xx
  
  DOUBLE PRECISION FUNCTION cylinder_soln_Grad_vel_yx(glob_i_in) RESULT(out)
! Diff of vel_x wrt y
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y,r,th
    
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
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y,r,th
    
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
    INTEGER, INTENT(IN) :: glob_i_in
    DOUBLE PRECISION :: x,y,r,th
    
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
  
! END FUNCTIONS FOR CARTESIAN FLOW PAST A CYLINDER

! UNIFORM FLOW PAST A CYLINDER: ALL are v_sphere, 0d0, or unknown - not required.

! BEGIN FUNCTIONS FOR 2-D POISEUILLE FLOW

! END FUNCTIONS FOR 2-D POISEUILLE FLOW

  
! END FUNCTIONS FOR STOKES MODEL SOLUTION

END MODULE functions_module