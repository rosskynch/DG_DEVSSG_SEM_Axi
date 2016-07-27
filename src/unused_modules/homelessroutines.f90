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

! Calc diff with respect to z of tau
      DO j=0,N
	DO i=0,N
	  ij=i+j*(N+1)
	  temp1=0d0
	  temp2=0d0
	  DO p=0,N
	    pj=p+j*(N+1)
	    ip=i+p*(N+1)
	    temp1 = temp1 + tempTxx(pj,el)*d(i,p)
	    temp2 = temp2 + tempTxx(ip,el)*d(j,p) ! switched p and q to save a new loop
	  ENDDO
	  strongGradTxx(ij,el) = (dyde(i,j,el)*temp1 - dydp(i,j,el)*temp2)/jac(i,j,el)
	  temp1=0d0
	  temp2=0d0
	  DO p=0,N
	    pj=p+j*(N+1)
	    ip=i+p*(N+1)
	    temp1 = temp1 + tempTxy(pj,el)*d(i,p)
	    temp2 = temp2 + tempTxy(ip,el)*d(j,p) ! switched p and q to save a new loop
	  ENDDO
	  strongGradTxy(ij,el) = (dyde(i,j,el)*temp1 - dydp(i,j,el)*temp2)/jac(i,j,el)
	  temp1=0d0
	  temp2=0d0
	  DO p=0,N
	    pj=p+j*(N+1)
	    ip=i+p*(N+1)
	    temp1 = temp1 + tempTyy(pj,el)*d(i,p)
	    temp2 = temp2 + tempTyy(ip,el)*d(j,p) ! switched p and q to save a new loop
	  ENDDO
	  strongGradTyy(ij,el) = (dyde(i,j,el)*temp1 - dydp(i,j,el)*temp2)/jac(i,j,el)
	  temp1=0d0
	  temp2=0d0
	  DO p=0,N
	    pj=p+j*(N+1)
	    ip=i+p*(N+1)
	    temp1 = temp1 + tempTzz(pj,el)*d(i,p)
	    temp2 = temp2 + tempTzz(ip,el)*d(j,p) ! switched p and q to save a new loop
	  ENDDO
	  strongGradTzz(ij,el) = (dyde(i,j,el)*temp1 - dydp(i,j,el)*temp2)/jac(i,j,el)
	ENDDO
      ENDDO




SUBROUTINE create_local_to_global_sparse_map

  IMPLICIT NONE
  
  ALLOCATE(sparse_rowIndex(1:glob_bd_dim+3*npint+1))
    
  non_zero_count=0
  Ax_local_to_sparse=0
  Ay_local_to_sparse=0
  Bx_local_to_sparse=0
  By_local_to_sparse=0
  
  
  seen_node=.false.
! boundary rows of x-part of global matrix:
  DO i=1,npedg
    IF (.not.bdflag(1,i)) THEN
      sparse_rowIndex(non_dir_bd_map_x(i))=non_zero_count+1
      DO j=i,npedg
	IF (.not.bdflag(1,j)) THEN  
	  DO el=1,numelm
	    DO ij=0,NP1SQ-1
	      IF (mapg(ij,el).eq.i) THEN
		DO kl=0,NP1SQ-1
		  IF (mapg(kl,el).eq.j) THEN
		    IF (i.eq.j) THEN ! Always store the diagonal entry, even if it is zero.
		      IF (seen_node(i,j)) THEN
			Ax_local_to_sparse(ij,kl,el)=non_zero_count
		      ELSE
			non_zero_count = non_zero_count+1
			Ax_local_to_sparse(ij,kl,el) = non_zero_count
			seen_node(i,j)=.true.
		      ENDIF
		    ELSEIF (abs(A_x(ij,kl,el)).gt.1d-15) THEN
		      IF (seen_node(i,j)) THEN
			Ax_local_to_sparse(ij,kl,el)=non_zero_count
		      ELSE
			non_zero_count = non_zero_count+1
			Ax_local_to_sparse(ij,kl,el) = non_zero_count
			seen_node(i,j)=.true.
		      ENDIF
		    ENDIF
		  ENDIF
		ENDDO
	      ENDIF
	    ENDDO
	  ENDDO
	ENDIF
      ENDDO
      A_x_interior_loop1: DO j=npedg+1,nptot
	DO el=1,numelm
	  DO ij=0,NP1SQ-1
	    IF (mapg(ij,el).eq.i) THEN
	      DO kl=0,NP1SQ-1
		IF (mapg(kl,el).eq.j) THEN
		  IF (abs(A_x(ij,kl,el)).gt.1d-15) THEN
		    non_zero_count=non_zero_count+1
		    Ax_local_to_sparse(ij,kl,el)=non_zero_count
		    CYCLE A_x_interior_loop1 ! Exit as we will only find 1 contribution to this node - it is interior
		  ENDIF
		ENDIF
	      ENDDO
	    ENDIF
	  ENDDO
	ENDDO
      ENDDO A_x_interior_loop1
      B_x_interior_loop1: DO j=1,npint
	DO el=1,numelm
	  DO kl=0,NP1SQ-1
	    IF (mapg(kl,el).eq.i) THEN
	      DO ii=1,NM1
		DO jj=1,NM1
		  ij = ii+jj*NP1
		  IF (mapg(ij,el).eq.j) THEN
		    intij=ii+(jj-1)*NM1
		    IF (abs(B_x(intij,kl,el)).gt.1d-15) THEN
		      non_zero_count=non_zero_count+1
		      Bx_local_to_sparse(intij,kl,el)=non_zero_count
		      CYCLE B_x_interior_loop1 ! Exit as we will only find 1 contribution to this node - it is partly interior
		    ENDIF
		  ENDIF
		ENDDO
	      ENDDO
	    ENDIF  
	  ENDDO
	ENDDO
      ENDDO B_x_interior_loop1
    ENDIF  
  ENDDO ! Now done all boundary nodes of the x-part of the global matrix.

! Interior rows of the x-part of the global matrix
  DO i=npedg+1,nptot
    sparse_rowIndex(i-npedg+bd_rows_x)=non_zero_count+1

    A_x_interior_loop2: DO j=i,nptot
	DO el=1,numelm
	  DO ij=0,NP1SQ-1
	    IF (mapg(ij,el).eq.i) THEN
	      DO kl=0,NP1SQ-1
		IF (mapg(kl,el).eq.j) THEN
		  IF (abs(A_x(ij,kl,el)).gt.1d-15) THEN
		    non_zero_count=non_zero_count+1
		    Ax_local_to_sparse(ij,kl,el)=non_zero_count
		    CYCLE A_x_interior_loop2 ! Exit as we will only find 1 contribution to this node - it is interior
		  ENDIF
		ENDIF
	      ENDDO
	    ENDIF
	  ENDDO
	ENDDO
      ENDDO A_x_interior_loop2
      B_x_interior_loop2: DO j=1,npint
	DO el=1,numelm
	  DO kl=0,NP1SQ-1
	    IF (mapg(kl,el).eq.i) THEN
	      DO ii=1,NM1
		DO jj=1,NM1
		  ij = ii+jj*NP1
		  IF (mapg(ij,el).eq.j) THEN
		    intij=ii+(jj-1)*NM1
		    IF (abs(B_x(intij,kl,el)).gt.1d-15) THEN
		      non_zero_count=non_zero_count+1
		      Bx_local_to_sparse(intij,kl,el)=non_zero_count
		      CYCLE B_x_interior_loop2 ! Exit as we will only find 1 contribution to this node - it is partly interior
		    ENDIF
		  ENDIF
		ENDDO
	      ENDDO
	    ENDIF  
	  ENDDO
	ENDDO
      ENDDO B_x_interior_loop2
    ENDIF  
  ENDDO ! Now done all interior nodes of the x-part of the global matrix.
		    
! boundary rows of y-part of global matrix:
  seen_node=.false.
  DO i=1,npedg
    IF (.not.bdflag(2,i)) THEN
      sparse_rowIndex(non_dir_bd_map_y(i))=non_zero_count+1
      DO j=i,npedg
	IF (.not.bdflag(2,j)) THEN  
	  DO el=1,numelm
	    DO ij=0,NP1SQ-1
	      IF (mapg(ij,el).eq.i) THEN
		DO kl=0,NP1SQ-1
		  IF (mapg(kl,el).eq.j) THEN
		    IF (i.eq.j) THEN ! Always store the diagonal entry, even if it is zero.
		      IF (seen_node(i,j)) THEN
			Ay_local_to_sparse(ij,kl,el )= non_zero_count
		      ELSE
			non_zero_count = non_zero_count+1
			Ay_local_to_sparse(ij,kl,el) = non_zero_count
			seen_node(i,j)=.true.
		      ENDIF
		    ELSEIF (abs(A_y(ij,kl,el)).gt.1d-15) THEN
		      IF (seen_node(i,j)) THEN
			Ay_local_to_sparse(ij,kl,el) = non_zero_count
		      ELSE
			non_zero_count = non_zero_count+1
			Ay_local_to_sparse(ij,kl,el) = non_zero_count
			seen_node(i,j)=.true.
		      ENDIF
		    ENDIF
		  ENDIF
		ENDDO
	      ENDIF
	    ENDDO
	  ENDDO
	ENDIF
      ENDDO
      A_y_interior_loop1: DO j=npedg+1,nptot
	DO el=1,numelm
	  DO ij=0,NP1SQ-1
	    IF (mapg(ij,el).eq.i) THEN
	      DO kl=0,NP1SQ-1
		IF (mapg(kl,el).eq.j) THEN
		  IF (abs(A_y(ij,kl,el)).gt.1d-15) THEN
		    non_zero_count=non_zero_count+1
		    Ay_local_to_sparse(ij,kl,el)=non_zero_count
		    CYCLE A_y_interior_loop1 ! Exit as we will only find 1 contribution to this node - it is interior
		  ENDIF
		ENDIF
	      ENDDO
	    ENDIF
	  ENDDO
	ENDDO
      ENDDO A_y_interior_loop1
      B_y_interior_loop1: DO j=1,npint
	DO el=1,numelm
	  DO kl=0,NP1SQ-1
	    IF (mapg(kl,el).eq.i) THEN
	      DO ii=1,NM1
		DO jj=1,NM1
		  ij = ii+jj*NP1
		  IF (mapg(ij,el).eq.j) THEN
		    intij=ii+(jj-1)*NM1
		    IF (abs(B_y(intij,kl,el)).gt.1d-15) THEN
		      non_zero_count=non_zero_count+1
		      By_local_to_sparse(intij,kl,el)=non_zero_count
		      CYCLE B_y_interior_loop1 ! Exit as we will only find 1 contribution to this node - it is partly interior
		    ENDIF
		  ENDIF
		ENDDO
	      ENDDO
	    ENDIF  
	  ENDDO
	ENDDO
      ENDDO B_x_interior_loop1
    ENDIF  
  ENDDO ! Now done all boundary nodes of the y-part of the global matrix.

! Interior rows of the y-part of the global matrix
  DO i=npedg+1,nptot
    sparse_rowIndex(i-npedg+glob_bd_dim+npint)=non_zero_count+1

      A_y_interior_loop2: DO j=i,nptot
	DO el=1,numelm
	  DO ij=0,NP1SQ-1
	    IF (mapg(ij,el).eq.i) THEN
	      DO kl=0,NP1SQ-1
		IF (mapg(kl,el).eq.j) THEN
		  IF (abs(A_y(ij,kl,el)).gt.1d-15) THEN
		    non_zero_count=non_zero_count+1
		    Ay_local_to_sparse(ij,kl,el)=non_zero_count
		    CYCLE A_y_interior_loop2 ! Exit as we will only find 1 contribution to this node - it is interior
		  ENDIF
		ENDIF
	      ENDDO
	    ENDIF
	  ENDDO
	ENDDO
      ENDDO A_y_interior_loop2
      B_y_interior_loop2: DO j=1,npint
	DO el=1,numelm
	  DO kl=0,NP1SQ-1
	    IF (mapg(kl,el).eq.i) THEN
	      DO ii=1,NM1
		DO jj=1,NM1
		  ij = ii+jj*NP1
		  IF (mapg(ij,el).eq.j) THEN
		    intij=ii+(jj-1)*NM1
		    IF (abs(B_y(intij,kl,el)).gt.1d-15) THEN
		      non_zero_count=non_zero_count+1
		      By_local_to_sparse(intij,kl,el)=non_zero_count
		      CYCLE B_y_interior_loop2 ! Exit as we will only find 1 contribution to this node - it is partly interior
		    ENDIF
		  ENDIF
		ENDDO
	      ENDDO
	    ENDIF  
	  ENDDO
	ENDDO
      ENDDO B_y_interior_loop2
    ENDIF  
  ENDDO ! Now done all interior nodes of the y-part of the global matrix.
  
! Interior rows of the pressure-part of the global matrix
! These are all zero, so we need only store the diagonal entries
  DO i=1,npint
    sparse_rowIndex(i-glob_bd_dim+2*npint)=non_zero_count+1
    non_zero_count=non_zero_count+1
  ENDDO
  sparse_rowIndex(glob_bd_dim+3*npint+1)=non_zero_count+1  


END SUBROUTINE create_local_to_global_sparse_map

SUBROUTINE update_global_sparse_storage
  IMPLICIT NONE
  INTEGER :: el,ij,kl,k,l,intkl
 
  IF (movingmeshflag.eq.1) THEN
    DO el=1,numelm
      DO ij=0,NP1SQ-1
	DO kl=0,NP1SQ-1
	  temp(Ax_local_to_sparse(kl,ij,el)) = temp(Ax_local_to_sparse(kl,ij,el)) + A_x(kl,ij,el)
	  temp(Ay_local_to_sparse(kl,ij,el)) = temp(Ay_local_to_sparse(kl,ij,el)) + A_y(kl,ij,el)
	ENDDO
	DO k=1,NM1
	  DO l=1,NM1
	    intkl=k+(l-1)*NM1
	    temp(Bx_local_to_sparse(intkl,ij,el)) = temp(Bx_local_to_sparse(intkl,ij,el)) + B_x(intkl,ij,el)
	    temp(By_local_to_sparse(intkl,ij,el)) = temp(By_local_to_sparse(intkl,ij,el)) + B_y(intkl,ij,el)
	  ENDDO
	ENDDO
      ENDDO
    ENDDO
  ENDIF
END SUBROUTINE update_global_sparse_storage










SUBROUTINE initial_SEM_setup
  IMPLICIT NONE
  
  
END SUBROUTINE initial_SEM_setup

! OLD VERSION OF ROUTINE - MAY BE USEFUL LATER!!
  SUBROUTINE construct_zero_pressure_matrix
    IMPLICIT NONE
    INTEGER :: i,j,ij,kl,el,k,l,edg,gij
!     DOUBLE PRECISION, DIMENSION(1:(N-1)**2) :: tempvec2
    DOUBLE PRECISION, DIMENSION(1:(N-1)**2,1:numelm) :: tempvec,tempvec2
    DOUBLE PRECISION, DIMENSION(1:npedg,1:nptot) :: temp

!
! Calculate value contribution of each pressure test function q=h_k h_l at every element edge node.
! E.g. If two elements are joined by 1 edge, there will be a contribution from both elements at each point on
!      that joining edge
!
!      If 4 elements meet at 1 point, then that point will contain contributions from all 4 elements.
!
    temp1=0d0
    DO el=1,numelm
      DO i=0,N
	j=0
	ij=i+j*(N+1)
	DO k=1,N-1
	DO l=1,N-1
	  kl=k+(l-1)*(N-1)
	  temp1(ij,kl,el)=temp1(ij,kl,el)+evalh(i,k)*evalh(j,l)
	ENDDO
	ENDDO
	j=N
	ij=i+j*(N+1)
	DO k=1,N-1
	DO l=1,N-1
	  kl=k+(l-1)*(N-1)
	  temp1(ij,kl,el)=temp1(ij,kl,el)+evalh(i,k)*evalh(j,l)
	ENDDO
	ENDDO
      DO j=0,N
	i=0
	ij=i+j*(N+1)
	DO k=1,N-1
	DO l=1,N-1
	  kl=k+(l-1)*(N-1)
	  temp1(ij,kl,el)=temp1(ij,kl,el)+evalh(i,k)*evalh(j,l)
	ENDDO
	ENDDO
	i=N
	ij=i+j*(N+1)
	DO k=1,N-1
	DO l=1,N-1
	  kl=k+(l-1)*(N-1)
	  temp1(ij,kl,el)=temp1(ij,kl,el)+evalh(i,k)*evalh(j,l)
	ENDDO
	ENDDO
      ENDDO
    ENDDO
    
    tempglob=0d0
    DO el=1,numelm
      DO ij=
	  
 
    temp=0d0
    DO edg=1,npedg
      DO el=1,numelm
	DO i=0,N
	  j=0
	  ij=i+j*(N+1)
	  IF (edg.eq.mapg(ij,el)) THEN
	    DO k=1,N-1
	    DO l=1,N-1
	      kl=k+l*(N+1)
	      temp(edg,mapg(kl,el))=temp(edg,mapg(kl,el))+evalh(k,i)*evalh(l,j)
	    ENDDO
	    ENDDO
	  ENDIF
	  j=N
	  ij=i+j*(N+1)
	  IF (edg.eq.mapg(ij,el)) THEN
	    DO k=1,N-1
	    DO l=1,N-1
	      kl=k+l*(N+1)
	      temp(edg,mapg(kl,el))=temp(edg,mapg(kl,el))+evalh(k,i)*evalh(l,j)
	    ENDDO
	    ENDDO
	  ENDIF
	ENDDO
	DO j=0,N
	  i=0
	  ij=i+j*(N+1)
	  IF (edg.eq.mapg(ij,el)) THEN
	    DO k=1,N-1
	    DO l=1,N-1
	      kl=k+l*(N+1)
	      temp(edg,mapg(kl,el))=temp(edg,mapg(kl,el))+evalh(k,i)*evalh(l,j)
	    ENDDO
	    ENDDO
	  ENDIF
	  i=N
	  ij=i+j*(N+1)
	  IF (edg.eq.mapg(ij,el)) THEN
	    DO k=1,N-1
	    DO l=1,N-1
	      kl=k+l*(N+1)
	      temp(edg,mapg(kl,el))=temp(edg,mapg(kl,el))+evalh(k,i)*evalh(l,j)
	    ENDDO
	    ENDDO
	  ENDIF
	ENDDO	
      ENDDO
    ENDDO
	  
	  
! probs wrong: 
!     temp=0d0
!     do el=1,numelm
!       do k=1,n-1
!       do l=1,n-1
! 	do i=0,N
! 	  j=0
! 	  ij=i+j*(N+1)
! 	  temp(mapg(ij,el))=temp(mapg(ij,el)) + evalh(k,i)*evalh(l,j)
! 	  j=N
! 	  ij=i+j*(N+1)
! 	  temp(mapg(ij,el))=temp(mapg(ij,el)) + evalh(k,i)*evalh(l,j)
! 	enddo
! 	do j=1,N-1
! 	  i=0
! 	  ij=i+j*(N+1)
! 	  temp(mapg(ij,el))=temp(mapg(ij,el)) + evalh(k,i)*evalh(l,j)
! 	  i=N
! 	  ij=i+j*(N+1)
! 	  temp(mapg(ij,el))=temp(mapg(ij,el)) + evalh(k,i)*evalh(l,j)
! 	enddo
!       enddo
!       enddo
!     enddo
    

    tempvec=0d0
    IF (coordflag.eq.0) THEN
      DO el=1,numelm
	DO j=1,N-1
	DO i=1,N-1
	  ij=i+(j-1)*(N-1)
	  tempvec(ij,el) = tempvec(ij,el) + jac(i,j,el)*w(i)*w(j) + &
				  evalh(i,0)*evalh(j,0)*jac(0,0,el)*w(0)*w(0) + &
				  evalh(i,0)*evalh(j,N)*jac(0,N,el)*w(0)*w(N) + &
				  evalh(i,N)*evalh(j,0)*jac(N,0,el)*w(N)*w(0) + &
				  evalh(i,N)*evalh(j,N)*jac(N,N,el)*w(N)*w(N)
	ENDDO
	ENDDO
      ENDDO
    ELSEIF (coordflag.eq.1) THEN 
      DO el=1,numelm
	DO j=1,N-1
	DO i=1,N-1
	  ij=i+(j-1)*(N-1)
	  tempvec(ij,el) = tempvec(ij,el) + jac(i,j,el)*w(i)*w(j)*nodeCoord(mapg(i+j*(N+1),el),2) + &
				  evalh(i,0)*evalh(j,0)*jac(0,0,el)*w(0)*w(0)*nodeCoord(mapg(0,el),2) + &
				  evalh(i,0)*evalh(j,N)*jac(0,N,el)*w(0)*w(N)*nodeCoord(mapg(N*(N+1),el),2) + &
				  evalh(i,N)*evalh(j,0)*jac(N,0,el)*w(N)*w(0)*nodeCoord(mapg(N,el),2) + &
				  evalh(i,N)*evalh(j,N)*jac(N,N,el)*w(N)*w(N)*nodeCoord(mapg(N+N*(N+1),el),2)
	ENDDO
	ENDDO
      ENDDO
    ELSE
      print*, 'Error: No co-ordinate system specified...'
      print*, 'Stopping'
      STOP
    ENDIF
    
    
    tempvec2=0d0
    IF (coordflag.eq.0) THEN
      DO el=1,numelm
	DO j=1,N-1
	DO i=1,N-1
	  ij=i+(j-1)*(N-1)
	  gij=i+j*(N+1)
	  tempvec2(ij,el) = tempvec2(ij,el) + jac(i,j,el)*w(i)*w(j) + &
				  temp(mapg(0,el),mapg(gij,el))*jac(0,N,el)*w(0)*w(0) + &
				  temp(mapg(N*(N+1),el),mapg(gij,el))*jac(0,N,el)*w(0)*w(N) + &
				  temp(mapg(N,el),mapg(gij,el))*jac(N,0,el)*w(N)*w(0) + &
				  temp(mapg(N+N*(N+1),el),mapg(gij,el))*jac(N,N,el)*w(N)*w(N)
	ENDDO
	ENDDO
      ENDDO
    ELSEIF (coordflag.eq.1) THEN 
      DO el=1,numelm
	DO j=1,N-1
	DO i=1,N-1
	  ij=i+(j-1)*(N-1)
	  tempvec(ij,el) = tempvec(ij,el) + jac(i,j,el)*w(i)*w(j)*nodeCoord(mapg(i+j*(N+1),el),2) + &
				  temp(mapg(0,el),mapg(gij,el))*jac(0,0,el)*w(0)*w(0)*nodeCoord(mapg(0,el),2) + &
				  temp(mapg(N*(N+1),el),mapg(gij,el))*jac(0,N,el)*w(0)*w(N)*nodeCoord(mapg(N*(N+1),el),2) + &
				  temp(mapg(N,el),mapg(gij,el))*jac(N,0,el)*w(N)*w(0)*nodeCoord(mapg(N,el),2) + &
				  temp(mapg(N+N*(N+1),el),mapg(gij,el))*jac(N,N,el)*w(N)*w(N)*nodeCoord(mapg(N+N*(N+1),el),2)
	ENDDO
	ENDDO
      ENDDO
    ELSE
      print*, 'Error: No co-ordinate system specified...'
      print*, 'Stopping'
      STOP
    ENDIF
    
    
! ! Calculate the integral over the whole domain for the pressure test function, q = h_k h_l
!     DO el=1,numelm
!       DO kl=1,(N-1)**2
!       	tempvec2(kl) = tempvec2(kl) + tempvec(kl,el)
!       ENDDO
!     ENDDO
    Z_p=0d0
    DO el=1,numelm
      DO kl=1,(N-1)**2
	DO ij=1,(N-1)**2
	  Z_p(ij,kl,el) = Z_p(ij,kl,el) + tempvec2(kl,el)*tempvec(ij,el)
	ENDDO
      ENDDO
    ENDDO
    Z_p=param_alphaZ*Z_p
    
!     do el=1,numelm
!     do ij=1,(N-1)**2
!     do kl=1,(N-1)**2
!       print*,Z_p(ij,kl,el)-Z_p(kl,ij,el)
!       enddo
!       enddo
!       enddo
!       stop

  END SUBROUTINE construct_zero_pressure_matrix
  
  
     
      
!       SUBROUTINE create_local_to_global_sparse_map
! 
!   IMPLICIT NONE
!   
!   seen_node=.false.
!     
!   non_zero_count=0
!   seen_node=.false.
! ! boundary rows of x-part of global matrix:
!   DO i=1,npedg
!     IF (.not.bdflag(1,i)) THEN
!       sparse_rowIndex(non_dir_bd_map_x(i))=non_zero_count+1
!       DO j=i,npedg
! 	IF (.not.bdflag(1,j)) THEN  
! 	  DO el=1,numelm
! 	    DO ij=0,NP1SQ-1
! 	      IF (mapg(ij,el).eq.i) THEN
! 		DO kl=0,NP1SQ-1
! 		  IF (mapg(kl,el).eq.j) THEN
! 		    IF (i.eq.j) THEN ! Always store the diagonal entry, even if it is zero.
! 		      IF (seen_node(i,j)) THEN
! 			Ax_local_to_sparse(ij,kl,el)=non_zero_count
! 			sparse_global_matrix(non_zero_count) = sparse_global_matrix(non_zero_count) + A_x(ij,kl,el)
! 		      ELSE
! 			non_zero_count = non_zero_count+1
! 			Ax_local_to_sparse(ij,kl,el) = non_zero_count
! 			sparse_global_matrix(non_zero_count) = A_x(ij,kl,el)
! 			sparse_column_count(non_zero_count) = non_dir_bd_map_x(j)
! 			seen_node(i,j)=.true.
! 		      ENDIF
! 		    ELSEIF (abs(A_x(ij,kl,el)).gt.1d-15) THEN
! 		      IF (seen_node(i,j)) THEN
! 			Ax_local_to_sparse(ij,kl,el)=non_zero_count
! 			sparse_global_matrix(non_zero_count) = sparse_global_matrix(non_zero_count) + A_x(ij,kl,el)
! 		      ELSE
! 			non_zero_count = non_zero_count+1
! 			Ax_local_to_sparse(ij,kl,el) = non_zero_count
! 			sparse_global_matrix(non_zero_count) = A_x(ij,kl,el)
! 			sparse_column_count(non_zero_count) = non_dir_bd_map_x(j)
! 			seen_node(i,j)=.true.
! 		      ENDIF
! 		    ENDIF
! 		  ENDIF
! 		ENDDO
! 	      ENDIF
! 	    ENDDO
! 	  ENDDO
! 	ENDIF
!       ENDDO
!       A_x_interior_loop1: DO j=npedg+1,nptot
! 	DO el=1,numelm
! 	  DO ij=0,NP1SQ-1
! 	    IF (mapg(ij,el).eq.i) THEN
! 	      DO kl=0,NP1SQ-1
! 		IF (mapg(kl,el).eq.j) THEN
! 		  IF (abs(A_x(ij,kl,el)).gt.1d-15) THEN
! 		    non_zero_count=non_zero_count+1
! 		    Ax_local_to_sparse(ij,kl,el)=non_zero_count
! 		    sparse_global_matrix(non_zero_count) = A_x(ij,kl,el)
! 		    sparse_column_count(non_zero_count)=j-npedg+bd_numrows_x
! 		    CYCLE A_x_interior_loop1 ! Exit as we will only find 1 contribution to this node - it is interior
! 		  ENDIF
! 		ENDIF
! 	      ENDDO
! 	    ENDIF
! 	  ENDDO
! 	ENDDO
!       ENDDO A_x_interior_loop1
!       B_x_interior_loop1: DO j=1,npint
! 	DO el=1,numelm
! 	  DO kl=0,NP1SQ-1
! 	    IF (mapg(kl,el).eq.i) THEN
! 	      DO ii=1,NM1
! 		DO jj=1,NM1
! 		  ij = ii+jj*NP1
! 		  IF (mapg(ij,el).eq.j) THEN
! 		    intij=ii+(jj-1)*NM1
! 		    IF (abs(B_x(intij,kl,el)).gt.1d-15) THEN
! 		      non_zero_count=non_zero_count+1
! 		      Bx_local_to_sparse(intij,kl,el)=non_zero_count
! 		      sparse_global_matrix(non_zero_count) = B_x(intij,kl,el)
! 		      sparse_column_count(non_zero_count)=j+glob_bd_dim+2*npint
! 		      CYCLE B_x_interior_loop1 ! Exit as we will only find 1 contribution to this node - it is partly interior
! 		    ENDIF
! 		  ENDIF
! 		ENDDO
! 	      ENDDO
! 	    ENDIF  
! 	  ENDDO
! 	ENDDO
!       ENDDO B_x_interior_loop1
!     ENDIF  
!   ENDDO ! Now done all boundary nodes of the x-part of the global matrix.
! 
! ! Interior rows of the x-part of the global matrix
!   DO i=npedg+1,nptot
!     sparse_rowIndex(i-npedg+bd_rows_x)=non_zero_count+1
! 
!     A_x_interior_loop2: DO j=i,nptot
! 	DO el=1,numelm
! 	  DO ij=0,NP1SQ-1
! 	    IF (mapg(ij,el).eq.i) THEN
! 	      DO kl=0,NP1SQ-1
! 		IF (mapg(kl,el).eq.j) THEN
! 		  IF (abs(A_x(ij,kl,el)).gt.1d-15) THEN
! 		    non_zero_count=non_zero_count+1
! 		    Ax_local_to_sparse(ij,kl,el)=non_zero_count
! 		    sparse_global_matrix(non_zero_count) = A_x(ij,kl,el)
! 		    sparse_column_count(non_zero_count) = j-npedg+bd_numrows_x
! 		    CYCLE A_x_interior_loop2 ! Exit as we will only find 1 contribution to this node - it is interior
! 		  ENDIF
! 		ENDIF
! 	      ENDDO
! 	    ENDIF
! 	  ENDDO
! 	ENDDO
!       ENDDO A_x_interior_loop2
!       B_x_interior_loop2: DO j=1,npint
! 	DO el=1,numelm
! 	  DO kl=0,NP1SQ-1
! 	    IF (mapg(kl,el).eq.i) THEN
! 	      DO ii=1,NM1
! 		DO jj=1,NM1
! 		  ij = ii+jj*NP1
! 		  IF (mapg(ij,el).eq.j) THEN
! 		    intij=ii+(jj-1)*NM1
! 		    IF (abs(B_x(intij,kl,el)).gt.1d-15) THEN
! 		      non_zero_count=non_zero_count+1
! 		      Bx_local_to_sparse(intij,kl,el)=non_zero_count
! 		      sparse_global_matrix(non_zero_count) = B_x(intij,kl,el)
! 		      sparse_column_count(non_zero_count)=j+glob_bd_dim+2*npint
! 		      CYCLE B_x_interior_loop2 ! Exit as we will only find 1 contribution to this node - it is partly interior
! 		    ENDIF
! 		  ENDIF
! 		ENDDO
! 	      ENDDO
! 	    ENDIF  
! 	  ENDDO
! 	ENDDO
!       ENDDO B_x_interior_loop2
!     ENDIF  
!   ENDDO ! Now done all interior nodes of the x-part of the global matrix.
! 		    
! ! boundary rows of y-part of global matrix:
!   seen_node=.false.
!   DO i=1,npedg
!     IF (.not.bdflag(2,i)) THEN
!       sparse_rowIndex(non_dir_bd_map_y(i))=non_zero_count+1
!       DO j=i,npedg
! 	IF (.not.bdflag(2,j)) THEN  
! 	  DO el=1,numelm
! 	    DO ij=0,NP1SQ-1
! 	      IF (mapg(ij,el).eq.i) THEN
! 		DO kl=0,NP1SQ-1
! 		  IF (mapg(kl,el).eq.j) THEN
! 		    IF (i.eq.j) THEN ! Always store the diagonal entry, even if it is zero.
! 		      IF (seen_node(i,j)) THEN
! 			Ay_local_to_sparse(ij,kl,el )= non_zero_count
! 			sparse_global_matrix(non_zero_count) = sparse_global_matrix(non_zero_count) + A_y(ij,kl,el)
! 		      ELSE
! 			non_zero_count = non_zero_count+1
! 			Ay_local_to_sparse(ij,kl,el) = non_zero_count
! 			sparse_global_matrix(non_zero_count) = A_y(ij,kl,el)
! 			sparse_column_count(non_zero_count) = non_dir_bd_map_y(j)
! 			seen_node(i,j)=.true.
! 		      ENDIF
! 		    ELSEIF (abs(A_y(ij,kl,el)).gt.1d-15) THEN
! 		      IF (seen_node(i,j)) THEN
! 			Ay_local_to_sparse(ij,kl,el) = non_zero_count
! 			sparse_global_matrix(non_zero_count) = sparse_global_matrix(non_zero_count) + A_y(ij,kl,el)
! 		      ELSE
! 			non_zero_count = non_zero_count+1
! 			Ay_local_to_sparse(ij,kl,el) = non_zero_count
! 			sparse_global_matrix(non_zero_count) = A_y(ij,kl,el)
! 			sparse_column_count(non_zero_count) = non_dir_bd_map_y(j)
! 			seen_node(i,j)=.true.
! 		      ENDIF
! 		    ENDIF
! 		  ENDIF
! 		ENDDO
! 	      ENDIF
! 	    ENDDO
! 	  ENDDO
! 	ENDIF
!       ENDDO
!       A_y_interior_loop1: DO j=npedg+1,nptot
! 	DO el=1,numelm
! 	  DO ij=0,NP1SQ-1
! 	    IF (mapg(ij,el).eq.i) THEN
! 	      DO kl=0,NP1SQ-1
! 		IF (mapg(kl,el).eq.j) THEN
! 		  IF (abs(A_y(ij,kl,el)).gt.1d-15) THEN
! 		    non_zero_count=non_zero_count+1
! 		    Ay_local_to_sparse(ij,kl,el)=non_zero_count
! 		    sparse_global_matrix(non_zero_count) = A_y(ij,kl,el)
! 		    sparse_column_count(non_zero_count)=j - npedg + glob_bd_dim + npint
! 		    CYCLE A_y_interior_loop1 ! Exit as we will only find 1 contribution to this node - it is interior
! 		  ENDIF
! 		ENDIF
! 	      ENDDO
! 	    ENDIF
! 	  ENDDO
! 	ENDDO
!       ENDDO A_y_interior_loop1
!       B_y_interior_loop1: DO j=1,npint
! 	DO el=1,numelm
! 	  DO kl=0,NP1SQ-1
! 	    IF (mapg(kl,el).eq.i) THEN
! 	      DO ii=1,NM1
! 		DO jj=1,NM1
! 		  ij = ii+jj*NP1
! 		  IF (mapg(ij,el).eq.j) THEN
! 		    intij=ii+(jj-1)*NM1
! 		    IF (abs(B_y(intij,kl,el)).gt.1d-15) THEN
! 		      non_zero_count=non_zero_count+1
! 		      By_local_to_sparse(intij,kl,el)=non_zero_count
! 		      sparse_global_matrix(non_zero_count) = B_y(intij,kl,el)
! 		      sparse_column_count(non_zero_count)=j+glob_bd_dim+2*npint
! 		      CYCLE B_y_interior_loop1 ! Exit as we will only find 1 contribution to this node - it is partly interior
! 		    ENDIF
! 		  ENDIF
! 		ENDDO
! 	      ENDDO
! 	    ENDIF  
! 	  ENDDO
! 	ENDDO
!       ENDDO B_x_interior_loop1
!     ENDIF  
!   ENDDO ! Now done all boundary nodes of the y-part of the global matrix.
! 
! ! Interior rows of the y-part of the global matrix
!   DO i=npedg+1,nptot
!     sparse_rowIndex(i-npedg+glob_bd_dim+npint)=non_zero_count+1
! 
!       A_y_interior_loop2: DO j=i,nptot
! 	DO el=1,numelm
! 	  DO ij=0,NP1SQ-1
! 	    IF (mapg(ij,el).eq.i) THEN
! 	      DO kl=0,NP1SQ-1
! 		IF (mapg(kl,el).eq.j) THEN
! 		  IF (abs(A_y(ij,kl,el)).gt.1d-15) THEN
! 		    non_zero_count=non_zero_count+1
! 		    Ay_local_to_sparse(ij,kl,el)=non_zero_count
! 		    sparse_global_matrix(non_zero_count) = A_y(ij,kl,el)
! 		    sparse_column_count(non_zero_count)=j - npedg + glob_bd_dim + npint
! 		    CYCLE A_y_interior_loop2 ! Exit as we will only find 1 contribution to this node - it is interior
! 		  ENDIF
! 		ENDIF
! 	      ENDDO
! 	    ENDIF
! 	  ENDDO
! 	ENDDO
!       ENDDO A_y_interior_loop2
!       B_y_interior_loop2: DO j=1,npint
! 	DO el=1,numelm
! 	  DO kl=0,NP1SQ-1
! 	    IF (mapg(kl,el).eq.i) THEN
! 	      DO ii=1,NM1
! 		DO jj=1,NM1
! 		  ij = ii+jj*NP1
! 		  IF (mapg(ij,el).eq.j) THEN
! 		    intij=ii+(jj-1)*NM1
! 		    IF (abs(B_y(intij,kl,el)).gt.1d-15) THEN
! 		      non_zero_count=non_zero_count+1
! 		      By_local_to_sparse(intij,kl,el)=non_zero_count
! 		      sparse_global_matrix(non_zero_count) = B_y(intij,kl,el)
! 		      sparse_column_count(non_zero_count)=j+glob_bd_dim+2*npint
! 		      CYCLE B_y_interior_loop2 ! Exit as we will only find 1 contribution to this node - it is partly interior
! 		    ENDIF
! 		  ENDIF
! 		ENDDO
! 	      ENDDO
! 	    ENDIF  
! 	  ENDDO
! 	ENDDO
!       ENDDO B_y_interior_loop2
!     ENDIF  
!   ENDDO ! Now done all interior nodes of the y-part of the global matrix.
!   
! ! Interior rows of the pressure-part of the global matrix
! ! These are all zero, so we need only store the diagonal entries
!   DO i=1,npint
!     sparse_rowIndex(i-glob_bd_dim+2*npint)=non_zero_count+1
!     non_zero_count=non_zero_count+1
!     sparse_global_matrix(non_zero_count) = 0d0
!     sparse_column_count(non_zero_count)=i+glob_bd_dim+2*npint
!   ENDDO
!   sparse_rowIndex(glob_bd_dim+3*npint+1)=non_zero_count+1  
! 
! 
! END SUBROUTINE create_local_to_global_sparse_map