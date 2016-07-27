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

MODULE lapack_solver_module
 ! Will require linking to blas/lapack at compilation.
 ! Cholesky fact: DPOTRF
 ! Solve: DPOTRS
 ! Invert: DPOTRI
 ! Estimate condition number: DPOCON
  USE shared_data
  USE functions_module
  
  IMPLICIT NONE
   
  INTEGER, ALLOCATABLE ::  schur_bd_local_to_boundary_map(:,:),&
			   schur_bd_local_map(:,:),&
			   normalise_from_el(:),&
			   normalise_from_local_vertex(:,:),&
			   normalise_to_local_vertex(:,:),&
			   normalise_pressure_element_order(:)

	     
  DOUBLE PRECISION, ALLOCATABLE :: schur_boundary_matrix(:,:),&
				   schur_crossover_matrix(:,:),&
				   schur_local_crossover_matrix(:,:,:),&
				   schur_internal_matrix(:,:,:),&
				   schur_fact_internal_matrix(:,:,:), &
				   schur_boundary_rhs(:),&
				   schur_internal_rhs(:,:),&
				   internal_diag_precon(:,:)


  CONTAINS
  
  SUBROUTINE lapack_release_memory
    IMPLICIT NONE
    DEALLOCATE( &
		schur_bd_local_map, &
		schur_bd_local_to_boundary_map, &
		schur_boundary_matrix, &
		schur_crossover_matrix, &
		schur_local_crossover_matrix, &
		schur_fact_internal_matrix,&
		schur_internal_matrix, &
		schur_boundary_rhs, &
		schur_internal_rhs, &
		internal_diag_precon,&
		normalise_from_el, &
		normalise_from_local_vertex, &
		normalise_to_local_vertex, &
		normalise_pressure_element_order &
		)
  END SUBROUTINE lapack_release_memory
  
  SUBROUTINE initialise_lapack_schur_comp
! This routine will construct the schur decomposition of our matrix.
! The 3 sub matrices are stored as:
! - schur_boundary_matrix - containing all boundary-boundary entries - stored globally
! - schur_crossover_matrix - containing the internal-boundary entries - stored globally (the transpose gives the other corner of the matrix).
! - schur_internal_matrix - containing the internal-internal entries - stored element-wise.
! 
! It will also construct the required RHS vectors:
! - schur_boundary_rhs - RHS boundary nodes - stored globally
! - schur_internal_rhs - RHS internal nodes - stored element-wise
!
! Finally, memory will be allocated here. memory should be released using:
!          CALL release_mem_lapack_schur_comp

    IMPLICIT NONE
    INTEGER :: globi,globj,el,info,i,jj,j,ij,intij,k,l,intkl,kl,minresconv,&
	       ipiv_int(1:threeNM1SQ,numelm),&
	       ipiv_bd(1:bd_numrows_x+bd_numrows_y),&
	       iwork_bd(1:bd_numrows_x+bd_numrows_y)
	       
	       
    DOUBLE PRECISION :: ferr_bd, berr_bd,test_pressure_const,&
			work_bd(3*(bd_numrows_x+bd_numrows_y)),&
			internal_vec1(1:threeNM1SQ,numelm), &
			internal_vec2(1:threeNM1SQ,numelm), &
			local_boundary_vec1(8*N,numelm), &
			local_boundary_vec2(8*N,numelm), &
			local_boundary_matrix1(1:8*N,1:8*N,numelm), &
			global_matrix1(1:bd_numrows_x+bd_numrows_y,1:bd_numrows_x+bd_numrows_y), &
			global_vec1(1:bd_numrows_x+bd_numrows_y),&
			local_crossover1(1:threeNM1SQ,1:8*N,numelm),&
			final_RHS_boundary(1:bd_numrows_x+bd_numrows_y),&
			final_RHS_internal(1:threeNM1SQ,numelm),&
			boundary_solution(1:bd_numrows_x+bd_numrows_y),&
			internal_solution(1:threeNM1SQ,numelm),&
			final_schur_comp_matrix(1:bd_numrows_x+bd_numrows_y,1:bd_numrows_x+bd_numrows_y),&
			fact_schur_comp_matrix(1:bd_numrows_x+bd_numrows_y,1:bd_numrows_x+bd_numrows_y),&
			boundary_rhs(1:bd_numrows_x+bd_numrows_y),&
			temp_local_vec1(1:threeNM1SQ,numelm),&
			temp_local_vec2(1:threeNM1SQ,numelm),&
			temp_local_bd_vec(1:8*N,numelm),&
			temp_boundary_vec1(1:bd_numrows_x+bd_numrows_y),&
			temp_boundary_vec2(1:bd_numrows_x+bd_numrows_y),&
			temp_local_crossover_matrix(1:threeNM1SQ,1:8*N,numelm),&
			temp_local_boundary_matrix(1:8*N,1:8*N,numelm),&
			fact_boundary_matrix(1:bd_numrows_x+bd_numrows_y,1:bd_numrows_x+bd_numrows_y),&
			test_pressure_sol(1:NM1SQ,numelm),&
			tempvec(1:bd_numrows_x+bd_numrows_y+1),&
			temp_global_crossover(1:3*npint,1:bd_numrows_x+bd_numrows_y+1),&
			global_internal(1:3*npint,1:3*npint),&
			global_int_vec(1:3*npint),&
			t1(1:threeNM1SQ),&
			t2(1:threeNM1SQ),tmat(1:threeNM1SQ,1:threeNM1SQ)
			
			


! Calculate required dimensions, assign memory and create our boundary map.
    
!     ALLOCATE(non_dir_bd_map_x(1:npedg))
!     ALLOCATE(non_dir_bd_map_y(1:npedg))
    ALLOCATE(schur_bd_local_map(0:NP1SQ-1,numelm)) 
    ALLOCATE(schur_bd_local_to_boundary_map(1:8*N,numelm))
    CALL create_schur_maps
    ALLOCATE(schur_boundary_matrix(1:global_bd_dim,1:global_bd_dim))
    ALLOCATE(schur_crossover_matrix(1:3*npint,1:global_bd_dim+1))
    ALLOCATE(schur_local_crossover_matrix(1:threeNM1SQ,1:8*N,numelm))
    ALLOCATE(schur_internal_matrix(1:threeNM1SQ,1:threeNM1SQ,numelm))
    ALLOCATE(schur_fact_internal_matrix(1:threeNM1SQ,1:threeNM1SQ,numelm))
    ALLOCATE(schur_boundary_rhs(1:global_bd_dim))
    ALLOCATE(schur_internal_rhs(1:threeNM1SQ,numelm))
    ALLOCATE(internal_diag_precon(1:threeNM1SQ,numelm))
    ALLOCATE(normalise_from_el(1:numelm))
    ALLOCATE(normalise_from_local_vertex(1:numelm,2))
    ALLOCATE(normalise_to_local_vertex(1:numelm,2))
    ALLOCATE(normalise_pressure_element_order(1:numelm))

! BEGIN check on maps
!     local_boundary_vec1=1
!     CALL assemble_local_to_global_boundary_vec(local_boundary_vec1,global_vec1)
!     tempvec(1:bd_numrows_x+bd_numrows_y) = global_vec1
!     DO i=1,npedg
!       IF (.not.bdflag(1,i)) print*,tempvec(non_dir_bd_map_x(i)),mult(i),bdflag(1,i)
!     ENDDO
!     DO i=1,npedg
!       IF (.not.bdflag(2,i)) print*,tempvec(non_dir_bd_map_y(i)),mult(i),bdflag(2,i)
!     ENDDO
!     
!     DO i=1,npedg
!       IF (.not.bdflag(1,i)) global_vec1(non_dir_bd_map_x(i))=mult(i)
!       IF (.not.bdflag(2,i)) global_vec1(non_dir_bd_map_y(i))=mult(i)
!     ENDDO
!     
!     CALL scatter_global_boundary_to_local_vector(global_vec1,local_boundary_vec1)
!     DO el=1,numelm
!     DO i=0,NP1SQ-1
!       IF (mapg(i,el).LE.npedg) THEN
!       IF (.not.bdflag(1,mapg(i,el))) print*,el,schur_bd_local_map(i,el),local_boundary_vec1(schur_bd_local_map(i,el),el)- &
! 	dfloat(mult(mapg(i,el))),bdflag(1,mapg(i,el))
!       ENDIF
!     ENDDO
!     DO i=0,NP1SQ-1
!       IF (mapg(i,el).LE.npedg) THEN
!       IF (.not.bdflag(2,mapg(i,el))) print*,el,schur_bd_local_map(i,el),local_boundary_vec1(4*N+schur_bd_local_map(i,el),el)-&
! 	      dfloat(mult(mapg(i,el))),bdflag(2,mapg(i,el))
!       ENDIF
!     ENDDO
!     ENDDO
!     
!     DO el=1,numelm
!       DO j=0,N
! 	DO i=0,N
! 	  print*,el,schur_bd_local_map(i+j*NP1,el),i,j,i+j*NP1
! 	ENDDO
!       ENDDO
!     ENDDO
!       
!     stop
! END check on maps
    CALL setup_pressure_normalisation_connectivity
    
    
    
! Build matrices:
    CALL build_boundary_matrix(schur_boundary_matrix)
    CALL build_local_crossover_matrices(schur_local_crossover_matrix)
    CALL build_crossover_matrix(schur_crossover_matrix)
    CALL build_local_internal_matrices(schur_internal_matrix)
! Diagonal preconditioner for internal minres routine
    CALL calc_internal_diag_preconditioner(schur_internal_matrix,internal_diag_precon)
    

    
    temp_global_crossover=0d0
    DO el=1,numelm
      DO i=1,NM1
	DO jj=1,NM1
	  ij=i+jj*NP1
	  intij=i+(jj-1)*NM1
	  globi = mapg(ij,el)-npedg
	  DO j=1,8*N
	    globj = schur_bd_local_to_boundary_map(j,el)
	    
	    temp_global_crossover(globi,globj) = temp_global_crossover(globi,globj) &
		    + schur_local_crossover_matrix(intij,j,el)
		    
	    temp_global_crossover(globi+npint,globj) = temp_global_crossover(globi+npint,globj) &
		  + schur_local_crossover_matrix(NM1SQ + intij,j,el)
		  
	    temp_global_crossover(globi+2*npint,globj) = temp_global_crossover(globi+2*npint,globj) &
		  + schur_local_crossover_matrix(2*NM1SQ + intij,j,el)
	  ENDDO
	ENDDO
      ENDDO
    ENDDO
    
    global_internal=0d0
    DO el=1,numelm
      DO i=1,NM1
      DO j=1,NM1
	ij=i+j*NP1
	intij=i+(j-1)*NM1
	globi=mapg(ij,el)-npedg
	DO k=1,NM1
	DO l=1,NM1
	  kl=k+l*NP1
	  intkl=k+(l-1)*NM1
	  globj=mapg(kl,el)-npedg
	  global_internal(globi,globj)= schur_internal_matrix(intij,intkl,el)

	  globj=globj+npint
	  intkl=intkl+NM1SQ
	  global_internal(globi,globj)= schur_internal_matrix(intij,intkl,el)

	  globj=globj+npint
	  intkl=intkl+NM1SQ
	  global_internal(globi,globj)= schur_internal_matrix(intij,intkl,el)

	ENDDO
	ENDDO
	globi=globi+npint
	intij=intij+NM1SQ
	DO k=1,NM1
	DO l=1,NM1
	  kl=k+l*NP1
	  intkl=k+(l-1)*NM1
	  globj=mapg(kl,el)-npedg
	  global_internal(globi,globj)= schur_internal_matrix(intij,intkl,el)

	  globj=globj+npint
	  intkl=intkl+NM1SQ
	  global_internal(globi,globj)= schur_internal_matrix(intij,intkl,el)

	  globj=globj+npint
	  intkl=intkl+NM1SQ
	  global_internal(globi,globj)= schur_internal_matrix(intij,intkl,el)

	ENDDO
	ENDDO
	globi=globi+npint
	intij=intij+NM1SQ
	DO k=1,NM1
	DO l=1,NM1
	  kl=k+l*NP1
	  intkl=k+(l-1)*NM1
	  globj=mapg(kl,el)-npedg
	  global_internal(globi,globj)= schur_internal_matrix(intij,intkl,el)

	  globj=globj+npint
	  intkl=intkl+NM1SQ
	  global_internal(globi,globj)= schur_internal_matrix(intij,intkl,el)

	  globj=globj+npint
	  intkl=intkl+NM1SQ
	  global_internal(globi,globj)= schur_internal_matrix(intij,intkl,el)

	ENDDO
	ENDDO
      ENDDO
      ENDDO
    ENDDO

global_internal=0d0
DO el=1,numelm
      DO l=1,NM1
      DO k=1,NM1
	kl = k+l*NP1
	intkl = k+(l-1)*NM1
	globj = mapg(kl,el)-npedg

	DO j=1,NM1
	DO i=1,NM1
	  ij = i+j*NP1
	  intij = i+(j-1)*NM1
	  globi = mapg(ij,el)-npedg

	  global_internal(globi, globj) = A_x(ij,kl,el) !+ const*Mv_x(ij,ij,el)
	  global_internal(globi +npint, globj +npint) = A_y(ij,kl,el)! + !const*Mv_y(ij,ij,el)
	  global_internal(globi +2*npint, globj) = B_x(intij,kl,el)
	  global_internal(globi +2*npint, globj +npint) = B_y(intij,kl,el)
	  global_internal(globi, globj +2*npint) = B_x(intkl,ij,el)
	  global_internal(globi + npint, globj +2*npint) = B_y(intkl,ij,el)
	  
	ENDDO
	ENDDO
      ENDDO
      ENDDO
enddo

!     stop
    
    

    
! Zero the additional row/column in each global matrix
!     schur_boundary_matrix(global_bd_dim+1,1:global_bd_dim+1)=0d0
!     schur_boundary_matrix(1:global_bd_dim,global_bd_dim+1)=0d0
!     schur_crossover_matrix(1:3*npint,global_bd_dim+1)=0d0

    
! Build RHS vectors:
    CALL build_boundary_rhs(schur_boundary_rhs)
    CALL build_internal_rhs(schur_internal_rhs)
    
        DO el=1,numelm
      DO i=1,NM1
      DO j=1,NM1
	ij=i+j*NP1
	intij=i+(j-1)*NM1
	globi=mapg(ij,el)-npedg
	global_int_vec(globi) = schur_internal_rhs(intij,el)
	global_int_vec(globi+npint) = schur_internal_rhs(intij+NM1SQ,el)
	global_int_vec(globi+2*npint) = schur_internal_rhs(intij+2*NM1SQ,el)
      ENDDO
      ENDDO
    ENDDO
	
     open(321,FILE='McG.dat')
    open(4321,FILE='fi.dat')
    open(54321,FILE='MiG.dat')   
open(3211,FILE='MbG.dat')
open(3212,FILE='fb.dat')
    !     DO i=1,3*npint
!       write(*,*) (schur_crossover_matrix(i,j)-temp_global_crossover(i,j),j=1,global_bd_dim)
!     ENDDO
    
    DO i=1,3*npint
      write(321,*) (schur_crossover_matrix(i,j),j=1,global_bd_dim)
    ENDDO
    DO i=1,3*npint
!       write(4321,*) (temp_global_crossover(i,j),j=1,global_bd_dim)
      write(4321,*) global_int_vec(i)
    ENDDO
    DO i=1,3*npint
      write(54321,*) (global_internal(i,j),j=1,3*npint)
    ENDDO
    DO i=1, global_bd_dim
      write(3211,*) (schur_boundary_matrix(i,j),j=1,global_bd_dim)
    ENDDO
    DO i=1, global_bd_dim
      write(3212,*) schur_boundary_rhs(i)
    ENDDO
    close(321)
    close(4321)
    close(54321)
    close(3211)
    close(3212)
!     stop
    
! Factorise the local internal matrices
! For a moving mesh this only needs to be done once for non-moving elements
! For a fixed mesh it is only needed once.

    schur_fact_internal_matrix = schur_internal_matrix
    CALL fact_all_internal(schur_fact_internal_matrix,ipiv_int)
!     
    open(123,FILE='Mi.dat')
    open(1234,FILE='Mc.dat')
    open(12345,FILE='C.dat')
    do el=1,numelm
    DO i=1,threeNM1SQ
      write(123,*) (schur_internal_matrix(i,j,el),j=1,threeNM1SQ)
    ENDDO
    enddo
    do el=1,numelm
    DO i=1,threeNM1SQ
      write(1234,*) (schur_local_crossover_matrix(i,j,el),j=1,8*N)
    ENDDO
    enddo
! ! Calculate local_crossover1 = inv(M_i)*M_c
    local_crossover1=schur_local_crossover_matrix ! local X_1 in workings
    CALL solve_all_int_mult_cross_matrix(schur_internal_matrix, &
			schur_fact_internal_matrix,local_crossover1,ipiv_int) ! STEP 1.

    
!       SUBROUTINE MINRES_internal_only(Mat_in,X_in,RHS_in,diag_precon,minresconv)
    

!     local_crossover1=0d0
!     DO el=1,numelm
!       DO j=1,8*N
! 	IF (schur_bd_local_to_boundary_map(j,el).gt.global_bd_dim) CYCLE
! 	print*,j,schur_bd_local_to_boundary_map(j,el)
! 	
! 	t1=0d0
! 	t2=schur_local_crossover_matrix(1:threeNM1SQ,j,el)
! 	tmat=schur_internal_matrix(1:threeNM1SQ,1:threeNM1SQ,el)
! 	CALL MINRES_internal_only(tmat,t1,t2, internal_diag_precon(1:threeNM1SQ,el), minresconv)
! 	print*,minresconv
! 	local_crossover1(1:threeNM1SQ,j,el) = t1
!       ENDDO
!     ENDDO
      
    do el=1,numelm
    DO i=1,threeNM1SQ
      write(12345,*) (local_crossover1(i,j,el),j=1,8*N)
    ENDDO
    enddo
    close(123)
    close(1234)
    close(12345)
!     stop

!     DO el=1,numelm
!     DO i=1,threeNM1SQ
!       print*, schur_internal_rhs(i,el)!,schur_internal_rhs(NM1SQ+i,el),schur_internal_rhs(2*NM1SQ+i,el)
!     enddo
!     enddo
! stop

    IF ( global_bd_dim.gt.0 ) THEN
! Calculate local_boundary_vec1 = transpose(local_crossover1)*f_i = transpose(M_c)*inv(M_i)*f_i
    local_boundary_vec1=0d0
    CALL element_wise_mv_mult('T', local_crossover1, threeNM1SQ, 8*N, 1d0, &
			      schur_internal_rhs, threeNM1SQ, 0d0, local_boundary_vec1, 8*N)
	
    DO el=1,numelm
    DO i=0,NP1SQ-1
    IF (mapg(i,el).le.npedg) then
      print*, el,local_boundary_vec1(schur_bd_local_map(i,el),el),local_boundary_vec1(4*N+schur_bd_local_map(i,el),el),&
      schur_bd_local_to_boundary_map(schur_bd_local_map(i,el),el), &
			      schur_bd_local_to_boundary_map(4*N+schur_bd_local_map(i,el),el),schur_bd_local_map(i,el),i,&
			      nodecoord(mapg(i,el),1),nodecoord(mapg(i,el),2)
    endif
    enddo
    enddo
!     stop

! Make local_boundary_vec1 global in global_vec1
!     do el=1,numelm
!     do ij=0,NP1SQ-1
!       if (mapg(ij,el).le.npedg) then
!       local_boundary_vec1(schur_bd_local_map(ij,el),el)=dfloat(mult(mapg(ij,el)))
!       local_boundary_vec1(4*N+schur_bd_local_map(ij,el),el)=dfloat(mult(mapg(ij,el)))
!       endif
!       enddo
!       enddo
    global_vec1=0d0
    CALL assemble_local_to_global_boundary_vec(local_boundary_vec1,global_vec1) ! Goes wrong here.
!     tempvec(1:global_bd_dim)=global_vec1
!     DO i=1,npedg
!       print*,tempvec(non_dir_bd_map_x(i)),tempvec(non_dir_bd_map_y(i)),mult(i),bdflag(1,i),bdflag(2,i)
!       enddo
!       stop
        DO i=1,bd_numrows_y
	  print*,global_vec1(bd_numrows_x+i),bd_numrows_x+i
	enddo

! Calculate final_RHS_boundary = f_b - global_vec1
    final_RHS_boundary = schur_boundary_rhs - global_vec1
    
    DO i=1,bd_numrows_x
      print*,'x',final_RHS_boundary(i)
    ENDDO
    DO i=1,bd_numrows_y
      print*,'y',final_RHS_boundary(bd_numrows_x+i)
    ENDDO
!     stop
! Calculate local_boundary_matrix1 = transpose(M_c)*local_crossover1 = transpose(M_c)*inv(M_i)*M_c
    local_boundary_matrix1=0d0
    CALL element_wise_mm_mult('T','N',schur_local_crossover_matrix, threeNM1SQ, 8*N,&
				      local_crossover1, threeNM1SQ, 8*N, &
				      1d0, 0d0, local_boundary_matrix1, 8*N, 8*N)

! Make local_boundary_matrix global in global_matrix1
    CALL assemble_local_to_global_boundary_matrix(local_boundary_matrix1,global_matrix1)
    
! Form final schur complement matrix. final_schur_comp_matrix = M_b - global_matrix1
    final_schur_comp_matrix = schur_boundary_matrix - global_matrix1
    
    open(321,FILE='schurcomp.dat')
    DO i=1,global_bd_dim
      write(321,*) (final_schur_comp_matrix(i,j),j=1,global_bd_dim)
    ENDDO
    close(321)
    
! Solve boundary system:
! Copy our final RHS and Schur comp matrices to solution memory:
    fact_boundary_matrix = final_schur_comp_matrix
    boundary_solution = final_RHS_boundary
 

! LU fact:
    CALL DGETRF( global_bd_dim, global_bd_dim, fact_boundary_matrix, &
		    global_bd_dim, ipiv_bd, info )

! Solve using LU fact:
    CALL DGETRS( 'N', global_bd_dim, 1, fact_boundary_matrix, &
		    global_bd_dim, ipiv_bd, boundary_solution, &
		      global_bd_dim, info )
		      
! Refine solution:
!     CALL DGERFS( 'N', global_bd_dim, 1, final_schur_comp_matrix, global_bd_dim,&
! 		    fact_boundary_matrix, global_bd_dim, ipiv_bd,&
! 		    final_RHS_boundary, global_bd_dim, &
! 		    boundary_solution, global_bd_dim, &
! 		      ferr_bd, berr_bd, work_bd, iwork_bd, info )
!     print*,ferr_bd,berr_bd

    ENDIF
!     tempvec(1:global_bd_dim)=boundary_solution 
!     tempvec(global_bd_dim+1)=-999
!     DO i=1,npedg
!       
!       print*,tempvec(non_dir_bd_map_x(i)),tempvec(non_dir_bd_map_y(i)),mult(i),bdflag(1,i),bdflag(2,i)
!     ENDDO
!           stop
    
!     DO i=1,npedg
!       IF (bdflag(1,i)) THEN
! 	print*,boundary_x(i)
!       ELSE
! 	print*,boundary_solution(non_dir_bd_map_x(i))
!       ENDIF
!     ENDDO
!     DO i=1,npedg
!       IF (bdflag(2,i)) THEN
! 	print*,boundary_y(i)
!       ELSE
! 	print*,boundary_solution(non_dir_bd_map_y(i))
!       ENDIF
!     ENDDO
!     stop
      
    
! Solve for the interior nodes:

! internal_vec1 = inv(M_i)*f_i
    internal_vec1 = schur_internal_rhs
    print*,'here?'
    CALL solve_all_internal(schur_internal_matrix,schur_fact_internal_matrix,internal_vec1,ipiv_int)
    
    
   
! Make x_b, the boundary solution, local in local_boundary_vec2
    CALL scatter_global_boundary_to_local_vector(boundary_solution,local_boundary_vec2)
    
! Calculate internal_vec2 = local_crossover1*local_boundary_vec2 = inv(M_i)*M_c*x_b
    CALL element_wise_mv_mult('N', local_crossover1, threeNM1SQ, 8*N, 1d0, &
			      local_boundary_vec2, 8*N, 0d0, internal_vec2, threeNM1SQ)
    
! Calculate the internal nodes solution, internal_solution = internal_vec1 - internal_vec2
    internal_solution = internal_vec1 - internal_vec2
    
    CALL integrate_pressure_on_domain(internal_solution(2*NM1SQ+1:threeNM1SQ,1:numelm),test_pressure_const)
    

    print*,test_pressure_const

!     
!     internal_vec1=0d0
!     CALL element_wise_mv_mult('N', schur_internal_matrix, threeNM1SQ, threeNM1SQ, 1d0, &
! 			      internal_solution, threeNM1SQ, 0d0, internal_vec1, threeNM1SQ)
!     
!     do el=1,numelm
!     do i=1,threeNM1SQ
!     print*, internal_solution(i,el) - schur_internal_rhs(i,el)
!     ENDDO
!     ENDDO
!     stop
    
    DO el=1,numelm
      DO i=1,NM1
	DO j=1,NM1
	  intij=i+(j-1)*NM1
	  ij=i+j*NP1
	  V_x(mapg(ij,el)) = internal_solution(intij,el)
	  V_y(mapg(ij,el)) = internal_solution(NM1SQ+intij,el)
	  pressure(mapg(ij,el)-npedg)= internal_solution(2*NM1SQ+intij,el)
	ENDDO
      ENDDO
    ENDDO
    DO i=1,npedg
      IF (bdflag(1,i)) THEN
	V_x(i) = boundary_x(i)
      ELSE
	V_x(i) = boundary_solution(non_dir_bd_map_x(i))
      ENDIF
      IF (bdflag(2,i)) THEN
	V_y(i) = boundary_y(i)
      ELSE
	V_y(i) = boundary_solution(non_dir_bd_map_y(i))
      ENDIF
    ENDDO
    
    DO i=1,nptot
      print*,V_x(i),V_y(i),mult(i),bdflag(1,i),bdflag(2,i)
    ENDDO
    DO i=1,npint
      print*,pressure(i)
    ENDDO
!     stop
    
    call lapack_release_memory
  END SUBROUTINE initialise_lapack_schur_comp
  
  

  
  SUBROUTINE create_schur_maps
! essentially removes the dirichlet nodes and creates a new boundary-only map to use in our schur decomposition.
! any dirichlet nodes are mapped to the last+1 node of the matrix (which will exist but will not be passed to the solvers)
    IMPLICIT NONE
    INTEGER :: i,j,k,l,kl,el,NtNP1,Nt4,local_counter

    
    
! a few useful numbers to save lots of wasted operations.
    NtNP1=N*NP1 !N times N plus 1
    Nt4=4*N
    
    schur_bd_local_map=-1
    DO el=1,numelm
      local_counter=0
      DO kl=0,NP1SQ-1
	IF (mapg(kl,el).LE.npedg) THEN
	  local_counter=local_counter+1
	  schur_bd_local_to_boundary_map(local_counter,el)=non_dir_bd_map_x(mapg(kl,el))
	  schur_bd_local_to_boundary_map(Nt4+local_counter,el)=non_dir_bd_map_y(mapg(kl,el))
	  schur_bd_local_map(kl,el)=local_counter
	ENDIF
      ENDDO
! ! loop over boundary nodes - ie around the edge of the element.
! !	  l=0
!       DO kl=0,N
! 	local_counter=local_counter+1
! 	schur_bd_local_to_boundary_map(local_counter,el)=non_dir_bd_map_x(mapg(kl,el))
! 	schur_bd_local_to_boundary_map(Nt4+local_counter,el)=non_dir_bd_map_y(mapg(kl,el))
! 	schur_bd_local_map(kl,el)=local_counter
!       ENDDO
! 
! 
! !	 k=N 
!       DO l=1,N
! 	kl=k+NtNP1 ! Already done k=0,l=N
! 	local_counter=local_counter+1
! 	schur_bd_local_to_boundary_map(local_counter,el)=non_dir_bd_map_x(mapg(kl,el))
! 	schur_bd_local_to_boundary_map(Nt4+local_counter,el)=non_dir_bd_map_y(mapg(kl,el))
! 	schur_bd_local_map(kl,el)=local_counter
!       ENDDO
!       
! !	  k=N
!       DO l=N,1,-1 ! Already done k=N,l=0
! 	kl=N+l*NP1
! 	local_counter=local_counter+1
! 	schur_bd_local_to_boundary_map(local_counter,el)=non_dir_bd_map_x(mapg(kl,el))
! 	schur_bd_local_to_boundary_map(Nt4+local_counter,el)=non_dir_bd_map_y(mapg(kl,el))
! 	schur_bd_local_map(kl,el)=local_counter
!       ENDDO
! 
! !	  l= 0
!       DO kl=1,N ! Already done k=0,l=0
! 	local_counter=local_counter+1
! 	schur_bd_local_to_boundary_map(local_counter,el)=non_dir_bd_map_x(mapg(kl,el))
! 	schur_bd_local_to_boundary_map(Nt4+local_counter,el)=non_dir_bd_map_y(mapg(kl,el))
! 	schur_bd_local_map(kl,el)=local_counter
! 	
!       ENDDO
    ENDDO

  END SUBROUTINE create_schur_maps
  
  SUBROUTINE assemble_local_to_global_boundary_vec(vec_in,vec_out)
    IMPLICIT NONE
    INTEGER :: i,el,globi
    DOUBLE PRECISION :: temp(1:bd_numrows_x+bd_numrows_y+1),&
			vec_in(1:8*N,numelm),&
			vec_out(1:bd_numrows_x+bd_numrows_y)
    
    vec_out=0d0
    DO el=1,numelm
    temp=0d0
      DO i=1,8*N
	globi = schur_bd_local_to_boundary_map(i,el)
	temp(globi) = temp(globi) + vec_in(i,el)
      ENDDO
      vec_out = vec_out + temp(1:bd_numrows_x+bd_numrows_y)
    ENDDO
    
    
  END SUBROUTINE assemble_local_to_global_boundary_vec
  
  SUBROUTINE assemble_local_to_global_internal_vec(vec_in,vec_out)
    IMPLICIT NONE
    INTEGER :: i,el,globi
    DOUBLE PRECISION :: temp(1:3*npint),&
			vec_in(1:threeNM1SQ,numelm),&
			vec_out(1:3*npint)
    temp=0d0
    DO el=1,numelm
      DO i=1,NM1SQ
	globi = mapg(i,el)-npedg
	temp(globi) = temp(globi) + vec_in(i,el)
	temp(npint+globi) = temp(npint+globi) + vec_in(NM1SQ+i,el)
	temp(2*npint+globi) = temp(2*npint+globi) + vec_in(2*NM1SQ+i,el)
      ENDDO
    ENDDO
    vec_out=temp(1:bd_numrows_x+bd_numrows_y)
    
  END SUBROUTINE assemble_local_to_global_internal_vec 
  
  SUBROUTINE assemble_local_to_global_boundary_matrix(mat_in,mat_out)
    IMPLICIT NONE
    INTEGER :: i,j,el,globi,globj
    DOUBLE PRECISION :: temp(1:bd_numrows_x+bd_numrows_y+1,1:bd_numrows_x+bd_numrows_y+1),&
			mat_in(1:8*N,1:8*N,numelm),&
			mat_out(1:bd_numrows_x+bd_numrows_y,1:bd_numrows_x+bd_numrows_y)
    
    mat_out=0d0
    DO el=1,numelm
      temp=0d0
      DO i=1,8*N
	globi = schur_bd_local_to_boundary_map(i,el)
	DO j=1,8*N
	  globj = schur_bd_local_to_boundary_map(j,el)
	  temp(globi,globj) = temp(globi,globj) + mat_in(i,j,el)
	ENDDO
      ENDDO
      mat_out = mat_out + temp(1:bd_numrows_x+bd_numrows_y,1:bd_numrows_x+bd_numrows_y)
    ENDDO

  END SUBROUTINE assemble_local_to_global_boundary_matrix  
  
  SUBROUTINE scatter_global_boundary_to_local_vector(vec_in,vec_out)
    IMPLICIT NONE
    INTEGER :: i,el
    DOUBLE PRECISION :: temp(1:bd_numrows_x+bd_numrows_y+1),&
			vec_out(1:8*N,numelm),&
			vec_in(1:bd_numrows_x+bd_numrows_y)
    
    temp(1:bd_numrows_x+bd_numrows_y)=vec_in
    temp(bd_numrows_x+bd_numrows_y+1)=0d0
    DO el=1,numelm
      DO i=1,8*N
	vec_out(i,el)=temp(schur_bd_local_to_boundary_map(i,el))
      ENDDO
    ENDDO
  
  END SUBROUTINE scatter_global_boundary_to_local_vector
  
  SUBROUTINE element_wise_mv_mult(trans,matrix,mat_rows,mat_cols,alpha,vec_in,vec_rows,beta,vec_out,vec_cols)
! Perform an element-wise matrix-vector multiplication using DGEMV from LAPACK
! the result is vec_out = alpha*matrix*vec_in + beta*vec_out
! ie, we may add/subtract the result from output vector if it is non-zero.
!
! trans controls whether we perform transpose(matrix)*vec_in or not.
! - 'N' is normal
! - 'T' is transpose
! NOTE:
!      - If we use normal then we need to match mat_cols to vec_rows and mat_rows to vec_cols
!      - If we use a transpose then we need to match mat_rows to vec_rows and mat_cols to vec_cols
    IMPLICIT NONE
    
    CHARACTER(1) :: trans
    INTEGER :: i,j,el, mat_rows, mat_cols, vec_rows, vec_cols
    DOUBLE PRECISION :: alpha, beta, &
			matrix(mat_rows,mat_cols,numelm),&
			vec_in(vec_rows,numelm),&
			vec_out(vec_cols,numelm),&
			tempmat(mat_rows,mat_cols),&
			tempvec_in(vec_rows),&
			tempvec_out(vec_cols)
    
    IF (trans.EQ.'N') THEN
      print*,'hi'
      IF (mat_cols.NE.vec_rows) THEN
	write(*,*) 'Error in element_wise_mv_mult: Input dimension mismatch!'
	STOP
      ENDIF
      IF (mat_rows.NE.vec_cols) THEN
	write(*,*) 'Error in element_wise_mv_mult: Output dimension mismatch!'
	STOP
      ENDIF
    ELSEIF (trans.EQ.'T') THEN
      IF (mat_rows.NE.vec_rows) THEN
	write(*,*) 'Error in element_wise_mv_mult: Input dimension mismatch!'
	STOP
      ENDIF
      IF (mat_cols.NE.vec_cols) THEN
	write(*,*) 'Error in element_wise_mv_mult: Output dimension mismatch!'
	STOP
      ENDIF
    ELSE
      write(*,*) 'Error in element_wise_mv_mult: Invalid matrix operator N or T only!!'
      STOP
    ENDIF
    
!     open(123,FILE='A.dat')
!     open(1234,FILE='x.dat')
!     open(12345,FILE='b.dat')
    DO el=1,numelm
      DO i=1,mat_rows
	DO j=1,mat_cols
	  tempmat(i,j) = matrix(i,j,el)
	ENDDO
! 	write(123,*) (tempmat(i,j),j=1,mat_cols)
      ENDDO
      DO i=1,vec_rows
	tempvec_in(i) = vec_in(i,el)
! 	write(1234,*) tempvec_in(i)
      ENDDO
      DO i=1,vec_cols
	tempvec_out = vec_out(i,el)
      ENDDO
      
      CALL DGEMV(trans, mat_rows, mat_cols, alpha, tempmat, mat_rows, &
			tempvec_in, 1, beta, tempvec_out, 1)
      DO i=1,vec_cols
	vec_out(i,el) = tempvec_out(i)
! 	write(12345,*) tempvec_out(i)
      ENDDO
            
    ENDDO
!     close(123)
!     close(1234)
!     close(12345)
  END SUBROUTINE element_wise_mv_mult  
  
  SUBROUTINE element_wise_mm_mult(trans1,trans2,matrix1,rows1,cols1,matrix2,rows2,cols2,alpha,beta,matrix_out,rows_out,cols_out)
! Perform an element-wise matrix-matrix multiplication using DGEMM from LAPACK
! the result is matrix_out = alpha*op(matrix1)*op(matrix2) + beta*matrix(out)
! ie, we may add/subtract the result from output matrix if the input is non-zero.
!
! Input for the rows/cols is as normal for the input matrices. We will sort out the transpose stuff
! within the routine.
! i.e. matrix1 is of shape rows1 by cols1
!      matrix2 is of shape rows2 by cols2
!      matrix_out is of shape rows_out by cols_out
!
! trans controls the op(matrix):
! - 'N' is normal
! - 'T' is transpose
    IMPLICIT NONE
    
    CHARACTER(1) :: trans1,trans2
    INTEGER :: i,j,el,rows,cols,col1_row2,lda,ldb,rows1,rows2,&
		cols1,cols2,rows_out,cols_out
    DOUBLE PRECISION :: alpha, beta, &
			matrix1(rows1,cols1,numelm),&
			matrix2(rows2,cols2,numelm),&
			matrix_out(rows_out,cols_out,numelm),&
			tempmat1(rows1,cols1),&
			tempmat2(rows2,cols2),&
			tempmat_out(rows_out,cols_out)
    
    IF (trans1.eq.'N') THEN
      rows = rows1
      col1_row2 = cols1
    ELSE
      rows = cols1
      col1_row2 = rows1
    ENDIF
    IF (trans2.eq.'N') THEN
      cols = cols2
      IF (col1_row2.ne.rows2) THEN
	write(*,*) 'Error in element_wise_mm_mult: Input dimension mismatch!'
	STOP
      ENDIF
    ELSE
      cols = rows2
      IF (col1_row2.ne.cols2) THEN
	write(*,*) 'Error in element_wise_mm_mult: Input dimension mismatch!'
	STOP
      ENDIF
    ENDIF
    IF (rows_out.ne.rows.OR.cols_out.ne.cols) THEN
      write(*,*) 'Error in element_wise_mm_mult: Output dimension mismatch!'
      STOP
    ENDIF

    lda=rows1
    ldb=rows2
    
    open(123,FILE='X1.dat')
    open(1234,FILE='Y.dat')
    open(12345,FILE='Z.dat')
    DO el=1,numelm
    
      tempmat1 = matrix1(1:rows1,1:cols1,el)
      tempmat2 = matrix2(1:rows2,1:cols2,el)
      tempmat_out = 0d0! matrix_out(1:rows_out,1:cols_out,el)
      
      DO i=1,rows1
	write(123,*) (tempmat1(i,j),j=1,cols1)
      ENDDO
      DO i=1,rows2
	write(1234,*) (tempmat2(i,j),j=1,cols2)
      ENDDO
      
      

      CALL DGEMM(trans1, trans2, rows, cols, col1_row2, alpha, &
		  tempmat1, lda, tempmat2, ldb, &
		  beta, tempmat_out, rows_out)


      matrix_out(1:rows,1:cols,el) = tempmat_out
      DO i=1,rows_out
	write(12345,*) (tempmat_out(i,j),j=1,cols_out)
      ENDDO
      
    ENDDO
  close(123)
  close(1234)
  close(12345)
  END SUBROUTINE element_wise_mm_mult
  
!     SUBROUTINE element_wise_transM_M_mult(rows,cols,col1_row2,alpha,matrix1,matrix2,beta,matrix_out)
! ! Perform an element-wise matrix-matrix multiplication using DGEMM from LAPACK
! ! the result is matrix_out = alpha*transpose(matrix1)*matrix2 + beta*matrix(out)
! ! ie, we may add/subtract the result from output matrix if the input is non-zero.
! !
! ! Input for the rows/cols is as normal for the input matrices. We will sort out the transpose stuff
! ! within the routine.
! ! i.e. matrix1 is of shape rows by col1_row2
! !      matrix2 is of shape col1_row2 by cols
! !      matrix_out is of shape rows by cols
! 
!     IMPLICIT NONE
!     
!     CHARACTER(1) :: trans1,trans2
!     INTEGER :: i,j,el, rows,cols,col1_row2,lda,ldb
!     DOUBLE PRECISION :: alpha, beta, &
! 			matrix1(col1_row2,rows,numelm),&
! 			matrix2(col1_row2,cols,numelm),&
! 			matrix_out(rows,cols,numelm),&
! 			tempmat1(rows,col1_row2),&
! 			tempmat2(col1_row2,cols),&
! 			tempmat_out(rows,cols)
!     
!     
!     DO el=1,numelm
!     
!       tempmat1 = matrix1(1:rows,1:col1_row2,el)
!       tempmat2 = matrix2(1:col1_row2,1:cols,el)
!       tempmat_out = matrix_out(1:rows,1:cols,el)
!       
!       CALL DGEMM( 'T', 'N', rows, cols, col1_row2, alpha, &
! 		  tempmat1, lda, tempmat2, ldb, &
! 		  beta, tempmat_out, rows )
! 
!       print*,'1'
!       do i=1,col1_row2
! 	write(*,*) (tempmat1(j,i),j=1,rows)
!       enddo
!       print*,'2'
!       do i=1,col1_row2
! 	write(*,*) (tempmat2(i,j),j=1,cols)
!       enddo	
!       print*,'3'
!       do i=1,rows
! 	write(*,*) (tempmat_out(i,j),j=1,cols)
!       enddo
!       stop
!       matrix_out(1:rows,1:cols,el) = tempmat_out
!       
!     ENDDO
! 
!   END SUBROUTINE element_wise_transM_M_mult
  
  SUBROUTINE solve_all_internal(matrix,matrixfact,rhs,ipiv)
    IMPLICIT NONE
  
    INTEGER :: i,el,INFO,matrix_order,&
	       iwork(1:threeNM1SQ),&
	       ipiv(1:threeNM1SQ,numelm), &
	       tempipiv(threeNM1SQ)
    DOUBLE PRECISION :: ferr,berr,temp_const,&
			work(1:9*NM1SQ),&
			matrix(threeNM1SQ,threeNM1SQ,numelm),&
			matrixfact(threeNM1SQ,threeNM1SQ,numelm),&
			rhs(threeNM1SQ,numelm),&
			temp(threeNM1SQ,threeNM1SQ),&
			temprhs(threeNM1SQ),&
			temp2(threeNM1SQ,threeNM1SQ),&
			temprhs2(threeNM1SQ),&
			local_pressure_only(1:NM1SQ,numelm)
			
  
    matrix_order = threeNM1SQ
! matrix will now contain the lower cholesky factorisation.
    DO el=1,numelm
      temp=matrixfact(1:matrix_order,1:matrix_order,el)
      temp2=matrix(1:matrix_order,1:matrix_order,el)
      temprhs=rhs(1:matrix_order,el)
      
!       temprhs(2*NM1SQ+1)=0d0
      
      temprhs2=temprhs
      tempipiv=ipiv(1:matrix_order,el)

! Solve
!     CALL DPOTRS('L', matrix_order, 1, temp, matrix_order, temprhs, matrix_order, INFO)
      CALL DGETRS( 'N', matrix_order, 1, temp, matrix_order, tempipiv, temprhs, matrix_order, info )
      if (info.lt.0) print*,'info is',info
! Refine solution
!       CALL DGERFS( 'N', matrix_order, 1, temp2, matrix_order, temp, matrix_order, tempipiv, &
! 			temprhs2, matrix_order, temprhs, matrix_order, ferr, berr, work, iwork, info )

      if (info.lt.0) print*,'info is',info
!       temprhs(2*NM1SQ+1)=0d0    
      rhs(1:matrix_order,el)=temprhs
    ENDDO
    
!     local_pressure_only = rhs(2*NM1SQ+1:matrix_order,1:numelm)
!     CALL normalise_local_pressure_soln(local_pressure_only)
!     CALL integrate_pressure_on_domain(local_pressure_only,temp_const)
! !     local_pressure_only=local_pressure_only-temp_const/area_of_domain
! 
!     rhs(2*NM1SQ+1:matrix_order,1:numelm)=local_pressure_only-temp_const/area_of_domain
  
  END SUBROUTINE solve_all_internal
  
  SUBROUTINE solve_all_int_mult_cross_matrix(matrix,matrixfact,matrixrhs,ipiv)
    IMPLICIT NONE
  
    INTEGER :: i,j,el,INFO,matrix_order,&
	       ipiv(1:threeNM1SQ,numelm),&
	       tempipiv(1:threeNM1SQ),&
	       iwork(1:threeNM1SQ)

    DOUBLE PRECISION :: temp_const,&
			ferr(1:8*N),&
			berr(1:8*N),&
			matrix(threeNM1SQ,threeNM1SQ,numelm),&
			matrixfact(threeNM1SQ,threeNM1SQ,numelm),&
			work(1:9*NM1SQ),&
			matrixrhs(threeNM1SQ,8*N,numelm),&
			tempmatfact(threeNM1SQ,threeNM1SQ),&
			tempsol(threeNM1SQ,8*N),&
			tempmat(threeNM1SQ,threeNM1SQ),&
			temprhs(threeNM1SQ,8*N),&
			local_pressure_only(1:NM1SQ,numelm)

    
    matrix_order = threeNM1SQ

    DO el=1,numelm
      tempmatfact = matrixfact(1:matrix_order,1:matrix_order,el)
            
      tempsol = matrixrhs(1:matrix_order,1:8*N,el)
      
!       tempsol(2*NM1SQ+1,1:8*N)=0d0
      
      tempmat = matrix(1:matrix_order,1:matrix_order,el)
!       tempmat(1:matrix_order,2*NM1SQ+1)=0d0
!       tempmat(2*NM1SQ+1,1:matrix_order)=0d0
!       tempmat(2*NM1SQ+1,2*NM1SQ+1)=1d0
      
      temprhs = tempsol
      tempipiv = ipiv(1:matrix_order,el)

! Solve:
!     CALL DPOTRS( 'L', matrix_order, 8*N, temp, matrix_order, temprhs, matrix_order, INFO )
      CALL DGETRS( 'N', matrix_order, 8*N, tempmatfact, matrix_order, tempipiv, tempsol, matrix_order, INFO )
      if (info.lt.0) print*,'info is',info
      
! Refine:
!       CALL DGERFS( 'N', matrix_order, 8*N, tempmat, matrix_order, tempmatfact, matrix_order, tempipiv,&
! 		   temprhs, matrix_order, tempsol, matrix_order, ferr, berr, work, iwork, info )
      if (info.lt.0) print*,'info is',info
!       tempsol(2*NM1SQ+1,1:8*N) = 0d0
      matrixrhs(1:matrix_order,1:8*N,el) = tempsol
    ENDDO
    
!     DO j=1,8*N
!       local_pressure_only = matrixrhs(2*NM1SQ+1:matrix_order,j,1:numelm)
!       
!       CALL normalise_local_pressure_soln(local_pressure_only)
!       CALL integrate_pressure_on_domain(local_pressure_only,temp_const)
! !       local_pressure_only=local_pressure_only-temp_const/area_of_domain
!       matrixrhs(2*NM1SQ+1:matrix_order,j,1:numelm) = local_pressure_only!-temp_const/area_of_domain
!       
!     ENDDO

  END SUBROUTINE solve_all_int_mult_cross_matrix
  
  
  SUBROUTINE fact_all_internal(matrix,ipiv)
    IMPLICIT NONE
  
    INTEGER :: info,el,matrix_order,&
		ipiv(1:threeNM1SQ,numelm), &
		tempipiv(threeNM1SQ)
    DOUBLE PRECISION :: matrix(threeNM1SQ,threeNM1SQ,numelm), &
			temp(threeNM1SQ,threeNM1SQ)
  
    matrix_order = threeNM1SQ
! matrix will now contain the lower cholesky factorisation.
    DO el=1,numelm
      temp=matrix(1:matrix_order,1:matrix_order,el)
!       temp(1:matrix_order,2*NM1SQ+1)=0d0
!       temp(2*NM1SQ+1,1:matrix_order)=0d0
!       temp(2*NM1SQ+1,2*NM1SQ+1)=1d0

!     CALL DPOTRF( 'L', matrix_order, temp, matrix_order, INFO)
      CALL DGETRF( matrix_order, matrix_order, temp, matrix_order, tempipiv, info )
      
      if (info.lt.0) print*,'info is',info
      
      matrix(1:matrix_order,1:matrix_order,el)=temp
      ipiv(1:matrix_order,el)=tempipiv
    ENDDO
  
  END SUBROUTINE fact_all_internal
  
  SUBROUTINE build_local_internal_matrices(local_matrix)
  
    IMPLICIT NONE
    
    INTEGER :: i,j,k,l,ij,kl,intij,intkl,el,NM14
    DOUBLE PRECISION :: const
    DOUBLE PRECISION, DIMENSION(1:threeNM1SQ,1:threeNM1SQ,numelm) :: local_matrix
    
! store (N-1)^2 and (N-1)^4 for later - saves recomputing
    NM14 = NM1SQ**2
! initialise storage matrix
    local_matrix=0d0
    const = 3d0*Re/(2d0*deltat)

! Build the local "stokes" matrix (including the additional terms for OIFS, etc)
! will need to modify for other schemes, eg the crank nicolson thingys!
    DO el=1,numelm
      DO l=1,NM1
      DO k=1,NM1
	kl = k+l*NP1
	intkl = k+(l-1)*NM1
	DO j=1,NM1
	DO i=1,NM1
	  ij = i+j*NP1
	  intij = i+(j-1)*NM1
	  
	  local_matrix(intij,intkl,el) = A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
	  local_matrix(intij+NM1SQ,intkl+NM1SQ,el) = A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
	  local_matrix(intij+2*NM1SQ,intkl,el) = B_x(intij,kl,el)
	  local_matrix(intij+2*NM1SQ,intkl+NM1SQ,el) = B_y(intij,kl,el)
	  local_matrix(intij,intkl+2*NM1SQ,el) = B_x(intkl,ij,el)
	  local_matrix(intij+NM1SQ,intkl+2*NM1SQ,el) = B_y(intkl,ij,el)
	  
	ENDDO
	ENDDO
      ENDDO
      ENDDO
! Zero the row/column of the first local pressure node and its test function, this is the row/column 2*NM1SQ+1
!      local_matrix(1+2*NM1SQ,1:3*NM1SQ,el) = 0d0
!      local_matrix(1:threeNM1SQ,1+2*NM1SQ,el) = 0d0
! ! Set the diagonal entry to be 1
!      local_matrix(1+2*NM1SQ,1+2*NM1SQ,el) = 1d0

      
      
    ENDDO
  
  END SUBROUTINE build_local_internal_matrices
  
  SUBROUTINE build_internal_rhs(local_rhs)
! Builds the local RHS vector for internal nodes only
! This is for our OIFS scheme. Other schemes will have a different RHS.
  
    IMPLICIT NONE
    
    INTEGER :: i,j,intij,ij,el
    DOUBLE PRECISION, DIMENSION(1:threeNM1SQ,numelm) :: local_rhs
    
    DO el=1,numelm
      DO j=1,NM1
	DO i=1,NM1
	  ij=i+j*NP1
	  intij=i+(j-1)*NM1
	  
	  local_rhs(intij,el) = f_x(mapg(ij,el))
	  local_rhs(NM1SQ+intij,el) = f_y(mapg(ij,el))
	  local_rhs(2*NM1SQ+intij,el) = g(mapg(ij,el)-npedg)
	ENDDO
      ENDDO
!       local_rhs(2*NM1SQ+1,el) = 0d0
    ENDDO
    
  
  END SUBROUTINE build_internal_rhs
  
  SUBROUTINE build_local_crossover_matrices(local_matrix)
  
    IMPLICIT NONE
    
    INTEGER :: i,j,k,l,ij,kl,intij,intkl,el
    DOUBLE PRECISION :: const
    DOUBLE PRECISION, DIMENSION(1:threeNM1SQ,1:8*N,numelm) :: local_matrix
    
! initialise storage matrix
    local_matrix=0d0
    const = 3d0*Re/(2d0*deltat)

! Build the local "stokes" matrix (including the additional terms for OIFS, etc)
! will need to modify for other schemes, eg the crank nicolson thingys!
    DO el=1,numelm
      DO kl=0,NP1SQ-1
	IF (mapg(kl,el).LE.npedg) THEN
! 	  IF (.not.bdflag(1,mapg(kl,el))) THEN
	    DO j=1,NM1
	      DO i=1,NM1
		ij = i+j*NP1
		intij = i+(j-1)*NM1
	      
		local_matrix(intij,schur_bd_local_map(kl,el),el) = A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
		local_matrix(intij+2*NM1SQ,schur_bd_local_map(kl,el),el) = B_x(intij,kl,el)
		
		local_matrix(intij+NM1SQ,4*N+schur_bd_local_map(kl,el),el) = A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
	        local_matrix(intij+2*NM1SQ,4*N+schur_bd_local_map(kl,el),el) = B_y(intij,kl,el)
	      ENDDO
	    ENDDO
! 	  ENDIF
! 	  IF (.not.bdflag(2,mapg(kl,el))) THEN
! 	    DO j=1,NM1
! 	      DO i=1,NM1
! 		ij = i+j*NP1
! 		intij = i+(j-1)*NM1
! 		
! 
! 	      ENDDO
! 	    ENDDO
! 	  ENDIF
	ENDIF
      ENDDO
! Zero the row for the first pressurenode 1+2*NM1SQ:
!       local_matrix(1+2*NM1SQ,1:8*N,el)=0d0
    ENDDO

  END SUBROUTINE build_local_crossover_matrices
  
  SUBROUTINE build_crossover_matrix(crossover_matrix)
! builds the global crossover matrix for the schur complement decomposition
! For the columns we only want the entries associated with the edge points of each element, so we go along
! each edge collecting the entries for kl.
! For the rows, we want only the internal entires, so we restrict to internal parts of each element (1 to N-1)
    IMPLICIT NONE
    INTEGER :: el,i,j,ij,intij,k,l,kl,NtNP1
    DOUBLE PRECISION :: const,&
			crossover_matrix(1:3*npint,1:bd_numrows_x+bd_numrows_y+1),&
			temp(1:3*npint,1:bd_numrows_x+bd_numrows_y+1)

! a few useful numbers to save lots of wasted operations.
    NtNP1=N*NP1 !N times N plus 1
    const = 3d0*Re/(2d0*deltat)
    
    crossover_matrix = 0d0
    DO el=1,numelm
      DO kl=0,NP1SQ-1
	IF (mapg(kl,el).LE.npedg) THEN
! loop over internal nodes
	  DO j=1,NM1
	    DO i=1,NM1
	      ij=i+j*NP1
	      intij=i+(j-1)*NM1
	      

	  
	      crossover_matrix(mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) = &
						crossover_matrix(mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) &
						+ A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
						
	      crossover_matrix(npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) = &
						crossover_matrix(npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) &
						+ A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
						
	      crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) = &
						crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el)))&
						+ B_x(intij,kl,el)
	    
	      crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) = &
						crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el)))&
						+ B_y(intij,kl,el)
	    
	    ENDDO
	  ENDDO
	ENDIF					
      ENDDO
    ENDDO	 

! loop over internal nodes
!       DO j=1,NM1
! 	DO i=1,NM1
! 	  ij=i+j*NP1
! 	  intij=i+(j-1)*NM1
! Now loop over boundary nodes - ie around the edge of the element.
! 	  k=0
! 	  DO l=0,N
! 	    kl=l*NP1
! 	    crossover_matrix(mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) = &
! 						crossover_matrix(mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	    crossover_matrix(npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) = &
! 						crossover_matrix(npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 						
! 	    crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) = &
! 						crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) + B_x(intij,kl,el)
! 	    
! 	    crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) = &
! 						crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) + B_y(intij,kl,el)
! 	  ENDDO
! !	  k=N
! 	  DO l=0,N
! 	    kl=N+l*NP1
! 	    crossover_matrix(mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) = &
! 						crossover_matrix(mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	    crossover_matrix(npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) = &
! 						crossover_matrix(npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 						
! 	    crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) = &
! 						crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) + B_x(intij,kl,el)
! 	    
! 	    crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) = &
! 						crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) + B_y(intij,kl,el)
! 	  ENDDO
! !	  l=0 | note we've already done k=0,l=0 and k=N,l=0 above
! 	  DO kl=1,NM1
! 	    
! 	    crossover_matrix(mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) = &
! 						crossover_matrix(mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	    crossover_matrix(npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) = &
! 						crossover_matrix(npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 						
! 	    crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) = &
! 						crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) + B_x(intij,kl,el)
! 	    
! 	    crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) = &
! 						crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) + B_y(intij,kl,el)
! 	  ENDDO
! !	  l=N | note we've already done k=0,l=N and k=N,l=N above.
! 	  DO k=1,NM1
! 	    kl=k+NtNP1
! 	    crossover_matrix(mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) = &
! 						crossover_matrix(mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	    crossover_matrix(npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) = &
! 						crossover_matrix(npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 						
! 	    crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) = &
! 						crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_x(mapg(kl,el))) + B_x(intij,kl,el)
! 	    
! 	    crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) = &
! 						crossover_matrix(2*npint + mapg(ij,el)-npedg,non_dir_bd_map_y(mapg(kl,el))) + B_y(intij,kl,el)
! 	  ENDDO
! 	ENDDO
!       ENDDO
!     ENDDO
	
  END SUBROUTINE build_crossover_matrix
  

  
  SUBROUTINE build_boundary_matrix(boundary_matrix)
! builds the global boundary matrix for the schur complement decomposition
! We only want the entries associated with the edge points of each element, so we go along
! each edge collecting the entries for both ij and kl.
!
! Note that the boundary matrix only contains the Laplace matrix (made up of A_x and A_y on the diagonal)

    IMPLICIT NONE
    INTEGER :: el,i,j,ij,k,l,kl,NtNP1
    DOUBLE PRECISION :: const,&
			temp(1:bd_numrows_x+bd_numrows_y+1,1:bd_numrows_x+bd_numrows_y+1),&
			boundary_matrix(1:bd_numrows_x+bd_numrows_y,1:bd_numrows_x+bd_numrows_y)
			

! a few useful numbers to save lots of wasted operations.
    NtNP1=N*NP1 !N times N plus 1
    const = 3d0*Re/(2d0*deltat)
    
    boundary_matrix=0d0
    DO el=1,numelm
    
      DO ij=0,NP1SQ-1
	IF (mapg(ij,el).LE.npedg) THEN
	  DO kl=0,NP1SQ-1
	    IF (mapg(kl,el).LE.npedg) THEN
	      temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
						
	      temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
	    ENDIF
	  ENDDO
	ENDIF
      ENDDO
! !   i=0      
!       DO j=0,N
! 	ij=j*NP1
! !     k=0
! 	DO l=0,N
! 	  kl=l*NP1
! 	  temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	  temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 	ENDDO
! !     k=N
! 	DO l=0,N
! 	  kl=N+l*NP1
! 	  temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	  temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 	ENDDO
! !     l=0 | note: we've already done k=0,l=0 and k=N,l=0 above.
! 	DO kl=1,NM1
! 	  temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	  temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 	ENDDO
! !     l=N | note: we've already done k=0,l=N and k=N,l=N above.
! 	DO k=1,NM1
! 	  kl=k+NtNP1
! 	  temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	  temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 	ENDDO
!       ENDDO
!       
! !   i=N     
!       DO j=0,N
! 	ij=N+j*NP1
! !     k=0
! 	DO l=0,N
! 	  kl=l*NP1
! 	  temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	  temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 	ENDDO
! !     k=N
! 	DO l=0,N
! 	  kl=N+l*NP1
! 	  temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	  temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 	ENDDO
! !     l=0 | note: we've already done k=0,l=0 and k=N,l=0 above.
! 	DO kl=1,NM1
! 	  temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	  temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 	ENDDO
! !     l=N | note: we've already done k=0,l=N and k=N,l=N above.
! 	DO k=1,NM1
! 	  kl=k+NtNP1
! 	  temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	  temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 	ENDDO
!       ENDDO
!       
!       
! !   j=0 | already done i=0,j=0 and i=N,j=0 above.  
!       DO ij=1,NM1
! 	
! !     k=0
! 	DO l=0,N
! 	  kl=l*NP1
! 	  temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	  temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 	ENDDO
! !     k=N
! 	DO l=0,N
! 	  kl=N+l*NP1
! 	  temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	  temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 	ENDDO
! !     l=0 | note: we've already done k=0,l=0 and k=N,l=0 above.
! 	DO kl=1,NM1
! 	  temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	  temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 	ENDDO
! !     l=N | note: we've already done k=0,l=N and k=N,l=N above.
! 	DO k=1,NM1
! 	  kl=k+NtNP1
! 	  temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	  temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 	ENDDO
!       ENDDO
!       
! !   j=N | already done i=0,j=N and i=N,j=N above.  
!       DO i=1,NM1
! 	ij=i+NtNP1
! !     k=0
! 	DO l=0,N
! 	  kl=l*NP1
! 	  temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	  temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 	ENDDO
! !     k=N
! 	DO l=0,N
! 	  kl=N+l*NP1
! 	  temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	  temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 	ENDDO
! !     l=0 | note: we've already done k=0,l=0 and k=N,l=0 above.
! 	DO kl=1,NM1
! 	  temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	  temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 	ENDDO
! !     l=N | note: we've already done k=0,l=N and k=N,l=N above.
! 	DO k=1,NM1
! 	  kl=k+NtNP1
! 	  temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_x(mapg(ij,el)),non_dir_bd_map_x(mapg(kl,el))) + A_x(ij,kl,el) + const*Mv_x(ij,ij,el)
! 						
! 	  temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) = &
! 			temp(non_dir_bd_map_y(mapg(ij,el)),non_dir_bd_map_y(mapg(kl,el))) + A_y(ij,kl,el) + const*Mv_y(ij,ij,el)
! 	ENDDO
!       ENDDO
      
    ENDDO
    
    DO i=1,bd_numrows_x+bd_numrows_y
      DO j=1,bd_numrows_x+bd_numrows_y
	boundary_matrix(i,j)=temp(i,j)
      ENDDO
    ENDDO
    
!     DO i=1,npedg
!       IF (bdflag(1,i)) THEN
! 	boundary_matrix(i,1:npedg)=0d0
! 	boundary_matrix(1:npedg,i)=0d0
! 	boundary_matrix(i,i)=1d0
!       ENDIF
!       IF (bdflag(2,i)) THEN
! 	boundary_matrix(npedg+i,npedg+1:2*npedg)=0d0
! 	boundary_matrix(npedg+1:2*npedg,npedg+i)=0d0
! 	boundary_matrix(npedg+i,npedg+i)=1d0
!       ENDIF
!     ENDDO
    
  END SUBROUTINE build_boundary_matrix
  
  SUBROUTINE build_boundary_rhs(boundary_rhs)
! Builds the global rhs vector. Any dirichlet nodes are mapped to the last entry by the schur complement map
! These will not be included when solving the system.
! This is for our OIFS scheme. Other schemes will have a different RHS.
  
    IMPLICIT NONE
    
    INTEGER :: i
    DOUBLE PRECISION :: temp(1:bd_numrows_x+bd_numrows_y+1),&
			boundary_rhs(1:bd_numrows_x+bd_numrows_y)
    
    DO i=1,npedg
      temp(non_dir_bd_map_x(i)) = f_x(i)
      temp(non_dir_bd_map_y(i)) = f_y(i)
    ENDDO
    DO i=1,bd_numrows_x+bd_numrows_y
      boundary_rhs(i)=temp(i)
    ENDDO
    
  END SUBROUTINE build_boundary_rhs
  
  SUBROUTINE reduce_global_crossover_to_locals(matrix_in,matrix_out)
! Takes a crossover matrix and splits the rows into seperate elements.
! Does not change the columns.
! The idea is that this will reduce the size of the matrix when multiplying with an internal vector (which may be stored element-wise too)
    IMPLICIT NONE
    INTEGER :: i,j,ij,intij,k,l,kl,intkl,el
    DOUBLE PRECISION :: matrix_in(1:3*npint,1:bd_numrows_x+bd_numrows_y),&
			matrix_out(1:threeNM1SQ,1:bd_numrows_x+bd_numrows_y,numelm)
    
    matrix_out=0d0
    DO el=1,numelm
      DO i=1,NM1
	DO j=1,NM1
	  ij=i+j*NP1
	  intij=i+(j-1)*NM1
	  DO kl=1,bd_numrows_x+bd_numrows_y
  
	    matrix_out(intij,kl,el) = matrix_in(mapg(ij,el),kl)
	    matrix_out(NM1SQ + intij,kl,el) = matrix_in(mapg(ij,el),kl)
	    matrix_out(2*NM1SQ + intij,kl,el) = matrix_in(mapg(ij,el),kl)

	  ENDDO
	  
	ENDDO
      ENDDO
    ENDDO
  
  END SUBROUTINE reduce_global_crossover_to_locals
  
  
  SUBROUTINE Mc_trans_mult_vec(matrix,vec,vec_out)
    IMPLICIT NONE
    INTEGER :: el,LDA,cols
    DOUBLE PRECISION :: matrix(1:threeNM1SQ,1:bd_numrows_x+bd_numrows_y,numelm),&
			vec(1:threeNM1SQ,numelm),&
			vec_out(1:bd_numrows_x+bd_numrows_y)
! 			y(1:bd_numrows_x+bd_numrows_y)
    
    vec_out=0d0
    LDA=threeNM1SQ
    cols=bd_numrows_x+bd_numrows_y

    DO el=1,numelm
      CALL DGEMV('T', LDA, cols, 1d0, matrix, LDA, vec, 1, 1d0, vec_out,1) ! Adds result to vec_out on each call.
    ENDDO
  
  END SUBROUTINE Mc_trans_mult_vec
  
  SUBROUTINE normalise_local_pressure_soln(vec)
    INTEGER :: els,el
    DOUBLE PRECISION :: p1,p2,local_pressure_const,&
			vec(1:NM1SQ,numelm), &
			tempvec_to(1:NM1SQ), &
			tempvec_from(1:NM1SQ)
    
! Sweep over order of elements to normalise
    DO els=1,numelm-1
! Pick the next current element in line:
      el=normalise_pressure_element_order(els)
      tempvec_to = vec(1:NM1SQ,el)
      tempvec_from = vec(1:NM1SQ,normalise_from_el(el))
! extrapolate pressure value at the current element's local node
      CALL extrapolate_pressure_in_element(tempvec_to,el,&
		    normalise_to_local_vertex(el,1),normalise_to_local_vertex(el,2),p1)
! extrapolate pressure value at the connected element's local node - this is the "correct" value
      CALL extrapolate_pressure_in_element(tempvec_from,normalise_from_el(el),&
		    normalise_from_local_vertex(el,1),normalise_from_local_vertex(el,2),p2)
      
      local_pressure_const = p2-p1
! Add constant to all internal local nodes in the current element.
      vec(1:NM1SQ,el) = tempvec_to + local_pressure_const
      
    ENDDO

  END SUBROUTINE normalise_local_pressure_soln
  
  
  SUBROUTINE setup_pressure_normalisation_connectivity
    IMPLICIT NONE
    LOGICAL :: done_element(1:numelm)
    INTEGER :: el,i,j,k,l,order_counter,curr_el,&
		vertex_to_local_sep(0:5,1:2),&
		vertex_to_local(0:5)

   
    done_element(1)=.true.
    done_element(2:numelm)=.false.
    
    vertex_to_local_sep(0,1)=0
    vertex_to_local_sep(0,2)=N
    vertex_to_local_sep(1,1)=0
    vertex_to_local_sep(1,2)=0
    vertex_to_local_sep(2,1)=N
    vertex_to_local_sep(2,2)=0
    vertex_to_local_sep(3,1)=N
    vertex_to_local_sep(3,2)=N
    vertex_to_local_sep(4,1)=0
    vertex_to_local_sep(4,2)=N
    vertex_to_local_sep(5,1)=0 ! repeat back to 1 - makes sense given the routine below.
    vertex_to_local_sep(5,2)=0
    
    vertex_to_local(0)=N*NP1
    vertex_to_local(1)=0    
    vertex_to_local(2)=N
    vertex_to_local(3)=NP1SQ -1
    vertex_to_local(4)=N*NP1 ! repeat back to 1 - makes sense given the routine below.
    vertex_to_local(5)=0
    
    
! Assume that element 1 is normalised.
! Go through each element and note which it is connected to.
! We will choose elements normalise using the current element, and log which local node in each is required.
! We log the order in which to perform normalisation too.
! i.e. we use element 1 to normalise all elements it is connected to, then we use those elements to normalise
!      elements that they are connected to, and so on.
! We also make sure that the same element is not normalised twice.

! conelm(i,j) is the element which element i, face j is touching
! conedg(i,j) is the face which element i, face j is touching on
!             the element identified in conelm(i,j)

    order_counter=0
    DO el=1,numelm
      DO i=1,4
	curr_el=conelm(el,i)
	IF (curr_el.ne.0) THEN    
	  IF (.not.done_element(curr_el)) THEN

	    normalise_from_el(curr_el) = el
	    IF (mapg(vertex_to_local(i),el).eq.mapg(vertex_to_local(conedg(el,i)),curr_el)) THEN
	      normalise_from_local_vertex(curr_el,1) = vertex_to_local_sep(i,1)
	      normalise_from_local_vertex(curr_el,2) = vertex_to_local_sep(i,2)
	      normalise_to_local_vertex(curr_el,1) = vertex_to_local_sep(conedg(el,i),1)
	      normalise_to_local_vertex(curr_el,2) = vertex_to_local_sep(conedg(el,i),2)
	    ELSE! (mapg(vertex_to_local(i+1),curr_el).eq.mapg(vertex_to_local(conedg(el,i)),curr_el)) THEN
	      k=i+1
	      normalise_from_local_vertex(curr_el,1) = vertex_to_local_sep(k,1)
	      normalise_from_local_vertex(curr_el,2) = vertex_to_local_sep(k,2)
	      normalise_to_local_vertex(curr_el,1) = vertex_to_local_sep(conedg(el,i),1)
	      normalise_to_local_vertex(curr_el,2) = vertex_to_local_sep(conedg(el,i),2)
	    ENDIF
	      
	    done_element(curr_el)=.true.
	    order_counter=order_counter+1
	    normalise_pressure_element_order(order_counter)=curr_el
	  ENDIF	  
	ENDIF
      ENDDO
    ENDDO
!     DO el=1,numelm-1
!     curr_el=normalise_pressure_element_order(el)
!     print*,curr_el,normalise_to_local_vertex(curr_el,1),&
! 	    normalise_to_local_vertex(curr_el,2),normalise_from_el(curr_el),&
! 	    normalise_from_local_vertex(curr_el,1),normalise_from_local_vertex(curr_el,2)
! 	    enddo
! 	    
! 	    stop
  
  END SUBROUTINE setup_pressure_normalisation_connectivity
  
  SUBROUTINE integrate_pressure_on_domain(p_vec,out)
    INTEGER :: el,i,j,intij
    DOUBLE PRECISION :: out,temp,&
			p_vec(1:NM1SQ,numelm)
    
    out=0d0
    DO el=1,numelm
    temp=0d0
      DO i=1,NM1
	DO j=1,NM1
	  intij=i+(j-1)*NM1
	  temp = temp + p_vec(intij,el)*( &
		  w(i)*(abs(jac(i,j,el))*w(j) + evalh(j,0)*abs(jac(i,0,el))*w(0) + evalh(j,N)*abs(jac(i,N,el))*w(N)) + &
       w(0)*evalh(i,0)*(abs(jac(0,j,el))*w(j) + evalh(j,0)*abs(jac(0,0,el))*w(0) + evalh(j,N)*abs(jac(0,N,el))*w(N)) + &
       w(N)*evalh(i,N)*(abs(jac(N,j,el))*w(j) + evalh(j,0)*abs(jac(N,0,el))*w(0) + evalh(j,N)*abs(jac(N,N,el))*w(N))   &
					)
	ENDDO
      ENDDO
      out=out+temp
    ENDDO
    out=temp
  END SUBROUTINE integrate_pressure_on_domain
  
  SUBROUTINE extrapolate_pressure_in_element(p_vec,el,point_i,point_j,out)
! Returns the value of pressure on the boundary point, point_ij, of an element for a local pressure-sized vector.
    IMPLICIT NONE
    INTEGER :: point_i,point_j,el,i,j,ij
    DOUBLE PRECISION :: out,temp,&
			p_vec(1:NM1SQ,numelm)

    temp=0d0
    DO i=1,NM1
      DO j=1,NM1
	ij=i+(j-1)*NM1
	temp = temp + p_vec(ij,el)*evalh(i,point_i)*evalh(j,point_j)
      ENDDO
    ENDDO
    out=temp
    
  END SUBROUTINE extrapolate_pressure_in_element
  
  SUBROUTINE MINRES_internal_only(Mat_in,X_in,RHS_in,diag_precon,minresconv)
! solving Mx = f where M is a matrix decomposed into the stokes problem with A(laplacian) and B(divergence),
! x is decomposed into velocity and pressure. f is our vector of known terms.
! System is as follows:
!			[ A_x  0  B_x' ] [ U_in ]   [ RHSx_in ]
!			[  0  A_y B_y' ] [ V_in ] = [ RHSy_in ]
!			[ B_x B_y  C   ] [ p_in ]   [ RHSp_in ]
!
  
    IMPLICIT NONE
    INTEGER :: i,j,k,l,counter,stag,restartcount,minresconv,el,ij,kl,intij
    
    DOUBLE PRECISION :: eta, sold, snew, cold, cnew, alpha0,alpha1,alpha2,alpha3,&
		      gammaold,gammanew,delta,temp,normr


    DOUBLE PRECISION, DIMENSION(1:NM1SQ) :: RHSx,RHSy,RHSp,&
					    v1xold,v1yold,v1xnew,v1ynew,&
					    z1xold,z1yold,z1xnew,z1ynew,&
					    w1xold,w1yold,w1xnew,w1ynew,&
					    tempx1,tempy1,Ax,Ay,x1,y1,&
					    xstore1,ystore1,&
					    v2old,v2new,z2old,z2new,w2old,w2new,temp2,x2,y2,Bx,By,Zp,&
					    X_x,X_y,X_p

    DOUBLE PRECISION :: Mat_in(1:threeNM1SQ,1:threeNM1SQ),&
			X_in(1:threeNM1SQ),&
			RHS_in(1:threeNM1SQ),&
			diag_precon(1:threeNM1SQ),&
			Matxx(1:NM1SQ,1:NM1SQ),Matyy(1:NM1SQ,1:NM1SQ),&
			Matpx(1:NM1SQ,1:NM1SQ),Matpy(1:NM1SQ,1:NM1SQ)
			
    DOUBLE PRECISION, EXTERNAL :: ddot
			
! Split input matrix and vector into seperate parts:
    Matxx = Mat_in(1:NM1SQ,1:NM1SQ)
    Matyy = Mat_in(NM1SQ+1:2*NM1SQ,NM1SQ+1:2*NM1SQ)
    Matpx = Mat_in(2*NM1SQ+1:threeNM1SQ,1:NM1SQ)
    Matpy = Mat_in(2*NM1SQ+1:threeNM1SQ,NM1SQ+1:2*NM1SQ)
    
    RHSx = RHS_in(1:NM1SQ)
    RHSy = RHS_in(NM1SQ+1:2*NM1SQ)
    RHSp = RHS_in(2*NM1SQ+1:threeNM1SQ)

    X_x = X_in(1:NM1SQ)
    X_y = X_in(NM1SQ+1:2*NM1SQ)
    X_p = X_in(2*NM1SQ+1:threeNM1SQ)
    
    
  
    minresconv=0 ! parameter to feedback if it converged.
 
!  SUBROUTINE minres_internal_matrix_vector_mult(matxx,matyy,matpx,matpy,vecx,vecy,vecp,outx,outy,outp)
 
! Use w1xold, v1yold and v2old to store the values of our matrix times the initial guess for velocity(x & y) and pressure
!   restartcount=0

    v1xold=0d0
    v1yold=0d0
    v2old=0d0
    CALL minres_internal_matrix_vector_mult(Matxx, Matyy, Matpx, Matpy, &
					    X_x, X_y, X_p, &
					    v1xold, v1yold, v2old)

    
    v1xnew = RHSx - v1xold
    v1ynew = RHSy - v1yold
    v2new = RHSp - v2old


    v1xold=0d0
    v1yold=0d0
    v2old=0d0
  
  
    w1xold=0d0
    w1yold=0d0
    w2old=0d0
    w1xnew=0d0
    w1ynew=0d0
    w2new=0d0


    z1xold = 0d0
    z1yold = 0d0
    z2old = 0d0
! PRECON STEP !
    IF (preconflag.ne.0) THEN
! Divide by 1/diag_precon
      DO i=1,NM1SQ
	z1xnew(i)=v1xnew(i)/diag_precon(i)
	z1ynew(i)=v1ynew(i)/diag_precon(NM1SQ+i)
	z2new(i)=v2new(i)/diag_precon(2*NM1SQ+i)
      ENDDO
! END PRECON STEP !
    ELSE
      z1xnew = v1xnew
      z1ynew = v1ynew
      z2new = v2new
    ENDIF
  
! scalar data setup
    gammanew = SQRT(ddot(NM1SQ,z1xnew,1,v1xnew,1) + ddot(NM1SQ,z1ynew,1,v1ynew,1) + ddot(NM1SQ,z2new,1,v2new,1))
    IF (abs(gammanew).lt.1d-16) THEN ! Initial guess was correct. Stop.
      RETURN
    ENDIF
    gammaold = 1d0
    eta=gammanew
    sold=0d0
    snew=0d0
    cold=1d0
    cnew=1d0
    counter=0
    stag=0
  
     
    DO
      counter = counter + 1
      z1xnew=z1xnew/gammanew
      z1ynew=z1ynew/gammanew
      z2new=z2new/gammanew
      tempx1 = v1xnew
      tempy1 = v1ynew
      temp2 = v2new


    
      CALL minres_internal_matrix_vector_mult(Matxx, Matyy, Matpx, Matpy, &
					    z1xnew, z1ynew, z2new, &
					    xstore1, ystore1, x2)
					    
! Essentially we have 	[Ax + x1] = [ A 0 B'][z1xnew] (=xstore1)
!		        [Ay + y1]   [ 0 A B'][z1ynew] (=ystore1)
!			[Bx + By]   [ B B C][z2new ] (=x2)

      delta = ddot(NM1SQ,z1xnew,1,xstore1,1) + ddot(NM1SQ,z1ynew,1,ystore1,1) + ddot(NM1SQ,z2new,1,x2,1)
      v1xnew = xstore1 - delta*v1xnew/gammanew - gammanew*v1xold/gammaold
      v1ynew = ystore1 - delta*v1ynew/gammanew - gammanew*v1yold/gammaold
      v2new = x2 - delta*v2new/gammanew - gammanew*v2old/gammaold

      v1xold=tempx1
      v1yold=tempy1
      v2old=temp2

      z1xold = z1xnew
      z1yold = z1ynew
      z2old = z2new
      

    
      IF (preconflag.ne.0) THEN
! PRECON STEP !
	DO i=1,NM1SQ
	  z1xnew(i)=v1xnew(i)/diag_precon(i)
	  z1ynew(i)=v1ynew(i)/diag_precon(NM1SQ+i)
	  z2new(i)=v2new(i)/diag_precon(2*NM1SQ+i)
	ENDDO
    
      ELSE
	z1xnew = v1xnew
	z1ynew = v1ynew
	z2new = v2new
      ENDIF
    
      gammaold=gammanew
      gammanew=SQRT(ddot(NM1SQ,z1xnew,1,v1xnew,1) + ddot(NM1SQ,z1ynew,1,v1ynew,1) + ddot(NM1SQ,z2new,1,v2new,1))
      alpha0 = cnew*delta - cold*snew*gammaold
      alpha1 = SQRT(alpha0**2 + gammanew**2)
      alpha2 = snew*delta + cold*cnew*gammaold
      alpha3 = sold*gammaold
      cold = cnew
      cnew = alpha0/alpha1
      sold = snew
      snew = gammanew/alpha1
      tempx1 = w1xnew
      tempy1 = w1ynew
      temp2 = w2new

      w1xnew = (z1xold - alpha3*w1xold - alpha2*w1xnew)/alpha1
      w1ynew = (z1yold - alpha3*w1yold - alpha2*w1ynew)/alpha1
      w2new = (z2old - alpha3*w2old - alpha2*w2new)/alpha1
    
      w1xold = tempx1
      w1yold = tempy1
      w2old = temp2
    
    
      X_in(1:NM1SQ) = X_in(1:NM1SQ) + cnew*eta*w1xnew
      X_in(NM1SQ+1:2*NM1SQ) = X_in(NM1SQ+1:2*NM1SQ) + cnew*eta*w1ynew
      X_in(2*NM1SQ+1:threeNM1SQ) = X_in(2*NM1SQ+1:threeNM1SQ) + cnew*eta*w2new
    
      eta=-snew*eta
    
! calc norm of residual for current step
      normr = abs(eta)

      IF (abs(normr).lt.1d-10) THEN
	minresconv=counter
	EXIT
      ENDIF
! Check if method is stagnating
      IF (abs(alpha0).lt.1d-16) THEN
	stag = stag + 1
	IF (stag.gt.3) THEN ! 3 iterations are same - stop!
	  print *, 'MINRES method is stagnant after',counter,' iterations'
	  STOP
	ENDIF
      ELSE
	stag = 0
      ENDIF
!     IF (counter.eq.10) THEN
! !       print*,'restarting minres normr is ', abs(normr)
!       restartcount = restartcount+1
!       IF (restartcount.gt.100) THEN
! 	print *,'MINRES failed after ', counter,' iterations. Normr was ',abs(normr)
! 	STOP
!       ELSE
!  	GO TO 1001
!       ENDIF
!     ENDIF

      IF (counter.gt.500) THEN
!         print *,'MINRES failed after ', counter,' iterations. Normr was ',abs(normr)
	minresconv=counter
	EXIT
      ENDIF
    
      DO i=1,NM1SQ
	IF (IsNaN(X_in(i))) THEN
	  print*, 'ERROR! X_x(',i,') is NaN on iteration', counter
	  STOP
	ELSEIF (IsNaN(X_in(NM1SQ+i))) THEN
	  print*, 'ERROR! X_y(',i,') is NaN on iteration', counter
	  STOP
	ELSEIF (IsNaN(X_in(2*NM1SQ+i))) THEN
	  print*, 'ERROR! X_p(',i,') is NaN on iteration', counter
	  STOP
	ENDIF
      ENDDO

    ENDDO
  
  END SUBROUTINE MINRES_internal_only

  SUBROUTINE minres_internal_matrix_vector_mult(matxx,matyy,matpx,matpy,vecx,vecy,vecp,outx,outy,outp)
! All matrices in this routine are square, so no row/col issues to worry about
    IMPLICIT NONE
    CHARACTER(1) :: trans
    INTEGER :: rowscols
    DOUBLE PRECISION :: matxx(1:NM1SQ,1:NM1SQ),&
			matyy(1:NM1SQ,1:NM1SQ),&
			matpx(1:NM1SQ,1:NM1SQ),&
			matpy(1:NM1SQ,1:NM1SQ),&
			vecx(1:NM1SQ),&
			vecy(1:NM1SQ),&
			vecp(1:NM1SQ),&
			outx(1:NM1SQ),&
			outy(1:NM1SQ),&
			outp(1:NM1SQ)
			

    outx=0d0
    outy=0d0
    outp=0d0
! outx = Ax*v_x
    CALL DGEMV('N', NM1SQ, NM1SQ, 1d0, matxx, NM1SQ, &
			vecx, 1, 0d0, outx, 1)
! outx = Ax*v_x + Bx'*p
    CALL DGEMV('T', NM1SQ, NM1SQ, 1d0, matpx, NM1SQ, &
			vecp, 1, 1d0, outx, 1)
			
! outy = Ay*v_y
    CALL DGEMV('N', NM1SQ, NM1SQ, 1d0, matyy, NM1SQ, &
			vecy, 1, 0d0, outy, 1)
! outy = Ay*v_y + By'*p
    CALL DGEMV('T', NM1SQ, NM1SQ, 1d0, matpy, NM1SQ, &
			vecp, 1, 1d0, outy, 1)
			
! outp = Bx*v_x
    CALL DGEMV('N', NM1SQ, NM1SQ, 1d0, matpx, NM1SQ, &
			vecx, 1, 0d0, outp, 1)
! outp = Bx*v_x + By*v_y
    CALL DGEMV('N', NM1SQ, NM1SQ, 1d0, matpy, NM1SQ, &
			vecy, 1, 1d0, outp, 1)

! Zero integral of pressure bit can be put in here, will require 2 level 1 BLAS routines for v-v multiplication.
  
  END SUBROUTINE minres_internal_matrix_vector_mult
  
  SUBROUTINE calc_internal_diag_preconditioner(Mat_in,diag)
    IMPLICIT NONE
    INTEGER :: el,ij,i,j,intij
    DOUBLE PRECISION, DIMENSION(1:threeNM1SQ,numelm) :: diag
    DOUBLE PRECISION, DIMENSION(1:threeNM1SQ,1:threeNM1SQ,numelm) :: Mat_in
! Calulates the matrix with entries A(i,i) (avoiding the zeros from deleted dirichlet nodes)
    diag=0d0
    DO el=1,numelm
      DO ij=1,2*NM1SQ
	diag(ij,el) = Mat_in(ij,ij,el)
      ENDDO
    ENDDO
    DO el=1,numelm
      DO ij=2*NM1SQ+1,threeNM1SQ
	diag(ij,el) = M_pressure(ij-2*NM1SQ,ij-2*NM1SQ,el)
      ENDDO
    ENDDO
    DO el=1,numelm
      DO ij=1,threeNM1SQ
	IF (abs(diag(ij,el)).lt.1d-16) THEN
	  diag(ij,el)=1d0
	  print*, 'Warning: Diagonal preconditioner had zero value! Set to 1 instead.'
	ENDIF
      ENDDO
    ENDDO
    
  END SUBROUTINE calc_internal_diag_preconditioner
	
	
END MODULE lapack_solver_module

