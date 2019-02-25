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

MODULE pardiso_solver
  USE shared_data
  USE constants
  IMPLICIT NONE
!   INCLUDE 'mkl_pardiso.f90'

!   INTEGER*4 :: pt_stokes(64) !, ALLOCATABLE :: ! CHANGE TO INTEGER*8 for 64bit (ie merlin) 
  INTEGER*8 :: pt_stokes(64) !, ALLOCATABLE :: ! CHANGE TO INTEGER*8 for 64bit (ie merlin)
  INTEGER :: iparm_stokes(64)
  
!   INTEGER :: global_non_zeros
!          perm_stokes,&

  INTEGER, ALLOCATABLE :: Ax_local_to_sparse(:,:,:),&
    Ay_local_to_sparse(:,:,:),&
    Bx_local_to_sparse(:,:,:),&
    By_local_to_sparse(:,:,:),&
    Cp_local_to_sparse(:,:,:,:),&
    sparse_rowIndex_stokes(:),&
    sparse_columns_stokes(:)

  DOUBLE PRECISION, ALLOCATABLE :: sparse_values_stokes(:),&
    globalRHS_stokes(:),&
    globalSOL_stokes(:)

! Using MKL_ALLOC FOR MEMORY:
!   DOUBLE PRECISION :: sparse_values_stokes,&
!           globalRHS_stokes,&
!           globalSOL_stokes
          
!   POINTER    (sparse_values_stokes_PTR,sparse_values_stokes(1)),&
!        (globalRHS_stokes_PTR,globalRHS_stokes(1)),&
!        (globalSOL_stokes_PTR,globalSOL_stokes(1))
!                (perm_stokes_PTR,perm_stokes(1)),&

!   INTEGER MKL_MALLOC
!   INTEGER*8 MKL_MALLOC ! Use INTEGER*8 for 64bit (ie merlin)
!   EXTERNAL MKL_MALLOC,MKL_FREE,MKL_MEM_STAT
  
  CONTAINS
  
  SUBROUTINE assign_pardiso_memory
    IMPLICIT NONE
    
    IF (.not.ALLOCATED (Ax_local_to_sparse)) ALLOCATE(Ax_local_to_sparse(0:NP1SQM1,0:NP1SQM1,numelm))
    IF (.not.ALLOCATED (Ay_local_to_sparse)) ALLOCATE(Ay_local_to_sparse(0:NP1SQM1,0:NP1SQM1,numelm))
    IF (.not.ALLOCATED (Bx_local_to_sparse)) ALLOCATE(Bx_local_to_sparse(1:NM1SQ,0:NP1SQM1,numelm))
    IF (.not.ALLOCATED (By_local_to_sparse)) ALLOCATE(By_local_to_sparse(1:NM1SQ,0:NP1SQM1,numelm))
    IF (.not.ALLOCATED (Cp_local_to_sparse)) ALLOCATE(Cp_local_to_sparse(1:NM1SQ,1:NM1SQ,numelm,numelm))
    IF (.not.ALLOCATED (sparse_rowIndex_stokes)) ALLOCATE(sparse_rowIndex_stokes(1:global_dim+1))
    IF (.not.ALLOCATED (sparse_columns_stokes)) ALLOCATE(sparse_columns_stokes(1:global_non_zeros))
    
    IF (.not.ALLOCATED (sparse_values_stokes)) ALLOCATE(sparse_values_stokes(1:global_non_zeros))
    IF (.not.ALLOCATED (globalRHS_stokes)) ALLOCATE(globalRHS_stokes(1:global_dim))
    IF (.not.ALLOCATED (globalSOL_stokes)) ALLOCATE(globalSOL_stokes(1:global_dim))
!     IF (.not.ALLOCATED (perm_stokes)) ALLOCATE(perm_stokes(1:global_dim))

! USING MKL_ALLOC:
!     perm_stokes_PTR = MKL_MALLOC(4*global_dim,64)
!     sparse_values_stokes_PTR = MKL_MALLOC(8*global_non_zeros,64)
!     globalSOL_stokes_PTR = MKL_MALLOC(8*global_dim,64)
!     globalRHS_stokes_PTR = MKL_MALLOC(8*global_dim,64)
!     print*,sparse_values_stokes_PTR,globalRHS_stokes_PTR,globalSOL_stokes_PTR
    
    
  END SUBROUTINE assign_pardiso_memory

  
  SUBROUTINE release_pardiso_memory
    IMPLICIT NONE
    
    CALL Pardiso_Release_Memory(pt_stokes, iparm_stokes, global_dim)
    
    IF (ALLOCATED (Ax_local_to_sparse)) DEALLOCATE(Ax_local_to_sparse)
    IF (ALLOCATED (Ay_local_to_sparse)) DEALLOCATE(Ay_local_to_sparse)
    IF (ALLOCATED (Bx_local_to_sparse)) DEALLOCATE(Bx_local_to_sparse)
    IF (ALLOCATED (By_local_to_sparse)) DEALLOCATE(By_local_to_sparse)
    IF (ALLOCATED (Cp_local_to_sparse)) DEALLOCATE(Cp_local_to_sparse)
    IF (ALLOCATED (sparse_rowIndex_stokes)) DEALLOCATE(sparse_rowIndex_stokes)
    IF (ALLOCATED (sparse_columns_stokes)) DEALLOCATE(sparse_columns_stokes)
    IF (ALLOCATED (sparse_values_stokes)) DEALLOCATE(sparse_values_stokes)
    IF (ALLOCATED (globalRHS_stokes)) DEALLOCATE(globalRHS_stokes)
    IF (ALLOCATED (globalSOL_stokes)) DEALLOCATE(globalSOL_stokes)
!     IF (ALLOCATED (perm_stokes)) DEALLOCATE(perm_stokes)

! USING MKL_ALLOC TYPE WAY:
!     CALL MKL_FREE_BUFFERS
!     CALL MKL_FREE(perm_stokes_PTR)
!     CALL MKL_FREE(sparse_values_stokes_PTR)
!     CALL MKL_FREE(globalRHS_stokes_PTR)
!     CALL MKL_FREE(globalSOL_stokes_PTR)
    
    
  
  END SUBROUTINE release_pardiso_memory
  
  SUBROUTINE setup_iparm(iparm_in)
  
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: iparm_in(64)

    iparm_in = 0 ! set default
! INPUT ARGS:
    iparm_in(1) = 1 ! 1 = non-default values 
    iparm_in(2) = 0 ! 0=minimum degree algo | 2 = nested dissection(METIS) DEFAULT | 3 = openMP nested dissection
    iparm_in(4) = 0 ! Preconned CGS/CG. 0 = don't perform | 10*L + 1 w/ stop crit=10^(-L)  | 10*L+2 = same but for symmetric matrix.
    iparm_in(5) = 0 ! permutation array | 0 = use default.
    iparm_in(6) = 0 ! write solution to RHS 0 = no | 1 = yes.
    iparm_in(8) = 0 ! maximum number of iterative refinement steps. Default is 0 => 2 performed.
    iparm_in(10) = 8 ! perturbation for zero pivots, uses 10^(value) ! Default for symmetric is 8
    iparm_in(11) = 1 ! scaling, default off for symmetric indef. on for nonsymmetric.
    iparm_in(12) = 0 ! solver with transposed/conjugate transposed matrix A. ! default 0 (off) | 1 = conjugate tranpose | 2 = transpose.
    iparm_in(13) = 1 ! improved accuracy using non-symmetric weighted matching. !0 off (default for symmetric indef matrices). 1 on. Recommended to use both scaling and this on for highly indefinite symm matrices.
    iparm_in(18) = -1 ! report number of non-zero elements in factors (default -1)
    iparm_in(19) = 0 ! report mflops to factor A, 0 = off (default 0)
    iparm_in(21) = 1 ! pivoting for symmetric indef matrices. ! default 1
    iparm_in(24) = 0 ! Parallel fact control. default 0.
    iparm_in(25) = 0 !Parallel forward/backward solve. default 0 (on). 1 = off (sequential)
    iparm_in(27) = 0 ! check sparse representation.
    iparm_in(28) = 0 ! single or double precision. 0 = double.
    iparm_in(31) = 0 ! Partial solve and computing selected components of solution vectors. Makes fact slower, but solver fast. 0 = off | see docs for more info. MIGHT BE USEFUL
    iparm_in(34) = 0 ! Optimal number of threads for conitional numberical reproducibility. default 0 = automatic.
    iparm_in(35) = 0 ! fortran or C style indexing. 0 = fortran (start at 1)
    iparm_in(60) = 0 ! 0 = In-core | 2 = Out-of-core (reduces memory required, but slower)
    
! OUTPUT ARGS:
!     iparm_in(7) = ! iterative refinement steps reqd
!     iparm_in(14) = ! number of perturbed pivots
!     iparm_in(15) = ! peak memory for symbolic fact.
!     iparm_in(16) = ! permanent memory for symbolic fact.
!     iparm_in(17) = ! peak memory for numerical fact and soln
!     iparm_in(20) = ! CG/CGS diagnostics
!     iparm_in(22) = ! number of positive eigenvalues
!     iparm_in(23) = ! number of negative eigenvalues
!     iparm_in(30) = ! number of zero/negative pivots.

! ALL OTHERS ARE ZERO.

  END SUBROUTINE setup_iparm
  
  SUBROUTINE initialise_pardiso_stokes
    IMPLICIT NONE
    INTEGER :: i, nrhs, mnum, maxfct, msglvl, mtype, phase, idum, error
    DOUBLE PRECISION :: ddum
    
    CALL mkl_set_num_threads(1)
    CALL calc_global_symmetric_non_zeroes(global_non_zeros)
    
    CALL assign_pardiso_memory      
        
    CALL initialise_local_to_global_sparse(sparse_values_stokes, &
             sparse_rowIndex_stokes, &
             sparse_columns_stokes)
             
             




    nrhs = 1 ! number of rhs vectors
    mnum = 1 ! number of matrices
    maxfct = 1 ! max number of matrices (i think)
    msglvl = 0 !print no (=0) statistical information on screen
    ddum = 0d0 ! dummy variables for things not used at this stage.
    idum = 0
    
!    mtype = 11 !real unsymmetric matrix 
    mtype = -2 !for real symmetric INDEFINITE matrix (2 = positive definite)

! initialise pardiso.
    iparm_stokes=0
    pt_stokes=0
    CALL pardisoinit(pt_stokes, mtype, iparm_stokes)
    
    CALL setup_iparm(iparm_stokes) ! use non-default
!     iparm_stokes(1)=0 ! use default
    
! Reordering step - only required once.    
    phase = 11 ! 11 perform analysis and numerical factorisation, 22 just numerical factorisation

    CALL pardiso (pt_stokes, maxfct, mnum, mtype, phase, global_dim, sparse_values_stokes,&
      sparse_rowIndex_stokes, sparse_columns_stokes,& 
      idum, nrhs, iparm_stokes, msglvl, ddum, ddum, error)

!     WRITE(*,*) 'Reordering completed ... '
    IF (error.ne.0) THEN
      WRITE(*,*) 'At phase 11 of PARDISO: The following ERROR was detected: ', error
      STOP
    ENDIF
!     WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
!     WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)

! If the mesh is non-moving then this is only required once, so we do it now.
! If the mesh is moving we will need to update the values vector and update the factorisation.
    CALL Pardiso_Fact(pt_stokes, iparm_stokes, global_non_zeros, global_dim, &
      sparse_values_stokes, sparse_columns_stokes, sparse_rowIndex_stokes)

  END SUBROUTINE initialise_pardiso_stokes
  
  SUBROUTINE pardiso_solve_stokes
    IMPLICIT NONE
    INTEGER :: i
    
    IF (movingmeshflag.eq.1) THEN !.or.param_delta_a.gt.0d0

      CALL update_global_sparse_storage(sparse_values_stokes,global_non_zeros)

      CALL Pardiso_Fact(pt_stokes, iparm_stokes, global_non_zeros, global_dim, &
      sparse_values_stokes, sparse_columns_stokes, sparse_rowIndex_stokes)

    ENDIF
    
!     call cpu_time(cputime1)
    CALL update_RHS_stokes(globalRHS_stokes)
    
    CALL Pardiso_Solve(pt_stokes, iparm_stokes, global_non_zeros, global_dim,&
      sparse_values_stokes, sparse_columns_stokes, sparse_rowIndex_stokes,&
      globalRHS_stokes, globalSOL_stokes)

    CALL PARDISO_UPDATE_SOLUTION(globalSOL_stokes,V_x,V_y,pressure)
!     call cpu_time(cputime2)
!     print*,'solve',cputime2-cputime1
  END SUBROUTINE pardiso_solve_stokes
  
  SUBROUTINE Pardiso_Fact(pt, iparm, nnz, mat_dim, values, columns, rowindex)
    IMPLICIT NONE
    
!     INTEGER*4 :: pt(64)
    INTEGER*8 :: pt(64)
!     TYPE(MKL_PARDISO_HANDLE) :: pt(64)

    INTEGER :: maxfct, mnum, mtype, phase, nrhs, error, msglvl, i,&
      nnz, mat_dim, idum,&
!       perm(mat_dim),&
      rowindex(mat_dim+1),&
      columns(nnz),&
      iparm(64)

    DOUBLE PRECISION :: ddum,&
      values(nnz)
    
    nrhs = 1 ! number of rhs vectors
    mnum = 1 ! number of matrices
    maxfct = 1 ! max number of matrices (i think)
    msglvl = 0 !print no (=0) statistical information on screen
    ddum = 0d0 ! double dummy variable for things not used at this stage.
    idum = 0   ! integer dummy variable for things not used at this stage.

!    mtype = 11 !real unsymmetric matrix 
    mtype = -2 !for real symmetric INDEFINITE matrix (2 = positive definite)
    
!     phase = 11 ! 11 perform analysis and numerical factorisation, 22 just numerical factorisation
!     
!     CALL pardiso (pt, maxfct, mnum, mtype, phase, mat_dim, values, rowindex, columns,& 
!                  idum, nrhs, iparm, msglvl, ddum, ddum, error) !perm
! 
! !     WRITE(*,*) 'Reordering completed ... '
!     IF (error.ne.0) THEN
!       WRITE(*,*) 'At phase 11 of PARDISO: The following ERROR was detected: ', error
!       STOP
!     END IF
!     WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
!     WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)
    
    !.. Factorization.
    phase = 22 ! only factorization
    CALL pardiso (pt, maxfct, mnum, mtype, phase, mat_dim, values, rowindex, columns,& 
      idum, nrhs, iparm, msglvl, ddum, ddum, error) !perm

    IF (error.ne.0) THEN
      WRITE(*,*) 'At phase 22 of PARDISO: The following ERROR was detected: ', error
      STOP
    END IF
  END SUBROUTINE Pardiso_Fact

  SUBROUTINE Pardiso_Solve(pt, iparm, nnz, mat_dim, values, columns, rowindex, rhs_in, sol_out)
    IMPLICIT NONE

!     INTEGER*4, INTENT(INOUT) :: pt(64)
    INTEGER*8, INTENT(INOUT) :: pt(64)
!     TYPE(MKL_PARDISO_HANDLE) :: pt(64)
    INTEGER, INTENT (INOUT) :: iparm(64)
    INTEGER, INTENT (IN) :: nnz,mat_dim, &
!       perm(mat_dim), &
      rowindex(mat_dim+1), &
      columns(nnz)

    INTEGER :: i,maxfct, mnum, mtype, phase, nrhs, error, msglvl, idum
    DOUBLE PRECISION, INTENT(IN) :: values(nnz),&
      rhs_in(mat_dim)

    DOUBLE PRECISION, INTENT(OUT) :: sol_out(mat_dim)
    DOUBLE PRECISION :: ddum
      
      

    
    nrhs = 1 ! number of rhs vectors
    mnum = 1 ! number of matrices
    maxfct = 1 ! max number of matrices (i think)

!     mtype = 11 !real unsymmetric matrix 
    mtype = -2 !for real symmetric INDEFINITE matrix (2 = positive definite)

    phase = 33 !solves equation

    msglvl = 0 !print no (=0) statistical information on screen 

    ddum = 0d0
    idum = 0

   
    do i = 1, mat_dim
      if((rhs_in(i)+1d0).eq.rhs_in(i).or.isNaN(rhs_in(i)))then
        write(*,*) 'globalRHS_stokes is NaN, before solve sym',i,rhs_in(i)
        stop
      endif
    enddo


    call pardiso (pt, maxfct, mnum, mtype, phase, mat_dim, &
      values, rowindex, columns,& 
      idum, nrhs, iparm, msglvl, &
      rhs_in, sol_out, error)


    do i = 1, mat_dim
      if((sol_out(i)+1d0).eq.sol_out(i)) then
        write(*,*) 'sol_out is NaN, after solve',i
        stop
      endif
    enddo
    
    IF (error.ne.0) THEN
      WRITE(0,*) 'At phase 33 of PARDISO: The following ERROR was detected: ', error
!       STOP
    END IF
  END SUBROUTINE Pardiso_Solve
  
  SUBROUTINE Pardiso_Release_Memory(pt, iparm, mat_dim)
    IMPLICIT NONE

!     INTEGER*4 :: pt(64)
    INTEGER*8 :: pt(64)
!     TYPE(MKL_PARDISO_HANDLE) :: pt(64)
    INTEGER :: maxfct, mnum, mtype, phase, nrhs, error, msglvl, idum, mat_dim,&
      iparm(64)

    DOUBLE PRECISION :: ddum

    nrhs = 1 ! number of rhs vectors
    mnum = 1 ! number of matrices
    maxfct = 1 ! max number of matrices (i think)

!     mtype = 11 !real unsymmetric matrix 
    mtype = -2 !for real symmetric INDEFINITE matrix (2 = positive definite)

    phase = -1 !release memory

    msglvl = 0 !print no (=0) statistical information on screen 


    idum = 0 !dummy variables which aren't used for this particular phase
    ddum = 0

    call pardiso(pt, maxfct, mnum, mtype, phase, mat_dim, ddum, idum, idum,& 
      idum, nrhs, iparm, msglvl, ddum, ddum, error) 

  END SUBROUTINE Pardiso_Release_Memory


  SUBROUTINE PARDISO_UPDATE_SOLUTION(sol_in,out_x,out_y,out_p)
    IMPLICIT NONE
    INTEGER :: rowcount,i
    DOUBLE PRECISION, INTENT(IN) :: sol_in(1:global_bd_dim+3*npint)
    DOUBLE PRECISION, INTENT(OUT) :: out_x(1:nptot),&
      out_y(1:nptot),&
      out_p(1:npint)
! Transfer global solution to V_x, V_y, p
! release global memory.
    
    rowcount=0
    DO i=1,nptot
      IF (bdflag(1,i)) THEN
        out_x(i) = boundary_x(i)
      ELSE
        rowcount=rowcount+1
        out_x(i)= sol_in(rowcount)
      ENDIF
    ENDDO
    DO i=1,nptot
      IF (bdflag(2,i)) THEN
        out_y(i) = boundary_y(i)
      ELSE
        rowcount=rowcount+1
        out_y(i)= sol_in(rowcount)
      ENDIF
    ENDDO
    DO i=1,npint
      rowcount=rowcount+1
      out_p(i)= sol_in(rowcount)
    ENDDO

  END SUBROUTINE PARDISO_UPDATE_SOLUTION


  SUBROUTINE update_global_sparse_storage(values,nnz)
! Only needed if the mesh is moving, otherwise our matrix does not change.
    IMPLICIT NONE
    INTEGER :: el,ij,intij,kl,intkl,el1,el2,nnz
    DOUBLE PRECISION :: values(1:nnz)
    DOUBLE PRECISION, ALLOCATABLE :: temp(:)
    
    ALLOCATE(temp(0:nnz))
    temp=0d0

!     temp(1:global_non_zeros)=sparse_values_stokes
!       print*,global_non_zeros
!       DO i=0,global_non_zeros
!   temp=0d0
!       ENDDO
!DEC$ NOPARALLEL
    DO el=1,numelm
      DO ij=0,NP1SQM1
        DO kl=0,NP1SQM1
!       IF (Ax_local_to_sparse(kl,ij,el).lt.0) print*,'ALERT: Ax_local_to_sparse DODGY VALUE',Ax_local_to_sparse(kl,ij,el)
!       IF (Ay_local_to_sparse(kl,ij,el).lt.0) print*,'ALERT: Ay_local_to_sparse DODGY VALUE',Ay_local_to_sparse(kl,ij,el)
          temp(Ax_local_to_sparse(kl,ij,el)) = temp(Ax_local_to_sparse(kl,ij,el)) + A_x(kl,ij,el)
          temp(Ay_local_to_sparse(kl,ij,el)) = temp(Ay_local_to_sparse(kl,ij,el)) + A_y(kl,ij,el)
        ENDDO
        DO intkl=1,NM1SQ
!           IF (Bx_local_to_sparse(intkl,ij,el).lt.0) print*,'ALERT: Ax_local_to_sparse DODGY VALUE',Bx_local_to_sparse(intkl,ij,el)
!           IF (By_local_to_sparse(intkl,ij,el).lt.0) print*,'ALERT: Ay_local_to_sparse DODGY VALUE',By_local_to_sparse(intkl,ij,el)
          temp(Bx_local_to_sparse(intkl,ij,el)) = temp(Bx_local_to_sparse(intkl,ij,el)) + B_x(intkl,ij,el)
          temp(By_local_to_sparse(intkl,ij,el)) = temp(By_local_to_sparse(intkl,ij,el)) + B_y(intkl,ij,el)
        ENDDO
      ENDDO
    ENDDO

    IF (param_alphaZ.gt.0d0) THEN
!DEC$ NOPARALLEL    
      DO el2=1,numelm
        DO el1=1,numelm
          DO intkl=1,NM1SQ
            DO intij=1,NM1SQ
              temp(Cp_local_to_sparse(intij,intkl,el1,el2)) = param_alphaZ*Z_p(intij,el1)*Z_p(intkl,el2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    values = temp(1:global_non_zeros)
    DEALLOCATE(temp)
  END SUBROUTINE update_global_sparse_storage
  
  SUBROUTINE update_RHS_stokes(rhs_out)
    IMPLICIT NONE
    INTEGER :: i,fixed1,fixed2
    DOUBLE PRECISION :: rhs_out(global_dim)
    
!     IF ((global_bd_dim+1+npint).gt.global_dim) THEN
!       print*,'Memory leak in update_RHS_stokes'
!       stop
!     ENDIF
! This will write some values to global_bd_dim+1 in the rhs_out vector.
! However, the next loop will copy over them with the 1st interior value of f_x/f_y
    DO i=1,npedg
      rhs_out(non_dir_bd_map_x(i)) = f_x(i)
    ENDDO
    DO i=1,npedg
      rhs_out(npint+non_dir_bd_map_y(i)) = f_y(i)
    ENDDO
    fixed1=global_bd_dim+npint
    fixed2=fixed1+npint
    DO i=1,npint
      rhs_out(bd_numrows_x+i) = f_x(npedg+i)
      rhs_out(fixed1+i) = f_y(npedg+i)
      rhs_out(fixed2+i) = g(i)
    ENDDO
  END SUBROUTINE update_RHS_stokes
  
  SUBROUTINE calc_global_symmetric_non_zeroes(nnz)
    IMPLICIT NONE
    LOGICAL,ALLOCATABLE :: seen_nodex(:,:),seen_nodey(:,:)

    INTEGER :: el,el1,el2,i,j,ij,kl,nnz

    nnz=0
    ALLOCATE(seen_nodex(npedg,npedg),seen_nodey(npedg,npedg))
    seen_nodex=.false.
    seen_nodey=.false.

    DO el=1,numelm
      DO ij=0,NP1SQM1
! scan through A_x
        i=mapg(ij,el)
        IF (.not.bdflag(1,i)) THEN
          DO kl=0,NP1SQM1
            j=mapg(kl,el)
            IF (.not.bdflag(1,j)) THEN
              IF (i.le.npedg.and.j.le.npedg) THEN        
                IF (seen_nodex(i,j)) CYCLE    
                IF (i.eq.j) THEN
                  nnz=nnz+1
                  seen_nodex(i,j)=.true.        
                ELSEIF (abs(A_x(ij,kl,el)).gt.1d-15) THEN
                  nnz=nnz+1
                  seen_nodex(i,j)=.true.
                ENDIF
              ELSE
                IF (i.eq.j) THEN
                  nnz=nnz+1
                ELSEIF (abs(A_x(ij,kl,el)).gt.1d-15) THEN
                  nnz=nnz+1      
                ENDIF
              ENDIF
            ENDIF
          ENDDO
! scan through B_x (add 2 because of B_x transpose too)
          DO kl=1,NM1SQ
            IF (abs(B_x(kl,ij,el)).gt.1d-15) THEN
              nnz=nnz+2  
            ENDIF
          ENDDO
        ENDIF
! scan through A_y
        IF (.not.bdflag(2,i)) THEN
          DO kl=0,NP1SQM1
            j=mapg(kl,el)
            IF (.not.bdflag(2,j)) THEN
              IF (i.le.npedg.and.j.le.npedg) THEN        
                IF (seen_nodey(i,j)) CYCLE    
                IF (i.eq.j) THEN
                  nnz=nnz+1
                  seen_nodey(i,j)=.true.        
                ELSEIF (abs(A_y(ij,kl,el)).gt.1d-15) THEN
                  nnz=nnz+1
                  seen_nodey(i,j)=.true.
                ENDIF
              ELSE
                IF (i.eq.j) THEN
                  nnz=nnz+1      
                ELSEIF (abs(A_y(ij,kl,el)).gt.1d-15) THEN
                  nnz=nnz+1      
                ENDIF
              ENDIF
            ENDIF
          ENDDO
! scan through B_y (add 2 because of B_y transpose too)
          DO kl=1,NM1SQ
            IF (abs(B_y(kl,ij,el)).gt.1d-15) THEN
              nnz=nnz+2        
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDDO

    IF (param_alphaZ.gt.0d0) THEN
      DO el1=1,numelm
        DO ij=1,NM1SQ
          DO el2=1,numelm
            DO kl=1,NM1SQ
              IF (mapg_pressure(ij,el1).eq.mapg_pressure(kl,el2)) THEN
                nnz=nnz+1
              ELSEIF (abs(param_alphaZ*Z_p(ij,el1)*Z_p(kl,el2)).gt.1d-15) THEN
                nnz=nnz+1
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ELSE     
! zero-entries on the diagonal of lower block diagonal of zeros.
      nnz=nnz+npint
    ENDIF

! Now have non-zero entries for the whole matrix.
! Want the number of entries in the upper triangle:
    IF (MOD(nnz-global_dim,2).ne.0) THEN
      print*,'ERROR: Odd number of off-diagonal entries in calc_global_symmetric_non_zeroes!'
      print*, 'nnz: ',nnz,'global_dim', global_dim
      STOP
    ENDIF
    nnz=int((nnz-global_dim)/int(2)) + global_dim ! non-diagonal entries must be divisible by 2 as the matrix is symmetric.

    DEALLOCATE(seen_nodex,seen_nodey)
  END SUBROUTINE calc_global_symmetric_non_zeroes
  
  SUBROUTINE initialise_local_to_global_sparse(values,rowIndex,columns)

    IMPLICIT NONE
    LOGICAL :: seen_node
    INTEGER :: i,j,el,el1,el2,ij,intij,kl,intkl,non_zero_count
    INTEGER, INTENT(OUT) :: rowindex(1:global_bd_dim+3*npint+1),&
          columns(1:global_non_zeros)
    DOUBLE PRECISION :: temp
    DOUBLE PRECISION, INTENT(OUT) :: values(1:global_non_zeros)

    
    Ax_local_to_sparse=0
    Ay_local_to_sparse=0
    Bx_local_to_sparse=0
    By_local_to_sparse=0
    Cp_local_to_sparse=0
    
    non_zero_count=0
!     seen_node=.false.
! boundary rows of x-part of global matrix:
    
    DO i=1,npedg
      IF (.not.bdflag(1,i)) THEN
        rowindex(non_dir_bd_map_x(i))=non_zero_count+1
        DO j=i,npedg
          seen_node=.false.
          IF (.not.bdflag(1,j)) THEN
            DO el=1,numelm
              ij = global_to_local_map(i,el)
              IF (ij.lt.0) CYCLE
              kl = global_to_local_map(j,el)
              IF (kl.lt.0) CYCLE
              IF (i.eq.j) THEN ! Always store the diagonal entry, even if it is zero.
                IF (seen_node) THEN
                  Ax_local_to_sparse(ij,kl,el)=non_zero_count
                  values(non_zero_count) = values(non_zero_count) + A_x(ij,kl,el)
                ELSE
                  non_zero_count = non_zero_count+1
                  Ax_local_to_sparse(ij,kl,el) = non_zero_count
                  values(non_zero_count) = A_x(ij,kl,el)
                  columns(non_zero_count) = non_dir_bd_map_x(j)
                  seen_node=.true.
                ENDIF
              ELSEIF (abs(A_x(ij,kl,el)).gt.1d-15) THEN
                IF (seen_node) THEN
                  Ax_local_to_sparse(ij,kl,el)=non_zero_count
                  values(non_zero_count) = values(non_zero_count) + A_x(ij,kl,el)
                ELSE
                  non_zero_count = non_zero_count+1
                  Ax_local_to_sparse(ij,kl,el) = non_zero_count
                  values(non_zero_count) = A_x(ij,kl,el)
                  columns(non_zero_count) = non_dir_bd_map_x(j)
                  seen_node=.true.
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        A_x_interior_loop1: DO j=npedg+1,nptot
          DO el=1,numelm
            ij = global_to_local_map(i,el)
            IF (ij.lt.0) CYCLE
            kl = global_to_local_map(j,el)
            IF (kl.lt.0) CYCLE
            IF (abs(A_x(ij,kl,el)).gt.1d-15) THEN
              non_zero_count=non_zero_count+1
              Ax_local_to_sparse(ij,kl,el)=non_zero_count
              values(non_zero_count) = A_x(ij,kl,el)
              columns(non_zero_count)=j-npedg+bd_numrows_x
      
              CYCLE A_x_interior_loop1 ! Exit as we will only find 1 contribution to this node - it is interior
            ENDIF
          ENDDO
        ENDDO A_x_interior_loop1
        B_x_interior_loop1: DO j=1,npint
          DO el=1,numelm
            ij = global_to_local_map(npedg+j,el)
            IF (ij.lt.0) CYCLE
            intij = local_to_interior_node(ij)
            kl = global_to_local_map(i,el)
            IF (kl.lt.0) CYCLE
            IF (abs(B_x(intij,kl,el)).gt.1d-15) THEN
              non_zero_count=non_zero_count+1
              Bx_local_to_sparse(intij,kl,el)=non_zero_count
              values(non_zero_count) = B_x(intij,kl,el)
              columns(non_zero_count)=j+global_bd_dim+2*npint
              
              CYCLE B_x_interior_loop1 ! Exit as we will only find 1 contribution to this node - it is partly interior
            ENDIF
          ENDDO
        ENDDO B_x_interior_loop1
      ENDIF  
    ENDDO ! Now done all boundary nodes of the x-part of the global matrix.
    
    

! Interior rows of the x-part of the global matrix
    DO i=npedg+1,nptot
      rowindex(i-npedg+bd_numrows_x)=non_zero_count+1

      A_x_interior_loop2: DO j=i,nptot
        DO el=1,numelm
          ij = global_to_local_map(i,el)
          IF (ij.lt.0) CYCLE
          kl = global_to_local_map(j,el)
          IF (kl.lt.0) CYCLE
          IF (i.eq.j) THEN ! Always store the diagonal entry, even if it is zero.
            non_zero_count=non_zero_count+1
            Ax_local_to_sparse(ij,kl,el)=non_zero_count
            values(non_zero_count) = A_x(ij,kl,el)
            columns(non_zero_count) = j - npedg + bd_numrows_x

            CYCLE A_x_interior_loop2 ! Exit as we will only find 1 contribution to this node - it is interior
          ELSEIF (abs(A_x(ij,kl,el)).gt.1d-15) THEN
            non_zero_count=non_zero_count+1
            Ax_local_to_sparse(ij,kl,el)=non_zero_count
            values(non_zero_count) = A_x(ij,kl,el)
            columns(non_zero_count) = j-npedg+bd_numrows_x

            CYCLE A_x_interior_loop2 ! Exit as we will only find 1 contribution to this node - it is interior
          ENDIF
        ENDDO
      ENDDO A_x_interior_loop2
      B_x_interior_loop2: DO j=1,npint
        DO el=1,numelm
          ij = global_to_local_map(npedg+j,el)
          IF (ij.lt.0) CYCLE
          intij = local_to_interior_node(ij)
          kl =  global_to_local_map(i,el)
          IF (kl.lt.0) CYCLE
          IF (abs(B_x(intij,kl,el)).gt.1d-15) THEN
            non_zero_count=non_zero_count+1
            Bx_local_to_sparse(intij,kl,el)=non_zero_count
            values(non_zero_count) = B_x(intij,kl,el)
            columns(non_zero_count)=j+global_bd_dim+2*npint
            
            CYCLE B_x_interior_loop2 ! Exit as we will only find 1 contribution to this node - it is partly interior
          ENDIF  
        ENDDO
      ENDDO B_x_interior_loop2
    ENDDO ! Now done all interior nodes of the x-part of the global matrix.

! boundary rows of y-part of global matrix:
    DO i=1,npedg
      IF (.not.bdflag(2,i)) THEN
        rowindex(npint+non_dir_bd_map_y(i))=non_zero_count+1
        DO j=i,npedg
          seen_node=.false.
          IF (.not.bdflag(2,j)) THEN
            DO el=1,numelm
              ij = global_to_local_map(i,el)
              IF (ij.lt.0) CYCLE
              kl =  global_to_local_map(j,el)
              IF (kl.lt.0) CYCLE
              IF (i.eq.j) THEN ! Always store the diagonal entry, even if it is zero.
                IF (seen_node) THEN
                  Ay_local_to_sparse(ij,kl,el )= non_zero_count
                  values(non_zero_count) = values(non_zero_count) + A_y(ij,kl,el)
                ELSE
                  non_zero_count = non_zero_count+1
                  Ay_local_to_sparse(ij,kl,el) = non_zero_count
                  values(non_zero_count) = A_y(ij,kl,el)
                  columns(non_zero_count) = npint + non_dir_bd_map_y(j)
                  seen_node=.true.
                  
                ENDIF
              ELSEIF (abs(A_y(ij,kl,el)).gt.1d-15) THEN
                IF (seen_node) THEN
                  Ay_local_to_sparse(ij,kl,el) = non_zero_count  
                  values(non_zero_count) = values(non_zero_count) + A_y(ij,kl,el)
                ELSE
                  non_zero_count = non_zero_count+1
                  Ay_local_to_sparse(ij,kl,el) = non_zero_count
                  values(non_zero_count) = A_y(ij,kl,el)
                  columns(non_zero_count) = npint + non_dir_bd_map_y(j)
                  seen_node=.true.
                  
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        A_y_interior_loop1: DO j=npedg+1,nptot
          DO el=1,numelm
            ij = global_to_local_map(i,el)
            IF (ij.lt.0) CYCLE
            kl =  global_to_local_map(j,el)
            IF (kl.lt.0) CYCLE
            IF (abs(A_y(ij,kl,el)).gt.1d-15) THEN
              non_zero_count=non_zero_count+1
              Ay_local_to_sparse(ij,kl,el)=non_zero_count
              values(non_zero_count) = A_y(ij,kl,el)
              columns(non_zero_count)=j - npedg + global_bd_dim + npint           
              
              CYCLE A_y_interior_loop1 ! Exit as we will only find 1 contribution to this node - it is interior
            ENDIF
          ENDDO
        ENDDO A_y_interior_loop1
        B_y_interior_loop1: DO j=1,npint
          DO el=1,numelm
            ij = global_to_local_map(npedg+j,el)
            IF (ij.lt.0) CYCLE
            intij = local_to_interior_node(ij)
            kl =  global_to_local_map(i,el)
            IF (kl.lt.0) CYCLE
            IF (abs(B_y(intij,kl,el)).gt.1d-15) THEN
              non_zero_count=non_zero_count+1
              By_local_to_sparse(intij,kl,el)=non_zero_count
              values(non_zero_count) = B_y(intij,kl,el)
              columns(non_zero_count)=j+global_bd_dim+2*npint        

              CYCLE B_y_interior_loop1 ! Exit as we will only find 1 contribution to this node - it is partly interior
            ENDIF  
          ENDDO
        ENDDO B_y_interior_loop1
      ENDIF  
    ENDDO ! Now done all boundary nodes of the y-part of the global matrix.

! Interior rows of the y-part of the global matrix
    DO i=npedg+1,nptot
      rowindex(i-npedg+global_bd_dim+npint)=non_zero_count+1
      A_y_interior_loop2: DO j=i,nptot
        DO el=1,numelm
          ij = global_to_local_map(i,el)
          IF (ij.lt.0) CYCLE
          kl =  global_to_local_map(j,el)
          IF (kl.lt.0) CYCLE
          IF (i.eq.j) THEN ! Always store the diagonal entry, even if it is zero.
            non_zero_count=non_zero_count+1
            Ay_local_to_sparse(ij,kl,el)=non_zero_count
            values(non_zero_count) = A_y(ij,kl,el)
            columns(non_zero_count)= j - npedg + global_bd_dim + npint      
            
            CYCLE A_y_interior_loop2 ! Exit as we will only find 1 contribution to this node - it is interior
          ELSEIF (abs(A_y(ij,kl,el)).gt.1d-15) THEN
            non_zero_count=non_zero_count+1
            Ay_local_to_sparse(ij,kl,el)=non_zero_count
            values(non_zero_count) = A_y(ij,kl,el)
            columns(non_zero_count)=j - npedg + global_bd_dim + npint      
            
            CYCLE A_y_interior_loop2 ! Exit as we will only find 1 contribution to this node - it is interior
          ENDIF
        ENDDO
      ENDDO A_y_interior_loop2

      B_y_interior_loop2: DO j=1,npint
        DO el=1,numelm
          ij = global_to_local_map(npedg+j,el)
          IF (ij.lt.0) CYCLE
          intij = local_to_interior_node(ij)
          kl =  global_to_local_map(i,el)
          IF (kl.lt.0) CYCLE
          IF (abs(B_y(intij,kl,el)).gt.1d-15) THEN
            non_zero_count=non_zero_count+1
            By_local_to_sparse(intij,kl,el)=non_zero_count
            values(non_zero_count) = B_y(intij,kl,el)
            columns(non_zero_count)=j+global_bd_dim+2*npint      
            
            CYCLE B_y_interior_loop2 ! Exit as we will only find 1 contribution to this node - it is partly interior
          ENDIF
        ENDDO
      ENDDO B_y_interior_loop2
    ENDDO ! Now done all interior nodes of the y-part of the global matrix.

! Interior rows of the pressure-part of the global matrix
! Alternative formulation:
    IF (param_alphaZ.gt.0d0) THEN
      DO i=1,npint
        rowindex(i+global_bd_dim+2*npint)=non_zero_count+1
! find correct element/local point for global i. There must be one and only 1 (internal points)
        DO el1=1,numelm
          ij = global_to_local_map(npedg+i,el1)
          IF (ij.ge.0) THEN
            intij = local_to_interior_node(ij)
            EXIT
          ENDIF
        ENDDO
        DO j=i,npint  
! find correct element/local point for global j. There must be one and only 1 (internal points)
          DO el2=1,numelm
            kl =  global_to_local_map(npedg+j,el2)
            IF (kl.ge.0) THEN
              intkl=local_to_interior_node(kl)
              EXIT
            ENDIF
          ENDDO
          temp=param_alphaZ*Z_p(intij,el1)*Z_p(intkl,el2)
          IF (i.eq.j) THEN
            non_zero_count=non_zero_count+1
            Cp_local_to_sparse(intij,intkl,el1,el2)=non_zero_count
            values(non_zero_count) = temp!param_alphaZ*Z_p(intij,el1)*Z_p(intkl,el2) !C_p(intij,intkl,el)
            columns(non_zero_count)=j+global_bd_dim+2*npint
          ELSEIF (abs(temp).gt.1d-15) THEN
            non_zero_count=non_zero_count+1
            Cp_local_to_sparse(intij,intkl,el1,el2)=non_zero_count
            values(non_zero_count) = temp!param_alphaZ*Z_p(intij,el1)*Z_p(intkl,el2)
            columns(non_zero_count)=j+global_bd_dim+2*npint      
          ENDIF
        ENDDO
      ENDDO
    ELSE
! Normal formulation:    
! These are all zero, so we need only store the diagonal entries
      DO i=1,npint
        non_zero_count=non_zero_count+1
        rowindex(i+global_bd_dim+2*npint)=non_zero_count
        values(non_zero_count) = 0d0!1d-16
        columns(non_zero_count)=i+global_bd_dim+2*npint
      ENDDO
    ENDIF
    
    rowindex(global_bd_dim+3*npint+1)=non_zero_count+1  
    
    IF (global_non_zeros.ne.non_zero_count) THEN
      write(*,*) 'Error in sparse matrix construction. Inconsistent non-zero count! global_non_zeros = ', &
      global_non_zeros,'internal_count = ',non_zero_count
      STOP
    ENDIF
!     print *, huge(global_non_zeros)
!     print*,non_zero_count,global_non_zeros
    
!     DO el=1,numelm
!     DO intij=1,NP1SQ
!     DO kl=0,NP1SQM1
!       IF Bx_local_to_sparse(intij,kl,el).gt.
  END SUBROUTINE initialise_local_to_global_sparse

END MODULE pardiso_solver
