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

 MODULE scalapack_solver_module
 ! Will require linking to blas/blacs/scalapack at compilation.
 
!   USE shared_data
!   USE functions_module
  
  IMPLICIT NONE
  
  CONTAINS
  
  SUBROUTINE scalapack_calc_internal_cholesky(A,param_N)
    IMPLICIT NONE

! SEE http://netlib.org/scalapack/slug/node30.html#secvariables FOR INFO ON THE VARIABLES !

! Parameters:
    INTEGER :: DLEN_, IA, JA, IB, JB, M, N, MB, NB, RSRC, &
		CSRC, MXLLDA, MXLLDB, NRHS, NBRHS, &
		MXLOCR, MXLOCC, MXRHSC, LOCr,LOCc, param_N
    PARAMETER :: ( DLEN_ = 9, M = 3*(param_N-1)**2, N = M )
    
    CHARACTER*1 :: uplo

! Local memory:
    INTEGER :: NP,NQ,ICTXT, INFO, MYCOL, MYROW, NPCOL, NPROW, &
		DESCA( DLEN_ ), DESCB( DLEN_ )
    DOUBLE PRECISION :: matrix_to_solve(1:(param_N-1)**4,numelm), RHS_to_solve(1:(param_N-1)**2,numelm)
    
! Set parameters:
    
    IA = 1
    JA = 1
    IB = 1
    JB = 1
    
    MB = 32
    NB = MB
    RSRC = 0
    CSRC = 0
    NP = NUMROC( M , NB, MYROW, 0, NPROW )
    NQ = NUMROC( N , NB, MYCOL, 0, NPCOL )

    MXLLDA = max( 1, NP )
    MXLLDB = MXLLDA
    NRHS = 1
    NBRHS = 1
    MXLOCR = MXLLDA
    MXLOCC = MXLLDA
    MXRHSC = 1
    
! INITIALIZE THE PROCESS GRID
    CALL SL_INIT( ICTXT,NPROW,NPCOL )
    CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
    
!     DISTRIBUTE THE MATRIX ON THE PROCESS GRID
!     Initialize the array descriptors for the matrices A and B
    CALL DESCINIT( DESCA, M, N, MB, NB, RSRC, CSRC, ICTXT, MXLLDA, dINFO )
    CALL DESCINIT( DESCB, N, NRHS, NB, NBRHS, RSRC, CSRC, ICTXT, MXLLDB, INFO )    
      
    
      
  END SUBROUTINE scalapack_solve
  
  SUBROUTINE build_internal_local_matrices(local_matrix,N,numelm,A_x,B_x,Mv_x,A_y,B_y,Mv_y)
  
    IMPLICIT NONE
    
    INTEGER :: N,numelm,ij,kl,el,NM1
    DOUBLE PRECISION :: const
    DOUBLE PRECISION, DIMENSION(1:(N-1)**4,numelm) :: local_matrix
    DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1,0:(N+1)**2,numelm) :: A_x,B_x,Mv_x,A_y,B_y,Mv_y
    
! store (N-1)^2 and (N-1)^4 for later - saves recomputing
    NM1SQ = (N-1)**2
    NM14 = NM1SQ**2
! initialise storage matrix
    local_matrix=0d0
    const = 3d0*Re/(2d0*deltat)

! Build the local "stokes" matrix (including the additional terms for OIFS, etc)
! will need to modify for other schemes, eg the crank nicolson thingys!
    DO el=1,numelm
      DO kl=1,NM1SQ
	DO ij=1,NM1SQ
	  local_matrix(ij+(kl-1)*NM1SQ) = A_x(ij,kl,el) + const*Mv_x(ij,kl,el)
	  local_matrix(ij+NM1SQ+(kl+NM1SQ-1)*NM1SQ) = A_y(ij,kl,el) + const*Mv_y(ij,kl,el)
	  local_matrix(ij+3*NM1SQ+(kl-1)*NM1SQ) = B_x(kl,ij,el)
	  local_matrix(ij+3*NM1SQ+(kl+NM1SQ-1)*NM1SQ) = B_y(kl,ij,el)
	  local_matrix(ij+NM1SQ+(kl+3*NM1SQ-1)*NM1SQ) = B_x(ij,kl,el)
	  local_matrix(ij+2*NM1SQ+(kl+3*NM1SQ-1)*NM1SQ) = B_y(ij,kl,el)
	ENDDO
      ENDDO
    ENDDO
       
  
  END SUBROUTINE build_internal_local_matrices
  
  SUBROUTINE build_internal_local_rhs(local_rhs)
  
  END SUBROUTINE build_internal_local_rhs
  
 
 END MODULE scalapack_solver_module