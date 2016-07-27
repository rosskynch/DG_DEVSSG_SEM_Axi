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

MODULE schur_complement_module
  USE shared_data
  USE constants
  IMPLICIT NONE
  
  CONTAINS
  
  SUBROUTINE schur_solve_internal
    IMPLICIT NONE
    INTEGER :: i,j,ij
    DOUBLE PRECISION ::
    DOUBLE PRECISION, DIMENSION(npint) ::
    DOUBLE PRECISION, DIMENSION(N-1,N-1,el)
    
    cholesky_fact_out
    
    
  END SUBROUTINE schur_solve_internal
  
  SUBROUTINE build_internal_local_matrices
    IMPLICIT NONE
    INTEGER :: i,j,el
    DOUBLE PRECISION ::
    DOUBLE PRECISION, DIMENSION(3*(N-1),3*(N-1),numelm) :: local_out
    
    DO el=1,numelm
      DO i=1,
      
  SUBROUTINE create_schur_maps
! essentially removes the dirichlet nodes and creates a new boundary-only map to use in our schur decomposition.
! any dirichlet nodes are mapped to the last+1 node of the matrix (which will exist but will not be passed to the solvers)
    IMPLICIT NONE
    INTEGER :: i,counter_x,counter_y
    counter_x=1
    counter_y=1
    DO i=1,npedg
      IF (bdflag(1,i)) THEN
	schur_bd_map_x(i) = numrows_x+numrows_y+1
      ELSE
	schur_bd_map_x(i) = counter_x
	counter_x=counter_x+1
      ENDIF
      IF (bdflag(2,i)) THEN
	schur_bd_map_y(i) = numrows_x+numrows_y+1
      ELSE
	schur_bd_map_y(i) = counter_y
	counter_y=counter_y+1
      ENDIF
    ENDDO
      
  END SUBROUTINE create_schur_maps
      
      
    

END MODULE schur_complement_module
