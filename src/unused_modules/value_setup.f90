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

MODULE value_setup
  USE constants
  USE shared_data
  USE utils
  USE known_terms
  USE functions_module

  IMPLICIT NONE
  CONTAINS

! Old module, kept only for the known solution testing part - will not be used at compile time in most cases.




  

  




! THIS FUNCTION IS NOW REDUNDANT.. may want to use for validation later
  SUBROUTINE construct_velocity_ACT

  IMPLICIT NONE
  INTEGER i,j,k,ij
  DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1) :: uactloc
  DOUBLE PRECISION, DIMENSION(nptot) :: temp

  V_xACT=0d0
    DO k=1,numelm
      temp=0
      DO j=0,N
      DO i=0,N
	ij=i+j*(N+1)
        uactloc(ij)=uterm_x(nodeCoord(mapg(ij,k),1),nodeCoord(mapg(ij,k),2))
      ENDDO
      ENDDO
      temp=0
      CALL vecglobalprolongationUACT(uactloc,k,temp)
      V_xACT = V_xACT + temp
    ENDDO
    
    V_yACT=0d0
    DO k=1,numelm
      temp=0
      DO j=0,N
      DO i=0,N
	ij=i+j*(N+1)
        uactloc(ij)=uterm_y(nodeCoord(mapg(ij,k),1),nodeCoord(mapg(ij,k),2))	
      ENDDO
      ENDDO
      temp=0d0
      CALL vecglobalprolongationUACT(uactloc,k,temp)
      V_yACT = V_yACT + temp
    ENDDO

  END SUBROUTINE construct_velocity_act

! THIS FUNCTION IS NOW REDUNDANT.. may want to use for validation later
  SUBROUTINE construct_pressure_ACT

  IMPLICIT NONE
  INTEGER i,j,k,ij
  DOUBLE PRECISION, DIMENSION((N-1)**2) :: pactloc
  DOUBLE PRECISION, DIMENSION(npint) :: temp

  pressureACT=0d0
    DO k=1,numelm
      temp=0
      DO i=1,N-1
      DO j=1,N-1
	ij=i+j*(N+1)
        pactloc(i+(j-1)*(N-1))=pressureterm(nodeCoord(mapg(ij,k),1),nodeCoord(mapg(ij,k),2))!*w(i)*w(j)*jac(k,i,j)   
      ENDDO
      ENDDO
      temp=0d0
      CALL vecglobalprolongation_internal_nodes(pactloc,k,temp)
      pressureACT = pressureACT + temp
    ENDDO

  END SUBROUTINE construct_pressure_act


END MODULE value_setup