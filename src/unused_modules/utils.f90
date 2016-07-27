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

MODULE utils
  IMPLICIT NONE 

  CONTAINS

  SUBROUTINE printvec(vec,title,M)
    IMPLICIT NONE
    INTEGER :: i,M
    CHARACTER*(*) :: title
    DOUBLE PRECISION, DIMENSION(0:M) :: vec
  
    write(*,*) title
    DO i=0,M
      write(*,*) i,vec(i)
    ENDDO
    
  END SUBROUTINE printvec
  ! CHECK - needs formatting in write statement.
  SUBROUTINE printmat(mat,title,M1,M2)
    IMPLICIT NONE
    INTEGER :: i,j,M1,M2
    CHARACTER*(*) :: title
    DOUBLE PRECISION, DIMENSION(0:M1,0:M2) :: mat
    
    write(*,*) title
    DO i=0,M1
      write(*,*) i,(mat(i,j),j=0,M2)
    ENDDO
  END SUBROUTINE printmat

END MODULE utils