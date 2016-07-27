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

MODULE mkl_solver
  USE shared_data
  IMPLICIT NONE
  CONTAINS
  
  SUBROUTINE assignGlobalDimension
    IMPLICIT NONE
    INTEGER :: i,nodecount
! Count global dimensions for x and y components
! note we have a square matrix, so rows/cols are same
    numrows_x=nptot
    numrows_y=nptot
    DO i=1,npedg
      IF (bdflag(1,i)) numrows_x=numrows_x-1
      IF (bdflag(2,i)) numrows_y=numrows_y-1
    ENDDO
    global_dim = numrows_x+numrows_y+npint
! Allocate arrays for global matrix, RHS and solution
! these are the correct size for solving - ie they have all boundary nodes removed
    ALLOCATE(global_matrix(1:global_dim,1:global_dim))
    ALLOCATE(globalRHS(global_dim),globalSOL(global_dim))
    globalRHS=0d0
    globalSOL=0d0
    global_matrix=0d0
    
!     ALLOCATE(globalNodeNumberX(1:numrows_x),globalNodeNumberY(1:numrows_Y))
!     globalNodeNumberX=0
!     globalNodeNumberY=0
!     nodecount=0
!     DO i=1,npedg
!       IF (bdflag(1,i)) CYCLE
!       nodecount=nodecount+1
!       globalNodeNumberX(nodecount)=i
!     ENDDO
!     DO i=npedg+1,nptot
!       nodecount=nodecount+1
!       globalNodeNumberX(nodecount)=i
!     ENDDO 
!     nodecount=0
!     DO i=1,npedg
!       IF (bdflag(2,i)) CYCLE
!       nodecount=nodecount+1
!       globalNodeNumberY(nodecount)=i
!     ENDDO
!     DO i=npedg+1,nptot
!       nodecount=nodecount+1
!       globalNodeNumberY(nodecount)=i
!     ENDDO
   
    
  END SUBROUTINE assignGlobalDimension
  
  SUBROUTINE buildGlobalMatrix
  IMPLICIT NONE
  INTEGER :: el,i,j,ij,k,l,kl,intij,rowcount,colcount,tempint
  DOUBLE PRECISION, DIMENSION(2*nptot+npint,2*nptot+npint):: tempglob
  
  ! Assuming global matrix is already allocated

! Fill in temporary global matrix including the boundary nodes
  tempglob=0d0
  DO el=1,numelm
    DO kl=0,(N+1)**2-1
      DO ij=0,(N+1)**2-1
	tempglob(mapg(ij,el),mapg(kl,el)) = tempglob(mapg(ij,el),mapg(kl,el)) + A_x(ij,kl,el)
     
	tempglob(nptot+mapg(ij,el),nptot+mapg(kl,el)) = tempglob(nptot+mapg(ij,el),nptot+mapg(kl,el)) +A_y(ij,kl,el)
      ENDDO
    ENDDO
    DO j=1,N-1
    DO i=1,N-1
      ij=i+j*(N+1)
      intij=i+(j-1)*(N-1)
      DO kl=0,(N+1)**2-1
	tempglob(mapg(ij,el)+nptot+npint,mapg(kl,el)) = tempglob(mapg(ij,el)+nptot+npint,mapg(kl,el)) + B_x(intij,kl,el)

	tempglob(mapg(ij,el)+nptot+npint,nptot+mapg(kl,el)) = tempglob(mapg(ij,el)+nptot+npint,nptot+mapg(kl,el)) + B_y(intij,kl,el)
	  
	tempglob(mapg(kl,el),mapg(ij,el)+nptot+npint) = tempglob(mapg(kl,el),mapg(ij,el)+nptot+npint) + B_x(intij,kl,el)
	  
	tempglob(nptot+mapg(kl,el),mapg(ij,el)+nptot+npint) = tempglob(nptot+mapg(kl,el),mapg(ij,el)+nptot+npint) + B_y(intij,kl,el)
      ENDDO
    ENDDO
    ENDDO
  ENDDO
 
!
! Now fill the final global matrix excluding boundary nodes
! This is done by keeping a count of rows copied into the final global array
! By going through each section, first A_x, then A_y and removing rows, we can delete the required rows, and for each row we also delete the required column
!
  rowcount=0
  
! x-component dirichlet nodes removed 
  DO i=1,nptot
    IF (bdflag(1,i)) CYCLE
    rowcount=rowcount+1
    colcount=0
! x-component dirichlet nodes removed
    DO j=1,nptot
      IF (bdflag(1,j)) CYCLE
      colcount=colcount+1
      global_matrix(rowcount,colcount)=tempglob(i,j)
    ENDDO
! y-component dirichlet nodes removed
    DO j=nptot+1,2*nptot
      IF (bdflag(2,j-nptot)) CYCLE
      colcount=colcount+1
      global_matrix(rowcount,colcount)=tempglob(i,j)
    ENDDO
! include all internal nodes
    DO j=2*nptot+1,2*nptot+npint
      colcount=colcount+1
      global_matrix(rowcount,colcount)=tempglob(i,j)
    ENDDO
    ! Check for correct number of columns
    IF (colcount.ne.global_dim) THEN
    print*,'ERROR 1: Dimension error in global matrix columns!. global dimension:',global_dim,'column count is:',colcount
    STOP
    ENDIF
  ENDDO
! y-component dirichlet nodes removed
  DO i=nptot+1,2*nptot
    IF (bdflag(2,i-nptot)) CYCLE
    rowcount=rowcount+1
    colcount=0
! x-component dirichlet nodes removed
    DO j=1,nptot
      IF (bdflag(1,j)) CYCLE
      colcount=colcount+1
      global_matrix(rowcount,colcount)=tempglob(i,j)
    ENDDO
! y-component dirichlet nodes removed
    DO j=nptot+1,2*nptot
      IF (bdflag(2,j-nptot)) CYCLE
      colcount=colcount+1
      global_matrix(rowcount,colcount)=tempglob(i,j)
    ENDDO
    DO j=2*nptot+1,2*nptot+npint
      colcount=colcount+1
      global_matrix(rowcount,colcount)=tempglob(i,j)
    ENDDO
    ! Check for correct number of columns
    IF (colcount.ne.global_dim) THEN
    print*,'ERROR 2: Dimension error in global matrix columns!. global dimension:',global_dim,'column count is:',colcount
    STOP
    ENDIF
  ENDDO
  DO i=2*nptot+1,2*nptot+npint
!     IF (bdflag(1,i)) CYCLE
    rowcount=rowcount+1
    colcount=0
    DO j=1,nptot
      IF (bdflag(1,j)) CYCLE
      colcount=colcount+1
      global_matrix(rowcount,colcount)=tempglob(i,j)
    ENDDO      
    DO j=nptot+1,2*nptot
      IF (bdflag(2,j-nptot)) CYCLE
      colcount=colcount+1
      global_matrix(rowcount,colcount)=tempglob(i,j)
    ENDDO
    DO j=2*nptot+1,2*nptot+npint
      colcount=colcount+1
      global_matrix(rowcount,colcount)=tempglob(i,j)
    ENDDO
! Check for correct number of columns
    IF (colcount.ne.global_dim) THEN
    	print*,'ERROR 3: Dimension error in global matrix columns!. global dimension:',global_dim,'column count is:',colcount
    	print*,nptot,npint,i,rowcount
    STOP
    ENDIF
  ENDDO
! Check for correct number of rows
    IF (rowcount.ne.global_dim) THEN
    	print*,'ERROR: Dimension error in global matrix rows! global dimension:',global_dim,'row count is:',rowcount
    STOP
    ENDIF
	
  END SUBROUTINE buildGlobalMatrix
  
  SUBROUTINE build_global_sparse
!
! Builds global sparse matrix
!
  IMPLICIT NONE
  INTEGER :: i,j,k,l,rowcount,first

! Assuming we have our final global matrix.
! Next need to count the number of non-zero elements
  nonz=0
  DO i=1,global_dim
    DO j=1,global_dim ! <--- CHANGE FOR SYMMETRIC MATRIX !
      IF (abs(global_matrix(i,j)).gt.1d-15) nonz=nonz+1
    ENDDO
  ENDDO

! Assign memory for sparse storage:
  ALLOCATE(values(nonz),columns(nonz),rowIndex(global_dim+1))

! Putting into sparse form (NON-SYMMETRIC):
  i = 1
  DO k = 1, global_dim
    first = 0
    DO l = 1, global_dim
      IF(abs(global_matrix(k,l)).gt.1d-15) THEN
	values(i) = global_matrix(k,l) 
	columns(i) = l
	IF(first.eq.0) THEN
	  rowIndex(k) = i
	  first = 1
	ENDIF
	i = i + 1
      ENDIF
    ENDDO
  ENDDO
  
! Putting into sparse form (SYMMETRIC):
!   i = 1
!   DO k = 1, global_dim
!     first = 0
!     DO l = k, global_dim
!       IF(abs(global_matrix(k,l)).gt.1d-15) THEN
! 	values(i) = global_matrix(k,l) 
! 	columns(i) = l
! 	IF(first.eq.0) THEN
! 	  rowIndex(k) = i
! 	  first = 1
! 	ENDIF
! 	i = i + 1
!       ENDIF
!     ENDDO
!   ENDDO  

  rowIndex(global_dim+1) = nonz + 1
  
! Now have global matrix stored in sparse form
! Fill globalRHS - the known vector in Ax=b

  rowcount=0
  DO i=1,nptot
    IF (bdflag(1,i)) CYCLE
    rowcount=rowcount+1
    globalRHS(rowcount)=f_x(i)
  ENDDO
  DO i=1,nptot
    IF (bdflag(2,i)) CYCLE
    rowcount=rowcount+1
    globalRHS(rowcount)=f_y(i)
  ENDDO
  DO i=1,npint
    globalRHS(rowcount+i)=g(i)
  ENDDO
  
!   OPEN(99,FILE='values.dat')
!   DO i=1,nonz
!     write(99,*) values(i)
!   ENDDO
!   close(99)
!   OPEN(99,FILE='columns.dat')
!   DO i=1,nonz
!     write(99,*) columns(i)
!   ENDDO
!   close(99)
!   OPEN(99,FILE='rowIndex.dat')
!   DO i=1,global_dim+1
!     write(99,*) rowIndex(i)
!   ENDDO
!   close(99)
!   
!   open(99,FILE='globalmatrix.dat')
!   DO i=1,global_dim
!     write(99,*) (global_matrix(i,j),j=1,global_dim)
!   ENDDO
!   close(99)
!   open(99,FILE='globalRHS.dat')
!   DO i=1,global_dim
!     write(99,*) globalRHS(i)
!   ENDDO
!   close(99)
!   
!   OPEN(99,FILE='xdata.dat')
!   DO i=1,nptot
!   	IF (bdflag(1,i)) CYCLE
!   	write(99,*) nodeCoord(i,1)
!   ENDDO
!   CLOSE(99)
!   
!   OPEN(99,FILE='ydata.dat')
!   DO i=1,nptot	
!   	IF (bdflag(2,i)) CYCLE
!   	write(99,*) nodeCoord(i,2)
!   ENDDO
!   CLOSE(99)
!   stop
!   print*,'sparse matrix done'

  END SUBROUTINE build_global_sparse  
  
  SUBROUTINE Pardiso_Fact(pt,iparm)!, Nonz, Gmax, values,columns,rowIndex)
    implicit none

    integer*8  pt(64)
    integer :: maxfct, mnum, mtype, phase, nrhs, error, msglvl, i
    integer :: solver, idum, Gmax!, Nonz
    integer :: iparm(64)!, rowIndex(Gmax+1),columns(Nonz)
 
!  double precision :: values(Nonz)
    double precision :: dparm(64), ddum
    Gmax=global_dim
    pt = 0 !internal pointer - MUST NOT BE CHANGED FROM HERE ON IN
    nrhs = 1 ! number of rhs vectors
    mnum = 1 ! number of matrices
    maxfct = 1 ! max number of matrices (i think)

   mtype = 11 !real unsymmetric matrix 
!     mtype = -2 !for real symmetric INDEFINITE matrix (2 = positive definite)

    phase = 12 !perform analysis and numerical factorisation

    msglvl = 0 !print no (=0) statistical information on screen 

    iparm(1) = 0 !use default solver settings
!     iparm(3) = 1
!     iparm(27)=1 !think this checks matrix - check the pardiso table

!   idum = 0 !dummy variables which aren't used for this particular phase=12
    ddum = 0

!   print*,'entering paradiso_fact stage'
    CALL pardiso (pt, maxfct, mnum, mtype, phase, Gmax, values, rowIndex, columns,& 
                 idum, nrhs, iparm, msglvl, ddum, ddum, error) 

!   write(*,*) 'error from pardiso_fact',error
    !STOP
  END SUBROUTINE Pardiso_Fact

  SUBROUTINE Pardiso_Solve(pt,iparm)!, Nonz, Gmax, values,columns,rowIndex,VGlobal)
    implicit none
    integer*8 ::  pt(64)
    integer :: maxfct, mnum, mtype, phase, nrhs, error, msglvl, i
    integer :: solver, idum, Gmax!, Nonz
    integer :: iparm(64)!, rowIndex(Gmax+1),columns(Nonz)
    double precision :: dparm(64), ddum

    Gmax=global_dim

    nrhs = 1 ! number of rhs vectors
    mnum = 1 ! number of matrices
    maxfct = 1 ! max number of matrices (i think)

    mtype = 11 !real unsymmetric matrix 
!      mtype = -2 !for real symmetric INDEFINITE matrix (2 = positive definite)

    phase = 33 !solves equation

    msglvl = 0 !print no (=0) statistical information on screen 


    iparm(1) = 0 !use default solver settings
!   iparm(3) = 1
!   iparm(6) = 0 !write solution onto vglobal

!   idum = 0 !dummy variables which aren't used for this particular phase=12
    ddum = 0


    do i = 1, Gmax
      if((globalRHS(i)+1d0).eq.globalRHS(i))then
	write(*,*) 'globalRHS is NaN, before solve sym',i
	stop
      endif
    enddo


! print*,'entering paradiso_solve stage'
    call pardiso (pt, maxfct, mnum, mtype, phase, Gmax, values, rowIndex, columns,& 
                 idum, nrhs, iparm, msglvl, globalRHS, globalSOL, error) 

    do i = 1, Gmax
      if((globalRHS(i)+1d0).eq.globalRHS(i))then
	write(*,*) 'Vglobal is NaN, after solve sym',i
	stop
      endif
    enddo


!    write(*,*) 'error from pardiso_solve',error

  END SUBROUTINE Pardiso_Solve
  
  SUBROUTINE Pardiso_Release_Memory(pt,iparm)!,Gmax)
    implicit none

    integer*8 :: pt(64)
    integer :: maxfct, mnum, mtype, phase, nrhs, error, msglvl, i
    integer :: solver, idum, Gmax, Nonz
    integer :: iparm(64)
 
    double precision :: dparm(64), ddum
    
    Gmax=global_dim
    nrhs = 1 ! number of rhs vectors
    mnum = 1 ! number of matrices
    maxfct = 1 ! max number of matrices (i think)

    mtype = 11 !real unsymmetric matrix 
!      mtype = -2 !for real symmetric INDEFINITE matrix (2 = positive definite)

    phase = -1 !release memory

    msglvl = 0 !print no (=0) statistical information on screen 

    iparm(1) = 0 !use default solver settings
   iparm(3) = 0 ! was 1

   idum = 0 !dummy variables which aren't used for this particular phase
    ddum = 0

! print*,'entering paradiso memory release stage'
   call pardiso(pt, maxfct, mnum, mtype, phase, Gmax, ddum, idum, idum,& 
                idum, nrhs, iparm, msglvl, ddum, ddum, error) 


!    write(*,*) 'error from pardiso_release_memory',error

  END SUBROUTINE Pardiso_Release_Memory
 
  SUBROUTINE RUN_PARDISO
    IMPLICIT NONE
  
    integer*8 ::  pt(64)
    integer :: iparm(64),i
! factorisation
    CALL Pardiso_Fact(pt,iparm)!nonz,global_dim,values,columns,rowIndex)
! solve
    CALL Pardiso_Solve(pt,iparm)!, Nonz, global_dim, values,columns,rowIndex,VGlobal)
! Release memory
    CALL Pardiso_Release_Memory(pt,iparm)
    DEALLOCATE(values,columns,rowIndex)
!   print*,'Paradiso solver done'

  END SUBROUTINE RUN_PARDISO
  
  SUBROUTINE PARDISO_UPDATE_SOLUTION
    IMPLICIT NONE
    INTEGER :: rowcount,i
! Transfer global solution to V_x, V_y, p
! release global memory.
    
    rowcount=0
    DO i=1,nptot
      IF (bdflag(1,i)) THEN
	V_x(i) = boundary_x(i)
      ELSE
	rowcount=rowcount+1
	V_x(i)= globalSOL(rowcount)
      ENDIF
    ENDDO
    DO i=1,nptot
      IF (bdflag(2,i)) THEN
	V_y(i) = boundary_y(i)
      ELSE
	rowcount=rowcount+1
	V_y(i)= globalSOL(rowcount)
      ENDIF
    ENDDO
    DO i=1,npint
      rowcount=rowcount+1
      pressure(i)= globalSOL(rowcount)
    ENDDO
! Release memory: not needed right now, as global dimension remains constant from beginning.
!    DEALLOCATE(global_matrix,globalRHS,globalSOL)
    
    
  END SUBROUTINE PARDISO_UPDATE_SOLUTION
  
  SUBROUTINE PARDISO_UPDATE_RHS
    IMPLICIT NONE
    INTEGER :: rowcount,i
    
! Transfer to the globalRHS vector for pardiso:
    rowcount=0
    DO i=1,nptot
      IF (bdflag(1,i)) CYCLE
      rowcount=rowcount+1
      globalRHS(rowcount)=f_x(i)
    ENDDO
    DO i=1,nptot
      IF (bdflag(2,i)) CYCLE
      rowcount=rowcount+1
      globalRHS(rowcount)=f_y(i)
    ENDDO
    DO i=1,npint
      globalRHS(rowcount+i)=g(i)
    ENDDO
  END SUBROUTINE PARDISO_UPDATE_RHS

END MODULE mkl_solver

