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

MODULE solver
  USE shared_data
!  USE utils
  USE functions_module

  IMPLICIT NONE

  CONTAINS


  SUBROUTINE minressolve(U_in,V_in,p_in,RHSx_in,RHSy_in,RHSp_in,minresconv)
! solving Mx = f where M is a matrix decomposed into the stokes problem with A(laplacian) and B(gradient),
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
! vectors with 1 represent the "velocity part" of the vector (of dimenion nptot)
! vectors with 2 represent the "pressure part" of the vector (of dimension numelm*(N-1)**2=npint
    DOUBLE PRECISION, DIMENSION(nptot) :: U_in,V_in,RHSx_in,RHSy_in,&
					  v1xold,v1yold,v1xnew,v1ynew,z1xold,z1yold,z1xnew,z1ynew,&
					  w1xold,w1yold,w1xnew,w1ynew,tempx1,tempy1,Ax,Ay,x1,y1,&
					  xstore1,ystore1
    DOUBLE PRECISION, DIMENSION(npint) :: p_in,RHSp_in,&
					  v2old,v2new,z2old,z2new,w2old,w2new,temp2,x2,y2,Bx,By,Zp
    DOUBLE PRECISION, DIMENSION(2*nptot+npint) :: tempdiag
  
    minresconv=0 ! parameter to feedback if it converged.
! x1,y1 and x2,y2 are used for storing matrix vector multiplications - saves on operations.

 
! Use v1xold, v1yold and v2old to store the values of our matrix times the initial guess for velocity(x & y) and pressure
!   restartcount=0
    CALL Bmvmult(B_x,U_in,Bx,1) ! Bx = B*velocityguess (x)
    CALL Bmvmult(B_y,V_in,By,2) ! By = B*velocityguess (y)
    CALL Amvmult(A_x,U_in,Ax,1) ! Ax = A*velocityguess (x)
    CALL Amvmult(A_y,V_in,Ay,2) ! Ay = A*velocityguess (y)
    CALL B_transpose_mvmult(B_x,p_in,x1) !x1 = B'*pressureguess
    CALL B_transpose_mvmult(B_y,p_in,y1) !y1 = B'*pressureguess
    CALL Z_pmvmult(p_in,Zp) !Zp = Z*pressureguess
    v1xold = Ax + x1
    v1yold = Ay + y1
    v2old = Bx + By! - Zp

    v1xnew=0d0
    v1ynew=0d0
    v2new=0d0
    DO i=1,nptot
!       IF (bdflag(1,i)) CYCLE
	v1xnew(i)=RHSx_in(i) - v1xold(i)
    ENDDO
    DO i=1,nptot
!       IF (bdflag(2,i)) CYCLE
      v1ynew(i)=RHSy_in(i) - v1yold(i)
    ENDDO
    v2new = RHSp_in - v2old
  
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
!   CALL preconcgsolve(v1xnew,v1ynew,v2new,z1xnew,z1ynew,z2new)

! Divide by reciprical of diag(A_x/y)
      tempdiag=0d0
      CALL calc_reciprical_of_diagA(tempdiag)
      DO i=1,nptot
! 	IF (bdflag(1,i)) CYCLE
	z1xnew(i)=v1xnew(i)/tempdiag(i)
      ENDDO
      DO i=1,nptot
! 	IF (bdflag(2,i)) CYCLE
	z1ynew(i)=v1ynew(i)/tempdiag(nptot+i)
      ENDDO
      DO i=1,npint
	z2new(i)=v2new(i)/tempdiag(2*nptot+i)
      ENDDO
! END PRECON STEP !
    ELSE
!   print*, 'No Precon'
      z1xnew = v1xnew
      z1ynew = v1ynew
      z2new = v2new
    ENDIF
  

! scalar data setup
  gammanew = SQRT(bddotp(z1xnew,v1xnew,1) + bddotp(z1ynew,v1ynew,2) + dotp(z2new,v2new,npint-1))
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
    
    CALL Bmvmult(B_x,z1xnew,Bx,1) ! Bx = B*z1xnew
    CALL Bmvmult(B_y,z1ynew,By,2) ! By = B*z1ynew
    CALL Amvmult(A_x,z1xnew,Ax,1) ! Ax = A*z1xnew
    CALL Amvmult(A_y,z1ynew,Ay,2) ! Ay = A*z1ynew
    CALL B_transpose_mvmult(B_x,z2new,x1) !x1 = B'*z2new
    CALL B_transpose_mvmult(B_y,z2new,y1) !y1 = B'*z2new
    CALL Z_pmvmult(z2new,Zp) ! Zp = Z*z2new


    xstore1 = Ax + x1
    ystore1 = Ay + y1
    x2 = Bx + By! - Zp
    

! Essentially we have 	[Ax + x1] = [ A 0 B'][z1xnew] (=xstore1)
!		        [Ay + y1]   [ 0 A B'][z1ynew] (=ystore1)
!			[Bx + By]   [ B B C][z2new ] (=x2)

    delta = bddotp(z1xnew,xstore1,1) + bddotp(z1ynew,ystore1,2) + dotp(z2new,x2,npint-1) ! Delta = v'Mv where M is the whole stokes matrix, and v the vector for both vel and pressure.
!     delta = dotp(v1xnew,xstore1,nptot-1) + dotp(v1ynew,ystore1,nptot-1) + dotp(v2new,x2,npint-1)
    DO i=1,nptot
!       IF (bdflag(1,i)) CYCLE
      v1xnew(i) = xstore1(i) - delta*v1xnew(i)/gammanew - gammanew*v1xold(i)/gammaold
    ENDDO
    DO i=1,nptot
!       IF (bdflag(2,i)) CYCLE
      v1ynew(i) = ystore1(i) - delta*v1ynew(i)/gammanew - gammanew*v1yold(i)/gammaold
    ENDDO
    DO i=1,npint
      v2new(i) = x2(i) - delta*v2new(i)/gammanew - gammanew*v2old(i)/gammaold
    ENDDO
    v1xold=tempx1
    v1yold=tempy1
    v2old=temp2

    z1xold = z1xnew
    z1yold = z1ynew
    z2old = z2new
    IF (preconflag.ne.0) THEN
! PRECON STEP !
!       CALL preconcgsolve(v1xnew,v1ynew,v2new,z1xnew,z1ynew,z2new)
    
    DO i=1,nptot
!       IF (bdflag(1,i)) CYCLE
      z1xnew(i)=v1xnew(i)/tempdiag(i)
    ENDDO
    DO i=1,nptot
!       IF (bdflag(2,i)) CYCLE
      z1ynew(i)=v1ynew(i)/tempdiag(nptot+i)
    ENDDO
    DO i=1,npint
      z2new(i)=v2new(i)/tempdiag(2*nptot+i)
    ENDDO
! END PRECON STEP !
    ELSE
      z1xnew = v1xnew
      z1ynew = v1ynew
      z2new = v2new
    ENDIF
    
    gammaold=gammanew
    gammanew=SQRT(bddotp(z1xnew,v1xnew,1) + bddotp(z1ynew,v1ynew,2) + dotp(z2new,v2new,npint-1))
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
    DO i=1,nptot
!       IF (bdflag(1,i)) CYCLE
      w1xnew(i) = (z1xold(i) - alpha3*w1xold(i) - alpha2*w1xnew(i))/alpha1
    ENDDO
    DO i=1,nptot
!       IF (bdflag(2,i)) CYCLE
      w1ynew(i) = (z1yold(i) - alpha3*w1yold(i) - alpha2*w1ynew(i))/alpha1
    ENDDO
    w2new = (z2old - alpha3*w2old - alpha2*w2new)/alpha1
    
    w1xold = tempx1
    w1yold = tempy1
    w2old = temp2
    DO i=1,nptot
!       IF (bdflag(1,i)) CYCLE
      U_in(i) = U_in(i) + cnew*eta*w1xnew(i)
    ENDDO
    DO i=1,nptot
!       IF (bdflag(2,i)) CYCLE
      V_in(i) = V_in(i) + cnew*eta*w1ynew(i)
    ENDDO
    p_in=p_in+cnew*eta*w2new
    eta=-snew*eta
    
! calc norm of residual for current step
    normr = abs(eta)
	!print*,gammanew,normr
!     print*,'normr is',normr,'alpha0 is',alpha0
! Check for convergence of residual
    IF (abs(normr).lt.1d-14) THEN

!       print*,'converged in ',counter,'iterations with normr',abs(normr)
      minresconv=counter
      EXIT
    ENDIF
! Check if method is stagnating
    IF (abs(alpha0).lt.1d-14) THEN
!      print*,alpha0,'warning stagnation',stag
      stag = stag + 1
      IF (stag.gt.3) THEN ! 3 iterations are same - stop!
	print *, 'MINRES method is stagnant after',counter,' iterations'
	STOP
      ENDIF
    ELSE
      stag = 0
    ENDIF
!     IF (counter.eq.20) THEN
!       print*,'restarting minres normr is ', abs(normr)
!       restartcount = restartcount+1
!       IF (restartcount.gt.1000) THEN
!             print *,'MINRES failed after ', counter,' iterations. Normr was ',abs(normr)
!       STOP
!      ELSE
!  	GO TO 1001
!      ENDIF
!     ENDIF

    IF (counter.gt.50000) THEN
!         print *,'MINRES failed after ', counter,' iterations. Normr was ',abs(normr)
      minresconv=counter
      EXIT
    ENDIF
    
    DO i=1,nptot
      IF (IsNaN(V_x(i))) THEN
	print*, 'ERROR! V_x(',i,') is NaN'
	EXIT
      ELSEIF (IsNaN(V_y(i))) THEN
	print*, 'ERROR! V_y(',i,') is NaN'
	EXIT
      ENDIF
    ENDDO
    DO i=1,npint
      IF (IsNaN(pressure(i))) THEN
	print*, 'ERROR! pressure(',i,') is NaN'
	EXIT
      ENDIF
    ENDDO
  ENDDO
  
END SUBROUTINE minressolve


 
  
  
  SUBROUTINE mmmult(A,M1,N1,B,M2,N2,C)
  
! multiplies matrices A and B, output is C
! A is of size (M1 + 1) x (N1 + 1)
! B is of size (M2 + 1) x (N2 + 1)
! C is of size (M1 + 1) x (N2 + 1)
  
    IMPLICIT NONE
    INTEGER :: i,j,k,M1,M2,N1,N2,ij,kj,ik
    DOUBLE PRECISION :: temp
    DOUBLE PRECISION, DIMENSION(0:M1,0:N1) :: A ! was (0:M1,0:N1)
    DOUBLE PRECISION, DIMENSION(0:M2,0:N2) :: B ! was (0:M2,0:N2)
    DOUBLE PRECISION, DIMENSION(0:M1,0:N2) :: C ! was (0:M1,0:N2)
    
    IF (N1.ne.M2) THEN
      write(*,*) 'ERROR: Size of matrices not consistent, returned C=0'
      C = 0d0
    ELSE
      DO i=0,M1
        DO j=0,N2
          temp = 0d0
          DO k=0,N1
            temp = temp + A(i,k)*B(k,j)
          ENDDO
          C(i,j) = temp
        ENDDO
      ENDDO
    ENDIF
  END SUBROUTINE mmmult
  
  SUBROUTINE Amvmult(Amatrix,x,y,c)
!
! Takes a global vector x, and multiplies it by the global stiffness matrix
! However, all multiplication is done at a local level, using the prolongation and restriction operators
! Dirichlet boundaries are removed via the bdflag array.
! c determines if we are using the x or y component
  
  IMPLICIT NONE
  INTEGER :: ij,kl,el,c
  DOUBLE PRECISION, DIMENSION(nptot) :: x,y,tempvec
  DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1) :: xk,yk
  DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1,0:(N+1)**2-1,1:numelm) :: Amatrix

  y=0d0
  DO el=1,numelm
!     xk=0d0
    CALL veclocalrestrict(x,el,xk)
! perform local matrix-vector multiplication
    yk=0d0
    DO kl=0,(N+1)**2-1
!      IF (bdflag(c,mapg(ij,el))) CYCLE
      DO ij=0,(N+1)**2-1
!	IF (bdflag(c,mapg(kl,el))) CYCLE
	  yk(ij) = yk(ij) + Amatrix(ij,kl,el)*xk(kl)
      ENDDO
    ENDDO
    tempvec=0d0
    CALL vecglobalprolongation(yk,el,tempvec)
! y = Ak*yk +s
    y = y + tempvec
  ENDDO
  END SUBROUTINE Amvmult
  
  SUBROUTINE Bmvmult(Bmatrix,x,y,c)
    
    IMPLICIT NONE
    INTEGER :: el,ij,kl,c
    DOUBLE PRECISION, DIMENSION(npint) :: y,tempvec
    DOUBLE PRECISION, DIMENSION(nptot) :: x
    DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1) :: xk
    DOUBLE PRECISION, DIMENSION((N-1)**2) :: yk 
    DOUBLE PRECISION, DIMENSION(1:(N-1)**2,0:(N+1)**2-1,1:numelm) :: Bmatrix

    y=0d0
    DO el=1,numelm
      ! pull out local elements from vector
!       xk=0d0
      CALL veclocalrestrict(x,el,xk)
! perform local matrix-vector multiplication      
      yk=0d0
      DO kl=0,(N+1)**2-1
	DO ij=1,(N-1)**2
!	  IF (bdflag(c,mapg(kl,el))) CYCLE !needed here.
	  yk(ij) = yk(ij) + Bmatrix(ij,kl,el)*xk(kl)
	ENDDO
      ENDDO
      tempvec=0d0
      CALL vecglobalprolongation_internal_nodes(yk,el,tempvec)
! y = Ak*yk +s
      y = y + tempvec
    ENDDO
  END SUBROUTINE Bmvmult
  
  SUBROUTINE B_transpose_mvmult(Bmatrix,x,y)
    
    IMPLICIT NONE
    INTEGER :: el,ij,kl
    DOUBLE PRECISION, DIMENSION(npint) :: x
    DOUBLE PRECISION, DIMENSION(nptot) :: y,tempvec
    DOUBLE PRECISION, DIMENSION((N-1)**2) :: xk
    DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1) :: yk 
    DOUBLE PRECISION, DIMENSION(1:(N-1)**2,0:(N+1)**2-1,1:numelm) :: Bmatrix

    y=0d0
    DO el=1,numelm
! pull out local elements from vector
!      xk=0d0
      CALL veclocalrestrict_internal_nodes(x,el,xk)
! perform local matrix-vector multiplication
      yk=0d0
      DO kl=0,(N+1)**2-1
	DO ij=1,(N-1)**2
	  yk(kl) = yk(kl) + Bmatrix(ij,kl,el)*xk(ij) 
	ENDDO
      ENDDO
      tempvec=0d0
      CALL vecglobalprolongation(yk,el,tempvec)
! y = Ak*yk +s
      y = y + tempvec
    ENDDO
  END SUBROUTINE B_transpose_mvmult
  
  SUBROUTINE Z_pmvmult(x,y)
! Multiplies a vector x, of global size npint (pressure), by the Zero pressure integral matrix.
! This involves integrating the local part of x times the pressure test function, over the whole domain.
! We then integrate the result times the pressure test function.
! 
! Note: Z_p(kl,el) gives the integral of pressure test function q_k,l = h_k*h_l (tilde) over the element el.
! SEE SEM_module for the routine where this is calculated.
    
    IMPLICIT NONE
    INTEGER :: el,ij,kl
    DOUBLE PRECISION :: temp
    DOUBLE PRECISION, DIMENSION(npint) :: x,y,tempvec
    DOUBLE PRECISION, DIMENSION(1:(N-1)**2) :: xk,yk
    

    temp=0d0
    DO el=1,numelm
! pull out local elements from vector
!       xk=0d0
      CALL veclocalrestrict_internal_nodes(x,el,xk)
! perform local matrix-vector multiplication
      DO ij=1,(N-1)**2
	temp=temp+xk(ij)*Z_p(ij,el)
      ENDDO
    ENDDO
    
    y=0d0
    DO el=1,numelm
!       yk=0d0
      DO kl=1,(N-1)**2
	yk(kl) = Z_p(kl,el)*temp
      ENDDO
      tempvec=0d0
      CALL vecglobalprolongation_internal_nodes(yk,el,tempvec)
! y = Ak*yk +s
      y = y + tempvec
    ENDDO
  END SUBROUTINE Z_pmvmult
  
!   SUBROUTINE Z_pmvmult(x,y)
!     
!     IMPLICIT NONE
!     INTEGER :: el,ij,kl
!     DOUBLE PRECISION, DIMENSION(npint) :: x,y,tempvec
!     DOUBLE PRECISION, DIMENSION(1:(N-1)**2) :: xk,yk
!     
! 
!     y=0d0
!     DO el=1,numelm
! ! pull out local elements from vector
!       xk=0d0
!       CALL veclocalrestrict_internal_nodes(x,el,xk)
! ! perform local matrix-vector multiplication
!       yk=0d0
!       DO kl=1,(N-1)**2
! 	DO ij=1,(N-1)**2
! 	  yk(kl) = yk(kl) + Z_p(kl,ij,el)*xk(ij)
! 	ENDDO
!       ENDDO
!       tempvec=0d0
!       CALL vecglobalprolongation_internal_nodes(yk,el,tempvec)
! ! y = Ak*yk +s
!       y = y + tempvec
!     ENDDO
!   END SUBROUTINE Z_pmvmult
  
  
  
  SUBROUTINE mvmult(A,x,M,y)
! performs y = A*x
! where A is a (M+1)x(M+1) matrix
! and x is a vector of length N+1

    IMPLICIT NONE
    INTEGER :: i,j,M
    DOUBLE PRECISION, DIMENSION(0:M) :: x,y
    DOUBLE PRECISION, DIMENSION(0:M,0:M) :: A
    
    DO i=0,M
      y(i) = 0d0
      DO j=0,M
	y(i) = y(i) + A(i,j)*x(j)
      ENDDO
    ENDDO
  END SUBROUTINE mvmult
  
  DOUBLE PRECISION FUNCTION dotp(x,y,M) RESULT(answer)
  
    IMPLICIT NONE
    INTEGER :: i,M
    DOUBLE PRECISION, DIMENSION(0:M) :: x,y
    
    answer=0d0
    DO i=0,M
      answer = answer + x(i)*y(i)
    ENDDO
    RETURN
  END FUNCTION dotp
  
  DOUBLE PRECISION FUNCTION bddotp(x,y,c) RESULT(answer)
! dot product of x and y restricted on the dirichlet boundary (ie, answer does not include dirichlet nodes)  
! c = component, 1=x 2=y
  
    IMPLICIT NONE
    INTEGER :: i,c
    DOUBLE PRECISION, DIMENSION(nptot) :: x,y
    
    answer=0d0
    DO i=1,nptot
!      IF (bdflag(c,i)) CYCLE
	answer = answer + x(i)*y(i)
    ENDDO
    RETURN
  END FUNCTION bddotp
  
  
  DOUBLE PRECISION FUNCTION A_IP(Amatrix,x,y,c) RESULT(answer)
! calculates Matrix inner product using stiffness matrix (restricted on dirichlet boundary)
! that is: answer = x'Ay
    IMPLICIT NONE
    INTEGER :: c
    DOUBLE PRECISION, DIMENSION(nptot) :: x,y,temp
    DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1,0:(N+1)**2-1,1:numelm) :: Amatrix
    
    CALL Amvmult(Amatrix,y,temp,c) ! temp=Ay
    answer = bddotp(x,temp,c) !x'temp
    RETURN
  END FUNCTION A_IP
  
  
!   DOUBLE PRECISION FUNCTION StokesIP(x1,x2,y1,y2) RESULT(answer)
! ! Input is in two parts of the vector you wish to mult: x1,y1 (dimension nptot) and x2,y2 (dimension npint)
!   
!     IMPLICIT NONE
!     DOUBLE PRECISION, DIMENSION(nptot) :: x1,y1,temp1,temp2
!     DOUBLE PRECISION, DIMENSION(npint) :: x2,y2,temp3
!     
! ! these two make up the nptot sized part of the vector results from the stokes matrix multiplied by [x1;y1]
!  
!     CALL stiffmvmult(y1,temp1) !temp1 = A*y1
!     print *,'STIFF MATRIX DONE'
!     CALL Bmatrix_transpose_mvmult(y2,temp2) !temp2 = B'*y2
!     print *,'B TRANS Matrix done'
! ! This makes up the npint sized part of the vector
!     CALL Bmatrixmvmult(y1,temp3) !temp3 = B*y1
!     print *,'BMatrix done'
! 
!     
! ! Now form the dot product of the two
!     answer = bddotp(x1,temp1+temp3) + dotp(x2,temp2,npint)
!     RETURN
!   END FUNCTION StokesIP 
    
  
  DOUBLE PRECISION FUNCTION Aip(A,x,y,M) RESULT(answer)
! calculates the A-Inner Product <x,y> = x'Ay
! where x and y are vectors, and A is some matrix
    IMPLICIT NONE
    INTEGER :: M
    DOUBLE PRECISION, DIMENSION(0:M) :: x,y,temp
    DOUBLE PRECISION, DIMENSION(0:M,0:M) :: A

!calc Ay
    CALL mvmult(A,y,M,temp)
! calc x'(Ay) = x'temp
    answer = dotp(x,temp,M)
    RETURN
  END FUNCTION Aip
  
  
  SUBROUTINE calc_reciprical_of_diagA(diag)
    IMPLICIT NONE
    INTEGER :: el,ij,i,j,intij
    DOUBLE PRECISION, DIMENSION(1:2*nptot+npint) :: diag
! Calulates the matrix with entries 1/A(i,i) (avoiding the zeros from deleted dirichlet nodes)
    diag=0d0
    DO el=1,numelm
      DO ij=0,(N+1)**2-1
	IF (bdflag(1,mapg(ij,el))) THEN
	  diag(mapg(ij,el))=1d0
	ELSE
	  diag(mapg(ij,el))=diag(mapg(ij,el))+A_x(ij,ij,el)
	ENDIF
      ENDDO
    ENDDO
    DO el=1,numelm
      DO ij=0,(N+1)**2-1
	IF (bdflag(2,mapg(ij,el))) THEN 
	  diag(nptot+mapg(ij,el))=1d0
	ELSE
	  diag(nptot+mapg(ij,el))=diag(nptot+mapg(ij,el))+A_y(ij,ij,el)
	ENDIF
      ENDDO
    ENDDO
    DO el=1,numelm
      DO j=1,N-1
      DO i=1,N-1
	ij=i+j*(N+1)
	intij=i+(j-1)*(N-1)
	diag(2*nptot+mapg(ij,el)-(nptot-npint))=&!1d0!
	    diag(2*nptot+mapg(ij,el)-(nptot-npint)) + M_pressure(intij,intij,el)
      ENDDO
      ENDDO
    ENDDO
    DO i=1,2*nptot+npint
      IF (abs(diag(i)).lt.1d-16) THEN
	diag(i)=1d0
	print*, 'Warning: Diagonal preconditioner had zero value! Set to 1 instead.'
      ENDIF
    ENDDO
    
  END SUBROUTINE calc_reciprical_of_diagA
  
!   SUBROUTINE element_matrix_vector_mult(trans,mat_in,vec_in,vec_out,rows,cols
  
END MODULE solver
