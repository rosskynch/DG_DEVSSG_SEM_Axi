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
      print*,normr
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
	diag(ij,el) = pressuremassmatrix(ij-2*NM1SQ,ij-2*NM1SQ,el)
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