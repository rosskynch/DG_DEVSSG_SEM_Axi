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

MODULE devss_module
  
  USE shared_data
  USE functions_module
  IMPLICIT NONE

  
  
  CONTAINS
  
  SUBROUTINE applyDEVSS
    IMPLICIT NONE
    INTEGER :: globi,el,ij

! Save old solution
!       devss_contrib_xNm2=devss_contrib_xNm1
!       devss_contrib_yNm2=devss_contrib_yNm1
      devss_contrib_xNm1=devss_contrib_x
      devss_contrib_yNm1=devss_contrib_y

! USE DEVSS:
!       CALL calc_devss_solution
!       devss_soln_xx=0d0
!       devss_soln_xy=0d0
!       devss_soln_yy=0d0
!       devss_soln_zz=0d0
!       DO el=1,numelm
!   DO ij=0,NP1SQM1
!     globi=mapg(ij,el)
!     devss_soln_xx(ij,el)=devss_soln_xx(ij,el) + localGradUxx(ij,el)/dfloat(mult(globi))
!     devss_soln_xy(ij,el)=devss_soln_xy(ij,el) + 5d-1*(localGradUxy(ij,el)+localGradUyx(ij,el))/dfloat(mult(globi))
!     devss_soln_yy(ij,el)=devss_soln_yy(ij,el) + localGradUyy(ij,el)/dfloat(mult(globi))
!     devss_soln_zz(ij,el)=devss_soln_zz(ij,el) + localGradUzz(ij,el)/dfloat(mult(globi))
!   ENDDO
!       ENDDO
    
  
!       devss_soln_xx=localGradUxx
!       devss_soln_xy=5d-1*(localGradUxy+localGradUyx)
!       devss_soln_yy=localGradUyy
!       devss_soln_zz=localGradUzz
!       CALL calc_devss_contribution
      
! USE DEVSS-G
    CALL calc_devssg_solution
!     devss_soln_xx=localGradUxx
!     devss_soln_xy=localGradUxy
!     devss_soln_yx=localGradUyx
!     devss_soln_yy=localGradUyy
!     devss_soln_zz=localGradUzz
    CALL calc_devssg_contribution


! DEVSS version: (using (beta_s-beta)Div(1.2( G + G^T )))
!       f_x = f_x - 2d0*(param_beta_s-param_beta)*(time_beta_0*devss_contrib_x + time_beta_1*devss_contrib_xNm1)! + time_beta_2*devss_contrib_xNm2)
!       f_y = f_y - 2d0*(param_beta_s-param_beta)*(time_beta_0*devss_contrib_y + time_beta_1*devss_contrib_yNm1)! + time_beta_2*devss_contrib_yNm2)
!DEVSS-G version: (using (beta_s-beta)Div(G))
      f_x = f_x - (param_beta_s-param_beta)*(time_beta_0*devss_contrib_x + time_beta_1*devss_contrib_xNm1)! + time_beta_2*devss_contrib_xNm2)
      f_y = f_y - (param_beta_s-param_beta)*(time_beta_0*devss_contrib_y + time_beta_1*devss_contrib_yNm1)! + time_beta_2*devss_contrib_yNm2)

  END SUBROUTINE applyDEVSS
 
 
  SUBROUTINE initialise_devss
    IMPLICIT NONE
    INTEGER :: el,info
    DOUBLE PRECISION :: temp_matrix(NM1SQ,NM1SQ)
    
    devss_matrix = M_pressure
    devss_soln_xx=0d0
    devss_soln_xy=0d0
    devss_soln_yx=0d0
    devss_soln_yy=0d0
    devss_soln_zz=0d0
    devss_contrib_x=0d0
    devss_contrib_y=0d0
    devss_contrib_xNm1=0d0
    devss_contrib_yNm1=0d0
! Factorise the matrix
    DO el=1,numelm
      temp_matrix=devss_matrix(1:NM1SQ,1:NM1SQ,el)
      CALL dpotrf( 'U', NM1SQ, temp_matrix, NM1SQ, info )
      IF (info.ne.0) THEN
        write(*,*) 'Error in initialise_devss:',el,info
        STOP
      ENDIF
      devss_matrix(1:NM1SQ,1:NM1SQ,el)=temp_matrix
    ENDDO
    
  END SUBROUTINE initialise_devss
  
!   SUBROUTINE calc_adaptive_viscosity
!     IMPLICIT NONE
!     INTEGER :: el,ij
!     DOUBLE PRECISION :: a,b
! 
!     a=0d0
!     b=1d-3
!     DO el=1,numelm
!       DO ij=0,NP1SQM1
!   IF (abs(localTxx(ij,el)).gt.a) a=abs(localTxx(ij,el))
!   IF (abs(localTxy(ij,el)).gt.a) a=abs(localTxy(ij,el))
!   IF (abs(localTyy(ij,el)).gt.a) a=abs(localTyy(ij,el))
!   IF (abs(localTzz(ij,el)).gt.a) a=abs(localTzz(ij,el))
!   IF (abs(2d0*devss_soln_xx(ij,el)).gt.b) b=abs(2d0*devss_soln_xx(ij,el))
!   IF (abs(devss_soln_xy(ij,el)+devss_soln_yx(ij,el)).gt.b) b=abs(devss_soln_xy(ij,el)+devss_soln_yx(ij,el))
!   IF (abs(2d0*devss_soln_yy(ij,el)).gt.b) b=abs(2d0*devss_soln_yy(ij,el))
!   IF (abs(2d0*devss_soln_zz(ij,el)).gt.b) b=abs(2d0*devss_soln_zz(ij,el))
! 
! !       param_beta_a(el) =sqrt( 1d0 + 0.5d0*param_delta_a*(localTxx(ij,el)**2 + 2d0*localTxy(ij,el)**2 + localTyy(ij,el)**2 + localTzz(ij,el)**2) )&
! !                 / &
! !       sqrt( 1d0 + 0.5d0*param_delta_a*( 4d0*(devss_soln_xx(ij,el)**2 + devss_soln_yy(ij,el)**2 + devss_soln_zz(ij,el)**2) + &
! !               2d0*(devss_soln_xy(ij,el) + devss_soln_yx(ij,el))**2 ) ) 
! 
! 
!       ENDDO
! !         print*, param_beta_a(el)  
!     ENDDO
!     param_beta_a = sqrt(1d0 + param_delta_a*a/b)
! 
!   END SUBROUTINE calc_adaptive_viscosity
  
! DEVSS-G VERSION: i.e. we calculate G = gradU

  SUBROUTINE calc_devssg_rhs
    IMPLICIT NONE
    INTEGER :: el,ij,kl,i,j,k,l,jj,ll
    DOUBLE PRECISION :: temp11,temp12,temp21,temp22,temp33,&
      temp1,temp2,temp3,temp4,temp5,tempmult

    IF (coordflag.eq.0) THEN
! CARTESIAN CASE !    
    DO el=1,numelm
      DO j=1,NM1
        jj=(j-1)*NM1
        DO i=1,NM1
          ij=i+jj
          temp11=0d0
          temp12=0d0
          temp21=0d0
          temp22=0d0  
          DO l=0,N
            ll=l*NP1
            tempmult=evalh(j,l)*w(l)
            temp1=0d0
            temp2=0d0
            temp3=0d0
            temp4=0d0
            DO k=0,N
              kl=k+ll
              temp1 = temp1 + localGradUxx(kl,el)*evalh(i,k)*jac(k,l,el)*w(k)
              temp2 = temp2 + localGradUxy(kl,el)*evalh(i,k)*jac(k,l,el)*w(k)
              temp3 = temp3 + localGradUyx(kl,el)*evalh(i,k)*jac(k,l,el)*w(k)
              temp4 = temp4 + localGradUyy(kl,el)*evalh(i,k)*jac(k,l,el)*w(k)
            ENDDO
            temp11 = temp11 + tempmult*temp1
            temp12 = temp12 + tempmult*temp2
            temp21 = temp21 + tempmult*temp3
            temp22 = temp22 + tempmult*temp4
          ENDDO

          devss_rhs_xx(ij,el) = temp11
          devss_rhs_xy(ij,el) = temp12
          devss_rhs_yx(ij,el) = temp21
          devss_rhs_yy(ij,el) = temp22
        ENDDO
      ENDDO
    ENDDO
    ELSEIF (coordflag.eq.1) THEN
      DO el=1,numelm
        DO j=1,NM1
          jj=(j-1)*NM1
          DO i=1,NM1
            ij=i+jj
            temp11=0d0
            temp12=0d0
            temp21=0d0
            temp22=0d0
            temp33=0d0
            DO l=0,N
              ll=l*NP1
              tempmult=evalh(j,l)*w(l)
              temp1=0d0
              temp2=0d0
              temp3=0d0
              temp4=0d0
              temp5=0d0
              DO k=0,N
                kl=k+ll
                temp1 = temp1 + localGradUxx(kl,el)*evalh(i,k)*jac(k,l,el)*w(k)*nodeCoord(mapg(kl,el),2)
                temp2 = temp2 + localGradUxy(kl,el)*evalh(i,k)*jac(k,l,el)*w(k)*nodeCoord(mapg(kl,el),2)
                temp3 = temp3 + localGradUyx(kl,el)*evalh(i,k)*jac(k,l,el)*w(k)*nodeCoord(mapg(kl,el),2)
                temp4 = temp4 + localGradUyy(kl,el)*evalh(i,k)*jac(k,l,el)*w(k)*nodeCoord(mapg(kl,el),2)
                temp5 = temp5 + V_y(mapg(kl,el))*evalh(i,k)*jac(k,l,el)*w(k)
              ENDDO
              temp11 = temp11 + tempmult*temp1
              temp12 = temp12 + tempmult*temp2
              temp21 = temp21 + tempmult*temp3
              temp22 = temp22 + tempmult*temp4
              temp33 = temp33 + tempmult*temp5
            ENDDO
            devss_rhs_xx(ij,el) = temp11
            devss_rhs_xy(ij,el) = temp12
            devss_rhs_yx(ij,el) = temp21
            devss_rhs_yy(ij,el) = temp22
            devss_rhs_zz(ij,el) = temp33
          ENDDO
        ENDDO
      ENDDO

    ELSE
      print*, 'Error: No co-ordinate system specified...'
      print*, 'Stopping'
      STOP
    ENDIF

  END SUBROUTINE calc_devssg_rhs
  
  SUBROUTINE calc_devssg_solution
    IMPLICIT NONE
    INTEGER :: el,info,edge,i,j,jj,k,l,p,q,pq,intij,ij,nrhs
    DOUBLE PRECISION :: temp_matrix(NM1SQ,NM1SQ),temp,&
      tempxx,tempxy,tempyx,tempyy,tempzz
      
    DOUBLE PRECISION, ALLOCATABLE :: temp_rhs(:,:)
    
    IF (coordflag.eq.1) THEN
      nrhs=5
      ALLOCATE(temp_rhs(NM1SQ,5))
    ELSE
      nrhs=4
      ALLOCATE(temp_rhs(NM1SQ,4))
    ENDIF
    
    IF (movingmeshflag.eq.1) THEN
      CALL initialise_devss
    ENDIF
    
    CALL calc_devssg_rhs

    DO el=1,numelm
      temp_rhs(1:NM1SQ,1) = devss_rhs_xx(1:NM1SQ,el)
      temp_rhs(1:NM1SQ,2) = devss_rhs_xy(1:NM1SQ,el)
      temp_rhs(1:NM1SQ,3) = devss_rhs_yx(1:NM1SQ,el)
      temp_rhs(1:NM1SQ,4) = devss_rhs_yy(1:NM1SQ,el)
      IF (coordflag.eq.1) THEN 
        temp_rhs(1:NM1SQ,5) = devss_rhs_zz(1:NM1SQ,el)
      ENDIF
      temp_matrix = devss_matrix(1:NM1SQ,1:NM1SQ,el)
      CALL dpotrs( 'U', NM1SQ, nrhs, temp_matrix, NM1SQ, temp_rhs, NM1SQ, info )
!       CALL dgetrs('N', NM1SQ, nrhs, temp_matrix, NM1SQ, devss_ipiv, temp_rhs, NM1SQ, info )
      IF (info.ne.0) THEN
        write(*,*) 'Error in calc_devssg_solution:',el,info
        STOP
      ENDIF
      devss_rhs_xx(1:NM1SQ,el) = temp_rhs(1:NM1SQ,1)
      devss_rhs_xy(1:NM1SQ,el) = temp_rhs(1:NM1SQ,2)
      devss_rhs_yx(1:NM1SQ,el) = temp_rhs(1:NM1SQ,3)
      devss_rhs_yy(1:NM1SQ,el) = temp_rhs(1:NM1SQ,4)
      IF (coordflag.eq.1) THEN 
        devss_rhs_zz(1:NM1SQ,el) = temp_rhs(1:NM1SQ,5)
      ENDIF
    ENDDO
    DEALLOCATE(temp_rhs)

! Copy to solution storage:
    DO el=1,numelm
! Extrapolate boundary values...
      DO edge=1,4
        DO k=0,N
          pq = local_edge_node(k,edge)
          p=local_ij_to_i(pq)
          q=local_ij_to_j(pq)
          tempxx=0d0
          tempxy=0d0
          tempyx=0d0
          tempyy=0d0
          DO j=1,NM1
!             jj=j*NP1
            jj=(j-1)*NM1
            DO i=1,NM1
              ij=i+jj
              tempxx = tempxx + devss_rhs_xx(ij,el)*evalh(i,p)*evalh(j,q)
              tempxy = tempxy + devss_rhs_xy(ij,el)*evalh(i,p)*evalh(j,q)
              tempyx = tempyx + devss_rhs_yx(ij,el)*evalh(i,p)*evalh(j,q)
              tempyy = tempyy + devss_rhs_yy(ij,el)*evalh(i,p)*evalh(j,q)

            ENDDO
          ENDDO

          devss_soln_xx(pq,el) = tempxx
          devss_soln_xy(pq,el) = tempxy
          devss_soln_yx(pq,el) = tempyx
          devss_soln_yy(pq,el) = tempyy
          IF (coordflag.eq.1) THEN
            tempzz=0d0
            DO i=1,NM1
              DO j=1,NM1
                ij=i+(j-1)*NM1
                tempzz = tempzz + devss_rhs_zz(ij,el)*evalh(i,p)*evalh(j,q)
              ENDDO
            ENDDO
            devss_soln_zz(pq,el) = tempzz
          ENDIF
        ENDDO
      ENDDO
! Fill in interior values
      DO intij=1,NM1SQ
        ij=interior_to_local_node(intij)
        devss_soln_xx(ij,el)=devss_rhs_xx(intij,el)
        devss_soln_xy(ij,el)=devss_rhs_xy(intij,el)
        devss_soln_yx(ij,el)=devss_rhs_yx(intij,el)
        devss_soln_yy(ij,el)=devss_rhs_yy(intij,el)
      ENDDO
      IF (coordflag.eq.1) THEN
        DO intij=1,NM1SQ
          ij=interior_to_local_node(intij)
          devss_soln_zz(ij,el)=devss_rhs_zz(intij,el)
        ENDDO
      ENDIF
    ENDDO  
  END SUBROUTINE calc_devssg_solution
  
  SUBROUTINE calc_devssg_contribution
    IMPLICIT NONE
    INTEGER :: el,i,j,k,l,kl,il,ki,ij,edge,ll
    DOUBLE PRECISION :: temp1,temp2,temp3,temp4,normx,normy
    DOUBLE PRECISION, DIMENSION(nptot) :: tempglob
    DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1) :: templocal_OUTx,templocal_OUTy
     
    
! Integrate the solutions in each element
    devss_contrib_x=0d0
    devss_contrib_y=0d0

    IF (coordflag.eq.0) THEN
! CARTESIAN CASE !
      DO el=1,numelm
        templocal_OUTx=0d0
        templocal_OUTy=0d0
        DO l=0,N
          ll=l*NP1
          DO k=0,N
            kl=k+ll
            temp1=0d0
            temp2=0d0
            temp3=0d0
            temp4=0d0
            DO i=0,N
              ki=k+i*NP1
              il=i+ll

! Integrating Div(G)
! z-component:
              temp1 = temp1 + w(i)*d(i,k)*( devss_soln_xx(il,el)*dyde(i,l,el) - devss_soln_yx(il,el)*dxde(i,l,el) ) 
              temp2 = temp2 + w(i)*d(i,l)*( devss_soln_yx(ki,el)*dxdp(k,i,el) - devss_soln_xx(ki,el)*dydp(k,i,el) )
! r-component:
              temp3 = temp3 + w(i)*d(i,k)*( devss_soln_xy(il,el)*dyde(i,l,el) - devss_soln_yy(il,el)*dxde(i,l,el) )
              temp4 = temp4 + w(i)*d(i,l)*( devss_soln_yy(ki,el)*dxdp(k,i,el) - devss_soln_xy(ki,el)*dydp(k,i,el) )

! Integrating Div( D ) | D = 1/2(G + G^T)
! x-comp sums:
! used i as dummy variable, workings may differ.
!       temp1 = temp1 + w(i)*d(i,k)*( 2d0*devss_soln_xx(il,el)*dyde(i,l,el) - (devss_soln_xy(il,el) + devss_soln_yx(il,el))*dxde(i,l,el) ) 
!       temp2 = temp2 + w(i)*d(i,l)*( (devss_soln_xy(ki,el) + devss_soln_yx(ki,el))*dxdp(k,i,el) - 2d0*devss_soln_xx(ki,el)*dydp(k,i,el) )
!       temp1 = temp1 + w(i)*d(i,k)*( devss_soln_xx(il,el)*dyde(i,l,el) - 5d-1*(devss_soln_xy(il,el) + devss_soln_yx(il,el))*dxde(i,l,el) ) 
!       temp2 = temp2 + w(i)*d(i,l)*( 5d-1*(devss_soln_xy(ki,el) + devss_soln_yx(ki,el))*dxdp(k,i,el) - devss_soln_xx(ki,el)*dydp(k,i,el) )
! y-comp sums:
!       temp3 = temp3 + w(i)*d(i,k)*( (devss_soln_xy(il,el) + devss_soln_yx(il,el))*dyde(i,l,el) - 2d0*devss_soln_yy(il,el)*dxde(i,l,el) )
!       temp4 = temp4 + w(i)*d(i,l)*( 2d0*devss_soln_yy(ki,el)*dxdp(k,i,el) - (devss_soln_xy(ki,el) + devss_soln_yx(ki,el))*dydp(k,i,el) )
!       temp3 = temp3 + w(i)*d(i,k)*( 5d-1*(devss_soln_xy(il,el) + devss_soln_yx(il,el))*dyde(i,l,el) - devss_soln_yy(il,el)*dxde(i,l,el) )
!       temp4 = temp4 + w(i)*d(i,l)*( devss_soln_yy(ki,el)*dxdp(k,i,el) - 5d-1*(devss_soln_xy(ki,el) + devss_soln_yx(ki,el))*dydp(k,i,el) )
      
            ENDDO
            IF (.not.bdflag(1,mapg(kl,el))) THEN 
              templocal_OUTx(kl) = (w(l)*temp1 + w(k)*temp2)
            ENDIF
            IF (.not.bdflag(2,mapg(kl,el))) THEN
              templocal_OUTy(kl) = (w(l)*temp3 + w(k)*temp4)
            ENDIF
          ENDDO
        ENDDO
        templocal_OUTx=templocal_OUTx
        templocal_OUTy=templocal_OUTy
        tempglob=0d0
        CALL vecglobalprolongation(templocal_OUTx,el,tempglob)
        devss_contrib_x = devss_contrib_x - tempglob
        tempglob=0d0
        CALL vecglobalprolongation(templocal_OUTy,el,tempglob)
        devss_contrib_y = devss_contrib_y - tempglob
      ENDDO
    ELSEIF (coordflag.eq.1) THEN 
! CYLINDERICAL POLAR (AXISYMMETRIC) CASE !
!
! NOTE:
! The x co-ordinate in the cartesian case is equivalent to the z component
! in the cylindrical polar coordinate case.
!      
      DO el=1,numelm
      templocal_OUTx=0d0
      templocal_OUTy=0d0
        DO l=0,N
          ll=l*NP1
          DO k=0,N
            kl=k+ll
! Contribution from this element for the x-component:
            temp1=0d0
            temp2=0d0
            temp3=0d0
            temp4=0d0
            DO i=0,N
              ki=k+i*NP1
              il=i+ll
      
!Integrating Div(G)
! z-component:
              temp1 = temp1 + nodeCoord(mapg(il,el),2)*w(i)*d(i,k)*( devss_soln_xx(il,el)*dyde(i,l,el) - devss_soln_yx(il,el)*dxde(i,l,el) ) 
              temp2 = temp2 + nodeCoord(mapg(ki,el),2)*w(i)*d(i,l)*( devss_soln_yx(ki,el)*dxdp(k,i,el) - devss_soln_xx(ki,el)*dydp(k,i,el) )
! r-component:
              temp3 = temp3 + nodeCoord(mapg(il,el),2)*w(i)*d(i,k)*( devss_soln_xy(il,el)*dyde(i,l,el) - devss_soln_yy(il,el)*dxde(i,l,el) )
              temp4 = temp4 + nodeCoord(mapg(ki,el),2)*w(i)*d(i,l)*( devss_soln_yy(ki,el)*dxdp(k,i,el) - devss_soln_xy(ki,el)*dydp(k,i,el) )

! Integrating Div( D ) | D = 1/2(G + G^T)
! z-comp sums:
! used i as dummy variable, workings may differ.
!       temp1 = temp1 + nodeCoord(mapg(il,el),2)*w(i)*d(i,k)*( 2d0*devss_soln_xx(il,el)*dyde(i,l,el) - (devss_soln_xy(il,el) + devss_soln_yx(il,el))*dxde(i,l,el) ) 
!       temp2 = temp2 + nodeCoord(mapg(ki,el),2)*w(i)*d(i,l)*( (devss_soln_xy(ki,el) + devss_soln_yx(ki,el))*dxdp(k,i,el) - 2d0*devss_soln_xx(ki,el)*dydp(k,i,el) )
!       temp1 = temp1 + nodeCoord(mapg(il,el),2)*w(i)*d(i,k)*( devss_soln_xx(il,el)*dyde(i,l,el) - 5d-1*(devss_soln_xy(il,el) + devss_soln_yx(il,el))*dxde(i,l,el) ) 
!       temp2 = temp2 + nodeCoord(mapg(ki,el),2)*w(i)*d(i,l)*( 5d-1*(devss_soln_xy(ki,el) + devss_soln_yx(ki,el))*dxdp(k,i,el) - devss_soln_xx(ki,el)*dydp(k,i,el) )
! r-comp sums:
!       temp3 = temp3 + nodeCoord(mapg(il,el),2)*w(i)*d(i,k)*( (devss_soln_xy(il,el) + devss_soln_yx(il,el))*dyde(i,l,el) - 2d0*devss_soln_yy(il,el)*dxde(i,l,el) )
!       temp4 = temp4 + nodeCoord(mapg(ki,el),2)*w(i)*d(i,l)*( 2d0*devss_soln_yy(ki,el)*dxdp(k,i,el) - (devss_soln_xy(ki,el) + devss_soln_yx(ki,el))*dydp(k,i,el) )
!       temp3 = temp3 + nodeCoord(mapg(il,el),2)*w(i)*d(i,k)*( 5d-1*(devss_soln_xy(il,el) + devss_soln_yx(il,el))*dyde(i,l,el) - devss_soln_yy(il,el)*dxde(i,l,el) )
!       temp4 = temp4 + nodeCoord(mapg(ki,el),2)*w(i)*d(i,l)*( devss_soln_yy(ki,el)*dxdp(k,i,el) - 5d-1*(devss_soln_xy(ki,el) + devss_soln_yx(ki,el))*dydp(k,i,el) )
            ENDDO

            IF (.not.bdflag(1,mapg(kl,el))) THEN
              templocal_OUTx(kl) = w(l)*temp1 + w(k)*temp2
            ENDIF
! Add the extra rr and theta,theta contribution to divergence in y(r)-comp case:
            IF (.not.bdflag(2,mapg(kl,el))) THEN
!               templocal_OUTy(kl) = 2d0*devss_soln_zz(kl,el)*jac(k,l,el)*w(k)*w(l) + w(l)*temp3 + w(k)*temp4
              templocal_OUTy(kl) = devss_soln_zz(kl,el)*jac(k,l,el)*w(k)*w(l) + w(l)*temp3 + w(k)*temp4
            ENDIF
          ENDDO
        ENDDO
        templocal_OUTx=templocal_OUTx
        templocal_OUTy=templocal_OUTy
        tempglob=0d0
        CALL vecglobalprolongation(templocal_OUTx,el,tempglob)
        devss_contrib_x = devss_contrib_x - tempglob
        tempglob=0d0
        CALL vecglobalprolongation(templocal_OUTy,el,tempglob)
        devss_contrib_y = devss_contrib_y - tempglob
      ENDDO
    ELSE
      print*, 'Error: No co-ordinate system specified...'
      print*, 'Stopping'
      STOP
    ENDIF

! Boundary Contribution from non-dirichlet edge nodes (i.e. the axis of symmetry)
! For the axis of symmetry, only the x/z-component has a neumann BC. The test function is zero in other cases
!     DO el=1,numelm
!       edge=axisymm_edge(el)
! !       temp_beta=param_beta_a(el) - param_beta
!       IF (edge.ne.0) THEN
!   DO i=0,N
!     ij=local_edge_node(i,edge)
!     k=mapg(ij,el)
!     normx=norm_to_edge_node(1,i,edge,el)
!     normy=norm_to_edge_node(2,i,edge,el)
!     IF (.not.bdflag(1,k)) THEN
!       devss_contrib_x(k) = devss_contrib_x(k) + &
!                 ( devss_soln_xx(ij,el)*normx + 5d-1*(devss_soln_xy(ij,el) + devss_soln_yx(ij,el))*normy )*jac_on_edge(edge,el)*w(i)
!     ENDIF
!     IF (.not.bdflag(2,k)) THEN
!       devss_contrib_y(k) = devss_contrib_y(k) + &
!                 ( 5d-1*(devss_soln_xy(ij,el) + devss_soln_yx(ij,el))*normx + devss_soln_yy(ij,el)*normy )*jac_on_edge(edge,el)*w(i)
!     ENDIF
! ! WRONG ???
! !     IF (.not.bdflag(1,k)) THEN
! !       devss_contrib_x(k) = devss_contrib_x(k) + &
! !                 ( 2d0*devss_soln_xx(ij,el)*normx + (devss_soln_xy(ij,el) + devss_soln_yx(ij,el))*normy )*jac_on_edge(edge,el)*w(i)
! !     ENDIF
! !     IF (.not.bdflag(2,k)) THEN
! !       devss_contrib_y(k) = devss_contrib_y(k) + &
! !                 ( (devss_soln_xy(ij,el) + devss_soln_yx(ij,el))*normx + 2d0*devss_soln_yy(ij,el)*normy )*jac_on_edge(edge,el)*w(i)
! !     ENDIF
!   ENDDO
!       ENDIF
!     ENDDO

  END SUBROUTINE calc_devssg_contribution


! NORMAL DEVSS (i.e. We calculate a projection of the rate of strain tensors: gamma = gradU + gradU^T
  SUBROUTINE calc_devss_rhs
    IMPLICIT NONE
    INTEGER :: el,ij,kl,i,j,k,l,jj,ll
    DOUBLE PRECISION :: tempxx,tempxy,tempyy,tempzz,&
      temp1,temp2,temp3,temp4,tempmult
!       tempBx(NM1SQ,0:NP1SQ-1),tempBy(NM1SQ,0:NP1SQ-1)
    IF (coordflag.eq.0) THEN
! CARTESIAN CASE !    
      DO el=1,numelm
         DO j=1,NM1
          jj=(j-1)*NM1
          DO i=1,NM1
            ij=i+jj
            tempxx=0d0
            tempxy=0d0
            tempyy=0d0
            DO l=0,N
              ll=l*NP1
              tempmult=evalh(j,l)*w(l)
              temp1=0d0
              temp2=0d0
              temp3=0d0
              DO k=0,N
                kl=k+ll
                temp1 = temp1 + localGradUxx(kl,el)*evalh(i,k)*jac(k,l,el)*w(k)
                temp2 = temp2 + 5d-1*(localGradUxy(kl,el) + localGradUyx(kl,el))*evalh(i,k)*jac(k,l,el)*w(k)
                temp3 = temp3 + localGradUyy(kl,el)*evalh(i,k)*jac(k,l,el)*w(k)
!         temp1 = temp1 + 2d0*localGradUxx(kl,el)*evalh(i,k)*jac(k,l,el)*w(k)
!         temp2 = temp2 + (localGradUxy(kl,el) + localGradUyx(kl,el))*evalh(i,k)*jac(k,l,el)*w(k)
!         temp3 = temp3 + 2d0*localGradUyy(kl,el)*evalh(i,k)*jac(k,l,el)*w(k)
      ENDDO
      tempxx = tempxx + tempmult*temp1
      tempxy = tempxy + tempmult*temp2
      tempyy = tempyy + tempmult*temp3

    ENDDO

    devss_rhs_xx(ij,el) = tempxx
    devss_rhs_xy(ij,el) = tempxy
    devss_rhs_yy(ij,el) = tempyy
  ENDDO
      ENDDO
    ENDDO
    ELSEIF (coordflag.eq.1) THEN
      DO el=1,numelm
  DO j=1,NM1
    jj=(j-1)*NM1
    DO i=1,NM1
      ij=i+jj
      tempxx=0d0
      tempxy=0d0
      tempyy=0d0
      tempzz=0d0
      DO l=0,N
        ll=l*NP1
        tempmult=evalh(j,l)*w(l)
        temp1=0d0
        temp2=0d0
        temp3=0d0
        temp4=0d0
        DO k=0,N
    kl=k+ll
    temp1 = temp1 + localGradUxx(kl,el)*evalh(i,k)*jac(k,l,el)*w(k)*nodeCoord(mapg(kl,el),2)
    temp2 = temp2 + 5d-1*(localGradUyx(kl,el) + localGradUxy(kl,el))*evalh(i,k)*jac(k,l,el)*w(k)*nodeCoord(mapg(kl,el),2)
    temp3 = temp3 + localGradUyy(kl,el)*evalh(i,k)*jac(k,l,el)*w(k)*nodeCoord(mapg(kl,el),2)
    temp4 = temp4 + V_y(mapg(kl,el))*evalh(i,k)*jac(k,l,el)*w(k)
!     temp1 = temp1 + 2d0*localGradUxx(kl,el)*evalh(i,k)*jac(k,l,el)*w(k)*nodeCoord(mapg(kl,el),2)
!     temp2 = temp2 + (localGradUyx(kl,el) + localGradUxy(kl,el))*evalh(i,k)*jac(k,l,el)*w(k)*nodeCoord(mapg(kl,el),2)
!     temp3 = temp3 + 2d0*localGradUyy(kl,el)*evalh(i,k)*jac(k,l,el)*w(k)*nodeCoord(mapg(kl,el),2)
!     temp4 = temp4 + 2d0*V_y(mapg(kl,el))*evalh(i,k)*jac(k,l,el)*w(k)
        ENDDO
        tempxx = tempxx + tempmult*temp1
        tempxy = tempxy + tempmult*temp2
        tempyy = tempyy + tempmult*temp3
        tempzz = tempzz + tempmult*temp4
      ENDDO
      devss_rhs_xx(ij,el) = tempxx
      devss_rhs_xy(ij,el) = tempxy
      devss_rhs_yy(ij,el) = tempyy
      devss_rhs_zz(ij,el) = tempzz
    ENDDO
  ENDDO
      ENDDO

    ELSE
      print*, 'Error: No co-ordinate system specified...'
      print*, 'Stopping'
      STOP
    ENDIF

  END SUBROUTINE calc_devss_rhs
  
  SUBROUTINE calc_devss_solution
    IMPLICIT NONE
    INTEGER :: el,info,edge,i,j,jj,k,l,p,q,pq,intij,ij,nrhs
    DOUBLE PRECISION :: temp_matrix(NM1SQ,NM1SQ),temp,&
      tempxx,tempxy,tempyx,tempyy,tempzz
      
    DOUBLE PRECISION, ALLOCATABLE :: temp_rhs(:,:)
    
    IF (coordflag.eq.1) THEN
      nrhs=4
      ALLOCATE(temp_rhs(NM1SQ,5))
    ELSE
      nrhs=3
      ALLOCATE(temp_rhs(NM1SQ,4))
    ENDIF
    
    IF (movingmeshflag.eq.1) THEN
      CALL initialise_devss
    ENDIF
    
    CALL calc_devss_rhs

    DO el=1,numelm
      temp_rhs(1:NM1SQ,1) = devss_rhs_xx(1:NM1SQ,el)
      temp_rhs(1:NM1SQ,2) = devss_rhs_xy(1:NM1SQ,el)
      temp_rhs(1:NM1SQ,3) = devss_rhs_yy(1:NM1SQ,el)
      IF (coordflag.eq.1) THEN 
  temp_rhs(1:NM1SQ,4) = devss_rhs_zz(1:NM1SQ,el)
      ENDIF
      temp_matrix = devss_matrix(1:NM1SQ,1:NM1SQ,el)
      CALL dpotrs( 'U', NM1SQ, nrhs, temp_matrix, NM1SQ, temp_rhs, NM1SQ, info )
!       CALL dgetrs('N', NM1SQ, nrhs, temp_matrix, NM1SQ, devss_ipiv, temp_rhs, NM1SQ, info )
      IF (info.ne.0) THEN
  write(*,*) 'Error in calc_devss_solution:',el,info
  STOP
      ENDIF
      devss_rhs_xx(1:NM1SQ,el) = temp_rhs(1:NM1SQ,1)
      devss_rhs_xy(1:NM1SQ,el) = temp_rhs(1:NM1SQ,2)
      devss_rhs_yy(1:NM1SQ,el) = temp_rhs(1:NM1SQ,3)
      IF (coordflag.eq.1) THEN 
  devss_rhs_zz(1:NM1SQ,el) = temp_rhs(1:NM1SQ,4)
      ENDIF
    ENDDO
    DEALLOCATE(temp_rhs)

! Copy to solution storage:
    DO el=1,numelm
! Extrapolate boundary values...
      DO edge=1,4
  DO k=0,N
    pq = local_edge_node(k,edge)
    p=local_ij_to_i(pq)
    q=local_ij_to_j(pq)
    tempxx=0d0
    tempxy=0d0
    tempyy=0d0
    DO j=1,NM1
      jj=(j-1)*NM1
      DO i=1,NM1
        ij=i+jj
        tempxx = tempxx + devss_rhs_xx(ij,el)*evalh(i,p)*evalh(j,q)
        tempxy = tempxy + devss_rhs_xy(ij,el)*evalh(i,p)*evalh(j,q)
        tempyy = tempyy + devss_rhs_yy(ij,el)*evalh(i,p)*evalh(j,q)
        
      ENDDO
    ENDDO

    devss_soln_xx(pq,el) = tempxx
    devss_soln_xy(pq,el) = tempxy
    devss_soln_yy(pq,el) = tempyy
    IF (coordflag.eq.1) THEN
      tempzz=0d0
      DO j=1,NM1
      jj = (j-1)*NM1
        DO i=1,NM1
    ij=i+jj
    tempzz = tempzz + devss_rhs_zz(ij,el)*evalh(i,p)*evalh(j,q)
        ENDDO
      ENDDO
      devss_soln_zz(pq,el) = tempzz
    ENDIF
  ENDDO
      ENDDO
! Fill in interior values
      DO intij=1,NM1SQ
  ij=interior_to_local_node(intij)
  devss_soln_xx(ij,el)=devss_rhs_xx(intij,el)
  devss_soln_xy(ij,el)=devss_rhs_xy(intij,el)
  devss_soln_yy(ij,el)=devss_rhs_yy(intij,el)
      ENDDO
      IF (coordflag.eq.1) THEN
  DO intij=1,NM1SQ
    ij=interior_to_local_node(intij)
    devss_soln_zz(ij,el)=devss_rhs_zz(intij,el)
  ENDDO
      ENDIF
      
    ENDDO
    
   
  END SUBROUTINE calc_devss_solution
  
  SUBROUTINE calc_devss_contribution
    IMPLICIT NONE
    INTEGER :: el,i,j,k,l,kl,il,ki,ij,edge,ll
    DOUBLE PRECISION :: temp1,temp2,temp3,temp4,normx,normy,temp_beta
    DOUBLE PRECISION, DIMENSION(nptot) :: tempglob
    DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1) :: templocal_OUTx,templocal_OUTy
    
! Integrate the solutions in each element
    devss_contrib_x=0d0
    devss_contrib_y=0d0

    IF (coordflag.eq.0) THEN
! CARTESIAN CASE !      
    DO el=1,numelm
  templocal_OUTx=0d0
  templocal_OUTy=0d0
!   temp_beta=param_beta_a(el) - param_beta
  DO l=0,N
    ll=l*NP1
  DO k=0,N
    kl=k+ll
    temp1=0d0
    temp2=0d0
    temp3=0d0
    temp4=0d0
    DO i=0,N
      ki=k+i*NP1
      il=i+ll
      
      
! integrating Div(G) 
! x-comp sums:
! used i as dummy variable, workings may differ.
      temp1 = temp1 + w(i)*d(i,k)*( devss_soln_xx(il,el)*dyde(i,l,el) - devss_soln_xy(il,el)*dxde(i,l,el) ) 
      temp2 = temp2 + w(i)*d(i,l)*( devss_soln_xy(ki,el)*dxdp(k,i,el) - devss_soln_xx(ki,el)*dydp(k,i,el) )
! y-comp sums:
      temp3 = temp3 + w(i)*d(i,k)*( devss_soln_xy(il,el)*dyde(i,l,el) - devss_soln_yy(il,el)*dxde(i,l,el) )
      temp4 = temp4 + w(i)*d(i,l)*( devss_soln_yy(ki,el)*dxdp(k,i,el) - devss_soln_xy(ki,el)*dydp(k,i,el) )
      
    ENDDO
    IF (.not.bdflag(1,mapg(kl,el))) THEN
      templocal_OUTx(kl) = (w(l)*temp1 + w(k)*temp2)
    ENDIF
    IF (.not.bdflag(2,mapg(kl,el))) THEN
      templocal_OUTy(kl) = (w(l)*temp3 + w(k)*temp4)
    ENDIF

  ENDDO
  ENDDO
  tempglob=0d0
  CALL vecglobalprolongation(templocal_OUTx,el,tempglob)
  devss_contrib_x = devss_contrib_x - tempglob
  tempglob=0d0
  CALL vecglobalprolongation(templocal_OUTy,el,tempglob)
  devss_contrib_y = devss_contrib_y - tempglob
      ENDDO
    ELSEIF (coordflag.eq.1) THEN 
! CYLINDERICAL POLAR (AXISYMMETRIC) CASE !
!
! NOTE:
! The x co-ordinate in the cartesian case is equivalent to the z component
! in the cylindrical polar coordinate case.
!      
      DO el=1,numelm
  templocal_OUTx=0d0
  templocal_OUTy=0d0
!   temp_beta=param_beta_a(el) - param_beta
  DO l=0,N
    ll=l*NP1
  DO k=0,N
    kl=k+ll

    temp1=0d0
    temp2=0d0
    temp3=0d0
    temp4=0d0

    DO i=0,N
      ki=k+i*NP1
      il=i+ll

! integrating Div(G) 
! x-comp sums:
! used i as dummy variable, workings may differ.
      temp1 = temp1 + nodeCoord(mapg(il,el),2)*w(i)*d(i,k)*(devss_soln_xx(il,el)*dyde(i,l,el) - devss_soln_xy(il,el)*dxde(i,l,el)) 
      temp2 = temp2 + nodeCoord(mapg(ki,el),2)*w(i)*d(i,l)*(devss_soln_xy(ki,el)*dxdp(k,i,el) - devss_soln_xx(ki,el)*dydp(k,i,el))
! y-comp sums:
      temp3 = temp3 + nodeCoord(mapg(il,el),2)*w(i)*d(i,k)*(devss_soln_xy(il,el)*dyde(i,l,el) - devss_soln_yy(il,el)*dxde(i,l,el))
      temp4 = temp4 + nodeCoord(mapg(ki,el),2)*w(i)*d(i,l)*(devss_soln_yy(ki,el)*dxdp(k,i,el) - devss_soln_xy(ki,el)*dydp(k,i,el))
            
    ENDDO

    IF (.not.bdflag(1,mapg(kl,el))) THEN
    templocal_OUTx(kl) = (w(l)*temp1 + w(k)*temp2)
    ENDIF
! Add the extra rr and theta,theta contribution to divergence in y(r)-comp case:
    IF (.not.bdflag(2,mapg(kl,el))) THEN
    templocal_OUTy(kl) = (devss_soln_zz(kl,el)*jac(k,l,el)*w(k)*w(l) + w(l)*temp3 + w(k)*temp4)
    ENDIF

  ENDDO
  ENDDO
  tempglob=0d0
  CALL vecglobalprolongation(templocal_OUTx,el,tempglob)
  devss_contrib_x = devss_contrib_x - tempglob
  tempglob=0d0
  CALL vecglobalprolongation(templocal_OUTy,el,tempglob)
  devss_contrib_y = devss_contrib_y - tempglob
      ENDDO
    ELSE
      print*, 'Error: No co-ordinate system specified...'
      print*, 'Stopping'
      STOP
    ENDIF
    
! Boundary Contribution from non-dirichlet edge nodes (i.e. the axis of symmetry)
! For the axis of symmetry, only the x/z-component has a neumann BC. The test function is zero in other cases
!     DO el=1,numelm
!       edge=axisymm_edge(el)
!       IF (edge.ne.0) THEN
!   DO i=0,N
!     ij=local_edge_node(i,edge)
!     k=mapg(ij,el)
!     normx=norm_to_edge_node(1,i,edge,el)
!     normy=norm_to_edge_node(2,i,edge,el)
!     IF (.not.bdflag(1,k)) THEN  
!       devss_contrib_x(k) = devss_contrib_x(k) + ( devss_soln_xx(ij,el)*normx + (devss_soln_xy(ij,el) + devss_soln_xy(ij,el))*normy )*jac_on_edge(edge,el)*w(i)
!     ENDIF
!     IF (.not.bdflag(2,k)) THEN
!       devss_contrib_y(k) = devss_contrib_y(k) + ( (devss_soln_xy(ij,el) + devss_soln_xy(ij,el))*normx + devss_soln_yy(ij,el)*normy )*jac_on_edge(edge,el)*w(i)
!     ENDIF
!   ENDDO
!       ENDIF
!     ENDDO
      
  END SUBROUTINE calc_devss_contribution


END MODULE devss_module