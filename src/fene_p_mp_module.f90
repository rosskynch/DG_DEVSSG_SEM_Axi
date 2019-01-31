!    DG_DEVSSG_SEM
!    Copyright (C) 2008-2019 Ross M Kynch
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

MODULE fene_p_mp_module
  USE shared_data
  USE functions_module
  USE IO_module
  USE viscoelastic_module

  IMPLICIT NONE
  CONTAINS

  SUBROUTINE applyElasticStress_FENE_PMP
    IMPLICIT NONE
    stress_cont_to_stokes_xNm1=stress_cont_to_stokes_x
    stress_cont_to_stokes_yNm1=stress_cont_to_stokes_y

    CALL calcStress_weakform_FENE_PMP
    CALL integrate_divergence_of_local_stress

! EXJ scheme:
    f_x = f_x + time_beta_0*stress_cont_to_stokes_x + time_beta_1*stress_cont_to_stokes_xNm1
    f_y = f_y + time_beta_0*stress_cont_to_stokes_y + time_beta_1*stress_cont_to_stokes_yNm1

  END SUBROUTINE applyElasticStress_FENE_PMP

  SUBROUTINE calcStress_weakform_FENE_PMP
    IMPLICIT NONE
    INTEGER :: i,j,k,l,ij,p,pj,ip,el,kl,meh,fromel,edge,fromij,localp,&
! LAPACK bits:
    info,ipiv(4),iwork(4),ii,jj,sum_counter

    DOUBLE PRECISION :: a11,a22,a33,a44,a12,a21,a23,a32,idetA,&
      r1,r2,r3,r4,&
      i11,i12,i13,i21,i22,i23,i31,i32,i33,i44,&
      sum_temp,total_points,sumxx,sumxy,sumyy,sumzz,&
      temp12sq,temp_giesekus_const,tempVx(1:nptot),&
! LAPACK bits:
      temp_matrix(4,4),temp_rhs(4),temp1,temp2,temp3,temp4,&
      work(16),anorm,rcond

    DOUBLE PRECISION :: convective_contrib_xx(0:NP1SQM1,numelm),&
      convective_contrib_xy(0:NP1SQM1,numelm),&
      convective_contrib_yy(0:NP1SQM1,numelm),&
      convective_contrib_zz(0:NP1SQM1,numelm),&
      temp_contrib_xx(0:NP1SQM1,numelm),&
      temp_contrib_xy(0:NP1SQM1,numelm),&
      temp_contrib_yy(0:NP1SQM1,numelm),&
      temp_contrib_zz(0:NP1SQM1,numelm),&
      temp_contrib_xxNm1(0:NP1SQM1,numelm),&
      temp_contrib_xyNm1(0:NP1SQM1,numelm),&
      temp_contrib_yyNm1(0:NP1SQM1,numelm),&
      temp_contrib_zzNm1(0:NP1SQM1,numelm),&
      tempTxx(0:NP1SQM1,numelm),&
      tempTxy(0:NP1SQM1,numelm),&
      tempTyy(0:NP1SQM1,numelm),&
      tempTzz(0:NP1SQM1,numelm),&
      tempTxxNext(0:NP1SQM1,numelm),&
      tempTxyNext(0:NP1SQM1,numelm),&
      tempTyyNext(0:NP1SQM1,numelm),&
      tempTzzNext(0:NP1SQM1,numelm)

    temp_giesekus_const=param_giesekus/(1d0-param_beta)
    sum_temp=1d0
    temp_contrib_xx=0d0
    temp_contrib_xy=0d0
    temp_contrib_yy=0d0
    temp_contrib_zz=0d0
    temp_contrib_xxNm1=0d0
    temp_contrib_xyNm1=0d0
    temp_contrib_yyNm1=0d0


! Will iterate with tempTpq as the solution at time N+1, and localTpq as Tpq at time N, and localTpqNm1 as Txx at time N-1
! We then update Txx with the converged value.
    tempTxx = localTxx
    tempTxy = localTxy
    tempTyy = localTyy
    tempTzz = localTzz
    tempTxxNext = localTxx
    tempTxyNext = localTxy
    tempTyyNext = localTyy
    tempTzzNext = localTzz


    IF (coordflag.eq.0) THEN
! CARTESIAN CASE !
! Not implemented. Will need to base it on the OldroydB/Giesekus version within the main viscoelastic_module.
      print*, 'Error: The FENE-P-MP model is not yet implemented for Cartesian co-ordinates...'
      print*, 'Stopping'
      STOP
    ELSEIF (coordflag.eq.1) THEN 
! CYLINDERICAL POLAR (AXISYMMETRIC) CASE !
!
! NOTE:
! The x co-ordinate in the cartesian case is equivalent to the z component
! and y is equivalent to r component.
! in the cylindrical polar coordinate case.
!
! Semi-Implicit iteration.
      IF (param_iterative_convection) THEN
        sum_counter=0
        DO WHILE (sum_temp.gt.1d-9) 
          sum_counter=sum_counter+1
          tempVx=V_x-mesh_velocity
          CALL calc_axisymm_convective_term(tempVx,V_y,tempTxx,tempTxy,tempTyy, tempTzz, &
            convective_contrib_xx,convective_contrib_xy, &
            convective_contrib_yy,convective_contrib_zz)

          DO el=1,numelm
! Now loop over all stress nodes and compute solutions or apply boundary conditions, etc.
            DO ij=0,NP1SQM1
              i=mapg(ij,el)
              IF (inflowflag(i)) THEN 
!  Apply inflow boundary conditions on stress.
                tempTxxNext(ij,el) = boundary_stress_xx(i)
                tempTxyNext(ij,el) = boundary_stress_xy(i)
                tempTyyNext(ij,el) = boundary_stress_yy(i)
                tempTzzNext(ij,el) = boundary_stress_yy(i)
                CYCLE
              ENDIF

! Calculate entries of the matrix for OldroydB:

! DEVSS-G:
! Where G tensor is used for deformation terms as well as Strain Rate.

              a11 = (1d0 + Wetime_constant1 - 2d0*We*localGradUxx(ij,el))
              a22 = (1d0 + Wetime_constant1 - We*(localGradUxx(ij,el) + localGradUyy(ij,el)))
              a33 = (1d0 + Wetime_constant1 - 2d0*We*localGradUyy(ij,el)) 
              a44 = (1d0 + Wetime_constant1 - 2d0*We*localGradUzz(ij,el))! - 2d0*We*V_y(k)*jac(i,j,el)*w(i)*w(j) !

              a12 = -2d0*We*localGradUyx(ij,el)
              a21 = -We*localGradUxy(ij,el)
              a23 = -We*localGradUyx(ij,el)
              a32 = -2d0*We*localGradUxy(ij,el)

! Calculate RHS entries
! BDFJ:
              temp12sq = tempTxy(ij,el)**2

              r1 = 2d0*(1d0-param_beta)*localGradUxx(ij,el) &
                + Wetime_constant2*( time_alpha_0*localTxx(ij,el) + time_alpha_1*localTxxNm1(ij,el) ) &! + time_alpha_2*localTxxNm2(ij,el) ) &
                - We*convective_contrib_xx(ij,el) &
                - temp_giesekus_const*(tempTxx(ij,el)**2 + temp12sq)

              r2 = (1d0-param_beta)*( localGradUyx(ij,el) + localGradUxy(ij,el) ) &
                + Wetime_constant2*( time_alpha_0*localTxy(ij,el) + time_alpha_1*localTxyNm1(ij,el) ) &! + time_alpha_2*localTxyNm2(ij,el) ) &
                - We*convective_contrib_xy(ij,el) &
                - temp_giesekus_const*(tempTxx(ij,el)*tempTxy(ij,el) + tempTxy(ij,el)*tempTyy(ij,el))

              r3 = 2d0*(1d0-param_beta)*localGradUyy(ij,el) &
                + Wetime_constant2*( time_alpha_0*localTyy(ij,el) + time_alpha_1*localTyyNm1(ij,el) ) &! + time_alpha_2*localTyyNm2(ij,el) ) &
                - We*convective_contrib_yy(ij,el) &
                - temp_giesekus_const*(temp12sq + tempTyy(ij,el)**2)

              r4 =  2d0*(1d0-param_beta)*localGradUzz(ij,el) &
                + Wetime_constant2*( time_alpha_0*localTzz(ij,el) + time_alpha_1*localTzzNm1(ij,el) ) &! + time_alpha_2*localTzzNm2(ij,el) ) &
                - We*convective_contrib_zz(ij,el) &
                - temp_giesekus_const*tempTzz(ij,el)**2

! Using LAPACK to solve 4x4:
              temp_matrix=0d0
              temp_matrix(1,1)=a11
              temp_matrix(1,2)=a12
              temp_matrix(2,1)=a21
              temp_matrix(2,2)=a22
              temp_matrix(2,3)=a23
              temp_matrix(3,2)=a32
              temp_matrix(3,3)=a33
              temp_matrix(4,4)=a44
              temp_rhs(1)=r1
              temp_rhs(2)=r2
              temp_rhs(3)=r3
              temp_rhs(4)=r4

              call dgetrf( 4, 4, temp_matrix, 4, ipiv, info )
              IF (info.ne.0) THEN
                write(*,*) 'Error in calcStress_weakform:',el,info
                STOP
              ENDIF

              call dgetrs( 'N', 4, 1, temp_matrix, 4, ipiv, temp_rhs, 4, info )  
              IF (info.ne.0) THEN
                write(*,*) 'Error in calcStress_weakform:',el,info
                STOP
              ENDIF

! Copy solution from lapack:
              tempTxxNext(ij,el) = temp_rhs(1)
              tempTxyNext(ij,el) = temp_rhs(2)
              tempTyyNext(ij,el) = temp_rhs(3)
              tempTzzNext(ij,el) = temp_rhs(4)

            ENDDO
          ENDDO

          sum_temp=0d0
          sumxx=0d0
          sumxy=0d0
          sumyy=0d0
          sumzz=0d0
          DO el=1,numelm
            DO ij=0,NP1SQM1
              sumxx=sumxx+abs(tempTxxNext(ij,el)-tempTxx(ij,el))
              sumxy=sumxy+abs(tempTxyNext(ij,el)-tempTxy(ij,el))
              sumyy=sumyy+abs(tempTyyNext(ij,el)-tempTyy(ij,el))
              sumzz=sumzz+abs(tempTzzNext(ij,el)-tempTzz(ij,el))
            ENDDO
          ENDDO
          sumxx=sumxx/(numelm*NP1SQ)
          sumxy=sumxy/(numelm*NP1SQ)
          sumyy=sumyy/(numelm*NP1SQ)
          sumzz=sumzz/(numelm*NP1SQ)

          sum_temp=(sumxx+sumxy+sumyy+sumzz)/4d0

          IF(sum_temp.gt.1d10.or.sum_counter.gt.1000) THEN
            print*,'Iterative scheme in constitutive equation failed to converge! Sum counter = ',sum_counter
            STOP
          ENDIF

          tempTxx = tempTxxNext
          tempTxy = tempTxyNext
          tempTyy = tempTyyNext
          tempTzz = tempTzzNext

        ENDDO
      ELSE
! EXJ VERSION
        tempVx=V_x-mesh_velocity
        CALL calc_axisymm_convective_term(tempVx,V_y,localTxx,localTxy,localTyy, localTzz, &
          temp_contrib_xx,temp_contrib_xy, &
          temp_contrib_yy,temp_contrib_zz)
        IF (param_time_order.eq.2) THEN
          tempVx=V_xNm1-mesh_velocityNm1
          CALL calc_axisymm_convective_term(tempVx,V_yNm1,localTxxNm1,localTxyNm1,localTyyNm1, localTzzNm1, &
            temp_contrib_xxNm1,temp_contrib_xyNm1, &
            temp_contrib_yyNm1,temp_contrib_zzNm1)
        ENDIF
              
        convective_contrib_xx = time_beta_0*temp_contrib_xx + time_beta_1*temp_contrib_xxNm1
        convective_contrib_xy = time_beta_0*temp_contrib_xy + time_beta_1*temp_contrib_xyNm1
        convective_contrib_yy = time_beta_0*temp_contrib_yy + time_beta_1*temp_contrib_yyNm1
        convective_contrib_zz = time_beta_0*temp_contrib_zz + time_beta_1*temp_contrib_zzNm1

        DO el=1,numelm

! Now loop over all stress nodes and compute solutions or apply boundary conditions, etc.
          DO ij=0,NP1SQM1
            i=mapg(ij,el)
            IF (inflowflag(i)) THEN 
!  Apply inflow boundary conditions on stress.
              tempTxxNext(ij,el) = boundary_stress_xx(i)
              tempTxyNext(ij,el) = boundary_stress_xy(i)
              tempTyyNext(ij,el) = boundary_stress_yy(i)
              tempTzzNext(ij,el) = boundary_stress_yy(i)
              CYCLE
            ENDIF
! Calculate entries of the matrix for OldroydB:

! DEVSS-G:
! Where G tensor is used for deformation terms as well as Strain Rate.

            a11 = (1d0 + Wetime_constant1 - 2d0*We*localGradUxx(ij,el))
            a22 = (1d0 + Wetime_constant1 - We*(localGradUxx(ij,el) + localGradUyy(ij,el)))
            a33 = (1d0 + Wetime_constant1 - 2d0*We*localGradUyy(ij,el)) 
            a44 = (1d0 + Wetime_constant1 - 2d0*We*localGradUzz(ij,el))

            a12 = -2d0*We*localGradUyx(ij,el)
            a21 = -We*localGradUxy(ij,el)
            a23 = -We*localGradUyx(ij,el)
            a32 = -2d0*We*localGradUxy(ij,el)

! Calculate RHS entries
! BDFJ:
            temp12sq = time_beta_0*localTxy(ij,el)**2 + time_beta_1*localTxyNm1(ij,el)**2

            r1 = 2d0*(1d0-param_beta)*localGradUxx(ij,el) &
              + Wetime_constant2*( time_alpha_0*localTxx(ij,el) + time_alpha_1*localTxxNm1(ij,el) ) &
              - We*convective_contrib_xx(ij,el) &
              - temp_giesekus_const*(time_beta_0*localTxx(ij,el)**2 + time_beta_1*localTxxNm1(ij,el)**2 + temp12sq)

            r2 = (1d0-param_beta)*( localGradUyx(ij,el) + localGradUxy(ij,el) ) &
              + Wetime_constant2*( time_alpha_0*localTxy(ij,el) + time_alpha_1*localTxyNm1(ij,el) ) &
              - We*convective_contrib_xy(ij,el) &
              - temp_giesekus_const*(time_beta_0*(localTxx(ij,el)*localTxy(ij,el) + localTxy(ij,el)*localTyy(ij,el)) &
              + time_beta_1*(localTxxNm1(ij,el)*localTxyNm1(ij,el) + localTxyNm1(ij,el)*localTyyNm1(ij,el)) )

            r3 = 2d0*(1d0-param_beta)*localGradUyy(ij,el) &
              + Wetime_constant2*( time_alpha_0*localTyy(ij,el) + time_alpha_1*localTyyNm1(ij,el) ) &
              - We*convective_contrib_yy(ij,el) &
              - temp_giesekus_const*(temp12sq + time_beta_0*localTyy(ij,el)**2 + time_beta_1*localTyyNm1(ij,el)**2)

            r4 =  2d0*(1d0-param_beta)*localGradUzz(ij,el) &
              + Wetime_constant2*( time_alpha_0*localTzz(ij,el) + time_alpha_1*localTzzNm1(ij,el) ) &
              - We*convective_contrib_zz(ij,el) &
              - temp_giesekus_const*(time_beta_0*localTzz(ij,el)**2 + time_beta_1*localTzzNm1(ij,el)**2)


! Using LAPACK to solve 4x4:
            temp_matrix=0d0
            temp_matrix(1,1)=a11
            temp_matrix(1,2)=a12
            temp_matrix(2,1)=a21
            temp_matrix(2,2)=a22
            temp_matrix(2,3)=a23
            temp_matrix(3,2)=a32
            temp_matrix(3,3)=a33
            temp_matrix(4,4)=a44
            temp_rhs(1)=r1
            temp_rhs(2)=r2
            temp_rhs(3)=r3
            temp_rhs(4)=r4


            call dgetrf( 4, 4, temp_matrix, 4, ipiv, info )
            IF (info.ne.0) THEN
              write(*,*) 'Error in calcStress_weakform:',el,info
              STOP
            ENDIF
    
            call dgetrs( 'N', 4, 1, temp_matrix, 4, ipiv, temp_rhs, 4, info )  
            IF (info.ne.0) THEN
              write(*,*) 'Error in calcStress_weakform:',el,info
              STOP
            ENDIF

! Copy solution from lapack:
            tempTxxNext(ij,el) = temp_rhs(1)
            tempTxyNext(ij,el) = temp_rhs(2)
            tempTyyNext(ij,el) = temp_rhs(3)
            tempTzzNext(ij,el) = temp_rhs(4)
  
          ENDDO
        ENDDO
        tempTxx=tempTxxNext
        tempTxy=tempTxyNext
        tempTyy=tempTyyNext
        tempTzz=tempTzzNext
      ENDIF

! Update stored valeus of Tpq 
! ie move time along by one for these, so that localTpq now hold the current stress values
      localTxxNm2=localTxxNm1
      localTxyNm2=localTxyNm1
      localTyyNm2=localTyyNm1
      localTzzNm2=localTzzNm1
      
      localTxxNm1=localTxx
      localTxyNm1=localTxy
      localTyyNm1=localTyy
      localTzzNm1=localTzz
      
      localTxx=tempTxx
      localTxy=tempTxy
      localTyy=tempTyy
      localTzz=tempTzz

    ELSE
      print*, 'Error: No co-ordinate system specified...'
      print*, 'Stopping'
      STOP
    ENDIF
  END SUBROUTINE calcStress_weakform_FENE_PMP