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

MODULE viscoelastic_module
  USE shared_data
  USE functions_module
  USE IO_module
!   USE result_analysis
  
  IMPLICIT NONE
  CONTAINS
  
  SUBROUTINE applyElasticStress!(in_gradUxx,in_gradUyx,in_gradUxy,in_gradUyy,in_gradUzz)
    IMPLICIT NONE
!     DOUBLE PRECISION, DIMENSION(0:NP1SQM1,numelm), INTENT(IN) :: in_gradUxx,in_gradUyx,in_gradUxy,in_gradUyy,in_gradUzz

!     stress_cont_to_stokes_xNm2=stress_cont_to_stokes_xNm1
!     stress_cont_to_stokes_yNm2=stress_cont_to_stokes_yNm1 
    stress_cont_to_stokes_xNm1=stress_cont_to_stokes_x
    stress_cont_to_stokes_yNm1=stress_cont_to_stokes_y
    
    CALL calcStress_weakform
    CALL integrate_divergence_of_local_stress

! EXJ scheme:
    f_x = f_x + time_beta_0*stress_cont_to_stokes_x + time_beta_1*stress_cont_to_stokes_xNm1! + time_beta_2*stress_cont_to_stokes_xNm2
    f_y = f_y + time_beta_0*stress_cont_to_stokes_y + time_beta_1*stress_cont_to_stokes_yNm1! + time_beta_2*stress_cont_to_stokes_yNm2

  END SUBROUTINE applyElasticStress

 
  SUBROUTINE calc_2d_convective_term(u_in,v_in,in_xx,in_xy,in_yy,out_xx,out_xy,out_yy)
    IMPLICIT NONE
    INTEGER :: p,q,i,k,ij,kl,el,edge,fromel,fromij
    DOUBLE PRECISION :: diffx_xx,diffy_xx, &
      diffx_xy,diffy_xy, &
      diffx_yy,diffy_yy, &
      temp
    DOUBLE PRECISION, INTENT(IN) :: u_in(nptot),v_in(nptot),in_xx(0:NP1SQM1,numelm),in_xy(0:NP1SQM1,numelm),in_yy(0:NP1SQM1,numelm)
    DOUBLE PRECISION, INTENT(OUT) :: out_xx(0:NP1SQM1,numelm),out_xy(0:NP1SQM1,numelm),out_yy(0:NP1SQM1,numelm)
    
    out_xx=0d0
    out_xy=0d0
    out_yy=0d0
    DO el=1,numelm
! Integral for convective contribution.
      DO kl=0,NP1SQM1
        k=mapg(kl,el)
        
        diffx_xx=0d0
        diffy_xx=0d0
        diffx_xy=0d0
        diffy_xy=0d0
        diffx_yy=0d0
        diffy_yy=0d0

        DO ij=0,NP1SQM1
          diffx_xx = diffx_xx + diff_x(kl,ij,el)*in_xx(ij,el)
          diffy_xx = diffy_xx + diff_y(kl,ij,el)*in_xx(ij,el)
          
          diffx_xy = diffx_xy + diff_x(kl,ij,el)*in_xy(ij,el)
          diffy_xy = diffy_xy + diff_y(kl,ij,el)*in_xy(ij,el)
          
          diffx_yy = diffx_yy + diff_x(kl,ij,el)*in_yy(ij,el)
          diffy_yy = diffy_yy + diff_y(kl,ij,el)*in_yy(ij,el)
        ENDDO
        out_xx(kl,el) = u_in(k)*diffx_xx + v_in(k)*diffy_xx
        out_xy(kl,el) = u_in(k)*diffx_xy + v_in(k)*diffy_xy
        out_yy(kl,el) = u_in(k)*diffx_yy + v_in(k)*diffy_yy
      ENDDO
! Boundary integral for DG contribution:
      DO edge=1,4
        DO i=0,N
          ij=local_edge_node(i,edge)
          p=local_ij_to_i(ij)
          q=local_ij_to_j(ij)  
          k=mapg(ij,el)
          fromel=upwind_local_edge_node(i,edge,el)
          fromij=global_to_local_map(k,fromel)
          IF (fromel.ne.el) THEN
            temp = (u_in(k)*norm_to_edge_node(1,i,edge,el) + v_in(k)*norm_to_edge_node(2,i,edge,el))*w(i)*jac_on_edge(edge,el)/(w(p)*w(q)*jac(p,q,el))
            
            out_xx(ij,el) = out_xx(ij,el) + temp*(in_xx(fromij,fromel) - in_xx(ij,el))
            out_xy(ij,el) = out_xy(ij,el) + temp*(in_xy(fromij,fromel) - in_xy(ij,el))
            out_yy(ij,el) = out_yy(ij,el) + temp*(in_yy(fromij,fromel) - in_yy(ij,el))
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE calc_2d_convective_term
   
  
  SUBROUTINE calc_axisymm_convective_term(u_in,v_in,in_xx,in_xy,in_yy,in_zz,out_xx,out_xy,out_yy,out_zz)
    IMPLICIT NONE
    INTEGER :: p,q,i,k,ij,kl,el,edge,fromel,fromij
    DOUBLE PRECISION :: diffx_xx,diffy_xx, &
      diffx_xy,diffy_xy, &
      diffx_yy,diffy_yy, &
      diffx_zz,diffy_zz, &
      temp
    
    DOUBLE PRECISION, INTENT(IN) :: u_in(nptot),v_in(nptot),&
            in_xx(0:NP1SQM1,numelm),in_xy(0:NP1SQM1,numelm),&
            in_yy(0:NP1SQM1,numelm),in_zz(0:NP1SQM1,numelm)
    DOUBLE PRECISION, INTENT(OUT) :: out_xx(0:NP1SQM1,numelm),out_xy(0:NP1SQM1,numelm),&
             out_yy(0:NP1SQM1,numelm),out_zz(0:NP1SQM1,numelm)
    
    out_xx=0d0
    out_xy=0d0
    out_yy=0d0
    out_zz=0d0
    DO el=1,numelm
! Integral for convective term.
      DO kl=0,NP1SQM1
        k=mapg(kl,el)
  
        diffx_xx=0d0
        diffy_xx=0d0
        diffx_xy=0d0
        diffy_xy=0d0
        diffx_yy=0d0
        diffy_yy=0d0
        diffx_zz=0d0
        diffy_zz=0d0
      
        DO ij=0,NP1SQM1
          diffx_xx = diffx_xx + diff_x(kl,ij,el)*in_xx(ij,el)
          diffy_xx = diffy_xx + diff_y(kl,ij,el)*in_xx(ij,el)
          
          diffx_xy = diffx_xy + diff_x(kl,ij,el)*in_xy(ij,el)
          diffy_xy = diffy_xy + diff_y(kl,ij,el)*in_xy(ij,el)
          
          diffx_yy = diffx_yy + diff_x(kl,ij,el)*in_yy(ij,el)
          diffy_yy = diffy_yy + diff_y(kl,ij,el)*in_yy(ij,el)
          
          diffx_zz = diffx_zz + diff_x(kl,ij,el)*in_zz(ij,el)
          diffy_zz = diffy_zz + diff_y(kl,ij,el)*in_zz(ij,el)
        ENDDO
        out_xx(kl,el) = u_in(k)*diffx_xx + v_in(k)*diffy_xx
        out_xy(kl,el) = u_in(k)*diffx_xy + v_in(k)*diffy_xy
        out_yy(kl,el) = u_in(k)*diffx_yy + v_in(k)*diffy_yy
        out_zz(kl,el) = u_in(k)*diffx_zz + v_in(k)*diffy_zz
      ENDDO

! Boundary integral for DG contribution:
      DO edge=1,4
        DO i=0,N
          ij=local_edge_node(i,edge)
          p=local_ij_to_i(ij)
          q=local_ij_to_j(ij)  
          k=mapg(ij,el)
          fromel=upwind_local_edge_node(i,edge,el)
          fromij=global_to_local_map(k,fromel)
          IF (fromel.ne.el) THEN
            temp = (u_in(k)*norm_to_edge_node(1,i,edge,el) + v_in(k)*norm_to_edge_node(2,i,edge,el))*w(i)*jac_on_edge(edge,el)/(w(p)*w(q)*jac(p,q,el))
            
            
            out_xx(ij,el) = out_xx(ij,el) + temp*(in_xx(fromij,fromel) - in_xx(ij,el))
            out_xy(ij,el) = out_xy(ij,el) + temp*(in_xy(fromij,fromel) - in_xy(ij,el))
            out_yy(ij,el) = out_yy(ij,el) + temp*(in_yy(fromij,fromel) - in_yy(ij,el))
            out_zz(ij,el) = out_zz(ij,el) + temp*(in_zz(fromij,fromel) - in_zz(ij,el)) 
      
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE calc_axisymm_convective_term
  
  SUBROUTINE calc_1d_convective_term(u_in,v_in,in_xx,in_xy,in_yy,in_zz,out_xx,out_xy,out_yy,out_zz)
    INTEGER :: el,edge,i,ij,k,p,localp,fromij,fromel
    DOUBLE PRECISION :: temp1,temp2,temp3,temp4,temp_sign=1d0
    DOUBLE PRECISION, INTENT(IN) :: u_in(nptot),v_in(nptot),&
            in_xx(0:NP1SQM1,numelm),in_xy(0:NP1SQM1,numelm),&
            in_yy(0:NP1SQM1,numelm),in_zz(0:NP1SQM1,numelm)
    DOUBLE PRECISION, INTENT(OUT) :: out_xx(0:N,numelm),out_xy(0:N,numelm),&
             out_yy(0:N,numelm),out_zz(0:N,numelm)

! temp_sign=1d0
    out_xx=0d0
    out_xy=0d0
    out_yy=0d0
    out_zz=0d0
! Deal with axisymmetric wall using a 1-D treatment of the integrals. I.e. We write over the previous values with the
! 1-D integral.
! Note: These nodes may be solved independently of the rest of the stress nodes.
    DO el=1,numelm
      edge=axisymm_edge(el)
      IF (edge.eq.4.or.edge.eq.3) then
        temp_sign=-1d0
      ELSE
        temp_sign=1d0
      ENDIF
! If the edge is an axisymmetric wall, then we deal with it as an independent 1-D integral, although still using DG:
      IF (edge.ne.0) THEN
    
! As it is an edge integral.. include the 1-D boundary integral for the convection term
! I think i=0,N is ok here. It may end up going backwards in some cases.. but I think for our meshes its ok!
        DO i=0,N
          temp1=0d0
          temp2=0d0
          temp3=0d0
          temp4=0d0
          ij = local_edge_node(i,edge)
          k = mapg(ij,el)
          
          DO p=0,N
            localp = local_edge_node(p,edge)
            temp1 = temp1 + in_xx(localp,el)*d(i,p)
            temp2 = temp2 + in_xy(localp,el)*d(i,p)
            temp3 = temp3 + in_yy(localp,el)*d(i,p)
            temp4 = temp4 + in_zz(localp,el)*d(i,p)
          ENDDO
          out_xx(i,el) = u_in(k)*temp1/jac_wallsymm_edge_1d(edge,el)
          out_xy(i,el) = u_in(k)*temp2/jac_wallsymm_edge_1d(edge,el)
          out_yy(i,el) = u_in(k)*temp3/jac_wallsymm_edge_1d(edge,el)
          out_zz(i,el) = u_in(k)*temp4/jac_wallsymm_edge_1d(edge,el)
        ENDDO
! Next include the DG treatment of the boundary terms, which only occur at the element edges (2 points per element).
        k = node(el,edge)
        ij = global_to_local_map(k,el)
        i=local_ij_to_edge_node(ij,edge)
        fromel = upwind_wallsymm_1d(1,edge,el)
        fromij = global_to_local_map(k,fromel)
    
        IF (fromel.ne.el) THEN
          out_xx(i,el) = out_xx(i,el) - temp_sign*u_in(k)*(in_xx(fromij,fromel)-in_xx(ij,el))/(w(i)*jac_wallsymm_edge_1d(edge,el))
          out_xy(i,el) = out_xy(i,el) - temp_sign*u_in(k)*(in_xy(fromij,fromel)-in_xy(ij,el))/(w(i)*jac_wallsymm_edge_1d(edge,el))
          out_yy(i,el) = out_yy(i,el) - temp_sign*u_in(k)*(in_yy(fromij,fromel)-in_yy(ij,el))/(w(i)*jac_wallsymm_edge_1d(edge,el))
          out_zz(i,el) = out_zz(i,el) - temp_sign*u_in(k)*(in_zz(fromij,fromel)-in_zz(ij,el))/(w(i)*jac_wallsymm_edge_1d(edge,el))
        ENDIF
        k = node(el,edge+1)
        ij = global_to_local_map(k,el)
        i=local_ij_to_edge_node(ij,edge+1)
        fromel = upwind_wallsymm_1d(2,edge,el)
        fromij = global_to_local_map(k,fromel)
        IF (fromel.ne.el) THEN
          out_xx(i,el) = out_xx(i,el) + temp_sign*u_in(k)*(in_xx(fromij,fromel)-in_xx(ij,el))/(w(i)*jac_wallsymm_edge_1d(edge,el))
          out_xy(i,el) = out_xy(i,el) + temp_sign*u_in(k)*(in_xy(fromij,fromel)-in_xy(ij,el))/(w(i)*jac_wallsymm_edge_1d(edge,el))
          out_yy(i,el) = out_yy(i,el) + temp_sign*u_in(k)*(in_yy(fromij,fromel)-in_yy(ij,el))/(w(i)*jac_wallsymm_edge_1d(edge,el))
          out_zz(i,el) = out_zz(i,el) + temp_sign*u_in(k)*(in_zz(fromij,fromel)-in_zz(ij,el))/(w(i)*jac_wallsymm_edge_1d(edge,el))
        ENDIF
      ENDIF
    ENDDO
  END SUBROUTINE calc_1d_convective_term
 
  SUBROUTINE calcStress_weakform
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
!       temp_contrib_xxNm2(0:NP1SQM1,numelm),&
!       temp_contrib_xyNm2(0:NP1SQM1,numelm),&
!       temp_contrib_yyNm2(0:NP1SQM1,numelm),&
!       temp_contrib_zzNm2(0:NP1SQM1,numelm),&
!       conv_1d_xx(0:N,numelm),&
!       conv_1d_xy(0:N,numelm),&
!       conv_1d_yy(0:N,numelm),&
!       conv_1d_zz(0:N,numelm),&
!       conv_1d_xxNm1(0:N,numelm),&
!       conv_1d_xyNm1(0:N,numelm),&
!       conv_1d_yyNm1(0:N,numelm),&
!       conv_1d_zzNm1(0:N,numelm),&
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
    temp_contrib_zzNm1=0d0
!     temp_contrib_xxNm2=0d0
!     temp_contrib_xyNm2=0d0
!     temp_contrib_yyNm2=0d0
!     temp_contrib_zzNm2=0d0



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

      IF (param_iterative_convection) THEN
! ITERATIVE VERSION
        sum_counter=0
        DO WHILE (sum_temp.gt.1d-9) ! Semi-Implicit iteration.
          sum_counter=sum_counter+1
          CALL calc_2d_convective_term(V_x,V_y,tempTxx,tempTxy,tempTyy, &
            convective_contrib_xx,convective_contrib_xy,convective_contrib_yy)

          DO el=1,numelm
            DO ij=0,NP1SQM1
              i=mapg(ij,el)

! Set inflow boundary conditions on stress.
!!! MUST BE CHANGED ACCORDING TO THE PROBLEM !!!
              IF (inflowflag(i)) THEN
                tempTxxNext(ij,el) = boundary_stress_xx(i)
                tempTxyNext(ij,el) = boundary_stress_xy(i)
                tempTyyNext(ij,el) = boundary_stress_yy(i)
                CYCLE
              ENDIF

! Calculate entries of the matrix arriving from OldroydB:
! DEVSS-G:
! Where G tensor is used for deformation terms as well as Strain Rate.
              a11 = (1d0 + Wetime_constant1 - 2d0*We*localGradUxx(ij,el))
              a22 = (1d0 + Wetime_constant1 - We*(localGradUxx(ij,el) + localGradUyy(ij,el)))
              a33 = (1d0 + Wetime_constant1 - 2d0*We*localGradUyy(ij,el))

              a12 = -2d0*We*localGradUyx(ij,el)
              a21 = -We*localGradUxy(ij,el)
              a23 = -We*localGradUyx(ij,el)
              a32 = -2d0*We*localGradUxy(ij,el)

! Calculate entries of the inverse:
              idetA=1d0/(a11*(a22*a33-a23*a32)-a12*a21*a33)
              i11 = idetA*(a22*a33-a23*a32)
              i12 = -idetA*a12*a33
              i13 = idetA*a12*a23
              i21 = -idetA*a21*a33
              i22 = idetA*a11*a33
              i23 = -idetA*a11*a23
              i31 = idetA*a21*a32
              i32 = -idetA*a11*a32
              i33 = idetA*(a11*a22-a12*a21)

! Calculate RHS entries
! BDFJ:
              temp12sq = tempTxy(ij,el)**2

              r1 = 2d0*(1d0-param_beta)*localGradUxx(ij,el) &
                + Wetime_constant2*(time_alpha_0*localTxx(ij,el) + time_alpha_1*localTxxNm1(ij,el) ) &! + time_alpha_2*localTxxNm2(ij,el) ) &
                - We*convective_contrib_xx(ij,el) &
                - temp_giesekus_const*(tempTxx(ij,el)**2 + temp12sq)

              r2 = (1d0-param_beta)*(localGradUyx(ij,el) + localGradUxy(ij,el)) &
                + Wetime_constant2*(time_alpha_0*localTxy(ij,el) + time_alpha_1*localTxyNm1(ij,el) ) &! + time_alpha_2*localTxyNm2(ij,el) ) &
                - We*convective_contrib_xy(ij,el) &
                - temp_giesekus_const*(tempTxx(ij,el)*tempTxy(ij,el) + tempTxy(ij,el)*tempTyy(ij,el))

              r3 = 2d0*(1d0-param_beta)*localGradUyy(ij,el) &
                + Wetime_constant2*(time_alpha_0*localTyy(ij,el) + time_alpha_1*localTyyNm1(ij,el) ) &! + time_alpha_2*localTyyNm2(ij,el) ) &
                - We*convective_contrib_yy(ij,el) &
                - temp_giesekus_const*(temp12sq + tempTyy(ij,el)**2)
 
! Calculate stress components from Inverse*RHS
              tempTxxNext(ij,el) = i11*r1 + i12*r2 + i13*r3 
              tempTxyNext(ij,el) = i21*r1 + i22*r2 + i23*r3
              tempTyyNext(ij,el) = i31*r1 + i32*r2 + i33*r3
            ENDDO
          ENDDO
          
          sum_temp=0d0
          DO el=1,numelm
            DO ij=0,NP1SQM1
              sum_temp=sum_temp + abs(tempTxxNext(ij,el)-tempTxx(ij,el)) + &
                abs(tempTxyNext(ij,el)-tempTxy(ij,el)) + &
                abs(tempTyyNext(ij,el)-tempTyy(ij,el))
            ENDDO
          ENDDO
          sum_temp=sum_temp/(numelm*3d0*NP1SQ)
          
          IF(sum_temp.gt.1d10.or.sum_counter.gt.1000) THEN
            print*,'Iterative scheme in constitutive equation failed to converge! Sum counter = ',sum_counter
            STOP
          ENDIF
          
          tempTxx = tempTxxNext
          tempTxy = tempTxyNext
          tempTyy = tempTyyNext
        ENDDO
      ELSE
! EXJ VERSION

!     tempTxx=time_beta_0*localTxx + time_beta_1*localTxxNm1! + time_beta_2*localTxxNm2
!     tempTxy=time_beta_0*localTxy + time_beta_1*localTxyNm1! + time_beta_2*localTxyNm2
!     tempTyy=time_beta_0*localTyy + time_beta_1*localTyyNm1! + time_beta_2*localTyyNm2
!     tempVx=time_beta_0*V_x + time_beta_1*V_xNm1! + time_beta_2*Vx_Nm2
!     tempVy=time_beta_0*V_y + time_beta_1*V_yNm1! + time_beta_2*Vy_Nm2
    
        CALL calc_2d_convective_term(V_x,V_y,localTxx,localTxy,localTyy, &
          temp_contrib_xx,temp_contrib_xy,temp_contrib_yy)
        IF (param_time_order.eq.2) THEN
      
! Should we use V_x/y or V_x/yNm1 here???
          CALL calc_2d_convective_term(V_xNm1,V_yNm1,localTxxNm1,localTxyNm1,localTyyNm1, &
            temp_contrib_xxNm1,temp_contrib_xyNm1,temp_contrib_yyNm1)
! REMOVED for now:
!     ELSEIF (param_time_order.eq.3) THEN
!       CALL calc_2d_convective_term(V_xNm1,V_yNm1,localTxxNm1,localTxyNm1,localTyyNm1, &
!       temp_contrib_xxNm1,temp_contrib_xyNm1,temp_contrib_yyNm1)
!       
!       CALL calc_2d_convective_term(V_xNm2,V_yNm2,localTxxNm2,localTxyNm2,localTyyNm2, &
!       temp_contrib_xxNm2,temp_contrib_xyNm2,temp_contrib_yyNm2)
        ENDIF
      
        convective_contrib_xx = time_beta_0*temp_contrib_xx + time_beta_1*temp_contrib_xxNm1! + time_beta_2*temp_contrib_xxNm2
        convective_contrib_xy = time_beta_0*temp_contrib_xy + time_beta_1*temp_contrib_xyNm1! + time_beta_2*temp_contrib_xyNm2
        convective_contrib_yy = time_beta_0*temp_contrib_yy + time_beta_1*temp_contrib_yyNm1! + time_beta_2*temp_contrib_yyNm2
      
        DO el=1,numelm
          DO ij=0,NP1SQM1
            i=mapg(ij,el)
      
! Set inflow boundary conditions on stress.
!!! MUST BE CHANGED ACCORDING TO THE PROBLEM !!!
            IF (inflowflag(i)) THEN
              tempTxxNext(ij,el) = boundary_stress_xx(i)
              tempTxyNext(ij,el) = boundary_stress_xy(i)
              tempTyyNext(ij,el) = boundary_stress_yy(i)
              CYCLE
            ENDIF
  
! Calculate entries of the matrix arriving from OldroydB:
! DEVSS-G:
! Where G tensor is used for deformation terms as well as Strain Rate.
            a11 = (1d0 + Wetime_constant1 - 2d0*We*localGradUxx(ij,el))
            a22 = (1d0 + Wetime_constant1 - We*(localGradUxx(ij,el) + localGradUyy(ij,el)))
            a33 = (1d0 + Wetime_constant1 - 2d0*We*localGradUyy(ij,el))
    
            a12 = -2d0*We*localGradUyx(ij,el)
            a21 = -We*localGradUxy(ij,el)
            a23 = -We*localGradUyx(ij,el)
            a32 = -2d0*We*localGradUxy(ij,el)
   
! Calculate entries of the inverse:
            idetA=1d0/(a11*(a22*a33-a23*a32)-a12*a21*a33)
            i11 = idetA*(a22*a33-a23*a32)
            i12 = -idetA*a12*a33
            i13 = idetA*a12*a23
            i21 = -idetA*a21*a33
            i22 = idetA*a11*a33
            i23 = -idetA*a11*a23
            i31 = idetA*a21*a32
            i32 = -idetA*a11*a32
            i33 = idetA*(a11*a22-a12*a21)

! Calculate RHS entries
! BDFJ:
            temp12sq = time_beta_0*localTxy(ij,el)**2 + time_beta_1*localTxyNm1(ij,el)**2
      
            r1 = 2d0*(1d0-param_beta)*localGradUxx(ij,el) &
              + Wetime_constant2*(time_alpha_0*localTxx(ij,el) + time_alpha_1*localTxxNm1(ij,el) ) &! + time_alpha_2*localTxxNm2(ij,el) ) &
              - We*convective_contrib_xx(ij,el) &
              - temp_giesekus_const*(time_beta_0*localTxx(ij,el)**2 + time_beta_1*localTxxNm1(ij,el)**2 + temp12sq)
                    
            r2 = (1d0-param_beta)*(localGradUyx(ij,el) + localGradUxy(ij,el)) &
              + Wetime_constant2*(time_alpha_0*localTxy(ij,el) + time_alpha_1*localTxyNm1(ij,el) ) &! + time_alpha_2*localTxyNm2(ij,el) ) &
              - We*convective_contrib_xy(ij,el) &
              - temp_giesekus_const*(time_beta_0*(localTxx(ij,el)*localTxy(ij,el) + localTxy(ij,el)*localTyy(ij,el)) &
              + time_beta_1*(localTxxNm1(ij,el)*localTxyNm1(ij,el) + localTxyNm1(ij,el)*localTyyNm1(ij,el)) )
            
            
            r3 = 2d0*(1d0-param_beta)*localGradUyy(ij,el) &
              + Wetime_constant2*(time_alpha_0*localTyy(ij,el) + time_alpha_1*localTyyNm1(ij,el) ) &! + time_alpha_2*localTyyNm2(ij,el) ) &
              - We*convective_contrib_yy(ij,el) &
              - temp_giesekus_const*(temp12sq + time_beta_0*localTyy(ij,el)**2 + time_beta_1*localTyyNm1(ij,el)**2)
 
! Calculate stress components from Inverse*RHS
            tempTxxNext(ij,el) = i11*r1 + i12*r2 + i13*r3 
            tempTxyNext(ij,el) = i21*r1 + i22*r2 + i23*r3
            tempTyyNext(ij,el) = i31*r1 + i32*r2 + i33*r3
          ENDDO
        ENDDO

        tempTxx=tempTxxNext
        tempTxy=tempTxyNext
        tempTyy=tempTyyNext

      ENDIF
      localTxxNm2 = localTxxNm1
      localTxyNm2 = localTxyNm1
      localTyyNm2 = localTyyNm1
    
      localTxxNm1 = localTxx
      localTxyNm1 = localTxy
      localTyyNm1 = localTyy
    
      localTxx = tempTxx
      localTxy = tempTxy
      localTyy = tempTyy
    
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
     
! Calculate entries of the inverse:
!               idetA=1d0/(a11*(a22*a33-a23*a32)-a12*a21*a33)
!               i11 = (a22*a33-a23*a32)
!               i12 = -a12*a33
!               i13 = a12*a23
!               i21 = -a21*a33
!               i22 = a11*a33
!               i23 = -a11*a23
!               i31 = a21*a32
!               i32 = -a11*a32
!               i33 = (a11*a22-a12*a21)

! Using LAPACK instead:
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
  
!               anorm=0d0
!               do ii=1,4
!                 temp1=temp_matrix(ii,1)
!                 do jj=2,4
!                   if (temp_matrix(ii,jj).gt.temp1) temp1=temp_matrix(ii,jj)
!                 enddo
!                 anorm=anorm+temp1
!               enddo
  
              call dgetrf( 4, 4, temp_matrix, 4, ipiv, info )
              IF (info.ne.0) THEN
                write(*,*) 'Error in calcStress_weakform:',el,info
                STOP
              ENDIF
  
    
!    call dgecon( '1', 4, temp_matrix, 4, anorm, rcond, work, iwork, info )
!    IF (info.ne.0) THEN
!   write(*,*) 'Error in calcStress_weakform:',el,info
!   STOP
!   ENDIF
!   IF (el.eq.2.and.ij.eq.273) print*,'rcond: ',rcond
  
              call dgetrs( 'N', 4, 1, temp_matrix, 4, ipiv, temp_rhs, 4, info )  
              IF (info.ne.0) THEN
                write(*,*) 'Error in calcStress_weakform:',el,info
                STOP
              ENDIF
! Calculate stress components from Inverse*RHS
!               tempTxxNext(ij,el) = idetA*(i11*r1 + i12*r2 + i13*r3)
!               tempTxyNext(ij,el) = idetA*(i21*r1 + i22*r2 + i23*r3)
!               tempTyyNext(ij,el) = idetA*(i31*r1 + i32*r2 + i33*r3)
!               tempTzzNext(ij,el) = r4/a44  

! LAPACK instead:
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
! REMOVED for now:
!         ELSEIF (param_time_order.eq.3) THEN 
!           CALL calc_axisymm_convective_term(V_xNm1,V_yNm1,localTxxNm1,localTxyNm1,localTyyNm1, localTzzNm1, &
!             temp_contrib_xxNm1,temp_contrib_xyNm1, &
!             temp_contrib_yyNm1,temp_contrib_zzNm1)
!           CALL calc_axisymm_convective_term(V_xNm2,V_yNm2,localTxxNm2,localTxyNm2,localTyyNm2, localTzzNm2, &
!             temp_contrib_xxNm2,temp_contrib_xyNm2, &
!             temp_contrib_yyNm2,temp_contrib_zzNm2)
        ENDIF
              
        convective_contrib_xx = time_beta_0*temp_contrib_xx + time_beta_1*temp_contrib_xxNm1! + time_beta_2*temp_contrib_xxNm2
        convective_contrib_xy = time_beta_0*temp_contrib_xy + time_beta_1*temp_contrib_xyNm1! + time_beta_2*temp_contrib_xyNm2
        convective_contrib_yy = time_beta_0*temp_contrib_yy + time_beta_1*temp_contrib_yyNm1! + time_beta_2*temp_contrib_yyNm2
        convective_contrib_zz = time_beta_0*temp_contrib_zz + time_beta_1*temp_contrib_zzNm1! + time_beta_2*temp_contrib_zzNm2

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
            temp12sq = time_beta_0*localTxy(ij,el)**2 + time_beta_1*localTxyNm1(ij,el)**2

            r1 = 2d0*(1d0-param_beta)*localGradUxx(ij,el) &
              + Wetime_constant2*( time_alpha_0*localTxx(ij,el) + time_alpha_1*localTxxNm1(ij,el) ) &! + time_alpha_2*localTxxNm2(ij,el) ) &
              - We*convective_contrib_xx(ij,el) &
              - temp_giesekus_const*(time_beta_0*localTxx(ij,el)**2 + time_beta_1*localTxxNm1(ij,el)**2 + temp12sq)

            r2 = (1d0-param_beta)*( localGradUyx(ij,el) + localGradUxy(ij,el) ) &
              + Wetime_constant2*( time_alpha_0*localTxy(ij,el) + time_alpha_1*localTxyNm1(ij,el) ) &! + time_alpha_2*localTxyNm2(ij,el) ) &
              - We*convective_contrib_xy(ij,el) &
              - temp_giesekus_const*(time_beta_0*(localTxx(ij,el)*localTxy(ij,el) + localTxy(ij,el)*localTyy(ij,el)) &
              + time_beta_1*(localTxxNm1(ij,el)*localTxyNm1(ij,el) + localTxyNm1(ij,el)*localTyyNm1(ij,el)) )

            r3 = 2d0*(1d0-param_beta)*localGradUyy(ij,el) &
              + Wetime_constant2*( time_alpha_0*localTyy(ij,el) + time_alpha_1*localTyyNm1(ij,el) ) &! + time_alpha_2*localTyyNm2(ij,el) ) &
              - We*convective_contrib_yy(ij,el) &
              - temp_giesekus_const*(temp12sq + time_beta_0*localTyy(ij,el)**2 + time_beta_1*localTyyNm1(ij,el)**2)

            r4 =  2d0*(1d0-param_beta)*localGradUzz(ij,el) &
              + Wetime_constant2*( time_alpha_0*localTzz(ij,el) + time_alpha_1*localTzzNm1(ij,el) ) &! + time_alpha_2*localTzzNm2(ij,el) ) &
              - We*convective_contrib_zz(ij,el) &
              - temp_giesekus_const*(time_beta_0*localTzz(ij,el)**2 + time_beta_1*localTzzNm1(ij,el)**2)


! Using LAPACK instead:
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

!             anorm=0d0
!            do ii=1,4
!               temp1=temp_matrix(ii,1)
!               do jj=2,4
!                 if (temp_matrix(ii,jj).gt.temp1) temp1=temp_matrix(ii,jj)
!              enddo
!               anorm=anorm+temp1
!             enddo

            call dgetrf( 4, 4, temp_matrix, 4, ipiv, info )
            IF (info.ne.0) THEN
              write(*,*) 'Error in calcStress_weakform:',el,info
              STOP
            ENDIF


!            call dgecon( '1', 4, temp_matrix, 4, anorm, rcond, work, iwork, info )
!            IF (info.ne.0) THEN
!           write(*,*) 'Error in calcStress_weakform:',el,info
!           STOP
!           ENDIF
!           IF (el.eq.2.and.ij.eq.273) print*,'rcond: ',rcond

            call dgetrs( 'N', 4, 1, temp_matrix, 4, ipiv, temp_rhs, 4, info )  
            IF (info.ne.0) THEN
              write(*,*) 'Error in calcStress_weakform:',el,info
              STOP
            ENDIF

! Calculate stress components from Inverse*RHS
!            tempTxxNext(ij,el) = idetA*(i11*r1 + i12*r2 + i13*r3)
!            tempTxyNext(ij,el) = idetA*(i21*r1 + i22*r2 + i23*r3)
!            tempTyyNext(ij,el) = idetA*(i31*r1 + i32*r2 + i33*r3)
!            tempTzzNext(ij,el) = r4/a44  

! LAPACK instead:
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
  END SUBROUTINE calcStress_weakform


  SUBROUTINE integrate_divergence_of_local_stress
! New version where we transfer the derivative onto the test function:
    IMPLICIT NONE
    INTEGER :: el,i,j,k,l,kl,il,ki,edge,ij
    DOUBLE PRECISION :: temp1,temp2,temp3,temp4,normx,normy
    DOUBLE PRECISION, DIMENSION(nptot) :: tempglob
    DOUBLE PRECISION, DIMENSION(0:(N+1)**2-1) :: templocal_INxx,templocal_INxy,templocal_INyy,templocal_INzz,&
             templocal_OUTx,templocal_OUTy
! Integrate the final solutions in each element
    stress_cont_to_stokes_x=0d0
    stress_cont_to_stokes_y=0d0

    IF (coordflag.eq.0) THEN
! CARTESIAN CASE !      
      DO el=1,numelm
        templocal_OUTx=0d0
        templocal_OUTy=0d0
        DO l=0,N
          DO k=0,N
            kl=k+l*(N+1)
            temp1=0d0
            temp2=0d0
            temp3=0d0
            temp4=0d0
            DO i=0,N
              ki=k+i*(N+1)
              il=i+l*(N+1)
 
! x-comp sums:
! used i as dummy variable, workings may differ.
              temp1 = temp1 + w(i)*d(i,k)*(localTxx(il,el)*dyde(i,l,el) - localTxy(il,el)*dxde(i,l,el))
              
              temp2 = temp2 + w(i)*d(i,l)*(localTxy(ki,el)*dxdp(k,i,el) - localTxx(ki,el)*dydp(k,i,el))
            
! y-comp sums:
              temp3 = temp3 + w(i)*d(i,k)*(localTxy(il,el)*dyde(i,l,el) - localTyy(il,el)*dxde(i,l,el))
      
              temp4 = temp4 + w(i)*d(i,l)*(localTyy(ki,el)*dxdp(k,i,el) - localTxy(ki,el)*dydp(k,i,el))
      
            ENDDO
!              IF (.not.bdflag(1,k))
            templocal_OUTx(kl) = w(l)*temp1 + w(k)*temp2
!             IF (.not.bdflag(2,k))
            templocal_OUTy(kl) = w(l)*temp3 + w(k)*temp4

          ENDDO
        ENDDO

        tempglob=0d0
        CALL vecglobalprolongation(templocal_OUTx,el,tempglob)
        stress_cont_to_stokes_x = stress_cont_to_stokes_x - tempglob
        tempglob=0d0
        CALL vecglobalprolongation(templocal_OUTy,el,tempglob)
        stress_cont_to_stokes_y = stress_cont_to_stokes_y - tempglob
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
          DO k=0,N
            kl=k+l*(N+1)

            temp1=0d0
            temp2=0d0
            temp3=0d0
            temp4=0d0

            DO i=0,N
              ki=k+i*(N+1)
              il=i+l*(N+1)
      
! z-comp sums:
! used i as dummy variable, workings may differ.
              temp1 = temp1 + nodeCoord(mapg(il,el),2)*w(i)*d(i,k)*(localTxx(il,el)*dyde(i,l,el) - localTxy(il,el)*dxde(i,l,el))
      
              temp2 = temp2 + nodeCoord(mapg(ki,el),2)*w(i)*d(i,l)*(localTxy(ki,el)*dxdp(k,i,el) - localTxx(ki,el)*dydp(k,i,el))
      
! r-comp sums:
              temp3 = temp3 + nodeCoord(mapg(il,el),2)*w(i)*d(i,k)*(localTxy(il,el)*dyde(i,l,el) - localTyy(il,el)*dxde(i,l,el))
      
              temp4 = temp4 + nodeCoord(mapg(ki,el),2)*w(i)*d(i,l)*(localTyy(ki,el)*dxdp(k,i,el) - localTxy(ki,el)*dydp(k,i,el))
            ENDDO
! x(z)-comp case:
!             IF (.not.bdflag(1,k))
            templocal_OUTx(kl) = w(l)*temp1 + w(k)*temp2
! Add the tau_theta,theta contribution to divergence in y(r)-comp case:
!              IF (.not.bdflag(2,k))
            templocal_OUTy(kl) = localTzz(kl,el)*jac(k,l,el)*w(k)*w(l) + w(l)*temp3 + w(k)*temp4
          ENDDO
        ENDDO

        tempglob=0d0
        CALL vecglobalprolongation(templocal_OUTx,el,tempglob)
        stress_cont_to_stokes_x = stress_cont_to_stokes_x - tempglob
        tempglob=0d0
        CALL vecglobalprolongation(templocal_OUTy,el,tempglob)
        stress_cont_to_stokes_y = stress_cont_to_stokes_y - tempglob
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
!         DO i=0,N
!           ij=local_edge_node(i,edge)
!           k=mapg(ij,el)
!           normx=norm_to_edge_node(1,i,edge,el)
!           normy=norm_to_edge_node(2,i,edge,el)
! !         IF (.not.bdflag(1,k))
!           stress_cont_to_stokes_x(k) = stress_cont_to_stokes_x(k) + ( localTxx(ij,el)*normx + localTxy(ij,el)*normy )*jac_on_edge(edge,el)*w(i)
! !         IF (.not.bdflag(2,k))
!           stress_cont_to_stokes_y(k) = stress_cont_to_stokes_y(k) + ( localTxy(ij,el)*normx + localTyy(ij,el)*normy )*jac_on_edge(edge,el)*w(i)
!         ENDDO
!       ENDIF
!     ENDDO  
  END SUBROUTINE integrate_divergence_of_local_stress
  
END MODULE viscoelastic_module
