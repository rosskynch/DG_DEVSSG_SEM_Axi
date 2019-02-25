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

MODULE SEM_module
  USE constants
  USE shared_data
  USE functions_module ! This may need to be replaced by the module which constants evalh( ) for now it is still functions_module.
  

!
! This module contains everything for the Spectral Element Method up to the construction of the Linear system for the Stokes problem.
!
  IMPLICIT NONE
  CONTAINS
  
  SUBROUTINE initialiseSEM
! Calls functions in SEM that will only be needed once
! In general, these are those which are not dependent upon the mesh or time.
! Anything which is mesh-dependent is dealt with in updateSEM.
    IMPLICIT NONE
    INTEGER :: i
!Mesh independent:
    CALL calcGLpoints
    CALL calcL
    CALL calcL1
    CALL calc_d
    CALL GLweights
    CALL create_evalh_array

  END SUBROUTINE initialiseSEM
  
  SUBROUTINE initialiseSEMarrays
! Calls functions in SEM that are mesh dependent.
! In general these are the values within SEM which change over time, due to mesh change.
! It is assumed that the geometry has already been updated when using this routine.
! Note boundary conditions must be applied via routine in the boundary_module.
    IMPLICIT NONE
    INTEGER :: el,ij,kl
    
    Retime_constant1 = Re*time_gamma_0/deltat
    Retime_constant2 = Re/deltat
    Wetime_constant1 = We*time_gamma_0/deltat
    Wetime_constant2 = We/deltat
    
    storeA_x=0d0
    storeA_y=0d0
    storeB_x=0d0
    storeB_y=0d0
!     storeC_x=0d0
!     storeC_y=0d0
!     storeCb=0d0
    storef_x=0d0
    storef_y=0d0
    storeg_p=0d0
    storeM_pressure=0d0
!     storeM_stress=0d0
    storeMv_x=0d0
    storeMv_y=0d0
    storeZ_p=0d0
    diff_x=0d0
    diff_y=0d0
    

    DO el=1,numelm
    
      CALL construct_diff_arrays(el)
    
      CALL constructA_x(el)
      CALL constructA_y(el)
    
      CALL constructB_x(el)
      CALL constructB_y(el)
    
      CALL constructpressuremassmatrix(el)
      CALL constructvelocitymassmatrix(el)
      CALL construct_zero_pressure_matrix(el)
    
!       CALL constructC_x(el)
!       CALL constructC_y(el)
!       CALL constructCb(el)
!       CALL constructStressMassMatrix(el)
    
    ENDDO
    CALL constructF

!!! REMOVED
!     IF (param_delta_a.gt.0d0) THEN
!       DO el=1,numelm
!   DO ij=0,NP1SQM1
!   DO kl=0,NP1SQM1
!     A_x(ij,kl,el) = param_beta_a(el)*storeA_x(ij,kl,el)
!     A_y(ij,kl,el) = param_beta_a(el)*storeA_y(ij,kl,el)
!   ENDDO
!   ENDDO
!       ENDDO
!     ELSEIF (param_delta_a.eq.0d0) THEN
    IF (param_beta_s.gt.0d0) THEN
      A_x = param_beta_s*storeA_x
      A_y = param_beta_s*storeA_y
    ELSE
      A_x = param_beta*storeA_x
      A_y = param_beta*storeA_y
    ENDIF

    B_x = storeB_x
    B_y = storeB_y
    Z_p = storeZ_p
    Mv_x = storeMv_x
    Mv_y = storeMv_y
    M_pressure = storeM_pressure
!     M_stress = storeM_stress
    f_x = storef_x
    f_y = storef_y
    g = storeg_p
!     C_x = storeC_x
!     C_y = storeC_y
!     Cb = storeCb

    IF (Re.gt.1d-5) THEN
      DO el=1,numelm
        DO ij=0,NP1SQM1
          A_x(ij,ij,el) = A_x(ij,ij,el) + Retime_constant1*Mv_x(ij,ij,el)
          A_y(ij,ij,el) = A_y(ij,ij,el) + Retime_constant1*Mv_y(ij,ij,el)
        ENDDO
      ENDDO
    ENDIF

!     DO el=1,numelm
!     DO ij=1,NM1SQ
!     write(*,*) (M_pressure(ij,kl,el),kl=1,NM1SQ)
!     ENDDO
!     ENDDO
!     stop

  END SUBROUTINE initialiseSEMarrays

  SUBROUTINE updateSEM
! Calls functions in SEM that are mesh dependent.
! In general these are the values within SEM which change over time, due to mesh change.
! It is assumed that the geometry has already been updated when using this routine.
! Note boundary conditions must be applied via routine in the boundary_module.
    IMPLICIT NONE
    INTEGER :: el,ij,kl

    IF (movingmeshflag.eq.1) THEN
      DO el=1,numelm
        IF (accordianflag(el)) THEN
          CALL construct_diff_arrays(el)
          CALL constructA_x(el)
          CALL constructA_y(el)    
          CALL constructB_x(el)
          CALL constructB_y(el)
          CALL constructpressuremassmatrix(el)
          CALL constructvelocitymassmatrix(el)
          CALL construct_zero_pressure_matrix(el)
!           CALL constructC_x(el)
!           CALL constructC_y(el)
!           CALL constructCb(el)
!           CALL constructStressMassMatrix(el)
        ENDIF
      ENDDO
      CALL constructF
    ENDIF

!     IF (param_delta_a.gt.0d0) THEN
!       DO el=1,numelm
!   DO ij=0,NP1SQM1
!   DO kl=0,NP1SQM1
!     A_x(ij,kl,el) = param_beta_a(el)*storeA_x(ij,kl,el)
!     A_y(ij,kl,el) = param_beta_a(el)*storeA_y(ij,kl,el)
!   ENDDO
!   ENDDO
!       ENDDO
!     ELSEIF (param_delta_a.eq.0d0) THEN
    IF (param_beta_s.gt.0d0) THEN
      A_x = param_beta_s*storeA_x
      A_y = param_beta_s*storeA_y
    ELSE
      A_x = param_beta*storeA_x
      A_y = param_beta*storeA_y
    ENDIF

    B_x = storeB_x
    B_y = storeB_y
    Z_p = storeZ_p
    Mv_x = storeMv_x
    Mv_y = storeMv_y
    M_pressure = storeM_pressure
!     M_stress = storeM_stress
    f_x = storef_x
    f_y = storef_y
    g = storeg_p
!     C_x = storeC_x
!     C_y = storeC_y
!     Cb = storeCb

    IF (Re.gt.1d-5) THEN
      DO el=1,numelm
        DO ij=0,NP1SQM1
          A_x(ij,ij,el) = A_x(ij,ij,el) + Retime_constant1*Mv_x(ij,ij,el)
          A_y(ij,ij,el) = A_y(ij,ij,el) + Retime_constant1*Mv_y(ij,ij,el)
        ENDDO
      ENDDO
    ENDIF
  END SUBROUTINE updateSEM
  
  SUBROUTINE add_local_vec_to_global_vec(loc_in,mult_by,glob_out)
    IMPLICIT NONE
    INTEGER :: el,ij
    DOUBLE PRECISION, INTENT(IN) :: loc_in(0:NP1SQM1,numelm)
    DOUBLE PRECISION, INTENT(INOUT) :: glob_out(1:nptot)
    DOUBLE PRECISION :: mult_by,temp_glob(1:nptot),temp_loc(0:NP1SQM1)
    
    DO el=1,numelm
      temp_loc=mult_by*loc_in(0:NP1SQM1,el)
      temp_glob = 0d0
      DO ij=0,NP1SQM1
        temp_glob(mapg(ij,el)) = temp_loc(ij)
      ENDDO
      glob_out=glob_out + temp_glob
    ENDDO
  END SUBROUTINE add_local_vec_to_global_vec
  
  
  SUBROUTINE calcGLpoints
!
!   computes the Gauss-Lobatto grid points from the roots of the degree N Lagrange
!   polynomial: (1-x^2)L' where L is the L' is the differential of the degree N Legendre polynomial. 
!
!   N:      degree of polynomial
!   alpha:  parameter in jacobi weight (set to zero in code)
!   beta:   parameter in jacobi weight (set to zero in code)
!   
!   gl:   output array with the Gauss-Lobatto roots
!           they are ordered from the smallest (-1.0) to the largest (1.0)

    IMPLICIT NONE
    INTEGER :: np,nh,i,j,jm,k,kstop

    DOUBLE PRECISION :: pnp1p,pdnp1p,pnp,pdnp, &
      pnm1p,pdnm1,pnp1m,pdnp1m,pnm,pnm1m, &
      pdnm,det,pnm1,cs,x,pnp1,pdnp1,pn,pdn, &
      rp,rm,ag,bg,dth,cd,sd,ss,poly,pder, &
      cssave,delx,epsg,recsum,alp,bet

    DOUBLE PRECISION, DIMENSION(0:N) :: hulpar
    alp=0d0
    bet=0d0
    kstop = 10
    epsg = 1.0d-25

    np = n+1
!
!   compute the parameters in the polynomial whose roots are desired
!
    CALL jacobf(alp,bet,np,pnp1p,pdnp1p,pnp,pdnp,pnm1p,pdnm1,1d0)
    CALL jacobf(alp,bet,np,pnp1m,pdnp1m,pnm,pdnm,pnm1m,pdnm1,-1d0)
    det = pnp*pnm1m-pnm*pnm1p
    rp = -pnp1p
    rm = -pnp1m
    ag = (rp*pnm1m-rm*pnm1p)/det
    bg = (rm*pnp-rp*pnm)/det
     
    gl(0) = 1d0
    nh = N
!
!   set-up recursion relation for initial guess for the roots
!
    dth = PI/(2d0*N+1d0)
    cd = dcos(2d0*dth)
    sd = dsin(2d0*dth)
    cs = dcos(dth)
    ss = dsin(dth)
!
!   compute the first half of the roots by polynomial deflation
!
    DO j=1,nh-1
      x=cs
      DO k=1,kstop
        CALL jacobf(alp,bet,np,pnp1,pdnp1,pn,pdn,pnm1,pdnm1,x)
        poly = pnp1+ag*pn+bg*pnm1
        pder = pdnp1+ag*pdn+bg*pdnm1
        recsum = 0d0
        jm = j-1
        DO i=0,jm
          recsum = recsum + 1d0/(x-gl(i))
        ENDDO
        delx = -poly/(pder-recsum*poly)
        x = x + delx
      ENDDO
      gl(j) = x
      cssave = cs*cd-ss*sd
      ss = cs*sd+ss*cd
      cs = cssave
    ENDDO
    gl(np-1) = -1d0
!
!   use symmetry for second half of the roots
!
    DO i=0,N
      hulpar(N-i) = gl(i)
    ENDDO
    DO i=0,N
      gl(i) =hulpar(i)
    ENDDO
    IF (N.eq.2*int(N/2d0)) gl(N/2)=0d0
  END SUBROUTINE calcGLpoints

  SUBROUTINE jacobf(alp,bet,m,poly,pder,polym1,pderm1,polym2,pderm2,x)
!
!   computes the jacobi polynomial (poly) and its derivative
!   (pder) of degree n at x

    INTEGER :: k,m
    DOUBLE PRECISION :: alp,bet,apb,x,&
      a1,a2,a3,a4,b3,  &
      poly,polylst,pder,&
      pderlst,polyn,pdern,psave,pdsave,&
      polym1,pderm1,polym2,pderm2
!      common /jacpar/alp,bet,rv

    apb = alp + bet
    poly = 1d0
    pder = 0d0
    IF (m.eq.0) RETURN
    polylst = poly
    pderlst = pder
    poly = 5d-1*(1d0+bet)*(x-1d0) + 5d-1*(1d0+alp)*(x+1d0)
    pder = 5d-1*(2d0+apb)
    IF (m.eq.1) RETURN
    DO k=2,m
      a1 = 2d0*k*(k+apb)*(2d0*k+apb-2d0)
      a2 = (2d0*k+apb-1d0)*(alp**2-bet**2)
      b3 = (2d0*k+apb-2d0)
      a3 = b3*(b3+1d0)*(b3+2d0)
      a4 = 2d0*(k+alp-1d0)*(k+bet-1d0)*(2d0*k+apb)
      polyn = ((a2+a3*x)*poly-a4*polylst)/a1
      pdern = ((a2+a3*x)*pder-a4*pderlst+a3*poly)/a1
      psave = polylst
      pdsave = pderlst
      polylst = poly
      poly = polyn
      pderlst = pder
      pder = pdern
    ENDDO
    polym1 = polylst
    pderm1 = pderlst
    polym2 = psave
    pderm2 = pdsave
  END SUBROUTINE jacobf
  
  SUBROUTINE calcL
! calculates values of legendre polynomials at gauss-lobatto points gl, using
! Bonnet’s recursion formula.
! outputs to leg.
  
    IMPLICIT NONE
    INTEGER :: i,j
    DOUBLE PRECISION, DIMENSION(0:N,0:N) :: temp
    
    DO i=0,N
      temp(i,0)=1d0
      temp(i,1)=gl(i)
    ENDDO
    
    DO j=1,NM1
      DO i=0,N
        temp(i,j+1) = ((2d0*dfloat(j)+1d0)*gl(i)*temp(i,j) - dfloat(j)*temp(i,j-1))/dfloat(j+1)
      ENDDO
    ENDDO
    
    DO i=0,N
      leg(i) = temp(i,N)
    ENDDO
  END SUBROUTINE calcL


  SUBROUTINE calcL1
    IMPLICIT NONE
    INTEGER :: i
    DOUBLE PRECISION :: F

    F=1d0
    IF (N.eq.2*int(N/2d0)) F = -1d0
    Leg1(N) = dfloat(N)*(dfloat(N)+1d0)/2d0
    Leg1(0) = F*Leg1(N)
    DO i=1,NM1
      Leg1(i) = 0d0
    ENDDO
  END SUBROUTINE calcL1


  SUBROUTINE GLweights
! calculates weight for Gauss-Lobatto integration
  
    IMPLICIT NONE
    INTEGER :: i,j,jj,ij
  
    DO i=0,N
      w(i) = 2d0/dfloat(N)/(dfloat(N)+1d0)/Leg(i)/Leg(i)
    ENDDO
    DO j=0,N
      jj=j*NP1
      DO i=0,N
        ij=i+jj
        weight_2d(ij)=w(i)*w(j)
      ENDDO
    ENDDO
    
  END SUBROUTINE GLweights


  SUBROUTINE calc_d

! code for differentials d(i,j) := h'_j(xi)
! xleg - GL points
! leg - legendre values at GL points
! d - ouput
    IMPLICIT NONE
    INTEGER :: i, j
    DOUBLE PRECISION :: const
  
    const = dfloat(N)*(dfloat(N)+1d0)/4d0
  
    DO i=0,N
      DO j=0,N
        IF (i.eq.j) THEN
          IF (i.eq.0) THEN
            d(i,j) = -const
          ELSE IF (i.eq.N) THEN
            d(i,j) = const
          ELSE
            d(i,j) = 0d0
          ENDIF
        ELSE
          d(i,j) = Leg(i)/(Leg(j)*(gl(i) - gl(j)))
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE calc_d
  

  
  SUBROUTINE create_evalh_array
    INTEGER :: i,j
    
    DO j=0,N
      DO i=1,NM1
        evalh(i,j) = calc_gl_htilde(i,j)
      ENDDO
    ENDDO
    
  END SUBROUTINE create_evalh_array
  
  SUBROUTINE construct_diff_arrays(el)
! Creates two arrays, diff_x and diff_y which represent differentiation of a test function(ij) in a local element w.r.t. x and y respectively.
!
! diff_x(pq,ij,el) represents the differential of velocity test function v_ij at the GLL point, pq in element el.
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: el
    INTEGER :: p,q,i,j,pq,ij,qq,jj
    DOUBLE PRECISION :: temp1,temp2
    
    DO j=0,N
      jj=j*NP1
      DO i=0,N
        ij = i + jj
        DO q=0,N
          qq=q*NP1
          DO p=0,N
            pq = p + qq
            temp1=0d0
            temp2=0d0
            IF (q.eq.j) temp1 = dyde(p,q,el)*d(p,i)
            IF (p.eq.i) temp2 = dydp(p,q,el)*d(q,j)
            diff_x(pq,ij,el) = (temp1-temp2)/jac(p,q,el)
            temp1=0d0
            temp2=0d0
            IF (p.eq.i) temp1 = dxdp(p,q,el)*d(q,j)
            IF (q.eq.j) temp2 = dxde(p,q,el)*d(p,i)
            diff_y(pq,ij,el) = (temp1-temp2)/jac(p,q,el)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE construct_diff_arrays


  SUBROUTINE constructF
  
    IMPLICIT NONE
    INTEGER :: i,j,jj,el,ij,globi

    DOUBLE PRECISION, DIMENSION(nptot) :: temp   ! NOTE: this may be the wrong counting style, might wanna change to 0:nptot-1
    DOUBLE PRECISION, DIMENSION(0:NP1SQM1) :: floc ! matrix for f on the local element
    
    storef_x=0d0
    IF (coordflag.eq.0) THEN
    DO el=1,numelm
      temp=0d0
      floc=0d0
      DO j=0,N
        jj=j*NP1
        DO i=0,N
          ij=i+jj
          globi=mapg(ij,el)
!     IF (.not.bdflag(1,globi)) THEN
          floc(ij)= fterm_x(globi)*w(i)*w(j)*jac(i,j,el)
!     ENDIF
        ENDDO
      ENDDO
      CALL vecglobalprolongation(floc,el,temp)
      storef_x = storef_x + temp
    ENDDO
! CYLINDRICAL (AXISYMMETRIC) CASE
    ELSEIF (coordflag.eq.1) THEN
    DO el=1,numelm
      temp=0d0
      floc=0d0
      DO j=0,N
        jj=j*NP1
        DO i=0,N
          ij=i+jj
          globi=mapg(ij,el)
!     IF (.not.bdflag(1,globi)) THEN
          floc(ij) = fterm_x(globi)*w(i)*w(j)*jac(i,j,el)*nodeCoord(globi,2)
!     ENDIF
        ENDDO
      ENDDO
      CALL vecglobalprolongation(floc,el,temp)
      storef_x = storef_x + temp
    ENDDO
    ENDIF

    storef_y=0d0
    IF (coordflag.eq.0) THEN
    DO el=1,numelm
      temp=0d0
      floc=0d0
      DO j=0,N
        jj=j*NP1
        DO i=0,N
          ij=i+jj
          globi=mapg(ij,el)
!     IF (.not.bdflag(2,globi)) THEN
          floc(ij)=fterm_y(globi)*w(i)*w(j)*jac(i,j,el)
!     ENDIF
        ENDDO
      ENDDO
      CALL vecglobalprolongation(floc,el,temp)
      storef_y = storef_y + temp
    ENDDO
! CYLINDRICAL (AXISYMMETRIC) CASE
    ELSEIF (coordflag.eq.1) THEN
    DO el=1,numelm
      temp=0d0
      floc=0d0
      DO j=0,N
        jj=j*NP1
        DO i=0,N
          ij=i+jj
          globi=mapg(ij,el)
!     IF (.not.bdflag(2,globi)) THEN
          floc(ij) = fterm_y(globi)*w(i)*w(j)*jac(i,j,el)*nodeCoord(globi,2)
!     ENDIF
        ENDDO
      ENDDO
      CALL vecglobalprolongation(floc,el,temp)
      storef_y = storef_y + temp
    ENDDO
    ENDIF
! For now, g, the pressure part of the f vector is zero - update if this ever changes
! (eg, for the integral pressure condition)

    storeg_p=0d0

  END SUBROUTINE constructF

  
  SUBROUTINE constructA_x(el)

! calculates stiffness matrix for all elements
! uses:
! d - the differentials d_j(x(i))
! w - the Gauss-Lobatto weights
! the differentials from the jacobian - dxdpsi etc)
  
    IMPLICIT NONE
  
    INTEGER, INTENT(IN) :: el
    INTEGER :: i,j,k,l,m,ij,kl
    DOUBLE PRECISION :: a
  
    IF (coordflag.eq.0) THEN 
! CARTESIAN CASE !
!     DO el=1,numelm
      DO l=0,N
        DO k=0,N
          kl=k+l*NP1
          DO j=0,N
            DO i=0,N
              ij=i+j*NP1
              IF (ij.ge.kl) THEN ! Ensure that the matrix is symmetric
                a = 0d0
                DO m=0,N 
                  IF (i.eq.k) THEN ! represents delta(k,i)
      
                    a = a + (w(k)*w(m)*d(m,j)*d(m,l)* &   !here switched n to m from my working as N is taken.
                      (dxdp(k,m,el)**2+dydp(k,m,el)**2))/jac(k,m,el)

                  ENDIF
                  IF (m.eq.l) THEN !represents delta(l,n) 
  
                    a = a - (w(i)*w(l)*d(l,j)*d(i,k)* & !here m is actually n
                            (dxdp(i,l,el)*dxde(i,l,el)+dydp(i,l,el)*dyde(i,l,el)))/jac(i,l,el)
                            
                  ENDIF
                  IF (j.eq.l) THEN !represents delta(l,j)
  
                    a = a + (w(m)*w(l)*d(m,i)*d(m,k)* &
                            (dxde(m,l,el)**2+dyde(m,l,el)**2))/jac(m,l,el)
                  ENDIF
                  IF (m.eq.k) THEN !represents delta(k,m)
  
                    a = a - (w(k)*w(j)*d(k,i)*d(j,l)* &
                            (dxdp(k,j,el)*dxde(k,j,el)+dydp(k,j,el)*dyde(k,j,el)))/jac(k,j,el)
                            
                  ENDIF
                ENDDO
                storeA_x(ij,kl,el) = a
                storeA_x(kl,ij,el) = a
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!     ENDDO  
    ELSEIF (coordflag.eq.1) THEN
! CYLINDERICAL POLAR (AXISYMMETRIC) CASE !
!
! NOTE:
! The x co-ordinate in the cartesian case is equivalent to the z component
! in the cylindrical polar coordinate case.
!
!     DO el=1,numelm
      DO l=0,N
        DO k=0,N
          kl=k+l*NP1
          DO j=0,N
            DO i=0,N
              ij=i+j*NP1
              IF (ij.ge.kl) THEN ! Ensure that the matrix is symmetric
                a = 0d0
                DO m=0,N 
              
                  IF (i.eq.k) THEN ! represents delta(k,i)

                    a = a + nodeCoord(mapg(k+m*NP1,el),2)*(w(k)*w(m)*d(m,j)*d(m,l)* &   !here switched n to m from my working as N is taken.
                      (dxdp(k,m,el)**2+dydp(k,m,el)**2))/jac(k,m,el)

                  ENDIF
                  IF (m.eq.l) THEN !represents delta(l,n) 

                    a = a - nodeCoord(mapg(i+l*NP1,el),2)*(w(i)*w(l)*d(l,j)*d(i,k)* & !here m is actually n
                      (dxdp(i,l,el)*dxde(i,l,el)+dydp(i,l,el)*dyde(i,l,el)))/jac(i,l,el)
                  ENDIF
                  IF (j.eq.l) THEN !represents delta(l,j)

                    a = a + nodeCoord(mapg(m+l*NP1,el),2)*(w(m)*w(l)*d(m,i)*d(m,k)* &
                      (dxde(m,l,el)**2+dyde(m,l,el)**2))/jac(m,l,el)
                  ENDIF
                  IF (m.eq.k) THEN !represents delta(k,m)

                    a = a - nodeCoord(mapg(k+j*NP1,el),2)*(w(k)*w(j)*d(k,i)*d(j,l)* &
                      (dxdp(k,j,el)*dxde(k,j,el)+dydp(k,j,el)*dyde(k,j,el)))/jac(k,j,el)
                  ENDIF
                ENDDO
                storeA_x(ij,kl,el) = a
                storeA_x(kl,ij,el) = a
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!     ENDDO
  ELSE
    print*, 'Error: No co-ordinate system specified...'
    print*, 'Stopping'
    STOP
  ENDIF
  
! OLD VERSION (Re/We):
!   A_x=param_beta*A_x
     
  END SUBROUTINE constructA_x
  
  SUBROUTINE constructA_y(el)

! calculates stiffness matrix for all elements
! uses:
! d - the differentials d_j(x(i))
! w - the Gauss-Lobatto weights
! the differentials from the jacobian - dxdpsi etc)
  
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: el
    INTEGER :: i,j,k,l,m,ij,kl
    DOUBLE PRECISION :: a,temp1,temp2,tempaxi
  
    IF (coordflag.eq.0) THEN 
! CARTESIAN CASE !
! IN THIS CASE A_x and A_y are the same, so no work is required
      storeA_y(0:NP1SQM1,0:NP1SQM1,el) = storeA_x(0:NP1SQM1,0:NP1SQM1,el)
    ELSEIF (coordflag.eq.1) THEN
! CYLINDERICAL POLAR (AXISYMMETRIC) CASE !
!
! NOTE:
! The y co-ordinate in the cartesian case is equivalent to the r component
! in the cylindrical polar coordinate case.
!
!     DO el=1,numelm
      DO l=0,N
        DO k=0,N
          kl=k+l*NP1
          DO j=0,N
            DO i=0,N
              ij=i+j*NP1
              IF (ij.ge.kl) THEN ! Ensure that the matrix is symmetric
                a = 0d0
                DO m=0,N 
                  IF (i.eq.k) THEN ! represents delta(k,i)

                    a = a + nodeCoord(mapg(k+m*NP1,el),2)*(w(k)*w(m)*d(m,j)*d(m,l)* &   !here switched n to m from my working as N is taken.
                      (dxdp(k,m,el)**2+dydp(k,m,el)**2))/jac(k,m,el)

                  ENDIF
                  IF (m.eq.l) THEN !represents delta(l,n) 

                    a = a - nodeCoord(mapg(i+l*NP1,el),2)*(w(i)*w(l)*d(l,j)*d(i,k)* & !here m is actually n
                      (dxdp(i,l,el)*dxde(i,l,el)+dydp(i,l,el)*dyde(i,l,el)))/jac(i,l,el)
                  ENDIF
                  IF (j.eq.l) THEN !represents delta(l,j)

                    a = a + nodeCoord(mapg(m+l*NP1,el),2)*(w(m)*w(l)*d(m,i)*d(m,k)* &
                      (dxde(m,l,el)**2+dyde(m,l,el)**2))/jac(m,l,el)
                  ENDIF
                  IF (m.eq.k) THEN !represents delta(k,m)

                    a = a - nodeCoord(mapg(k+j*NP1,el),2)*(w(k)*w(j)*d(k,i)*d(j,l)* &
                      (dxdp(k,j,el)*dxde(k,j,el)+dydp(k,j,el)*dyde(k,j,el)))/jac(k,j,el)
                  ENDIF
                ENDDO
! Additional terms arising in r-component
! Include IF STATEMENT to deal with the case that r=0, see workings for why we may ignore this term entirely (u_r = 0 when r=0)
     !nodeCoord(mapg(ij,el),2).gt.1d-6
                temp1=0d0
                temp2=0d0
                tempaxi=0d0
                IF ((i.eq.k).and.(j.eq.l)) THEN
                  IF (.not.wallsymmflag(mapg(ij,el))) THEN
                    tempaxi = tempaxi + jac(i,j,el)*w(i)*w(j)/nodeCoord(mapg(ij,el),2) 
                  ENDIF
                ENDIF
! Extra loop to generate terms in the quadrature sum s.t. r_pq=0
                IF (wallsymmflag(mapg(kl,el))) THEN
                  IF (k.eq.i) temp1 = temp1 + dxdp(k,l,el)*d(l,j)
                  IF (l.eq.j) temp1 = temp1 - dxde(k,l,el)*d(k,i)
                    temp1=temp1*w(k)*w(l)
                ENDIF
                IF (wallsymmflag(mapg(ij,el))) THEN
                  IF (i.eq.k) temp2 = temp2 + dxdp(i,j,el)*d(j,l)
                  IF (j.eq.l) temp2 = temp2 - dxde(i,j,el)*d(i,k)
                  temp2=temp2*w(i)*w(j)
                ENDIF
                a = a + tempaxi + temp1 + temp2
! End extra loop to generate the extra terms.
                storeA_y(ij,kl,el) = a
                storeA_y(kl,ij,el) = a
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!     ENDDO
    ELSE
      print*, 'Error: No co-ordinate system specified...'
      print*, 'Stopping'
      STOP
    ENDIF
  END SUBROUTINE constructA_y
    
  SUBROUTINE constructB_x(el)
! ! B_x(ij,kl,el) is the divergence operator acting on velocity space. The transpose gives the gradient operator acting on pressure space.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: el
    INTEGER :: i,j,k,l,ij,kl  
    DOUBLE PRECISION :: temp1,temp2

    IF (coordflag.eq.0) THEN 
! CARTESIAN CASE !
!     DO el=1,numelm

    DO l=0,N
      DO k=0,N
        kl=k+l*NP1
        DO j=1,NM1
          DO i=1,NM1
            ij=i+(j-1)*NM1
      
              temp1 = evalh(j,l)*w(l)*(dyde(i,l,el)*d(i,k)*w(i) + &
                evalh(i,0)*dyde(0,l,el)*d(0,k)*w(0) + &
                evalh(i,N)*dyde(N,l,el)*d(N,k)*w(N))
      
              temp2 = evalh(i,k)*w(k)*(dydp(k,j,el)*d(j,l)*w(j) + &
                evalh(j,0)*dydp(k,0,el)*d(0,l)*w(0) + &
                evalh(j,N)*dydp(k,N,el)*d(N,l)*w(N))
      
            storeB_x(ij,kl,el) = temp2-temp1
                
          ENDDO
        ENDDO
      ENDDO
    ENDDO 
!     ENDDO
    ELSEIF (coordflag.eq.1) THEN
! CYLINDERICAL POLAR (AXISYMMETRIC) CASE !
!
! NOTE:
! The x co-ordinate in the cartesian case is equivalent to the z component
! in the cylindrical polar coordinate case.
!
!     DO el=1,numelm
      DO l=0,N
        DO k=0,N
          kl=k+l*NP1
          DO j=1,NM1
            DO i=1,NM1
              ij=i+(j-1)*NM1
            
              temp1 = evalh(j,l)*w(l)*(nodeCoord(mapg(i+l*NP1,el),2)*&
                dyde(i,l,el)*d(i,k)*w(i) + &
                nodeCoord(mapg(l*NP1,el),2)*&
                evalh(i,0)*dyde(0,l,el)*d(0,k)*w(0) + &
                nodeCoord(mapg(N+l*NP1,el),2)*&
                evalh(i,N)*dyde(N,l,el)*d(N,k)*w(N))
          
              temp2 = evalh(i,k)*w(k)*(nodeCoord(mapg(k+j*NP1,el),2)*&
                dydp(k,j,el)*d(j,l)*w(j) + &
                nodeCoord(mapg(k,el),2)*&
                evalh(j,0)*dydp(k,0,el)*d(0,l)*w(0) + &
                nodeCoord(mapg(k+N*NP1,el),2)*&
                evalh(j,N)*dydp(k,N,el)*d(N,l)*w(N))
            
              storeB_x(ij,kl,el) = temp2-temp1
            ENDDO
          ENDDO
        ENDDO
      ENDDO 
!   ENDDO
    ELSE
      print*, 'Error: No co-ordinate system specified...'
      print*, 'Stopping'
      STOP
    ENDIF
  END SUBROUTINE constructB_x
  
  SUBROUTINE constructB_y(el)
! B_y(ij,kl,el) is the divergence operator acting on velocity space. The transpose gives the gradient operator acting on pressure space.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: el
    INTEGER :: i,j,k,l,ij,kl
    DOUBLE PRECISION :: temp1,temp2,temp3
    
    
    IF (coordflag.eq.0) THEN 
! CARTESIAN CASE !
!     DO el=1,numelm
      DO l=0,N
        DO k=0,N
          kl=k+l*NP1
          DO j=1,NM1
            DO i=1,NM1
              ij=i+(j-1)*NM1
              
              temp1 = evalh(i,k)*w(k)*(dxdp(k,j,el)*d(j,l)*w(j) + &
                evalh(j,0)*dxdp(k,0,el)*d(0,l)*w(0) + &    
                evalh(j,N)*dxdp(k,N,el)*d(N,l)*w(N))
          
              temp2 = evalh(j,l)*w(l)*(dxde(i,l,el)*d(i,k)*w(i) + &
                evalh(i,0)*dxde(0,l,el)*d(0,k)*w(0) + &
                evalh(i,N)*dxde(N,l,el)*d(N,k)*w(N))
                
              storeB_y(ij,kl,el) = temp2-temp1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!     ENDDO
    ELSEIF (coordflag.eq.1) THEN
! CYLINDERICAL POLAR (AXISYMMETRIC) CASE !
!
! NOTE:
! The y co-ordinate in the cartesian case is equivalent to the r component
! in the cylindrical polar coordinate case.
!
!     DO el=1,numelm
      DO l=0,N
        DO k=0,N
          kl=k+l*NP1
          DO j=1,NM1
            DO i=1,NM1
              ij=i+(j-1)*NM1
              temp1 = evalh(i,k)*w(k)*(nodeCoord(mapg(k+j*NP1,el),2)*& 
                dxdp(k,j,el)*d(j,l)*w(j) + &
                nodeCoord(mapg(k,el),2)*evalh(j,0)*&
                dxdp(k,0,el)*d(0,l)*w(0) + &
                nodeCoord(mapg(k+N*NP1,el),2)*evalh(j,N)*&
                dxdp(k,N,el)*d(N,l)*w(N))

              temp2 = evalh(j,l)*w(l)*(nodeCoord(mapg(i+l*NP1,el),2)*&
                dxde(i,l,el)*d(i,k)*w(i) + &
                nodeCoord(mapg(l*NP1,el),2)*evalh(i,0)*&
                dxde(0,l,el)*d(0,k)*w(0) + &
                nodeCoord(mapg(N+l*NP1,el),2)*evalh(i,N)*&
                dxde(N,l,el)*d(N,k)*w(N))

              temp3 = evalh(i,k)*evalh(j,l)*jac(k,l,el)*w(k)*w(l)

              storeB_y(ij,kl,el) = temp2-temp1-temp3
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!     ENDDO
    ELSE
      print*, 'Error: No co-ordinate system specified...'
      print*, 'Stopping'
      STOP
    ENDIF
  END SUBROUTINE constructB_y

  SUBROUTINE constructpressuremassmatrix(el)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: el
    INTEGER :: i,j,k,l,ij,kl,jj,ll!,p,q,pq,
    DOUBLE PRECISION :: temp

    IF (coordflag.eq.0) THEN 
! CARTESIAN CASE !
!     DO el=1,numelm
      DO l=1,NM1
        ll=(l-1)*NM1
        DO k=1,NM1
          kl=k+ll
          DO j=1,NM1
            jj=(j-1)*NM1
            DO i=1,NM1
              ij=i+jj
              IF (ij.ge.kl) THEN
! Add boundary terms which are always present in sum
                temp = evalh(k,0)*evalh(i,0)*w(0)*( evalh(l,0)*evalh(j,0)*jac(0,0,el)*w(0) + &
                  evalh(l,N)*evalh(j,N)*jac(0,N,el)*w(N) ) + &
                  evalh(k,N)*evalh(i,N)*w(N)*( evalh(l,0)*evalh(j,0)*jac(N,0,el)*w(0) + &
                  evalh(l,N)*evalh(j,N)*jac(N,N,el)*w(N) )
     
! Now include the terms depending upon kronecker-deltas...
                IF (k.eq.i) THEN
                  temp = temp + w(k)*( evalh(l,0)*evalh(j,0)*jac(k,0,el)*w(0) + &
                    evalh(l,N)*evalh(j,N)*jac(k,N,el)*w(N) )
                  IF (l.eq.j) THEN
                    temp = temp + w(k)*jac(k,l,el)*w(l)
                  ENDIF
                ENDIF
              IF (l.eq.j) THEN
                temp = temp + evalh(k,0)*evalh(i,0)*w(0)*jac(0,l,el)*w(l) + &
                  evalh(k,N)*evalh(i,N)*w(N)*jac(N,l,el)*w(l)
              ENDIF
              storeM_pressure(ij,kl,el) = temp
              storeM_pressure(kl,ij,el) = temp
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
!     ENDDO
    ELSEIF (coordflag.eq.1) THEN
! CYLINDERICAL POLAR (AXISYMMETRIC) CASE !
!     DO el=1,numelm
      DO l=1,NM1
        ll=(l-1)*NM1
        DO k=1,NM1
          kl=k+ll
          DO j=1,NM1
            jj=(j-1)*NM1
            DO i=1,NM1
              ij=i+jj
              temp=0d0
              IF (ij.ge.kl) THEN
! Add boundary terms which are always present in sum
                temp = evalh(k,0)*evalh(i,0)*w(0)*( evalh(l,0)*evalh(j,0)*jac(0,0,el)*w(0)*nodeCoord(mapg(0,el),2) + &
                  evalh(l,N)*evalh(j,N)*jac(0,N,el)*w(N)*nodeCoord(mapg(N*NP1,el),2) ) + &
                  evalh(k,N)*evalh(i,N)*w(N)*( evalh(l,0)*evalh(j,0)*jac(N,0,el)*w(0)*nodeCoord(mapg(N,el),2) + &
                  evalh(l,N)*evalh(j,N)*jac(N,N,el)*w(N)*nodeCoord(mapg(N+N*NP1,el),2) )

! Now include the terms depending upon kronecker-deltas...
                IF (k.eq.i) THEN
                  temp = temp + w(k)*( evalh(l,0)*evalh(j,0)*jac(k,0,el)*w(0)*nodeCoord(mapg(k,el),2) + &
                    evalh(l,N)*evalh(j,N)*jac(k,N,el)*w(N)*nodeCoord(mapg(k+N*NP1,el),2) )
                  IF (l.eq.j) THEN
                    temp = temp + w(k)*jac(k,l,el)*w(l)*nodeCoord(mapg(k+l*NP1,el),2)
                  ENDIF
                ENDIF
                IF (l.eq.j) THEN
                  temp = temp + evalh(k,0)*evalh(i,0)*w(0)*jac(0,l,el)*w(l)*nodeCoord(mapg(l*NP1,el),2) + &
                    evalh(k,N)*evalh(i,N)*w(N)*jac(N,l,el)*w(l)*nodeCoord(mapg(N+l*NP1,el),2)
                ENDIF
                storeM_pressure(ij,kl,el) = temp
                storeM_pressure(kl,ij,el) = temp
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!     ENDDO
    ELSE
      print*, 'Error: No co-ordinate system specified...'
      print*, 'Stopping'
      STOP
    ENDIF
  END SUBROUTINE constructpressuremassmatrix
  
  
  SUBROUTINE constructvelocitymassmatrix(el)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: el
    INTEGER :: i,j,ij
  
    IF (coordflag.eq.0) THEN 
!   DO el=1,numelm
      DO j=0,N
        DO i=0,N
          ij=i+j*NP1
          IF (.NOT.bdflag(1,mapg(ij,el))) THEN
            storeMv_x(ij,ij,el) = jac(i,j,el)*w(i)*w(j)
          ENDIF
          IF (.NOT.bdflag(2,mapg(ij,el))) THEN
            storeMv_y(ij,ij,el) = jac(i,j,el)*w(i)*w(j)
          ENDIF
        ENDDO
      ENDDO
!   ENDDO
    ELSEIF (coordflag.eq.1) THEN 
!   DO el=1,numelm
      DO j=0,N
        DO i=0,N
          ij=i+j*NP1
          IF (.NOT.bdflag(1,mapg(ij,el))) THEN
            storeMv_x(ij,ij,el) = jac(i,j,el)*w(i)*w(j)*nodeCoord(mapg(ij,el),2)
          ENDIF
          IF (.NOT.bdflag(2,mapg(ij,el))) THEN
            storeMv_y(ij,ij,el) = jac(i,j,el)*w(i)*w(j)*nodeCoord(mapg(ij,el),2)
          ENDIF
        ENDDO
      ENDDO
!   ENDDO
    ELSE
      print*, 'Error: No co-ordinate system specified...'
      print*, 'Stopping'
      STOP
    ENDIF
  END SUBROUTINE constructvelocitymassmatrix
  

! Redundant for now:
!   SUBROUTINE constructC_x(el)
!     IMPLICIT NONE
!     INTEGER, INTENT(IN) :: el
!     INTEGER :: i,j,k,l,ij,kl
!     DOUBLE PRECISION :: temp 
!     
!     IF (coordflag.eq.0) THEN 
! !       DO el=1,numelm
!   DO i=0,N
!     DO j=0,N
!       ij=i+j*NP1
!       DO k=0,N
!         DO l=0,N
!     kl=k+l*NP1
!     temp=0d0
! !     IF (j.eq.l) temp = dyde(i,j,el)*d(i,k)
! !     IF (i.eq.k) temp = temp - dydp(i,j,el)*d(j,l)
! !     storeC_x(kl,ij,el) = temp*w(i)*w(j)
! ! No integration by parts:
!     IF (.not.inflowflag(mapg(kl,el))) THEN
!       IF (j.eq.l) temp = dyde(k,l,el)*d(k,i)
!       IF (i.eq.k) temp = temp - dydp(k,l,el)*d(l,j)
!     ENDIF
!     storeC_x(kl,ij,el) = temp*w(k)*w(l)
!         ENDDO
!       ENDDO
!     ENDDO
!   ENDDO
! !       ENDDO
!     ELSEIF ( coordflag.eq.1 ) THEN
! !       DO el=1,numelm
!   DO i=0,N
!     DO j=0,N
!       ij=i+j*NP1
!       DO k=0,N
!         DO l=0,N
!     kl=k+l*NP1
!     temp=0d0
! !     IF (j.eq.l) temp = dyde(i,j,el)*d(i,k)
! !     IF (i.eq.k) temp = temp - dydp(i,j,el)*d(j,l)
! !     IF (abs(nodeCoord(mapg(ij,el),2)).gt.1d-5) THEN
! !       storeC_x(kl,ij,el) = temp*w(i)*w(j)*nodeCoord(mapg(ij,el),2)
! !     ELSE
! !       storeC_x(kl,ij,el) = temp*w(i)*w(j)
! !     ENDIF
! ! No integration by parts:
!     IF (.not.inflowflag(mapg(kl,el))) THEN
!       IF (j.eq.l) temp = dyde(k,l,el)*d(k,i)
!       IF (i.eq.k) temp = temp - dydp(k,l,el)*d(l,j)
!     ENDIF
! !       IF (wallsymmflag(mapg(kl,el))) THEN
! !         storeC_x(kl,ij,el) = temp*w(k)*w(l)
! !       ELSE
!         storeC_x(kl,ij,el) = temp*w(k)*w(l)!*nodeCoord(mapg(kl,el),2)
! !       ENDIF
!     
!         ENDDO
!       ENDDO
!     ENDDO
!   ENDDO
! !       ENDDO  
!     ENDIF
!   END SUBROUTINE constructC_x
!   
!   SUBROUTINE constructC_y(el)
!     IMPLICIT NONE
!     INTEGER, INTENT(IN) :: el
!     INTEGER :: i,j,k,l,ij,kl
!     DOUBLE PRECISION :: temp
!     
!   
!     IF (coordflag.eq.0) THEN 
! !       DO el=1,numelm
!   DO i=0,N
!     DO j=0,N
!       ij=i+j*NP1
!       DO k=0,N
!         DO l=0,N
!     kl=k+l*NP1
!     temp=0d0
! !     IF (i.eq.k) temp = dxdp(i,j,el)*d(j,l)
! !     IF (j.eq.l) temp = temp - dxde(i,j,el)*d(i,k)
! !     storeC_y(kl,ij,el) = temp*w(i)*w(j)
! ! No integration by parts:
!     IF(.not.inflowflag(mapg(kl,el))) THEN
!       IF (i.eq.k) temp = dxdp(k,l,el)*d(l,j)
!       IF (j.eq.l) temp = temp - dxde(k,l,el)*d(k,i)
!     ENDIF
!     storeC_y(kl,ij,el) = temp*w(k)*w(l)
!         ENDDO
!       ENDDO
!     ENDDO
!   ENDDO
! !       ENDDO
!     ELSEIF ( coordflag.eq.1 ) THEN
! !       DO el=1,numelm
!   DO i=0,N
!     DO j=0,N
!       ij=i+j*NP1
!       DO k=0,N
!         DO l=0,N
!     kl=k+l*NP1
!     temp=0d0
! !         IF (i.eq.k) temp = dxdp(i,j,el)*d(j,l)
! !     IF (j.eq.l) temp = temp - dxde(i,j,el)*d(i,k)
! !     IF (abs(nodeCoord(mapg(ij,el),2)).gt.1d-5) THEN
! !       storeC_y(kl,ij,el) = temp*w(i)*w(j)*nodeCoord(mapg(ij,el),2)
! !     ELSE
! !       storeC_y(kl,ij,el) = temp*w(i)*w(j)
! !     ENDIF
! ! No integration by parts:
!     IF(.not.inflowflag(mapg(kl,el))) THEN
!       IF (i.eq.k) temp = dxdp(k,l,el)*d(l,j)
!       IF (j.eq.l) temp = temp - dxde(k,l,el)*d(k,i)
!     ENDIF
! !       IF (wallsymmflag(mapg(kl,el))) THEN
! !         storeC_y(kl,ij,el) = temp*w(k)*w(l)
! !       ELSE
!         storeC_y(kl,ij,el) = temp*w(k)*w(l)!*nodeCoord(mapg(kl,el),2) 
! !       ENDIF
! 
!         ENDDO
!       ENDDO
!     ENDDO
!   ENDDO
! !       ENDDO  
!     ENDIF
!     
!   END SUBROUTINE constructC_y
  
!   SUBROUTINE constructCb(el)
!     IMPLICIT NONE
!     INTEGER, INTENT(IN) :: el
!     INTEGER :: i,edge,ij
! 
! ! Boundary integral for DG contribution:
! !     DO el=1,numelm
! !     IF (coordflag.eq.0) THEN 
!       DO edge=1,4
!   DO i=0,N
!     ij=local_edge_node(i,edge)
!     IF(.not.inflowflag(mapg(ij,el))) storeCb(i,i,edge,el) = jac_on_edge(edge,el)*w(i) 
!   ENDDO
!       ENDDO
! !     ELSEIF (coordflag.eq.0) THEN 
! !       DO edge=1,4
! !   DO i=0,N
! !     ij=local_edge_node(i,edge)
! !     IF(.not.inflowflag(mapg(ij,el))) storeCb(i,i,edge,el) = jac_on_edge(edge,el)*w(i)*nodeCoord(mapg(ij,el),2)
! !   ENDDO
! !       ENDDO
! !     ENDIF
! !     ENDDO
!       
!   END SUBROUTINE constructCb
  
!   SUBROUTINE constructStressMassMatrix(el)
!     IMPLICIT NONE
!     INTEGER, INTENT(IN) :: el
!     INTEGER :: i,j,ij
!     
!     IF (coordflag.eq.0) THEN 
! !       DO el=1,numelm
!   DO i=0,N
!     DO j=0,N
!       ij=i+j*NP1
!       IF (.NOT.inflowflag(mapg(ij,el))) THEN
!         storeM_stress(ij,el) = jac(i,j,el)*w(i)*w(j)
!       ENDIF
!     ENDDO
!   ENDDO
! !       ENDDO
!     ELSEIF ( coordflag.eq.1 ) THEN
! !       DO el=1,numelm
!   DO i=0,N
!     DO j=0,N
!       ij=i+j*NP1
!       IF (.NOT.inflowflag(mapg(ij,el))) THEN
! !         IF (wallsymmflag(mapg(ij,el))) THEN
! !     storeM_stress(ij,el) = jac(i,j,el)*w(i)*w(j)
! !         ELSE
!     storeM_stress(ij,el) = jac(i,j,el)*w(i)*w(j)!*nodeCoord(mapg(ij,el),2)
! !         ENDIF
!       ENDIF
!     ENDDO
!   ENDDO
! !       ENDDO  
!     ENDIF
!   END SUBROUTINE constructStressMassMatrix
  
  SUBROUTINE construct_zero_pressure_matrix(el)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: el
    INTEGER :: k,l,kl
    
    

!
! Z_p(kl,el) gives the integral of pressure test function q_k,l = h_k*h_l (tilde) over the element el.
! It is also multiplied by the zero integral of pressure constant, alpha (param_alphaZ)
!
!     storeZ_p=0d0
    IF (coordflag.eq.0) THEN
!       DO el=1,numelm
      DO l=1,NM1
        DO k=1,NM1
          kl=k+(l-1)*NM1
  
          storeZ_p(kl,el) = w(k)*( &
            jac(k,l,el)*w(l) + &
            evalh(l,0)*jac(k,0,el)*w(0) + &
            evalh(l,N)*jac(k,N,el)*w(N) ) + &
            evalh(k,0)*w(0)*( &
            jac(0,l,el)*w(l) + &
            evalh(l,0)*jac(0,0,el)*w(0) + &
            evalh(l,N)*jac(0,N,el)*w(N) ) + &
            evalh(k,N)*w(N)*( &
            jac(N,l,el)*w(l) + &
            evalh(l,0)*jac(N,0,el)*w(0) + &
            evalh(l,N)*jac(N,N,el)*w(N))
        ENDDO
      ENDDO
!       ENDDO
    ELSEIF (coordflag.eq.1) THEN 
!       DO el=1,numelm
      DO l=1,NM1
        DO k=1,NM1
          kl=k+(l-1)*NM1
    
          storeZ_p(kl,el) = w(k)*( &
            jac(k,l,el)*w(l)*nodeCoord(mapg(k+l*NP1,el),2) + &
            evalh(l,0)*jac(k,0,el)*w(0)*nodeCoord(mapg(k,el),2) + &
            evalh(l,N)*jac(k,N,el)*w(N)*nodeCoord(mapg(k+N*NP1,el),2) ) + &
            evalh(k,0)*w(0)*( &
            jac(0,l,el)*w(l)*nodeCoord(mapg(l*NP1,el),2) + &
            evalh(l,0)*jac(0,0,el)*w(0)*nodeCoord(mapg(0,el),2) + &
            evalh(l,N)*jac(0,N,el)*w(N)*nodeCoord(mapg(N*NP1,el),2) ) + &
            evalh(k,N)*w(N)*( &
            jac(N,l,el)*w(l)*nodeCoord(mapg(N+l*NP1,el),2) + &
            evalh(l,0)*jac(N,0,el)*w(0)*nodeCoord(mapg(N,el),2) + &
            evalh(l,N)*jac(N,N,el)*w(N)*nodeCoord(mapg(N+N*NP1,el),2) )
        ENDDO
      ENDDO
!       ENDDO
    ELSE
      print*, 'Error: No co-ordinate system specified...'
      print*, 'Stopping'
      STOP
    ENDIF
  END SUBROUTINE construct_zero_pressure_matrix


!   SUBROUTINE calc_gradient_of_U
!     IMPLICIT NONE
!     INTEGER :: el,i,j,p,ij,pj,ip,jN
!     DOUBLE PRECISION :: temp1,temp2,temp3,temp4
!     DOUBLE PRECISION, DIMENSION(0:NP1SQM1) :: localV_x,localV_y,localxx,localxy,localyx,localyy,localzz
!     DOUBLE PRECISION, DIMENSION(nptot) :: globalxx,globalxy,globalyx,globalyy,globalzz
! 
!     uxx=0d0
!     uxy=0d0
!     uyx=0d0
!     uyy=0d0
!     uzz=0d0
! ! calculate the derivatives of velocity at all node points:
! ! Cartesian Coords
!     IF (coordflag.eq.0) THEN
!       DO el=1,numelm
!       
! ! Restrict into local element
!   CALL veclocalrestrict(V_x,el,localV_x)
!   CALL veclocalrestrict(V_y,el,localV_y)
!   DO j=0,N
!     jN=j*NP1
!   DO i=0,N
!     ij=i+jN
!     temp1=0d0
!     temp2=0d0
!     temp3=0d0
!     temp4=0d0
! ! HAVE VERIFIED THIS WITH WORKINGS !
!     DO p=0,N
!       pj=p+jN
!       ip=i+p*NP1
!       temp1 = temp1 + localV_x(pj)*d(i,p)
!       temp2 = temp2 + localV_x(ip)*d(j,p) ! switched p and q to save a new loop
!       temp3 = temp3 + localV_y(pj)*d(i,p)
!       temp4 = temp4 + localV_y(ip)*d(j,p) ! switched p and q to save a new loop
!     ENDDO
!     
!     localxx(ij) = (dyde(i,j,el)*temp1 - dydp(i,j,el)*temp2)/jac(i,j,el)
!     localxy(ij) = (dxdp(i,j,el)*temp2 - dxde(i,j,el)*temp1)/jac(i,j,el)
!     localyx(ij) = (dyde(i,j,el)*temp3 - dydp(i,j,el)*temp4)/jac(i,j,el)
!     localyy(ij) = (dxdp(i,j,el)*temp4 - dxde(i,j,el)*temp3)/jac(i,j,el)
!   ENDDO
!   ENDDO
! ! Map back into global nodes
!   globalxx=0d0
!   CALL vecglobalprolongation(localxx,el,globalxx)
!   globalxy=0d0
!   CALL vecglobalprolongation(localxy,el,globalxy)
!   globalyx=0d0
!   CALL vecglobalprolongation(localyx,el,globalyx)  
!   globalyy=0d0
!   CALL vecglobalprolongation(localyy,el,globalyy)
! ! Split loop into edge and internal parts.
! ! internal parts will vectorise.
!   DO i=1,npedg
!     IF (upwinded_element(i).eq.el) THEN! May comment out if we wish to remove upwinding scheme.
!       uxx(i) = globalxx(i)
!       uxy(i) = globalxy(i)
!       uyx(i) = globalyx(i)
!       uyy(i) = globalyy(i)
!     ENDIF
!   ENDDO
!   DO i=npedg+1,nptot
!     uxx(i) = globalxx(i)
!     uxy(i) = globalxy(i)
!     uyx(i) = globalyx(i)
!     uyy(i) = globalyy(i)
!   ENDDO
! 
!       ENDDO
! ! Cylindrical Co-ords
!     ELSEIF (coordflag.eq.1) THEN
!       DO el=1,numelm
! ! Restrict into local element
!   CALL veclocalrestrict(V_x,el,localV_x)
!   CALL veclocalrestrict(V_y,el,localV_y)
!   
!   DO j=0,N
!   DO i=0,N
!     ij=i+j*NP1
!     temp1=0d0
!     temp2=0d0
!     temp3=0d0
!     temp4=0d0
! ! HAVE VERIFIED THIS WITH WORKINGS !
!     DO p=0,N
!       pj=p+j*NP1
!       ip=i+p*NP1
!       temp1 = temp1 + localV_x(mapg(pj,el))*d(i,p)
!       temp2 = temp2 + localV_x(mapg(ip,el))*d(j,p) ! switched p and q to save a new loop
!       temp3 = temp3 + localV_y(mapg(pj,el))*d(i,p)
!       temp4 = temp4 + localV_y(mapg(ip,el))*d(j,p) ! switched p and q to save a new loop
!     ENDDO
!     
!     localxx(ij) = (dyde(i,j,el)*temp1 - dydp(i,j,el)*temp2)/jac(i,j,el)
!     localxy(ij) = (dxdp(i,j,el)*temp2 - dxde(i,j,el)*temp1)/jac(i,j,el)
!     localyx(ij) = (dyde(i,j,el)*temp3 - dydp(i,j,el)*temp4)/jac(i,j,el)
!     localyy(ij) = (dxdp(i,j,el)*temp4 - dxde(i,j,el)*temp3)/jac(i,j,el)
!   
! ! Cylindrical Polar co-ordinate case has a theta,theta contribution even in axisymmetric:
!     IF (abs(nodeCoord(mapg(ij,el),2)).gt.1d-6) THEN
!       localzz(ij) = localV_y(ij)/nodeCoord(mapg(ij,el),2)
!     ELSE
!       localzz(ij) = 0d0
!     ENDIF
!   
!   ENDDO
!   ENDDO
!   
! ! Map back into global nodes
!   globalxx=0d0
!   CALL vecglobalprolongation(localxx,el,globalxx)
!   globalxy=0d0
!   CALL vecglobalprolongation(localxy,el,globalxy)
!   globalyx=0d0
!   CALL vecglobalprolongation(localyx,el,globalyx)  
!   globalyy=0d0
!   CALL vecglobalprolongation(localyy,el,globalyy)
!   globalzz=0d0
!   CALL vecglobalprolongation(localzz,el,globalzz)
! ! Split loop into edge and internal parts.
! ! internal parts will vectorise.
!   DO i=1,npedg
!     IF (upwinded_element(i).eq.el) THEN ! May comment out if we wish to remove upwinding scheme.
!       uxx(i) = globalxx(i)
!       uxy(i) = globalxy(i)
!       uyx(i) = globalyx(i)
!       uyy(i) = globalyy(i)
!       uzz(i) = globalzz(i)
!     ENDIF
!   ENDDO
!   DO i=npedg+1,nptot
!     uxx(i) = globalxx(i)
!     uxy(i) = globalxy(i)
!     uyx(i) = globalyx(i)
!     uyy(i) = globalyy(i)
!     uzz(i) = globalzz(i)
!   ENDDO
!   
!   
!       ENDDO
!     ENDIF
!     
!       
!   END SUBROUTINE calc_gradient_of_U
  
!   SUBROUTINE calc_gradient_of_U
!     IMPLICIT NONE
!     INTEGER :: el,i,j,p,ij,pj,ip
!     DOUBLE PRECISION :: temp1,temp2,temp3,temp4
!     
!     uxx=0d0
!     uxy=0d0
!     uyx=0d0
!     uyy=0d0
!     uzz=0d0
! ! calculate the derivatives of velocity at all node points:
! ! Cartesian Coords
!     IF (coordflag.eq.0) THEN
!       DO el=1,numelm
!   DO j=0,N
!   DO i=0,N
!     ij=i+j*NP1
!     temp1=0d0
!     temp2=0d0
!     temp3=0d0
!     temp4=0d0
! ! HAVE VERIFIED THIS WITH WORKINGS !
!     DO p=0,N
!       pj=p+j*NP1
!       ip=i+p*NP1
!       temp1 = temp1 + V_x(mapg(pj,el))*d(i,p)
!       temp2 = temp2 + V_x(mapg(ip,el))*d(j,p) ! switched p and q to save a new loop
!       temp3 = temp3 + V_y(mapg(pj,el))*d(i,p)
!       temp4 = temp4 + V_y(mapg(ip,el))*d(j,p) ! switched p and q to save a new loop
!     ENDDO
!     uxx(mapg(ij,el)) = uxx(mapg(ij,el)) + &
!           (dyde(i,j,el)*temp1 - dydp(i,j,el)*temp2)/jac(i,j,el)
!     uxy(mapg(ij,el)) = uxy(mapg(ij,el)) + &
!           (dxdp(i,j,el)*temp2 - dxde(i,j,el)*temp1)/jac(i,j,el)
!     uyx(mapg(ij,el)) = uyx(mapg(ij,el)) + &
!           (dyde(i,j,el)*temp3 - dydp(i,j,el)*temp4)/jac(i,j,el)
!     uyy(mapg(ij,el)) = uyy(mapg(ij,el)) + &
!           (dxdp(i,j,el)*temp4 - dxde(i,j,el)*temp3)/jac(i,j,el)
!   ENDDO
!   ENDDO
!       ENDDO
! ! Cylindrical Co-ords
!     ELSEIF (coordflag.eq.1) THEN
!       DO el=1,numelm
!   DO j=0,N
!   DO i=0,N
!     ij=i+j*NP1
!     temp1=0d0
!     temp2=0d0
!     temp3=0d0
!     temp4=0d0
! ! HAVE VERIFIED THIS WITH WORKINGS !
!     DO p=0,N
!       pj=p+j*NP1
!       ip=i+p*NP1
!       temp1 = temp1 + V_x(mapg(pj,el))*d(i,p)
!       temp2 = temp2 + V_x(mapg(ip,el))*d(j,p) ! switched p and q to save a new loop
!       temp3 = temp3 + V_y(mapg(pj,el))*d(i,p)
!       temp4 = temp4 + V_y(mapg(ip,el))*d(j,p) ! switched p and q to save a new loop
!     ENDDO
!     uxx(mapg(ij,el)) = uxx(mapg(ij,el)) + &
!           (dyde(i,j,el)*temp1 - dydp(i,j,el)*temp2)/jac(i,j,el)
!     uxy(mapg(ij,el)) = uxy(mapg(ij,el)) + &
!           (dxdp(i,j,el)*temp2 - dxde(i,j,el)*temp1)/jac(i,j,el)
!     uyx(mapg(ij,el)) = uyx(mapg(ij,el)) + &
!           (dyde(i,j,el)*temp3 - dydp(i,j,el)*temp4)/jac(i,j,el)
!     uyy(mapg(ij,el)) = uyy(mapg(ij,el)) + &
!           (dxdp(i,j,el)*temp4 - dxde(i,j,el)*temp3)/jac(i,j,el)
!   
! ! Cylindrical Polar co-ordinate case has a theta,theta contribution even in axisymmetric:
! !     IF (abs(nodeCoord(mapg(ij,el),2)).gt.1d-5) THEN
!     IF (.not.wallsymmflag(mapg(ij,el))) THEN
!       uzz(mapg(ij,el)) = uzz(mapg(ij,el)) + &
!          V_y(mapg(ij,el))/nodeCoord(mapg(ij,el),2)
!     ELSE
!       uzz(mapg(ij,el)) = uyy(mapg(ij,el))!0d0
!     ENDIF
!   
!   ENDDO
!   ENDDO
!       ENDDO
!     ENDIF
!     
!     
! ! Take into account doubling up on the edges of elements
! ! This takes equal contributions from all elements
! !
! ! May wish to change this to take the element in the direction of the streamline
! !
!     DO i=1,npedg
!       IF (mult(i).ne.1) THEN
!   uxx(i)=uxx(i)/dfloat(mult(i))
!   uxy(i)=uxy(i)/dfloat(mult(i))
!   uyx(i)=uyx(i)/dfloat(mult(i))
!   uyy(i)=uyy(i)/dfloat(mult(i))
!   uzz(i)=uzz(i)/dfloat(mult(i))
!       ENDIF
!     ENDDO
!       
!   END SUBROUTINE calc_gradient_of_U
  
  SUBROUTINE calc_local_gradient_of_U
! Calculates the gradient at each local point within an element.
! NOTE that the labelling convention is:
! GradUxy = dU_y/dx
! GradUyx = dU_x/dy, etc!
    IMPLICIT NONE
    INTEGER :: el,i,j,p,ij,pj,ip,jj
    DOUBLE PRECISION :: temp1,temp2,temp3,temp4
    
    localGradUxxNm1=localGradUxx
    localGradUyxNm1=localGradUyx
    localGradUxyNm1=localGradUxy
    localGradUyyNm1=localGradUyy
    localGradUzzNm1=localGradUzz    
    localGradUxx=0d0
    localGradUxy=0d0
    localGradUyx=0d0
    localGradUyy=0d0
    localGradUzz=0d0
! calculate the derivatives of velocity at all node points:
! Cartesian Coords
    IF (coordflag.eq.0) THEN
      DO el=1,numelm
        DO j=0,N
          jj=j*NP1
          DO i=0,N
            ij=i+jj
            temp1=0d0
            temp2=0d0
            temp3=0d0
            temp4=0d0
! HAVE VERIFIED THIS WITH WORKINGS !
            DO p=0,N
              pj=p+jj
              ip=i+p*NP1
!       temp1 = temp1 + V_x(mapg(pj,el))*d(i,p)
!       temp2 = temp2 + V_x(mapg(ip,el))*d(j,p) ! switched p and q to save a new loop
!       temp3 = temp3 + V_y(mapg(pj,el))*d(i,p)
!       temp4 = temp4 + V_y(mapg(ip,el))*d(j,p) ! switched p and q to save a new loop
              temp1 = temp1 + localV_x(pj,el)*d(i,p)
              temp2 = temp2 + localV_x(ip,el)*d(j,p) ! switched p and q to save a new loop
              temp3 = temp3 + localV_y(pj,el)*d(i,p)
              temp4 = temp4 + localV_y(ip,el)*d(j,p) ! switched p and q to save a new loop
            ENDDO
            localGradUxx(ij,el) = (dyde(i,j,el)*temp1 - dydp(i,j,el)*temp2)/jac(i,j,el)
            localGradUyx(ij,el) = (dxdp(i,j,el)*temp2 - dxde(i,j,el)*temp1)/jac(i,j,el)
            localGradUxy(ij,el) = (dyde(i,j,el)*temp3 - dydp(i,j,el)*temp4)/jac(i,j,el)
            localGradUyy(ij,el) = (dxdp(i,j,el)*temp4 - dxde(i,j,el)*temp3)/jac(i,j,el)
            IF (wallsymmflag(mapg(ij,el))) THEN
              localGradUyx(ij,el)=0d0
            ENDIF
          ENDDO
        ENDDO
      ENDDO
! Cylindrical Co-ords
    ELSEIF (coordflag.eq.1) THEN
      DO el=1,numelm
        DO j=0,N
          jj=j*NP1
          DO i=0,N
            ij=i+jj
            temp1=0d0
            temp2=0d0
            temp3=0d0
            temp4=0d0
! HAVE VERIFIED THIS WITH WORKINGS !
            DO p=0,N
              pj=p+jj
              ip=i+p*NP1
!       temp1 = temp1 + V_x(mapg(pj,el))*d(i,p)
!       temp2 = temp2 + V_x(mapg(ip,el))*d(j,p) ! switched p and q to save a new loop
!       temp3 = temp3 + V_y(mapg(pj,el))*d(i,p)
!       temp4 = temp4 + V_y(mapg(ip,el))*d(j,p) ! switched p and q to save a new loop
              temp1 = temp1 + localV_x(pj,el)*d(i,p)
              temp2 = temp2 + localV_x(ip,el)*d(j,p) ! switched p and q to save a new loop
              temp3 = temp3 + localV_y(pj,el)*d(i,p)
              temp4 = temp4 + localV_y(ip,el)*d(j,p) ! switched p and q to save a new loop
            ENDDO
            localGradUxx(ij,el) = (dyde(i,j,el)*temp1 - dydp(i,j,el)*temp2)/jac(i,j,el)
            localGradUyx(ij,el) = (dxdp(i,j,el)*temp2 - dxde(i,j,el)*temp1)/jac(i,j,el)
            localGradUxy(ij,el) = (dyde(i,j,el)*temp3 - dydp(i,j,el)*temp4)/jac(i,j,el)
            localGradUyy(ij,el) = (dxdp(i,j,el)*temp4 - dxde(i,j,el)*temp3)/jac(i,j,el)
! Cylindrical Polar co-ordinate case has a theta,theta contribution even in axisymmetric:
            IF (wallsymmflag(mapg(ij,el))) THEN ! in the r=0 case, L'Hopital's rules gives lim r->0 V_r/r = lim r->0 dV_r/dr
              localGradUzz(ij,el) = localGradUyy(ij,el)
!       localGradUyx(ij,el) = 0d0 ! shouldnt be setting this?!?!
            ELSE
!         IF (nodeCoord(mapg(ij,el),2).lt.1d-6) THEN
!     print*,'r=0 but wallsymmflag not true!',el,ij,mapg(ij,el)
!     STOP
!         ENDIF
              localGradUzz(ij,el) = localV_y(ij,el)/nodeCoord(mapg(ij,el),2)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF
!     DO i=1,npedg
!       IF (wallsymmflag(i)) THEN
!   DO el=1,numelm
!     ij=global_to_local_map(i,el)
!     IF (ij.ge.0d0) THEN
!       localGradUyx(ij,el)=0d0
!     ENDIF
!   ENDDO
!       ENDIF
!     ENDDO
  END SUBROUTINE calc_local_gradient_of_U
  
END MODULE SEM_module


