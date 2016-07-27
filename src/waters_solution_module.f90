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

MODULE waters_solution_module
  USE constants
  USE shared_data
  USE IFPORT
  IMPLICIT NONE
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: transient_cxx,transient_cxy,&
						  transient_txxNm1, transient_txxNm2, transient_txyNm1, transient_txyNm2

						  
  
  
  
  CONTAINS
  
  SUBROUTINE initialise_waters_solution
    IMPLICIT NONE
    INTEGER :: nterms,i
    DOUBLE PRECISION :: y_in(nptot),temp_const
    
    ALLOCATE(transient_u(nptot),transient_gradUyx(nptot)) 

    y_in = nodeCoord(1:nptot,2)
    
    IF (param_function_choice.eq.7) THEN ! 2-D Viscoelastic
      ALLOCATE( transient_txx(nptot), transient_txy(nptot), transient_cxx(nptot), transient_cxy(nptot))
      nterms=20
      CALL WatersIni(Re,We,nterms,(1d0-param_beta),param_beta,nptot,transient_cxx,transient_cxy,y_in)
      transient_u=0d0
      transient_txx=0d0
      transient_txy=0d0
    ELSEIF (param_function_choice.eq.8) THEN ! 3-D Viscoelastic.
      ALLOCATE(transient_txx(nptot), transient_txy(nptot), transient_cxx(nptot), transient_cxy(nptot))
      ALLOCATE(transient_txxNm1(nptot), transient_txyNm1(nptot), transient_txxNm2(nptot), transient_txyNm2(nptot))
      transient_txxNm2=0d0
      transient_txyNm2=0d0
      transient_txxNm1=0d0
      transient_txyNm1=0d0
      transient_u=0d0
      transient_txx=0d0
      transient_txy=0d0
! Approximates stress solution for 3-D Axisymmetric using a BDF2 treatment of the time derivative.
! Could use this, but it's zero anyway.
!       temp_const=(1d0+Wetime_constant1)
!       DO i=1,npedg
! 	transient_txy(i) = ((1d0-param_beta)*transient_gradUyx(i) +  Wetime_constant2*(time_alpha_0*transient_txyNm1(i) + time_alpha_1*transient_txyNm2(i)))/temp_const
! 	transient_txx(i) = (We*2d0*transient_gradUyx(i)*transient_txy(i) + Wetime_constant2*(time_alpha_0*transient_txxNm1(i) + time_alpha_1*transient_txxNm2(i)))/temp_const
!       ENDDO
    ENDIF
    
  END SUBROUTINE initialise_waters_solution
  
  
     
  SUBROUTINE calc_waters_solution!(y_in,t_in,u_out,txx_out,txy_out)
! Verified against paper by van Os & Phillips (2004)
! Produces the correct soln
! Requires the that initialise_waters_solution has been run beforehand (to calculate cxx/cxy)
    IMPLICIT NONE
    INTEGER :: hits,i,nterms
    DOUBLE PRECISION :: temp_const
    DOUBLE PRECISION, DIMENSION(nptot) :: y_in !, INTENT(IN)
!     DOUBLE PRECISION, DIMENSION(nptot), INTENT(OUT)  :: u_out,txx_out,txy_out

    nterms=20
    hits=nptot
    y_in = nodeCoord(1:nptot,2)
    
    IF (param_function_choice.eq.5) THEN ! 2-D Newtonian solution
      DO i=1,nptot
	transient_u(i) = newtonian_waters_V_x(i,timeN)
	transient_gradUyx(i) = newtonian_waters_GradV_yx(i,timeN)
      ENDDO
    ELSEIF (param_function_choice.eq.6) THEN ! 3-D Axisymmetric Newtonian
      DO i=1,nptot
	transient_u(i) = newtonian_waters_axisymmetric_V_x(i,timeN)
	transient_gradUyx(i) = newtonian_waters_axisymmetric_GradV_yx(i,timeN)
      ENDDO
      
    ELSEIF ( param_function_choice.eq.7 ) THEN 
! 2-D Viscoelastic solution, output v_x, txx,txy
! MUST HAVE CALLED INITALISE ROUTINE ABOVE ALREADY (done in initial conditions module)
      DO i=1,nptot
	transient_gradUyx(i) = waters_GradV_yx(i,timeN)
      ENDDO
      CALL Waters(Re, We, timeN, nterms, 1d0-param_beta, param_beta, nptot,&
		transient_cxx, transient_cxy, transient_u, transient_txx, transient_txy, y_in)

    ELSEIF (param_function_choice.eq.8) THEN 
! 3-D Axisymmetric Viscoelastic, note that we approximate the stress solution.
      DO i=1,nptot
	transient_u(i) = waters_axisymmetric_V_x(i,timeN)
	transient_gradUyx(i) = waters_axisymmetric_GradV_yx(i,timeN)
      ENDDO
  ! Approximates stress solution for 3-D Axisymmetric using a BDF2 treatment of the time derivative.
      transient_txxNm2=transient_txxNm1
      transient_txxNm1=transient_txx
      transient_txyNm2=transient_txyNm1
      transient_txyNm1=transient_txy      
      temp_const=(1d0+Wetime_constant1)
      DO i=1,nptot
	transient_txy(i) = ((1d0-param_beta)*transient_gradUyx(i) +  Wetime_constant2*(time_alpha_0*transient_txyNm1(i) + time_alpha_1*transient_txyNm2(i)))/temp_const
	transient_txx(i) = (We*2d0*transient_gradUyx(i)*transient_txy(i) + Wetime_constant2*(time_alpha_0*transient_txxNm1(i) + time_alpha_1*transient_txxNm2(i)))/temp_const
      ENDDO
    ELSE
      write(*,*) 'ERROR: called calc_waters_solution with inconsistent param_function_choice: ',param_function_choice
      STOP
    ENDIF


  END SUBROUTINE calc_waters_solution
  
  
! BEGIN FUNCTIONS FOR NEWTOWNIAN TIME-DEPENDENT POISEUILLE FLOW:

  DOUBLE PRECISION FUNCTION newtonian_waters_V_x(globi_in,t) RESULT(out)
! 2-D Time dependent newtonian POISEUILLE flow, from waters
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: globi_in
    INTEGER :: i,n_crit
    DOUBLE PRECISION, INTENT(IN) :: t
    DOUBLE PRECISION :: Ay,add_on_u,bigN,temp,crit,y

    y=nodeCoord(globi_in,2)
    
    crit=1d-6
    n_crit=0
    Ay = -4d0*y*(y - 1d0)
!     Ay = 0.25d0*(4d0 - y**2)
    
    i=0
    add_on_u = 0d0
    DO WHILE (n_crit.le.3)
      i=i+1
      bigN = (2d0*dfloat(i)-1d0)*pi
      temp = (sin(bigN*y)/bigN**3)*exp(-bigN**2*(t/Re))
      add_on_u = add_on_u + temp
      IF (abs(temp).lt.crit) THEN
	n_crit= n_crit + 1
      ELSE
	n_crit = 0
      ENDIF
      
    ENDDO
    
    
    out = Ay - 32d0*add_on_u
    
  END FUNCTION newtonian_waters_V_x
  
  
  DOUBLE PRECISION FUNCTION newtonian_waters_GradV_yx(globi_in,t) RESULT(out)
! 2-D Time dependent newtonian POISEUILLE flow, from waters
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: globi_in
    INTEGER :: i,n_crit
    DOUBLE PRECISION, INTENT(IN) :: t
    DOUBLE PRECISION :: Ayprime,add_on_u,bigN,temp,crit,y

    y=nodeCoord(globi_in,2)
    
    crit=1d-6
    n_crit=0
    Ayprime = 4d0-8d0*y
!     Ay = 0.25d0*(4d0 - y**2)
    
    i=0
    add_on_u = 0d0
    DO WHILE (n_crit.le.3)
      i=i+1
      bigN = (2d0*dfloat(i)-1d0)*pi
      temp = (cos(bigN*y)/bigN**2)*exp(-bigN**2*(t/Re))
      add_on_u = add_on_u + temp
      IF (abs(temp).lt.crit) THEN
	n_crit= n_crit + 1
      ELSE
	n_crit = 0
      ENDIF
      
    ENDDO
    
    
    out = Ayprime - 32d0*add_on_u
    
  END FUNCTION newtonian_waters_GradV_yx


! NOT IN USE:
! DOUBLE PRECISION FUNCTION calcWaters_analytical_U(y,t) RESULT(out_U)
!   
!     IMPLICIT NONE
!     INTEGER :: i,n_crit
!     DOUBLE PRECISION, INTENT(IN) :: y,t
!     DOUBLE PRECISION :: Ay,El,&
! 			add_on_u,temp,crit,&
! 			bigN,alphaN,betaN,temp_betaN,gammaN,alphaNStar,betaNstar,aN,bN,pN,qN,bigGN
!     
!     crit = 1d-6
!     n_crit=0
!     Ay = -4d0*y*(y - 1d0)
!     El = We/Re
!     
!     i=0
!     add_on_u = 0d0
! 
!     DO WHILE (n_crit.le.3)
!     
!       i=i+1
!       bigN = (2d0*dfloat(i) - 1d0)*pi
!       alphaN = 1d0 + param_beta*El*bigN**2
!       temp_betaN = alphaN**2 - 4d0*El*bigN**2
!       IF (temp_betaN.lt.0d0) THEN 
! 	betaN=sqrt(-temp_betaN)
!       ELSE ! betaN^2 >= 0 !
! 	betaN=sqrt(temp_betaN)
!       ENDIF
!       gammaN = 1d0 + bigN**2*El*(param_beta - 2d0)
!       alphaNStar = alphaN/2d0/El
!       betaNStar = betaN/2d0/El
!       aN = 1d0 + gammaN/betaN
!       bN = 1d0 - gammaN/betaN
!       pN = -alphaNStar + betaNStar
!       qN = -alphaNStar - betaNStar
! 
!       IF (temp_betaN.ge.0d0) THEN
! 	bigGN = 0.5d0*( aN*exp(pN*t) + bN*exp(qN*t) )
!       ELSE ! betaN^2 < 0 !
! 	bigGN = exp(-alphaNStar*t)*(cos(betaNStar*t)+gammaN/betaN*sin(betaNStar*t))
!       ENDIF
!       
!       temp = (sin(bigN*y)/bigN**3)*bigGN
!       add_on_u = add_on_u + temp
!       
!       IF (abs(temp).lt.crit) THEN
! 	n_crit= n_crit + 1
!       ELSE
! 	n_crit = 0
!       ENDIF
! 
!     ENDDO
!     
!     
!     out_U = Ay - 32d0*add_on_u
!     
!   END FUNCTION calcWaters_analytical_U
  
  DOUBLE PRECISION FUNCTION waters_GradV_yx(globi,t) RESULT(out_Uxy)
  
    IMPLICIT NONE
    INTEGER :: i,nterms
    INTEGER, INTENT(IN) :: globi
    DOUBLE PRECISION, INTENT(IN) :: t
    DOUBLE PRECISION :: Ay1,El,y,&
			add_on_uxy,temp,&
			bigN,alphaN,betaN,temp_betaN,gammaN,alphaNStar,betaNstar,aN,bN,pN,qN,bigGN    
!     INTEGER :: i,n_crit
!     DOUBLE PRECISION, INTENT(IN) ::y,t
!     DOUBLE PRECISION :: Ay1,El,&
! 			add_on_uxy,temp,crit,&
! 			bigN,alphaN,betaN,temp_betaN,gammaN,alphaNStar,betaNstar,aN,bN,pN,qN,bigGN
    
    
    y = nodeCoord(globi,2)    
    nterms=20
    Ay1= -8d0*y + 4d0
    El = We/Re
    
  
    add_on_uxy = 0d0
    
    DO i=1,nterms
    
      bigN = (2d0*dfloat(i) - 1d0)*pi
      alphaN = 1d0 + param_beta*El*bigN**2
      temp_betaN = alphaN**2 - 4d0*bigN**2*El
      IF (temp_betaN.lt.0d0) THEN 
	betaN=sqrt(-temp_betaN)
      ELSE ! betaN^2 >= 0 !
	betaN=sqrt(temp_betaN)
      ENDIF
      gammaN = 1d0 + bigN**2*El*(param_beta - 2d0)
      alphaNStar = alphaN/2d0/El
      betaNStar = betaN/2d0/El
      aN = 1d0 + gammaN/betaN
      bN = 1d0 - gammaN/betaN
      pN = -alphaNStar + betaNStar
      qN = -alphaNStar - betaNStar

      IF (temp_betaN.ge.0d0) THEN
	bigGN = 0.5d0*( aN*exp(pN*t) + bN*exp(qN*t) )
      ELSE ! betaN^2 < 0 !
	bigGN = exp(-alphaNStar*t)*(cos(betaNStar*t)+gammaN/betaN*sin(betaNStar*t))
      ENDIF
      
      temp = (cos(bigN*y)/bigN**2)*bigGN
      
      add_on_uxy = add_on_uxy + temp 
      
      
    ENDDO
    
    
    out_Uxy = Ay1 - 32d0*add_on_uxy
    
  END FUNCTION waters_GradV_yx
  
  DOUBLE PRECISION FUNCTION newtonian_waters_axisymmetric_V_x(globi,t) RESULT(out)

    IMPLICIT NONE
    INTEGER :: i,nterms
    INTEGER, INTENT(IN) :: globi
    DOUBLE PRECISION, INTENT(IN) :: t
    DOUBLE PRECISION :: Ay, add_on_u,temp,y,bigN
    DOUBLE PRECISION, DIMENSION(100) :: besj0_zeros


    besj0_zeros =(/ 2.4048255576957729, 5.5200781102863106, 8.6537279129110125, 11.791534439014281, 14.930917708487787,&
		    18.071063967910924, 21.211636629879258, 24.352471530749302, 27.493479132040253, 30.634606468431976,&
		    33.775820213573567, 36.917098353664045, 40.058425764628240, 43.199791713176730, 46.341188371661815,&
		    49.482609897397815, 52.624051841114998, 55.765510755019982, 58.906983926080940, 62.048469190227166,&
		    65.189964800206866, 68.331469329856802, 71.472981603593738, 74.614500643701831, 77.756025630388052,& 
		    80.897555871137627, 84.039090776938195, 87.180629843641157, 90.322172637210485, 93.463718781944777,&
		    96.605267950996264, 99.746819858680595, 102.88837425419480, 106.02993091645162, 109.17148964980538,&
		    112.31305028049491, 115.45461265366694, 118.59617663087253, 121.73774208795096, 124.87930891323295,&
		    128.02087700600833, 131.16244627521391, 134.30401663830546, 137.44558802028428, 140.58716035285428,&
		    143.72873357368974, 146.87030762579664, 150.01188245695477, 153.15345801922788, 156.29503426853353,&
		    159.43661116426316, 162.57818866894667, 165.71976674795502, 168.86134536923583, 172.00292450307819,&
		    175.14450412190274, 178.28608420007376, 181.42766471373105, 184.56924564063871, 187.71082696004936,&
		    190.85240865258152, 193.99399070010912, 197.13557308566141, 200.27715579333241, 203.41873880819864,&
		    206.56032211624446, 209.70190570429406, 212.84348955994949, 215.98507367153402, 219.12665802804057,&
		    222.26824261908430, 225.40982743485932, 228.55141246609881, 231.69299770403853, 234.83458314038324,&
		    237.97616876727565, 241.11775457726802, 244.25934056329569, 247.40092671865281, 250.54251303696995,&
		    253.68409951219309, 256.82568613856444, 259.96727291060449, 263.10885982309549, 266.25044687106589,&
		    269.39203404977604, 272.53362135470491, 275.67520878153744, 278.81679632615311, 281.95838398461490,&
		    285.09997175315954, 288.24155962818770, 291.38314760625519, 294.52473568406498, 297.66632385845895,&
		    300.80791212641111, 303.94950048502056, 307.09108893150506, 310.23267746319499, 313.37426607752786 /)
    
    y = nodeCoord(globi,2)!/rad_cylinder
    nterms=20
    Ay= 1d0-y**2
    
    add_on_u = 0d0
    DO i=1,nterms

      bigN = besj0_zeros(i)
      temp = (DBESJ0(bigN*y)/DBESJ1(bigN)/bigN**3)*exp(-(t/Re)*bigN**2)
      add_on_u = add_on_u + temp

    ENDDO
    
    out = Ay - 8d0*add_on_u
    
  END FUNCTION newtonian_waters_axisymmetric_V_x


  DOUBLE PRECISION FUNCTION newtonian_waters_axisymmetric_GradV_yx(globi,t) RESULT(out)

    IMPLICIT NONE
    INTEGER :: i,nterms
    INTEGER, INTENT(IN) :: globi
    DOUBLE PRECISION, INTENT(IN) :: t
    DOUBLE PRECISION :: Ayprime, add_on_u,temp,y,bigN
    DOUBLE PRECISION, DIMENSION(100) :: besj0_zeros


    besj0_zeros =(/ 2.4048255576957729, 5.5200781102863106, 8.6537279129110125, 11.791534439014281, 14.930917708487787,&
		    18.071063967910924, 21.211636629879258, 24.352471530749302, 27.493479132040253, 30.634606468431976,&
		    33.775820213573567, 36.917098353664045, 40.058425764628240, 43.199791713176730, 46.341188371661815,&
		    49.482609897397815, 52.624051841114998, 55.765510755019982, 58.906983926080940, 62.048469190227166,&
		    65.189964800206866, 68.331469329856802, 71.472981603593738, 74.614500643701831, 77.756025630388052,& 
		    80.897555871137627, 84.039090776938195, 87.180629843641157, 90.322172637210485, 93.463718781944777,&
		    96.605267950996264, 99.746819858680595, 102.88837425419480, 106.02993091645162, 109.17148964980538,&
		    112.31305028049491, 115.45461265366694, 118.59617663087253, 121.73774208795096, 124.87930891323295,&
		    128.02087700600833, 131.16244627521391, 134.30401663830546, 137.44558802028428, 140.58716035285428,&
		    143.72873357368974, 146.87030762579664, 150.01188245695477, 153.15345801922788, 156.29503426853353,&
		    159.43661116426316, 162.57818866894667, 165.71976674795502, 168.86134536923583, 172.00292450307819,&
		    175.14450412190274, 178.28608420007376, 181.42766471373105, 184.56924564063871, 187.71082696004936,&
		    190.85240865258152, 193.99399070010912, 197.13557308566141, 200.27715579333241, 203.41873880819864,&
		    206.56032211624446, 209.70190570429406, 212.84348955994949, 215.98507367153402, 219.12665802804057,&
		    222.26824261908430, 225.40982743485932, 228.55141246609881, 231.69299770403853, 234.83458314038324,&
		    237.97616876727565, 241.11775457726802, 244.25934056329569, 247.40092671865281, 250.54251303696995,&
		    253.68409951219309, 256.82568613856444, 259.96727291060449, 263.10885982309549, 266.25044687106589,&
		    269.39203404977604, 272.53362135470491, 275.67520878153744, 278.81679632615311, 281.95838398461490,&
		    285.09997175315954, 288.24155962818770, 291.38314760625519, 294.52473568406498, 297.66632385845895,&
		    300.80791212641111, 303.94950048502056, 307.09108893150506, 310.23267746319499, 313.37426607752786 /)
    
    y = nodeCoord(globi,2)
    nterms=20
    Ayprime = -2d0*y
    
    add_on_u = 0d0
    DO i=1,nterms

      bigN = besj0_zeros(i)
      temp = (DBESJ1(bigN*y)/DBESJ1(bigN)/bigN**2)*exp(-(t/Re)*bigN**2)
      add_on_u = add_on_u + temp
      

    ENDDO
    
    out = Ayprime + 8d0*add_on_u
    
  END FUNCTION newtonian_waters_axisymmetric_GradV_yx


  DOUBLE PRECISION FUNCTION waters_axisymmetric_V_x(globi,t) RESULT(out)
  
    IMPLICIT NONE
    INTEGER :: i,nterms
    INTEGER, INTENT(IN) :: globi
    DOUBLE PRECISION, INTENT(IN) :: t
    DOUBLE PRECISION :: Ay,El,y,&
			add_on_u,temp,crit,&
			bigN,alphaN,betaN,temp_betaN,gammaN,alphaNStar,betaNstar,aN,bN,pN,qN,bigGN
    DOUBLE PRECISION, DIMENSION(20) :: besj0_zeros


    besj0_zeros = (/ 2.4048255576957729, 5.5200781102863106, 8.6537279129110125, 11.791534439014281, 14.930917708487787, &
		    18.071063967910924, 21.211636629879258, 24.352471530749302, 27.493479132040253, 30.634606468431976, &
		    33.775820213573567, 36.917098353664045, 40.058425764628240, 43.199791713176730, 46.341188371661815, &
		    49.482609897397815, 52.624051841114998, 55.765510755019982, 58.906983926080940, 62.048469190227166 /)
    
    y = nodeCoord(globi,2)
    nterms=20
    Ay = 1d0-y*y
    El = We/Re
    
    i=0
    add_on_u = 0d0

    DO i=1,nterms
    

      bigN = besj0_zeros(i)
      alphaN = 1d0 + param_beta*El*bigN**2
      temp_betaN = alphaN**2 - 4d0*El*bigN**2
      IF (temp_betaN.lt.0d0) THEN 
	betaN=sqrt(-temp_betaN)
      ELSE ! betaN^2 >= 0 !
	betaN=sqrt(temp_betaN)
      ENDIF
      gammaN = 1d0 + bigN**2*El*(param_beta - 2d0)
      alphaNStar = alphaN/2d0/El
      betaNStar = betaN/2d0/El
      aN = 1d0 + gammaN/betaN
      bN = 1d0 - gammaN/betaN
      pN = -alphaNStar + betaNStar
      qN = -alphaNStar - betaNStar

      IF (temp_betaN.ge.0d0) THEN
	bigGN = 0.5d0*( aN*exp(pN*t) + bN*exp(qN*t) )
      ELSE ! betaN^2 < 0 !
	bigGN = exp(-alphaNStar*t)*(cos(betaNStar*t)+gammaN/betaN*sin(betaNStar*t))
      ENDIF
      
      temp = (DBESJ0(bigN*y)/DBESJ1(bigN)/bigN**3)*bigGN
      add_on_u = add_on_u + temp
      

    ENDDO
    
    
    out = Ay - 8d0*add_on_u
    
  END FUNCTION waters_axisymmetric_V_x
  
  
  DOUBLE PRECISION FUNCTION waters_axisymmetric_GradV_yx(globi,t) RESULT(out)
  
    IMPLICIT NONE
    INTEGER :: i,nterms
    INTEGER, INTENT(IN) :: globi
    DOUBLE PRECISION, INTENT(IN) :: t
    DOUBLE PRECISION :: Ayprime,El,y,&
			add_on_u,temp,crit,&
			bigN,alphaN,betaN,temp_betaN,gammaN,alphaNStar,betaNstar,aN,bN,pN,qN,bigGN
    DOUBLE PRECISION, DIMENSION(20) :: besj0_zeros


    besj0_zeros = (/ 2.4048255576957729, 5.5200781102863106, 8.6537279129110125, 11.791534439014281, 14.930917708487787, &
		    18.071063967910924, 21.211636629879258, 24.352471530749302, 27.493479132040253, 30.634606468431976, &
		    33.775820213573567, 36.917098353664045, 40.058425764628240, 43.199791713176730, 46.341188371661815, &
		    49.482609897397815, 52.624051841114998, 55.765510755019982, 58.906983926080940, 62.048469190227166 /)
    
    y = nodeCoord(globi,2)
    nterms=20
    Ayprime = -2d0*y
    El = We/Re
    
    i=0
    add_on_u = 0d0

    DO i=1,nterms
    

      bigN = besj0_zeros(i)
      alphaN = 1d0 + param_beta*El*bigN**2
      temp_betaN = alphaN**2 - 4d0*El*bigN**2
      IF (temp_betaN.lt.0d0) THEN 
	betaN=sqrt(-temp_betaN)
      ELSE ! betaN^2 >= 0 !
	betaN=sqrt(temp_betaN)
      ENDIF
      gammaN = 1d0 + bigN**2*El*(param_beta - 2d0)
      alphaNStar = alphaN/2d0/El
      betaNStar = betaN/2d0/El
      aN = 1d0 + gammaN/betaN
      bN = 1d0 - gammaN/betaN
      pN = -alphaNStar + betaNStar
      qN = -alphaNStar - betaNStar

      IF (temp_betaN.ge.0d0) THEN
	bigGN = 0.5d0*( aN*exp(pN*t) + bN*exp(qN*t) )
      ELSE ! betaN^2 < 0 !
	bigGN = exp(-alphaNStar*t)*(cos(betaNStar*t)+gammaN/betaN*sin(betaNStar*t))
      ENDIF
      
      temp = (DBESJ1(bigN*y)/DBESJ1(bigN)/bigN**2)*bigGN
      add_on_u = add_on_u + temp
      

    ENDDO
    
    
    out = Ayprime + 8d0*add_on_u
    
  END FUNCTION waters_axisymmetric_GradV_yx
  
  

  SUBROUTINE Waters (REYNLD,EL1,TIMEW,NTERMS,VISC1,VISC2,hits,cxx,cxy,utrans,t11trans,t12trans,yin)
    IMPLICIT NONE   

    INTEGER, INTENT(IN) :: NTERMS,hits

    DOUBLE PRECISION, INTENT(IN) ::REYNLD,EL1,TIMEW,VISC1,VISC2
    DOUBLE PRECISION :: REYNLD1
    DOUBLE PRECISION, DIMENSION(nptot), INTENT(IN) :: yin,cxx,cxy
    DOUBLE PRECISION, DIMENSION(nptot), INTENT(OUT) :: utrans,t11trans,t12trans

    INTEGER :: ro
    INTEGER NOD,NODDU,FSTEP
    INTEGER I, I1, N2, M, J, NT
    DOUBLE PRECISION SST1, SST2, SST11, CONST1
    DOUBLE PRECISION SN,SM,T1,TL
    DOUBLE PRECISION AP,BQ
    DOUBLE PRECISION YY,U1,U2,U2T,FAC
    DOUBLE PRECISION NPI,NP2,NP3,CNP,SNP,UMAX
    DOUBLE PRECISION MPI,MP2,CMP,SUMFA,SUMFA2,SUMF2,SUMF22
    DOUBLE PRECISION AN,BN,CNN,ASN,BSN,PNA,QN,PNL,QNL,PNT,QNT
    DOUBLE PRECISION AM,BM,CM,ASM,BSM,PM,QM,PML,QML,PMT,QMT
    DOUBLE PRECISION ANP,BNQ
    DOUBLE PRECISION WKN,WK0,WK1,WK2
    DOUBLE PRECISION PNPM ,PNQM ,QNPM ,QNQM
    DOUBLE PRECISION PNPMT,PNQMT,QNPMT,QNQMT
    DOUBLE PRECISION ANM,BNM,DNM,ENM,PMLAN 
    DOUBLE PRECISION ASMAN,ASMBN,BSMAN,BSMBN
    DOUBLE PRECISION UUU,TXY1,TXX1,TXX2,TXX3,TXX4,TXX5,TXX6
    DOUBLE PRECISION SUMUU,SUMXY,SUMXX1,SUMXX2,SUMXX3,SUMXX4
    
    REYNLD1=1d0 
      
    IF(REYNLD.ne.0d0) THEN
      T1=TIMEW/REYNLD
    ELSE
      T1=TIMEW
    ENDIF

    IF (EL1.ne.0d0) THEN
      CONST1 = -VISC1/(EL1*REYNLD1)
      IF(REYNLD.ne.0d0) THEN
	SST1 =EL1/REYNLD
      ELSE
	SST1 =EL1
      ENDIF
      SST2 =EL1*VISC2
      SST11=1.D0/SST1
      TL=T1*SST11
    ENDIF

    UMAX=1.D0

!      Starting of loop of inlet and outlet nodes
!      ------------------------------------------

    DO ro=1,hits
!       IF (coordflag.eq.1) THEN ! Axisymmetric
! 	YY=(yin(ro)+1.D0)/2.D0 !adjust the 2.D0 for the height of the channel.
!       ELSE 
	YY=yin(ro)
!       ENDIF
      
      U2=-4.D0*(2.D0*YY-1.D0)*UMAX
      U1=-4.D0*YY*(YY-1.D0)*UMAX
      SUMUU=0.D0
      SUMXY=0.D0
      SUMXX1=0.D0
      SUMXX2=0.D0
      SUMXX3=0.D0
      SUMXX4=0.D0
      SUMFA=32.D0
      
      SUMFA=SUMFA*UMAX
      SUMF2=SUMFA/2.D0
      SUMFA2=SUMFA*SUMFA
      SUMF22=SUMF2*SUMF2

      U2T=U2*T1

      DO N2=1,NTERMS
	NPI=(2.D0*DFLOAT(N2)-1.D0)*PI
        NP2=NPI*NPI
        NP3=NP2*NPI
        SNP=DSIN(NPI*YY)/NP3
        CNP=DCOS(NPI*YY)/NP2

        IF (EL1.ne.0.D0) THEN
	  AN=1.D0+SST2*NP2
          CNN=AN*AN-4.D0*SST1*NP2
          AN=0.5D0*AN*SST11
          IF(CNN.LT.0.D0) THEN
	    BN=DSQRT(-CNN)
            SN=1.D0+NP2*(SST2-2.D0*SST1)
            BSN=SN/BN
            BN=0.5D0*BN*SST11
            
            PNA=-AN
            QN=BN
            PNT=PNA*T1
            QNT=QN*T1
            AP=BSN
            BQ=1.D0
            FAC=1.D0
!*  evaluates trig. terms for u
            FAC=1.D0
            CALL Get_Gn_Sin(UUU,FAC,PNT,AP,BQ,QNT)
            SUMUU=SUMUU+SNP*UUU

            PNL=PNA+SST11
            AP=BN+PNL*BSN
            BQ=PNL-BN*BSN
            WK1=1.D0/(PNL*PNL+BN*BN)
            FAC=2.D0*WK1
!*  evaluates trig. terms for Txy/2nd trig. terms for Tyy
            CALL Get_Gn_Sin(TXY1,FAC,PNT,AP,BQ,QNT)
            SUMXY=SUMXY+CNP*TXY1

            AP=BN-AN*BSN
            BQ=-AN-BN*BSN
            PNT=PNT-TL
            WK0=1.D0/(AN*AN+BN*BN)
            FAC=2.D0*WK0
!*  evaluates 1st trig. terms for Txx
            CALL Get_Gn_Sin(TXX1,FAC,PNT,AP,BQ,QNT)
            SUMXX1=SUMXX1+CNP*TXX1

            ANP=AP
            BNQ=BQ
            WKN=WK0
          ELSE    !CNN
	    BN=DSQRT(CNN)
	    SN=(1.D0+NP2*(SST2-2.D0*SST1))/BN
	    BN=0.5D0*BN*SST11
	      
	    ASN=1.D0+SN
	    BSN=1.D0-SN
	    PNA=-AN+BN
	    QN=-AN-BN
	    AP=ASN/PNA
	    BQ=BSN/QN
	    PNT=PNA*T1
	    QNT=QN*T1
!*  evaluates exp. terms for u
	    FAC=0.5D0
	    CALL Get_Gn_Exp(UUU,FAC,ASN,PNT,BSN,QNT)
	    SUMUU=SUMUU+SNP*UUU
	      
	    PNL=PNA+SST11
	    QNL=QN+SST11
	    ASN=ASN/PNL
	    BSN=BSN/QNL
!*  evaluates exp. terms for Txy/2nd exponential terms for Tyy
	    FAC=1.D0
	    CALL Get_Gn_Exp(TXY1,FAC,ASN,PNT,BSN,QNT)
	    SUMXY=SUMXY+CNP*TXY1
	      
	    ASN=AP
	    BSN=BQ
	    PNT=PNT-TL
	    QNT=QNT-TL
!*  evaluates 1st exp. terms for Txx
	    FAC=1.D0
	    CALL Get_Gn_Exp(TXX1,FAC,ASN,PNT,BSN,QNT)
	    SUMXX1=SUMXX1+CNP*TXX1
	      
	    ANP=AP
	    BNQ=BQ
	    ASN=1.D0+SN
	    BSN=1.D0-SN
          ENDIF   !CNN
          DO M=1,NTERMS
	    MPI=(2.D0*DFLOAT(M)-1.D0)*PI
	    MP2=MPI*MPI
	    CMP=DCOS(MPI*YY)/MP2
	      
	    AM =1.D0+SST2*MP2
	    CM =AM*AM-4.D0*SST1*MP2
	    AM=0.5D0*AM*SST11
!*** CASE 1 ***
	    IF(CM.GT.0.D0) THEN
	      BM=DSQRT(CM)
	      SM=(1.D0+MP2*(SST2-2.D0*SST1))/BM
	      ASM=1.D0+SM
	      BSM=1.D0-SM
	      BM=0.5D0*BM*SST11
	      PM=-AM+BM
	      QM=-AM-BM
	      PML=PM+SST11
	      QML=QM+SST11
!** single sum Txx 
	      IF(N2.EQ.1) THEN
		PMT=-TL
		QMT=0.D0
		AP=1.D0
		BQ=0.D0
		FAC=(ASM/PML+BSM/QML)
		CALL Get_Gn_Exp(TXX2,FAC,AP,PMT,BQ,QMT)
		SUMXX2=SUMXX2+CMP*TXX2
		
		PMT=PM*T1
		QMT=QM*T1
		AP=ASM/(PML*PML)
		BQ=BSM/(QML*QML)
!*  evaluates 3rd exp. terms for Txx
		FAC=1.D0
		CALL Get_Gn_Exp(TXX3,FAC,AP,PMT,BQ,QMT)
		SUMXX3=SUMXX3+CMP*TXX3
	      ENDIF
	    ENDIF
!** double sum Txx
	    IF(CNN.GT.0.D0.AND.CM.GT.0.D0) THEN
	      FAC=(ASM/PML+BSM/QML)
	      CALL Get_Gn_Exp(TXX6,FAC,ANP,PNT,BNQ,QNT)
	      
	      PNPM=PNA+PM
	      PNQM=PNA+QM
	      QNPM=QN+PM
	      QNQM=QN+QM
	      PNPMT=PNPM*T1
	      PNQMT=PNQM*T1
	      QNPMT=QNPM*T1
	      QNQMT=QNQM*T1
	      ASMAN=ASM*ASN/(PML*(PML+PNA))
	      BSMAN=BSM*ASN/(QML*(QML+PNA))
	      ASMBN=ASM*BSN/(PML*(PML+QN))
	      BSMBN=BSM*BSN/(QML*(QML+QN))
!* evaluates 4th exp. terms for Txx
	      FAC=1.D0
	      CALL Get_Gn_Exp(TXX4,FAC,ASMAN,PNPMT,BSMAN,PNQMT)
	      CALL Get_Gn_Exp(TXX5,FAC,ASMBN,QNPMT,BSMBN,QNQMT)
	      SUMXX4=SUMXX4+CNP*CMP*(TXX4+TXX5-TXX6)
	    ENDIF
!*** CASE4 ***      IF(CNN.LT.0.D0.AND.CM.LT.0.D0)
	    IF(CM.LT.0.D0) THEN
	      BM=DSQRT(-CM)
	      SM=1.D0+MP2*(SST2-2.D0*SST1)
	      BSM=SM/BM
	      BM=0.5D0*BM*SST11
	      PML=-AM+SST11
	      BNM=PML-BM*BSM
	      WK0=PML*PML+BM*BM
!** single sum Txx
	      IF(N2.EQ.1) THEN
		PMT=-TL
		QMT=0.D0
		AP=1.D0
		BQ=0.D0
		FAC=2.D0*BNM/WK0
		CALL Get_Gn_Exp(TXX2,FAC,AP,PMT,BQ,QMT)
		SUMXX2=SUMXX2+CMP*TXX2
		
		PM=-AM
		QM=BM
		PMT=PM*T1
		QMT=QM*T1
		ANM=BM+PML*BSM
		AP=PML*ANM+BM*BNM
		BQ=-BM*ANM+PML*BNM
		WK1=1.D0/WK0
		WK2=WK1*WK1
		FAC=2.D0*WK2
!*  evaluates 3rd trig. terms for Txx
		CALL Get_Gn_Sin(TXX3,FAC,PMT,AP,BQ,QMT)
		SUMXX3=SUMXX3+CMP*TXX3
	      ENDIF
	    ENDIF
!** double sum Txx
	    IF(CNN.LT.0.D0.AND.CM.LT.0.D0) THEN
	      FAC=4.D0*BNM*WKN/WK0
	      CALL Get_Gn_Sin(TXX6,FAC,PNT,ANP,BNQ,QNT)
	      
	      PMLAN=PML-AN
	      PNPM=-AN-AM
	      QNPM=BN+BM
	      QNQM=BN-BM
	      WK1=1.D0/((PMLAN*PMLAN+QNPM*QNPM)*WK0)
	      WK2=1.D0/((PMLAN*PMLAN+QNQM*QNQM)*WK0)
	      PNPMT=PNPM*T1
	      QNPMT=QNPM*T1
	      QNQMT=QNQM*T1
	      ANM=0.5*( BM+PML*BSM+PML*BSN-BM*BSN*BSM)
	      BNM=0.5*(-BM-PML*BSM+PML*BSN-BM*BSN*BSM)
	      DNM=0.5*( PML-BM*BSM-BM*BSN-PML*BSN*BSM)
	      ENM=0.5*( PML-BM*BSM+BM*BSN+PML*BSN*BSM)
	      ASMAN= ANM*PMLAN+DNM*QNPM
	      BSMAN=-ANM*QNPM+DNM*PMLAN
	      ASMBN= BNM*PMLAN+ENM*QNQM
	      BSMBN=-BNM*QNQM+ENM*PMLAN
!*  evaluates 4th trig. terms for Txx
	      FAC=4.D0*WK1
	      CALL Get_Gn_Sin(TXX4,FAC,PNPMT,ASMAN,BSMAN,QNPMT)
	      FAC=4.D0*WK2
	      CALL Get_Gn_Sin(TXX5,FAC,PNPMT,ASMBN,BSMBN,QNQMT)
	      SUMXX4=SUMXX4+CNP*CMP*(TXX4+TXX5-TXX6)
	    ENDIF

	  ENDDO  !M,NTERMS

	ELSE  !EL1
	  PNT=-NP2*T1
	  CALL Get_Gn_Exp(UUU,1.D0,1.D0,PNT,0.D0,0.D0)
	  SUMUU=SUMUU+SNP*UUU
	ENDIF !EL1
      ENDDO  !N2,NTERMS

!     velocity x-component

      utrans(ro)=U1-SUMFA*SUMUU

      IF (EL1.ne.0d0) THEN
!     stress xx-component
	t11trans(ro) = (2.D0*VISC1*U2*&
			(U2*SST1-SUMF2*SUMXY-U2T*DEXP(-TL)+SUMF2*SUMXX1) &
			+SUMFA*CONST1*(U2*SUMXX3-U2T*SUMXX2) &
			-2.D0*SUMF22*CONST1*SUMXX4) &
			 +cxx(ro)*DEXP(-TL)
!     stress xy or rz-component
	t12trans(ro) = -CONST1*(U2*SST1-SUMF2*SUMXY) &
			+cxy(ro)*DEXP(-TL)
      ENDIF
    ENDDO
	    ! I Loop End of nodes

    RETURN
  END SUBROUTINE Waters
!=======================================================================
  SUBROUTINE Get_Gn_Exp (GN,FAC,AS,PNA,BS,QN)
    DOUBLE PRECISION GN,FAC,AS,BS,PNA,QN

   GN=FAC*(AS*DEXP(PNA)+BS*DEXP(QN))

   RETURN
  END SUBROUTINE Get_Gn_Exp
!=======================================================================
  SUBROUTINE Get_Gn_Sin (GN,FAC,PNA,AS,BS,QN)

    DOUBLE PRECISION GN,FAC,AS,BS,PNA,QN

    GN=FAC*DEXP(PNA)*(AS*DSIN(QN)+BS*DCOS(QN))

    RETURN
  END SUBROUTINE Get_Gn_Sin
!=======================================================================
  SUBROUTINE WatersIni(REYNLD,EL1,NTERMS,VISC1,VISC2,hits,cxx,cxy,yin)
    
    IMPLICIT NONE   
      
    INTEGER, INTENT(IN) :: NTERMS,hits
    DOUBLE PRECISION, INTENT(IN) :: REYNLD,EL1,VISC1,VISC2
    DOUBLE PRECISION, DIMENSION(nptot), INTENT(IN) :: yin
    DOUBLE PRECISION, DIMENSION(nptot), INTENT(OUT) :: cxx,cxy

    INTEGER :: ro
    INTEGER :: NOD, NODDU, FSTEP, I, I1, J, N2, M
    DOUBLE PRECISION :: REYNLD1
    DOUBLE PRECISION :: SN, SM
    DOUBLE PRECISION :: AP, BQ
    DOUBLE PRECISION :: YY, U1, U2, FAC
    DOUBLE PRECISION :: NPI, NP2, NP3, CNP, SNP, UMAX
    DOUBLE PRECISION :: MPI,MP2,CMP,SUMFA,SUMFA2,SUMF2,SUMF22
    DOUBLE PRECISION :: AN,BN,CNN,ASN,BSN,PNA,QN,PNL,QNL
    DOUBLE PRECISION :: AM,BM,CM,ASM,BSM,PM,QM,PML,QML,PMT,QMT
    DOUBLE PRECISION :: ANP,BNQ
    DOUBLE PRECISION :: WKN,WK0,WK1,WK2
    DOUBLE PRECISION :: PNPM ,PNQM ,QNPM ,QNQM
    DOUBLE PRECISION :: PNPMT,PNQMT,QNPMT,QNQMT
    DOUBLE PRECISION :: ANM,BNM,DNM,ENM,PMLAN
    DOUBLE PRECISION :: ASMAN,ASMBN,BSMAN,BSMBN
    DOUBLE PRECISION :: SST1, SST11, SST2, CONST1

    REYNLD1=1d0
	
    IF (EL1.ne.0d0) THEN
      CONST1 = -VISC1/(EL1*REYNLD1)
      IF(REYNLD.ne.0d0) THEN
	SST1 =EL1/REYNLD
      ELSE
	SST1 =EL1
      ENDIF
	SST2  = EL1*VISC2
	SST11 = 1.0D0/SST1
    ENDIF

    UMAX=1.D0
    
    DO ro=1,hits
!       IF (coordflag.eq.1) THEN
! 	YY=(yin(ro)+1.D0)/2.D0
!       ELSE
	YY=yin(ro)
!       ENDIF
      
      U2=-4.D0*(2.D0*YY-1.D0)*UMAX
      U1=-4.D0*YY*(YY-1.D0)*UMAX
      
      cxy(ro) = 0.0D0
      cxx(ro) = 0.0D0

      UMAX  = 1.0D0
      SUMFA = 32.0D0
      SUMFA  = SUMFA*UMAX
      SUMF2  = 0.5d0*SUMFA
      SUMFA2 = SUMFA*SUMFA
      SUMF22 = SUMF2*SUMF2
      
      DO N2=1,NTERMS
	NPI = (2.0D0*DFLOAT(N2)-1.0D0)*PI
	NP2 = NPI*NPI
	NP3 = NP2*NPI
	SNP = DSIN(NPI*YY) / NP3
	CNP = DCOS(NPI*YY) / NP2
	
        IF (EL1.NE.0.D0) THEN
	  AN  = 1.D0 + SST2*NP2
          CNN = AN*AN - 4.0D0*SST1*NP2
          AN  = 0.5D0*AN*SST11
          IF ( CNN.lt.0 ) THEN
	    BN  = DSQRT(-CNN)
	    SN  = 1.0D0 + NP2*( SST2 - 2.0D0*SST1 )
	    BSN = SN / BN
	    BN  = 0.5D0*BN*SST11
          
	    PNA =-AN
	    QN  = BN
	    AP  = BSN
	    BQ  = 1.0D0
	    FAC = 1.0D0
          
	    PNL = PNA+SST11
	    AP  = BN+PNL*BSN
	    BQ  = PNL-BN*BSN
	    WK1 = 1.0D0/(PNL*PNL+BN*BN)
	    FAC = 2.0D0*WK1
!*  sum for constant CXY
	    cxy(ro) = cxy(ro) - SUMF2*CNP*BQ*FAC
          
	    AP  = BN - AN*BSN
	    BQ  =-AN - BN*BSN
	    WK0 = 1.0D0 / ( AN*AN + BN*BN )
	    FAC = 2.0D0*WK0
!*  sum for constant CXX
	    cxx(ro) = cxx(ro) - SUMFA*VISC1*U2*CNP*BQ*FAC

	    ANP = AP
	    BNQ = BQ
	    WKN = WK0
	  ELSE  !CNN
	    BN = DSQRT(CNN)
	    SN = ( 1.0D0 + NP2*( SST2 - 2.0D0*SST1 ) ) / BN
	    BN = 0.5D0*BN*SST11

	    ASN = 1.0D0 + SN
	    BSN = 1.0D0 - SN
	    PNA =-AN + BN
	    QN  =-AN - BN
	    AP  = ASN / PNA
	    BQ  = BSN / QN
 
	    PNL = PNA + SST11
	    QNL = QN + SST11
	    ASN = ASN / PNL
	    BSN = BSN / QNL
!*  sum for constant CXY
	    cxy(ro) = cxy(ro) - SUMF2*CNP*( ASN + BSN )
!*  evaluates exp. terms for Txy/2nd exponential terms for Tyy
	    ASN = AP
	    BSN = BQ
!*  sum for constant CXX
	    cxx(ro) = cxx(ro) - SUMFA*VISC1*U2*CNP*( ASN + BSN )
          
	    ANP = AP
	    BNQ = BQ
	    ASN = 1.0D0 + SN
	    BSN = 1.0D0 - SN
	  ENDIF  !CNN

	  DO M=1,NTERMS
	    MPI = ( 2.0D0*DFLOAT(M) - 1.0D0 )*PI
	    MP2 = MPI*MPI
	    CMP = DCOS(MPI*YY) / MP2
	    AM = 1.0D0 + SST2*MP2
	    CM = AM*AM - 4.0D0*SST1*MP2
	    AM = 0.5D0*AM*SST11
!*** CASE 1 ***	
	    IF ( CM.gt.0.0D0 ) THEN 
	      BM  = DSQRT(CM)
	      SM  = ( 1.0D0 + MP2*( SST2 - 2.0D0*SST1 ) ) / BM
	      ASM = 1.0D0 + SM
	      BSM = 1.0D0 - SM
	      BM  = 0.5D0*BM*SST11
	      PM  =-AM + BM
	      QM  =-AM - BM
	      PML = PM + SST11
	      QML = QM + SST11
	    
	      IF ( N2.EQ.1 ) THEN
		AP  = ASM/(PML*PML)
		BQ  = BSM/(QML*QML)
!*sum for constant CXX
		cxx(ro) = cxx(ro) - SUMFA*CONST1*U2*CMP*( AP + BQ )
	      ENDIF
	    ENDIF

	    IF ((CNN.gt.0.0D0).AND.(CM.gt.0.0D0)) THEN
	      FAC  = (ASM/PML+BSM/QML)
	  !* sum for constant CXX
	      cxx(ro) = cxx(ro) - 2.0D0*SUMF22*CONST1*CNP*CMP*FAC*(ANP+BNQ)
	  
	      ASMAN = ASM*ASN/(PML*(PML+PNA))
	      BSMAN = BSM*ASN/(QML*(QML+PNA))
	      ASMBN = ASM*BSN/(PML*(PML+QN))
	      BSMBN = BSM*BSN/(QML*(QML+QN))
!* sum for constant CXX
	      cxx(ro) = cxx(ro) + 2.0D0*SUMF22*CONST1*CNP*CMP*(ASMAN+BSMAN+ASMBN+BSMBN)
	    ENDIF
!*** CASE4 ***
	    IF ( CM.lt.0.D0 ) THEN
	      BM  = DSQRT(-CM)
	      SM  = 1.0D0+MP2*(SST2-2.0D0*SST1)
	      BSM = SM/BM
	      BM  = 0.5D0*BM*SST11
	      PML =-AM+SST11
	      BNM = PML-BM*BSM
	      WK0 = PML*PML+BM*BM
!** single sum Txx
	      IF ( N2.EQ.1 ) THEN
		ANM = BM+PML*BSM
		BQ  =-BM*ANM+PML*BNM
		WK1 = 1.0D0/WK0
		WK2 = WK1*WK1
		FAC = 2.0D0*WK2
!*  sum for constant CXX
		cxx(ro) = cxx(ro) - SUMFA*CONST1*U2*CMP*BQ*FAC
	      ENDIF
	    ENDIF
	    IF ((CNN.lt.0.0D0).AND.(CM.lt.0.0D0)) THEN
	      FAC  = 4.0D0*BNM*WKN/WK0
!* sum for constant CXX
	      cxx(ro) = cxx(ro) - 2.0D0*SUMF22*CONST1*CNP*CMP*FAC*BNQ
	  
	      PMLAN = PML-AN
	      QNPM  = BN+BM
	      QNQM  = BN-BM
	      WK1   = 1.0D0/((PMLAN*PMLAN+QNPM*QNPM)*WK0)
	      WK2   = 1.0D0/((PMLAN*PMLAN+QNQM*QNQM)*WK0)
	      ANM   = 0.5d0*( BM+PML*BSM+PML*BSN-BM*BSN*BSM)
	      BNM   = 0.5d0*(-BM-PML*BSM+PML*BSN-BM*BSN*BSM)
	      DNM   = 0.5d0*( PML-BM*BSM-BM*BSN-PML*BSN*BSM)
	      ENM   = 0.5d0*( PML-BM*BSM+BM*BSN+PML*BSN*BSM)
	      BSMAN =-ANM*QNPM + DNM*PMLAN
	      BSMBN =-BNM*QNQM + ENM*PMLAN
!* sum for constant CXX
	      cxx(ro) = cxx(ro)+2.0D0*SUMFA2*CONST1*CNP*CMP*(BSMAN*WK1+BSMBN*WK2)
	    ENDIF
	  ENDDO  !M
	ENDIF !EL1
	
      ENDDO  !N2

!**   stress components
      IF (EL1.ne.0.D0) THEN
	cxx(ro) = cxx(ro) - 2.0D0*U2*VISC1*( U2*SST1 + cxy(ro) )
	cxy(ro) = CONST1*( U2*SST1 + cxy(ro) )
      ENDIF
       
    ENDDO    ! I Loop
    RETURN

  END SUBROUTINE WatersIni
!=======================================================================
!      call WatersIni (Re,We,20,1d0-beta,beta,hits,
!     $                 cxx,cxy,yin)
!      call Waters (Re,We,timew,20,1d0-beta,beta,hits,
!     $             cxx,cxy,utrans,t11trans,t12trans,yin)
!=======================================================================





! Axisymmetric version is fairly tricky - we instead just use the velocity and not bother with stress components. See function above
!  
!   SUBROUTINE Waters_Axisymmetric (REYNLD,EL1,TIMEW,NTERMS,VISC1,VISC2,hits,cxx,cxy,utrans,t11trans,t12trans,yin)
!     IMPLICIT NONE   
! 
!     INTEGER, INTENT(IN) :: NTERMS,hits
! 
!     DOUBLE PRECISION, INTENT(IN) ::REYNLD,EL1,TIMEW,VISC1,VISC2
!     DOUBLE PRECISION :: REYNLD1
!     DOUBLE PRECISION, DIMENSION(nptot), INTENT(IN) :: yin,cxx,cxy
!     DOUBLE PRECISION, DIMENSION(nptot), INTENT(OUT) :: utrans,t11trans,t12trans
! 
!     INTEGER :: ro
!     INTEGER NOD,NODDU,FSTEP
!     INTEGER I, I1, N2, M, J, NT
!     DOUBLE PRECISION SST1, SST2, SST11, CONST1
!     DOUBLE PRECISION SN,SM,T1,TL
!     DOUBLE PRECISION AP,BQ
!     DOUBLE PRECISION YY,U1,U2,U2T,FAC
!     DOUBLE PRECISION NPI,NP2,NP3,CNP,SNP,UMAX
!     DOUBLE PRECISION MPI,MP2,CMP,SUMFA,SUMFA2,SUMF2,SUMF22
!     DOUBLE PRECISION AN,BN,CNN,ASN,BSN,PNA,QN,PNL,QNL,PNT,QNT
!     DOUBLE PRECISION AM,BM,CM,ASM,BSM,PM,QM,PML,QML,PMT,QMT
!     DOUBLE PRECISION ANP,BNQ
!     DOUBLE PRECISION WKN,WK0,WK1,WK2
!     DOUBLE PRECISION PNPM ,PNQM ,QNPM ,QNQM
!     DOUBLE PRECISION PNPMT,PNQMT,QNPMT,QNQMT
!     DOUBLE PRECISION ANM,BNM,DNM,ENM,PMLAN
!     DOUBLE PRECISION ASMAN,ASMBN,BSMAN,BSMBN
!     DOUBLE PRECISION UUU,TXY1,TXX1,TXX2,TXX3,TXX4,TXX5,TXX6
!     DOUBLE PRECISION SUMUU,SUMXY,SUMXX1,SUMXX2,SUMXX3,SUMXX4
!     
!     DOUBLE PRECISION, DIMENSION(20) :: besj0_zeros
!     
! ! First 20 zeros of the bessel function J_0
! ! These are required for the sums and will form our set "Z_n" (see Waters & King 1971)
!     besj0_zeros = (/ 2.4048255576957729, 5.5200781102863106, 8.6537279129110125, 11.791534439014281, 14.930917708487787, &
! 				    18.071063967910924, 21.211636629879258, 24.352471530749302, 27.493479132040253, 30.634606468431976, &
! 				    33.775820213573567, 36.917098353664045, 40.058425764628240, 43.199791713176730, 46.341188371661815, &
! 				    49.482609897397815, 52.624051841114998, 55.765510755019982, 58.906983926080940, 62.048469190227166 /)
! 				    
! 				    
!     REYNLD1=1d0 
!       
!     IF(REYNLD.ne.0d0) THEN
!       T1=TIMEW/REYNLD
!     ELSE
!       T1=TIMEW
!     ENDIF
! 
!     IF (EL1.ne.0d0) THEN
!       CONST1 = -VISC1/(EL1*REYNLD1)
!       IF(REYNLD.ne.0d0) THEN
! 	SST1 =EL1/REYNLD
!       ELSE
! 	SST1 =EL1
!       ENDIF
!       SST2 =EL1*VISC2
!       SST11=1.D0/SST1
!       TL=T1*SST11
!     ENDIF
! 
!     UMAX=1.D0
! 
! !      Starting of loop of inlet and outlet nodes
! !      ------------------------------------------
! 
!     DO ro=1,hits
! ! Adjustment for axisymmetric case.
!       YY=yin(ro)
!       U2=(-2.D0*YY)*UMAX
!       U1=(1.D0-YY*YY)*UMAX
!       SUMUU=0.D0
!       SUMXY=0.D0
!       SUMXX1=0.D0
!       SUMXX2=0.D0
!       SUMXX3=0.D0
!       SUMXX4=0.D0
!       SUMFA=8.D0!32.D0 ! Change for Axisymmetric
!       
!       SUMFA=SUMFA*UMAX
!       SUMF2=SUMFA/2.D0
!       SUMFA2=SUMFA*SUMFA
!       SUMF22=SUMF2*SUMF2
! 
!       U2T=U2*T1
! 
!       DO N2=1,NTERMS
! ! 	NPI=besj0_zeros(N2)  !(2.D0*DFLOAT(N2)-1.D0)*PI
!         NP2=NPI*NPI
!         NP3=NP2*NPI
! !         SNP=BESJ !DSIN(NPI*YY)/NP3
!         CNP=DCOS(NPI*YY)/NP2
! 
!         IF (EL1.ne.0.D0) THEN
! 	  AN=1.D0+SST2*NP2
!           CNN=AN*AN-4.D0*SST1*NP2
!           AN=0.5D0*AN*SST11
!           IF(CNN.LT.0.D0) THEN
! 	    BN=DSQRT(-CNN)
!             SN=1.D0+NP2*(SST2-2.D0*SST1)
!             BSN=SN/BN
!             BN=0.5D0*BN*SST11
!             
!             PNA=-AN
!             QN=BN
!             PNT=PNA*T1
!             QNT=QN*T1
!             AP=BSN
!             BQ=1.D0
!             FAC=1.D0
! !*  evaluates trig. terms for u
!             FAC=1.D0
!             CALL Get_Gn_Sin(UUU,FAC,PNT,AP,BQ,QNT)
!             SUMUU=SUMUU+SNP*UUU
! 
!             PNL=PNA+SST11
!             AP=BN+PNL*BSN
!             BQ=PNL-BN*BSN
!             WK1=1.D0/(PNL*PNL+BN*BN)
!             FAC=2.D0*WK1
! !*  evaluates trig. terms for Txy/2nd trig. terms for Tyy
!             CALL Get_Gn_Sin(TXY1,FAC,PNT,AP,BQ,QNT)
!             SUMXY=SUMXY+CNP*TXY1
! 
!             AP=BN-AN*BSN
!             BQ=-AN-BN*BSN
!             PNT=PNT-TL
!             WK0=1.D0/(AN*AN+BN*BN)
!             FAC=2.D0*WK0
! !*  evaluates 1st trig. terms for Txx
!             CALL Get_Gn_Sin(TXX1,FAC,PNT,AP,BQ,QNT)
!             SUMXX1=SUMXX1+CNP*TXX1
! 
!             ANP=AP
!             BNQ=BQ
!             WKN=WK0
!           ELSE    !CNN
! 	    BN=DSQRT(CNN)
! 	    SN=(1.D0+NP2*(SST2-2.D0*SST1))/BN
! 	    BN=0.5D0*BN*SST11
! 	      
! 	    ASN=1.D0+SN
! 	    BSN=1.D0-SN
! 	    PNA=-AN+BN
! 	    QN=-AN-BN
! 	    AP=ASN/PNA
! 	    BQ=BSN/QN
! 	    PNT=PNA*T1
! 	    QNT=QN*T1
! !*  evaluates exp. terms for u
! 	    FAC=0.5D0
! 	    CALL Get_Gn_Exp(UUU,FAC,ASN,PNT,BSN,QNT)
! 	    SUMUU=SUMUU+SNP*UUU
! 	      
! 	    PNL=PNA+SST11
! 	    QNL=QN+SST11
! 	    ASN=ASN/PNL
! 	    BSN=BSN/QNL
! !*  evaluates exp. terms for Txy/2nd exponential terms for Tyy
! 	    FAC=1.D0
! 	    CALL Get_Gn_Exp(TXY1,FAC,ASN,PNT,BSN,QNT)
! 	    SUMXY=SUMXY+CNP*TXY1
! 	      
! 	    ASN=AP
! 	    BSN=BQ
! 	    PNT=PNT-TL
! 	    QNT=QNT-TL
! !*  evaluates 1st exp. terms for Txx
! 	    FAC=1.D0
! 	    CALL Get_Gn_Exp(TXX1,FAC,ASN,PNT,BSN,QNT)
! 	    SUMXX1=SUMXX1+CNP*TXX1
! 	      
! 	    ANP=AP
! 	    BNQ=BQ
! 	    ASN=1.D0+SN
! 	    BSN=1.D0-SN
!           ENDIF   !CNN
!           DO M=1,NTERMS
! 	    MPI=(2.D0*DFLOAT(M)-1.D0)*PI
! 	    MP2=MPI*MPI
! 	    CMP=DCOS(MPI*YY)/MP2
! 	      
! 	    AM =1.D0+SST2*MP2
! 	    CM =AM*AM-4.D0*SST1*MP2
! 	    AM=0.5D0*AM*SST11
! !*** CASE 1 ***
! 	    IF(CM.GT.0.D0) THEN
! 	      BM=DSQRT(CM)
! 	      SM=(1.D0+MP2*(SST2-2.D0*SST1))/BM
! 	      ASM=1.D0+SM
! 	      BSM=1.D0-SM
! 	      BM=0.5D0*BM*SST11
! 	      PM=-AM+BM
! 	      QM=-AM-BM
! 	      PML=PM+SST11
! 	      QML=QM+SST11
! !** single sum Txx 
! 	      IF(N2.EQ.1) THEN
! 		PMT=-TL
! 		QMT=0.D0
! 		AP=1.D0
! 		BQ=0.D0
! 		FAC=(ASM/PML+BSM/QML)
! 		CALL Get_Gn_Exp(TXX2,FAC,AP,PMT,BQ,QMT)
! 		SUMXX2=SUMXX2+CMP*TXX2
! 		
! 		PMT=PM*T1
! 		QMT=QM*T1
! 		AP=ASM/(PML*PML)
! 		BQ=BSM/(QML*QML)
! !*  evaluates 3rd exp. terms for Txx
! 		FAC=1.D0
! 		CALL Get_Gn_Exp(TXX3,FAC,AP,PMT,BQ,QMT)
! 		SUMXX3=SUMXX3+CMP*TXX3
! 	      ENDIF
! 	    ENDIF
! !** double sum Txx
! 	    IF(CNN.GT.0.D0.AND.CM.GT.0.D0) THEN
! 	      FAC=(ASM/PML+BSM/QML)
! 	      CALL Get_Gn_Exp(TXX6,FAC,ANP,PNT,BNQ,QNT)
! 	      
! 	      PNPM=PNA+PM
! 	      PNQM=PNA+QM
! 	      QNPM=QN+PM
! 	      QNQM=QN+QM
! 	      PNPMT=PNPM*T1
! 	      PNQMT=PNQM*T1
! 	      QNPMT=QNPM*T1
! 	      QNQMT=QNQM*T1
! 	      ASMAN=ASM*ASN/(PML*(PML+PNA))
! 	      BSMAN=BSM*ASN/(QML*(QML+PNA))
! 	      ASMBN=ASM*BSN/(PML*(PML+QN))
! 	      BSMBN=BSM*BSN/(QML*(QML+QN))
! !* evaluates 4th exp. terms for Txx
! 	      FAC=1.D0
! 	      CALL Get_Gn_Exp(TXX4,FAC,ASMAN,PNPMT,BSMAN,PNQMT)
! 	      CALL Get_Gn_Exp(TXX5,FAC,ASMBN,QNPMT,BSMBN,QNQMT)
! 	      SUMXX4=SUMXX4+CNP*CMP*(TXX4+TXX5-TXX6)
! 	    ENDIF
! !*** CASE4 ***      IF(CNN.LT.0.D0.AND.CM.LT.0.D0)
! 	    IF(CM.LT.0.D0) THEN
! 	      BM=DSQRT(-CM)
! 	      SM=1.D0+MP2*(SST2-2.D0*SST1)
! 	      BSM=SM/BM
! 	      BM=0.5D0*BM*SST11
! 	      PML=-AM+SST11
! 	      BNM=PML-BM*BSM
! 	      WK0=PML*PML+BM*BM
! !** single sum Txx
! 	      IF(N2.EQ.1) THEN
! 		PMT=-TL
! 		QMT=0.D0
! 		AP=1.D0
! 		BQ=0.D0
! 		FAC=2.D0*BNM/WK0
! 		CALL Get_Gn_Exp(TXX2,FAC,AP,PMT,BQ,QMT)
! 		SUMXX2=SUMXX2+CMP*TXX2
! 		
! 		PM=-AM
! 		QM=BM
! 		PMT=PM*T1
! 		QMT=QM*T1
! 		ANM=BM+PML*BSM
! 		AP=PML*ANM+BM*BNM
! 		BQ=-BM*ANM+PML*BNM
! 		WK1=1.D0/WK0
! 		WK2=WK1*WK1
! 		FAC=2.D0*WK2
! !*  evaluates 3rd trig. terms for Txx
! 		CALL Get_Gn_Sin(TXX3,FAC,PMT,AP,BQ,QMT)
! 		SUMXX3=SUMXX3+CMP*TXX3
! 	      ENDIF
! 	    ENDIF
! !** double sum Txx
! 	    IF(CNN.LT.0.D0.AND.CM.LT.0.D0) THEN
! 	      FAC=4.D0*BNM*WKN/WK0
! 	      CALL Get_Gn_Sin(TXX6,FAC,PNT,ANP,BNQ,QNT)
! 	      
! 	      PMLAN=PML-AN
! 	      PNPM=-AN-AM
! 	      QNPM=BN+BM
! 	      QNQM=BN-BM
! 	      WK1=1.D0/((PMLAN*PMLAN+QNPM*QNPM)*WK0)
! 	      WK2=1.D0/((PMLAN*PMLAN+QNQM*QNQM)*WK0)
! 	      PNPMT=PNPM*T1
! 	      QNPMT=QNPM*T1
! 	      QNQMT=QNQM*T1
! 	      ANM=0.5*( BM+PML*BSM+PML*BSN-BM*BSN*BSM)
! 	      BNM=0.5*(-BM-PML*BSM+PML*BSN-BM*BSN*BSM)
! 	      DNM=0.5*( PML-BM*BSM-BM*BSN-PML*BSN*BSM)
! 	      ENM=0.5*( PML-BM*BSM+BM*BSN+PML*BSN*BSM)
! 	      ASMAN= ANM*PMLAN+DNM*QNPM
! 	      BSMAN=-ANM*QNPM+DNM*PMLAN
! 	      ASMBN= BNM*PMLAN+ENM*QNQM
! 	      BSMBN=-BNM*QNQM+ENM*PMLAN
! !*  evaluates 4th trig. terms for Txx
! 	      FAC=4.D0*WK1
! 	      CALL Get_Gn_Sin(TXX4,FAC,PNPMT,ASMAN,BSMAN,QNPMT)
! 	      FAC=4.D0*WK2
! 	      CALL Get_Gn_Sin(TXX5,FAC,PNPMT,ASMBN,BSMBN,QNQMT)
! 	      SUMXX4=SUMXX4+CNP*CMP*(TXX4+TXX5-TXX6)
! 	    ENDIF
! 
! 	  ENDDO  !M,NTERMS
! 
! 	ELSE  !EL1
! 	  PNT=-NP2*T1
! 	  CALL Get_Gn_Exp(UUU,1.D0,1.D0,PNT,0.D0,0.D0)
! 	  SUMUU=SUMUU+SNP*UUU
! 	ENDIF !EL1
!       ENDDO  !N2,NTERMS
! 
! !     velocity x-component
! 
!       utrans(ro)=U1-SUMFA*SUMUU
! 
!       IF (EL1.ne.0d0) THEN
! !     stress xx-component
! 	t11trans(ro) = (2.D0*VISC1*U2*&
! 			(U2*SST1-SUMF2*SUMXY-U2T*DEXP(-TL)+SUMF2*SUMXX1) &
! 			+SUMFA*CONST1*(U2*SUMXX3-U2T*SUMXX2) &
! 			-2.D0*SUMF22*CONST1*SUMXX4) &
! 			 +cxx(ro)*DEXP(-TL)
! !     stress xy or rz-component
! 	t12trans(ro) = -CONST1*(U2*SST1-SUMF2*SUMXY) &
! 			+cxy(ro)*DEXP(-TL)
!       ENDIF
!     ENDDO
! 	    ! I Loop End of nodes
! 
!     RETURN
!   END SUBROUTINE Waters_Axisymmetric
!   
!     SUBROUTINE WatersIni_Axisymmetric(REYNLD,EL1,NTERMS,VISC1,VISC2,hits,cxx,cxy,yin)
!     
!     IMPLICIT NONE   
!       
!     INTEGER, INTENT(IN) :: NTERMS,hits
!     DOUBLE PRECISION, INTENT(IN) :: REYNLD,EL1,VISC1,VISC2
!     DOUBLE PRECISION, DIMENSION(nptot), INTENT(IN) :: yin
!     DOUBLE PRECISION, DIMENSION(nptot), INTENT(OUT) :: cxx,cxy
! 
!     INTEGER :: ro
!     INTEGER :: NOD, NODDU, FSTEP, I, I1, J, N2, M
!     DOUBLE PRECISION :: REYNLD1
!     DOUBLE PRECISION :: SN, SM
!     DOUBLE PRECISION :: AP, BQ
!     DOUBLE PRECISION :: YY, U1, U2, FAC
!     DOUBLE PRECISION :: NPI, NP2, NP3, CNP, SNP, UMAX
!     DOUBLE PRECISION :: MPI,MP2,CMP,SUMFA,SUMFA2,SUMF2,SUMF22
!     DOUBLE PRECISION :: AN,BN,CNN,ASN,BSN,PNA,QN,PNL,QNL
!     DOUBLE PRECISION :: AM,BM,CM,ASM,BSM,PM,QM,PML,QML,PMT,QMT
!     DOUBLE PRECISION :: ANP,BNQ
!     DOUBLE PRECISION :: WKN,WK0,WK1,WK2
!     DOUBLE PRECISION :: PNPM ,PNQM ,QNPM ,QNQM
!     DOUBLE PRECISION :: PNPMT,PNQMT,QNPMT,QNQMT
!     DOUBLE PRECISION :: ANM,BNM,DNM,ENM,PMLAN
!     DOUBLE PRECISION :: ASMAN,ASMBN,BSMAN,BSMBN
!     DOUBLE PRECISION :: SST1, SST11, SST2, CONST1
! 
!     REYNLD1=1d0
! 	
!     IF (EL1.ne.0d0) THEN
!       CONST1 = -VISC1/(EL1*REYNLD1)
!       IF(REYNLD.ne.0d0) THEN
! 	SST1 =EL1/REYNLD
!       ELSE
! 	SST1 =EL1
!       ENDIF
! 	SST2  = EL1*VISC2
! 	SST11 = 1.0D0/SST1
!     ENDIF
! 
!     UMAX=1.D0
!     
!     DO ro=1,hits
!      
! ! Adjustment for axisymmetric case.
!       YY=(yin(ro)+1.D0)/2.D0
! 
! 
!       U2=-4.D0*(2.D0*YY-1.D0)*UMAX
!       U1=-4.D0*YY*(YY-1.D0)*UMAX
!       
!       cxy(ro) = 0.0D0
!       cxx(ro) = 0.0D0
! 
!       UMAX  = 1.0D0
!       SUMFA = 32.D0 ! Change for Axisymmetric
!       SUMFA  = SUMFA*UMAX
!       SUMF2  = 0.5d0*SUMFA
!       SUMFA2 = SUMFA*SUMFA
!       SUMF22 = SUMF2*SUMF2
!       
!       DO N2=1,NTERMS
! 	NPI = (2.0D0*DFLOAT(N2)-1.0D0)*PI
! 	NP2 = NPI*NPI
! 	NP3 = NP2*NPI
! 	SNP = DSIN(NPI*YY) / NP3
! 	CNP = DCOS(NPI*YY) / NP2
! 	
!         IF (EL1.NE.0.D0) THEN
! 	  AN  = 1.D0 + SST2*NP2
!           CNN = AN*AN - 4.0D0*SST1*NP2
!           AN  = 0.5D0*AN*SST11
!           IF ( CNN.lt.0 ) THEN
! 	    BN  = DSQRT(-CNN)
! 	    SN  = 1.0D0 + NP2*( SST2 - 2.0D0*SST1 )
! 	    BSN = SN / BN
! 	    BN  = 0.5D0*BN*SST11
!           
! 	    PNA =-AN
! 	    QN  = BN
! 	    AP  = BSN
! 	    BQ  = 1.0D0
! 	    FAC = 1.0D0
!           
! 	    PNL = PNA+SST11
! 	    AP  = BN+PNL*BSN
! 	    BQ  = PNL-BN*BSN
! 	    WK1 = 1.0D0/(PNL*PNL+BN*BN)
! 	    FAC = 2.0D0*WK1
! !*  sum for constant CXY
! 	    cxy(ro) = cxy(ro) - SUMF2*CNP*BQ*FAC
!           
! 	    AP  = BN - AN*BSN
! 	    BQ  =-AN - BN*BSN
! 	    WK0 = 1.0D0 / ( AN*AN + BN*BN )
! 	    FAC = 2.0D0*WK0
! !*  sum for constant CXX
! 	    cxx(ro) = cxx(ro) - SUMFA*VISC1*U2*CNP*BQ*FAC
! 
! 	    ANP = AP
! 	    BNQ = BQ
! 	    WKN = WK0
! 	  ELSE  !CNN
! 	    BN = DSQRT(CNN)
! 	    SN = ( 1.0D0 + NP2*( SST2 - 2.0D0*SST1 ) ) / BN
! 	    BN = 0.5D0*BN*SST11
! 
! 	    ASN = 1.0D0 + SN
! 	    BSN = 1.0D0 - SN
! 	    PNA =-AN + BN
! 	    QN  =-AN - BN
! 	    AP  = ASN / PNA
! 	    BQ  = BSN / QN
!  
! 	    PNL = PNA + SST11
! 	    QNL = QN + SST11
! 	    ASN = ASN / PNL
! 	    BSN = BSN / QNL
! !*  sum for constant CXY
! 	    cxy(ro) = cxy(ro) - SUMF2*CNP*( ASN + BSN )
! !*  evaluates exp. terms for Txy/2nd exponential terms for Tyy
! 	    ASN = AP
! 	    BSN = BQ
! !*  sum for constant CXX
! 	    cxx(ro) = cxx(ro) - SUMFA*VISC1*U2*CNP*( ASN + BSN )
!           
! 	    ANP = AP
! 	    BNQ = BQ
! 	    ASN = 1.0D0 + SN
! 	    BSN = 1.0D0 - SN
! 	  ENDIF  !CNN
! 
! 	  DO M=1,NTERMS
! 	    MPI = ( 2.0D0*DFLOAT(M) - 1.0D0 )*PI
! 	    MP2 = MPI*MPI
! 	    CMP = DCOS(MPI*YY) / MP2
! 	    AM = 1.0D0 + SST2*MP2
! 	    CM = AM*AM - 4.0D0*SST1*MP2
! 	    AM = 0.5D0*AM*SST11
! !*** CASE 1 ***	
! 	    IF ( CM.gt.0.0D0 ) THEN 
! 	      BM  = DSQRT(CM)
! 	      SM  = ( 1.0D0 + MP2*( SST2 - 2.0D0*SST1 ) ) / BM
! 	      ASM = 1.0D0 + SM
! 	      BSM = 1.0D0 - SM
! 	      BM  = 0.5D0*BM*SST11
! 	      PM  =-AM + BM
! 	      QM  =-AM - BM
! 	      PML = PM + SST11
! 	      QML = QM + SST11
! 	    
! 	      IF ( N2.EQ.1 ) THEN
! 		AP  = ASM/(PML*PML)
! 		BQ  = BSM/(QML*QML)
! !*sum for constant CXX
! 		cxx(ro) = cxx(ro) - SUMFA*CONST1*U2*CMP*( AP + BQ )
! 	      ENDIF
! 	    ENDIF
! 
! 	    IF ((CNN.gt.0.0D0).AND.(CM.gt.0.0D0)) THEN
! 	      FAC  = (ASM/PML+BSM/QML)
! 	  !* sum for constant CXX
! 	      cxx(ro) = cxx(ro) - 2.0D0*SUMF22*CONST1*CNP*CMP*FAC*(ANP+BNQ)
! 	  
! 	      ASMAN = ASM*ASN/(PML*(PML+PNA))
! 	      BSMAN = BSM*ASN/(QML*(QML+PNA))
! 	      ASMBN = ASM*BSN/(PML*(PML+QN))
! 	      BSMBN = BSM*BSN/(QML*(QML+QN))
! !* sum for constant CXX
! 	      cxx(ro) = cxx(ro) + 2.0D0*SUMF22*CONST1*CNP*CMP*(ASMAN+BSMAN+ASMBN+BSMBN)
! 	    ENDIF
! !*** CASE4 ***
! 	    IF ( CM.lt.0.D0 ) THEN
! 	      BM  = DSQRT(-CM)
! 	      SM  = 1.0D0+MP2*(SST2-2.0D0*SST1)
! 	      BSM = SM/BM
! 	      BM  = 0.5D0*BM*SST11
! 	      PML =-AM+SST11
! 	      BNM = PML-BM*BSM
! 	      WK0 = PML*PML+BM*BM
! !** single sum Txx
! 	      IF ( N2.EQ.1 ) THEN
! 		ANM = BM+PML*BSM
! 		BQ  =-BM*ANM+PML*BNM
! 		WK1 = 1.0D0/WK0
! 		WK2 = WK1*WK1
! 		FAC = 2.0D0*WK2
! !*  sum for constant CXX
! 		cxx(ro) = cxx(ro) - SUMFA*CONST1*U2*CMP*BQ*FAC
! 	      ENDIF
! 	    ENDIF
! 	    IF ((CNN.lt.0.0D0).AND.(CM.lt.0.0D0)) THEN
! 	      FAC  = 4.0D0*BNM*WKN/WK0
! !* sum for constant CXX
! 	      cxx(ro) = cxx(ro) - 2.0D0*SUMF22*CONST1*CNP*CMP*FAC*BNQ
! 	  
! 	      PMLAN = PML-AN
! 	      QNPM  = BN+BM
! 	      QNQM  = BN-BM
! 	      WK1   = 1.0D0/((PMLAN*PMLAN+QNPM*QNPM)*WK0)
! 	      WK2   = 1.0D0/((PMLAN*PMLAN+QNQM*QNQM)*WK0)
! 	      ANM   = 0.5d0*( BM+PML*BSM+PML*BSN-BM*BSN*BSM)
! 	      BNM   = 0.5d0*(-BM-PML*BSM+PML*BSN-BM*BSN*BSM)
! 	      DNM   = 0.5d0*( PML-BM*BSM-BM*BSN-PML*BSN*BSM)
! 	      ENM   = 0.5d0*( PML-BM*BSM+BM*BSN+PML*BSN*BSM)
! 	      BSMAN =-ANM*QNPM + DNM*PMLAN
! 	      BSMBN =-BNM*QNQM + ENM*PMLAN
! !* sum for constant CXX
! 	      cxx(ro) = cxx(ro)+2.0D0*SUMFA2*CONST1*CNP*CMP*(BSMAN*WK1+BSMBN*WK2)
! 	    ENDIF
! 	  ENDDO  !M
! 	ENDIF !EL1
! 	
!       ENDDO  !N2
! 
! !**   stress components
!       IF (EL1.ne.0.D0) THEN
! 	cxx(ro) = cxx(ro) - 2.0D0*U2*VISC1*( U2*SST1 + cxy(ro) )
! 	cxy(ro) = CONST1*( U2*SST1 + cxy(ro) )
!       ENDIF
!        
!     ENDDO    ! I Loop
!     RETURN
! 
!   END SUBROUTINE WatersIni_Axisymmetric


      

END MODULE waters_solution_module