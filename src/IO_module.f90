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

MODULE IO_module

  USE constants
  USE shared_data
  USE functions_module
  IMPLICIT NONE
  CONTAINS
  
  SUBROUTINE readreport(ierror)
    IMPLICIT NONE
    INTEGER :: ierror
    
    IF (ierror.ne.0) THEN
      print*
      print*,'Error reading file! Value returned was ',ierror
      print*
      STOP
    ENDIF
  END SUBROUTINE readreport

  SUBROUTINE openreport(file,ierror)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ierror
    CHARACTER (LEN=100),INTENT(IN) :: file
    
    IF (ierror.ne.0) THEN
      print*
      print*,'Error opening file ',file,'! Value returned was',ierror
      print*
      STOP
    ENDIF
  END SUBROUTINE openreport  
  
! reads an input file containing the required geometry information.
! this routine will also assign memory for the main geometrical arrays, etc
! using an additional subroutine, assignMem

  SUBROUTINE read_input

  IMPLICIT NONE
  INTEGER :: ierror,i,j,k,l,el,temp,temp1,temp2,lp1,lm1,jp1
  DOUBLE PRECISION :: max_height
                        
  CHARACTER (LEN=64) :: scrap,scrap2
 
! READ INPUT FROM COMMAND LINE:
  CALL parse_command_line_arguments
  
! Once we have N, assign values to a number of constants associated with it:
  NM1=(N-1)
  NP1=(N+1)
  NM1SQ=NM1**2
  NP1SQ=NP1**2
  NP1SQM1=NP1SQ-1
  threeNM1SQ=3*NM1SQ
  
  OPEN(UNIT=10,FILE=input_filename,&
       STATUS='OLD',IOSTAT=ierror)
  REWIND(10)

  CALL openreport(input_filename,ierror)
    

! Skip first 2 lines of file
! will be used later on too:
! each section is seperated by 2 spaces then a header.
  DO i=1,2
    READ(10,'(A)') scrap
  ENDDO
! read in array sizes req'd
  READ(10,'(A8,I6)') scrap,numnp
  READ(10,'(A9,I6)') scrap,numelm
  READ(10,'(A13,D25.18)') scrap,rad_sphere
  READ(10,'(A10,D25.18)') scrap,centre_of_sphere(1)
  READ(10,'(A10,D25.18)') scrap,centre_of_sphere(2)
  READ(10,'(A12,I6)') scrap,coordflag
  READ(10,'(A13,I6)') scrap,preconflag
  READ(10,'(A13,I6)') scrap,movingmeshflag


! Check movingmeshflag and Re make sense:
  IF (Re.lt.1d-6.AND.movingmeshflag.eq.1) THEN
    print*,'ERROR Moving mesh with Re = 0?'
    STOP
  ENDIF
  
  DO i=1,3
    READ(10,'(A)') scrap
  ENDDO

! Allocate memory - any new memory for this stage required should be added to assignMem
  CALL assignMem

! Initialize flags
  bdflag=.false.
  neuflag=.false.
  circflag=.false.
  wallflag=.false.
  inflowflag=.false.
  outflowflag=.false.
  wallsymmflag=.false.
  circnodeflag=.false.
  fixed_node=.false.
  
! Read in node co-ords and log which is the largest value - store as cylinder height, rad_cyl
  max_height=-999d0
  DO
    READ(10,'(A)',IOSTAT=ierror) scrap
    CALL readreport(ierror)
    IF (scrap.eq.'END') EXIT
    READ(scrap,*) el,(vertx(el,k), k=1,2)
    IF (vertx(el,2).gt.max_height) max_height=vertx(el,2)
  ENDDO 
  rad_cylinder = max_height

  DO i=1,3
    READ(10,'(A)') scrap
  ENDDO    

! Read in each element's vertex nodes
  DO    
    READ(10,'(A)',IOSTAT=ierror) scrap
    CALL readreport(ierror)
    IF (scrap.eq.'END') EXIT
    READ(scrap,*) el,(node(el,k), k=1,4)
    node(el,5)=node(el,1)
    node(el,0)=node(el,4)
  ENDDO
  
  DO i=1,3
    READ(10,'(A)') scrap
  ENDDO
  
! Read in boundary nodes
  bdnodes=0
  DO
    READ(10,'(A)',IOSTAT=ierror) scrap
    CALL readreport(ierror)
    IF (scrap.eq.'END') EXIT
    bdnodes=bdnodes+1
    READ(scrap,'(I6)') temp
    bdflag(1,temp)=.true.
    bdflag(2,temp)=.true.
  ENDDO
  
  DO i=1,3
    READ(10,'(A)') scrap
  ENDDO
  
! Read in neumann nodes
  neumann=0
  DO
    READ(10,'(A)',IOSTAT=ierror) scrap
    CALL readreport(ierror)
    IF (scrap.eq.'END') EXIT
    neumann=neumann+1
    READ(scrap,*) temp,temp1,temp2
    IF (temp1.eq.1) THEN
      neuflag(1,temp)=.true.
!       bdflag(1,temp)=.false.
    ENDIF
    IF (temp2.eq.2) THEN
      neuflag(2,temp)=.true.
!       bdflag(2,temp)=.false.
    ENDIF
  ENDDO


  DO i=1,3
    READ(10,'(A)') scrap
  ENDDO
  
! Read in wall nodes
  DO
    READ(10,'(A)',IOSTAT=ierror) scrap
    CALL readreport(ierror)
    IF (scrap.eq.'END') EXIT
    READ(scrap,'(I6)') temp
    wallflag(temp)=.true.
  ENDDO
  
  DO i=1,3
    READ(10,'(A)') scrap
  ENDDO
  
! Read in inflow nodes
  DO
    READ(10,'(A)',IOSTAT=ierror) scrap
    CALL readreport(ierror)
    IF (scrap.eq.'END') EXIT
    READ(scrap,'(I6)') temp
    inflowflag(temp)=.true.
  ENDDO
  
  DO i=1,3
    READ(10,'(A)') scrap
  ENDDO
  
! Read in outflow nodes
  DO
    READ(10,'(A)',IOSTAT=ierror) scrap
    CALL readreport(ierror)
    IF (scrap.eq.'END') EXIT
    READ(scrap,'(I6)') temp
    outflowflag(temp)=.true.
  ENDDO
  
  DO i=1,3
    READ(10,'(A)') scrap
  ENDDO
  
! Read in wall of symmetry nodes
  DO
    READ(10,'(A)',IOSTAT=ierror) scrap
    CALL readreport(ierror)
    IF (scrap.eq.'END') EXIT
    READ(scrap,'(I6)') temp
    wallsymmflag(temp)=.true.
  ENDDO
  
  DO i=1,3
    READ(10,'(A)') scrap
  ENDDO

! Read in circular angle info
  circelm=0
  DO
    READ(10,'(A)',IOSTAT=ierror) scrap
    CALL readreport(ierror)
    IF (scrap.eq.'END') EXIT
    
    circelm=circelm+1
    READ(scrap,*) el,theta(1,el),theta(2,el)
    circflag(el)=.true.
    circnodeflag(node(el,1))=.true.
    circnodeflag(node(el,2))=.true.
    
!Read in theta values in degrees so convert into radians
    theta(1:2,el)=theta(1:2,el)*pi/180d0
    vertx(node(el,1),1)=centre_of_sphere(1) + rad_sphere*cos(theta(1,el))
    vertx(node(el,1),2)=centre_of_sphere(2) + rad_sphere*sin(theta(1,el))  
    vertx(node(el,2),1)=centre_of_sphere(1) + rad_sphere*cos(theta(2,el))
    vertx(node(el,2),2)=centre_of_sphere(2) + rad_sphere*sin(theta(2,el))
  ENDDO

    DO i=1,3
    READ(10,'(A)') scrap
  ENDDO
  
! Read in list of fixed (ie, non-moving) nodes
  DO
    READ(10,'(A)',IOSTAT=ierror) scrap
    CALL readreport(ierror)
    IF (scrap.eq.'END') EXIT
    READ(scrap,'(I6)') temp
    fixed_node(temp)=.true.
  ENDDO
  
  DO i=1,3
    READ(10,'(A)') scrap
  ENDDO
  
! Read in list of accordian (ie, deformable) elements
  accordianflag=.false.
  DO
    READ(10,'(A)',IOSTAT=ierror) scrap
    CALL readreport(ierror)
    IF (scrap.eq.'END') EXIT
    READ(scrap,'(I6)') temp
    accordianflag(temp)=.true.
  ENDDO
  ! should now have all vertices of elements
  CLOSE(10) ! data now read in, close file
  
!Setup the global map and multiplicity with the acquired data
  
!initialize multiplicity to 0s
  DO i=1,numelm*(N+1)**2
    mult(i)=0    
  ENDDO
  
! write known nodes into global map
  mapg=0
  DO i=1,numelm
    mapg(0,i)=node(i,1)
    mapg(N,i)=node(i,2)
    mapg((N+1)**2-1,i)=node(i,3)
    mapg(N*(N+1),i)=node(i,4)
  ENDDO
      
  
! Setup node connectivity, this is used to construct the global map in mapglobal2
  DO i=1,numelm
    DO j=1,4
      jp1=j+1
      IF (jp1.EQ.5) jp1=1
      mult(node(i,j))=1
      conedg(i,j)=0
      conelm(i,j)=0
      DO k=1,numelm
        DO l=1,4
          lp1=l+1
          IF (lp1.EQ.5) lp1=1
          lm1=l-1
          IF (lm1.EQ.0) lm1=4
          IF (i.NE.k) THEN
            IF (node(i,j).EQ.node(k,l)) THEN
              mult(node(i,j))=mult(node(i,j))+1
              IF (node(i,jp1).EQ.node(k,lm1)) THEN
                conedg(i,j)=lm1
                conelm(i,j)=k
              ELSEIF (node(i,jp1).EQ.node(k,lp1)) THEN
                WRITE(*,*) 'ERROR: Elements wrongly connected'
                		print *, i,jp1,k,lp1,node(i,jp1)
                STOP
              ENDIF
            ENDIF
          ENDIF
        ENDDO    
      ENDDO
    ENDDO
  ENDDO

!   do i=1,numelm
!   do j=1,4
!   write(*,*) i,j,node(i,j),' ',conelm(i,j),conedg(i,j),' ',&
!   mult(node(i,j))
!   enddo
!   enddo  
END SUBROUTINE read_input

  SUBROUTINE output_to_tecplot
    IMPLICIT NONE
    INTEGER :: el,ij

! UPDATE TO SET UP FROM STOKES TO NEWTONIAN TO VISCOELASTIC:

    IF (param_problem_choice.eq.11) THEN ! General 2-D Stokes (eg. Poiseuille flow past Cylinder)
    
      write(tecplot_output_fileid,*) 'VARIABLES = "x","y","V_x","V_y","P",', & ! 1 2 3 4 5
			'"GradUxx","GradUyx","GradUxy","GradUyy",', & ! 6 7 8 9 ! components of gradient of velocity
			'"Gxy"' ! 10 ! xy component of rate of strain (= 0.5*(GradU + GradU^(T)))
      DO el=1,numelm
	write(tecplot_output_fileid,*) 'ZONE I =', NP1,', J =', NP1,',F=POINT'
	DO ij=0,NP1SQM1
	  write(tecplot_output_fileid,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),&
		    localV_x(ij,el),localV_y(ij,el),localpressure(ij,el),&
! VELOCITY GRADIENT
		    localGradUxx(ij,el),localGradUyx(ij,el),localGradUxy(ij,el),localGradUyy(ij,el),&
! RATE OF STRAIN:
		    5d-1*(localGradUxy(ij,el) + localGradUyx(ij,el))
	ENDDO
      ENDDO
      
    ELSEIF (param_problem_choice.eq.12) THEN ! General 3-D Axisymm Stokes (eg. Uniform flow past Fixed sphere)
    
      write(tecplot_output_fileid,*) 'VARIABLES = "x","y","V_x","V_y","P",', & ! 1 2 3 4 5
			'"GradUxx","GradUyx","GradUxy","GradUyy","GradUzz",', & ! 6 7 8 9 10 ! components of gradient of velocity
			'"Gxy"' ! 11 ! xy component of rate of strain (= 0.5*(GradU + GradU^(T)))
      DO el=1,numelm
	write(tecplot_output_fileid,*) 'ZONE I =', NP1,', J =', NP1,',F=POINT'
	DO ij=0,NP1SQM1
	  write(tecplot_output_fileid,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),&
		    localV_x(ij,el),localV_y(ij,el),localpressure(ij,el),&
! VELOCITY GRADIENT
		    localGradUxx(ij,el),localGradUyx(ij,el),localGradUxy(ij,el),localGradUyy(ij,el),localGradUzz(ij,el),&
! RATE OF STRAIN:
		    5d-1*(localGradUxy(ij,el) + localGradUyx(ij,el))
	ENDDO
      ENDDO  


    ELSEIF (param_problem_choice.eq.13.or.param_problem_choice.eq.14) THEN ! Stokes 2-D Model solutions.
    
    
      write(tecplot_output_fileid,*) 'VARIABLES = "x","y","V_x","V_y","P",', & ! 1 2 3 4 5
			'"KnownV_x","KnownV_y","KnownPressure",',&
			'"GradUxx","GradUyx","GradUxy","GradUyy",', & ! 6 7 8 9 ! components of gradient of velocity
			'"Gxy",',& ! 10 ! xy component of rate of strain (= 0.5*(GradU + GradU^(T)))
			'"V_xERR","V_yERR","pressureERR",',& ! 11 12 13 ! Error of field variables
			'"GradUxxERR","GradUyxERR","GradUxyERR","GradUyyERR"' ! 14 15 16 17 ! Error of velocity gradient
      DO el=1,numelm
	write(tecplot_output_fileid,*) 'ZONE I =', NP1,', J =', NP1,',F=POINT'
	DO ij=0,NP1SQM1
	  write(tecplot_output_fileid,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),&		    
! 		    V_x(mapg(ij,el)),V_y(mapg(ij,el)),pressure(mapg(ij,el)-npedg),&
		    localV_x(ij,el),localV_y(ij,el),localpressure(ij,el),&
		    V_x_analytic(ij,el),V_y_analytic(ij,el),pressure_analytic(ij,el),&
! VELOCITY GRADIENT
		    localGradUxx(ij,el),localGradUyx(ij,el),localGradUxy(ij,el),localGradUyy(ij,el),&
! RATE OF STRAIN:
		    5d-1*(localGradUyx(ij,el) + localGradUyx(ij,el)),&
! ERRORS:
		    V_x_error(ij,el),V_y_error(ij,el),pressure_error(ij,el),&
		    gradUxx_error(ij,el),gradUyx_error(ij,el),gradUxy_error(ij,el),gradUyy_error(ij,el)

	ENDDO
      ENDDO
    
    ELSEIF (param_problem_choice.eq.21.or.param_problem_choice.eq.41) THEN ! General Newtonian 2-D (Moving and fixed)
    
      write(tecplot_output_fileid,*) 'VARIABLES = "x","y","V_x","V_y","P",', & ! 1 2 3 4 5
			'"GradUxx","GradUyx","GradUxy","GradUyy",', & ! 6 7 8 9 ! components of gradient of velocity
			'"Gxy"' ! 10 ! xy component of rate of strain (= 0.5*(GradU + GradU^(T)))
      DO el=1,numelm
	write(tecplot_output_fileid,*) 'ZONE T="time:',timeN,'" I =', NP1,', J =', NP1,',F=POINT'
	DO ij=0,NP1SQM1
	  write(tecplot_output_fileid,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),&
		    localV_x(ij,el),localV_y(ij,el),localpressure(ij,el),&
! VELOCITY GRADIENT
		    localGradUxx(ij,el),localGradUyx(ij,el),localGradUxy(ij,el),localGradUyy(ij,el),&
! RATE OF STRAIN:
		    5d-1*(localGradUxy(ij,el) + localGradUyx(ij,el))
	ENDDO
      ENDDO
      
    ELSEIF (param_problem_choice.eq.22.or.param_problem_choice.eq.42) THEN ! General Newtonian 3-D Axisymmetric
    
      write(tecplot_output_fileid,*) 'VARIABLES = "x","y","V_x","V_y","P",', & ! 1 2 3 4 5
			'"GradUxx","GradUyx","GradUxy","GradUyy","GradUzz",', & ! 6 7 8 9 10 ! components of gradient of velocity
			'"Gxy"' ! 11 ! xy component of rate of strain (= 0.5*(GradU + GradU^(T)))
      DO el=1,numelm
	write(tecplot_output_fileid,*) 'ZONE T="time:',timeN,'" I =', NP1,', J =', NP1,',F=POINT'
	DO ij=0,NP1SQM1
	  write(tecplot_output_fileid,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),&
		    localV_x(ij,el),localV_y(ij,el),localpressure(ij,el),&
! VELOCITY GRADIENT
		    localGradUxx(ij,el),localGradUyx(ij,el),localGradUxy(ij,el),localGradUyy(ij,el),localGradUzz(ij,el),&
! RATE OF STRAIN:
		    5d-1*(localGradUxy(ij,el) + localGradUyx(ij,el))
	ENDDO
      ENDDO  
    
    ELSEIF (param_problem_choice.eq.23) THEN ! Newtonian known soln (2D)
    
      write(tecplot_output_fileid,*) 'VARIABLES = "x","y","V_x","V_y","P",', & ! 1 2 3 4 5
			'"GradUxx","GradUyx","GradUxy","GradUyy",', & ! 6 7 8 9 ! components of gradient of velocity
			'"Gxy",',& ! 10 ! xy component of rate of strain (= 0.5*(GradU + GradU^(T)))
			'"V_xERR","V_yERR","pressureERR",',& ! 11 12 13 ! Error of field variables
			'"GradUxxERR","GradUyxERR","GradUxyERR","GradUyyERR"' ! 14 15 16 17 ! Error of velocity gradient
      DO el=1,numelm
	write(tecplot_output_fileid,*) 'ZONE T="time:',timeN,'" I =', NP1,', J =', NP1,',F=POINT'
	DO ij=0,NP1SQM1
	  write(tecplot_output_fileid,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),&
		    localV_x(ij,el),localV_y(ij,el),localpressure(ij,el),&
! VELOCITY GRADIENT
		    localGradUxx(ij,el),localGradUyx(ij,el),localGradUxy(ij,el),localGradUyy(ij,el),&
! RATE OF STRAIN:
		    5d-1*(localGradUyx(ij,el) + localGradUyx(ij,el)),&
! ERRORS:
		    V_x_error(ij,el),V_y_error(ij,el),pressure_error(ij,el),&
		    gradUxx_error(ij,el),gradUyx_error(ij,el),gradUxy_error(ij,el),gradUxx_error(ij,el)
	ENDDO
      ENDDO
    
    ELSEIF (param_problem_choice.eq.24) THEN ! Newtonian known soln (3D Axisymmetric)

      write(tecplot_output_fileid,*) 'VARIABLES = "x","y","V_x","V_y","P",', & ! 1 2 3 4 5
			'"GradUxx","GradUyx","GradUxy","GradUyy","GradUzz",', & ! 6 7 8 9 10 ! components of gradient of velocity
			'"Gxy",',& ! 11 ! xy component of rate of strain (= 0.5*(GradU + GradU^(T)))
			'"V_xERR","V_yERR","pressureERR",',& ! 12 13 14 ! Error of field variables
			'"GradUxxERR","GradUyxERR","GradUxyERR","GradUyyERR","GradUzzERR",' ! 15 16 17 18 19 ! Error of velocity gradient
      DO el=1,numelm
	write(tecplot_output_fileid,*) 'ZONE T="time:',timeN,'" I =', NP1,', J =', NP1,',F=POINT'
	DO ij=0,NP1SQM1
	  write(tecplot_output_fileid,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),&
		    localV_x(ij,el),localV_y(ij,el),localpressure(ij,el),&
! VELOCITY GRADIENT
		    localGradUxx(ij,el),localGradUyx(ij,el),localGradUxy(ij,el),localGradUyy(ij,el),localGradUzz(ij,el),&
! RATE OF STRAIN:
		    5d-1*(localGradUyx(ij,el) + localGradUyx(ij,el)),&
! ERRORS:
		    V_x_error(ij,el),V_y_error(ij,el),pressure_error(ij,el),&
		    gradUxx_error(ij,el),gradUyx_error(ij,el),gradUxy_error(ij,el),gradUyy_error(ij,el),gradUzz_error(ij,el)	

	ENDDO
      ENDDO

    ELSEIF (param_problem_choice.eq.31.or.param_problem_choice.eq.51) THEN ! General Viscoelastic 2-D  (Moving and fixed)
    
      write(tecplot_output_fileid,*) 'VARIABLES = "x","y","V_x","V_y","P",', & ! 1 2 3 4 5
			'"GradUxx","GradUyx","GradUxy","GradUyy",', & ! 6 7 8 9 ! components of gradient of velocity
			'"Gxy",',& ! 10 ! xy component of rate of strain (= 0.5*(GradU + GradU^(T)))
			'"Tao_xx","Tao_xy","Tao_yy"' ! 11 12 13 components of elastic stress
      DO el=1,numelm
	write(tecplot_output_fileid,*) 'ZONE T="time:',timeN,'" I =', NP1,', J =', NP1,',F=POINT'
	DO ij=0,NP1SQM1
	  write(tecplot_output_fileid,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),&
		    localV_x(ij,el),localV_y(ij,el),localpressure(ij,el),&
! VELOCITY GRADIENT
		    localGradUxx(ij,el),localGradUyx(ij,el),localGradUxy(ij,el),localGradUyy(ij,el),&
! RATE OF STRAIN:
		    5d-1*(localGradUxy(ij,el) + localGradUyx(ij,el)),&
! ELASTIC STRESS:
		    localTxx(ij,el),localTxy(ij,el),localTyy(ij,el)
	ENDDO
      ENDDO
      
    ELSEIF (param_problem_choice.eq.32.or.param_problem_choice.eq.52) THEN ! General Viscoelastic 3-D Axisymmetric (Moving and fixed)
    
      write(tecplot_output_fileid,*) 'VARIABLES = "x","y","V_x","V_y","P",', & ! 1 2 3 4 5
			'"GradUxx","GradUyx","GradUxy","GradUyy","GradUzz",', & ! 6 7 8 9 10 ! components of gradient of velocity
			'"Gxy",',& ! 11 ! xy component of rate of strain (= 0.5*(GradU + GradU^(T)))
			'"Tao_xx","Tao_xy","Tao_yy","Tao_zz"' ! 12 13 14 15 components of elastic stress
      DO el=1,numelm
	write(tecplot_output_fileid,*) 'ZONE T="time:',timeN,'" I =', NP1,', J =', NP1,',F=POINT'
	DO ij=0,NP1SQM1
	  write(tecplot_output_fileid,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),&
		    localV_x(ij,el),localV_y(ij,el),localpressure(ij,el),&
! VELOCITY GRADIENT
		    localGradUxx(ij,el),localGradUyx(ij,el),localGradUxy(ij,el),localGradUyy(ij,el),localGradUzz(ij,el),&
! RATE OF STRAIN:
		    5d-1*(localGradUxy(ij,el) + localGradUyx(ij,el)),&
! ELASTIC STRESS:
		    localTxx(ij,el),localTxy(ij,el),localTyy(ij,el),localTzz(ij,el)
	ENDDO
      ENDDO


    ELSEIF (param_problem_choice.eq.33) THEN ! Viscoelastic known soln (2D)
    
      write(tecplot_output_fileid,*) 'VARIABLES = "x","y","V_x","V_y","P",', & ! 1 2 3 4 5
			'"GradUxx","GradUyx","GradUxy","GradUyy",', & ! 6 7 8 9 ! components of gradient of velocity
			'"Gxy",',& ! 10 ! xy component of rate of strain (= 0.5*(GradU + GradU^(T)))
			'"Tao_xx","Tao_xy","Tao_yy",',& ! 11 12 13 components of elastic stress
			'"V_xERR","V_yERR","pressureERR",',& ! 14 15 16 ! Error of field variables
			'"GradUxxERR","GradUyxERR","GradUxyERR","GradUyyERR",',& ! 17 18 19 20 ! Error of velocity gradient
			'"Tao_xxERR","Tao_xyERR","Tao_yyERR"' ! 21 22 23 Error of components of elastic stress
      DO el=1,numelm
	write(tecplot_output_fileid,*) 'ZONE T="time:',timeN,'" I =', NP1,', J =', NP1,',F=POINT'
	DO ij=0,NP1SQM1
	  write(tecplot_output_fileid,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),&
		    localV_x(ij,el),localV_y(ij,el),localpressure(ij,el),&
! VELOCITY GRADIENT
		    localGradUxx(ij,el),localGradUyx(ij,el),localGradUxy(ij,el),localGradUyy(ij,el),&
! RATE OF STRAIN:
		    5d-1*(localGradUyx(ij,el) + localGradUyx(ij,el)),&
! ELASTIC STRESS:
		    localTxx(ij,el),localTxy(ij,el),localTyy(ij,el),&
! ERRORS:
		    V_x_error(ij,el),V_y_error(ij,el),pressure_error(ij,el),&
		    gradUxx_error(ij,el),gradUyx_error(ij,el),gradUxy_error(ij,el),gradUxx_error(ij,el),&
		    Txx_error(ij,el),Txy_error(ij,el),Tyy_error(ij,el)
 	ENDDO
      ENDDO
      
    ELSEIF (param_problem_choice.eq.34) THEN ! Viscoelastic known soln (3D Axisymmetric)

      write(tecplot_output_fileid,*) 'VARIABLES = "x","y","V_x","V_y","P",', & ! 1 2 3 4 5
			'"GradUxx","GradUyx","GradUxy","GradUyy","GradUzz",', & ! 6 7 8 9 10 ! components of gradient of velocity
			'"Gxy",',& ! 11 ! xy component of rate of strain (= 0.5*(GradU + GradU^(T)))
			'"Tao_xx","Tao_xy","Tao_yy","Tao_zz",',& ! 11 12 13 14 components of elastic stress
			'"V_xERR","V_yERR","pressureERR",',& ! 15 16 17 ! Error of field variables
			'"GradUxxERR","GradUyxERR","GradUxyERR","GradUyyERR","GradUzzERR",',& ! 18 19 20 21 22 ! Error of velocity gradient
			'"Tao_xxERR","Tao_xyERR","Tao_yyERR","Tao_zzERR"' ! 23 24 25 26 Error of components of elastic stress
      DO el=1,numelm
	write(tecplot_output_fileid,*) 'ZONE T="time:',timeN,'" I =', NP1,', J =', NP1,',F=POINT'
	DO ij=0,NP1SQM1
	  write(tecplot_output_fileid,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),&
		    localV_x(ij,el),localV_y(ij,el),localpressure(ij,el),&
! VELOCITY GRADIENT
		    localGradUxx(ij,el),localGradUyx(ij,el),localGradUxy(ij,el),localGradUyy(ij,el),localGradUzz(ij,el),&
! RATE OF STRAIN:
		    5d-1*(localGradUyx(ij,el) + localGradUyx(ij,el)),&
! ELASTIC STRESS:
		    localTxx(ij,el),localTxy(ij,el),localTyy(ij,el),localTzz(ij,el),&
! ERRORS:
		    V_x_error(ij,el),V_y_error(ij,el),pressure_error(ij,el),&
		    gradUxx_error(ij,el),gradUyx_error(ij,el),gradUxy_error(ij,el),gradUyy_error(ij,el),gradUzz_error(ij,el),&
		    Txx_error(ij,el),Txy_error(ij,el),Tyy_error(ij,el),Tzz_error(ij,el)
	ENDDO
      ENDDO
      
    ENDIF
  
  END SUBROUTINE output_to_tecplot
  
  SUBROUTINE output_nodetypes
    IMPLICIT NONE
    INTEGER :: el,ij
! UNCOMMENT TO PRINT OUT THE NEUMANN AND DIRICHLET NODES TO THE OUT PUT FILE!!
    OPEN(99,FILE='dirichletorneumanngrid.dat')

    write(99,*) 'VARIABLES = "x", "y","node_x","node_y"'
!,"T12","T21","T22"'
    DO el=1,numelm
    write(99,*) 'ZONE I =',NP1,', J =',NP1,',F=POINT'
! ZONE T = timestep |  I = points in y direction | J = points in x direction

!write(20,*) k
      DO ij=0,NP1SQM1
	  IF (bdflag(1,mapg(ij,el)).AND.bdflag(2,mapg(ij,el))) THEN
	    write(99,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),1d0,1d0
	  ELSEIF (neuflag(1,mapg(ij,el)).AND.bdflag(2,mapg(ij,el))) THEN
	    write(99,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),0d0,1d0
	  ELSEIF (neuflag(2,mapg(ij,el)).AND.bdflag(1,mapg(ij,el))) THEN
	    write(99,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),1d0,0d0
	  ELSEIF (neuflag(1,mapg(ij,el)).AND.neuflag(2,mapg(ij,el))) THEN
	    write(99,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),0d0,0d0
	  ELSE
	    write(99,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),-1d0,-1d0
! 	  IF (inflowflag(mapg(ij,el))) THEN
! 	    write(99,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),1d0
! 	  ELSEIF (wallflag(mapg(ij,el))) THEN
! 	    write(99,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),0d0
! 	  ELSEIF (outflowflag(mapg(ij,el))) THEN
! 	    write(99,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),-1d0
! 	  
! 	  ELSE
! 	    write(99,*) nodeCoord(mapg(ij,el),1),nodeCoord(mapg(ij,el),2),0d0
	  ENDIF
      ENDDO
    ENDDO
    CLOSE(99)
  END SUBROUTINE output_nodetypes



  SUBROUTINE output_to_tecplot_finegrid
    IMPLICIT NONE
    INTEGER :: el,ifine,jfine,temp
    
    
  
    write(tecplot_fine_output_fileid,*) 'VARIABLES = "x","y","V_x","V_y","P","GradUxx","GradUyx","GradUxy","GradUyy","GradUzz","Gxy","Tao_xx","Tao_xy","Tao_yy","Tao_zz"'
    DO el=1,numelm
    write(tecplot_fine_output_fileid,*) 'ZONE T="time:',timeN,'" I =', (Nfine + 1),', J =', (Nfine + 1),',F=POINT' 
      DO jfine=0,Nfine
	DO ifine=0,Nfine
	  write(tecplot_fine_output_fileid,*) fine_node_coords(1,ifine,jfine,el),fine_node_coords(2,ifine,jfine,el), &
		      fine_velocity_x(ifine,jfine,el),fine_velocity_y(ifine,jfine,el),fine_pressure(ifine,jfine,el),&   
		      fineGxx(ifine,jfine,el),fineGyx(ifine,jfine,el),fineGxy(ifine,jfine,el),fineGyy(ifine,jfine,el),fineGzz(ifine,jfine,el),&
		      5d-1*(fineGyx(ifine,jfine,el) + fineGxy(ifine,jfine,el)),&
		      fine_tao_xx(ifine,jfine,el),fine_tao_xy(ifine,jfine,el),fine_tao_yy(ifine,jfine,el),fine_tao_zz(ifine,jfine,el)

	ENDDO
      ENDDO
    ENDDO
    
  
  END SUBROUTINE output_to_tecplot_finegrid
  

  
  
  SUBROUTINE create_fine_solution!(Nf,grid_in,sol_outx,sol_outy)
! Interpolates the velocity onto a fine grid, element-wise
    IMPLICIT NONE
    INTEGER :: el,i,j,ij,jj,ifine,jfine,jjf,ijfine,globalij,intij
    DOUBLE PRECISION :: temp,tempx,tempy,tempxx,tempxy,tempyy,tempzz,tempp,&
			tempGxx,tempGxy,tempGyx,tempGyy,tempGzz
    
    DO el=1,numelm
      DO jfine=0,Nfine
	jjf=jfine*(Nfine+1)
	DO ifine=0,Nfine
	  ijfine=ifine+jjf
	  tempx=0d0
	  tempy=0d0
	  tempxx=0d0
	  tempxy=0d0
	  tempyy=0d0
	  tempzz=0d0
	  tempGxx=0d0
	  tempGxy=0d0
	  tempGyx=0d0
	  tempGyy=0d0
	  tempGzz=0d0
	  DO ij=0,NP1SQM1
! 	    globalij=mapg(ij,el)
	    temp = fine_velocitybasis(ij,ijfine)!fine_hbasis(i,ifine)*fine_hbasis(j,jfine)
	    tempx = tempx + localV_x(ij,el)*temp
	    tempy = tempy + localV_y(ij,el)*temp
	    tempxx = tempxx + localTxx(ij,el)*temp
	    tempxy = tempxy + localTxy(ij,el)*temp
	    tempyy = tempyy + localTyy(ij,el)*temp
	    tempzz = tempzz + localTzz(ij,el)*temp
	    tempGxx = tempGxx + localGradUxx(ij,el)*temp
	    tempGxy = tempGxy + localGradUxy(ij,el)*temp
	    tempGyx = tempGyx + localGradUyx(ij,el)*temp
	    tempGyy = tempGyy + localGradUyy(ij,el)*temp
	    tempGzz = tempGzz + localGradUzz(ij,el)*temp
	  ENDDO
	  fine_velocity_x(ifine,jfine,el) = tempx
	  fine_velocity_y(ifine,jfine,el) = tempy
	  fine_tao_xx(ifine,jfine,el) = tempxx
	  fine_tao_xy(ifine,jfine,el) = tempxy
	  fine_tao_yy(ifine,jfine,el) = tempyy
	  fine_tao_zz(ifine,jfine,el) = tempzz
	  fineGxx(ifine,jfine,el) = tempGxx
	  fineGxy(ifine,jfine,el) = tempGxy
	  fineGyx(ifine,jfine,el) = tempGyx
	  fineGyy(ifine,jfine,el) = tempGyy
	  fineGzz(ifine,jfine,el) = tempGzz
	  
	  tempp=0d0
	  DO intij=1,NM1SQ
	    ij=interior_to_local_node(intij)
	    tempp = tempp + localpressure(ij,el)*fine_pressurebasis(intij,ijfine)!fine_hbasis_tilde(i,ifine)*fine_hbasis_tilde(j,jfine)
	  ENDDO
	  fine_pressure(ifine,jfine,el) = tempp
	ENDDO
      ENDDO
    ENDDO
  
  END SUBROUTINE create_fine_solution
  
  
  SUBROUTINE output_along_wallsymm
    IMPLICIT NONE
    INTEGER :: el,edge,ii,jj
    
    DO el=1,numelm
      DO edge=1,4
	IF (is_wallsymm_edge(edge,el).or.is_circ_edge(edge,el)) THEN
	  IF (edge.eq.1) THEN
	    jj=0
	    DO ii=0,Nfine
	      write(wallsymm_output_fileid,*) timeN,fine_node_coords(1,ii,jj,el),& ! 1 2
			      fine_velocity_x(ii,jj,el),fine_velocity_y(ii,jj,el),fine_pressure(ii,jj,el),& ! 3 4 5
			      fine_tao_xx(ii,jj,el),fine_tao_xy(ii,jj,el),fine_tao_yy(ii,jj,el),fine_tao_zz(ii,jj,el),& ! 6 7 8 9
			      fineGxx(ii,jj,el),fineGxy(ii,jj,el),fineGyx(ii,jj,el),fineGyy(ii,jj,el),fineGzz(ii,jj,el) ! 10 11 12 13 14
	    ENDDO
	  ELSEIF (edge.eq.2) THEN
	    ii=Nfine
	    DO jj=0,Nfine
	      write(wallsymm_output_fileid,*) timeN,fine_node_coords(1,ii,jj,el),& ! 1 2
			      fine_velocity_x(ii,jj,el),fine_velocity_y(ii,jj,el),fine_pressure(ii,jj,el),& ! 3 4 5
			      fine_tao_xx(ii,jj,el),fine_tao_xy(ii,jj,el),fine_tao_yy(ii,jj,el),fine_tao_zz(ii,jj,el),& ! 6 7 8 9
			      fineGxx(ii,jj,el),fineGyx(ii,jj,el),fineGxy(ii,jj,el),fineGyy(ii,jj,el),fineGzz(ii,jj,el) ! 10 11 12 13 14
	    ENDDO
	  ELSEIF (edge.eq.3) THEN
	    jj=Nfine
	    DO ii=0,Nfine
	      write(wallsymm_output_fileid,*) timeN,fine_node_coords(1,ii,jj,el),& ! 1 2
			      fine_velocity_x(ii,jj,el),fine_velocity_y(ii,jj,el),fine_pressure(ii,jj,el),& ! 3 4 5
			      fine_tao_xx(ii,jj,el),fine_tao_xy(ii,jj,el),fine_tao_yy(ii,jj,el),fine_tao_zz(ii,jj,el),& ! 6 7 8 9
			      fineGxx(ii,jj,el),fineGxy(ii,jj,el),fineGyx(ii,jj,el),fineGyy(ii,jj,el),fineGzz(ii,jj,el) ! 10 11 12 13 14
	    ENDDO
	  ELSEIF (edge.eq.4) THEN
	    ii=0
	    DO jj=0,Nfine
	      write(wallsymm_output_fileid,*) timeN,fine_node_coords(1,ii,jj,el),& ! 1 2
			      fine_velocity_x(ii,jj,el),fine_velocity_y(ii,jj,el),fine_pressure(ii,jj,el),& ! 3 4 5
			      fine_tao_xx(ii,jj,el),fine_tao_xy(ii,jj,el),fine_tao_yy(ii,jj,el),fine_tao_zz(ii,jj,el),& ! 6 7 8 9
			      fineGxx(ii,jj,el),fineGxy(ii,jj,el),fineGyx(ii,jj,el),fineGyy(ii,jj,el),fineGzz(ii,jj,el) ! 10 11 12 13 14
	    ENDDO
	  ENDIF
	ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE output_along_wallsymm
  
!!! NEEDS UPDATING !!!
  SUBROUTINE output_to_screen
    IMPLICIT NONE
    INTEGER :: globi

! STOKES:
    IF (param_problem_choice.eq.11) THEN
! Stokes Poiseuille flow past cylinder.
! Want numelm, N, drag, cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim
	WRITE(*,*) numelm, N, drag, & 
		    cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim
    ELSEIF (param_problem_choice.eq.12) THEN
! Uniform Stokes flow past sphere.
! Want numelm, N, drag, dragstar, stopping_criteria, cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim
	WRITE(*,*) numelm, N, drag, drag_star, & 
		    cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim
		    
    ELSEIF (param_problem_choice.eq.13) THEN
! Model problem.
! Want numelm, N, H1norm(velocity), L2norm(pressure), cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim   
	WRITE(*,*) numelm, N, L2norm_vel_err, H1norm_vel_err, L2norm_press_err, & 
		    cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim

    ELSEIF (param_problem_choice.eq.14) THEN
! Cylinder known soln.
! Want numelm, N, drag, H1norm(velocity), L2norm(pressure), cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim
	WRITE(*,*) numelm, N, drag, L2norm_vel_err, H1norm_vel_err, L2norm_press_err, & 
		    cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim

! GENERAL NEWTONIAN/VISCOELASTIC:
    ELSEIF (param_problem_choice.eq.21.or.param_problem_choice.eq.31) THEN
! General 2-D output (usually poiseuille flow past cylinder).
! Want numelm, N, time, drag, drag_diff, stopping_criteria, cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim
	WRITE(*,*) numelm, N, timeN, drag, (drag-dragNm1)/drag, stopping_criteria, & 
		    cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim  
		    
    ELSEIF (param_problem_choice.eq.22.or.param_problem_choice.eq.32) THEN
! General 3-D Axisymmetric output (usually Uniform flow past sphere).
! Want numelm, N, time, drag, dragstar, drag_diff, stopping_criteria, cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim
	WRITE(*,*) numelm, N, timeN, drag, drag_star, (drag-dragNm1)/drag, stopping_criteria, & 
		    cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim

    ELSEIF (param_problem_choice.eq.41.or.param_problem_choice.eq.51) THEN
! Moving mesh, 2-D
! Want numelm, N, time, drag, V_sphere, centre of sphere, V_sphere diff, drag diff, stopping_criteria, cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim
      WRITE(*,*) numelm, N, timeN, drag, &
		  V_sphere, centre_of_sphere(1), (V_sphere-V_sphereNm1)/V_sphere,(drag-dragNm1)/drag, stopping_criteria, & 
		    cputime_initialise,cputime_setup,cputime_solve,cputime_total, global_dim
		    
    ELSEIF (param_problem_choice.eq.42.or.param_problem_choice.eq.52) THEN
! Moving mesh, 3-D Axisymmetric
! Want numelm, N, time, drag, drag_star, V_sphere, centre of sphere, V_sphere diff, drag diff, stopping_criteria, cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim
      WRITE(*,*) numelm, N, timeN, drag, drag_star, & ! 1 2 3 4 5
		  V_sphere, centre_of_sphere(1), (V_sphere-V_sphereNm1)/V_sphere,(drag-dragNm1)/drag, stopping_criteria, & ! 6 7 8 9 10
		    cputime_initialise,cputime_setup,cputime_solve,cputime_total, global_dim ! 11 12 13 14 15
		    
! NEWTONIAN TRANSIENT KNOWN:
    ELSEIF (param_problem_choice.eq.23.or.param_problem_choice.eq.24) THEN
      globi = mapg(transient_stress_testnode,transient_stress_testelement)
! NEWTONIAN TRANSIENT POISEUILLE 2-D / 3-D
! Want numelm, N, time, transient_U_error, transient_Txy_error, stopping_criteria, cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim
    	WRITE(*,*) numelm, N, timeN, & ! 1 2 3 
		V_x(transient_velocity_testnode), transient_u(transient_velocity_testnode), & ! 4 5
		localGradUyx(transient_stress_testnode,transient_stress_testelement),& ! 6
		transient_gradUyx(globi), & ! 7
		L2norm_vel_err, H1norm_vel_err, L2norm_press_err, transient_U_error, transient_gradUyx_error, & ! 8 9 10 11 12
		cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim

! VISCOELASTIC TRANSIENT KNOWN:
    ELSEIF (param_problem_choice.eq.33.or.param_problem_choice.eq.34) THEN
      globi = mapg(transient_stress_testnode,transient_stress_testelement)
! VISCOELASTIC TRANSIENT POISEUILLE 2-D / 3-D
! Want numelm, N, time, transient_U_error, transient_Txx_error, transient_Txy_error, stopping_criteria, cputime_initialise, cputime_setup, cputime_solve, cputime_total
    	WRITE(*,*) numelm, N, timeN, & ! 1 2 3
		V_x(transient_velocity_testnode), transient_u(transient_velocity_testnode), & ! 4 5
		localTxx(transient_stress_testnode,transient_stress_testelement),& ! 6 
		transient_txx(globi), & ! 7
		localTxy(transient_stress_testnode,transient_stress_testelement),& ! 8
		transient_txy(globi), & ! 9
		localGradUyx(transient_stress_testnode,transient_stress_testelement),& ! 10
		transient_gradUyx(globi), & ! 11
		L2norm_vel_err, H1norm_vel_err, L2norm_press_err, L2norm_stress_err, & ! 12 13 14 15
		transient_U_error, transient_gradUyx_error, transient_Txx_error, transient_Txy_error, & ! 16 17 18 19 20
		cputime_initialise, cputime_setup, cputime_solve, cputime_total, global_dim

    ENDIF
  END SUBROUTINE output_to_screen
  
  SUBROUTINE open_output_files
    IMPLICIT NONE
    INTEGER :: ierror

    OPEN(tecplot_output_fileid,FILE=tecplot_output_filename,IOSTAT=ierror)
    CALL openreport(tecplot_output_filename,ierror)
    OPEN(wallsymm_output_fileid,FILE=wallsymm_output_filename,IOSTAT=ierror)
    CALL openreport(wallsymm_output_filename,ierror)

  END SUBROUTINE open_output_files
  
  SUBROUTINE final_stage_output
! Prints out finegrid and wall of symmetry values.
! Must have generated the finegrid variables first.
    IMPLICIT NONE
    INTEGER :: ierror
    
! Print out final solution to both tecplot and tecplotfine grids. Only need the non-fine version if we've not got a fine grid
    IF (movingmeshflag.eq.0) THEN
      OPEN(tecplot_output_fileid,FILE=tecplot_output_filename,IOSTAT=ierror)
      CALL openreport(tecplot_output_filename,ierror)
      CALL output_to_tecplot
    ENDIF
    
    OPEN(tecplot_fine_output_fileid,FILE=tecplot_fine_output_filename,IOSTAT=ierror)
    CALL openreport(tecplot_fine_output_filename,ierror)
    CALL output_to_tecplot_finegrid
    CLOSE(tecplot_fine_output_fileid)
    
    IF ( MOD(param_problem_choice,10).eq.1.or.MOD(param_problem_choice,10).eq.2 ) THEN
 ! Only required for the general flow past cylinder or sphere or moving mesh
      OPEN(wallsymm_output_fileid,FILE=wallsymm_output_filename,IOSTAT=ierror)
      CALL openreport(wallsymm_output_filename,ierror)
      CALL output_along_wallsymm
      CLOSE(wallsymm_output_fileid)
      
    ENDIF
    
    
  END SUBROUTINE final_stage_output
  
  
  SUBROUTINE parse_command_line_arguments
! READ IN ALL COMMAND LINE ARGUMENTS, DEFAULTS ARE SET IN SHARED_DATA.
    IMPLICIT NONE
    CHARACTER (LEN=256) :: arg_in,temp
    INTEGER :: i,temp_int
    
    i=1
    DO 
      CALL GETARG(i,temp)
      IF (len_trim(temp) == 0) EXIT
      arg_in=trim(temp)
      IF (arg_in.eq.'-help') THEN
	print*, 'Usage, ./stokes [arg list]'
	print*, 'Available commandline arguments:'
	print*, ''
	print*, '-input [filename]  | input filename, may include path and must include suffix. '
	print*, '-N [integer]       | Order of spectral approximation'
	print*, '-Re [double]       | Reynolds number (default: 0.0 / OFF)'
	print*, '-We [double]       | Weissenberg number (default: 0.0 / OFF)'
	print*, '-beta [double]     | Viscosity Ratio (default: 1.0 / OFF)' 
	print*, '-giesekus [double] | Giesekus slip parameter, alpha (default: 0.0 / OFF)'
	print*, '-beta_s [double]   | DEVSS artificial viscosity parameter (default: 0.0 / OFF)'
	print*, '-prob [integer]    | problem choice - use -prob_help for more info'
	print*, '-output [filename] | output leading filename (3 will be produced), do not include suffix'
	print*, '-deltat [double]   | timestep size (default: 0.001)'
	print*, '-alphaz [double]   | Co-efficient of integral of pressure in conservation of mass equation (default: 0.0 / OFF)'
	print*, '-con_it [double]   | Choose the treatment of the convective term in the constituitive equation. 0 = EX2, 1=semi-iterative (default: 1 / semi-iterative)'
	print*, '-time_order [integer]   | Choose the time-order of the numerical scheme. Available options: 1=1st-order, 2=2nd-order (default), 3=3rd-order'
	print*, '-rho_f [double]   | density of fluid (default: 0.868)'
	print*, '-rho_s [double]   | density of sphere (default: 3.581)'
	STOP
      ELSEIF (arg_in.eq.'-prob_help'.or.arg_in.eq.'-probhelp') THEN
            
	print*,'parameter to set the boundary conditions and result analysis'
	print*,'GUIDE:'
	print*,''	
	print*,'1x = Stokes (Fixed mesh)'
	print*,'2x = Newtonian (Fixed mesh)'
	print*,'3x = Viscoelastic (Fixed mesh)'
	print*,'4x = Newtonian (Moving mesh)'
	print*,'5x = Viscoelastic (Moving mesh)'
	print*,''
	print*,'x1 = Poiseuille flow (past cylinder)'
	print*,'x2 = Uniform flow (past sphere)'
	print*,'x3 = Stokes: model solution | Newt/Visco: transient Poiseuille (2-D)'
	print*,'x4 = Stokes: known cylinder solution | Newt/Visco: transient Poiseuille (Axisymm 3-D)'
	STOP
      ENDIF
      i=i+1
      CALL GETARG(i,temp)
      IF (len_trim(temp) == 0) THEN
	print*,'ERROR: Mismatch in command line arguments: ',arg_in,' and ',trim(temp)
	STOP
      ENDIF
      IF (arg_in.eq.'-input') THEN
	input_filename=trim(temp)
	
      ELSEIF (arg_in.eq.'-n'.or.arg_in.eq.'-N') THEN
	READ(temp,*) N
	
      ELSEIF (arg_in.eq.'-Re'.or.arg_in.eq.'-re') THEN
	READ(temp,*) Re
	
      ELSEIF (arg_in.eq.'-We'.or.arg_in.eq.'-we') THEN
	READ(temp,*) We
	
      ELSEIF (arg_in.eq.'-beta'.or.arg_in.eq.'-Beta') THEN
	READ(temp,*) param_beta
	
      ELSEIF (arg_in.eq.'-giesekus'.or.arg_in.eq.'-Giesekus') THEN
	READ(temp,*) param_giesekus
	
      ELSEIF (arg_in.eq.'-beta_s'.or.arg_in.eq.'-Beta_s') THEN
	READ(temp,*) param_beta_s
	
      ELSEIF (arg_in.eq.'-prob'.or.arg_in.eq.'-Prob') THEN
	READ(temp,*) param_problem_choice
	
      ELSEIF (arg_in.eq.'-output') THEN      
	output_filename=trim(temp)
	tecplot_output_filename = trim(output_filename)//'_tec.dat'
	tecplot_fine_output_filename = trim(output_filename)//'_tec_fine.dat'
	wallsymm_output_filename = trim(output_filename)//'_wallsymm.txt'
	
      ELSEIF (arg_in.eq.'-deltat'.or.arg_in.eq.'-Deltat') THEN      
	READ(temp,*) deltat
      ELSEIF (arg_in.eq.'-alphaz'.or.arg_in.eq.'-alphaZ' &
	      .or.arg_in.eq.'-Alphaz'.or.arg_in.eq.'-AlphaZ' ) THEN
	READ(temp,*) param_alphaZ
      ELSEIF (arg_in.eq.'-cons_it'.or.arg_in.eq.'-Cons_it') THEN
	READ(temp,*) temp_int
	IF (temp_int.eq.1) THEN
	  param_iterative_convection=.true.
	ELSEIF (temp_int.eq.0) THEN
	  param_iterative_convection=.false.
	ELSE
	  print*,'ERROR: Input argument cons_it invalid: use values 0 (off) or 1 (on)'
	  STOP
	ENDIF
      ELSEIF (arg_in.eq.'-time_order'.or.arg_in.eq.'-Time_order') THEN
      	READ(temp,*) param_time_order
      	IF (param_time_order.gt.2.or.param_time_order.lt.1) THEN
	  print*,'ERROR: Input argument time_order invalid: use values 1, 2'! or 3' REMOVED FOR NOW.
	  STOP
	ENDIF
      ELSEIF (arg_in.eq.'-rho_s') THEN
	READ(temp,*) rho_s
      ELSEIF (arg_in.eq.'-rho_f') THEN
	READ(temp,*) rho_f
      ELSE
	print*,'ERROR: Unknown commandline arguments: ',arg_in,' and ',trim(temp)
	print*,'Try using the argument -help'
	STOP
      ENDIF
      i=i+1
      ENDDO

  END SUBROUTINE parse_command_line_arguments  
  
  SUBROUTINE create_cross_stream_points
    IMPLICIT NONE
    INTEGER :: i,j,el,k
    DOUBLE PRECISION :: fixed_x, y_coords(1:25),temp_x,temp_y,temp_xi,temp_eta,
  
    fixed_x=2d0
    DO i=1:25
      temp_x=fixed_x
      temp_y=y_coords(i)
! Find local element and xi/eta co-ordinate from x/y co-ord :
      CALL map_x_y_to_xi_eta(temp_x, temp_y, el_out, temp_xi, temp_eta)
      
            
      
    ENDDO
    
    
  
  END SUBROUTINE create_cross_stream_points
  
  
  
END MODULE IO_module
