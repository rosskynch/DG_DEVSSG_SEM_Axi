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

MODULE shared_data
!
! This module contains the first set of shared data and storage arrays which
! can be thrown between various modules/programs
!
! Explanation of shared data:

! CHARACTER DATA:
! glgrid_output_filename - name for Gauss-Lobatto grid point output file
!                          the default is glgrid.dat

! INTEGER DATA:
! N - the number of GL points
! numelm - the number of elements
! nptot - total number of nodes
! numnp - number of input nodes (ie, 1 element has 4)
! npedg - total number of boundary nodes
! npint - total number of internal nodes
! bdnodes - number of input boundary nodes
! circelm - number of circular elements
! conelm & conedg - connectivity of elements, elm tells which element the edge is touching
!                   and edg tells which edge of that element it is touching
! mapg - the global map
! locel -
! locnp - not clear on what these 2 are - do I need??
! mult - multiplicity of all nodes (including GL)
!
! LOGICAL DATA:
! bdflag - flag for each node declaring if dirichlet boundary or not
! neuflag - flag for each node declaring if neumann boundary or not
! circflag - flag for each ELEMENT declaring if the element has a circular boundary
!		NOTE: if = .true. => node 1 and node 2 are the circular nodes
! 
! GLOBAL node flags in priority order
! (eg, if a node is between a wall vertex and an inflow vertex. it will be flagged as wall!)
! wallflag - fixed wall boundary node - u=0
! inflowflag - inflow boundary node
! outflowflag - outflow boundary node
! wallsymmflag - boundary node on the line of symmetry (cylindrical co-ords)
! circnodeflag - boundary node on a circular boundary
!
!
! DOUBLE PRECISION DATA:
!
! rad_sphere - the radius of the sphere considered in the circular geometry case
! Re - Reynold's number
! sphereXleft, et al - max/min X/Y values of co-ords for sphere. Used for boundary purposes
! inflowX and outflowX,- X co-ord of the inflow and outflow boundaries.
! wallYtop/bottom - Y co-ord of the wall boundaries
!
! centre_of_sphere - X and Y co-ords of sphere's centre
!
! vertx - the co-ords of the vertex nodes <-- subject to change
! gl - the G-L points
! w - the G-L weights
! Leg - values of the legendre polynomials at gl points
! Leg1 - values of the differential of the legendre polynomials at gl points
! d - the differential operator matrix
! dxdp,dxde,dydp,dyde - all matrices storing these values for each element
! jac  - matrix storing value of det(J_k), for each gl point in each element k
! A_r & A_z - Laplacian operator (in respective co-ord) matrix for each element
! f_x & f_y - load (vector) at each node
! u - velocity (assumed scalar for now) at each node
! above will soon be ux, uy 
! pressure - vector of pressures at each node
!
!
!
! have added jac - dimensions (1:numelm,0:N,0:N) - stores jacobian matrix for each element
!
! enable_output_to_file - controls if output to file is allowed from the main program.

  IMPLICIT NONE
  
  LOGICAL :: enable_output_to_file = .true.
  
  CHARACTER(LEN=256) ::	output_filename='output',&
			tecplot_output_filename = 'output_tec.dat',&
			tecplot_fine_output_filename = 'output_tec_fine.dat',&
			input_filename = 'input.dat',&
			wallsymm_output_filename = 'output_wallsymm.txt'
  
  INTEGER ::	N=8, NM1, NM1SQ, threeNM1SQ, NP1, NP1SQ, NP1SQM1, &
		tecplot_output_fileid=30, wallsymm_output_fileid=40, tecplot_fine_output_fileid=50, &
		numelm,&
		nptot,&
		numnp,&
		npedg,&
		npint,&
		bdnodes,&
		neumann,&
		circelm,&
		preconflag=0,& !default = no precon
		coordflag=0,&    !default = cartesian co-ords
		movingmeshflag=0,& ! allow mesh to move? - 1 = yes
		global_dim,& !dimension of square matrix when dirichlet nodes are removed = bd_numrows_x+bd_numrows_y+3*npint
		global_non_zeros=0,& ! number of non-zeros in global matrix (for pardiso)
		bd_numrows_x,&
		bd_numrows_y,&
		global_bd_dim,&
		RK4_timesteps,&
		Nfine = 20,&
		num_wallsymm_nodes,&
		param_problem_choice=0,&
		param_function_choice=0,&
		param_time_order=2,&
		transient_velocity_testnode,&
		transient_stress_testnode,&
		transient_stress_testelement

  INTEGER, ALLOCATABLE, DIMENSION(:) :: mult, &
					non_dir_bd_map_x, &
					non_dir_bd_map_y

  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: conelm,&
                                          conedg,&
                                          node,&
                                          mapg,&
                                          mapg_pressure!,&
!                                           locel,&
!                                           locnp
                                          
  INTEGER, ALLOCATABLE :: upwind_local_edge_node(:,:,:),&
			  local_edge_node(:,:),&
			  global_to_local_map(:,:),&
			  local_to_interior_node(:),&
			  interior_to_local_node(:),&
			  local_ij_to_edge_node(:,:),&
			  local_ij_to_i(:),&
			  local_ij_to_j(:),&
			  upwind_wallsymm_1d(:,:,:),&
			  fine_edge_ifine(:,:),&
			  fine_edge_jfine(:,:),&
			  axisymm_edge(:)
                                          
					  
  LOGICAL :: param_waters=.false.,&
	      param_error=.false.,&
	      param_iterative_convection=.true.
  
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: circflag,&
  					wallflag,&
					inflowflag,&
					outflowflag,&
					wallsymmflag,&
					circnodeflag,&
					fixed_node,&
					accordianflag

  
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: bdflag,&
					  neuflag
					  
  LOGICAL, ALLOCATABLE :: is_wallsymm_edge(:,:),&
			  is_circ_edge(:,:),&
			  is_wall_edge(:,:),&
			  is_inflow_edge(:,:),&
			  is_outflow_edge(:,:)
			  

					  
  DOUBLE PRECISION :: rad_sphere,rad_cylinder,&
		      Re=0d0,&
		      We=0d0,&
		      param_giesekus=0d0,&
		      param_beta=1d0,&
		      param_beta_s=0d0,&
! 		      param_delta_a,&
		      param_alphaZ=0d0,&
		      area_of_domain,&
! 		      param_elasticity,&
! 		      param_lambda,&
		      gravity_const=9.8d0,&
		      rho_s=3.581d0,& ! sphere density
		      rho_f=0.868d0,& ! fluid density
		      sphereXleft,&
		      sphereXright,&
		      sphereYtop,&
		      sphereYbottom,&
		      inflowX,& ! could also be a fixed wall
		      outflowX,&
		      wallYtop,&
		      wallYbottom,& ! only needed if we don't have a sphere in the geometry
		      deltat=1d-3,&
		      timeN,& ! time NOW
		      timeNm1,& ! time LAST 1
		      timeNm2,& ! time LAST 2
		      drag=0d0,&
		      dragNm1=0d0,&
		      dragNm2=0d0,&
		      drag_star=0d0,&
		      H1norm_vel_err,&
		      L2norm_vel_err,&
		      L2norm_press_err,&
		      L2norm_stress_err,&
		      transient_U_error, transient_gradUyx_error, transient_Txx_error, transient_Txy_error,&
		      stopping_criteria,&
		      V_sphere=0d0,&
		      V_sphereNm1=0d0,&
		      V_sphereNm2=0d0,&
		      h_RK4,& !stepsize for RK4 solver
! constants for OIFS/BDF/EX schemes
		      time_gamma_0=0d0,&
		      time_alpha_0=0d0,&
		      time_alpha_1=0d0,&
		      time_alpha_2=0d0,&
		      time_beta_0=0d0,&
		      time_beta_1=0d0,&
		      time_beta_2=0d0,&
		      Retime_constant1=0d0,&
		      Retime_constant2=0d0,&
		      Wetime_constant1=0d0,&
		      Wetime_constant2=0d0,&
		      cputime_initialise,cputime_setup,cputime_solve,cputime_total,cputime_newtonian,cputime_viscoelastic

! NEW BITS
  DOUBLE PRECISION, ALLOCATABLE :: storeA_x(:,:,:),& !param_beta_a(:),&  
				   storeA_y(:,:,:),&				   
				   storeB_x(:,:,:),&
				   storeB_y(:,:,:),&
				   storeZ_p(:,:),&
				   storef_x(:),&
				   storef_y(:),&
				   storeg_p(:),&
! 				   storeM_stress(:,:),&
				   storeM_pressure(:,:,:),&
				   storeMv_x(:,:,:),&
				   storeMv_y(:,:,:),&
! 				   storeC_x(:,:,:),&
! 				   storeC_y(:,:,:),&
! 				   storeCb(:,:,:,:),&
! 				   M_stress(:,:),&
				   Mv_x(:,:,:),&
				   Mv_y(:,:,:),&
! 				   C_x(:,:,:),&
! 				   C_y(:,:,:),&
! 				   Cb(:,:,:,:),&
				   localTxx(:,:),&
				   localTxy(:,:),&
				   localTyy(:,:),&
				   localTzz(:,:),&
				   localTxxNm1(:,:),&
				   localTxyNm1(:,:),&
				   localTyyNm1(:,:),&
				   localTzzNm1(:,:),&
				   localTxxNm2(:,:),&
				   localTxyNm2(:,:),&
				   localTyyNm2(:,:),&
				   localTzzNm2(:,:),&
				   localConformationxx(:,:),&
				   localConformationxy(:,:),&
				   localConformationyy(:,:),&
				   localConformationzz(:,:),&
				   localConformationxxNm1(:,:),&
				   localConformationxyNm1(:,:),&
				   localConformationyyNm1(:,:),&
				   localConformationzzNm1(:,:),&
				   localConformationxxNm2(:,:),&
				   localConformationxyNm2(:,:),&
				   localConformationyyNm2(:,:),&
				   localConformationzzNm2(:,:),&
				   norm_to_edge_node(:,:,:,:),&
				   jac_on_edge(:,:),&
				   jac_wallsymm_edge_1d(:,:),&
				   localGradUxx(:,:),&
				   localGradUxy(:,:),&
				   localGradUyx(:,:),&
				   localGradUyy(:,:),&
				   localGradUzz(:,:),&
				   localGradUxxNm1(:,:),&
				   localGradUxyNm1(:,:),&
				   localGradUyxNm1(:,:),&
				   localGradUyyNm1(:,:),&
				   localGradUzzNm1(:,:),&
				   localGradUxxNm2(:,:),&
				   localGradUxyNm2(:,:),&
				   localGradUyxNm2(:,:),&
				   localGradUyyNm2(:,:),&
				   localGradUzzNm2(:,:),&
				   diff_x(:,:,:),&
				   diff_y(:,:,:),&
				   fine_uniform_points(:),&
				   fine_node_coords(:,:,:,:),&
				   fine_hbasis(:,:),&
				   fine_hbasis_tilde(:,:),&
				   fine_velocitybasis(:,:),&
				   fine_pressurebasis(:,:),&
				   fine_velocity_x(:,:,:),&
				   fine_velocity_y(:,:,:),&
				   fine_pressure(:,:,:),&
				   fine_tao_xx(:,:,:),&
				   fine_tao_xy(:,:,:),&
				   fine_tao_yy(:,:,:),&
				   fine_tao_zz(:,:,:),&
				   fineGxx(:,:,:),&
				   fineGxy(:,:,:),&
				   fineGyx(:,:,:),&
				   fineGyy(:,:,:),&
				   fineGzz(:,:,:),&
				   localpressure(:,:),&
				   localpressureNm1(:,:),&
				   localpressureNm2(:,:),&
				   localV_x(:,:),&
				   localV_y(:,:),&
				   localV_xNm1(:,:),&
				   localV_yNm1(:,:),&
				   localV_xNm2(:,:),&
				   localV_yNm2(:,:),&
! analytical storage:
				   V_x_analytic(:,:),&
				   V_y_analytic(:,:),&
				   gradUxx_analytic(:,:),&
				   gradUyx_analytic(:,:),&
				   gradUxy_analytic(:,:),&
				   gradUyy_analytic(:,:),&
				   pressure_analytic(:,:),&
! error storage:
				   V_x_error(:,:),&
				   V_y_error(:,:),&
				   gradUxx_error(:,:),&
				   gradUyx_error(:,:),&
				   gradUxy_error(:,:),&
				   gradUyy_error(:,:),&
				   gradUzz_error(:,:),&
				   pressure_error(:,:),&
				   Txx_error(:,:),&
				   Txy_error(:,:),&
				   Tyy_error(:,:),&
				   Tzz_error(:,:)

  DOUBLE PRECISION, DIMENSION(2) :: centre_of_sphere
  
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: gl,&
                                                 w,&
                                                 weight_2d,&
                                                 Leg,&
                                                 Leg1,&
                                                 f_x,&
                                                 f_y,&
                                                 V_x,&
                                                 V_y,&
                                                 V_xNm1,&!Stores velocity at previous timestep #1
                                                 V_yNm1,& 
                                                 V_xNm2,&!Stores velocity at previous timestep #2
                                                 V_yNm2,& 
                                                 pressure,&
                                                 g,&
                                                 boundaryContribution_x,&
                                                 boundaryContribution_y,&
                                                 boundaryContribution_p,&
                                                 boundary_x,&
                                                 boundary_y,&
                                                 boundary_stress_xx, boundary_stress_xy, boundary_stress_yy, boundary_stress_zz,&
!                                                  uxx,&
!                                                  uxy,&
!                                                  uyx,&
!                                                  uyy,&
!                                                  uzz,&
!                                                  values,&
!                                                  globalRHS,&
!                                                  globalSOL,&
                                                 OIFS_stokes_contrib_x,&
                                                 OIFS_stokes_contrib_y,&
                                                 stress_cont_to_stokes_x,&
                                                 stress_cont_to_stokes_y,&
                                                 stress_cont_to_stokes_xNm1,&
                                                 stress_cont_to_stokes_yNm1,&
                                                 stress_cont_to_stokes_xNm2,&
                                                 stress_cont_to_stokes_yNm2,&
                                                 OIFSstress_xx,&
                                                 OIFSstress_xy,&
                                                 OIFSstress_yy,&
                                                 OIFSstress_zz,&
                                                 upwinded_element,&
                                                 RHS_x,&
                                                 RHS_y,&
! 						 Txx,&
! 						 Txy,&
! 						 Tyy,&
! 						 Tzz,&
! 						 TxxNm1,&
! 						 TxyNm1,&
! 						 TyyNm1,&
! 						 TzzNm1,&
! 						 TxxNm2,&
! 						 TxyNm2,&
! 						 TyyNm2,&
! 						 TzzNm2,&
						 mesh_velocity,&
						 mesh_velocityNm1,&
						 mesh_velocityNm2,&
						 transient_u,transient_gradUyx,transient_txx,transient_txy
						 !transient_u,transient_txx,transient_txy,transient_cxx,transient_cxy
						 
  
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: vertx,&
                                                   d,&
                                                   evalh,&
                                                   theta,&
                                                   global_matrix,&
                                                   Z_p,&
						   nodeCoord,&
						   nodeCoordNm1,&
						   nodeCoordNm2

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: jac,&
                                                     dxdp,&
                                                     dxde,&
                                                     dydp,&
                                                     dyde,&
						     A_x,&
						     A_y,&
						     B_x,&
						     B_y,&
						     M_pressure,&  ! MAY NEED TO RENAME
						     norm_to_edge,&
						     angles_at_vertex
! DEVSS module memory:						     
  INTEGER, ALLOCATABLE :: devss_ipiv(:)
  DOUBLE PRECISION, ALLOCATABLE :: devss_contrib_x(:), devss_contrib_y(:), &
				   devss_contrib_xNm1(:), devss_contrib_yNm1(:), &
				   devss_contrib_xNm2(:), devss_contrib_yNm2(:), &
				   devss_matrix(:,:,:), devss_rhs_xx(:,:), &
				   devss_rhs_xy(:,:), devss_rhs_yx(:,:), &
				   devss_rhs_yy(:,:), devss_rhs_zz(:,:), &
				   devss_soln_xx(:,:),devss_soln_xy(:,:), &
				   devss_soln_yx(:,:),devss_soln_yy(:,:), &
				   devss_soln_zz(:,:)
  CONTAINS
  SUBROUTINE assignMem
  
    IMPLICIT NONE

! Integer & logical mem
    ALLOCATE( &
      node(1:numelm,0:5),& ! really 4, but want corner node 0 to equal node 4 and 5 to equal corner node 1
      bdflag(1:2,numelm*NP1SQ),& 
      neuflag(1:2,numelm*NP1SQ),&
      circflag(1:numelm),&
      wallflag(numelm*NP1SQ),&
      inflowflag(numelm*NP1SQ),&
      outflowflag(numelm*NP1SQ),&
      wallsymmflag(numelm*NP1SQ),&
      circnodeflag(numelm*NP1SQ),&
      fixed_node(numelm*NP1SQ),&
      accordianflag(1:numelm),&
      conedg(1:numelm,4),&
      conelm(1:numelm,4),&
      mult(numelm*NP1SQ),&
      mapg(0:NP1SQM1,1:numelm),&
      mapg_pressure(1:NM1SQ,1:numelm),&
!       locel(numelm*NP1SQ,1:4), &
!       locnp(numelm*NP1SQ,1:4),&
      local_edge_node(0:N,1:4),&
      local_ij_to_edge_node(0:NP1SQM1,0:5),&
      local_ij_to_i(0:NP1SQM1),&
      local_ij_to_j(0:NP1SQM1),&
      upwind_local_edge_node(0:N,1:4,1:numelm),&
      local_to_interior_node(0:NP1SQM1),&
      interior_to_local_node(1:NM1SQ),&
      is_wallsymm_edge(1:4,1:numelm),&
      is_circ_edge(1:4,1:numelm),&
      is_wall_edge(1:4,1:numelm),&
      is_inflow_edge(1:4,1:numelm),&
      is_outflow_edge(1:4,1:numelm),&
      axisymm_edge(numelm),&
      upwind_wallsymm_1d(1:2,1:4,1:numelm), &
      jac_wallsymm_edge_1d(1:4,numelm)&
    )
  
! double mem
    ALLOCATE( &
      vertx(numnp,2),&
      gl(0:N),&
      w(0:N),&
      weight_2d(0:NP1SQM1),&
      Leg(0:N),&
      Leg1(0:N),&
      d(0:N,0:N),&
      evalh(1:NM1,0:N),&
      jac(0:N,0:N,1:numelm),&
! may comment the below out depending on efficiency of passing them around
      dxdp(0:N,0:N,1:numelm),&
      dxde(0:N,0:N,1:numelm),&
      dydp(0:N,0:N,1:numelm),&
      dyde(0:N,0:N,1:numelm),&
      diff_x(0:NP1SQM1,0:NP1SQM1,numelm),&
      diff_y(0:NP1SQM1,0:NP1SQM1,numelm),&
      A_x(0:NP1SQM1,0:NP1SQM1,1:numelm),&
      A_y(0:NP1SQM1,0:NP1SQM1,1:numelm),&
! Size is, in matrix form, (N-1)^2 by (N+1)^2, the (N-1)^2 comes from pressure part
! and (N+1)^2 from the velocity trial function.
      B_x(1:NM1SQ,0:NP1SQM1,1:numelm),&
      B_y(1:NM1SQ,0:NP1SQM1,1:numelm),&
      Z_p(1:NM1SQ,1:numelm),&
      M_pressure(1:NM1SQ,1:NM1SQ,1:numelm),&
      Mv_x(0:NP1SQM1,0:NP1SQM1,1:numelm),&
      Mv_y(0:NP1SQM1,0:NP1SQM1,1:numelm),&
      theta(2,1:numelm),&
      norm_to_edge(2,4,numelm),&
      angles_at_vertex(2,4,numelm),&      
      storeA_x(0:NP1SQM1,0:NP1SQM1,1:numelm),&
      storeA_y(0:NP1SQM1,0:NP1SQM1,1:numelm),&
      storeB_x(1:NM1SQ,0:NP1SQM1,1:numelm),&
      storeB_y(1:NM1SQ,0:NP1SQM1,1:numelm),&
      storeZ_p(1:NM1SQ,1:numelm),&
!       storeM_stress(0:NP1SQM1,1:numelm),&
      storeM_pressure(1:NM1SQ,1:NM1SQ,1:numelm),&
      storeMv_x(0:NP1SQM1,0:NP1SQM1,1:numelm),&
      storeMv_y(0:NP1SQM1,0:NP1SQM1,1:numelm),&
!       storeC_x(0:NP1SQM1,0:NP1SQM1,1:numelm),&
!       storeC_y(0:NP1SQM1,0:NP1SQM1,1:numelm),&
!       storeCb(0:N,0:N,1:4,1:numelm),&
!       M_stress(0:NP1SQM1,1:numelm),&
!       C_x(0:NP1SQM1,0:NP1SQM1,1:numelm),&
!       C_y(0:NP1SQM1,0:NP1SQM1,1:numelm),&
!       Cb(0:N,0:N,1:4,1:numelm),&
      localTxx(0:NP1SQM1,1:numelm),&
      localTxy(0:NP1SQM1,1:numelm),&
      localTyy(0:NP1SQM1,1:numelm),&
      localTzz(0:NP1SQM1,1:numelm),&
      localTxxNm1(0:NP1SQM1,1:numelm),&
      localTxyNm1(0:NP1SQM1,1:numelm),&
      localTyyNm1(0:NP1SQM1,1:numelm),&
      localTzzNm1(0:NP1SQM1,1:numelm),&
      localTxxNm2(0:NP1SQM1,1:numelm),&
      localTxyNm2(0:NP1SQM1,1:numelm),&
      localTyyNm2(0:NP1SQM1,1:numelm),&
      localTzzNm2(0:NP1SQM1,1:numelm),&
      localConformationxx(0:NP1SQM1,1:numelm),&
      localConformationxy(0:NP1SQM1,1:numelm),&
      localConformationyy(0:NP1SQM1,1:numelm),&
      localConformationzz(0:NP1SQM1,1:numelm),&
      localConformationxxNm1(0:NP1SQM1,1:numelm),&
      localConformationxyNm1(0:NP1SQM1,1:numelm),&
      localConformationyyNm1(0:NP1SQM1,1:numelm),&
      localConformationzzNm1(0:NP1SQM1,1:numelm),&
      localConformationxxNm2(0:NP1SQM1,1:numelm),&
      localConformationxyNm2(0:NP1SQM1,1:numelm),&
      localConformationyyNm2(0:NP1SQM1,1:numelm),&
      localConformationzzNm2(0:NP1SQM1,1:numelm),&
      norm_to_edge_node(2,0:N,1:4,numelm),&
      jac_on_edge(4,numelm) &
!       param_beta_a(numelm)&
   ) 
! could move all non-mapglobal2 data to 2nd routine below...
  END SUBROUTINE assignMem

  SUBROUTINE assignMem2
    IMPLICIT NONE
    ALLOCATE(&
      boundaryContribution_x(nptot),&
      boundaryContribution_y(nptot),&
      boundaryContribution_p(npint),&
      non_dir_bd_map_x(npedg), &
      non_dir_bd_map_y(npedg), &
      nodeCoord(nptot,1:2),&
      nodeCoordNm1(nptot,1:2),&
      nodeCoordNm2(nptot,1:2),& ! needed for AB3
      boundary_x(npedg),&
      boundary_y(npedg),&
      boundary_stress_xx(npedg),&
      boundary_stress_xy(npedg),&
      boundary_stress_yy(npedg),&
      boundary_stress_zz(npedg),&      
      f_x(nptot),&
      f_y(nptot),&
      g(npint),&
      V_x(nptot),&
      V_y(nptot),&
      V_xNm1(nptot),& !Stores velocity at previous timestep #1
      V_yNm1(nptot),&
      V_xNm2(nptot),& !Stores velocity at previous timestep #2
      V_yNm2(nptot),&
      pressure(npint),&
!       uxx(nptot),&
!       uxy(nptot),&
!       uyx(nptot),&
!       uyy(nptot),&
!       uzz(nptot),&
      stress_cont_to_stokes_x(nptot),&
      stress_cont_to_stokes_y(nptot),&
      stress_cont_to_stokes_xNm1(nptot),&
      stress_cont_to_stokes_yNm1(nptot),&
      stress_cont_to_stokes_xNm2(nptot),&
      stress_cont_to_stokes_yNm2(nptot),&
      OIFS_stokes_contrib_x(nptot),&
      OIFS_stokes_contrib_y(nptot),&
      OIFSstress_xx(nptot),&
      OIFSstress_xy(nptot),&
      OIFSstress_yy(nptot),&
      OIFSstress_zz(nptot),&
      upwinded_element(nptot),&
      RHS_x(nptot),& ! no longer required ??
      RHS_y(nptot),&
      storef_x(nptot),&
      storef_y(nptot),&
      storeg_p(npint),&
!       Txx(nptot),&
!       Txy(nptot),&
!       Tyy(nptot),&
!       Tzz(nptot),&
!       TxxNm1(nptot),&
!       TxyNm1(nptot),&
!       TyyNm1(nptot),&
!       TzzNm1(nptot),&
!       TxxNm2(nptot),&
!       TxyNm2(nptot),&
!       TyyNm2(nptot),&
!       TzzNm2(nptot),&
      mesh_velocity(nptot),&
      mesh_velocityNm1(nptot),&
      mesh_velocityNm2(nptot),&
!       transient_u(nptot),&
!       transient_txx(nptot),&
!       transient_txy(nptot),&
!       transient_cxx(nptot),&
!       transient_cxy(nptot),&
      localGradUxx(0:NP1SQM1,1:numelm),&
      localGradUxy(0:NP1SQM1,1:numelm),&
      localGradUyx(0:NP1SQM1,1:numelm),&
      localGradUyy(0:NP1SQM1,1:numelm),&
      localGradUzz(0:NP1SQM1,1:numelm),&
      localGradUxxNm1(0:NP1SQM1,1:numelm),&
      localGradUxyNm1(0:NP1SQM1,1:numelm),&
      localGradUyxNm1(0:NP1SQM1,1:numelm),&
      localGradUyyNm1(0:NP1SQM1,1:numelm),&
      localGradUzzNm1(0:NP1SQM1,1:numelm),&
      localpressure(0:NP1SQM1,1:numelm),&
      localpressureNm1(0:NP1SQM1,1:numelm),&
      localV_x(0:NP1SQM1,1:numelm),&
      localV_y(0:NP1SQM1,1:numelm),&
      localV_xNm1(0:NP1SQM1,1:numelm),&
      localV_yNm1(0:NP1SQM1,1:numelm),&
      global_to_local_map(nptot,1:numelm),&
      fine_uniform_points(0:Nfine),&
      fine_node_coords(1:2,0:Nfine,0:Nfine,1:numelm),&
      fine_hbasis(0:N,0:Nfine),&
      fine_hbasis_tilde(1:N-1,0:Nfine),&
      fine_edge_ifine(0:Nfine,1:4),&
      fine_edge_jfine(0:Nfine,1:4),&
      fine_velocitybasis(0:NP1SQM1,0:(Nfine+1)**2-1),&
      fine_pressurebasis(0:NP1SQM1,0:(Nfine+1)**2-1),&      
      fine_velocity_x(0:Nfine,0:Nfine,1:numelm),&
      fine_velocity_y(0:Nfine,0:Nfine,1:numelm),&
      fine_pressure(0:Nfine,0:Nfine,1:numelm),&
      fine_tao_xx(0:Nfine,0:Nfine,1:numelm),&
      fine_tao_xy(0:Nfine,0:Nfine,1:numelm),&
      fine_tao_yy(0:Nfine,0:Nfine,1:numelm),&
      fine_tao_zz(0:Nfine,0:Nfine,1:numelm),&
      fineGxx(0:Nfine,0:Nfine,1:numelm),&
      fineGxy(0:Nfine,0:Nfine,1:numelm),&
      fineGyx(0:Nfine,0:Nfine,1:numelm),&
      fineGyy(0:Nfine,0:Nfine,1:numelm),&
      fineGzz(0:Nfine,0:Nfine,1:numelm),&
! Analytical solution storage:
      V_x_analytic(0:NP1SQM1,1:numelm),&
      V_y_analytic(0:NP1SQM1,1:numelm),&
      gradUxx_analytic(0:NP1SQM1,1:numelm),&
      gradUyx_analytic(0:NP1SQM1,1:numelm),&
      gradUxy_analytic(0:NP1SQM1,1:numelm),&
      gradUyy_analytic(0:NP1SQM1,1:numelm),&
      pressure_analytic(0:NP1SQM1,1:numelm),&
! error storage:
      V_x_error(0:NP1SQM1,1:numelm),&
      V_y_error(0:NP1SQM1,1:numelm),&
      gradUxx_error(0:NP1SQM1,1:numelm),&
      gradUyx_error(0:NP1SQM1,1:numelm),&
      gradUxy_error(0:NP1SQM1,1:numelm),&
      gradUyy_error(0:NP1SQM1,1:numelm),&
      gradUzz_error(0:NP1SQM1,1:numelm),&
      pressure_error(0:NP1SQM1,1:numelm),&
      Txx_error(0:NP1SQM1,1:numelm),&
      Txy_error(0:NP1SQM1,1:numelm),&
      Tyy_error(0:NP1SQM1,1:numelm),&
      Tzz_error(0:NP1SQM1,1:numelm) &
            )
            
      CALL assign_devss_memory
  END SUBROUTINE
  
  SUBROUTINE deassignMem
    IMPLICIT NONE
    DEALLOCATE(node)
    DEALLOCATE(bdflag)
    DEALLOCATE(neuflag)
    DEALLOCATE(circflag)
    DEALLOCATE(wallflag)
    DEALLOCATE(inflowflag)
    DEALLOCATE(outflowflag)
    DEALLOCATE(wallsymmflag)
    DEALLOCATE(circnodeflag)
    DEALLOCATE(fixed_node)
    DEALLOCATE(accordianflag)
    DEALLOCATE(conedg)
    DEALLOCATE(conelm)
    DEALLOCATE(mult)
    DEALLOCATE(non_dir_bd_map_x)
    DEALLOCATE(non_dir_bd_map_y)
    DEALLOCATE(mapg)
!     DEALLOCATE(locel)
!     DEALLOCATE(locnp)
    DEALLOCATE(local_edge_node)
    DEALLOCATE(upwind_local_edge_node)
    DEALLOCATE(global_to_local_map)
    DEALLOCATE(local_to_interior_node)
    DEALLOCATE(interior_to_local_node)
    DEALLOCATE(is_wallsymm_edge)
    DEALLOCATE(is_circ_edge)
    DEALLOCATE(is_wall_edge)
    DEALLOCATE(is_inflow_edge)
    DEALLOCATE(is_outflow_edge)
    DEALLOCATE(axisymm_edge)
    DEALLOCATE(local_ij_to_edge_node)
    DEALLOCATE(local_ij_to_i)
    DEALLOCATE(local_ij_to_j)
    DEALLOCATE(upwind_wallsymm_1d)
    DEALLOCATE(jac_wallsymm_edge_1d)
    
    DEALLOCATE(vertx)
    DEALLOCATE(gl)
    DEALLOCATE(w)
    DEALLOCATE(weight_2d)
    DEALLOCATE(Leg)
    DEALLOCATE(Leg1)
    DEALLOCATE(d)
    DEALLOCATE(evalh)
    DEALLOCATE(jac)
    DEALLOCATE(dxdp)
    DEALLOCATE(dxde)
    DEALLOCATE(dyde)
    DEALLOCATE(dydp)
    DEALLOCATE(A_x)
    DEALLOCATE(A_y)
    DEALLOCATE(B_x)
    DEALLOCATE(B_y)
    DEALLOCATE(Z_p)
    DEALLOCATE(M_pressure)
    DEALLOCATE(Mv_x,Mv_y)
    DEALLOCATE(theta)
!     DEALLOCATE(Txx,Txy,Tyy,Tzz)
!     DEALLOCATE(transient_u,transient_txx,transient_txy,transient_cxx,transient_cxy)
!     DEALLOCATE(uxx,uxy,uyx,uyy,uzz)
!     DEALLOCATE(TxxNm1,TxyNm1,TyyNm1,TzzNm1)
!     DEALLOCATE(TxxNm2,TxyNm2,TyyNm2,TzzNm2)
    DEALLOCATE(RHS_x)
    DEALLOCATE(RHS_y)
    DEALLOCATE(OIFS_stokes_contrib_x,OIFS_stokes_contrib_y)
    DEALLOCATE(stress_cont_to_stokes_x,stress_cont_to_stokes_y,stress_cont_to_stokes_xNm1,stress_cont_to_stokes_yNm1)
    DEALLOCATE(OIFSstress_xx,OIFSstress_xy,OIFSstress_yy,OIFSstress_zz)
    DEALLOCATE(upwinded_element)
    DEALLOCATE(mesh_velocity,mesh_velocityNm1,mesh_velocityNm2)
    DEALLOCATE(boundaryContribution_x,boundaryContribution_y)
    DEALLOCATE(boundaryContribution_p)
    DEALLOCATE(nodeCoord)
    DEALLOCATE(nodeCoordNm1)
    DEALLOCATE(nodeCoordNm2)
    DEALLOCATE(boundary_x)
    DEALLOCATE(boundary_y) 
    DEALLOCATE(boundary_stress_xx, boundary_stress_xy, boundary_stress_yy, boundary_stress_zz)
    DEALLOCATE(f_x)
    DEALLOCATE(f_y)
    DEALLOCATE(g)
    DEALLOCATE(V_x)
    DEALLOCATE(V_y)
    DEALLOCATE(V_xNm1)
    DEALLOCATE(V_yNm1)
    DEALLOCATE(V_xNm2)
    DEALLOCATE(V_yNm2)
    DEALLOCATE(pressure)
    DEALLOCATE(norm_to_edge_node)
    DEALLOCATE(jac_on_edge)
    DEALLOCATE(diff_x)
    DEALLOCATE(diff_y)
!     DEALLOCATE(param_beta_a)


    DEALLOCATE(norm_to_edge,angles_at_vertex)
    
!     DEALLOCATE(M_stress)
!     DEALLOCATE(C_x,C_y,Cb)
    DEALLOCATE(localTxx, localTxy, localTyy, localTzz,&
    localTxxNm1, localTxyNm1, localTyyNm1, localTzzNm1,&
    localTxxNm2, localTxyNm2, localTyyNm2, localTzzNm2)
    
    DEALLOCATE(localConformationxx, localConformationxy, localConformationyy, localConformationzz,&
    localConformationxxNm1, localConformationxyNm1, localConformationyyNm1, localConformationzzNm1,&
    localConformationxxNm2, localConformationxyNm2, localConformationyyNm2, localConformationzzNm2)

    DEALLOCATE(storeA_x, storeA_y, storeB_x, storeB_y, &
	      storef_x, storef_y, storeg_p, storeZ_p, &
	      storeM_pressure, storeMv_x, storeMv_y)
!     DEALLOCATE(storeM_stress,storeC_x, storeC_y, storeCb)
    
    DEALLOCATE(fine_uniform_points,fine_node_coords,fine_edge_ifine,fine_edge_jfine,&
		fine_hbasis,fine_hbasis_tilde,fine_velocitybasis,fine_pressurebasis,&
		fine_velocity_x,fine_velocity_y,fine_pressure,&
		fine_tao_xx,fine_tao_xy,fine_tao_yy,fine_tao_zz,&
		fineGxx,fineGxy,fineGyx,fineGyy,fineGzz)

    DEALLOCATE(localpressure, localpressureNm1, localV_x, localV_xNm1, localV_y, localV_yNm1, &
		V_x_analytic, V_y_analytic, &
		gradUxx_analytic, gradUyx_analytic, gradUxy_analytic, gradUyy_analytic, pressure_analytic, &
		V_x_error, V_y_error, pressure_error, &
		gradUxx_error, gradUyx_error, gradUxy_error, gradUyy_error, gradUzz_error, &
		Txx_error, Txy_error, Tyy_error, Tzz_error)
		
    CALL deassign_devss_memory
  END SUBROUTINE deassignMem
  
  
    SUBROUTINE assign_devss_memory
    IMPLICIT NONE
    
    IF (.not.ALLOCATED(devss_contrib_x)) ALLOCATE(devss_contrib_x(nptot))
    IF (.not.ALLOCATED(devss_contrib_y)) ALLOCATE(devss_contrib_y(nptot))
    IF (.not.ALLOCATED(devss_contrib_xNm1)) ALLOCATE(devss_contrib_xNm1(nptot))
    IF (.not.ALLOCATED(devss_contrib_yNm1)) ALLOCATE(devss_contrib_yNm1(nptot))
    IF (.not.ALLOCATED(devss_contrib_xNm2)) ALLOCATE(devss_contrib_xNm2(nptot))
    IF (.not.ALLOCATED(devss_contrib_yNm2)) ALLOCATE(devss_contrib_yNm2(nptot))
    IF (.not.ALLOCATED(devss_matrix)) ALLOCATE(devss_matrix(NM1SQ,NM1SQ,numelm))
    IF (.not.ALLOCATED(devss_rhs_xx)) ALLOCATE(devss_rhs_xx(NM1SQ,numelm))
    IF (.not.ALLOCATED(devss_rhs_xy)) ALLOCATE(devss_rhs_xy(NM1SQ,numelm))
    IF (.not.ALLOCATED(devss_rhs_yx)) ALLOCATE(devss_rhs_yx(NM1SQ,numelm))
    IF (.not.ALLOCATED(devss_rhs_yy)) ALLOCATE(devss_rhs_yy(NM1SQ,numelm))
    IF (.not.ALLOCATED(devss_rhs_zz)) ALLOCATE(devss_rhs_zz(NM1SQ,numelm))
    IF (.not.ALLOCATED(devss_soln_xx)) ALLOCATE(devss_soln_xx(0:NP1SQ-1,numelm))
    IF (.not.ALLOCATED(devss_soln_xy)) ALLOCATE(devss_soln_xy(0:NP1SQ-1,numelm))
    IF (.not.ALLOCATED(devss_soln_yx)) ALLOCATE(devss_soln_yx(0:NP1SQ-1,numelm))
    IF (.not.ALLOCATED(devss_soln_yy)) ALLOCATE(devss_soln_yy(0:NP1SQ-1,numelm))
    IF (.not.ALLOCATED(devss_soln_zz)) ALLOCATE(devss_soln_zz(0:NP1SQ-1,numelm))
    IF (.not.ALLOCATED(devss_ipiv)) ALLOCATE(devss_ipiv(NM1SQ))
  
  END SUBROUTINE assign_devss_memory
  
  SUBROUTINE deassign_devss_memory
    IMPLICIT NONE
    
    IF (ALLOCATED(devss_contrib_x)) DEALLOCATE(devss_contrib_x)
    IF (ALLOCATED(devss_contrib_y)) DEALLOCATE(devss_contrib_y)
    IF (ALLOCATED(devss_contrib_xNm1)) DEALLOCATE(devss_contrib_xNm1)
    IF (ALLOCATED(devss_contrib_yNm1)) DEALLOCATE(devss_contrib_yNm1)
    IF (ALLOCATED(devss_contrib_xNm2)) DEALLOCATE(devss_contrib_xNm2)
    IF (ALLOCATED(devss_contrib_yNm2)) DEALLOCATE(devss_contrib_yNm2)    
    IF (ALLOCATED(devss_matrix)) DEALLOCATE(devss_matrix)
    IF (ALLOCATED(devss_rhs_xx)) DEALLOCATE(devss_rhs_xx)
    IF (ALLOCATED(devss_rhs_xy)) DEALLOCATE(devss_rhs_xy)
    IF (ALLOCATED(devss_rhs_yx)) DEALLOCATE(devss_rhs_yx)
    IF (ALLOCATED(devss_rhs_yy)) DEALLOCATE(devss_rhs_yy)
    IF (ALLOCATED(devss_rhs_zz)) DEALLOCATE(devss_rhs_zz)
    IF (ALLOCATED(devss_soln_xx)) DEALLOCATE(devss_soln_xx)
    IF (ALLOCATED(devss_soln_xy)) DEALLOCATE(devss_soln_xy)
    IF (ALLOCATED(devss_soln_yx)) DEALLOCATE(devss_soln_yx)
    IF (ALLOCATED(devss_soln_yy)) DEALLOCATE(devss_soln_yy)
    IF (ALLOCATED(devss_soln_zz)) DEALLOCATE(devss_soln_zz)
    IF (ALLOCATED(devss_ipiv)) DEALLOCATE(devss_ipiv)
  
  END SUBROUTINE deassign_devss_memory

END MODULE shared_data