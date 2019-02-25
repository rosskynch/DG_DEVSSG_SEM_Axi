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

MODULE geometry_module
  USE constants
  USE shared_data
  USE IO_module

  IMPLICIT NONE

  CONTAINS
  
  SUBROUTINE initialise_geometry
    IMPLICIT NONE
    
    CALL calcNodeCoords
    nodeCoordNm1=nodeCoord
    nodeCoordNm2=nodeCoordNm1
    CALL calcJacobian
    
    CALL create_local_ij_to_i_and_j
    CALL create_edge_flags    
!     CALL calc_angle_at_edges
    CALL store_local_edge_nodes
    CALL create_global_to_local_map
    CALL create_local_to_interior_map
    
    CALL calc_normal_to_elements

    CALL calc_domain_area
    mesh_velocity=0d0
    
  END SUBROUTINE initialise_geometry
  
  SUBROUTINE update_geometry
! Updates all the geometry parts required for a new timestep, providing the vertex nodes have already been moved
    IMPLICIT NONE
    nodeCoordNm2=nodeCoordNm1
    nodeCoordNm1=nodeCoord
    CALL calcNodeCoords 
    CALL calcJacobian
    CALL calc_normal_to_elements
    CALL calc_angle_at_edges
! Check for mesh deformation in moving mesh case
    IF (movingmeshflag.eq.1) THEN
      CALL check_accordian_mesh
      CALL calcMeshVelocity
    ENDIF
  END SUBROUTINE update_geometry


  SUBROUTINE create_local_ij_to_i_and_j
    IMPLICIT NONE
    INTEGER :: i,j,ij

    DO i=0,N
      DO j=0,N
        ij=i+j*(N+1)
        local_ij_to_i(ij)=i
        local_ij_to_j(ij)=j
      ENDDO
    ENDDO
  END SUBROUTINE create_local_ij_to_i_and_j


  SUBROUTINE create_edge_flags
    IMPLICIT NONE
    INTEGER :: edge,el
    
    is_wallsymm_edge=.false.
    is_circ_edge=.false.
    is_wall_edge=.false.
    is_inflow_edge=.false.
    is_outflow_edge=.false.
    
    DO el=1,numelm
      DO edge=1,4
        IF (inflowflag(node(el,edge)).and.inflowflag(node(el,edge+1))) THEN
          is_inflow_edge(edge,el)=.true.
        ELSEIF (outflowflag(node(el,edge)).and.outflowflag(node(el,edge+1))) THEN
          is_outflow_edge(edge,el)=.true.
        ELSEIF (wallsymmflag(node(el,edge)).and.wallsymmflag(node(el,edge+1))) THEN
          is_wallsymm_edge(edge,el)=.true.
        ELSEIF (circnodeflag(node(el,edge)).and.circnodeflag(node(el,edge+1))) THEN
          is_circ_edge(edge,el)=.true.
        ELSEIF (wallflag(node(el,edge)).and.wallflag(node(el,edge+1))) THEN
          is_wall_edge(edge,el)=.true.
        ENDIF
      ENDDO
    ENDDO
    
    num_wallsymm_nodes=0
    DO el=1,numelm
      axisymm_edge(el)=0
      DO edge=1,4
        IF (is_wallsymm_edge(edge,el)) THEN
          axisymm_edge(el)=edge
          num_wallsymm_nodes=num_wallsymm_nodes+NP1
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE create_edge_flags
  
  SUBROUTINE create_local_to_interior_map
    IMPLICIT NONE
    INTEGER :: i,j,intij,ij,jj,intjj
    
    local_to_interior_node=-999
    interior_to_local_node=-1
    DO j=1,NM1
      jj=j*NP1
      intjj=(j-1)*NM1
      DO i=1,NM1
        ij=i+jj
        intij=i+intjj
        local_to_interior_node(ij)=intij
        interior_to_local_node(intij)=ij
      ENDDO
    ENDDO
  END SUBROUTINE create_local_to_interior_map


  SUBROUTINE create_global_to_local_map
    IMPLICIT NONE
    INTEGER :: i,ij,el
    
! Initialise to something impossible - will help with debugging
    global_to_local_map=-999
    
    DO el=1,numelm
      DO ij=0,NP1SQM1
        i=mapg(ij,el)
        global_to_local_map(i,el)=ij
      ENDDO
    ENDDO
  END SUBROUTINE create_global_to_local_map
  
  SUBROUTINE store_local_edge_nodes
    IMPLICIT NONE
    INTEGER :: i,j,ij
    
! Edge 1: i varies, j=0
    DO ij=0,N
      local_edge_node(ij,1)=ij
      local_ij_to_edge_node(ij,1)=ij
      local_ij_to_edge_node(ij,5)=ij
    ENDDO    
! Edge 2: i=N, j varies
    DO j=0,N
      ij=N+j*NP1
      local_edge_node(j,2)=ij
      local_ij_to_edge_node(ij,2)=j
    ENDDO
! Edge 3: i varies, j=N
    DO i=0,N
      ij=i+N*NP1
      local_edge_node(i,3)=ij
      local_ij_to_edge_node(ij,3)=i
    ENDDO
! Edge 4: i=0, j varies
    DO j=0,N
      ij = j*NP1
      local_edge_node(j,4)=ij
      local_ij_to_edge_node(ij,4)=j
      local_ij_to_edge_node(ij,0)=j
    ENDDO 
  END SUBROUTINE store_local_edge_nodes
  
  SUBROUTINE calc_domain_area
    IMPLICIT NONE
    INTEGER :: el,p,q
    DOUBLE PRECISION :: temp

    temp=0d0
    DO el=1,numelm
      DO p=0,N
        DO q=0,N
          temp = temp + jac(p,q,el)*w(p)*w(q)
        ENDDO
      ENDDO
    ENDDO
    area_of_domain = temp

  END SUBROUTINE calc_domain_area
  
  SUBROUTINE map_non_dirichlet_edge_nodes
! counts how many non-dirichlet edge nodes we have - doesnt include internal nodes!!!!
! note this is different to the routine in the mkl_solver
! Also creates a map from npedg -> glob_bd_dim for all nodes which will deliberately miss out dirichlet nodes.

  
    IMPLICIT NONE
    INTEGER :: i,counter_x,counter_y
    
    bd_numrows_x=npedg
    bd_numrows_y=npedg
    DO i=1,npedg
      IF (bdflag(1,i)) bd_numrows_x=bd_numrows_x-1
      IF (bdflag(2,i)) bd_numrows_y=bd_numrows_y-1
    ENDDO
    
    global_bd_dim = bd_numrows_x + bd_numrows_y
    global_dim = global_bd_dim + 3*npint
    
    counter_x=0
    counter_y=bd_numrows_x
    DO i=1,npedg
      IF (bdflag(1,i)) THEN
        non_dir_bd_map_x(i) = bd_numrows_x+bd_numrows_y+1
      ELSE
        counter_x=counter_x+1
        non_dir_bd_map_x(i) = counter_x
      ENDIF
      IF (bdflag(2,i)) THEN
        non_dir_bd_map_y(i) = bd_numrows_x+bd_numrows_y+1
      ELSE
        counter_y=counter_y+1
        non_dir_bd_map_y(i) = counter_y
      ENDIF
    ENDDO
    
  END SUBROUTINE map_non_dirichlet_edge_nodes

  
  SUBROUTINE calc_normal_to_elements
! Calculates the normal to each edge of every element.
! Since all element vertices are assigned in a counter-clockwise order, we simply need to take the outward
! normal direction through a clockwise rotation of 90 degrees of the vector from the "beginning" edge node
! to the "ending" edge node. Eg, edge1 is the vector from vertex1 to vertex2, so we rotate this through 270
! degrees, counterclock-wise.
    IMPLICIT NONE
    INTEGER :: el,i,ij,edge,edgep1
    DOUBLE PRECISION :: x1,x2,y1,y2,&
      r11,r12,r21,r22,&
      temp1,temp2,normfact

! Create elements of rotation matrix - counter-clockwise rotation of 270 degrees.
    r11 = 0d0
    r12 = 1d0
    r21 = -1d0
    r22 = 0d0

    DO el=1,numelm
      DO edge=1,4
        edgep1=edge+1
        x1 = vertx(node(el,edge),1)
        x2 = vertx(node(el,edgep1),1)
        y1 = vertx(node(el,edge),2)
        y2 = vertx(node(el,edgep1),2)
        temp1 = x2-x1
        temp2 = y2-y1
        normfact = SQRT(temp1**2+temp2**2)
        norm_to_edge(1,edge,el) = (r11*temp1 + r12*temp2)/normfact
        norm_to_edge(2,edge,el) = (r21*temp1 + r22*temp2)/normfact
      ENDDO
    ENDDO
! Array of normal components for all local edge nodes per element.
    norm_to_edge_node=0d0
    DO el=1,numelm
      DO edge=1,4
        IF (is_circ_edge(edge,el)) THEN
          DO i=0,N
            ij=local_edge_node(i,edge)
            x1=-nodeCoord(mapg(ij,el),1)
            y1=-nodeCoord(mapg(ij,el),2)
            normfact=2d0*rad_sphere!SQRT(x1**2+y1**2)
            norm_to_edge_node(1,i,edge,el)=-x1/normfact
            norm_to_edge_node(2,i,edge,el)=-y1/normfact
          ENDDO
        ELSE
          DO i=0,N
            norm_to_edge_node(1,i,edge,el)=norm_to_edge(1,edge,el)
            norm_to_edge_node(2,i,edge,el)=norm_to_edge(2,edge,el)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE calc_normal_to_elements


  SUBROUTINE calc_angle_at_edges
    IMPLICIT NONE
    INTEGER :: el,i
    DOUBLE PRECISION :: x1,x2,x3,x4,y1,y2,y3,y4,&
      edge1x,edge1y,edge2x,edge2y,&
      ang1,ang2,temp
    
    DO el=1,numelm
      x1 = vertx(node(el,1),1)
      x2 = vertx(node(el,2),1)
      x3 = vertx(node(el,3),1)
      x4 = vertx(node(el,4),1)
      y1 = vertx(node(el,1),2)
      y2 = vertx(node(el,2),2)
      y3 = vertx(node(el,3),2)
      y4 = vertx(node(el,4),2)

! Edge 1:
      edge1x = x4-x1
      edge1y = y4-y1
      edge2x = x2-x1
      edge2y = y2-y1
            
      IF (edge1x.gt.0d0) THEN
        ang1 = atan(edge1y/edge1x)
      ELSEIF (edge1x.lt.0d0.AND.edge1y.gt.0d0) THEN
        ang1 = atan(edge1y/edge1x) + PI
      ELSEIF (edge1x.lt.0d0.AND.edge1y.le.0d0) THEN
        ang1 = atan(edge1y/edge1x) - PI
      ELSEIF (edge1x.eq.0d0.AND.edge1y.gt.0d0) THEN
        ang1 = PI/2d0
      ELSEIF (edge1x.eq.0d0.AND.edge1y.lt.0d0) THEN
        ang1 = -PI/2d0
      ELSE
        ang1 = 0d0
      ENDIF
      
      IF (edge2x.gt.0d0) THEN
        ang2 = atan(edge2y/edge2x)
      ELSEIF (edge2x.lt.0d0.AND.edge2y.gt.0d0) THEN
        ang2 = atan(edge2y/edge2x) + PI
      ELSEIF (edge2x.lt.0d0.AND.edge2y.le.0d0) THEN
        ang2 = atan(edge2y/edge2x) - PI
      ELSEIF (edge2x.eq.0d0.AND.edge2y.gt.0d0) THEN
        ang2 = PI/2d0
      ELSEIF (edge2x.eq.0d0.AND.edge2y.lt.0d0) THEN
        ang2 = -PI/2d0
      ELSE
        ang2 = 0d0
      ENDIF

! Make sure that first angle is smaller than second:
      IF (ang1.gt.ang2) THEN
        temp=ang1
        ang1=ang2
        ang2=temp
      ENDIF

! Check if we have an angle of -PI or PI (i.e. same angle two possible values. Angular scheme always picks -PI.
! If so, match it to the sign of the other angle at that node:
      IF (ang1.eq.-PI) THEN
        IF (ang2.gt.0d0) THEN
! In this case, ang1 would now become the largest, but ang2 should be largest of the 2, so we switch:
          ang1=ang2
          ang2=PI
        ENDIF
      ENDIF

! Store in global memory:
      angles_at_vertex(1,1,el) = ang1
      angles_at_vertex(2,1,el) = ang2

! Edge 2:
      edge1x = x1-x2
      edge1y = y1-y2
      edge2x = x3-x2
      edge2y = y3-y2
            
      IF (edge1x.gt.0d0) THEN
        ang1 = atan(edge1y/edge1x)
      ELSEIF (edge1x.lt.0d0.AND.edge1y.gt.0d0) THEN
        ang1 = atan(edge1y/edge1x) + PI
      ELSEIF (edge1x.lt.0d0.AND.edge1y.le.0d0) THEN
        ang1 = atan(edge1y/edge1x) - PI
      ELSEIF (edge1x.eq.0d0.AND.edge1y.gt.0d0) THEN
        ang1 = PI/2d0
      ELSEIF (edge1x.eq.0d0.AND.edge1y.lt.0d0) THEN
        ang1 = -PI/2d0
      ELSE
        ang1 = 0d0
      ENDIF
      
      IF (edge2x.gt.0d0) THEN
        ang2 = atan(edge2y/edge2x)
      ELSEIF (edge2x.lt.0d0.AND.edge2y.gt.0d0) THEN
        ang2 = atan(edge2y/edge2x) + PI
      ELSEIF (edge2x.lt.0d0.AND.edge2y.le.0d0) THEN
        ang2 = atan(edge2y/edge2x) - PI
      ELSEIF (edge2x.eq.0d0.AND.edge2y.gt.0d0) THEN
        ang2 = PI/2d0
      ELSEIF (edge2x.eq.0d0.AND.edge2y.lt.0d0) THEN
        ang2 = -PI/2d0
      ELSE
        ang2 = 0d0
      ENDIF

! Make sure that first angle is smaller than second:
      IF (ang1.gt.ang2) THEN
        temp=ang1
        ang1=ang2
        ang2=temp
      ENDIF

! Check if we have an angle of -PI or PI (i.e. same angle two possible values. Angular scheme always picks -PI.
! If so, match it to the sign of the other angle at that node:
      IF (ang1.eq.-PI) THEN
        IF (ang2.gt.0d0) THEN
! In this case, ang1 would now become the largest, but ang2 should be largest of the 2, so we switch:
          ang1=ang2
          ang2=PI
        ENDIF
      ENDIF

! Store in global memory:
      angles_at_vertex(1,2,el) = ang1
      angles_at_vertex(2,2,el) = ang2
      
! Edge 3:
      edge1x = x2-x3
      edge1y = y2-y3
      edge2x = x4-x3
      edge2y = y4-y3
            
      IF (edge1x.gt.0d0) THEN
        ang1 = atan(edge1y/edge1x)
      ELSEIF (edge1x.lt.0d0.AND.edge1y.ge.0d0) THEN
        ang1 = atan(edge1y/edge1x) + PI
      ELSEIF (edge1x.lt.0d0.AND.edge1y.lt.0d0) THEN
        ang1 = atan(edge1y/edge1x) - PI
      ELSEIF (edge1x.eq.0d0.AND.edge1y.gt.0d0) THEN
        ang1 = PI/2d0
      ELSEIF (edge1x.eq.0d0.AND.edge1y.lt.0d0) THEN
        ang1 = -PI/2d0
      ELSE
        ang1 = 0d0
      ENDIF

      IF (edge2x.gt.0d0) THEN
        ang2 = atan(edge2y/edge2x)
      ELSEIF (edge2x.lt.0d0.AND.edge2y.gt.0d0) THEN
        ang2 = atan(edge2y/edge2x) + PI
      ELSEIF (edge2x.lt.0d0.AND.edge2y.le.0d0) THEN
        ang2 = atan(edge2y/edge2x) - PI
      ELSEIF (edge2x.eq.0d0.AND.edge2y.gt.0d0) THEN
        ang2 = PI/2d0
      ELSEIF (edge2x.eq.0d0.AND.edge2y.lt.0d0) THEN
        ang2 = -PI/2d0
      ELSE
        ang2 = 0d0
      ENDIF

! Make sure that first angle is smaller than second:
      IF (ang1.gt.ang2) THEN
        temp=ang1
        ang1=ang2
        ang2=temp
      ENDIF

! Check if we have an angle of -PI or PI (i.e. same angle two possible values. Angular scheme always picks -PI.
! If so, match it to the sign of the other angle at that node:
      IF (ang1.eq.-PI) THEN
        IF (ang2.gt.0d0) THEN
! In this case, ang1 would now become the largest, but ang2 should be largest of the 2, so we switch:
        ang1=ang2
        ang2=PI
        ENDIF
      ENDIF

! Store in global memory:
      angles_at_vertex(1,3,el) = ang1
      angles_at_vertex(2,3,el) = ang2

! Edge 4:
      edge1x = x3-x4
      edge1y = y3-y4
      edge2x = x1-x4
      edge2y = y1-y4

      IF (edge1x.gt.0d0) THEN
        ang1 = atan(edge1y/edge1x)
      ELSEIF (edge1x.lt.0d0.AND.edge1y.gt.0d0) THEN
        ang1 = atan(edge1y/edge1x) + PI
      ELSEIF (edge1x.lt.0d0.AND.edge1y.le.0d0) THEN
        ang1 = atan(edge1y/edge1x) - PI
      ELSEIF (edge1x.eq.0d0.AND.edge1y.gt.0d0) THEN
        ang1 = PI/2d0
      ELSEIF (edge1x.eq.0d0.AND.edge1y.lt.0d0) THEN
        ang1 = -PI/2d0
      ELSE
        ang1 = 0d0
      ENDIF
      
      IF (edge2x.gt.0d0) THEN
        ang2 = atan(edge2y/edge2x)
      ELSEIF (edge2x.lt.0d0.AND.edge2y.gt.0d0) THEN
        ang2 = atan(edge2y/edge2x) + PI
      ELSEIF (edge2x.lt.0d0.AND.edge2y.le.0d0) THEN
        ang2 = atan(edge2y/edge2x) - PI
      ELSEIF (edge2x.eq.0d0.AND.edge2y.gt.0d0) THEN
        ang2 = PI/2d0
      ELSEIF (edge2x.eq.0d0.AND.edge2y.lt.0d0) THEN
        ang2 = -PI/2d0
      ELSE
        ang2 = 0d0
      ENDIF

! Make sure that first angle is smaller than second:
      IF (ang1.gt.ang2) THEN
        temp=ang1
        ang1=ang2
        ang2=temp
      ENDIF

! Check if we have an angle of -PI or PI (i.e. same angle two possible values. Angular scheme always picks PI.
! If so, match it to the sign of the other angle at that node:
      IF (ang2.eq.PI.and.ang1.lt.0d0) THEN
! In this case, ang1 would now become the largest, but ang2 should be largest of the 2, so we switch:
        ang2 = ang1
        ang1 = -PI
      ENDIF

! This one is redundant now:
! Check if we have an angle of -PI or PI (i.e. same angle two possible values. Angular scheme always picks -PI.
! If so, match it to the sign of the other angle at that node:
!       IF (ang1.eq.-PI) THEN
!   IF (ang2.gt.0d0) THEN
! In this case, ang1 would now become the largest, but ang2 should be largest of the 2, so we switch:
!     ang1=ang2
!     ang2=PI
!   ENDIF
!       ENDIF

! Store in global memory:
      angles_at_vertex(1,4,el) = ang1
      angles_at_vertex(2,4,el) = ang2
    ENDDO

! Check that all angles lie within 180degrees of each other, otherwise will have the outer, rather than inner, angle.
    DO el=1,numelm
      DO i=1,4
        IF ((angles_at_vertex(2,i,el) - angles_at_vertex(1,i,el)).gt.PI) THEN
          ang1 = angles_at_vertex(1,i,el)
          angles_at_vertex(1,i,el) = angles_at_vertex(2,i,el)
          angles_at_vertex(2,i,el) = ang1
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE calc_angle_at_edges

  
  SUBROUTINE calcMeshVelocity
    IMPLICIT NONE
    INTEGER :: i
    
!     mesh_velocity=V_sphere
!Just use V_sphere ??
!     DO i=1,nptot
!       mesh_velocity(i)=0d0!V_sphere
!     ENDDO


! 1st order mesh velocity:
!     mesh_velocity = (nodeCoord(1:nptot,1) - nodeCoordNm1(1:nptot,1))/deltat
! 2nd order mesh velocity:
!     mesh_velocity = (time_gamma_0*nodeCoord(1:nptot,1) + time_alpha_0*nodeCoordNm1(1:nptot,1) + time_alpha_1*nodeCoordNm2(1:nptot,1))/deltat
    mesh_velocity(1:nptot) = 0.5d0*(3d0*nodeCoord(1:nptot,1) - 4d0*nodeCoordNm1(1:nptot,1) + nodeCoordNm2(1:nptot,1))/deltat
! 3rd order mesh velocity:
!     mesh_velocity = (1d0/12d0)*(23d0*nodeCoord(1:nptot,1) - 16d0*nodeCoordNm1(1:nptot,1) + nodeCoordNm2(1:nptot,1))/deltat
! AB3:
!     mesh_velocity = mesh_velocity + (deltat/12d0)*(23d0*nodeCoord(1:nptot,1) - 16d0*nodeCoordNm1(1:nptot,1) + 5d0*nodeCoordNm2(1:nptot,1))

  END SUBROUTINE calcMeshVelocity
  
  SUBROUTINE check_accordian_mesh
! Must check, in the case of an accordian mesh, that the elements which expand/contract do not get too distorted.
! The general rule is that the width never exceeds twice the height, and never drops below half the height.
    IMPLICIT NONE
    INTEGER :: el
    DOUBLE PRECISION :: dist1,dist2,dist3,dist4,ratio12,ratio23,ratio34,ratio41
    DOUBLE PRECISION :: thresh1,thresh2
    
    thresh1=120d0
    thresh2=0.0075d0
    DO el=1,numelm
      IF (accordianflag(el)) THEN
! Calc distances of each side of the element:
        dist1 = SQRT((5d-1*dxdp(0,0,el))**2 + (5d-1*dydp(0,0,el))**2)
        dist2 = SQRT((5d-1*dxde(N,0,el))**2 + (5d-1*dyde(N,0,el))**2)
        dist3 = SQRT((5d-1*dxdp(N,N,el))**2 + (5d-1*dydp(N,N,el))**2)
        dist4 = SQRT((5d-1*dxde(0,0,el))**2 + (5d-1*dyde(0,0,el))**2)
        ratio12 = dist1/dist2
        ratio23 = dist2/dist3
        ratio34 = dist3/dist4
        ratio41 = dist4/dist1
! Check the ratios of each side
! We allow ratios of 20:1 or less:
! Could do this in 1 IF statement, but we may wish to allow different ratios for
! horizontal and vertical sides in the future?
        IF (ratio12.lt.thresh2.OR.ratio12.gt.thresh1) THEN
          write(*,*) 'Accordian mesh element ',el,' becoming too deformed! Ratio12 is, ',ratio12
          STOP
        ELSEIF (ratio23.lt.thresh2.OR.ratio23.gt.thresh1) THEN
          write(*,*) 'Accordian mesh element ',el,' becoming too deformed! Ratio23 is, ',ratio23
          STOP
        ELSEIF (ratio34.lt.thresh2.OR.ratio34.gt.thresh1) THEN
          write(*,*) 'Accordian mesh element ',el,' becoming too deformed! Ratio34 is, ',ratio34
          STOP
        ELSEIF (ratio41.lt.thresh2.OR.ratio41.gt.thresh1) THEN
          write(*,*) 'Accordian mesh element ',el,' becoming too deformed! Ratio41 is, ',ratio41
          STOP
        ENDIF
      ENDIF
    ENDDO
  END SUBROUTINE check_accordian_mesh

! global/local mapping routine
! needs streamlining in terms of the variables
  SUBROUTINE create_global_map

! mapg(np,el) maps local node np=i+j*(N+1) to global node
! locel(np,mult)  global node np to local element, lowest
!                 element for lowest multiplicity
! locnp(np,mult)  global node np to local np=i+j*(N+1)
!
! Input is numnp, mapg, conelm and conedg
! mapg will already contain the first numnp nodes from the
! input file (ie the vertex nodes)
! conelm(i,j) is the element which element i, face j is touching
! conedg(i,j) is the face which element i, face j is touching on
!             the element identified in conelm(i,j)
 
    IMPLICIT NONE
    INTEGER :: i,j,k,l,el,ed,beg,s1,ij,intij,jj,intjj,&!,co,mul,el2,mm,
      vis(numelm*NP1SQ)
!   INTEGER, DIMENSION(numelm*(N+1)**2) :: vis
    LOGICAL, DIMENSION(numelm,1:4) :: visit
!   LOGICAL, DIMENSION(numelm*(N+1)**2) :: visnp
    INTEGER, DIMENSION(numelm,1:4) :: npbegin

!   INTEGER, DIMENSION(numelm*NP1SQM1,1:4) :: locel,locnp
!   INTEGER, DIMENSION(0:NP1SQM1) :: mapm2ln,mapln2m
  
! numnp = vertex nodes.
! numedg = total boundary nodes ie, number of edge nodes (without vertices is numedg-numnp)
! numint = internal nodes
! numtot = total nodes
  
! initialise visit array
    visit=.false.

! Assign map to boundary nodes, this is done using the face connectivity information
! At the same time, non-vertex boundary node flags (wall/inflow/outflow/neumann) are generated.
    npedg=numnp+1
    DO i=1,numelm
      IF (.not.visit(i,1)) THEN
        npbegin(i,1)=npedg
! Fill in boundary flags for non-vertex boundary nodes in each element
        DO k=1,N-1
          l=0
          mapg(k+l*(N+1),i)=npedg
          mult(npedg)=1
!  
! Cycle over x and y (j=1 and 2) for dirichlet & neumann node flags
! TO EXPAND TO 3-D, j will go to 3
! 
          DO j=1,2
            bdflag(j,npedg)=.true.
            IF (neuflag(j,node(i,1)).AND.neuflag(j,node(i,2))) THEN
              bdflag(j,npedg)=.false.
              neuflag(j,npedg)=.true.
            ENDIF
          ENDDO
  
! Now fill node flags for boundary types  
          IF (inflowflag(node(i,1)).AND.inflowflag(node(i,2))) THEN
            inflowflag(npedg)=.true.
          ELSEIF (outflowflag(node(i,1)).AND.outflowflag(node(i,2))) THEN
            outflowflag(npedg)=.true.
          ELSEIF (wallflag(node(i,1)).AND.wallflag(node(i,2))) THEN
            wallflag(npedg)=.true.
          ELSEIF (circnodeflag(node(i,1)).AND.circnodeflag(node(i,2))) THEN
            circnodeflag(npedg)=.true.
          ELSEIF (wallsymmflag(node(i,1)).AND.wallsymmflag(node(i,2))) THEN
            wallsymmflag(npedg)=.true.
!  ELSE internal boundary node
          ENDIF
          npedg=npedg+1
        ENDDO
      ELSE ! Not a boundary node - a shared internal boundary !
        beg=npbegin(conelm(i,1),conedg(i,1))+(N-2)
        DO k=1,N-1
          l=0
          mapg(k+l*(N+1),i)=beg
          mult(beg)=2
          DO j=1,2
            bdflag(j,beg)=.false.
            neuflag(j,beg)=.false.
          ENDDO
          wallflag(beg)=.false.
          inflowflag(beg)=.false.
          outflowflag(beg)=.false.
          wallsymmflag(beg)=.false.
          circnodeflag(beg)=.false.
  
          beg=beg-1
        ENDDO
      ENDIF
      IF (.not.visit(i,2)) THEN
        npbegin(i,2)=npedg
        DO l=1,N-1
          k=N
          mapg(k+l*(N+1),i)=npedg
          mult(npedg)=1
!  
! Cycle over x and y (j=1 and 2) for dirichlet & neumann node flags
! TO EXPAND TO 3-D, j will go to 3
! 
          DO j=1,2
            bdflag(j,npedg)=.true.
            IF (neuflag(j,node(i,2)).AND.neuflag(j,node(i,3))) THEN
              bdflag(j,npedg)=.false.
              neuflag(j,npedg)=.true.
            ENDIF
          ENDDO
  
! Now fill node flags for boundary types  
          IF (inflowflag(node(i,2)).AND.inflowflag(node(i,3))) THEN
            inflowflag(npedg)=.true.
          ELSEIF (outflowflag(node(i,2)).AND.outflowflag(node(i,3))) THEN
            outflowflag(npedg)=.true.
          ELSEIF (wallflag(node(i,2)).AND.wallflag(node(i,3))) THEN
            wallflag(npedg)=.true.
          ELSEIF (circnodeflag(node(i,2)).AND.circnodeflag(node(i,3))) THEN
            circnodeflag(npedg)=.true.  
          ELSEIF (wallsymmflag(node(i,2)).AND.wallsymmflag(node(i,3))) THEN
            wallsymmflag(npedg)=.true.
!          ELSE internal boundary node
          ENDIF

          npedg=npedg+1
        ENDDO
      ELSE
        beg=npbegin(conelm(i,2),conedg(i,2))+(N-2)
        DO l=1,N-1
          k=N
          mapg(k+l*(N+1),i)=beg
          mult(beg)=2
          DO j=1,2
            bdflag(j,beg)=.false.
            neuflag(j,beg)=.false.
          ENDDO
          wallflag(beg)=.false.
          inflowflag(beg)=.false.
          outflowflag(beg)=.false.
          wallsymmflag(beg)=.false.
          circnodeflag(beg)=.false.
  
          beg=beg-1
        ENDDO
      ENDIF
      IF (.not.visit(i,3)) THEN
        npbegin(i,3)=npedg
        DO k=N-1,1,-1
          l=N
          mapg(k+l*(N+1),i)=npedg
          mult(npedg)=1
!  
! Cycle over x and y (j=1 and 2) for dirichlet & neumann node flags
! TO EXPAND TO 3-D, j will go to 3
! 
          DO j=1,2
            bdflag(j,npedg)=.true.
            IF (neuflag(j,node(i,3)).AND.neuflag(j,node(i,4))) THEN
              bdflag(j,npedg)=.false.
              neuflag(j,npedg)=.true.
            ENDIF
          ENDDO

! Now fill node flags for boundary types  
          IF (inflowflag(node(i,3)).AND.inflowflag(node(i,4))) THEN
            inflowflag(npedg)=.true.
          ELSEIF (outflowflag(node(i,3)).AND.outflowflag(node(i,4))) THEN
            outflowflag(npedg)=.true.
          ELSEIF (wallflag(node(i,3)).AND.wallflag(node(i,4))) THEN
            wallflag(npedg)=.true.
          ELSEIF (circnodeflag(node(i,3)).AND.circnodeflag(node(i,4))) THEN
            circnodeflag(npedg)=.true.
          ELSEIF (wallsymmflag(node(i,3)).AND.wallsymmflag(node(i,4))) THEN
            wallsymmflag(npedg)=.true.
!  ELSE internal boundary node
          ENDIF

          npedg=npedg+1
        ENDDO
      ELSE
        beg=npbegin(conelm(i,3),conedg(i,3))+(N-2)
        DO k=N-1,1,-1
          l=N
          mapg(k+l*(N+1),i)=beg
          mult(beg)=2
          DO j=1,2
            bdflag(j,beg)=.false.
            neuflag(j,beg)=.false.
          ENDDO
          wallflag(beg)=.false.
          inflowflag(beg)=.false.
          outflowflag(beg)=.false.
          wallsymmflag(beg)=.false.
          circnodeflag(beg)=.false.
  
          beg=beg-1
        ENDDO
      ENDIF 
      IF (.not.visit(i,4)) THEN
        npbegin(i,4)=npedg
        DO l=N-1,1,-1
          k=0
          mapg(k+l*(N+1),i)=npedg
          mult(npedg)=1
  
!  
! Cycle over x and y (j=1 and 2) for dirichlet & neumann node flags
! TO EXPAND TO 3-D, j will go to 3
! 
          DO j=1,2
            bdflag(j,npedg)=.true.
            IF (neuflag(j,node(i,4)).AND.neuflag(j,node(i,1))) THEN
              bdflag(j,npedg)=.false.
              neuflag(j,npedg)=.true.
            ENDIF
          ENDDO

! Now fill node flags for boundary types
          IF (inflowflag(node(i,4)).AND.inflowflag(node(i,1))) THEN
            inflowflag(npedg)=.true.
          ELSEIF (outflowflag(node(i,4)).AND.outflowflag(node(i,1))) THEN
            outflowflag(npedg)=.true.
          ELSEIF (wallflag(node(i,4)).AND.wallflag(node(i,1))) THEN
            wallflag(npedg)=.true.
          ELSEIF (circnodeflag(node(i,4)).AND.circnodeflag(node(i,1))) THEN
            circnodeflag(npedg)=.true.
          ELSEIF (wallsymmflag(node(i,4)).AND.wallsymmflag(node(i,1))) THEN
            wallsymmflag(npedg)=.true.
!  ELSE internal boundary node
          ENDIF

          npedg=npedg+1
        ENDDO
      ELSE
        beg=npbegin(conelm(i,4),conedg(i,4))+(N-2)
        DO l=N-1,1,-1
          k=0
          mapg(k+l*(N+1),i)=beg
          mult(beg)=2
          DO j=1,2
            bdflag(j,beg)=.false.
            neuflag(j,beg)=.false.
          ENDDO
          wallflag(beg)=.false.
          inflowflag(beg)=.false.
          outflowflag(beg)=.false.
          wallsymmflag(beg)=.false.
          circnodeflag(beg)=.false.
          beg=beg-1
        ENDDO
      ENDIF

      IF (i.le.numelm-1) THEN
        DO el=i+1,numelm
          DO ed=1,4
            IF (conelm(el,ed).eq.i) THEN
              visit(el,ed)=.true.
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    npedg=npedg-1
    nptot=npedg+1

! Finally map the internal nodes, labelled as normal in "kl" format
    DO i=1,numelm
      DO l=1,N-1
        DO k=1,N-1
          mapg(k+l*(N+1),i)=nptot
          mult(nptot)=1
          bdflag(1,nptot)=.false.
          bdflag(2,nptot)=.false.
          neuflag(1,nptot)=.false.
          neuflag(2,nptot)=.false.
          wallflag(nptot)=.false.
          inflowflag(nptot)=.false.
          outflowflag(nptot)=.false.
          wallsymmflag(nptot)=.false.
          circnodeflag(nptot)=.false.
          nptot=nptot+1
        ENDDO
      ENDDO
    ENDDO
    nptot=nptot-1
    npint=nptot-npedg
!   write(*,*) 'nr. of vertices           =',numnp
!   write(*,*) 'nr. of vertices+edgenodes =',npedg
!   write(*,*) 'nr. of edgenodes          =',npedg-numnp
!   write(*,*) 'nr. of total global nodes =',nptot
!   write(*,*) 'nr. of internal nodes     =',npint


! do i=numnp+1,nptot-npint
!   print*,i,wallflag(i),inflowflag(i),outflowflag(i),wallsymmflag(i),circnodeflag(i)
!   enddo

    vis=0
    DO j=1,numelm
      DO i=0,NP1SQM1
        vis(mapg(i,j))=vis(mapg(i,j))+1
!        locel(mapg(i,j),vis(mapg(i,j)))=j
!        locnp(mapg(i,j),vis(mapg(i,j)))=i
      ENDDO
    ENDDO
  
! Not convinced these are setup right
! If there are any errors, the i,j parts may not actually help??
  
    DO j=1,numelm
      DO i=0,NP1SQM1
        IF (vis(mapg(i,j)).ne.mult(mapg(i,j))) THEN
          write(*,*) 'ERROR in global to local mapping'
          write(*,*) j,i,mapg(i,j),vis(mapg(i,j)),mult(mapg(i,j))
          STOP
        ENDIF
      ENDDO
    ENDDO

    s1=0
    DO i=1,nptot
!       DO j=1,mult(i)
!        write(*,*) 'gn to m',i,'-->',locel(i,j),locnp(i,j)
!       ENDDO
      s1=s1+mult(i)
    ENDDO
    IF (s1.ne.numelm*(N+1)**2) THEN
      write(*,*) 'ERROR in multiplicity global to local'
      write(*,*) s1,i,j,mapg(j,i),vis(mapg(j,i)),mult(mapg(j,i))
      STOP
    ENDIF


! Finally, check if any nodes have been simultaneously labelled both neumann and dirichlet.
! For example, if a node is on the corner of a boundary which changes from an inflow edge, to a symmetric wall edge.
! Here we wish to keep the node as inflow (ie, dirichlet), and remove its neumann status.
! Similarly, if 
    DO i=1,numnp
      DO j=1,2
        IF ((outflowflag(i).OR.inflowflag(i).OR.circnodeflag(i).OR.wallflag(i)).AND.neuflag(j,i)) THEN
    neuflag(j,i)=.false.
    bdflag(j,i)=.true.
        ENDIF
        IF (.NOT.(outflowflag(i).OR.inflowflag(i).OR.circnodeflag(i).OR.wallflag(i)).AND.(wallsymmflag(i).AND.neuflag(j,i))) THEN
    neuflag(j,i)=.true.
    bdflag(j,i)=.false.
        ENDIF      
      ENDDO
    ENDDO
  
  
! Create pressure space mapg. This saves a lot of messing about later:
    DO el=1,numelm
      DO j=1,N-1
        jj=j*NP1
        intjj=(j-1)*NM1
        DO i=1,N-1
    ij=i+jj
    intij=i+intjj
    mapg_pressure(intij,el) = mapg(ij,el) - npedg
        ENDDO
      ENDDO
    ENDDO

! Modified Basis adjustment - the neuman nodes on a symmetric edge become dirichlet
!     IF (modified_basis_flag) THEN
!       DO i=1,npedg
!   IF (wallsymmflag(i)) THEN
!     bdflag(1,i) = .true.
!     neuflag(1,i) = .false.
! 
! ! Not req'd, node is dirichlet in either modified or normal basis.
! !     bdflag(2,i) = .true.
! !     neuflag(2,i) = .false.
!   ENDIF
!       ENDDO
!     ENDIF
  END SUBROUTINE create_global_map
  
 

  DOUBLE PRECISION FUNCTION trnsfnmap(x,y,el,k) RESULT(answer)
!
! local co-ords (ie, gauss-lobatto points)
! x - the psi co-ordinate 
! y - the eta co-ordinate
! el - the local element we're mapping to
! k - takes value 1 or 2, 1 for x co-ordinate, 2 for y co-ordinate
!
    IMPLICIT NONE
    INTEGER :: el,k
    DOUBLE PRECISION :: x,y

    IF (k.gt.2.or.k.lt.1) THEN
      write(*,*) 'ERROR: called incorrect co-ordinate in transfinite mapping'
      RETURN
    ENDIF
   
! CONDITION TO DEAL WITH THE CIRCULAR BOUNDARY CASE
    IF (circflag(el)) THEN
      IF (k.eq.1) THEN
        answer = 5d-1*(1d0 - y)*(centre_of_sphere(1) + rad_sphere*cos(5d-1*((1d0 - x)*theta(1,el) + &
          (1d0 + x)*theta(2,el)))) + 25d-2*((1d0 + y)*((1d0 + x)*vertx(node(el,3),k) + (1d0 - x)*vertx(node(el,4),k)))
      ELSE 
        answer = 5d-1*(1d0 - y)*(centre_of_sphere(2) + rad_sphere*sin(5d-1*((1d0 - x)*theta(1,el) + (1d0 + x)*theta(2,el)))) + &
          25d-2*(1d0 + y)*((1d0 + x)*vertx(node(el,3),k) + (1d0 - x)*vertx(node(el,4),k))
      ENDIF
    ELSE
      answer = 25d-2*( &
        (1d0 - x)*(1d0 - y)*vertx(node(el,1),k) + &
        (1d0 + x)*(1d0 - y)*vertx(node(el,2),k) + &
        (1d0 + x)*(1d0 + y)*vertx(node(el,3),k) + &
        (1d0 - x)*(1d0 + y)*vertx(node(el,4),k) )
    ENDIF
!     IF (abs(answer).lt.1d-14) answer=0d0
  END FUNCTION trnsfnmap
  
  SUBROUTINE calcNodeCoords
! Calculates the co-ordinates in physical space using the transfinite map
    IMPLICIT NONE
    INTEGER :: el,i,j,ij,k,jj
    LOGICAL,DIMENSION(nptot) :: visit
    visit=.true.
    DO el=1,numelm
      DO j=0,N
        jj=j*NP1
        DO i=0,N
          ij=i+jj
          IF (visit(mapg(ij,el))) THEN
            DO k=1,2
              nodeCoord(mapg(ij,el),k) = trnsfnmap(gl(i),gl(j),el,k)
            ENDDO
            visit(mapg(ij,el))=.false.
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE calcNodeCoords
  
  SUBROUTINE calcJacobian

! calculates value of components of the jacobian at G-L points:
! gl - vector of Gauss-Lobatto points
! vertx - vector containing x1,..,x4 - the corners of the quadrilateral element - use node in conjunction with this.
! output - dxdp, dydp, dxde, dyde, J (det(J) - the jacobian)
    IMPLICIT NONE
   
    INTEGER :: i,j,el,k
    DOUBLE PRECISION :: temp
!    DOUBLE PRECISION, DIMENSION(0:N,0:N) :: dxdp,dydp,dxde,dyde,jac


    DO el=1,numelm
! Circular Boundary case:
      IF (circflag(el)) THEN
        DO j=0,N
          DO i=0,N
! NOTE: REMOVED -1 from the transfinite mapping in y direction - circle must now centre @ (0,0) to work!
            dxdp(i,j,el) = 0.25d0*(rad_sphere*(1d0 - gl(j))*(theta(1,el)-theta(2,el))* &
              sin(0.5*((1d0 - gl(i))*theta(1,el) + (1d0+gl(i))*theta(2,el))) + &
              (1d0 + gl(j))*(vertx(node(el,3),1)-vertx(node(el,4),1)))
          
            dxde(i,j,el) = -0.5d0*(centre_of_sphere(1) + rad_sphere*cos(0.5d0*((1d0-gl(i))*theta(1,el) +&
              (1d0 + gl(i))*theta(2,el)))) + &
              0.25d0*((1d0+gl(i))*vertx(node(el,3),1) + (1d0 - gl(i))*vertx(node(el,4),1))
        
            dydp(i,j,el) = 0.25d0*(rad_sphere*(1d0 - gl(j))*(theta(2,el)-theta(1,el))* &
              cos(0.5d0*((1d0 - gl(i))*theta(1,el) + (1d0 + gl(i))*theta(2,el))) + &
              (1d0 + gl(j))*(vertx(node(el,3),2)-vertx(node(el,4),2)))
      
            dyde(i,j,el) = -0.5d0*(centre_of_sphere(2) + rad_sphere*sin(0.5d0*((1d0 - gl(i))*theta(1,el) +&
              (1d0 + gl(i))*theta(2,el)))) + & 
              0.25d0*((1d0 + gl(i))*vertx(node(el,3),2) + (1d0 - gl(i))*vertx(node(el,4),2))
          ENDDO
        ENDDO
      ELSE
! Non-Circular boundary case:
        DO j=0,N
          DO i=0,N
            dxdp(i,j,el) = 0.25d0*( &
              (1d0-gl(j))*(vertx(node(el,2),1)-vertx(node(el,1),1)) + &
              (1d0+gl(j))*(vertx(node(el,3),1)-vertx(node(el,4),1)) )

            dydp(i,j,el) = 0.25d0*( &
              (1d0-gl(j))*(vertx(node(el,2),2)-vertx(node(el,1),2)) + & 
              (1d0+gl(j))*(vertx(node(el,3),2)-vertx(node(el,4),2)) )

            dxde(i,j,el) = 0.25d0*( &
              (1d0-gl(i))*(vertx(node(el,4),1)-vertx(node(el,1),1)) + & 
              (1d0+gl(i))*(vertx(node(el,3),1)-vertx(node(el,2),1)) )

            dyde(i,j,el) = 0.25d0*( &
              (1d0-gl(i))*(vertx(node(el,4),2)-vertx(node(el,1),2)) + &
              (1d0+gl(i))*(vertx(node(el,3),2)-vertx(node(el,2),2)) )
          ENDDO
        ENDDO
      ENDIF
! calculate the jacobian for the current element, el
      DO i=0,N
        DO j=0,N
          temp = dxdp(i,j,el)*dyde(i,j,el)-dydp(i,j,el)*dxde(i,j,el)
          jac(i,j,el) = temp
          IF (temp.lt.0d0) THEN
            print*,'ERROR: Negative Jacobian!!!',i,j,el,temp, 'PRINTED GRID TO FILE'
!             print*,'more info ',dxdp(i,j,el),dyde(i,j,el),dydp(i,j,el),dxde(i,j,el)
!             print*,circflag(el),gl(i),gl(j)
!             write(*,*) (vertx(node(el,k),1),k=1,4)
!             write(*,*) (vertx(node(el,k),2),k=1,4)
!             write(*,*) (nodeCoord(node(el,k),1),k=1,4)
!             write(*,*) (nodeCoord(node(el,k),2),k=1,4)
!             OPEN(tecplot_output_fileid,FILE=tecplot_output_filename)
            OPEN(tecplot_output_fileid,FILE=tecplot_output_filename)
            CALL output_to_tecplot
            CLOSE(tecplot_output_fileid)
            STOP
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    
! calculate the jacobian for a surface integral on the edge of each element:
    DO el=1,numelm

! edge 1 - note, at the moment we only allow the first edge to be circular.
      IF (circflag(el)) THEN
        jac_on_edge(1,el) = 0.5d0*rad_sphere*(theta(2,el)-theta(1,el))
      ELSE
        jac_on_edge(1,el) = 0.5d0*SQRT((vertx(node(el,2),1)-vertx(node(el,1),1))**2 + (vertx(node(el,2),2)-vertx(node(el,1),2))**2)
      ENDIF
! edge 2
      jac_on_edge(2,el) = 0.5d0*SQRT((vertx(node(el,3),1)-vertx(node(el,2),1))**2 + (vertx(node(el,3),2)-vertx(node(el,2),2))**2)
! edge 3
      jac_on_edge(3,el) = 0.5d0*SQRT((vertx(node(el,3),1)-vertx(node(el,4),1))**2 + (vertx(node(el,3),2)-vertx(node(el,4),2))**2)
! edge 4
       jac_on_edge(4,el) = 0.5d0*SQRT((vertx(node(el,4),1)-vertx(node(el,1),1))**2 + (vertx(node(el,4),2)-vertx(node(el,1),2))**2)

      jac_wallsymm_edge_1d(1,el) = 0.5d0*(vertx(node(el,2),1)-vertx(node(el,1),1))
      jac_wallsymm_edge_1d(2,el) = 0.5d0*(vertx(node(el,3),1)-vertx(node(el,2),1))
      jac_wallsymm_edge_1d(3,el) = 0.5d0*(vertx(node(el,3),1)-vertx(node(el,4),1))
      jac_wallsymm_edge_1d(4,el) = 0.5d0*(vertx(node(el,4),1)-vertx(node(el,1),1))
    ENDDO
  END SUBROUTINE calcJacobian
  
  SUBROUTINE create_local_uniform_finegrid!(Nf,grid_out)
! Creates a uniform grid in the parent element (ie, from -1 to 1)
! It has Nf+1 points, with 2 being on the boundary points.
    IMPLICIT NONE
    INTEGER :: i
    DOUBLE PRECISION :: inc


    inc=2d0/dfloat(Nfine)
    fine_uniform_points(0)=-1d0
    fine_uniform_points(Nfine)=1d0

    DO i=1,Nfine-1
      fine_uniform_points(i) = dfloat(i)*inc - 1d0
    ENDDO
  END SUBROUTINE create_local_uniform_finegrid
  
  SUBROUTINE create_fine_basis_array!(Nf,grid_in,basis_out,tilde_basis_out)
! Pre-evaluates the value of the basis functions, h_i(x) on the
! refined grid (i.e. Non-GLL nodes) which is given as input
    IMPLICIT NONE
    INTEGER :: i,ifine,j,jj,ij,jfine,jjf,ijfine
        
    DO ifine=0,Nfine
      DO i=0,N
        fine_hbasis(i,ifine) = hbasis(i,fine_uniform_points(ifine))
      ENDDO
      DO i=1,NM1
        fine_hbasis_tilde(i,ifine) = hbasis_tilde(i,fine_uniform_points(ifine))
      ENDDO
    ENDDO
    
    DO jfine=0,Nfine
      jjf=jfine*(Nfine+1)
      DO ifine=0,Nfine
        ijfine=ifine + jjf
        DO j=0,N
          jj=j*NP1
          DO i=0,N
            ij=i+jj
            fine_velocitybasis(ij,ijfine) = fine_hbasis(i,ifine)*fine_hbasis(j,jfine)
          ENDDO
        ENDDO
        DO j=1,NM1
          jj= (j-1)*NM1
          DO i=1,NM1
            ij=i+jj
            fine_pressurebasis(ij,ijfine) = fine_hbasis_tilde(i,ifine)*fine_hbasis_tilde(j,jfine)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE create_fine_basis_array
  
  SUBROUTINE create_fine_node_points
! Creates an array of the (x,y) co-ordinates for the fine uniform grid. Stored element-wise
    IMPLICIT NONE
    INTEGER :: el,ifine,jfine
    
    DO el=1,numelm
      DO jfine=0,Nfine
        DO ifine=0,Nfine
          fine_node_coords(1,ifine,jfine,el) = trnsfnmap(fine_uniform_points(ifine),fine_uniform_points(jfine),el,1)
          fine_node_coords(2,ifine,jfine,el) = trnsfnmap(fine_uniform_points(ifine),fine_uniform_points(jfine),el,2)
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE create_fine_node_points
  
  SUBROUTINE create_fine_edge_points
    IMPLICIT NONE
    INTEGER :: ifine,jfine,edge
    
    DO edge=1,4
      IF (edge.eq.1) THEN
        fine_edge_ifine(0:Nfine,edge) = 0
        DO jfine=0,Nfine
          fine_edge_jfine(jfine,edge) = jfine
        ENDDO
      ELSEIF (edge.eq.2) THEN
        DO ifine=0,Nfine 
          fine_edge_ifine(ifine,edge) = ifine
        ENDDO
        fine_edge_jfine(0:Nfine,edge) = Nfine
      ELSEIF (edge.eq.3) THEN
        fine_edge_ifine(0:Nfine,edge) = Nfine
        DO jfine=0,Nfine
          fine_edge_jfine(jfine,edge) = jfine
        ENDDO
      ELSEIF (edge.eq.4) THEN
        DO ifine=0,Nfine
          fine_edge_ifine(ifine,edge) = ifine
        ENDDO
        fine_edge_jfine(0:Nfine,edge) = 0
      ENDIF
    ENDDO
  END SUBROUTINE create_fine_edge_points

  SUBROUTINE initialise_fine_grid
    IMPLICIT NONE

    CALL create_local_uniform_finegrid
    CALL create_fine_basis_array
    CALL create_fine_node_points

  END SUBROUTINE initialise_fine_grid
  
  SUBROUTINE update_fine_grid
    IMPLICIT NONE

    IF (movingmeshflag.eq.1) THEN
      CALL create_fine_node_points
    ENDIF

  END SUBROUTINE update_fine_grid


END MODULE geometry_module
