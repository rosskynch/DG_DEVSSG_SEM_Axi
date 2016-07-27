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

! MODULE NO LONGER USED !
MODULE dim_data
!
! This module contains the small amount of information on the Fdimensions required
! for storage array setup.
! for example:
! N - the number of GL points
! numelm - the number of elements
!
! Further dimensional data can be added as required
!
  IMPLICIT NONE
 INTEGER :: N=3,numelm,preconflag=0
 
END MODULE dim_data
