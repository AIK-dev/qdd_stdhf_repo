!This file is a part of PW-TELEMAN project.
!PW-TELEMAN is a Time-Dependent Electronic Dynamics in Molecules And Nanosystems library.
!Copyright (C) 2011-2015  Paul-Gerhard Reinhard, Eric Suraud, Florent Calvayrac,
!Phuong Mai Dinh, David Brusson, Philipp Wopperer, José María Escartín Esteban.
!
!PW-Teleman is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!PW-Teleman is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with PW-Teleman.  If not, see <http://www.gnu.org/licenses/>.


MODULE cuparams
USE cudafor
USE params
IMPLICIT NONE
SAVE
!
! general settings for dimensions, part.numbers etc
!

REAL(DP) , CONSTANT :: dvol_d
INTEGER,DEVICE,ALLOCATABLE :: ispin_d(:)                 !  spin of states
INTEGER,DEVICE,ALLOCATABLE :: nrel2abs_d(:),nabs2rel_d(:)  !  pointer to wfs
REAL(DP),DEVICE,ALLOCATABLE :: sphermask_d(:)  

#include "cupseudo.F90"

contains

SUBROUTINE cuinit_fields()

! fields for PsP projectors in device memory
IF(ipsptyp /=0) THEN
  ALLOCATE(ifin_d(0:ng),icount_d(knl,0:ng))
  ALLOCATE(p0_1_d(knl,0:ng),p0_2_d(knl,0:ng),p1_1_d(knl,0:ng),p1_1x_d(knl,0:ng))
  ALLOCATE(p1_1y_d(knl,0:ng),p1_1z_d(knl,0:ng) )
  ALLOCATE(p0_3_d(knl,0:ng),p1_2_d(knl,0:ng),p1_2x_d(knl,0:ng) )
  ALLOCATE(p1_2y_d(knl,0:ng),p1_2z_d(knl,0:ng) )
  ALLOCATE(p2_1_d(knl,0:ng),p2_xy_d(knl,0:ng),p2_xz_d(knl,0:ng) )
  ALLOCATE(p2_yz_d(knl,0:ng),p2_xy2_d(knl,0:ng),p2_z2_d(knl,0:ng))
  allocate(h0_11g_d(-99:99),h0_22g_d(-99:99),h0_33g_d(-99:99),h1_11g_d(-99:99),h1_22g_d(-99:99))
  allocate(h2_11g_d(-99:99),h0_12g_d(-99:99))
ENDIF
RETURN

END SUBROUTINE cuinit_fields

END MODULE cuparams
