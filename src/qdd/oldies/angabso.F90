! was part of abso_bc.F90
! works only for rough mesh of angular points

!------------------------------------------------------------

SUBROUTINE angabso

! Collects information on particles lost IN angular cones,
! for later analysis of photo-electron angular distribution (PAD).

  USE params
  IMPLICIT NONE

  INTEGER :: ind, iphi, itheta, ix, iy, iz, jj
  REAL(DP) :: dia2, dist2, tau, pp, tt, x, y, z, x1, y1, z1, xn, yn, zn
  REAL(DP) :: absocone(maxnang)
  INTEGER, SAVE :: ifilecount = 0 ! counter for output times

  DO jj = 1, maxnang
    absocone(jj) = 0D0
  END DO

  IF (ifilecount == 0) WRITE (47, '(a)') &
    '# time theta phi emission yield'
  WRITE (47, '(a,i8)') '# printed slot nr.', ifilecount

  jj = 1
  DO itheta = 1, nangtheta
    tt = (angthetah - angthetal)/(nangtheta - 1)*(itheta - 1) + angthetal
    DO iphi = 1, nangphi
      jj = jj + 1
      pp = (angphih - angphil)/(nangphi - 1)*(iphi - 1) + angphil
      xn = COS(pp)*SIN(tt)
      yn = SIN(pp)*SIN(tt)
      zn = COS(tt)
      ind = 0
      DO iz = minz, maxz
        z1 = (iz - nzsh)*dz
        DO iy = miny, maxy
          y1 = (iy - nysh)*dy
          DO ix = minx, maxx
            x1 = (ix - nxsh)*dx
            ind = ind + 1
            tau = xn*(x1 - xango) + yn*(y1 - yango) + &
                  zn*(z1 - zango)
            IF (tau >= 0D0) THEN
              x = xango + tau*xn
              y = yango + tau*yn
              z = zango + tau*zn
              dist2 = (x - x1)*(x - x1) + (y - y1)*(y - y1) + (z - z1)*(z - z1)
              dia2 = TAN(delomega/2)
              dia2 = dia2*dia2
              dia2 = dia2*((x - xango)*(x - xango) + &
                           (y - yango)*(y - yango) + (z - zango)*(z - zango))
              IF (dist2 <= dia2) THEN ! Hit! Grid point is IN cone
                absocone(jj) = absocone(jj) + rhoabso(ind)
              END IF
            END IF
          END DO
        END DO
      END DO
      WRITE (47, '(1f14.5,2f12.3,1e17.5,i5)') tfs, tt/pi*180.0D0, &
        pp/pi*180.0D0, absocone(jj),jj
      CALL flush (47)
    END DO
    WRITE (47, '(1x)') ! empty space TO enable gnuplot 3D plotting
  END DO
  WRITE (47, '(1x)') ! empty space TO enable gnuplot 3D plotting
  ifilecount = 1 + ifilecount
  RETURN
END SUBROUTINE angabso
!------------------------------------------------------------
