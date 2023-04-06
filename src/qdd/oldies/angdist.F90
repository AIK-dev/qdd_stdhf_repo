#include "define.h"

PROGRAM angdist
  USE params
  IMPLICIT REAL(DP) (A - H, O - Z)
! plots the angular distribution of ionization as density plot;
! requires pescmask FILE

! parameters for plotting:
  phi1 = 0. ! IN degrees
  phi2 = 360.
  theta1 = 0.
  theta2 = 180.
  nstepsphi = 40
  nstepstheta = 80.
  deltaomega = 10.
! deltaomega=180./nstepstheta
  WRITE (6, *) 'dOmega: ', deltaomega

  xo = 0.
  yo = 0.
  zo = 0.

  dx = 0.8
  dy = 0.8
  dz = 0.8

! READ IN mask FILE

  DO ind = 1, kdfull2
    READ (5, *) dum, dum, dum, rhoabso(ind)
  END DO

  WRITE (6, *) kdfull2
  WRITE (6, *) minx, miny, minz
  WRITE (6, *) maxx, maxy, maxz
  WRITE (6, *) nxsh, nysh, nzsh
  WRITE (6, *) dx, dy, dz

  OPEN (28, STATUS='unknown', FILE='pangdist.res')

  DO nth = 1, nstepstheta + 1
    theta = (nth - 1)*(theta2 - theta1)/nstepstheta/180.*pi

    DO nph = 1, nstepsphi + 1
      phi = (nph - 1)*(phi2 - phi1)/nstepsphi/180.*pi

      xni = COS(phi)*SIN(theta)
      yni = SIN(phi)*SIN(theta)
      zni = COS(theta)

      acc = 0.0
      ind = 0
      DO iz = minz, maxz
        z1 = (iz - nzsh)*dz

        DO iy = miny, maxy
          y1 = (iy - nysh)*dy

          DO ix = minx, maxx
            x1 = (ix - nxsh)*dx

            ind = ind + 1

            rxrel = x1 - xo
            ryrel = y1 - xo
            rzrel = z1 - xo

            para = xni*rxrel + yni*ryrel + zni*rzrel

            dist2 = (rxrel - para*xni)**2 + (ryrel - para*yni)**2 + &
                    (rzrel - para*zni)**2

            rkegel2 = para*para*TAN(deltaomega)**2

! direc = (xni*rxrel+yni*ryrel+zni*rzrel)

            IF (dist2 < rkegel2 .AND. para >= 0.) THEN
              acc = acc + rhoabso(ind)
            END IF

          END DO
        END DO
      END DO

      acc = acc*dx*dy*dz

! WRITE(6,*) '__________________________________'
! WRITE(6,'(2f10.3)') phi,theta
! WRITE(6,'(3f10.3,1e17.7)') xni,yni,zni,acc

      WRITE (28, '(2f17.4,1e17.7)') phi/pi*180., theta/pi*180., acc

    END DO

! WRITE(28,*)
! one free line TO separate blocks for pm3d IN gnuplot
    WRITE (28, *)

  END DO

  CLOSE (28)

  STOP " finished"
END PROGRAM angdist
