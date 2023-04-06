
!-----savings----------------------------------------------------

SUBROUTINE savings(psi, it)

! Check STATUS of computing resources and optionally SAVE wavefunctions.
!
! Input:
! psi = set of COMPLEX s.p. wavefunctions
! it = time step IN calling routine

  USE params
  IMPLICIT NONE
#if(mpi)
  INCLUDE 'mpif.h'
#endif

  COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)
  INTEGER, INTENT(IN) :: it

  REAL(DP) :: time_act, time_elapse

! IF walltime has almost expired, SAVE all relevant DATA TO
! prepare for restart and STOP:

  IF (trequest > 0D0) THEN
    IF (myn == 0) THEN ! cPW
      CALL cpu_time(time_act)
      time_elapse = time_act - time_absinit
      WRITE (6, '(a,5(1pg13.5))') &
        ' cpu_time: trequest,timefrac,t_elapse=', trequest, timefrac, time_elapse
      CALL FLUSH (6)
    END IF
#if(mpi)
    CALL pi_scatter(time_elapse)
#endif
    IF (time_elapse > trequest*timefrac) THEN
      CALL SAVE(psi, it, outnam)
      IF (myn == 0) THEN
        OPEN (660, STATUS='unknown', FILE='progstatus')
        WRITE (660, *) '0'
        CLOSE (660)
      END IF
#if(mpi)
      CALL mpi_barrier(mpi_comm_world, mpi_ierror)
      !CALL mpi_finalize (mpi_ierror)
#endif
      STOP ' finish at TREQUEST'
    END IF
  END IF ! cPW
! IF PROGRAM has received user message TO shut down PROGRAM, make it so!
  IF (ishutdown /= 0) THEN
    CALL SAVE(psi, it, outnam)
    STOP 'Terminated by user message.'
  END IF

  RETURN
END SUBROUTINE savings

SUBROUTINE escmask(it)

! PRINT collected information on escaping electrons.
!
! Input:
! it = nr. of time step IN calling routine.

  USE params
  USE util, ONLY: inttostring, printfield
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: it

  INTEGER :: nbe, nbeabs

  IF (jescmaskorb /= 0 .AND. MOD(it, jescmaskorb) == 0) THEN
    DO nbe = 1, nstate
      nbeabs = nrel2abs(nbe)
      IF (nbeabs < 1000) THEN
        OPEN (588, STATUS='unknown', FILE='pescmaskOrb.'//trim(adjustl(inttostring(nbeabs)))//'.'//outnam)
      ELSE
        STOP 'ERROR: Too many states for jescmaskOrb'
      END IF
      CALL printfield(588, rhoabsoorb(1, nbe), 'x')
      CLOSE (588)
    END DO
  END IF

  IF (myn == 0) THEN
    IF (jescmask .NE. 0 .AND. MOD(it, jescmask) == 0) THEN
      IF (it < 1000000000) THEN
        OPEN (589, STATUS='unknown', FILE='pescmask.'//trim(adjustl(inttostring(it)))//'.'//outnam)
      ELSE
        STOP '::too many time steps::'
      END IF
      CALL printfield(589, rhoabso, 'x')
      CLOSE (589)
    END IF

    IF (it == 0) THEN
      OPEN (589, STATUS='unknown', FILE='pescmask.0.'//outnam) ! Why DO it again ??
      CALL printfield(589, rhoabso, 'x')
      CLOSE (589)
    END IF

  END IF

  RETURN
END SUBROUTINE escmask

!-----init_scattel-------------------------------------------------

SUBROUTINE init_scattel(psi)

! Adds one electron IN scattering state as Gaussian wavepacket
! with a certain velocity given by 'scatterelectronv...'.
! The wavefunctions bookkeeping fields are extended accordingly.

  USE params
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)

  INTEGER :: i, ind, ix, iy, iz
  REAL(DP) :: arg, fac, fr, pnorm, rr, x1, y1, z1
!------------------------------------------------------------------

  pnorm = scatterelectronvxn**2 + scatterelectronvyn**2 + scatterelectronvzn**2
  pnorm = SQRT(pnorm)

  IF (pnorm == 0D0) STOP 'Momentum of Scatt. Electron vanishes!'

  scatterelectronvxn = scatterelectronvxn/pnorm
  scatterelectronvyn = scatterelectronvyn/pnorm
  scatterelectronvzn = scatterelectronvzn/pnorm

  pnorm = SQRT(scatterelectronenergy*2D0*ame)

  scatterelectronvxn = scatterelectronvxn*pnorm
  scatterelectronvyn = scatterelectronvyn*pnorm
  scatterelectronvzn = scatterelectronvzn*pnorm

  nstate = nstate + 1
  IF (nstate > kstate) STOP ' insufficient KSTATE IN INIT_SCATTEL'
  occup(nstate) = 1D0
  nrel2abs(nstate) = nstate
  nabs2rel(nstate) = nstate
  ispin(nstate) = 1

! test actual NUMBER of active states, ONLY occcupied states allowed
  DO i = 1, nstate
    IF (occup(i) < 0.5D0) STOP "ONLY occupied states allowed IN CASE of attachement"
  END DO

  fac = 1D0/SQRT(pi**1.5D0*scatterelectronw**3)

  ind = 0

  DO iz = minz, maxz
    z1 = (iz - nzsh)*dz
    DO iy = miny, maxy
      y1 = (iy - nysh)*dy
      DO ix = minx, maxx
        x1 = (ix - nxsh)*dx
        ind = ind + 1
        arg = (x1 - scatterelectronx)*scatterelectronvxn + &
              (y1 - scatterelectrony)*scatterelectronvyn + &
              (z1 - scatterelectronz)*scatterelectronvzn
        rr = (x1 - scatterelectronx)**2 &
             + (y1 - scatterelectrony)**2 &
             + (z1 - scatterelectronz)**2
        fr = fac*EXP(-rr/2D0/scatterelectronw**2)
        psi(ind, nstate) = CMPLX(fr, 0D0, DP)*CMPLX(COS(arg), SIN(arg), DP)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE init_scattel

!-----init_projwf-------------------------------------------------

SUBROUTINE init_projwf(psi)

! Boosts the electronic wavefunctions of the projectile
! by the same velocity, whose norm is given by 'eproj' and
! the direction by 'vpx', 'vpy',and 'vpz'.
! This routine is similar TO 'init_velocel', but applied
! ONLY TO projectile wavefunctions, distinguished by an
! ENTRY IN the array 'proj_states'.

  USE params
  IMPLICIT NONE

  COMPLEX(DP), INTENT(OUT) :: psi(kdfull2, kstate)

  INTEGER :: ind, ix, iy, iz, kk, nbe, nbee
  REAL(DP) :: arg, rnorm, v0, x1, y1, z1
  COMPLEX(DP) :: cfac

!------------------------------------------------------------------
! lionel : np(nion) ==> np(nproj)

  v0 = SQRT(2D0*eproj/(amu(np(nproj))*1836.0D0*ame))
  rnorm = vpx**2 + vpy**2 + vpz**2
  rnorm = SQRT(rnorm)

  IF (rnorm == 0) STOP 'Velocity vector not normalizable'

  vpx = vpx/rnorm*v0*ame
  vpy = vpy/rnorm*v0*ame
  vpz = vpz/rnorm*v0*ame

  IF (taccel > 0D0) RETURN ! boost is done adiabatically

  IF (.NOT. init_ao) THEN
    WRITE (*, *) ' instantaneous acceleration ONLY for INIT_AO==.TRUE.'
    STOP ' instantaneous acceleration ONLY for INIT_AO==.TRUE.'
  END IF
  IF (nproj_states == 0) THEN
    WRITE (*, *) ' CAUTION : atomic projectile without wf TO boost'
    WRITE (*, *) ' IF there are electrons on the projectile, please USE'
    WRITE (*, *) ' nproj_states IN GLOBAL, proj_states and nproj IN DYNAMIC'
  ELSE
    WRITE (*, *) 'Input states of the projectile', proj_states(:)
    ind = 0
    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz
      DO iy = miny, maxy
        y1 = (iy - nysh)*dy
        DO ix = minx, maxx
          x1 = (ix - nxsh)*dx
          ind = ind + 1
          arg = x1*vpx + y1*vpy + z1*vpz
          cfac = CMPLX(COS(arg), SIN(arg), DP)
          DO nbe = 1, nstate
            nbee = nrel2abs(nbe)
            DO kk = 1, nproj_states
              IF (nbee == proj_states(kk)) psi(ind, nbe) = cfac*psi(ind, nbe)
            END DO
          END DO
        END DO
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE init_projwf

!-----projmoms-------------------------------------------------------projmoms
! REAL version
!-----------------------------------------------------------------------
SUBROUTINE r_projmoms(rho, psi)

! Multipole moments relative TO c.m. on 'qetarget' and relative TO
! projectile coordinate 'cz' on 'qeproj'.

  USE params
  USE util, ONLY: getcm
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: rho(2*kdfull2)
  REAL(DP), INTENT(IN) :: psi(kdfull2, kstate)
#if(mpi)
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)
#endif

  INTEGER :: ind, ix, iy, iz, k
  REAL(DP) :: sproj, starget
  REAL(DP) :: x1, y1, z1, x1t, y1t, z1t, x1p, y1p, z1p

#if(nompi)
  INTEGER :: ik, ikk
#else
  INTEGER :: kk, nbe, nbee
  REAL(DP) :: sprojec
#endif
!----------------------------------------------------------------------
  nrmom = 35
  IF (nrmom > kmom) STOP ' too many moments IN projmoms'

  DO k = 1, nrmom
    qetarget(k) = 0D0
    qeproj(k) = 0D0
  END DO

! switch for calculating moments relative TO center of mass (1)
! or center of box (0)
  rvectmp = 0D0
  IF (tmoms_rel_cm .AND. nion2 > 0) CALL getcm(1, 0, 0)

  ind = 0
  DO iz = minz, maxz
    z1 = (iz - nzsh)*dz
    z1t = z1 - rvectmp(3)
    z1p = z1 - cz(nproj)
    DO iy = miny, maxy
      y1 = (iy - nysh)*dy
      y1t = y1 - rvectmp(2)
      y1p = y1 - cy(nproj)
      DO ix = minx, maxx
        ind = ind + 1
        IF ((ix <= nx2) .AND. (iy <= ny2) .AND. (iz <= nz2)) THEN
          x1 = (ix - nxsh)*dx
          x1t = x1 - rvectmp(1)
          x1p = x1 - cx(nproj)
          sproj = 0D0
#if(nompi)
          DO ik = 1, nproj_states
            ikk = proj_states(ik)
            sproj = sproj + psi(ind, ikk)*psi(ind, ikk)
          END DO
#else
          sprojec = 0D0
          DO nbe = 1, nstate
            nbee = nrel2abs(nbe)
            DO kk = 1, nproj_states
              IF (nbee == proj_states(kk)) THEN
                sprojec = sprojec + psi(ind, nbe)*psi(ind, nbe)
              END IF
            END DO
          END DO
          CALL mpi_barrier(mpi_comm_world, mpi_ierror)
          CALL mpi_allreduce(sprojec, sproj, 1, mpi_double_precision, &
                             mpi_sum, mpi_comm_world, mpi_ierror)
          CALL mpi_barrier(mpi_comm_world, mpi_ierror)
#endif
          starget = rho(ind) - sproj
! monopole
          qetarget(1) = qetarget(1) + starget
          qeproj(1) = qeproj(1) + sproj
! dipole
          qetarget(2) = qetarget(2) + starget*x1t
          qetarget(3) = qetarget(3) + starget*y1t
          qetarget(4) = qetarget(4) + starget*z1t

          qeproj(2) = qeproj(2) + sproj*x1p
          qeproj(3) = qeproj(3) + sproj*y1p
          qeproj(4) = qeproj(4) + sproj*z1p
        END IF
      END DO
    END DO
  END DO

  DO k = 1, nrmom
    qetarget(k) = qetarget(k)*dvol
    qeproj(k) = qeproj(k)*dvol
  END DO

  DO k = 2, nrmom
    qetarget(k) = qetarget(k)/qetarget(1) !normalization
    qeproj(k) = qeproj(k)/qeproj(1) !normalization
  END DO

  RETURN
END SUBROUTINE r_projmoms

!-----------------------------------------------------------------------
! COMPLEX version
!-----------------------------------------------------------------------
SUBROUTINE c_projmoms(rho, psi)
  USE params
  USE util, ONLY: getcm
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: rho(2*kdfull2)
  COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)
#if(mpi)
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)
#endif

  INTEGER :: ind, ix, iy, iz, k
  REAL(DP) :: sproj, starget
  REAL(DP) :: x1, y1, z1, x1t, y1t, z1t, x1p, y1p, z1p

#if(nompi)
  INTEGER :: ik, ikk
#else
  INTEGER :: kk, nbe, nbee
  REAL(DP) :: sprojec
#endif
!----------------------------------------------------------------------
  nrmom = 35
  IF (nrmom > kmom) STOP ' too many moments IN projmoms'

  DO k = 1, nrmom
    qetarget(k) = 0D0
    qeproj(k) = 0D0
  END DO

! switch for calculating moments relative TO center of mass (1)
! or center of box (0)
  rvectmp = 0D0
  IF (tmoms_rel_cm .AND. nion2 > 0) CALL getcm(1, 0, 0)

  ind = 0
  DO iz = minz, maxz
    z1 = (iz - nzsh)*dz
    z1t = z1 - rvectmp(3)
    z1p = z1 - cz(nproj)
    DO iy = miny, maxy
      y1 = (iy - nysh)*dy
      y1t = y1 - rvectmp(2)
      y1p = y1 - cy(nproj)
      DO ix = minx, maxx
        ind = ind + 1
        IF ((ix <= nx2) .AND. (iy <= ny2) .AND. (iz <= nz2)) THEN
          x1 = (ix - nxsh)*dx
          x1t = x1 - rvectmp(1)
          x1p = x1 - cx(nproj)
          sproj = 0D0
#if(nompi)
          DO ik = 1, nproj_states
            ikk = proj_states(ik)
            sproj = sproj + REAL(CONJG(psi(ind, ikk))*psi(ind, ikk), DP)
          END DO
#else
          sprojec = 0D0
          DO nbe = 1, nstate
            nbee = nrel2abs(nbe)
            DO kk = 1, nproj_states
              IF (nbee == proj_states(kk)) THEN
                sprojec = sprojec + REAL(CONJG(psi(ind, nbe))*psi(ind, nbe), DP)
              END IF
            END DO
          END DO
          CALL mpi_barrier(mpi_comm_world, mpi_ierror)
          CALL mpi_allreduce(sprojec, sproj, 1, mpi_double_precision, &
                             mpi_sum, mpi_comm_world, mpi_ierror)
          CALL mpi_barrier(mpi_comm_world, mpi_ierror)
#endif
          starget = rho(ind) - sproj
! monopole
          qetarget(1) = qetarget(1) + starget
          qeproj(1) = qeproj(1) + sproj
! dipole
          qetarget(2) = qetarget(2) + starget*x1t
          qetarget(3) = qetarget(3) + starget*y1t
          qetarget(4) = qetarget(4) + starget*z1t

          qeproj(2) = qeproj(2) + sproj*x1p
          qeproj(3) = qeproj(3) + sproj*y1p
          qeproj(4) = qeproj(4) + sproj*z1p
        END IF
      END DO
    END DO
  END DO

  DO k = 1, nrmom
    qetarget(k) = qetarget(k)*dvol
    qeproj(k) = qeproj(k)*dvol
  END DO

  DO k = 2, nrmom
    qetarget(k) = qetarget(k)/qetarget(1) !normalization
    qeproj(k) = qeproj(k)/qeproj(1) !normalization
  END DO

  RETURN
END SUBROUTINE c_projmoms

! **************************

SUBROUTINE mtv_fld(field, i)

! Prepares plotting of 'field' with MTV software.

  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i ! 1=static, 2=dynamic (ONLY change is filenaming)
  REAL(DP), INTENT(IN) :: field(kdfull2)

  INTEGER :: ind, jx, jy, num

! PRINT 2d-field

  num = nz
  ind = nxyf*num
  IF (i == 1) THEN
    OPEN (70, FILE='denstat.'//outnam, STATUS='unknown')
  ELSE IF (i == 2) THEN
    OPEN (70, FILE='densdyn.'//outnam, STATUS='unknown')
  END IF
  WRITE (70, *) '$ DATA=contour'
  WRITE (70, '(a,i4,a,f7.2,a,f7.2)') &
    '% nx =', nx2, ' xmin=', (minx - nxsh)*dx, ' xmax=', (maxx - nxsh)*dx
  WRITE (70, '(a,i4,a,f7.2,a,f7.2)') &
    '% ny =', ny2, ' ymin=', (miny - nysh)*dy, ' ymax=', (maxy - nysh)*dy
  WRITE (70, *) '% contstyle=2'
  WRITE (70, *) '% nsteps=50'
  WRITE (70, *) '% interp=2'
  WRITE (70, *) '% xlabel="x"'
  WRITE (70, *) '% ylabel="y"'
  WRITE (70, *) '% toplabel="mg2"'
  WRITE (70, *) '% equalscale=false'
  WRITE (70, *) '% fitpage=false'
  WRITE (70, *) '% xyratio=1.0'

  DO jy = miny, maxy
    DO jx = minx, maxx
      ind = ind + 1
      WRITE (70, '(g13.5)') field(ind)
    END DO
    WRITE (70, *) ' '
  END DO

  WRITE (70, *) '$ END'
  CLOSE (70)

  RETURN
END SUBROUTINE mtv_fld

! ****************************

SUBROUTINE fmtv_fld(psi, field, i)

! Prepares plotting of 'psi**2' with MTV software.

  USE params
  IMPLICIT NONE
  CHARACTER(LEN=6) :: ext
  CHARACTER(LEN=8) :: ext1

  REAL(DP), INTENT(IN) :: field(kdfull2)
  COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)
  INTEGER, INTENT(IN) :: i

  INTEGER :: ind, ishift, j, jx, jy, jz, k, nb
  REAL(DP) :: rhoz, xmin1, ymin1, zmin1, xmax1, ymax1, zmax1
  REAL(DP), ALLOCATABLE :: rho(:)

  ALLOCATE (rho(2*kdfull2))

! PRINT integrated 3d-field (DIMENSION TO integrate on depends of values of i3dx, i3dz)

  xmin1 = -nx*dx
  ymin1 = -ny*dy
  zmin1 = -nz*dz
  xmax1 = nx*dx
  ymax1 = ny*dy
  zmax1 = nz*dz
  WRITE (ext, '(a,i5.5)') '.', i

  IF (i3dz == 1) THEN
    OPEN (70, FILE='zdensdyn'//ext//'.'//outnam, STATUS='unknown')
    WRITE (70, '(a)') '$ DATA=contour'
    WRITE (70, '(a,i2,a,f7.3,a,f5.2)') '% nx = ', nx2, &
      ' xmin=', xmin1, ' xmax=', xmax1
    WRITE (70, '(a,i2,a,f7.3,a,f5.2)') '% ny = ', ny2, &
      ' ymin=', ymin1, ' ymax=', ymax1
    WRITE (70, '(a)') '% contstyle=2'
    WRITE (70, '(a)') '% nsteps=50'
    WRITE (70, '(a)') '% interp=2'
    WRITE (70, '(a)') '% xlabel="x"'
    WRITE (70, '(a)') '% ylabel="y"'
    WRITE (70, '(a,a3,a)') '% toplabel="', outnam, '"'
    WRITE (70, '(a)') '% equalscale=false'
    WRITE (70, '(a)') '% fitpage=false'
    WRITE (70, '(a)') '% xyratio=1.0'
    ind = 0
    DO jy = miny, maxy
      DO jx = minx, maxx
        ind = ind + 1
        j = 0
        rhoz = 0D0
        DO jz = minz, maxz
          j = j + 1
          k = (j - 1)*nxyf + ind
          rhoz = rhoz + field(k)
        END DO
        WRITE (70, '(g13.5)') rhoz
      END DO
      WRITE (70, *) ' '
    END DO
    WRITE (70, '(a)') '$ END'
    CLOSE (70)
  END IF

  IF (i3dx == 1) THEN
    OPEN (73, FILE='xdensdyn'//ext//'.'//outnam, STATUS='unknown')
    WRITE (73, '(a)') '$ DATA=contour'
    WRITE (73, '(a,i2,a,f7.3,a,f5.2)') '% nx = ', ny2, &
      ' xmin=', ymin1, ' xmax=', ymax1
    WRITE (73, '(a,i2,a,f7.3,a,f5.2)') '% ny = ', nz2, &
      ' ymin=', zmin1, ' ymax=', zmax1
    WRITE (73, '(a)') '% contstyle=2'
    WRITE (73, '(a)') '% nsteps=50'
    WRITE (73, '(a)') '% interp=2'
    WRITE (73, '(a)') '% xlabel="y"'
    WRITE (73, '(a)') '% ylabel="z"'
    WRITE (73, '(a,a3,a)') '% toplabel="', outnam, '"'
    WRITE (73, '(a)') '% equalscale=false'
    WRITE (73, '(a)') '% fitpage=false'
    WRITE (73, '(a)') '% xyratio=1.0'

    ind = 0
    DO jz = minz, maxz
      DO jy = miny, maxy
        rhoz = 0D0
        DO jx = minx, maxx
          ind = ind + 1
          rhoz = rhoz + field(ind)
        END DO
        WRITE (73, '(g13.5)') rhoz
      END DO
      WRITE (73, *) ' '
    END DO
    WRITE (73, '(a)') '$ END'
    CLOSE (73)
  END IF

  IF (i3dstate == 1) THEN
    DO nb = 1, nstate
      ishift = (ispin(nrel2abs(nb)) - 1)*nxyz
      DO ind = 1, nxyz
        rho(ind + ishift) = occup(nb)*(CONJG(psi(ind, nb)))*psi(ind, nb)
      END DO
      WRITE (ext1, '(a1,i1,a1,i5.5)') '.', nb, '.', i
      OPEN (71, FILE='zdenspsi'//ext1//'.'//outnam, STATUS='unknown')
      WRITE (71, '(a)') '$ DATA=contour'
      WRITE (71, '(a,i2,a,f7.3,a,f5.2)') '% nx = ', nx2, &
        ' xmin=', xmin1, ' xmax=', xmax1
      WRITE (71, '(a,i2,a,f7.3,a,f5.2)') '% ny = ', ny2, &
        ' ymin=', ymin1, ' ymax=', ymax1
      WRITE (71, '(a)') '% contstyle=2'
      WRITE (71, '(a)') '% nsteps=50'
      WRITE (71, '(a)') '% interp=2'
      WRITE (71, '(a)') '% xlabel="x"'
      WRITE (71, '(a)') '% ylabel="y"'
      WRITE (71, '(a,a3,a,i1,a)') '% toplabel="', outnam, '_', nb, '"'
      WRITE (71, '(a)') '% equalscale=false'
      WRITE (71, '(a)') '% fitpage=false'
      WRITE (71, '(a)') '% xyratio=1.0'
      ind = 0
      DO jy = miny, maxy
        DO jx = minx, maxx
          ind = ind + 1
          j = 0
          rhoz = 0D0
          DO jz = minz, maxz
            j = j + 1
            k = (j - 1)*nxyf + ind
            rhoz = rhoz + rho(k)
          END DO
          WRITE (71, '(g13.5)') rhoz
        END DO
        WRITE (71, *) ' '
      END DO
      WRITE (71, '(a)') '$ END'
      CLOSE (71)

      OPEN (72, FILE='xdenspsi'//ext1//'.'//outnam, STATUS='unknown')
      WRITE (72, '(a)') '$ DATA=contour'
      WRITE (72, '(a,i2,a,f7.3,a,f5.2)') '% nx = ', ny2, &
        ' xmin=', ymin1, ' xmax=', ymax1
      WRITE (72, '(a,i2,a,f7.3,a,f5.2)') '% ny = ', nz2, &
        ' ymin=', zmin1, ' ymax=', zmax1
      WRITE (72, '(a)') '% contstyle=2'
      WRITE (72, '(a)') '% nsteps=50'
      WRITE (72, '(a)') '% interp=2'
      WRITE (72, '(a)') '% xlabel="y"'
      WRITE (72, '(a)') '% ylabel="z"'
      WRITE (72, '(a,a3,a,i1,a)') '% toplabel="', outnam, '_', nb, '"'
      WRITE (72, '(a)') '% equalscale=false'
      WRITE (72, '(a)') '% fitpage=false'
      WRITE (72, '(a)') '% xyratio=1.0'
      ind = 0
      DO jz = minz, maxz
        DO jy = miny, maxy
          rhoz = 0D0
          DO jx = minx, maxx
            ind = ind + 1
            rhoz = rhoz + rho(ind)
          END DO
          WRITE (72, '(g13.5)') rhoz
        END DO
        WRITE (72, *) ' '
      END DO
      WRITE (72, '(a)') '$ END'
      CLOSE (72)
    END DO
  END IF

  DEALLOCATE (rho)

  RETURN
END SUBROUTINE fmtv_fld
! **************************************************

SUBROUTINE mtv_3dfld(rho, iframe)

! Prepares 3D plotting of 'rho' with MTV software.

  USE params

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iframe
  REAL(DP), INTENT(IN) :: rho(2*kdfull2)

  CHARACTER(LEN=5) :: ext
  INTEGER :: ind, ix, iy, iz
!REAL(DP) :: sumr

!sumr=0
  WRITE (ext, '(a,i4.4)') '.', iframe
! OPEN(UNIT=80,STATUS='UNKNOWN',FILE='xdens'//ext//'.'//outnam)
! OPEN(UNIT=81,STATUS='UNKNOWN',FILE='ydens'//ext//'.'//outnam)
  OPEN (UNIT=82, STATUS='UNKNOWN', FILE='dens'//ext//'.'//outnam)
! OPEN(UNIT=83,STATUS='UNKNOWN',FILE='xydens'//ext//'.'//outnam)
  ind = 0
  DO iz = minz, maxz
! sumr=0
    DO iy = miny, maxy
      DO ix = minx, maxx
        ind = ind + 1
! IF (ix.eq.24) THEN
! WRITE(80,*)Rho(ind)
! END IF
! IF (iy.eq.24) THEN
! WRITE(81,*)Rho(ind)
! END IF
! IF (iz.eq.24) THEN
        WRITE (82, *) rho(ind)
! END IF
! sumr=sumr+rho(ind)
      END DO
    END DO
! WRITE(83,*)sumr
  END DO
! CLOSE(80)
! CLOSE(81)
  CLOSE (82)
! CLOSE(83)

  RETURN
END SUBROUTINE mtv_3dfld

SUBROUTINE hpsi_boostinv(qact, aloc, current, rho, nbe)

! ACTION of boost-invariant Hamiltonian on one s.p. wavefunction:
! qact = wavefunction on which H acts and resulting w.f.
! aloc = local potential for the actual spin component
! current = average local momentum IN x,y, and z-direction
! nbe = NUMBER of state
! The routine requires that 'current' has been accumulated before.

  USE params
  USE kinetic
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT) :: qact(kdfull2)
  REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
  REAL(DP), INTENT(IN) :: current(kdfull2, 3)
  REAL(DP), INTENT(IN) :: rho(2*kdfull2)
  INTEGER, INTENT(IN) :: nbe

! workspaces
  COMPLEX(DP), ALLOCATABLE :: q1(:), q2(:)

  COMPLEX(DP) :: cf

!----------------------------------------------------------------------

! check availability of FFT
  IF (.NOT. ALLOCATED(akv)) STOP "HPSI_BOOSTINVARIANT requires FFT"
  ALLOCATE (q1(kdfull2), q2(kdfull2))

  q1 = qact
  CALL hpsi(qact, aloc, nbe, 0)

  CALL xgradient_rspace(q1, q2)
  qact = qact - h2m*current(:, 1)*q2
  q2 = current(:, 1)*q1
  CALL xgradient_rspace(q2, q2)
  qact = qact - h2m*q2

  CALL ygradient_rspace(q1, q2)
  qact = qact - h2m*current(:, 2)*q2
  q2 = current(:, 2)*q1
  CALL ygradient_rspace(q2, q2)
  qact = qact - h2m*q2

  CALL zgradient_rspace(q1, q2)
  qact = qact - eye*h2m*current(:, 3)*q2
  q2 = current(:, 3)*q1
  CALL zgradient_rspace(q2, q2)
  qact = qact - eye*h2m*q2

  qact = qact + h2m* &
         (current(:, 1)**2 + current(:, 2)**2 + current(:, 3)**2)

  WRITE (*, *) ' HPSI_BOOSTINV: E_coll=', dvol*h2m*SUM(rho(:)* &
                                                       (current(:, 1)**2 + current(:, 2)**2 + current(:, 3)**2))

  RETURN
END SUBROUTINE hpsi_boostinv

!-----afterburn---------------------------------------------------------

SUBROUTINE afterburn(psir, rho, aloc)

! Improving static solution by doing a couple of steps with
! imaginary-time propatation.
!
! Input/Output:
! psir = REAL wavefunctions
! rho = electron density
! aloc = local mean-field potential

  USE params

#if(fsic)
  USE twost, ONLY: tnearest, init_fsic, init_vecs, end_fsic, expdabvol_rotate_init
#endif

  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: psir(kdfull2, kstate)
  REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)
  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)

  COMPLEX(DP), ALLOCATABLE :: psiw(:, :), psi(:, :)

  INTEGER :: it, ion

! ****************************************************
! imaginary-time iteration (TO improve statics)
! ****************************************************

  WRITE (*, *) ' start afterburn. isitmax=', isitmax
  CALL flush (6)
  ALLOCATE (psi(kdfull2, kstate))
#if(fsic)
  IF (ifsicp > 7 .AND. nelect > 0) CALL init_fsic()
  IF (ifsicp >= 8) CALL expdabvol_rotate_init
#endif
  CALL restart2(psi, outnam, .true.) ! READ static wf's
#if(fsic)
  IF (ifsicp > 7) CALL init_vecs()
#endif
  CALL calcpseudo()
  CALL calclocal(rho, aloc) ! ??
  IF (ifsicp > 0) CALL calc_sic(rho, aloc, psi) ! USE some TYPE of SIC
  IF (ipsptyp == 1) THEN ! full Goedecker pseudopotentials require projectors
  DO ion = 1, nion
  IF (tinterpol) THEN
    CALL calc_projFine(cx(ion), cy(ion), cz(ion), cx(ion), cy(ion), cz(ion), ion)
! CALL calc_proj(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion) !???
  ELSE
    CALL calc_proj(cx(ion), cy(ion), cz(ion), cx(ion), cy(ion), cz(ion), ion)
  END IF
  END DO
  IF (tinterpol) CALL mergetabs
  END IF
  CALL info(psi, rho, aloc, -1)

  IF (ifexpevol) ALLOCATE (psiw(kdfull2, kstate))

  ! Start imaginary time iterations
  DO it = 1, isitmax
    WRITE (*, *) ' afterburn. iteration=', it
    IF (ifexpevol) THEN
      CALL tstep_exp(psi, aloc, rho, it, psiw, .TRUE.)
      CALL cschmidt(psi) ! schmidt normalization of results
    ELSE
      STOP 'imaginary-time step requires exponential evolution'
      CALL tstep(psi, aloc, rho, it)
      IF (it == isitmax) THEN
        psi = CMPLX(REAL(psi, DP), 0D0, DP)
        CALL cschmidt(psi) ! schmidt normalization of results
      END IF
    END IF
    IF (MOD(it, istinf) == 0) CALL info(psi, rho, aloc, it)
  END DO
  DEALLOCATE (psiw)
  CALL info(psi, rho, aloc, -1)
  CALL SAVE(psi, -1, outnam)
  DEALLOCATE (psi)
#if(fsic)
  IF (nelect > 0 .AND. ifsicp > 7) CALL end_fsic()
#endif
  IF (itmax <= 0) STOP ' terminate with afterburn '

  RETURN

END SUBROUTINE afterburn

!-----hstate---------------------------------------------------------

SUBROUTINE hstate()

! Cuts a hole into one occupied state 'nhstate' by resetting
! 'occup(nhstate).

  USE params

  IMPLICIT NONE

  REAL(DP) :: phangpi, ca, sa, cb, sb

  WRITE (*, *) 'HSTATE entered:', npstate, nhstate, occup(nhstate)

  IF (npstate .NE. 0) STOP ' HSTATE: particle state must be inactive (=0)'
  IF (nhstate > nstate) STOP ' HSTATE: hole state OUT of range'
  IF (occup(nhstate) < 0.999D0) &
    STOP 'HSTATE: hole state not sufficiently occupied'

  occup(nhstate) = 0D0

! reset particle numbers
  nelect = nelect - 1
  IF (ispin(nhstate) == 2) nspdw = nspdw - 1
  nhstate = 0 ! hinder second CALL

  WRITE (6, '(a,2i4,6(1pg13.5))') ' h state defined: nhstate,nelect,nspdw', &
    nhstate, nelect, nspdw
  WRITE (6, '(5(1pg13.5))') occup(1:nstate)
  CALL flush (6)

END SUBROUTINE hstate

!-----phstate---------------------------------------------------------

SUBROUTINE phstate(psi)

! Picks one particular 1ph state from the given wavefunctions
! 'psi' and mixes this TO a rotated ph configuration. The state
! is rotated by 'phangle'.
! The particle state is 'npstate' and the hole state 'nhstate'.
! The old hole wavefunction is stored and returned on 'oldhole'.
! The parameters 'phangle', 'npstate' and 'nhstate' are communicated
! through MODULE PARAMS.
! 'oldhole' is also communicated through PARAMS and allocated here.

  USE params
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)

  REAL(DP) :: phangpi, ca, phphasepi
  COMPLEX(DP) :: sa

 WRITE(*,*) 'PHSTATE entered:',npstate,nhstate,occup(npstate),occup(nhstate),phangle,phphase

  IF (npstate > nstate) STOP ' PHSTATE: particle state OUT of range'
  IF (nhstate > nstate) STOP ' PHSTATE: hole state OUT of range'
  IF (occup(npstate) > 0.5D0) STOP 'PHSTATE: particle state already occupied'
  IF (occup(nhstate) < 0.5D0) STOP 'PHSTATE: hole state not occupied'

  IF(.NOT.ALLOCATED(oldhole)) ALLOCATE (oldhole(kdfull2))
  IF(.NOT.ALLOCATED(newhole)) ALLOCATE (newhole(kdfull2))

  oldhole = psi(:, nhstate)
  phangpi = phangle/180D0*PI
  phphasepi = phphase/180D0*PI
  ca = cos(phangpi)
  sa = sin(phangpi)*EXP(CMPLX(0D0, phphasepi, DP))
!  cb = cos(phphase/180D0*PI)
!  sb = sin(phphase/180D0*PI)

!  psi(:, nhstate) = (ca*cb + CMPLX(0D0, sa*sb, DP))*psi(:, nhstate) &
!                    + (sa*cb - CMPLX(0D0, ca*sb, DP))*psi(:, npstate)
!  psi(:, npstate) = -(sa*cb + CMPLX(0D0, ca*sb, DP))*oldhole &
!                    + (ca*cb - CMPLX(0D0, sa*sb, DP))*psi(:, npstate)
  psi(:, nhstate) = ca*psi(:, nhstate) + sa*psi(:, npstate)
  psi(:, npstate) = CONJG(sa)*oldhole + ca*psi(:, npstate)

  newhole = psi(:, nhstate)       !! ???????????

  WRITE (6, '(a,2i4,6(1pg13.5))') ' 1ph state mixed:', nhstate, npstate, &
    phangle, phphase, ca, sa
  CALL flush (6)

END SUBROUTINE phstate

!-----phstate_static---------------------------------------------------------

SUBROUTINE phstate_static()

! Creates a 1ph state by exchanging the occupation numbers between
! the states 'npstate' <--> 'nhstate' in array 'occup' (in 'params').

  USE params
  IMPLICIT NONE


  REAL(DP) :: sav

  IF (npstate > nstate) STOP ' PHSTATE_STATIC: particle state OUT of range'
  IF (nhstate > nstate) STOP ' PHSTATE_STATIC: hole state OUT of range'
  IF (occup(npstate) > 0.5D0) STOP 'PHSTATE: particle state already occupied'
  IF (occup(nhstate) < 0.5D0) STOP 'PHSTATE: hole state not occupied'

  sav = occup(npstate)
  occup(npstate) = occup(nhstate)
  occup(nhstate) = sav

  WRITE (6, '(a,2i4)') ' static 1ph state generated:',nhstate, npstate

  CALL flush (6)

END SUBROUTINE phstate_static

!-----phoverl---------------------------------------------------------

SUBROUTINE phoverl(psi)

! Computes overlap of original Slater state with dynamically rotated
! 1ph state. This routine makes sense ONLY after having called
! 'phstate'.

  USE params
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)
  LOGICAL, SAVE :: tfirst = .true.
  COMPLEX(DP) :: covo, covn

  IF (TFIRST) THEN
    OPEN (381, FILE='phoverlap.'//outnam)
    WRITE (381, '(a/a,2f10.2,2i5/a)') &
      '# overlap of ph-rotated state with orinigal hole state:', &
      '# phangle,phphase,npstate,nhstate=', phangle, phphase, npstate, nhstate, &
      '# time (fs) REAL(ov) imag(ov) abs(ov) phase(ov)'
    tfirst = .false.
  END IF

  covo = SUM(CONJG(oldhole)*psi(:, nhstate))*dvol
  covn = SUM(CONJG(newhole)*psi(:, nhstate))*dvol
  WRITE (381, '(f10.3,8(1pg13.5))') tfs, &
    covo, SQRT(REAL(covo, DP)**2 + AIMAG(covo)**2), atan2(AIMAG(covo), REAL(covo, DP)), &
    covn, SQRT(REAL(covn, DP)**2 + AIMAG(covn)**2), atan2(AIMAG(covn), REAL(covn, DP))
  CALL FLUSH (381)

  RETURN

END SUBROUTINE phoverl

!-----view3D------------------------------------------------------

  SUBROUTINE view3d()

! Prepares FILE for viewing ionic structure with 'xbs'.

    USE params
    IMPLICIT NONE
    CHARACTER(LEN=3) :: ext
    INTEGER :: ion
!-----------------------------------------------------------------

    ext = outnam

    OPEN (39, FILE=ext//'.bs', STATUS='UNKNOWN')
    DO ion = 1, nion
      WRITE (39, '(a4,a6,a1,a6,3f10.4)') 'atom', ' ', 'H', ' ', cx(ion), cy(ion), cz(ion)
    END DO
    WRITE (39, *) ' '
    WRITE (39, '(a4,a6,a1,2(a6,a3))') 'spec', ' ', 'H', ' ', '0.5', ' ', '1.0'
    WRITE (39, *) ' '
    WRITE (39, '(a5,a5,a1,a6,a1,4(a6,a3))') 'bonds', ' ', 'H', ' ', &
      'H', ' ', '0.0', ' ', '9.0', ' ', '0.1', ' ', '1.0'
    WRITE (39, *) ' '
    WRITE (39, '(a4,9(a2,a3))') 'tmat', ' ', '1.0', ' ', '0.0', ' ', '0.0', &
      ' ', '0.0', ' ', '1.0', ' ', '0.0', ' ', '0.0', ' ', '0.0', ' ', '1.0'
    WRITE (39, '(a4,a6,a4)') 'dist', ' ', '12.0'
    WRITE (39, '(a3,a8,a3)') 'inc', ' ', '5.0'
    WRITE (39, '(a5,a5,a4)') 'scale', ' ', '20.0'
    WRITE (39, '(a4,a7,a3)') 'rfac', ' ', '1.0'
    WRITE (39, '(a4,a7,a3)') 'bfac', ' ', '1.0'
    WRITE (39, '(a3,2(a6,a3))') 'pos', ' ', '0.0', ' ', '0.0'
    WRITE (39, '(a8,9(a2,a1))') 'switches', ' ', '1', ' ', '0', ' ', '1', ' ', '0', &
      ' ', '0', ' ', '1', ' ', '1', ' ', '0', ' ', '0'
    CLOSE (39)

    RETURN
  END SUBROUTINE view3d
!-----omega_mieplasmon--------------------------------------------------

  REAL(DP) FUNCTION omega_mieplasmon(rho)

! Estimates the Mie plasmon frequency
! (? seems TO be correct ONLY for Na with soft local PsP ?)

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: rho(kdfull2)

    INTEGER :: i, ii, ix, iy, iz
    REAL(DP) :: psrho1, psrho2, acc, omegam
    REAL(DP) :: rx, ry, rz, r2, x1, y1, z1

    REAL(DP), DIMENSION(:), ALLOCATABLE :: rhops

    REAL(DP), PARAMETER:: prho1_data = -0.46073D0, prho2_data = 0.13287D0

!------------------------------------------------------------

    ALLOCATE (rhops(kdfull2))

    IF (nion2 == 2) RETURN

    IF (ipsptyp == 0) THEN

! total electronic charge from given density 'rho'

      apnum = 0D0
      DO i = 1, nxyz
        apnum = apnum + rho(i)
      END DO
      apnum = apnum*dvol

! positive background charge distribution for different cases:
! jellium or local Gaussian pseudo-densities

      IF (nion2 == 0) THEN
      DO i = 1, nxyz
        rhops(i) = rhojel(i)
      END DO
      ELSE
      DO i = 1, nxyz
        rhops(i) = 0D0
      END DO

      DO ii = 1, nion
        i = 0
        DO iz = minz, maxz
          z1 = (iz - nzsh)*dz
          rz = z1 - cz(ii)
          DO iy = miny, maxy
            y1 = (iy - nysh)*dy
            ry = y1 - cy(ii)
            DO ix = minx, maxx
              x1 = (ix - nxsh)*dx
              i = i + 1
              rx = x1 - cx(ii)
              r2 = rx*rx + ry*ry + rz*rz
              psrho1 = prho1_data*EXP(-r2/(2.0D0*sgm1(11)*sgm1(11)))
              psrho2 = prho2_data*EXP(-r2/(2.0D0*sgm2(11)*sgm2(11)))
              rhops(i) = rhops(i) + psrho1 + psrho2
            END DO
          END DO
        END DO
      END DO
      END IF

      rhomix = 0D0
      acc = 0D0
      DO i = 1, nxyz
        acc = acc + rhops(i)
        rhomix = rhomix + rho(i)*rhops(i)
      END DO
      rhopss = acc*dvol
      rhomix = rhomix*dvol

! 16*pi = 4*pi * e^2 * hbar^2/4m

      omegam = SQRT(16.0D0*pi/3.0D0*rhomix/apnum)
    ELSE
      omegam = 0D0 ! Mie plasmon not installed for Goedecker
      apnum = 0D0
      acc = 0D0
      rhomix = 0D0
    END IF

    DEALLOCATE (rhops)

    omega_mieplasmon = omegam

    RETURN
  END FUNCTION omega_mieplasmon
 