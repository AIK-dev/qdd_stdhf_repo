
!------------------------------------------------------------

SUBROUTINE evaluate(rho, aloc, psi)

! Not clear what this routine is doing. ????

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)
  REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)
  COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)
#if(raregas)
  INTEGER :: i
#endif
  REAL(DP) :: dummy, enpol, enpol0
  REAL(DP), EXTERNAL :: energ_ions

  WRITE (6, *) 'Doing some postrun evaluation ONLY...'

  CALL calcrho(rho, psi)
  CALL calclocal(rho, aloc)
  CALL calcpseudo()
  CALL calcrho(rho, psi)
  IF (ifsicp > 0) CALL calc_sic(rho, aloc, psi)
  jenergy = 1

  CALL info(psi, rho, aloc, 10)

  nion = 0
! CALL energ_ions() ??? BF
  dummy = energ_ions()

  enpol0 = 0D0
#if(raregas)
  OPEN (834, STATUS='unknown', FILE='penerinfty')
  enerinfty = engg
  WRITE (834, *) enerinfty
  CLOSE (834)

  DO i = 1, nc
    enpol0 = enpol0 + (xc(i) - xe(i))**2 + (yc(i) - ye(i))**2 + (zc(i) - ze(i))**2
  END DO
  enpol0 = 0.5D0*cspr*enpol0
#endif

  CALL restart2(psi, outnam, .false.)
  WRITE (6, *) 'Done.'

  CLOSE (163)

  CALL calcrho(rho, psi)
  CALL calclocal(rho, aloc)
  CALL calcpseudo()
  CALL calcrho(rho, psi)
  IF (ifsicp > 0) CALL calc_sic(rho, aloc, psi)
  jenergy = 1
  CALL info(psi, rho, aloc, 10)

  enpol = 0D0
#if(raregas)
  DO i = 1, nc
    enpol = enpol + (xc(i) - xe(i))**2 + (yc(i) - ye(i))**2 + (zc(i) - ze(i))**2
  END DO
  enpol = 0.5D0*cspr*enpol
  WRITE (6, '(a,3f17.8)') 'enii,enig,engg', enii, enig, engg - enerinfty
  WRITE (6, '(a,2f17.8)') 'enpol0,enpol', enpol0, enpol
#endif

  WRITE (6, *) 'Postrun evaluation done.'

  STOP 'Postrun evaluation done.'

  RETURN
END SUBROUTINE evaluate
!------------------------------------------------------------

SUBROUTINE evalprops
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! evaluates properties at a given iteration, such like
! kinetic energy of ions, cores, clouds and cations
! as well as coupling energies etc.

  INTEGER :: i
  REAL(DP) :: kinenergy
#if(raregas)
  REAL(DP) :: kinenergyc, kinenergye, kinenergyk
#endif
  REAL(DP) :: vsum, vvsum, vsumc, vvsumc, vsume, vvsume, vsumk, vvsumk
  REAL(DP) :: vx, vy, vz, xm

  vsum = 0D0
  vvsum = 0D0
  vsumc = 0D0
  vvsumc = 0D0
  vsume = 0D0
  vvsume = 0D0
  vsumk = 0D0
  vvsumk = 0D0

! kinetic energies first

  xm = amu(np(i)) ! is that correct??? check with ionmd.F

  DO i = 1, nion

! reset velocities by half a time step because of
! the USE of leapfrog
    vx = (cpx(i) - 0.5D0*fx(i)*dt1)
    vy = (cpy(i) - 0.5D0*fy(i)*dt1)
    vz = (cpz(i) - 0.5D0*fz(i)*dt1)
    vsum = vsum + vx + vy + vz
    vvsum = vvsum + vx*vx + vy*vy + vz*vz
  END DO

  kinenergy = 0.5D0*vvsum/xm

#if(raregas)
  xm = mion ! is that correct??? check with ionmd.F

  DO i = 1, nc

! reset velocities by half a time step because of
! the USE of leapfrog
    vx = (pxc(i) - 0.5D0*fxc(i)*dt1)
    vy = (pyc(i) - 0.5D0*fyc(i)*dt1)
    vz = (pzc(i) - 0.5D0*fzc(i)*dt1)
    vsumc = vsumc + vx + vy + vz
    vvsumc = vvsumc + vx*vx + vy*vy + vz*vz
  END DO

  kinenergyc = 0.5D0*vvsumc/xm/nc

  xm = me ! is that correct??? check with ionmd.F

  DO i = 1, NE

! reset velocities by half a time step because of
! the USE of leapfrog
    vx = (pxe(i) - 0.5D0*fxe(i)*dt1)
    vy = (pye(i) - 0.5D0*fye(i)*dt1)
    vz = (pze(i) - 0.5D0*fze(i)*dt1)
    vsume = vsume + vx + vy + vz
    vvsume = vvsume + vx*vx + vy*vy + vz*vz
  END DO

  kinenergye = 0.5D0*vvsume/xm/NE

  xm = mkat ! is that correct??? check with ionmd.F

  DO i = 1, nk

! reset velocities by half a time step because of
! the USE of leapfrog
    vx = (pxk(i) - 0.5D0*fxk(i)*dt1)
    vy = (pyk(i) - 0.5D0*fyk(i)*dt1)
    vz = (pzk(i) - 0.5D0*fzk(i)*dt1)
    vsumk = vsumk + vx + vy + vz
    vvsumk = vvsumk + vx*vx + vy*vy + vz*vz
  END DO

  kinenergyk = 0.5D0*vvsumk/xm/nk
#endif

  RETURN
END SUBROUTINE evalprops
!------------------------------------------------------------

