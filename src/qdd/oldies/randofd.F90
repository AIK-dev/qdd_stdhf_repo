
! **************************

SUBROUTINE ithion

! **************************

  USE params
!USE kinetic
  IMPLICIT REAL(DP) (A - H, O - Z)
! initial thermalization of the ions (ONLY executed IF needed)

! we ASSIGN 0.5*kb*t per degree of freedom

  ekin = 0.5*bk*tempion
  WRITE (7, '(a,f9.4)') 'wanted temp', tempion
  WRITE (7, '(a,f9.4)') 'wanted energ per degree of freedom', ekin
  xm = 0.5*1836.0*amu(np(1))*ame
  WRITE (7, *) 'warning : ithion not yet able TO treat unhomogeneous systems'

! attention TO masses !

  ekin = ekin/3/nion*(3*nion - 6.0)
  vmoy = SQRT(3.0*2.0*ekin/xm)
  WRITE (7, '(a,f9.4)') 'corresponding speed', vmoy
  pmoy = xm*vmoy

! we choose random directions of the speeds

  CALL spheric(pmoy, cpx, cpy, cpz, nion)
  CALL conslw(cx, cy, cz, cpx, cpy, cpz, nion)

! we rescale the speeds

  DO ion = 1, nion
    ek = cpx(ion)*cpx(ion)
    ek = ek + cpy(ion)*cpy(ion)
    ek = ek + cpz(ion)*cpz(ion)
    pn = SQRT(ek)
    cpx(ion) = cpx(ion)/pn*pmoy
    cpy(ion) = cpy(ion)/pn*pmoy
    cpz(ion) = cpz(ion)/pn*pmoy
  END DO
  CALL conslw(cx, cy, cz, cpx, cpy, cpz, nion)
  ekion = 0.0
  DO ion = 1, nion
    ek = 0.0
    ek = ek + cpx(ion)*cpx(ion)
    ek = ek + cpy(ion)*cpy(ion)
    ek = ek + cpz(ion)*cpz(ion)
    xm = 0.5*1836.0*amu(np(1))*ame
    ek = ek/2.0/xm
    WRITE (7, '(a,i2,a,f9.4)') 'ion', ion, 'ek=', ek
    ekion = ekion + ek
    WRITE (7, '(a,i2,3f9.4)') 'ion', ion, cpx(ion), cpy(ion), cpz(ion)
  END DO
  WRITE (7, '(a,f9.4)') 'kin.energ after renormalization', ekion
  ekion = ekion/(3*nion - 6.0)
  WRITE (7, '(a,f9.4)') 'av.kin.energ per net degree of freedom', ekion
  ekion = ekion*2.0/bk
  WRITE (7, '(a,f9.4)') 'corresponding temperature', ekion

  RETURN
END SUBROUTINE ithion
