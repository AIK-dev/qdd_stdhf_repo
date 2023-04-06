
!!-----model-------------------------------------------------------------

! Code converted using TO_F90 by Alan Miller
! Date: 2010-04-16 Time: 12:37:44

SUBROUTINE model(a, nparmx, yfit, ndatmx)

  REAL, INTENT(IN OUT) :: a(kparmx)
  INTEGER, INTENT(IN) :: nparmx
  REAL, INTENT(IN OUT) :: yfit(kdatmx)
  INTEGER, INTENT(IN) :: ndatmx
  IMPLICIT REAL*8(a - h, o - z)

! compute chi**2 for Ar atom IN SIC with local Goedecker

  INTEGER, PARAMETER :: kparmx = 5
  INTEGER, PARAMETER :: kdatmx = 10

  CHARACTER(LEN=71) :: command

!-----------------------------------------------------------------------

  IF (nparmx /= 3) STOP 'locgoe.f requires NPARMX=2'
  IF (ndatmx /= 4) STOP 'locgoe.f requires NDATMX=4'
  OPEN (13, STATUS='scratch', FORM='FORMATTED')
! WRITE(6,'(a,4f6.1)')
! & 'model: ',(A(i),i=1,nparmx)
  WRITE (13, '(a,4f15.10)') 'fitpsp.bat ', 0.8, (a(i), i=1, nparmx)
! WRITE(13,'(a,4f15.10)')
! & 'fitpsp.bat ',0.6,A(1),0.0,A(2)
  REWIND 13
  READ (13, '(a)') command
  CLOSE (13)
! WRITE(6,'(a)') command
  CALL system(command)
  OPEN (13, FILE='energies.res')
  READ (13, '(1x)')
  READ (13, '(1x)')
  DO i = 1, ndatmx
    READ (13, *) id1, id2, yfit(i)
  END DO
  CLOSE (13)
  WRITE (6, '(a,4f7.2)') 'energies:', (yfit(i), i=1, ndatmx)
  WRITE (7, '(a,4f7.2)') 'energies:', (yfit(i), i=1, ndatmx)
! CALL system('tail -n 1 sumchi')

  RETURN
END SUBROUTINE model

