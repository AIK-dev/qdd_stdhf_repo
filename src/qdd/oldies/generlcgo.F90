
SUBROUTINE genermowf(psiom, nmxst)

! generates orbital molecular wave FUNCTION: psiom

  USE params
  IMPLICIT NONE
#if(mpi)
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)
#endif

  REAL(DP), INTENT(OUT) :: psiom(kdfull2, kstate)
  INTEGER, INTENT(OUT) :: nmxst(1:ng) ! maximum nr. for each atom

  INTEGER :: i, ind, ion, ix, iy, iz
  INTEGER :: nadd, nadded, natlevel, nbr, ncy, ncycle, ndiff, nmaxval, nstx, nsty, nstz, numstate, numlocstate
  REAL(DP) :: x1, y1, z1
  INTEGER, ALLOCATABLE :: nactst(:) ! keep track of actual atomic state
  INTEGER, ALLOCATABLE :: ncount(:) ! keep track of actual atomic state
  INTEGER :: nnodes(1:3, 1:10)
  DATA nnodes/0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, &
    1, 1, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 2/

!---------------------------------------------------------------------

! reset wavefunctions 'psiom'

  ALLOCATE (nactst(1:ng))
  ALLOCATE (ncount(0:knode))
  nactst = 0
  ncount = 0

  WRITE (6, *) ' GENERLCGO entered. NION,NSTATE=', nion, nstate
  WRITE (6, *) ' CH(ion):', (ch(np(ion)), ion=1, nion)
  CALL flush (6)

  DO nbr = 1, nstate
    DO i = 1, kdfull2
      psiom(i, nbr) = 0D0
    END DO
  END DO

! book-keeping of NUMBER of states

  nmaxval = 0
  DO ion = 1, nion
    nmxst(ion) = ch(np(ion))
    nmaxval = nmaxval + ch(np(ion))
  END DO
  ndiff = nmaxval - nelect
  ncycle = ndiff/nion + 1
  IF (ndiff > 0) THEN
    nadd = +1
  ELSE
    nadd = -1
  END IF
!nadd=0 ! ?? check the initialization
  IF (ndiff /= 0) THEN
    nadded = 0
    cycleloop: DO ncy = 1, ncycle
      DO ion = nion, 1, -1
        nmxst(ion) = nmxst(ion) - nadd
        nadded = nadded + nadd
        IF (nadded == ndiff) EXIT cycleloop
      END DO
    END DO cycleloop
  END IF
  WRITE (6, '(a,5i5)') &
    ' ndiff,nmaxval,nstate,ncycle,nadd=', ndiff, nmaxval, nstate, ncycle, nadd
  WRITE (6, '(a,100i3)') ' nmxst:', (nmxst(ion), ion=1, nion)
  WRITE (6, '(a,100i3)') ' ipol:', (ipol(ion), ion=1, nion)
  WRITE (6, *) ' ipolcheck:', (ipol(ion), ion=1, nion)

! loop through ions and fill electron states successively

!nmaxact = nstate/nion+1 ! max. states per atom

  numstate = 0
  numlocstate = 0 ! IN CASE of parallelism
  ionloop: DO ion = 1, nion
    DO natlevel = 1, nmxst(ion)
      IF (numstate == ksttot) THEN
        WRITE (6, *) 'numstate = ksttot, increase kstate or nproc'
        WRITE (7, *) 'numstate = ksttot, increase kstate or nproc'
        EXIT ionloop
      END IF
      numstate = 1 + numstate
      ncount(nhome(numstate)) = ncount(nhome(numstate)) + 1
      IF (numspin == 2) THEN
        ispin(numstate) = 2 - MOD(numstate, 2)
        IF (ipol(ion) .eq. -1) ispin(numstate) = 2 - mod(numstate + 1, 2)
        nactst(ion) = nactst(ion) + MOD(natlevel, 2)
      ELSE
        nactst(ion) = nactst(ion) + 1
      END IF
#if(mpi)
      ispin_node(ncount(nhome(numstate)), nhome(numstate)) = ispin(numstate)
#endif
      IF (nhome(numstate) == myn) THEN
        numlocstate = numlocstate + 1

        WRITE (6, *) 'ion,natlev,numst,ispin', &
          ion, natlevel, numstate, ispin(numstate)
        WRITE (7, *) 'ksttot,nstate,ion,natlev,numst,ispin', &
          ksttot, nstate, ion, natlevel, numstate, ispin(numstate)
        IF (nactst(ion) > 10) STOP 'GENERMOWF: too high atomic level'
        IF (numlocstate > nstate) THEN
          WRITE (6, *) 'numlocstate > nstate, increase kstate or nproc'
          WRITE (7, *) 'numlocstate > nstate, increase kstate or nproc'
          EXIT ionloop
        END IF
! SELECT nodes
        nstx = nnodes(initord(1, ion), nactst(ion))
        nsty = nnodes(initord(2, ion), nactst(ion))
        nstz = nnodes(initord(3, ion), nactst(ion))
! compute raw wf
        ind = 0
        DO iz = 1, nz2
          z1 = ((iz - nzsh)*dz - cz(ion))/radini(ion)
          DO iy = 1, ny2
            y1 = ((iy - nysh)*dy - cy(ion))/radini(ion)
            DO ix = 1, nx2
              x1 = ((ix - nxsh)*dx - cx(ion))/radini(ion)
              ind = ind + 1
              psiom(ind, numlocstate) = (x1**nstx - 0.5D0*MAX(0, nstx - 1)) &
                                        *(y1**nsty - 0.5D0*MAX(0, nsty - 1))*(z1**nstz - 0.5D0*MAX(0, nstz - 1)) &
                                        *EXP(-(x1*x1 + y1*y1 + z1*z1)*0.5D0)
            END DO
          END DO
        END DO

        WRITE (6, '(a,i3,a,a,5i3,1pg13.5)') &
          ' electron state nr.', numstate, ': ', &
          ' ion,nstx,nsty,nstz=', ion, nstx, nsty, nstz, ispin(numstate), &
          SUM(psiom(:, numlocstate)**2)*dvol
      END IF
    END DO
  END DO ionloop

  DEALLOCATE (nactst)

!nstate=numlocstate

! Schmidt ortho-normalisation

!CALL schmidt(psiom)

!DO nbe=1,nstate
!
! DO IN=1,nbe
! IF((ispin(nbe) == ispin(IN))) THEN
! sum=0D0
! DO i=1,nxyz
! sum=sum+psiom(i,nbe)*psiom(i,IN)
! END DO
! sum=sum*dvol
! cs = SUM(psiom(:,nbe)*psiom(:,IN))*dvol
! IF(IN == nbe) THEN
! cs = 1D0/SQRT(cs)
! DO i=1,nxyz
! psiom(:,nbe)=psiom(:,nbe)*cs
! END DO
! ELSE
! DO i=1,nxyz
! psiom(:,nbe)=psiom(:,nbe)-cs*psiom(:,IN)
! END DO
! END IF
! END IF
! END DO
!END DO

!DO n=1,numstate
! WRITE(*,*) ' INITLCGO after Schmidt: state nr. norm=', &
! n,SUM(psiom(:,n)**2)*dvol
!END DO
  CALL flush (6)

  RETURN
END SUBROUTINE genermowf

