!-----spectr-------------------------------------------------------------------

PROGRAM spectr

  PARAMETER(klinmx=199999, kcolmx=9)
  PARAMETER(kcoly=kcolmx - 1)
! PARAMETER (timefc=197.32) ! input DATA IN fm/c output MeV
  PARAMETER(timefc=0.0484) ! input DATA IN fs output Ry
! PARAMETER (timefc=1.0) ! input DATA IN 1/Ry (or 1/MeV)
  PARAMETER(omax=67.09) ! maximum plotted frequency IN eV
! PARAMETER (omax=3.0) ! maximum plotted frequency IN units
  PARAMETER(escl=13.6) ! conversion energy from Ry TO eV
! PARAMETER (escl=1.0) ! conversion energy from Ry TO ...
! e.g. from Ry TO eV --> 13.6
! set 1 TO leave output as is
  PARAMETER(ncosfilt=4) ! order of cosine filtering
! ! 0 ==> Gaussian filtering with 'twidth'
  PARAMETER(twidth=5.0/0.048) ! time width for filtering (input units)
! PARAMETER (twidth=12.0/0.048) ! time width for filtering (input units)

! Fourier transforms the columns of input DATA on FILE 5
! The input FILE is scanned for the first occurence of a
! line beginning with 'H:'. This Header-line CONTAINS
! information about the interpretation of the subsequent
! columns: 'X' marks the equidistant points of the axis,
! 'Yi' marks the columns TO be transformed, and 'N' marks
! columns TO be overridden.
! The columns are THEN READ line by line until an empty
! line or EoF is encountered.
! The spacing IN Fourier-space is determined from the largest
! 'X'-value provided.

! Example: the line
! "H: X N Y"
! WHERE "H:" has TO be IN the FIRST TWO columns
! determines that from the immediately subsequent DATA
! the first column is READ as times (TO be equidistant!),
! the second column is not READ,
! and the third column is READ as DATA TO be analysed.
! Note: ONLY one X column is possible, but up TO 12 Y columns
! can be handled simultanously.

  CHARACTER(LEN=50) :: comment
  CHARACTER(LEN=256) :: a
  CHARACTER(LEN=1) :: ytype(kcoly)
  DIMENSION xin(klinmx), yin(klinmx, kcoly), dlinin(kcolmx)
  DIMENSION ytrans(0:klinmx, kcoly), ynorm(kcoly)
  DIMENSION nycol(kcoly)
  DATA ynorm/kcoly*0.0/
! COMPLEX cfac,cprod
  COMPLEX :: cysum(kcoly), c
  COMPLEX :: cysav(0:klinmx, kcoly)
  COMPLEX :: cftfil, crdfil
  EXTERNAL cftfil, crdfil

  PARAMETER(pi=3.1415926)

  DATA itimfi/1/ ! switch TO exp. filtering IN time domain
  DATA iftimpo/0/ ! switch TO PRINT filtered power spectrum
  DATA istd/0/ ! switch TO PRINT standard unfiltered
  DATA iphase/0/ ! switch TO printing of phase
  DATA ifcutm/0/ ! cut at maximum of DATA
  DATA ifaver/0/ ! subtract average signal
  DATA ifrel0/1/ ! scale TO zero at start
  DATA ifauto/0/ ! switch TO compute autocorrelation
  DATA ifnorm/0/ ! switch TO PRINT power normalized
  DATA iftims/0/ ! PRINT protocol of filtered time signal

! internal statement FUNCTION for COMPLEX abs. value squared

!cabs2(c) = REAL(c*CONJG(c))

!-------------------------------------------------------------------------------

! OPEN scratch FILE

  OPEN (UNIT=19, STATUS='scratch')

  ncount = 0

  READ (5, '(a)', ERR=991, END=999) comment
  DO iread = 1, 9999
    READ (5, '(a)', ERR=991, END=999) a
! WRITE(6,'(2a)') 'N:',a
    IF (a(1:2) == 'H:') THEN
! analysing column structure
      maxcol = 0
      nxcol = 0
      nycolc = 0

      i = 3
10    CONTINUE
      IF (a(i:i) /= ' ') THEN
        maxcol = 1 + maxcol
        IF (a(i:i) == 'X' .OR. a(i:i) == 'x') THEN
          nxcol = maxcol
        ELSE IF (a(i:i) == 'Y') THEN
          nycolc = 1 + nycolc
          nycol(nycolc) = maxcol
          IF (a(i + 1:i + 1) /= ' ') THEN
            i = 1 + i
            ytype(nycolc) = a(i:i)
          END IF
        ELSE IF (a(i:i) /= 'N') THEN
          WRITE (6, '(a,i5)') ' wrong specifier IN Header line, column=', i
          STOP ' wrong specifier IN Header line.'
        END IF
      END IF
      i = 1 + i
      IF (i <= 80) GO TO 10

      DO i = 1, 80
        a(i:i) = ' '
      END DO
      DO i = 1, maxcol
        i1 = 7*i - 4
        a(i1:i1 + 4) = ' N '
      END DO
      DO i = 1, nycolc
        i1 = 7*nycol(i) - 5
        a(i1:i1 + 3) = ' Y1,'
        i1 = i1 + 4
        a(i1:i1) = ytype(i)
      END DO
      i1 = 7*nxcol - 2
      a(i1:i1) = 'X'
      a(1:2) = 'H:'
! WRITE(6,'(a)') a

      IF (nxcol == 0) STOP ' missing X-column'
      IF (nycolc == 0) STOP ' missing Y-column'
      GO TO 1009
    END IF
  END DO
  STOP ' no H: header line'
1009 CONTINUE

! reading DATA

  maxdat = 0
  DO idat = 1, klinmx
    READ (5, '(a)', ERR=993, END=29) a
! IF(idat.le.1001) WRITE(6,'(a)') a(1:LEN(a))
    IF (LEN(a) == 0) THEN
      GO TO 29
    ELSE
      maxdat = 1 + maxdat
      WRITE (19, '(a)') a
      REWIND 19
      READ (19, *) (dlinin(j), j=1, maxcol)
      REWIND 19
      xin(maxdat) = dlinin(nxcol)
      DO i = 1, nycolc
        yin(maxdat, i) = dlinin(nycol(i))
      END DO
      IF (maxdat > 1) THEN
        IF (ABS(xin(maxdat) - xin(maxdat - 1)) < 1.e-6) THEN
          maxdat = maxdat - 1
        END IF
      END IF
    END IF
  END DO
  READ (5, '(a)', ERR=29, END=29) a
  IF (LEN(a) /= 0) STOP ' too much DATA'
29 CONTINUE

! optionally cut DATA at last maximum of the first DATA set

  IF (ifcutm == 1) THEN
    ymax = yin(maxdat, 1)
    slope = ymax - yin(maxdat - 1, 1)
    DO ix = maxdat - 1, 1, -1
      maxnew = ix
      IF (ymax > yin(ix, 1) .AND. slope < 0.0) GO TO 353
      ymax = yin(ix, 1)
      slope = yin(ix + 1, 1) - ymax
    END DO
    STOP ' backward search for maximum failed'
353 CONTINUE
    maxdat = maxnew
  END IF

! compute transformation parameters

  delx = (xin(maxdat) - xin(1))/(maxdat - 1)
  xmin = xin(1) - delx ! ???
  delp = (pi + pi)/(maxdat*delx)

! supress average value of observables

  IF (ifaver == 1) THEN
    DO i = 1, nycolc
      ynorm(i) = 1.0E-30 ! set accumulator for max. value
      acc = 0.0
      DO ix = 1, maxdat
        acc = yin(ix, i) + acc
      END DO
      acc = acc/maxdat
      DO ix = 1, maxdat
        yin(ix, i) = yin(ix, i) - acc
      END DO
    END DO
  END IF

! scale signal TO start with zero

  IF (ifrel0 == 1) THEN
    DO i = 1, nycolc
      acc = yin(1, i)
      DO ix = 1, maxdat
        yin(ix, i) = yin(ix, i) - acc
      END DO
    END DO
  END IF

! determine maximum frequency bin

  DO ip = 1, maxdat/2 + 1
    pact = ip*delp - delp
    IF (pact*escl*timefc > omax) GO TO 119
    maxplt = ip + 1
  END DO
119 CONTINUE
  maxplt = MIN(maxdat/2 + 1, maxplt)
! WRITE(6,'(/a,i5)') ' now maxplt=',maxplt

! Fourier transformation and plot

  IF (istd /= 0) THEN
    DO i = 1, nycolc
      CALL slowft(yin(1, i), cysav(1, i), maxdat, maxplt, delx)
    END DO
    DO ip = 1, maxplt
      DO i = 1, nycolc
        ytrans(ip, i) = (cysav(ip, i))*CONJG(cysav(ip, i))
!? ytrans(maxdat+1-ip,i) = ytrans(ip,i)
        ynorm(i) = MAX(ytrans(ip, i), ynorm(i))
      END DO
    END DO
    IF (ifnorm /= 1) THEN
      DO i = 1, nycolc
        ynorm(i) = 1.0
      END DO
    END IF
! OPEN(7,FILE='strength.res')
    WRITE (6, '(a)') '#'//comment, &
      '# nr omega REAL(strength) imag(strength) power'
! OPEN(8,FILE='power.res')
! WRITE(8,'(a)')
! & '#'//comment,
! & '# nr omega power'
    DO ip = 1, maxplt
      IF (ip > 1) THEN
        pact = ip*delp - delp
        WRITE (6, '(i5,18(g12.4))') ip, pact*escl*timefc, &
          ((cysav(ip, i)/SQRT(ynorm(i))), (ytrans(ip, i)/ynorm(i)), i=1, nycolc)
! WRITE(8,'(i5,10(g12.4))')
! & ip,pact*escl*timefc,
! & (ytrans(ip,i)/ynorm(i),i=1,nycolc)
      END IF
    END DO
! CLOSE(7)
! CLOSE(8)
    WRITE (6, '(a)') ' ', ' '
  END IF
!? maxplt = min(maxplt,maxdat/2-2)

! optionally PRINT the phase

  IF (iphase == 1) THEN
! OPEN(7,FILE='phase.res')
    WRITE (6, '(a)') '#'//comment, &
      '# nr omega phase of strength'
    pih = pi/2.0
    DO ip = 3, maxplt
      WRITE (6, '(i5,10(g12.4))') ip, (ip*delp - delp)*escl*timefc, &
        (AIMAG(LOG(cysav(ip, i))), i=1, nycolc), pih, -pih
    END DO
! CLOSE(7)
    WRITE (6, '(a)') ' ', ' '
  END IF

! correlation FUNCTION directly

  IF (ifauto == 1) THEN
! OPEN(7,FILE='autocor.res')
    WRITE (6, '(a)') '#'//comment, &
      '# nr time autocorrelation-FUNCTION'
    DO ix = 1, maxdat/2 + 1
      DO i = 1, nycolc
        cysum(i) = CMPLX(0.0, 0.0)
      END DO
      DO ix2 = 1, maxdat
        ishift = MOD(ix + ix2 - 2, maxdat) + 1
        DO i = 1, nycolc
          cysum(i) = yin(ix2, i)*yin(ishift, i) + cysum(i)
        END DO
      END DO
      IF (ix == 1) THEN
        DO i = 1, nycolc
          ynorm(i) = 1.0/REAL(cysum(i))
        END DO
      END IF
      WRITE (6, '(i5,8(g12.4))') &
        ix, (ix - 1)*delx, (REAL(cysum(i))*ynorm(i), i=1, nycolc)
    END DO
! CLOSE(7)
    WRITE (6, '(a)') ' ', ' '
  END IF

  IF (itimfi /= 1) STOP

! filtering IN time domain

  xmax = (maxdat)*delx
  DO ix = 1, maxdat
    xact = (ix - 1)*delx
    IF (ncosfilt > 0) THEN
      winfac = (0.5 + 0.5*COS(pi*xact/xmax))**ncosfilt
    ELSE
      winfac = EXP(-(xact/twidth)**2)
    END IF
    DO i = 1, nycolc
      yin(ix, i) = yin(ix, i)*winfac
    END DO
  END DO

! Fourier transformation of filtered DATA

  DO i = 1, nycolc
    CALL slowft(yin(1, i), cysav(1, i), maxdat, maxplt, delx)
  END DO
  DO ip = 1, maxplt
    DO i = 1, nycolc
      ytrans(ip, i) = (cysav(ip, i))*CONJG(cysav(ip, i))
!? ytrans(maxdat+1-ip,i) = ytrans(ip,i)
      ynorm(i) = MAX(ytrans(ip, i), ynorm(i))
    END DO
  END DO
  IF (ifnorm /= 1) THEN
    DO i = 1, nycolc
      ynorm(i) = 1.0
    END DO
  END IF
! OPEN(7,FILE='fstrength.res')
  WRITE (6, '(a)') '#'//comment, &
    '# with filtering IN time domain', &
    '# nr omega REAL(strength) imag(strength) power'
! OPEN(8,FILE='fpower.res')
! WRITE(8,'(a)')
! & '#'//comment,
! & '# nr omega power'
  DO ip = 1, maxplt
    IF (ip > 1) THEN
      pact = ip*delp - delp
      WRITE (6, '(i5,18(g12.4))') ip, pact*escl*timefc, &
        ((cysav(ip, i)/SQRT(ynorm(i))), (ytrans(ip, i)/ynorm(i)), i=1, nycolc)
! WRITE(8,'(i5,10(g12.4))')
! & ip,pact*escl*timefc,
! & (ytrans(ip,i)/ynorm(i),i=1,nycolc)
    END IF
  END DO

! optionally PRINT time filtered time signal

  IF (iftims == 1) THEN
    WRITE (6, '(a)') ' ', ' ', '# filtered time signal', '# time signals'
    DO ix = 1, maxdat
      xact = (ix - 1)*delx
      WRITE (6, '(f12.3,18(g12.4))') xact, (yin(ix, i), i=1, nycolc)
    END DO
  END IF
! CLOSE(7)
! CLOSE(8)

  STOP ' normal END'

991 STOP ' error IN input before H:'
993 STOP ' error IN input columns'
999 STOP ' EoF IN input before H:'
END PROGRAM spectr
!-----sftfil-----------------------------------------------------------

FUNCTION sftfil(field, ip)

! soft filtering IN frequency space.
! field = array with DATA TO be filtered
! ip = POSITION at which the filtering is TO be extracted

  INTEGER, PARAMETER :: klinmx = 29999
  REAL, INTENT(IN) :: field(klinmx)
  INTEGER, INTENT(IN) :: ip

!----------------------------------------------------------------------

  IF (ip < 2) STOP ' SFTFIL called with IP<2'

  IF (ip >= 5) THEN
    sftfil = 0.4697*field(ip) + 0.2349*(field(ip - 1) + field(ip + 1)) &
             + 0.0294*(field(ip - 2) + field(ip + 2)) + 0.00092*(field(ip - 3) + field(ip + 3))
  ELSE IF (ip == 4) THEN
    sftfil = 0.4697*field(ip) + 0.2349*(field(ip - 1) + field(ip + 1)) &
             + 0.0294*(field(ip - 2) + field(ip + 2))
  ELSE
    sftfil = 0.4697*field(ip) + 0.2349*(field(ip - 1) + field(ip + 1))
  END IF

  RETURN
END FUNCTION sftfil
!-----hrdfil-----------------------------------------------------------

FUNCTION hrdfil(field, ip)

! simple nearest neighbour filtering IN frequency space.
! field = array with DATA TO be filtered
! ip = POSITION at which the filtering is TO be extracted

  INTEGER, PARAMETER :: klinmx = 29999

  REAL, INTENT(IN) :: field(klinmx)
  INTEGER, INTENT(IN) :: ip

!----------------------------------------------------------------------

  IF (ip < 2) STOP ' HRDFIL called with IP<2'

  hrdfil = 0.5*field(ip) + 0.25*(field(ip - 1) + field(ip + 1))

  RETURN
END FUNCTION hrdfil
!-----cftfil-----------------------------------------------------------

COMPLEX FUNCTION cftfil(field, ip)

! soft filtering IN frequency space.
! field = array with DATA TO be filtered
! ip = POSITION at which the filtering is TO be extracted
  INTEGER, PARAMETER :: klinmx = 29999

  COMPLEX, INTENT(IN) :: field(klinmx)
  INTEGER, INTENT(IN) :: ip

!----------------------------------------------------------------------

  IF (ip < 2) STOP ' CFTFIL called with IP<2'

  IF (ip >= 5) THEN
    cftfil = 0.4697*field(ip) + 0.2349*(field(ip - 1) + field(ip + 1)) &
             + 0.0294*(field(ip - 2) + field(ip + 2)) + 0.00092*(field(ip - 3) + field(ip + 3))
  ELSE IF (ip == 4) THEN
    cftfil = 0.4697*field(ip) + 0.2349*(field(ip - 1) + field(ip + 1)) &
             + 0.0294*(field(ip - 2) + field(ip + 2))
  ELSE
    cftfil = 0.4697*field(ip) + 0.2349*(field(ip - 1) + field(ip + 1))
  END IF

  RETURN
END FUNCTION cftfil
!-----hrdfil-----------------------------------------------------------

COMPLEX FUNCTION crdfil(field, ip)

! simple nearest neighbour filtering IN frequency space.
! field = array with DATA TO be filtered
! ip = POSITION at which the filtering is TO be extracted

  INTEGER, PARAMETER :: klinmx = 29999

  COMPLEX, INTENT(IN) :: field(klinmx)
  INTEGER, INTENT(IN) :: ip

!----------------------------------------------------------------------

  IF (ip < 2) STOP ' CRDFIL called with IP<2'

  crdfil = 0.5*field(ip) + 0.25*(field(ip - 1) + field(ip + 1))

  RETURN
END FUNCTION crdfil
!-----slowft-----------------------------------------------------------

SUBROUTINE slowft(ytime, cyomeg, maxdat, maxplt, delt)

  INTEGER, PARAMETER :: klinmx = 29999

  REAL, INTENT(IN OUT) :: ytime(klinmx)
  COMPLEX, INTENT(OUT) :: cyomeg(klinmx)
  INTEGER, INTENT(IN) :: maxdat
  INTEGER, INTENT(IN) :: maxplt
  REAL, INTENT(IN) :: delt
  REAL, PARAMETER :: omax = 9.89

! (slow) Fourier transformation
! ytime = input of REAL array of signal IN time domain;
! it is assumed that the first point corresponds TO time=0.
! cyomeg = output of COMPLEX array of transform infrequency domain
! ONLY positive frequencies are computed and returnd
! (maxdat/2+1 points)
! maxdat = NUMBER of DATA points IN time domain
! maxplt = NUMBER of DATA points TO be computed (plus a marigin
! of 4 extra points)
! delt = step SIZE IN time domain

  COMPLEX :: cysum ! accumulator
  COMPLEX :: cprod, cfac ! auxiliary

  DATA pi/3.1415926/

!----------------------------------------------------------------------

  delp = (pi + pi)/(maxdat*delt)
  DO ip = 1, MIN(maxdat/2 + 1, maxplt + 4)
    pact = ip*delp - delp
    cfac = CEXP(CMPLX(0.0, pact*delt))
!old cprod = cexp(cmplx(0.0,pact*xin(1)))/cfac
    cprod = CMPLX(1.0, 0.0)/cfac
    cysum = CMPLX(0.0, 0.0)
    DO ix = 1, maxdat
      cprod = cprod*cfac
      cysum = cprod*CMPLX(ytime(ix), 0.0) + cysum
    END DO
    cyomeg(ip) = delt*cysum
    cyomeg(maxdat + 1 - ip) = cyomeg(ip)
  END DO

  RETURN
END SUBROUTINE slowft

