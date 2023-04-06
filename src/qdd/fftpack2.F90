! pieces from the FFT routines from NETLIB, converted to F90

!References
!
! P. N. Swarztrauber, Vectorizing the FFTs, IN Parallel
! Computations (G. Rodrigue, ed.), Academic Press,
! 1982, pp. 51-83.
!

SUBROUTINE cosqb(n, x, wsave, ifac)
! SUBROUTINE COSQB computes the fast Fourier transform of quarter
! wave DATA. That is, COSQB computes a SEQUENCE from its
! representation IN terms of a cosine series with odd wave numbers.
! The transform is defined below at output PARAMETER X.
!
! COSQB is the unnormalized inverse of COSQF since a CALL of COSQB
! followed by a CALL of COSQF will multiply the input SEQUENCE X
! by 4*N.
!
! The WSAVE and IFAC arrays which are used by SUBROUTINE COSQB must be
! initialized by calling SUBROUTINE COSQI(N,WSAVE).
!
! Input Parameters
!
! N the length of the array X TO be transformed. The method
! is most efficient when N is a product of small primes.
!
! X an array which CONTAINS the SEQUENCE TO be transformed
!
! WSAVE a work array which must be dimensioned at least 3*N+15
! IN the PROGRAM that calls COSQB. The WSAVE array must be
! initialized by calling SUBROUTINE COSQI(N,WSAVE), and a
! different WSAVE array must be used for each different
! value of N. This initialization does not have TO be
! repeated so long as N remains unchanged. Thus subsequent
! transforms can be obtained faster than the first.
!
! IFAC an INTEGER work array which must be dimensioned at least 15.
! Condtions applying on WSAVE also apply on IFAC.
!
! Output Parameters
!
! X For I=1,...,N
!
! X(I)= the sum from K=1 TO K=N of
!
! 2*X(K)*COS((2*K-1)*(I-1)*PI/(2*N))
!
! A CALL of COSQB followed by a CALL of
! COSQF will multiply the SEQUENCE X by 4*N.
! Therefore COSQF is the unnormalized inverse
! of COSQB.
!
! WSAVE and IFAC contain initialization calculations which must not
! be destroyed between calls of COSQB or COSQF.
  USE params, ONLY: DP, tsqrt2
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(OUT) :: x(n)
  REAL(DP), INTENT(IN OUT) :: wsave(*)
  INTEGER, INTENT(IN OUT) :: ifac(10)

  REAL(DP)::x1

  IF (n < 2) THEN
    x(1) = 4D0*x(1)
    RETURN
  ELSE IF (n == 2) THEN
    x1 = 4D0*(x(1) + x(2))
    x(2) = tsqrt2*(x(1) - x(2))
    x(1) = x1
    RETURN
  ELSE
    CALL cosqb1(n, x, wsave, wsave(n + 1), ifac)
    RETURN
  END IF

END SUBROUTINE cosqb

SUBROUTINE cosqb1(n, x, w, xh, ifac)
! SUBROUTINE COSQB1 computes the fast Fourier transform of quarter
! wave DATA. That is, COSQB1 computes a SEQUENCE from its
! representation IN terms of a cosine series with odd wave numbers.
! The transform is defined below at output PARAMETER X.
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN OUT) :: x(n)
  REAL(DP), INTENT(IN) :: w(*)
  REAL(DP), INTENT(OUT) :: xh(*)
  INTEGER, INTENT(OUT) :: ifac(10)

  INTEGER :: i, k, kc, modn, np2, ns2, xim1

  ns2 = (n + 1)/2
  np2 = n + 2
  DO i = 3, n, 2
    xim1 = x(i - 1) + x(i)
    x(i) = x(i) - x(i - 1)
    x(i - 1) = xim1
  END DO
  x(1) = x(1) + x(1)
  modn = MOD(n, 2)
  IF (modn == 0) x(n) = x(n) + x(n)
  CALL rfftb(n, x, xh, ifac)
  DO k = 2, ns2
    kc = np2 - k
    xh(k) = w(k - 1)*x(kc) + w(kc - 1)*x(k)
    xh(kc) = w(k - 1)*x(k) - w(kc - 1)*x(kc)
  END DO
  IF (modn == 0) x(ns2 + 1) = w(ns2)*(x(ns2 + 1) + x(ns2 + 1))
  DO k = 2, ns2
    kc = np2 - k
    x(k) = xh(k) + xh(kc)
    x(kc) = xh(k) - xh(kc)
  END DO
  x(1) = x(1) + x(1)
  RETURN
END SUBROUTINE cosqb1

SUBROUTINE radb2(ido, l1, cc, ch, wa1)
! Calculate the fast Fourier transform of subvectors of
! length two.
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, 2, l1)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, 2)
  REAL(DP), INTENT(IN) :: wa1(*)

  INTEGER :: i, ic, idp2, k
  REAL(DP) :: ti2, tr2

  DO k = 1, l1
    ch(1, k, 1) = cc(1, 1, k) + cc(ido, 2, k)
    ch(1, k, 2) = cc(1, 1, k) - cc(ido, 2, k)
  END DO

  IF (ido < 2) THEN
    RETURN
  ELSE IF (ido > 2) THEN
    idp2 = ido + 2
    DO k = 1, l1
      DO i = 3, ido, 2
        ic = idp2 - i
        ch(i - 1, k, 1) = cc(i - 1, 1, k) + cc(ic - 1, 2, k)
        tr2 = cc(i - 1, 1, k) - cc(ic - 1, 2, k)
        ch(i, k, 1) = cc(i, 1, k) - cc(ic, 2, k)
        ti2 = cc(i, 1, k) + cc(ic, 2, k)
        ch(i - 1, k, 2) = wa1(i - 2)*tr2 - wa1(i - 1)*ti2
        ch(i, k, 2) = wa1(i - 2)*ti2 + wa1(i - 1)*tr2
      END DO
    END DO
    IF (MOD(ido, 2) == 1) RETURN
  END IF

  DO k = 1, l1
    ch(ido, k, 1) = cc(ido, 1, k) + cc(ido, 1, k)
    ch(ido, k, 2) = -(cc(1, 2, k) + cc(1, 2, k))
  END DO

  RETURN
END SUBROUTINE radb2

SUBROUTINE radb3(ido, l1, cc, ch, wa1, wa2)
! Calculate the fast Fourier transform of subvectors of
! length three.
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, 3, l1)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, 3)
  REAL(DP), INTENT(IN) :: wa1(*)
  REAL(DP), INTENT(IN) :: wa2(*)

  INTEGER:: i, ic, idp2, k
  REAL(DP):: ci2, ci3, cr2, cr3
  REAL(DP):: di2, di3, dr2, dr3
  REAL(DP):: taui, taur, ti2, tr2

  DATA taur, taui/-.5D0, .866025403784439D0/

  DO k = 1, l1
    tr2 = cc(ido, 2, k) + cc(ido, 2, k)
    cr2 = cc(1, 1, k) + taur*tr2
    ch(1, k, 1) = cc(1, 1, k) + tr2
    ci3 = taui*(cc(1, 3, k) + cc(1, 3, k))
    ch(1, k, 2) = cr2 - ci3
    ch(1, k, 3) = cr2 + ci3
  END DO

  IF (ido == 1) RETURN
  idp2 = ido + 2
  DO k = 1, l1
    DO i = 3, ido, 2
      ic = idp2 - i
      tr2 = cc(i - 1, 3, k) + cc(ic - 1, 2, k)
      cr2 = cc(i - 1, 1, k) + taur*tr2
      ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2
      ti2 = cc(i, 3, k) - cc(ic, 2, k)
      ci2 = cc(i, 1, k) + taur*ti2
      ch(i, k, 1) = cc(i, 1, k) + ti2
      cr3 = taui*(cc(i - 1, 3, k) - cc(ic - 1, 2, k))
      ci3 = taui*(cc(i, 3, k) + cc(ic, 2, k))
      dr2 = cr2 - ci3
      dr3 = cr2 + ci3
      di2 = ci2 + cr3
      di3 = ci2 - cr3
      ch(i - 1, k, 2) = wa1(i - 2)*dr2 - wa1(i - 1)*di2
      ch(i, k, 2) = wa1(i - 2)*di2 + wa1(i - 1)*dr2
      ch(i - 1, k, 3) = wa2(i - 2)*dr3 - wa2(i - 1)*di3
      ch(i, k, 3) = wa2(i - 2)*di3 + wa2(i - 1)*dr3
    END DO
  END DO
  RETURN
END SUBROUTINE radb3

SUBROUTINE radb4(ido, l1, cc, ch, wa1, wa2, wa3)
! Calculate the fast Fourier transform of subvectors of
! length four.
  USE params, ONLY: DP, sq2
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, 4, l1)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, 4)
  REAL(DP), INTENT(IN) :: wa1(*)
  REAL(DP), INTENT(IN) :: wa2(*)
  REAL(DP), INTENT(IN) :: wa3(*)

  INTEGER:: i, ic, idp2, k
  REAL(DP):: ci2, ci3, ci4, cr2, cr3, cr4
  REAL(DP):: ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4

  DO k = 1, l1
    tr1 = cc(1, 1, k) - cc(ido, 4, k)
    tr2 = cc(1, 1, k) + cc(ido, 4, k)
    tr3 = cc(ido, 2, k) + cc(ido, 2, k)
    tr4 = cc(1, 3, k) + cc(1, 3, k)
    ch(1, k, 1) = tr2 + tr3
    ch(1, k, 2) = tr1 - tr4
    ch(1, k, 3) = tr2 - tr3
    ch(1, k, 4) = tr1 + tr4
  END DO
  IF (ido < 2) THEN
    RETURN
  ELSE IF (ido > 2) THEN
    idp2 = ido + 2
    DO k = 1, l1
      DO i = 3, ido, 2
        ic = idp2 - i
        ti1 = cc(i, 1, k) + cc(ic, 4, k)
        ti2 = cc(i, 1, k) - cc(ic, 4, k)
        ti3 = cc(i, 3, k) - cc(ic, 2, k)
        tr4 = cc(i, 3, k) + cc(ic, 2, k)
        tr1 = cc(i - 1, 1, k) - cc(ic - 1, 4, k)
        tr2 = cc(i - 1, 1, k) + cc(ic - 1, 4, k)
        ti4 = cc(i - 1, 3, k) - cc(ic - 1, 2, k)
        tr3 = cc(i - 1, 3, k) + cc(ic - 1, 2, k)
        ch(i - 1, k, 1) = tr2 + tr3
        cr3 = tr2 - tr3
        ch(i, k, 1) = ti2 + ti3
        ci3 = ti2 - ti3
        cr2 = tr1 - tr4
        cr4 = tr1 + tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4
        ch(i - 1, k, 2) = wa1(i - 2)*cr2 - wa1(i - 1)*ci2
        ch(i, k, 2) = wa1(i - 2)*ci2 + wa1(i - 1)*cr2
        ch(i - 1, k, 3) = wa2(i - 2)*cr3 - wa2(i - 1)*ci3
        ch(i, k, 3) = wa2(i - 2)*ci3 + wa2(i - 1)*cr3
        ch(i - 1, k, 4) = wa3(i - 2)*cr4 - wa3(i - 1)*ci4
        ch(i, k, 4) = wa3(i - 2)*ci4 + wa3(i - 1)*cr4
      END DO
    END DO
    IF (MOD(ido, 2) == 1) RETURN
  END IF

  DO k = 1, l1
    ti1 = cc(1, 2, k) + cc(1, 4, k)
    ti2 = cc(1, 4, k) - cc(1, 2, k)
    tr1 = cc(ido, 1, k) - cc(ido, 3, k)
    tr2 = cc(ido, 1, k) + cc(ido, 3, k)
    ch(ido, k, 1) = tr2 + tr2
    ch(ido, k, 2) = sq2*(tr1 - ti1)
    ch(ido, k, 3) = ti2 + ti2
    ch(ido, k, 4) = -sq2*(tr1 + ti1)
  END DO
  RETURN
END SUBROUTINE radb4

SUBROUTINE radb5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
! Calculate the fast Fourier transform of subvectors of
! length five.
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, 5, l1)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, 5)
  REAL(DP), INTENT(IN) :: wa1(*)
  REAL(DP), INTENT(IN) :: wa2(*)
  REAL(DP), INTENT(IN) :: wa3(*)
  REAL(DP), INTENT(IN) :: wa4(*)

  INTEGER:: i, ic, idp2, k
  REAL(DP):: ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5
  REAL(DP):: di2, di3, di4, di5, dr2, dr3, dr4, dr5
  REAL(DP):: ti2, ti3, ti4, ti5, tr2, tr3, tr4, tr5
  REAL(DP):: tr11, ti11, tr12, ti12

  DATA tr11, ti11, tr12, ti12/.309016994374947D0, .951056516295154D0, &
    -.809016994374947D0, .587785252292473D0/

  DO k = 1, l1
    ti5 = cc(1, 3, k) + cc(1, 3, k)
    ti4 = cc(1, 5, k) + cc(1, 5, k)
    tr2 = cc(ido, 2, k) + cc(ido, 2, k)
    tr3 = cc(ido, 4, k) + cc(ido, 4, k)
    ch(1, k, 1) = cc(1, 1, k) + tr2 + tr3
    cr2 = cc(1, 1, k) + tr11*tr2 + tr12*tr3
    cr3 = cc(1, 1, k) + tr12*tr2 + tr11*tr3
    ci5 = ti11*ti5 + ti12*ti4
    ci4 = ti12*ti5 - ti11*ti4
    ch(1, k, 2) = cr2 - ci5
    ch(1, k, 3) = cr3 - ci4
    ch(1, k, 4) = cr3 + ci4
    ch(1, k, 5) = cr2 + ci5
  END DO
  IF (ido == 1) RETURN
  idp2 = ido + 2
  DO k = 1, l1
    DO i = 3, ido, 2
      ic = idp2 - i
      ti5 = cc(i, 3, k) + cc(ic, 2, k)
      ti2 = cc(i, 3, k) - cc(ic, 2, k)
      ti4 = cc(i, 5, k) + cc(ic, 4, k)
      ti3 = cc(i, 5, k) - cc(ic, 4, k)
      tr5 = cc(i - 1, 3, k) - cc(ic - 1, 2, k)
      tr2 = cc(i - 1, 3, k) + cc(ic - 1, 2, k)
      tr4 = cc(i - 1, 5, k) - cc(ic - 1, 4, k)
      tr3 = cc(i - 1, 5, k) + cc(ic - 1, 4, k)
      ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2 + tr3
      ch(i, k, 1) = cc(i, 1, k) + ti2 + ti3
      cr2 = cc(i - 1, 1, k) + tr11*tr2 + tr12*tr3
      ci2 = cc(i, 1, k) + tr11*ti2 + tr12*ti3
      cr3 = cc(i - 1, 1, k) + tr12*tr2 + tr11*tr3
      ci3 = cc(i, 1, k) + tr12*ti2 + tr11*ti3
      cr5 = ti11*tr5 + ti12*tr4
      ci5 = ti11*ti5 + ti12*ti4
      cr4 = ti12*tr5 - ti11*tr4
      ci4 = ti12*ti5 - ti11*ti4
      dr3 = cr3 - ci4
      dr4 = cr3 + ci4
      di3 = ci3 + cr4
      di4 = ci3 - cr4
      dr5 = cr2 + ci5
      dr2 = cr2 - ci5
      di5 = ci2 - cr5
      di2 = ci2 + cr5
      ch(i - 1, k, 2) = wa1(i - 2)*dr2 - wa1(i - 1)*di2
      ch(i, k, 2) = wa1(i - 2)*di2 + wa1(i - 1)*dr2
      ch(i - 1, k, 3) = wa2(i - 2)*dr3 - wa2(i - 1)*di3
      ch(i, k, 3) = wa2(i - 2)*di3 + wa2(i - 1)*dr3
      ch(i - 1, k, 4) = wa3(i - 2)*dr4 - wa3(i - 1)*di4
      ch(i, k, 4) = wa3(i - 2)*di4 + wa3(i - 1)*dr4
      ch(i - 1, k, 5) = wa4(i - 2)*dr5 - wa4(i - 1)*di5
      ch(i, k, 5) = wa4(i - 2)*di5 + wa4(i - 1)*dr5
    END DO
  END DO
  RETURN
END SUBROUTINE radb5

SUBROUTINE radbg(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
! Calculate the fast Fourier transform of subvectors of
! arbitrary length.
  USE params, ONLY: DP, tPI
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: ip
  INTEGER, INTENT(IN) :: l1
  INTEGER, INTENT(IN) :: idl1
  REAL(DP), INTENT(IN) :: cc(ido, ip, l1)
  REAL(DP), INTENT(IN OUT) :: c1(ido, l1, ip)
  REAL(DP), INTENT(OUT) :: c2(idl1, ip)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, ip)
  REAL(DP), INTENT(IN OUT) :: ch2(idl1, ip)
  REAL(DP), INTENT(IN) :: wa(*)

  INTEGER:: i, ic, idij, idp2, ik, ipp2, ipph, is
  INTEGER:: j, j2, jc, k, l, lc, nbd
  REAL(DP):: ai1, ai2, ar1, ar1h, ar2, ar2h, arg, dc2, dcp, ds2, dsp

  arg = tpi/REAL(ip, DP)
  dcp = COS(arg)
  dsp = SIN(arg)
  idp2 = ido + 2
  nbd = (ido - 1)/2
  ipp2 = ip + 2
  ipph = (ip + 1)/2
  IF (ido < l1) THEN
    DO i = 1, ido
      DO k = 1, l1
        ch(i, k, 1) = cc(i, 1, k)
      END DO
    END DO
  ELSE
    DO k = 1, l1
      DO i = 1, ido
        ch(i, k, 1) = cc(i, 1, k)
      END DO
    END DO
  END IF

  DO j = 2, ipph
    jc = ipp2 - j
    j2 = j + j
    DO k = 1, l1
      ch(1, k, j) = cc(ido, j2 - 2, k) + cc(ido, j2 - 2, k)
      ch(1, k, jc) = cc(1, j2 - 1, k) + cc(1, j2 - 1, k)
    END DO
  END DO
  IF (ido /= 1) THEN
    IF (nbd < l1) THEN
      DO j = 2, ipph
        jc = ipp2 - j
        DO i = 3, ido, 2
          ic = idp2 - i
          DO k = 1, l1
            ch(i - 1, k, j) = cc(i - 1, 2*j - 1, k) + cc(ic - 1, 2*j - 2, k)
            ch(i - 1, k, jc) = cc(i - 1, 2*j - 1, k) - cc(ic - 1, 2*j - 2, k)
            ch(i, k, j) = cc(i, 2*j - 1, k) - cc(ic, 2*j - 2, k)
            ch(i, k, jc) = cc(i, 2*j - 1, k) + cc(ic, 2*j - 2, k)
          END DO
        END DO
      END DO
    ELSE
      DO j = 2, ipph
        jc = ipp2 - j
        DO k = 1, l1
          DO i = 3, ido, 2
            ic = idp2 - i
            ch(i - 1, k, j) = cc(i - 1, 2*j - 1, k) + cc(ic - 1, 2*j - 2, k)
            ch(i - 1, k, jc) = cc(i - 1, 2*j - 1, k) - cc(ic - 1, 2*j - 2, k)
            ch(i, k, j) = cc(i, 2*j - 1, k) - cc(ic, 2*j - 2, k)
            ch(i, k, jc) = cc(i, 2*j - 1, k) + cc(ic, 2*j - 2, k)
          END DO
        END DO
      END DO
    END IF
  END IF
  ar1 = 1D0
  ai1 = 0D0
  DO l = 2, ipph
    lc = ipp2 - l
    ar1h = dcp*ar1 - dsp*ai1
    ai1 = dcp*ai1 + dsp*ar1
    ar1 = ar1h
    DO ik = 1, idl1
      c2(ik, l) = ch2(ik, 1) + ar1*ch2(ik, 2)
      c2(ik, lc) = ai1*ch2(ik, ip)
    END DO
    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1
    DO j = 3, ipph
      jc = ipp2 - j
      ar2h = dc2*ar2 - ds2*ai2
      ai2 = dc2*ai2 + ds2*ar2
      ar2 = ar2h
      DO ik = 1, idl1
        c2(ik, l) = c2(ik, l) + ar2*ch2(ik, j)
        c2(ik, lc) = c2(ik, lc) + ai2*ch2(ik, jc)
      END DO
    END DO
  END DO
  DO j = 2, ipph
    DO ik = 1, idl1
      ch2(ik, 1) = ch2(ik, 1) + ch2(ik, j)
    END DO
  END DO
  DO j = 2, ipph
    jc = ipp2 - j
    DO k = 1, l1
      ch(1, k, j) = c1(1, k, j) - c1(1, k, jc)
      ch(1, k, jc) = c1(1, k, j) + c1(1, k, jc)
    END DO
  END DO
  IF (ido /= 1) THEN
    IF (nbd < l1) THEN
      DO j = 2, ipph
        jc = ipp2 - j
        DO i = 3, ido, 2
          DO k = 1, l1
            ch(i - 1, k, j) = c1(i - 1, k, j) - c1(i, k, jc)
            ch(i - 1, k, jc) = c1(i - 1, k, j) + c1(i, k, jc)
            ch(i, k, j) = c1(i, k, j) + c1(i - 1, k, jc)
            ch(i, k, jc) = c1(i, k, j) - c1(i - 1, k, jc)
          END DO
        END DO
      END DO
    ELSE
      DO j = 2, ipph
        jc = ipp2 - j
        DO k = 1, l1
          DO i = 3, ido, 2
            ch(i - 1, k, j) = c1(i - 1, k, j) - c1(i, k, jc)
            ch(i - 1, k, jc) = c1(i - 1, k, j) + c1(i, k, jc)
            ch(i, k, j) = c1(i, k, j) + c1(i - 1, k, jc)
            ch(i, k, jc) = c1(i, k, j) - c1(i - 1, k, jc)
          END DO
        END DO
      END DO
    END IF
  END IF

  IF (ido == 1) RETURN

  DO ik = 1, idl1
    c2(ik, 1) = ch2(ik, 1)
  END DO
  DO j = 2, ip
    DO k = 1, l1
      c1(1, k, j) = ch(1, k, j)
    END DO
  END DO
  IF (nbd > l1) THEN
    is = -ido
    DO j = 2, ip
      is = is + ido
      DO k = 1, l1
        idij = is
        DO i = 3, ido, 2
          idij = idij + 2
          c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) - wa(idij)*ch(i, k, j)
          c1(i, k, j) = wa(idij - 1)*ch(i, k, j) + wa(idij)*ch(i - 1, k, j)
        END DO
      END DO
    END DO
  ELSE
    is = -ido
    DO j = 2, ip
      is = is + ido
      idij = is
      DO i = 3, ido, 2
        idij = idij + 2
        DO k = 1, l1
          c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) - wa(idij)*ch(i, k, j)
          c1(i, k, j) = wa(idij - 1)*ch(i, k, j) + wa(idij)*ch(i - 1, k, j)
        END DO
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE radbg

SUBROUTINE rfftb(n, r, wsave, ifac)
!Calls for rfftb1
! Compute the backward fast Fourier transform of a REAL
! coefficient array
! SUBROUTINE RFFTB computes the REAL periodic SEQUENCE from its
! Fourier coefficients (Fourier synthesis). The transform is defined
! below at output PARAMETER R.
!
! Input Arguments
!
! N the length of the array R TO be transformed. The method
! is most efficient when N is a product of small primes.
! N may change so long as different work arrays are provided.
!
! R a REAL array of length N which CONTAINS the SEQUENCE
! TO be transformed.
!
! WSAVE a work array which must be dimensioned at least 2*N
! IN the PROGRAM that calls RFFTB. The WSAVE array must be
! initialized by calling SUBROUTINE RFFTI, and a different
! WSAVE array must be used for each different value of N.
! This initialization does not have TO be repeated so long as
! remains unchanged. Thus subsequent transforms can be
! obtained faster than the first. Moreover, the same WSAVE
! array can be used by RFFTF and RFFTB as long as N remains
! unchanged.
!
! IFAC an INTEGER work array which must be dimensioned at least 15.
! Condtions applying on WSAVE also apply on IFAC.
!
! Output Argument
!
! R For N even and for I = 1,...,N
!
! R(I) = R(1)+(-1)**(I-1)*R(N)
!
! plus the sum from K=2 TO K=N/2 of
!
! 2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
!
! -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
!
! For N odd and for I = 1,...,N
!
! R(I) = R(1) plus the sum from K=2 TO K=(N+1)/2 of
!
! 2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
!
! -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
!
! Note: This transform is unnormalized since a CALL of RFFTF
! followed by a CALL of RFFTB will multiply the input
! SEQUENCE by N.
!
! WSAVE and IFAC contain results which must not be destroyed between
! calls of RFFTB or RFFTF.
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN OUT) :: r(n)
  REAL(DP), INTENT(IN OUT) :: wsave(*)
  INTEGER, INTENT(IN OUT) :: ifac(10)

  IF (n == 1) RETURN
  CALL rfftb1(n, r, wsave, wsave(n + 1), ifac)
  RETURN
END SUBROUTINE rfftb

SUBROUTINE rfftb1(n, c, ch, wa, ifac)

  USE params, ONLY: DP
  IMPLICIT NONE
!See rfftb1 above
  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(OUT) :: c(n)
  REAL(DP), INTENT(IN OUT) :: ch(*)
  REAL(DP), INTENT(IN OUT) :: wa(*)
  INTEGER, INTENT(IN) :: ifac(10)

  INTEGER:: i, idl1, ido, ip, iw, ix2, ix3, ix4
  INTEGER:: k1, l1, l2, na, nf

  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1
  DO k1 = 1, nf
    ip = ifac(k1 + 2)
    l2 = ip*l1
    ido = n/l2
    idl1 = ido*l1
    !-------
    SELECT CASE (ip)
      !-------
    CASE (4)
      ix2 = iw + ido
      ix3 = ix2 + ido
      IF (na /= 0) THEN
        CALL radb4(ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3))
      ELSE
        CALL radb4(ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3))
      END IF
      na = 1 - na
      !-------
    CASE (2)
      IF (na /= 0) THEN
        CALL radb2(ido, l1, ch, c, wa(iw))
      ELSE
        CALL radb2(ido, l1, c, ch, wa(iw))
      END IF
      na = 1 - na
      !-------
    CASE (3)
      ix2 = iw + ido
      IF (na /= 0) THEN
        CALL radb3(ido, l1, ch, c, wa(iw), wa(ix2))
      ELSE
        CALL radb3(ido, l1, c, ch, wa(iw), wa(ix2))
      END IF
      na = 1 - na
      !-------
    CASE (5)
      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido
      IF (na /= 0) THEN
        CALL radb5(ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4))
      ELSE
        CALL radb5(ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4))
      END IF
      na = 1 - na
      !-------
    CASE DEFAULT
      IF (na /= 0) THEN
        CALL radbg(ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw))
      ELSE
        CALL radbg(ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw))
      END IF
      IF (ido == 1) na = 1 - na
      !-------
    END SELECT
    !-------
    l1 = l2
    iw = iw + (ip - 1)*ido
  END DO
  IF (na == 0) RETURN
  DO i = 1, n
    c(i) = ch(i)
  END DO
  RETURN
END SUBROUTINE rfftb1

SUBROUTINE cosqf(n, x, wsave, ifac)
! SUBROUTINE COSQF computes the fast Fourier transform of quarter
! wave DATA. That is, COSQF computes the coefficients IN a cosine
! series representation with ONLY odd wave numbers. The transform
! is defined below at Output PARAMETER X
!
! COSQF is the unnormalized inverse of COSQB since a CALL of COSQF
! followed by a CALL of COSQB will multiply the input SEQUENCE X
! by 4*N.
!
! The WSAVE and IFAC arrays which are used by SUBROUTINE COSQF must be
! initialized by calling SUBROUTINE COSQI(N,WSAVE,IFAC).
!
! Input Parameters
!
! N the length of the array X TO be transformed. The method
! is most efficient when N is a product of small primes.
!
! X an array which CONTAINS the SEQUENCE TO be transformed
!
! WSAVE a work array which must be dimensioned at least 3*N+15
! IN the PROGRAM that calls COSQF. The WSAVE array must be
! initialized by calling SUBROUTINE COSQI(N,WSAVE), and a
! different WSAVE array must be used for each different
! value of N. This initialization does not have TO be
! repeated so long as N remains unchanged. Thus subsequent
! transforms can be obtained faster than the first.
!
! IFAC an INTEGER work array which must be dimensioned at least 15.
! Condtions applying on WSAVE also apply on IFAC.
!
! Output Parameters
!
! X For I=1,...,N
!
! X(I) = X(1) plus the sum from K=2 TO K=N of
!
! 2*X(K)*COS((2*I-1)*(K-1)*PI/(2*N))
!
! A CALL of COSQF followed by a CALL of
! COSQB will multiply the SEQUENCE X by 4*N.
! Therefore COSQB is the unnormalized inverse
! of COSQF.
!
! WSAVE and IFAC contain initialization calculations which must not
! be destroyed between calls of COSQF or COSQB.
  USE params, ONLY: DP, sq2
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN OUT) :: x(n)
  REAL(DP), INTENT(IN OUT) :: wsave(*)
  INTEGER, INTENT(IN OUT) :: ifac(10)

  REAL(DP):: tsqx

  IF (n == 2) THEN
    tsqx = sq2*x(2)
    x(2) = x(1) - tsqx
    x(1) = x(1) + tsqx
  ELSE IF (n > 2) THEN
    CALL cosqf1(n, x, wsave, wsave(n + 1), ifac)
  END IF

  RETURN
END SUBROUTINE cosqf

SUBROUTINE cosqf1(n, x, w, xh, ifac)
! SUBROUTINE COSQF1 computes the fast Fourier transform of quarter
! wave DATA. That is, COSQF1 computes the coefficients IN a cosine
! series representation with ONLY odd wave numbers. The transform
! is defined below at Output PARAMETER X
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN OUT) :: x(n)
  REAL(DP), INTENT(IN) :: w(*)
  REAL(DP), INTENT(OUT) :: xh(*)
  INTEGER, INTENT(IN) :: ifac(10)

  INTEGER:: i, k, kc, modn, np2, ns2
  REAL(DP):: xim1

  ns2 = (n + 1)/2
  np2 = n + 2
  DO k = 2, ns2
    kc = np2 - k
    xh(k) = x(k) + x(kc)
    xh(kc) = x(k) - x(kc)
  END DO
  modn = MOD(n, 2)
  IF (modn == 0) xh(ns2 + 1) = x(ns2 + 1) + x(ns2 + 1)
  DO k = 2, ns2
    kc = np2 - k
    x(k) = w(k - 1)*xh(kc) + w(kc - 1)*xh(k)
    x(kc) = w(k - 1)*xh(k) - w(kc - 1)*xh(kc)
  END DO
  IF (modn == 0) x(ns2 + 1) = w(ns2)*xh(ns2 + 1)
  CALL rfftf(n, x, xh, ifac)
  DO i = 3, n, 2
    xim1 = x(i - 1) - x(i)
    x(i) = x(i - 1) + x(i)
    x(i - 1) = xim1
  END DO
  RETURN
END SUBROUTINE cosqf1

SUBROUTINE radf2(ido, l1, cc, ch, wa1)
! Calculate the fast Fourier transform of subvectors of
! length two.
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, l1, 2)
  REAL(DP), INTENT(OUT) :: ch(ido, 2, l1)
  REAL(DP), INTENT(IN) :: wa1(*)

  INTEGER:: i, ic, idp2, k
  REAL(DP):: ti2, tr2

  DO k = 1, l1
    ch(1, 1, k) = cc(1, k, 1) + cc(1, k, 2)
    ch(ido, 2, k) = cc(1, k, 1) - cc(1, k, 2)
  END DO
  IF (ido < 2) THEN
    RETURN
  ELSE IF (ido >= 2) THEN
    idp2 = ido + 2
    DO k = 1, l1
      DO i = 3, ido, 2
        ic = idp2 - i
        tr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
        ti2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
        ch(i, 1, k) = cc(i, k, 1) + ti2
        ch(ic, 2, k) = ti2 - cc(i, k, 1)
        ch(i - 1, 1, k) = cc(i - 1, k, 1) + tr2
        ch(ic - 1, 2, k) = cc(i - 1, k, 1) - tr2
      END DO
    END DO
    IF (MOD(ido, 2) == 1) RETURN
  END IF

  DO k = 1, l1
    ch(1, 2, k) = -cc(ido, k, 2)
    ch(ido, 1, k) = cc(ido, k, 1)
  END DO
  RETURN
END SUBROUTINE radf2

SUBROUTINE radf3(ido, l1, cc, ch, wa1, wa2)
! Calculate the fast Fourier transform of subvectors of
! length three.
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, l1, 3)
  REAL(DP), INTENT(OUT) :: ch(ido, 3, l1)
  REAL(DP), INTENT(IN) :: wa1(*)
  REAL(DP), INTENT(IN) :: wa2(*)

  REAL(DP) :: ci2, cr2, di2, di3, dr2, dr3
  REAL(DP) :: taur, taui, ti2, ti3, tr2, tr3
  INTEGER:: i, ic, idp2, k
  DATA taur, taui/-.5D0, .866025403784439D0/

  DO k = 1, l1
    cr2 = cc(1, k, 2) + cc(1, k, 3)
    ch(1, 1, k) = cc(1, k, 1) + cr2
    ch(1, 3, k) = taui*(cc(1, k, 3) - cc(1, k, 2))
    ch(ido, 2, k) = cc(1, k, 1) + taur*cr2
  END DO
  IF (ido == 1) RETURN
  idp2 = ido + 2
  DO k = 1, l1
    DO i = 3, ido, 2
      ic = idp2 - i
      dr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
      di2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
      dr3 = wa2(i - 2)*cc(i - 1, k, 3) + wa2(i - 1)*cc(i, k, 3)
      di3 = wa2(i - 2)*cc(i, k, 3) - wa2(i - 1)*cc(i - 1, k, 3)
      cr2 = dr2 + dr3
      ci2 = di2 + di3
      ch(i - 1, 1, k) = cc(i - 1, k, 1) + cr2
      ch(i, 1, k) = cc(i, k, 1) + ci2
      tr2 = cc(i - 1, k, 1) + taur*cr2
      ti2 = cc(i, k, 1) + taur*ci2
      tr3 = taui*(di2 - di3)
      ti3 = taui*(dr3 - dr2)
      ch(i - 1, 3, k) = tr2 + tr3
      ch(ic - 1, 2, k) = tr2 - tr3
      ch(i, 3, k) = ti2 + ti3
      ch(ic, 2, k) = ti3 - ti2
    END DO
  END DO
  RETURN
END SUBROUTINE radf3

SUBROUTINE radf4(ido, l1, cc, ch, wa1, wa2, wa3)
! Calculate the fast Fourier transform of subvectors of
! length four.
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, l1, 4)
  REAL(DP), INTENT(OUT) :: ch(ido, 4, l1)
  REAL(DP), INTENT(IN) :: wa1(*)
  REAL(DP), INTENT(IN) :: wa2(*)
  REAL(DP), INTENT(IN) :: wa3(*)

  INTEGER:: i, ic, idp2, k
  REAL(DP):: ci2, ci3, ci4, cr2, cr3, cr4
  REAL(DP):: hsqt2
  REAL(DP):: ti1, ti2, ti3, ti4
  REAL(DP):: tr1, tr2, tr3, tr4

  DATA hsqt2/.7071067811865475D0/

  DO k = 1, l1
    tr1 = cc(1, k, 2) + cc(1, k, 4)
    tr2 = cc(1, k, 1) + cc(1, k, 3)
    ch(1, 1, k) = tr1 + tr2
    ch(ido, 4, k) = tr2 - tr1
    ch(ido, 2, k) = cc(1, k, 1) - cc(1, k, 3)
    ch(1, 3, k) = cc(1, k, 4) - cc(1, k, 2)
  END DO

  IF (ido < 2) THEN
    RETURN
  ELSE IF (ido > 2) THEN
    idp2 = ido + 2
    DO k = 1, l1
      DO i = 3, ido, 2
        ic = idp2 - i
        cr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
        ci2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
        cr3 = wa2(i - 2)*cc(i - 1, k, 3) + wa2(i - 1)*cc(i, k, 3)
        ci3 = wa2(i - 2)*cc(i, k, 3) - wa2(i - 1)*cc(i - 1, k, 3)
        cr4 = wa3(i - 2)*cc(i - 1, k, 4) + wa3(i - 1)*cc(i, k, 4)
        ci4 = wa3(i - 2)*cc(i, k, 4) - wa3(i - 1)*cc(i - 1, k, 4)
        tr1 = cr2 + cr4
        tr4 = cr4 - cr2
        ti1 = ci2 + ci4
        ti4 = ci2 - ci4
        ti2 = cc(i, k, 1) + ci3
        ti3 = cc(i, k, 1) - ci3
        tr2 = cc(i - 1, k, 1) + cr3
        tr3 = cc(i - 1, k, 1) - cr3
        ch(i - 1, 1, k) = tr1 + tr2
        ch(ic - 1, 4, k) = tr2 - tr1
        ch(i, 1, k) = ti1 + ti2
        ch(ic, 4, k) = ti1 - ti2
        ch(i - 1, 3, k) = ti4 + tr3
        ch(ic - 1, 2, k) = tr3 - ti4
        ch(i, 3, k) = tr4 + ti3
        ch(ic, 2, k) = tr4 - ti3
      END DO
    END DO
    IF (MOD(ido, 2) == 1) RETURN
  END IF

  DO k = 1, l1
    ti1 = -hsqt2*(cc(ido, k, 2) + cc(ido, k, 4))
    tr1 = hsqt2*(cc(ido, k, 2) - cc(ido, k, 4))
    ch(ido, 1, k) = tr1 + cc(ido, k, 1)
    ch(ido, 3, k) = cc(ido, k, 1) - tr1
    ch(1, 2, k) = ti1 - cc(ido, k, 3)
    ch(1, 4, k) = ti1 + cc(ido, k, 3)
  END DO
  RETURN
END SUBROUTINE radf4

SUBROUTINE radf5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
! Calculate the fast Fourier transform of subvectors of
! length five.
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, l1, 5)
  REAL(DP), INTENT(OUT) :: ch(ido, 5, l1)
  REAL(DP), INTENT(IN) :: wa1(*)
  REAL(DP), INTENT(IN) :: wa2(*)
  REAL(DP), INTENT(IN) :: wa3(*)
  REAL(DP), INTENT(IN) :: wa4(*)

  INTEGER:: i, ic, idp2, k
  REAL(DP):: ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5
  REAL(DP):: di2, di3, di4, di5, dr2, dr3, dr4, dr5
  REAL(DP):: ti11, ti12, ti2, ti3, ti4, ti5
  REAL(DP):: tr11, tr12, tr2, tr3, tr4, tr5

  DATA tr11, ti11, tr12, ti12/.309016994374947D0, .951056516295154D0, &
    -.809016994374947D0, .587785252292473D0/

  DO k = 1, l1
    cr2 = cc(1, k, 5) + cc(1, k, 2)
    ci5 = cc(1, k, 5) - cc(1, k, 2)
    cr3 = cc(1, k, 4) + cc(1, k, 3)
    ci4 = cc(1, k, 4) - cc(1, k, 3)
    ch(1, 1, k) = cc(1, k, 1) + cr2 + cr3
    ch(ido, 2, k) = cc(1, k, 1) + tr11*cr2 + tr12*cr3
    ch(1, 3, k) = ti11*ci5 + ti12*ci4
    ch(ido, 4, k) = cc(1, k, 1) + tr12*cr2 + tr11*cr3
    ch(1, 5, k) = ti12*ci5 - ti11*ci4
  END DO
  IF (ido == 1) RETURN
  idp2 = ido + 2
  DO k = 1, l1
    DO i = 3, ido, 2
      ic = idp2 - i
      dr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
      di2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
      dr3 = wa2(i - 2)*cc(i - 1, k, 3) + wa2(i - 1)*cc(i, k, 3)
      di3 = wa2(i - 2)*cc(i, k, 3) - wa2(i - 1)*cc(i - 1, k, 3)
      dr4 = wa3(i - 2)*cc(i - 1, k, 4) + wa3(i - 1)*cc(i, k, 4)
      di4 = wa3(i - 2)*cc(i, k, 4) - wa3(i - 1)*cc(i - 1, k, 4)
      dr5 = wa4(i - 2)*cc(i - 1, k, 5) + wa4(i - 1)*cc(i, k, 5)
      di5 = wa4(i - 2)*cc(i, k, 5) - wa4(i - 1)*cc(i - 1, k, 5)
      cr2 = dr2 + dr5
      ci5 = dr5 - dr2
      cr5 = di2 - di5
      ci2 = di2 + di5
      cr3 = dr3 + dr4
      ci4 = dr4 - dr3
      cr4 = di3 - di4
      ci3 = di3 + di4
      ch(i - 1, 1, k) = cc(i - 1, k, 1) + cr2 + cr3
      ch(i, 1, k) = cc(i, k, 1) + ci2 + ci3
      tr2 = cc(i - 1, k, 1) + tr11*cr2 + tr12*cr3
      ti2 = cc(i, k, 1) + tr11*ci2 + tr12*ci3
      tr3 = cc(i - 1, k, 1) + tr12*cr2 + tr11*cr3
      ti3 = cc(i, k, 1) + tr12*ci2 + tr11*ci3
      tr5 = ti11*cr5 + ti12*cr4
      ti5 = ti11*ci5 + ti12*ci4
      tr4 = ti12*cr5 - ti11*cr4
      ti4 = ti12*ci5 - ti11*ci4
      ch(i - 1, 3, k) = tr2 + tr5
      ch(ic - 1, 2, k) = tr2 - tr5
      ch(i, 3, k) = ti2 + ti5
      ch(ic, 2, k) = ti5 - ti2
      ch(i - 1, 5, k) = tr3 + tr4
      ch(ic - 1, 4, k) = tr3 - tr4
      ch(i, 5, k) = ti3 + ti4
      ch(ic, 4, k) = ti4 - ti3
    END DO
  END DO
  RETURN
END SUBROUTINE radf5

SUBROUTINE radfg(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
! Calculate the fast Fourier transform of subvectors of
! arbitrary length.
  USE params, ONLY: DP, tPI
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: ip
  INTEGER, INTENT(IN) :: l1
  INTEGER, INTENT(IN) :: idl1
  REAL(DP), INTENT(OUT) :: cc(ido, ip, l1)
  REAL(DP), INTENT(IN OUT) :: c1(ido, l1, ip)
  REAL(DP), INTENT(IN OUT) :: c2(idl1, ip)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, ip)
  REAL(DP), INTENT(OUT) :: ch2(idl1, ip)
  REAL(DP), INTENT(IN) :: wa(*)

  INTEGER:: i, ic, idij, idp2, ik, ipp2, ipph, is
  INTEGER:: j, j2, jc, k, l, lc, nbd
  REAL(DP):: ai1, ai2, ar1, ar1h, ar2, ar2h, arg, dc2, dcp, dsp, ds2

  arg = tpi/REAL(ip, DP)
  dcp = COS(arg)
  dsp = SIN(arg)
  ipph = (ip + 1)/2
  ipp2 = ip + 2
  idp2 = ido + 2
  nbd = (ido - 1)/2
  IF (ido == 1) THEN
    DO ik = 1, idl1
      c2(ik, 1) = ch2(ik, 1)
    END DO
  ELSE
    DO ik = 1, idl1
      ch2(ik, 1) = c2(ik, 1)
    END DO
    DO j = 2, ip
      DO k = 1, l1
        ch(1, k, j) = c1(1, k, j)
      END DO
    END DO
    IF (nbd > l1) THEN
      is = -ido
      DO j = 2, ip
        is = is + ido
        DO k = 1, l1
          idij = is
          DO i = 3, ido, 2
            idij = idij + 2
            ch(i - 1, k, j) = wa(idij - 1)*c1(i - 1, k, j) + wa(idij)*c1(i, k, j)
            ch(i, k, j) = wa(idij - 1)*c1(i, k, j) - wa(idij)*c1(i - 1, k, j)
          END DO
        END DO
      END DO
    ELSE
      is = -ido
      DO j = 2, ip
        is = is + ido
        idij = is
        DO i = 3, ido, 2
          idij = idij + 2
          DO k = 1, l1
            ch(i - 1, k, j) = wa(idij - 1)*c1(i - 1, k, j) + wa(idij)*c1(i, k, j)
            ch(i, k, j) = wa(idij - 1)*c1(i, k, j) - wa(idij)*c1(i - 1, k, j)
          END DO
        END DO
      END DO
    END IF

    IF (nbd < l1) THEN
      DO j = 2, ipph
        jc = ipp2 - j
        DO i = 3, ido, 2
          DO k = 1, l1
            c1(i - 1, k, j) = ch(i - 1, k, j) + ch(i - 1, k, jc)
            c1(i - 1, k, jc) = ch(i, k, j) - ch(i, k, jc)
            c1(i, k, j) = ch(i, k, j) + ch(i, k, jc)
            c1(i, k, jc) = ch(i - 1, k, jc) - ch(i - 1, k, j)
          END DO
        END DO
      END DO
    ELSE
      DO j = 2, ipph
        jc = ipp2 - j
        DO k = 1, l1
          DO i = 3, ido, 2
            c1(i - 1, k, j) = ch(i - 1, k, j) + ch(i - 1, k, jc)
            c1(i - 1, k, jc) = ch(i, k, j) - ch(i, k, jc)
            c1(i, k, j) = ch(i, k, j) + ch(i, k, jc)
            c1(i, k, jc) = ch(i - 1, k, jc) - ch(i - 1, k, j)
          END DO
        END DO
      END DO
    END IF
  END IF

  DO j = 2, ipph
    jc = ipp2 - j
    DO k = 1, l1
      c1(1, k, j) = ch(1, k, j) + ch(1, k, jc)
      c1(1, k, jc) = ch(1, k, jc) - ch(1, k, j)
    END DO
  END DO

  ar1 = 1D0
  ai1 = 0D0
  DO l = 2, ipph
    lc = ipp2 - l
    ar1h = dcp*ar1 - dsp*ai1
    ai1 = dcp*ai1 + dsp*ar1
    ar1 = ar1h
    DO ik = 1, idl1
      ch2(ik, l) = c2(ik, 1) + ar1*c2(ik, 2)
      ch2(ik, lc) = ai1*c2(ik, ip)
    END DO
    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1
    DO j = 3, ipph
      jc = ipp2 - j
      ar2h = dc2*ar2 - ds2*ai2
      ai2 = dc2*ai2 + ds2*ar2
      ar2 = ar2h
      DO ik = 1, idl1
        ch2(ik, l) = ch2(ik, l) + ar2*c2(ik, j)
        ch2(ik, lc) = ch2(ik, lc) + ai2*c2(ik, jc)
      END DO
    END DO
  END DO

  DO j = 2, ipph
    DO ik = 1, idl1
      ch2(ik, 1) = ch2(ik, 1) + c2(ik, j)
    END DO
  END DO

  IF (ido < l1) THEN
    DO i = 1, ido
      DO k = 1, l1
        cc(i, 1, k) = ch(i, k, 1)
      END DO
    END DO
  ELSE
    DO k = 1, l1
      DO i = 1, ido
        cc(i, 1, k) = ch(i, k, 1)
      END DO
    END DO
  END IF

  DO j = 2, ipph
    jc = ipp2 - j
    j2 = j + j
    DO k = 1, l1
      cc(ido, j2 - 2, k) = ch(1, k, j)
      cc(1, j2 - 1, k) = ch(1, k, jc)
    END DO
  END DO

  IF (ido == 1) RETURN

  IF (nbd < l1) THEN
    DO j = 2, ipph
      jc = ipp2 - j
      j2 = j + j
      DO i = 3, ido, 2
        ic = idp2 - i
        DO k = 1, l1
          cc(i - 1, j2 - 1, k) = ch(i - 1, k, j) + ch(i - 1, k, jc)
          cc(ic - 1, j2 - 2, k) = ch(i - 1, k, j) - ch(i - 1, k, jc)
          cc(i, j2 - 1, k) = ch(i, k, j) + ch(i, k, jc)
          cc(ic, j2 - 2, k) = ch(i, k, jc) - ch(i, k, j)
        END DO
      END DO
    END DO
  ELSE
    DO j = 2, ipph
      jc = ipp2 - j
      j2 = j + j
      DO k = 1, l1
        DO i = 3, ido, 2
          ic = idp2 - i
          cc(i - 1, j2 - 1, k) = ch(i - 1, k, j) + ch(i - 1, k, jc)
          cc(ic - 1, j2 - 2, k) = ch(i - 1, k, j) - ch(i - 1, k, jc)
          cc(i, j2 - 1, k) = ch(i, k, j) + ch(i, k, jc)
          cc(ic, j2 - 2, k) = ch(i, k, jc) - ch(i, k, j)
        END DO
      END DO
    END DO
  END IF

  RETURN

END SUBROUTINE radfg

SUBROUTINE rfftf(n, r, wsave, ifac)
! Calls for rfftf1 IF n is not 1
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN OUT) :: r(*)
  REAL(DP), INTENT(IN OUT) :: wsave(*)
  INTEGER, INTENT(IN) ::ifac(10)

  IF (n /= 1) CALL rfftf1(n, r, wsave, wsave(n + 1), ifac)

  RETURN
END SUBROUTINE rfftf

SUBROUTINE rfftf1(n, c, ch, wa, ifac)
! SUBROUTINE RFFTF1 computes the Fourier coefficients of a REAL
! periodic SEQUENCE (Fourier analysis). The transform is defined
! below at output PARAMETER C.
!
! The WA and IFAC arrays which are used by SUBROUTINE RFFTB1 must be
! initialized by calling SUBROUTINE RFFTI1.
!
! Input Arguments
!
! N the length of the array C TO be transformed. The method
! is most efficient when N is a product of small primes.
! N may change so long as different work arrays are provided.
!
! C a REAL array of length N which CONTAINS the SEQUENCE
! TO be transformed.
!
! CH a REAL work array of length at least N.
!
! WA a REAL work array which must be dimensioned at least N.
!
! IFAC an INTEGER work array which must be dimensioned at least 15.
!
! The WA and IFAC arrays must be initialized by calling
! SUBROUTINE RFFTI1, and different WA and IFAC arrays must be
! used for each different value of N. This initialization
! does not have TO be repeated so long as N remains unchanged.
! Thus subsequent transforms can be obtained faster than the
! first. The same WA and IFAC arrays can be used by RFFTF1
! and RFFTB1.
!
! Output Argument
!
! C C(1) = the sum from I=1 TO I=N of R(I)
!
! IF N is even set L = N/2; IF N is odd set L = (N+1)/2
!
! THEN for K = 2,...,L
!
! C(2*K-2) = the sum from I = 1 TO I = N of
!
! C(I)*COS((K-1)*(I-1)*2*PI/N)
!
! C(2*K-1) = the sum from I = 1 TO I = N of
!
! -C(I)*SIN((K-1)*(I-1)*2*PI/N)
!
! IF N is even
!
! C(N) = the sum from I = 1 TO I = N of
!
! (-1)**(I-1)*C(I)
!
! Notes: This transform is unnormalized since a CALL of RFFTF1
! followed by a CALL of RFFTB1 will multiply the input
! SEQUENCE by N.
!
! WA and IFAC contain initialization calculations which must
! not be destroyed between calls of SUBROUTINE RFFTF1 or
! RFFTB1.
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(OUT) :: c(n)
  REAL(DP), INTENT(IN OUT) :: ch(*)
  REAL(DP), INTENT(IN) :: wa(*)
  INTEGER, INTENT(IN) :: ifac(10)

  INTEGER:: i, idl1, ido, ip, iw, ix2, ix3, ix4
  INTEGER:: k1, kh, l1, l2, na, nf

  nf = ifac(2)
  na = 1
  l2 = n
  iw = n
  DO k1 = 1, nf
    kh = nf - k1
    ip = ifac(kh + 3)
    l1 = l2/ip
    ido = n/l2
    idl1 = ido*l1
    iw = iw - (ip - 1)*ido
    na = 1 - na
    !-------
    SELECT CASE (ip)
      !-------
    CASE (4)
      ix2 = iw + ido
      ix3 = ix2 + ido
      IF (na /= 0) THEN
        CALL radf4(ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3))
      ELSE
        CALL radf4(ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3))
      END IF
      !-------
    CASE (2)
      IF (na /= 0) THEN
        CALL radf2(ido, l1, ch, c, wa(iw))
      ELSE
        CALL radf2(ido, l1, c, ch, wa(iw))
      END IF
      !-------
    CASE (3)
      ix2 = iw + ido
      IF (na /= 0) THEN
        CALL radf3(ido, l1, ch, c, wa(iw), wa(ix2))
      ELSE
        CALL radf3(ido, l1, c, ch, wa(iw), wa(ix2))
      END IF
      !-------
    CASE (5)
      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido
      IF (na /= 0) THEN
        CALL radf5(ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4))
      ELSE
        CALL radf5(ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4))
      END IF
      !-------
    CASE DEFAULT
      IF (ido == 1) na = 1 - na
      IF (na /= 0) THEN
        CALL radfg(ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw))
        na = 0
      ELSE
        CALL radfg(ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw))
        na = 1
      END IF
      !-------
    END SELECT
    !-------
    l2 = l1
  END DO
  IF (na == 1) RETURN
  DO i = 1, n
    c(i) = ch(i)
  END DO
  RETURN
END SUBROUTINE rfftf1

SUBROUTINE cosqi(n, wsave, ifac)
! SUBROUTINE COSQI initializes the work array WSAVE which is used IN
! both COSQF1 and COSQB1. The prime factorization of N together with
! a tabulation of the trigonometric functions are computed and
! stored IN WSAVE and IFAC.
!
! Input PARAMETER
!
! N the length of the array TO be transformed. The method
! is most efficient when N is a product of small primes.
!
! Output PARAMETER
!
! WSAVE a work array which must be dimensioned at least 3*N.
! The same work array can be used for both COSQF1 and COSQB1
! as long as N remains unchanged. Different WSAVE arrays
! are required for different values of N. The contents of
! WSAVE must not be changed between calls of COSQF1 or COSQB1.
!
! IFAC an INTEGER work array which must be dimensioned at least 10.
! Condtions applying on WSAVE also apply on IFAC.
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(OUT) :: wsave(*)
  INTEGER, INTENT(OUT) :: ifac(10)

  INTEGER:: k
  REAL(DP):: dt, fk, pih

  DATA pih/1.57079632679491D0/

  dt = pih/REAL(n, DP)
  fk = 0D0
  DO k = 1, n
    fk = fk + 1.0D0
    wsave(k) = COS(fk*dt)
  END DO
  CALL rffti(n, wsave(n + 1), ifac)
  RETURN
END SUBROUTINE cosqi

SUBROUTINE cost(n, x, wsave, ifac)
  USE params, ONLY: DP
  IMPLICIT NONE

! SUBROUTINE COST computes the discrete Fourier cosine transform
! of an even SEQUENCE X(I). The transform is defined below at output
! PARAMETER X.
!
! COST is the unnormalized inverse of itself since a CALL of COST
! followed by another CALL of COST will multiply the input SEQUENCE
! X by 2*(N-1). The transform is defined below at output PARAMETER X.
!
! The WSAVE and IFAC arrays which are used by SUBROUTINE COST must be
! initialized by calling SUBROUTINE COSTI(N,WSAVE,IFAC).
!
! Input Parameters
!
! N the length of the SEQUENCE X. N must be greater than 1.
! The method is most efficient when N-1 is a product of
! small primes.
!
! X an array which CONTAINS the SEQUENCE TO be transformed
!
! WSAVE a work array which must be dimensioned at least 3*N
! IN the PROGRAM that calls COST. The WSAVE array must be
! initialized by calling SUBROUTINE COSTI(N,WSAVE), and a
! different WSAVE array must be used for each different
! value of N. This initialization does not have TO be
! repeated so long as N remains unchanged. Thus subsequent
! transforms can be obtained faster than the first.
!
! IFAC an INTEGER work array which must be dimensioned at least 15.
! Condtions applying on WSAVE also apply on IFAC.

! Output Parameters
!
! X For I=1,...,N
!
! X(I) = X(1)+(-1)**(I-1)*X(N)
!
! + the sum from K=2 TO K=N-1
!
! 2*X(K)*COS((K-1)*(I-1)*PI/(N-1))
!
! A CALL of COST followed by another CALL of
! COST will multiply the SEQUENCE X by 2*(N-1).
! Hence COST is the unnormalized inverse
! of itself.
!
! WSAVE and IFAC contain initialization calculations which must not be
! destroyed between calls of COST.
!

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN OUT) :: x(n)
  REAL(DP), INTENT(IN OUT) :: wsave(*)
  INTEGER, INTENT(IN) :: ifac(10)

  INTEGER::i, k, kc, modn, nm1, np1, ns2
  REAL(DP):: c1, t1, t2, tx2, x1h, x1p3, xi, xim2

  nm1 = n - 1
  np1 = n + 1
  ns2 = n/2
!-------
  SELECT CASE (n)
!-------
  CASE (:1)
    STOP 'IN cost : vector of SIZE n too small for Fourier transform'
!-------
  CASE (2)
    x1h = x(1) + x(2)
    x(2) = x(1) - x(2)
    x(1) = x1h
!-------
  CASE (3)
    x1p3 = x(1) + x(3)
    tx2 = x(2) + x(2)
    x(2) = x(1) - x(3)
    x(1) = x1p3 + tx2
    x(3) = x1p3 - tx2
!-------
  CASE DEFAULT
    c1 = x(1) - x(n)
    x(1) = x(1) + x(n)
    DO k = 2, ns2
      kc = np1 - k
      t1 = x(k) + x(kc)
      t2 = x(k) - x(kc)
      c1 = c1 + wsave(kc)*t2
      t2 = wsave(k)*t2
      x(k) = t1 - t2
      x(kc) = t1 + t2
    END DO
    modn = MOD(n, 2)
    IF (modn /= 0) x(ns2 + 1) = x(ns2 + 1) + x(ns2 + 1)
    CALL rfftf(nm1, x, wsave(n + 1), ifac)
    xim2 = x(2)
    x(2) = c1
    DO i = 4, n, 2
      xi = x(i)
      x(i) = x(i - 2) - x(i - 1)
      x(i - 1) = xim2
      xim2 = xi
    END DO
    IF (modn /= 0) x(n) = xim2
!-------
  END SELECT
!-------
  RETURN
END SUBROUTINE cost

SUBROUTINE costi(n, wsave, ifac)
! SUBROUTINE COSTI initializes array WSAVE and IFAC which is used IN
! SUBROUTINE COST. The prime factorization of N together with
! a tabulation of the trigonometric functions are computed and
! stored IN WSAVE and IFAC.
!
! Input PARAMETER
!
! N the length of the SEQUENCE TO be transformed. The method
! is most efficient when N-1 is a product of small primes.
!
! Output PARAMETER
!
! WSAVE a work array which must be dimensioned at least 3*N.
! Different WSAVE arrays are required for different values
! of N. The contents of WSAVE must not be changed between
! calls of COST.
!
! IFAC an INTEGER work array which must be dimensioned at least 15.
! Different IFAC arrays are required for different values
! of N. The contents of IFAC must not be changed between
! calls of COST.
  USE params, ONLY: DP, PI
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(OUT) :: wsave(*)
  INTEGER, INTENT(OUT) :: ifac(10)

  INTEGER:: k, kc, nm1, np1, ns2
  REAL(DP):: dt, fk

  IF (n <= 3) RETURN
  nm1 = n - 1
  np1 = n + 1
  ns2 = n/2
  dt = pi/REAL(nm1, DP)
  fk = 0D0
  DO k = 2, ns2
    kc = np1 - k
    fk = fk + 1D0
    wsave(k) = 2D0*SIN(fk*dt)
    wsave(kc) = 2D0*COS(fk*dt)
  END DO
  CALL rffti(nm1, wsave(n + 1), ifac)
  RETURN
END SUBROUTINE costi

SUBROUTINE rffti(n, wsave, ifac)
! Calls rffti1 IF n is not 1
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(OUT) :: wsave(*)
  INTEGER, INTENT(OUT) :: ifac(10)

  IF (n == 1) RETURN
  CALL rffti1(n, wsave(n + 1), ifac)
  RETURN
END SUBROUTINE rffti

SUBROUTINE sinqb(n, x, wsave, ifac)
! SUBROUTINE SINQB computes the fast Fourier transform of quarter
! wave DATA. That is, SINQB computes a SEQUENCE from its
! representation IN terms of a sine series with odd wave numbers.
! the transform is defined below at output PARAMETER X.
!
! SINQF is the unnormalized inverse of SINQB since a CALL of SINQB
! followed by a CALL of SINQF will multiply the input SEQUENCE X
! by 4*N.
!
! The WSAVE and IFAC arrays which are used by SUBROUTINE SINQB must be
! initialized by calling SUBROUTINE SINQI(N,WSAVE,IFAC).
!
! Input Parameters
!
! N the length of the array X TO be transformed. The method
! is most efficient when N is a product of small primes.
!
! X an array which CONTAINS the SEQUENCE TO be transformed
!
! WSAVE a work array which must be dimensioned at least 3*N+15
! IN the PROGRAM that calls SINQB. The WSAVE array must be
! initialized by calling SUBROUTINE SINQI(N,WSAVE), and a
! different WSAVE array must be used for each different
! value of N. This initialization does not have TO be
! repeated so long as N remains unchanged. Thus subsequent
! transforms can be obtained faster than the first.
!
! IFAC an INTEGER work array which must be dimensioned at least 15.
! Condtions applying on WSAVE also apply on IFAC.
!
! Output Parameters
!
! X For I=1,...,N
!
! X(I)= the sum from K=1 TO K=N of
!
! 4*X(K)*SIN((2*K-1)*I*PI/(2*N))
!
! a CALL of SINQB followed by a CALL of
! SINQF will multiply the SEQUENCE X by 4*N.
! Therefore SINQF is the unnormalized inverse
! of SINQB.
!
! WSAVE and IFAC contain initialization calculations which must not
! be destroyed between calls of SINQB or SINQF.
!
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(OUT) :: x(n)
  REAL(DP), INTENT(IN OUT) :: wsave(*)
  INTEGER, INTENT(IN OUT) :: ifac(10)

  INTEGER:: k, kc, ns2
  REAL(DP):: xhold

  IF (n > 1) THEN
    ns2 = n/2
    DO k = 2, n, 2
      x(k) = -x(k)
    END DO
    CALL cosqb(n, x, wsave, ifac)
    DO k = 1, ns2
      kc = n - k
      xhold = x(k)
      x(k) = x(kc + 1)
      x(kc + 1) = xhold
    END DO
  ELSE
    x(1) = 4D0*x(1)
  END IF

  RETURN

END SUBROUTINE sinqb

SUBROUTINE sinqf(n, x, wsave, ifac)
! SUBROUTINE SINQF computes the fast Fourier transform of quarter
! wave DATA. That is, SINQF computes the coefficients IN a sine
! series representation with ONLY odd wave numbers. The transform
! is defined below at output PARAMETER X.
!
! SINQB is the unnormalized inverse of SINQF since a CALL of SINQF
! followed by a CALL of SINQB will multiply the input SEQUENCE X
! by 4*N.
!
! The WSAVE and IFAC arrays which are used by SUBROUTINE SINQF must be
! initialized by calling SUBROUTINE SINQI(N,WSAVE,IFAC).
!
! Input Parameters
!
! N the length of the array X TO be transformed. The method
! is most efficient when N is a product of small primes.
!
! X an array which CONTAINS the SEQUENCE TO be transformed
!
! WSAVE a work array which must be dimensioned at least 3*N+15
! IN the PROGRAM that calls SINQF. The WSAVE array must be
! initialized by calling SUBROUTINE SINQI(N,WSAVE), and a
! different WSAVE array must be used for each different
! value of N. This initialization does not have TO be
! repeated so long as N remains unchanged. Thus subsequent
! transforms can be obtained faster than the first.
!
!
! IFAC an INTEGER work array which must be dimensioned at least 15.
! Condtions applying on WSAVE also apply on IFAC.
!
! Output Parameters
!
! X For I=1,...,N
!
! X(I) = (-1)**(I-1)*X(N)
!
! + the sum from K=1 TO K=N-1 of
!
! 2*X(K)*SIN((2*I-1)*K*PI/(2*N))
!
! A CALL of SINQF followed by a CALL of
! SINQB will multiply the SEQUENCE X by 4*N.
! Therefore SINQB is the unnormalized inverse
! of SINQF.
!
! WSAVE and IFAC CONTAINS initialization calculations which must not
! be destroyed between calls of SINQF or SINQB.
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN OUT) :: x(n)
  REAL(DP), INTENT(IN OUT) :: wsave(*)
  INTEGER, INTENT(IN OUT) :: ifac(10)

  INTEGER:: k, kc, ns2
  REAL(DP):: xhold

  IF (n == 1) RETURN
  ns2 = n/2
  DO k = 1, ns2
    kc = n - k
    xhold = x(k)
    x(k) = x(kc + 1)
    x(kc + 1) = xhold
  END DO
  CALL cosqf(n, x, wsave, ifac)
  DO k = 2, n, 2
    x(k) = -x(k)
  END DO
  RETURN
END SUBROUTINE sinqf

SUBROUTINE sinqi(n, wsave, ifac)
! Initialize work arrays for SINQF and SINQB
! Same as COSQI
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN OUT) :: n
  REAL(DP), INTENT(IN OUT) :: wsave(*)
  INTEGER, INTENT(IN OUT) :: ifac(10)

  CALL cosqi(n, wsave, ifac)
  RETURN
END SUBROUTINE sinqi

SUBROUTINE sint(n, x, wsave, ifac)
! SUBROUTINE SINT computes the discrete Fourier sine transform
! of an odd SEQUENCE X(I). The transform is defined below at
! output PARAMETER X.
!
! SINT is the unnormalized inverse of itself since a CALL of SINT
! followed by another CALL of SINT will multiply the input SEQUENCE
! X by 2*(N+1).
!
! The WSAVE and IFAC arrays which are used by SUBROUTINE SINT must be
! initialized by calling SUBROUTINE SINTI(N,WSAVE,IFAC).
!
! Input Parameters
!
! N the length of the SEQUENCE TO be transformed. The method
! is most efficient when N+1 is the product of small primes.
!
! X an array which CONTAINS the SEQUENCE TO be transformed
!
!
! WSAVE a work array with DIMENSION at least INT(3.5*N+16)
! IN the PROGRAM that calls SINT. The WSAVE array must be
! initialized by calling SUBROUTINE SINTI(N,WSAVE), and a
! different WSAVE array must be used for each different
! value of N. This initialization does not have TO be
! repeated so long as N remains unchanged. Thus subsequent
! transforms can be obtained faster than the first.
!
! IFAC an INTEGER work array which must be dimensioned at least 15.
! Condtions applying on WSAVE also apply on IFAC.
!
! Output Parameters
!
! X For I=1,...,N
!
! X(I)= the sum from K=1 TO K=N
!
! 2*X(K)*SIN(K*I*PI/(N+1))
!
! A CALL of SINT followed by another CALL of
! SINT will multiply the SEQUENCE X by 2*(N+1).
! Hence SINT is the unnormalized inverse
! of itself.
!
! WSAVE and IFAC contain initialization calculations which must not be
! destroyed between calls of SINT.
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN OUT) :: x(n)
  REAL(DP), INTENT(IN OUT) :: wsave(*)
  INTEGER, INTENT(IN OUT) :: ifac(10)

  INTEGER:: iw1, iw2, iw3, np1

  np1 = n + 1
  iw1 = n/2 + 1
  iw2 = iw1 + np1
  iw3 = iw2 + np1
  CALL sint1(n, x, wsave, wsave(iw1), wsave(iw2), ifac)
  RETURN
END SUBROUTINE sint

SUBROUTINE sint1(n, war, was, xh, x, ifac)
! See sint
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN OUT) :: war(*)
  REAL(DP), INTENT(IN) :: was(*)
  REAL(DP), INTENT(OUT) :: xh(*)
  REAL(DP), INTENT(IN OUT) :: x(*)
  INTEGER, INTENT(IN OUT) :: ifac(10)

  INTEGER:: i, k, kc, modn, np1, ns2
  REAL(DP):: sqrt3, t1, t2, xhold

  DATA sqrt3/1.73205080756888D0/

  DO i = 1, n
    xh(i) = war(i)
    war(i) = x(i)
  END DO

  IF (n < 2) THEN
    xh(1) = xh(1) + xh(1)
  ELSE IF (n == 2) THEN
    xhold = sqrt3*(xh(1) + xh(2))
    xh(2) = sqrt3*(xh(1) - xh(2))
    xh(1) = xhold
  ELSE
    np1 = n + 1
    ns2 = n/2
    x(1) = 0D0
    DO k = 1, ns2
      kc = np1 - k
      t1 = xh(k) - xh(kc)
      t2 = was(k)*(xh(k) + xh(kc))
      x(k + 1) = t1 + t2
      x(kc + 1) = t2 - t1
    END DO
    modn = MOD(n, 2)
    IF (modn /= 0) x(ns2 + 2) = 4.*xh(ns2 + 1)
    CALL rfftf1(np1, x, xh, war, ifac)
    xh(1) = .5D0*x(1)
    DO i = 3, n, 2
      xh(i - 1) = -x(i)
      xh(i) = xh(i - 2) + x(i - 1)
    END DO
    IF (modn == 0) xh(n) = -x(n + 1)
  END IF

  DO i = 1, n
    x(i) = war(i)
    war(i) = xh(i)
  END DO

  RETURN
END SUBROUTINE sint1

SUBROUTINE sinti(n, wsave, ifac)
  USE params, ONLY: DP, PI
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(OUT) :: wsave(*)
  INTEGER, INTENT(OUT) :: ifac(10)

  INTEGER:: k, np1, ns2
  REAL(DP):: dt

  IF (n <= 1) RETURN
  ns2 = n/2
  np1 = n + 1
  dt = PI/REAL(np1, DP)
  DO k = 1, ns2
    wsave(k) = 2D0*SIN(k*dt)
  END DO
  CALL rffti(np1, wsave(ns2 + 1), ifac)
  RETURN
END SUBROUTINE sinti

SUBROUTINE rffti1(n, wa, ifac)
! SUBROUTINE RFFTI1 initializes the work arrays WA and IFAC which are
! used IN both RFFTF1 and RFFTB1. The prime factorization of N and a
! tabulation of the trigonometric functions are computed and stored IN
! IFAC and WA, respectively.
!
! Input Argument
!
! N the length of the SEQUENCE TO be transformed.
!
! Output Arguments
!
! WA a REAL work array which must be dimensioned at least N.
!
! IFAC an INTEGER work array which must be dimensioned at least 15.
!
! The same work arrays can be used for both RFFTF1 and RFFTB1 as long
! as N remains unchanged. Different WA and IFAC arrays are required
! for different values of N. The contents of WA and IFAC must not be
! changed between calls of RFFTF1 or RFFTB1.
  USE params, ONLY: DP, tPI
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(OUT) :: wa(*)
  INTEGER, INTENT(OUT) :: ifac(10)
  INTEGER :: ntryh(4)

  INTEGER:: i, ib, ido, ii, ip, ipm, is
  INTEGER:: j, k1, l1, l2, ld, nf, nfm1, nl, nq, nr, ntry

  REAL(DP):: arg, argh, argld, fi

  DATA ntryh(1), ntryh(2), ntryh(3), ntryh(4)/4, 2, 3, 5/

  nl = n
  nf = 0
  j = 1
  ntry = ntryh(1)

  DO WHILE (nl /= 1)
    nq = nl/ntry
    nr = nl - ntry*nq
    IF (nr /= 0) THEN
      j = j + 1
      IF (j <= 4) THEN
        ntry = ntryh(j)
      ELSE
        ntry = ntry + 2
      END IF
      CYCLE
    END IF

    nf = nf + 1
    ifac(nf + 2) = ntry
    nl = nq
    IF (ntry == 2 .AND. nf /= 1) THEN
      DO i = 2, nf
        ib = nf - i + 2
        ifac(ib + 2) = ifac(ib + 1)
      END DO
      ifac(3) = 2
    END IF
  END DO

  ifac(1) = n
  ifac(2) = nf
  argh = tpi/REAL(n, DP)
  is = 0
  nfm1 = nf - 1
  l1 = 1
  IF (nfm1 == 0) RETURN
  DO k1 = 1, nfm1
    ip = ifac(k1 + 2)
    ld = 0
    l2 = l1*ip
    ido = n/l2
    ipm = ip - 1
    DO j = 1, ipm
      ld = ld + l1
      i = is
      argld = REAL(ld, DP)*argh
      fi = 0D0
      DO ii = 3, ido, 2
        i = i + 2
        fi = fi + 1D0
        arg = fi*argld
        wa(i - 1) = COS(arg)
        wa(i) = SIN(arg)
      END DO
      is = is + ido
    END DO
    l1 = l2
  END DO
  RETURN
END SUBROUTINE rffti1

