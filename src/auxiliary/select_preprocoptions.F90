PROGRAM select_preprocoptions

! Routine to extract wanted subset of the code in QDD.
! The preprocessor option to be evaluated are set here in the code.
! See the lines where 'option(...)' is set.
! The routine eevaluates one option after another using a temporary
! work-FILE 'xxxtmp'.

! The files for which the preprocessor options are evaluated are
! also given here in the code.

! How to USE it:
! 1. Go to the subdirectory 'src/qdd' of the 'QDD' repository.
! 2. Make sure that a subdirectory 'public' exists.
! 3. Compile 'gfortran ../auxiliary/select_proprocoptions.F90.
! 4. Issue './a.out'.
! 5. Go to subdirectory 'public'. All *.f90 and *.F90 should be there.
! 6. Make sure that 'fftw3.f03' is accessible in this directory.
! 7. Try 'make -f Makefile.gnu'. You should finally find a working 'qdd'
! in this direcctory.

  IMPLICIT NONE

  INTEGER, PARAMETER :: maxpreproc = 20 ! maximum NUMBER of preprocessor options
  CHARACTER(20) :: option(maxpreproc), actoption
  INTEGER :: length_option(maxpreproc), len_actoption, numpreproc, n
  LOGICAL :: value_option(maxpreproc), tactvalue

  INTEGER, PARAMETER :: maxfiles = 100 ! maximum NUMBER of files to be treated
  CHARACTER(64) :: filename(maxfiles), actfile
  INTEGER :: len_filename(maxfiles), len_actfile, numfiles, nf

  INTEGER, PARAMETER :: maxcodelines = 99999
  CHARACTER(256) :: codeline
  INTEGER :: len_codeline, line, level

! copy *.f files into public directory because they DO not need preprocessing
  CALL system('cp *.f public/.')

! define files: input, output and work (preliminary)
  filename = ''
  filename(01) = 'abso_bc.F90'
  filename(02) = 'fft.F90'
  filename(03) = 'fftpack2.F90'
  filename(04) = 'fftpack.F90'
  filename(05) = 'fftw.F90'
  filename(06) = 'forces.F90'
  filename(07) = 'givens.F90'
  filename(08) = 'init.F90'
  filename(09) = 'ionmd.F90'
  filename(10) = 'kinetic.F90'
  filename(11) = 'lda.F90'
  filename(12) = 'localize.F90'
  filename(13) = 'loc_mfield.F90'
  filename(14) = 'main.F90'
  filename(15) = 'nonloc.F90'
  filename(16) = 'orthmat.F90'
  filename(17) = 'parallele.F90'
  filename(18) = 'params.F90'
  filename(19) = 'pseudo.F90'
  filename(20) = 'pseudogoed.F90'
  filename(21) = 'pseudosoft.F90'
!  filename(22) = 'QDD_info.F90'
  filename(22) = 'falr.F90'
  filename(23) = 'restart.F90'
  filename(24) = 'rho.F90'
  filename(25) = 'rta.F90'
  filename(26) = 'schmid.F90'
  filename(27) = 'sicnew.F90'
  filename(28) = 'static.F90'
  filename(29) = 'subgrids.F90'
  filename(30) = 'util.F90'
!  filename(31) = 'zeroforce.F90'
  filename(31) = 'carlo.F90'
  filename(32) = 'coulex.F90'
  filename(33) = 'coulsolv.F90'
  filename(34) = 'dynamic.F90'
  DO nf = 1, maxfiles
    IF (LEN_TRIM(filename(nf)) > 0) THEN
      numfiles = nf
      len_filename(nf) = LEN_TRIM(filename(nf))
    ELSE
      EXIT
    END IF
  END DO

  OPEN (3, FILE='xxxtmp')
  DO nf = 1, numfiles
    actfile = filename(nf)
    len_actfile = len_filename(nf)
    OPEN (1, FILE=actfile, ACTION='READ')
    OPEN (2, FILE='public/'//actfile)
    CALL system('cp '//actfile//' public/'//actfile)
    CLOSE (1)

    ! enter here preprocesseroptions to be selected
    option = ''
    option(1) = 'fftw_cpu'; value_option(1) = .TRUE.
    option(2) = 'netlib'; value_option(2) = .FALSE.
    option(3) = 'extended'; value_option(3) = .FALSE.
    option(4) = 'fsic'; value_option(4) = .FALSE.
    option(5) = 'raregas'; value_option(5) = .FALSE.
    option(6) = 'findiff'; value_option(6) = .FALSE.
    option(7) = 'numerov'; value_option(7) = .FALSE.
    option(8) = 'twostsic'; value_option(8) = .FALSE.
    option(9) = 'nompi'; value_option(9) = .TRUE.
    option(10) = 'mpi'; value_option(10) = .FALSE.    ! put 'mpi' last

    ! determine NUMBER of active options and their length
    DO n = 1, maxpreproc
      IF (LEN_TRIM(option(n)) > 0) THEN
        numpreproc = n
        length_option(n) = LEN_TRIM(option(n))
      ELSE
        EXIT
      END IF
    END DO

    ! run filter successively, USE scratch files for intermediate storage
    DO n = 1, numpreproc
      actoption = option(n)
      len_actoption = LEN_TRIM(actoption)
      tactvalue = .TRUE.
      level = 0
      REWIND (2)
      REWIND (3)
      WRITE (*, *) 'start run nr.', n, ' OPTION=', actoption
      DO line = 1, maxcodelines
        READ (2, '(a)', ERR=99, END=19) codeline
        len_codeline = LEN_TRIM(codeline)

        IF (len_codeline > 2 .AND. codeline(1:3) == '#if') THEN
          IF (INDEX(codeline(1:len_codeline), '&') > 0) &
            STOP 'code not suited to evaluate .AND. condition'
          ! WRITE(*,*) codeline(1:len_codeline),actoption,INDEX(codeline,actoption(1:len_actoption))
          IF (INDEX(codeline(1:len_codeline), '!'//actoption(1:len_actoption)) > 0) THEN
            tactvalue = .NOT. value_option(n)
            level = 1
            WRITE (*, *) 'NO in line:'//codeline(1:len_codeline)
          ELSE IF (INDEX(codeline(1:len_codeline), actoption(1:len_actoption)) > 0) THEN
            tactvalue = value_option(n)
            level = 1
            WRITE (*, *) 'YES in line:'//codeline(1:len_codeline)
          ELSE IF (level > 0) THEN
            level = level + 1
            IF (tactvalue) WRITE (3, '(a)') codeline(1:len_codeline)
          ELSE
            IF (tactvalue) WRITE (3, '(a)') codeline(1:len_codeline)
          END IF
        ELSE IF (level == 1 .AND. len_codeline > 4 .AND. codeline(1:5) == '#else') THEN
          tactvalue = .NOT. tactvalue
        ELSE IF (level > 0 .AND. len_codeline > 5 .AND. codeline(1:6) == '#endif') THEN
          level = level - 1
          IF (level == 0) THEN
            tactvalue = .TRUE.
          ELSE
            IF (tactvalue) WRITE (3, '(a)') codeline(1:len_codeline)
          END IF
        ELSE
          IF (tactvalue) WRITE (3, '(a)') codeline(1:len_codeline)
        END IF
      END DO
19    CONTINUE
      WRITE (*, *) 'EXIT loop after codeline nr. ', line, ' level=', level
      CALL system('cp xxxtmp'//' public/'//actfile)
    END DO

    CLOSE (2)
  END DO

  STOP 'regular END'
99 STOP ' error in input'

END PROGRAM select_preprocoptions
