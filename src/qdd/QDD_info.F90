MODULE QDDInfo
  USE iso_c_binding
#ifdef __INTEL_COMPILER
  USE ifport
#endif /* __INTEL_COMPILER */

  IMPLICIT NONE

  INTERFACE
    FUNCTION readlink(path, buf, bufsize) bind(C, NAME='readlink')
      IMPORT
      INTEGER(C_SIZE_T) :: readlink
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: path(*)
      CHARACTER(KIND=C_CHAR) :: buf(*)
      INTEGER(C_SIZE_T), value :: bufsize
    END FUNCTION
  END INTERFACE

  INTEGER :: pid, i, j, idx, nFilesChanged
  INTEGER(C_SIZE_T) :: szret
  CHARACTER(256) :: path, filesChanged, filename = 'QDD.info'
  CHARACTER(256), ALLOCATABLE :: files(:)
  CHARACTER(2) :: snFilesChanged
  CHARACTER(KIND=C_CHAR) :: cbuf(256)

CONTAINS

  SUBROUTINE writeQDDVersionInfo()
    filesChanged = FILES_CHANGED
    snFilesChanged = NFILES_CHANGED
    READ (snFilesChanged, *) nFilesChanged
    nFilesChanged = nFilesChanged*2
    ALLOCATE (files(nFilesChanged))
    READ (filesChanged, *) files(1:nFilesChanged)

    OPEN (UNIT=100, FILE=trim(filename))
    WRITE (100, "('********************************************************************************')")
    WRITE (100, "('*** QDD VERSION AND STATUS AT COMPILATION TIME ***')")
    WRITE (100, "('********************************************************************************')")
    WRITE (100, "('Built on : ', T22, A23)") COMPILE_TIME
    WRITE (100, "('GIT commit hash : ', T22, A40)") GIT_HASH
    WRITE (100, "('GIT branch : ', T22, A20)") GIT_BRANCH
    WRITE (100, "('Working tree state : ', T22, A5)") GIT_STATUS
    WRITE (100, *)
    IF (nFilesChanged /= 0) THEN
      WRITE (100, "('Changed files')")
      WRITE (100, "('-------------------------------------------')")
      DO i = 1, nFilesChanged, 2
        IF (files(i) == 'M') files(i) = '(modified)'
        IF (files(i) == 'D') files(i) = '(deleted)'
        IF (files(i) == 'A') files(i) = '(added)'
        IF (files(i) == 'C') files(i) = '(copied)'
        IF (files(i) == 'R') files(i) = '(renamed)'
        IF (files(i) == 'U') files(i) = '(updated but unmerged)'
        IF (files(i) == 'T') files(i) = '(typechange)'
        DO j = 1, 256
          IF (files(i + 1) (j:j) == '>') files(i + 1) (j:j) = '/' ! String-array substring replacement
        END DO
        WRITE (100, "(A32,1X,A10)") files(i + 1), files(i)
      END DO
      WRITE (100, "('-------------------------------------------')")
      WRITE (100, *)
    END IF
    CLOSE (100)

    DEALLOCATE (files)

#ifdef __linux
    CALL writeQDDLibInfo()
#endif /* __linux */

  END SUBROUTINE

  SUBROUTINE writeQDDLibInfo()
    CHARACTER(256) :: command
    CALL getQDDPath()
    OPEN (UNIT=100, FILE=trim(filename), STATUS='unknown', ACCESS='append')
    WRITE (100, "('QDD executable location, linked library versions and locations : ', A64)") path
    CLOSE (100)
    command = 'ldd '//trim(path)//' >> '//trim(filename)
    CALL execute_command_line(command)
    OPEN (UNIT=100, FILE=trim(filename), STATUS='unknown', ACCESS='append')
    WRITE (100, "('********************************************************************************')")
    CLOSE (100)
    WRITE (*, "('********************************************************************************')")
    WRITE (*, *)
  END SUBROUTINE

  SUBROUTINE getQDDPath()
    pid = GETPID()

    WRITE (path, '(i0)') pid
    path = '/proc/'//TRIM(path)//'/exe'
    szret = readlink(TRIM(path)//C_NULL_CHAR, cbuf, SIZE(cbuf, KIND=C_SIZE_T))
    IF (szret == -1) STOP 'Error reading link'

    path = ''
    DO i = 1, SIZE(cbuf)
      IF (cbuf(i) == C_NULL_CHAR) EXIT
      path(i:i) = cbuf(i)
    END DO
    path = TRIM(path)
  END SUBROUTINE
END MODULE
