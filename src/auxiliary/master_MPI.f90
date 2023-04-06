PROGRAM master_MPI

  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)
  INTEGER, PARAMETER :: knode = 2

  WRITE (*, *) ' before mpi_init'
  CALL mpi_init(icode)

  WRITE (*, *) ' before mpi_comm_size'
  CALL mpi_comm_size(mpi_comm_world, nprocs, icode)

  WRITE (*, *) ' before mpi_comm_rank'
  CALL mpi_comm_rank(mpi_comm_world, myn, icode)
  WRITE (*, *) ' nprocs,myn=', nprocs, myn

  IF (nprocs < knode) STOP " not enough nodes"

  IF (myn == 0) THEN
    WRITE (*, *) ' submitted myn=', myn
    CALL system('./testsub.bat 001')
  ELSE IF (myn == 1) THEN
    WRITE (*, *) ' submitted myn=', myn
    CALL system('./testsub.bat 010')
  ELSE IF (myn == 2) THEN
    WRITE (*, *) ' submitted myn=', myn
    CALL system('./testsub.bat 011')
  ELSE IF (myn == 3) THEN
    WRITE (*, *) ' submitted myn=', myn
    CALL system('./testsub.bat 100')
  ELSE IF (myn == 4) THEN
    WRITE (*, *) ' submitted myn=', myn
    CALL system('./testsub.bat 101')
  ELSE IF (myn == 5) THEN
    WRITE (*, *) ' submitted myn=', myn
    CALL system('./testsub.bat 110')
  END IF

  CALL mpi_barrier(mpi_comm_world, mpi_ierror)
  CALL mpi_finalize(icode)

END PROGRAM master_MPI
