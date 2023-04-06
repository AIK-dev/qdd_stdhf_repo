!------------------------------------------------------------
!
 PROGRAM analyze_MP
!
! TO analyze 'pMP' FILE for PES
!
  implicit REAL*8(a-h,o-z)
  INTEGER, PARAMETER :: ktimes=165400            ! Maximum number of time steps read in
  INTEGER, PARAMETER :: kpoints=92                ! number of measuring points
  INTEGER, PARAMETER :: kstate=8                 ! number of s.p. states
  CHARACTER*(*), PARAMETER :: name = 'H2O'       ! set up the name
  LOGICAL, PARAMETER :: tphasecorr = .FALSE.      ! switch to include phase correction
  INTEGER, PARAMETER :: komega=2000              ! maximum number of frequencies
  REAL(8), PARAMETER    :: delomega=0.005        ! step in frequency alias kinetic energies, in Ry
  CHARACTER*1000 command,cin
  CHARACTER*6 asci6
  CHARACTER*8 asci8

  complex(8) :: q0(ktimes,kpoints,kstate)
  complex(8) :: accum(kpoints,kstate)
  REAL(8)    :: str(kpoints)
  REAL(8)    :: fieldcomponent(2,kpoints)        ! projection of E-field on outgoing direction
  REAL(8)    :: fpulseintegs(4,ktimes)
  complex(8) :: cfac
  REAL(8)    :: time,time2,omega,timfin,deltime,dummy(3),pact
  INTEGER    :: iomega,i,j,nbe,ntimes,itimes

  WRITE(cin,*) kpoints+2
  cin = adjustl(cin)
  command = 'head --lines='// trim(cin) // ' pMP.'// name // ' > pPES.' // name
  call system(command)

  open(01,file='pMP.'//name)
  open(02,position='append',file='pPES.'//name)

  IF(tphasecorr) THEN
    read(1,*)
    read(1,*)
    DO j=1,kpoints
      read(1,'(a,i6,a,3f12.1,25x,2e14.5)',err=98) asci8,idum,asci6,dummy,fieldcomponent(:,j)
      WRITE(*,*) 'init:',j,fieldcomponent(:,j)
    END DO
  ELSE
    DO j=1,kpoints+2
      read(01,*,err=98) 
    END DO
  END IF

  inline = kpoints + 2
  timfin = -1.0

  DO itimes=0,ktimes
    inline = inline+1
    IF(tphasecorr) THEN
!      WRITE(*,*) 'itimes:',itimes
      read(01,'(1f14.5,1000e18.8)',err=98,end=99) time,q0(itimes,:,1),fpulseintegs(:,itimes) 
!      WRITE(*,*) 'time etc:',time,fpulseintegs(:,itimes) 
    ELSE
!      WRITE(*,*) 'itimesold:',itimes
      read(01,'(1f14.5,1000e18.8)',err=98,end=99) time,q0(itimes,:,1)
    END IF
    if(time.lt.timfin) goto 98
    ntimes = itimes
    timfin = time
 
    WRITE(*,*) 'line= ',inline,'read time= ',time

    do j=2,kstate
      inline=inline+1
!      IF(tphasecorr) THEN
!        read(01,'(1f14.5,1000e18.8)',err=98,end=98) time2,q0(itimes,:,j),fpulseintegs(:,itimes)
!      ELSE
        read(01,'(1f14.5,1000e18.8)',err=98,end=98) time2,q0(itimes,:,j)
!      END IF
      if(time2.ne.time) goto 98
    END DO
  END DO
  stop ' file count exhausted'

 98   continue
      WRITE(*,*) 'error in line= ',inline,' time= ',time
      stop ' inconsistent input file'
    
 99   continue

  close(01)
  timfin = timfin/0.0484D0                                                 !  time in 1/Ry
  deltim = timfin/ntimes
  DO iomega=1,komega
    omega = iomega*delomega                                                !  omega in Ry
    IF(tphasecorr) pact = SQRT(omega)

    accum = dcmplx(0.0D0,0.0D0)
    DO itimes=0,ntimes
      IF(tphasecorr) THEN
        cfac   = cdexp(dcmplx(0.0D0,omega*itimes*deltim))        
        DO j=1,kpoints
          accum(j,:) = accum(j,:) + cfac*q0(itimes,j,:)* &
                    EXP(CMPLX(0D0, &
                    +2D0*pact*fpulseintegs(1,itimes)*fieldcomponent(1,j) &
                    +2D0*pact*fpulseintegs(3,itimes)*fieldcomponent(2,j) &
                     +fpulseintegs(2,itimes)+fpulseintegs(4,itimes),8))
        END DO
      ELSE
        cfac   = cdexp(dcmplx(0.0D0,omega*itimes*deltim))          
        accum = accum + cfac*q0(itimes,:,:)
      END IF
    END DO

    str = 0.0
    DO j=1,kstate
      str = str + accum(:,j)*conjg(accum(:,j))
    END DO    

    WRITE(02,'(f14.5,1000(1pg18.8))') omega,str                            !  cPW
!    WRITE(6,'(f8.3,1000(1pg12.5))') omega,str


  END DO

  STOP
  END

