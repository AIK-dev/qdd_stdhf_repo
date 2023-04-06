!This file is a part of PW-TELEMAN project.
!PW-TELEMAN is a Time-Dependent Electronic Dynamics in Molecules And Nanosystems library.
!Copyright (C) 2011-2015  Paul-Gerhard Reinhard, Eric Suraud, Florent Calvayrac,
!Phuong Mai Dinh, David Brusson, Philipp Wopperer, José María Escartín Esteban.
!
!PW-Teleman is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!PW-Teleman is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with PW-Teleman.  If not, see <http://www.gnu.org/licenses/>.

 
!-----zeroforce-----------------------------------------------------------------

SUBROUTINE zeroforce(aloc,rho)

!     Corrects SIC-Slater and SIC-KLI potentials to meet the
!     zero-force condition.
!     
!     Input:
!       rho   = local density
!     Input/Output:
!       aloc  = local KS potential before and after correction

USE params
USE util, ONLY:wfovlp
USE kinetic
IMPLICIT NONE


REAL(DP), INTENT(IN OUT)   :: aloc(2*kdfull2)
REAL(DP), INTENT(IN)       :: rho(2*kdfull2)

INTEGER :: i, is, ishift
REAL(DP) :: denominator, counter, xlambda, ylambda, zlambda
COMPLEX(DP), ALLOCATABLE :: potk(:),dervk(:)     ! for Fourier transformed potentials
REAL(DP), ALLOCATABLE :: potwork(:)

!-----------------------------------------------------------------------------

IF(.NOT.ALLOCATED(akv)) STOP "ZEROFORCE requires FFT"

ALLOCATE(potk(kdfull2),dervk(kdfull2),potwork(kdfull2))

!     loop over spin-up and spin-down

DO is=1,numspin
    ishift = (is-1)*nxyz
    
!       Fourier transformation
    CALL rftf(aloc(1+ishift),potk)
    
!       Laplacian and integral
    
    DO  i=1,nxyz
      dervk(i) = potk(i)*akv(i)
    END DO
    CALL rfftback(dervk,potwork)
    denominator = wfovlp(rho(1+ishift : kdfull2+ishift),potwork)
    
!   x-derivative
    dervk = -akx*potk
    CALL rfftback(dervk,potwork)
    counter = wfovlp(rho(1+ishift : kdfull2+ishift),potwork)
    xlambda = counter/denominator
    DO i=1,nxyz
      aloc(i+ishift) = aloc(i+ishift)-xlambda*potwork(i)
    END DO
    
!   y-derivative
    dervk = -aky*potk
    CALL rfftback(dervk,potwork)
    counter = wfovlp(rho(1+ishift : kdfull2+ishift),potwork)
    ylambda = counter/denominator
    DO i=1,nxyz
      aloc(i+ishift) = aloc(i+ishift)-ylambda*potwork(i)
    END DO
    
!   z-derivative
    dervk = -akz*potk
    CALL rfftback(dervk,potwork)
    counter = wfovlp(rho(1+ishift : kdfull2+ishift),potwork)
    zlambda = counter/denominator
    DO i=1,nxyz
      aloc(i+ishift) = aloc(i+ishift)-zlambda*potwork(i)
    END DO
    
  END DO
  DEALLOCATE(potk,dervk,potwork)
  
  RETURN
END SUBROUTINE zeroforce

!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE checkzeroforce(rho,aloc)
!------------------------------------------------------------

!     Checks whether SIC-Slater or SIC-KLI potentials fulfill the
!     zero-force condition.
!     
!     Input:
!       rho   = local density
!       aloc  = local KS potential 

USE params
USE util, ONLY:wfovlp
USE kinetic
IMPLICIT NONE

REAL(DP), INTENT(IN)  :: rho(2*kdfull2)
REAL(DP), INTENT(IN)  :: aloc(2*kdfull2)

INTEGER :: is, ishift
REAL(DP) :: zforcex, zforcey, zforcez
REAL(DP) :: test1,test2,test3

!       workspaces

COMPLEX(DP), ALLOCATABLE :: potk(:),dervk(:)     ! for Fourier transformed potentials
REAL(DP), ALLOCATABLE :: potwork(:)

IF(.NOT.ALLOCATED(akv)) STOP "ZEROFORCE requires FFT"

ALLOCATE(potk(kdfull2),dervk(kdfull2),potwork(kdfull2))

DO is=1,1 !numspin
    ishift = (is-1)*nxyz
    
!       Fourier transformation
    CALL rftf(aloc(1+ishift),potk)
    
!   x-derivative
    dervk = -akx*potk
    CALL rfftback(dervk,potwork)
    zforcex = wfovlp(rho(1+ishift : kdfull2+ishift),potwork)
    
!   y-derivative
    dervk = -aky*potk
    CALL rfftback(dervk,potwork)
    zforcey = wfovlp(rho(1+ishift : kdfull2+ishift),potwork)
    
!   z-derivative
    dervk = -akz*potk
    CALL rfftback(dervk,potwork)
    zforcez = wfovlp(rho(1+ishift : kdfull2+ishift),potwork)

!    IF(is==1) THEN
    test1 = wfovlp(rho(1+ishift : kdfull2+ishift),rho(1+ishift : kdfull2+ishift))
    test2 = wfovlp(aloc(1+ishift : kdfull2+ishift),aloc(1+ishift : kdfull2+ishift))
    test3 = SUM(rho(1+ishift : kdfull2+ishift))
!    END IF
  WRITE(500,'(1f10.3,6(1pg13.5))') tfs,zforcex,zforcey,zforcez, test1,test2,test3
    
  END DO
  
  DEALLOCATE(potk,dervk,potwork)
  WRITE(6,'(a,3e17.7)') 'Checking Zero-Force Theorem: ',  &
      zforcex,zforcey,zforcez
  
!  WRITE(500,'(1f17.5,3e17.7)') tfs,zforcex,zforcey,zforcez
!  WRITE(500,'(a,4e17.7)') 'tests:',test1,test2,test3,SUM(rho)
  
  RETURN
END SUBROUTINE checkzeroforce
