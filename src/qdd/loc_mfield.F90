
!-----coul_mfield------------------------------------------------

SUBROUTINE coul_mfield(rho)

! The Coulomb part of the mean field.

! Input:
! rho = electron density
! Output:
! Coulomb field via 'chpcoul' in module 'params'

  USE params
  USE coulsolv, ONLY: solv_poisson_f, solv_poisson_e, tcoulfalr

  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)

#if(raregas)
  INTEGER :: ii
  REAL(DP), ALLOCATABLE :: rhotmp(:)
#endif

!----------------------------------------------------------------

! copy density on intermediate storage
! add images for dielectric contribution
#if(raregas)
  IF (idielec == 1) THEN

    ALLOCATE (rhotmp(2*kdfull2))
    DO ii = 1, 2*kdfull2
      rhotmp(ii) = rho(ii)
    END DO

    CALL addimage(rho, 1)

  END IF
#endif

! Coulomb of the electronic density

  IF (tcoulfalr) THEN
    CALL solv_poisson_f(rho, chpcoul, kdfull2)
  ELSE
    CALL solv_poisson_e(rho, chpcoul, kdfull2)
  END IF

! computation of the Coulomb potential from the electronic density
! before adjustdip or vstep for pseudodensity description
! restore electron density, SAVE Coulomb fields,
! add images and solve Coulomb next round
#if(raregas)
  IF (idielec == 1) THEN

    DO ii = 1, 2*kdfull2
      rho(ii) = rhotmp(ii)
    END DO
    DO ii = 1, kdfull2
      rfieldtmp(ii) = chpcoul(ii)
    END DO

    CALL addimage(rho, 0)

    IF (tcoulfalr) THEN
      CALL solv_poisson_f(rho, chpcoul, kdfull2)
    ELSE
      CALL solv_poisson_e(rho, chpcoul, kdfull2)
    END IF

    DO ii = 1, kdfull2
      CALL conv1to3(ii)
      IF (iindtmp(1) >= nint(xdielec/dx) + nxsh) THEN
        chpcoul(ii) = rfieldtmp(ii)
      END IF
    END DO
    DO ii = 1, 2*kdfull2
      rho(ii) = rhotmp(ii)
    END DO
    DEALLOCATE (rhotmp)

  END IF
#endif

  RETURN
END SUBROUTINE coul_mfield

!-----ldapsp_mfield------------------------------------------------

SUBROUTINE ldapsp_mfield(rho, aloc)

! The LDA and pseudopotential part of the mean field.
! Input:
! rho = electron density
! ionic positions for PsP come in via COMMON
! Output:
! aloc = local mean field

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)
  REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)

!------------------------------------------------------------------

  CALL calcpseudo() ! update pseudo-potentials

  CALL calclocal(rho, aloc) ! LDA part of the potential

  RETURN
END SUBROUTINE ldapsp_mfield

