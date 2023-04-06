
!-----V_Ar_Ar-----------------------------------------------------

FUNCTION v_ar_ar(r)

  USE params
!USE kinetic
  IMPLICIT REAL(DP) (A - H, O - Z)

! Ar-Ar potential

  REAL, INTENT(IN OUT) :: r
  DATA alpha_ar, beta_ar, c6_ar, c8_ar, a_ar/1.7301, 1.7966, 55.465, 3672.9, 794.21/

!----------------------------------------------------------------------

  rabs = ABS(r)
  core = a_ar*EXP(-alpha_ar*rabs)/rabs - 2.0/(1.0 + EXP((beta_ar/rabs)**2)) &
         *(c6_ar/rabs**6 + c8_ar/rabs**8)

  v_ar_ar = core + ch(18)**2*v_soft(r, 2.0*sgm1(18))

  RETURN
END FUNCTION v_ar_ar
!-----V_Ar_Na-----------------------------------------------------

FUNCTION v_ar_na(r)

  USE params
!USE kinetic
  IMPLICIT REAL(DP) (A - H, O - Z)

! Ar-Na potential (repulsive core yet TO be found)

!----------------------------------------------------------------------

  rabs = ABS(r)

  v_ar_na = ch(18)*ch(11)*v_soft(r, sq2*sgm1(18))

  RETURN
END FUNCTION v_ar_na
!-----V_Ar_el-------------------------------------------------------

FUNCTION v_ar_el(r)

  USE params
!USE kinetic
  IMPLICIT REAL(DP) (A - H, O - Z)

  REAL, INTENT(IN) :: r
!DATA c_dipmod/3.378/ ! spring constant of dipole model IN Ha
! DATA sigma_Ar/1.43/
  DATA rstep/1.0E-3/
  DATA vdip/11.08/ ! effective polar.pot. IN Ha
  DATA vrep/52.8/ ! repulsive core IN Ha
  DATA rcut/3.1/ ! cut radius for pol.pot.
  DATA rcor/1.35/ ! inverse radius**2 for core
! DIMENSION ch(18:18)
! DATA ch(18)/6.119/
! DATA e2/2.0/
  DATA espil/1.0E-6/

!-------------------------------------------------------------------

  sigma_ar = sgm1(18)*sq2
  effch = ch(18)*e2
  effc = c_dipmod*e2 ! IN Ry
  rpmin = -effch*(v_soft(r + rstep, sigma_ar) - v_soft(r - rstep, sigma_ar)) &
          /(2.0*rstep*effc)
! WRITE(6,'(1x,i4,2g14.5)') 0,rpmin,abs(rpold-rpmin)
  DO iter = 1, 100
    rpold = rpmin
    rpmin = -effch*(v_soft(r + rstep + rpmin, sigma_ar) &
                    - v_soft(r + rpmin - rstep, sigma_ar))/(2.0*rstep*effc)
! WRITE(6,'(1x,i4,2g14.5)') iter,rpmin,abs(rpold-rpmin)
    IF (ABS(rpold - rpmin) < epsil) GO TO 19
  END DO
19 CONTINUE
  v_dipa = effch*(v_soft(r + rpmin, sigma_ar) - v_soft(r, sigma_ar)) &
           + 0.5*effc*rpmin*rpmin
  v_ar_el = vrep*e2*EXP(-rcor*r*r) - vdip*e2/(1.0 + EXP((3.1/r)**8))/r**4 &
            - v_dipa
!test WRITE(6,'(1x,f8.2,3g13.5)') r,V_dipa,V_Ar_el,V_Ar_el+V_dipa

  RETURN
END FUNCTION v_ar_el
