
#if(raregas)
! This FILE CONTAINS tabulated versions of some COMMON
! functions

!------------------------------------------------------------
!
SUBROUTINE initfunctions
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
  REAL(DP) :: add0, add1, add1t
  REAL(DP) :: coreheight, coreradius, corewidth
  REAL(DP) :: r
  REAL(DP) :: sigvdw, vdwh
  INTEGER :: ind

  REAL(DP), EXTERNAL :: v_soft, dvsdr

  WRITE (6, *) 'INITIALISING TABULATED FUNCTIONS...'

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Initialisation of the Fermi FUNCTION
! (argon core-electron potential)

  IF (ivdw == 0 .OR. ivdw == 1) THEN

    coreheight = 0.47D0
    coreradius = 2.2D0
    corewidth = 1.6941D0

  ELSE IF (ivdw == 2) THEN

    coreheight = 0.14D0
    coreradius = 2.2D0
    corewidth = 1.515502D0

    sigvdw = 6.0D0
    vdwh = 201.25D0

  END IF
#if(raregas)
  IF (isurf /= 0) THEN
    IF (ivdw == 0 .OR. ivdw == 1) THEN
      DO ind = 1, kfermi
        r = ind*dx/1733D0
        varelcore0(ind) = -e2*coreheight/(1D0 + EXP(corewidth*(r - coreradius)))
        varelcore1(ind) = &
          e2*coreheight*corewidth*EXP(corewidth*(r - coreradius))/ &
          (1D0 + EXP(corewidth*(r - coreradius)))**2
      END DO
    ELSE IF (ivdw == 2) THEN
! WRITE(7,*) 'vArElCore0:'
      DO ind = 1, kfermi
        r = ind*dx/1733D0

        varelcore0(ind) = -e2*coreheight/(1.+EXP(corewidth*(r - coreradius)))

        varelcore1(ind) = &
          e2*coreheight*corewidth*EXP(corewidth*(r - coreradius))/ &
          (1D0 + EXP(corewidth*(r - coreradius)))**2

        add0 = v_soft(r, sigvdw)
        add0 = add0**8*e2*r*r*vdwh

        add1t = v_soft(r, sigvdw)
        add1 = -(-add1t**8*r*2 - r*r*8*add1t**7*dvsdr(r, sigvdw))
        add1 = add1*e2*vdwh

        varelcore0(ind) = varelcore0(ind) + add0
        varelcore1(ind) = varelcore1(ind) + add1
! WRITE(7,'(2(1pg13.5))') r,varelcore0(ind)
      END DO
    END IF
! determine cutoff for V_elArcore
    DO ind = 1, kfermi
      r = ind*dx/1733D0
      IF (ABS(varelcore0(ind)) < varelcorelimit) THEN
        nsg_arelcore = INT(r/MIN(dx, dy, dz)) + 1
        EXIT
      END IF
    END DO
    WRITE (6, '(a,i5)') ' nsg_arelcore=', nsg_arelcore
    WRITE (7, '(a,i5)') ' nsg_arelcore=', nsg_arelcore
  END IF
#endif
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! general Fermi FUNCTION
  DO ind = 1, kfermi
    r = ind*dx/1733D0

    vfermi0(ind) = -e2*fermia/(1D0 + EXP(fermib*(r - fermic)))

    vfermi1(ind) = e2*fermia*fermib*EXP(fermib*(r - fermic))/ &
                   (1D0 + EXP(fermib*(r - fermic)))**2
  END DO

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! MORE TABULATED FUNCTIONS CAN BE IMPLEMENTED HERE:

!old DO ind=1,kfermi
!old r=ind*dx/1733.
!old
!old V_soft_tab0(ind,1)=V_soft(r,sigmac*SQ2)
!old V_soft_tab0(ind,2)=V_soft(r,sigmav*SQ2)
!old V_soft_tab0(ind,3)=V_soft(r,sigmak*SQ2)
!old
!old V_soft_tab1(ind,1)=
!old & 2.*exp(-r**2/2./sigmac**2)/sqrt(2.*PI)/sigmac/r
!old & -V_soft(r,sigmac*SQ2)/r
!old V_soft_tab1(ind,2)=
!old & 2.*exp(-r**2/2./sigmav**2)/sqrt(2.*PI)/sigmav/r
!old & -V_soft(r,sigmav*SQ2)/r
!old V_soft_tab1(ind,3)=
!old & 2.*exp(-r**2/2./sigmak**2)/sqrt(2.*PI)/sigmak/r
!old & -V_soft(r,sigmak*SQ2)/r
!old
!old V_soft_tab2(ind,1)=2./(sqrt(PI)*SQ2**3*sigmac**3*r**3)*
!old & (-2.*exp(-r**2/2./sigmac**2)*r**3-2.*exp(-r**2/2./sigmac**2)
!old & +V_soft(r,sigmac*SQ2)*sqrt(PI)*sigmac**3*SQ2**3)
!old
!old V_soft_tab2(ind,2)=2./(sqrt(PI)*SQ2**3*sigmav**3*r**3)*
!old & (-2.*exp(-r**2/2./sigmav**2)*r**3-2.*exp(-r**2/2./sigmav**2)
!old & +V_soft(r,sigmav*SQ2)*sqrt(PI)*sigmav**3*SQ2**3)
!old
!old V_soft_tab2(ind,3)=2./(sqrt(PI)*SQ2**3*sigmak**3*r**3)*
!old & (-2.*exp(-r**2/2./sigmak**2)*r**3-2.*exp(-r**2/2./sigmak**2)
!old & +V_soft(r,sigmak*SQ2)*sqrt(PI)*sigmak**3*SQ2**3)
!old
!old END DO

!tab DO ind=1,kfermi
!tab r=ind*dx/1733.
!tab
!tab V_soft_tabc0(ind)=V_soft(r,sigmac*SQ2)
!tab V_soft_tabe0(ind)=V_soft(r,sigmav*SQ2)
!tab V_soft_tabk0(ind)=V_soft(r,sigmak*SQ2)
!tab
!tab V_soft_tabc1(ind)=
!tab & 2.*exp(-r**2/2./sigmac**2)/sqrt(2.*PI)/sigmac/r
!tab & -V_soft(r,sigmac*SQ2)/r
!tab V_soft_tabe1(ind)=
!tab & 2.*exp(-r**2/2./sigmav**2)/sqrt(2.*PI)/sigmav/r
!tab & -V_soft(r,sigmav*SQ2)/r
!tab V_soft_tabk1(ind)=
!tab & 2.*exp(-r**2/2./sigmak**2)/sqrt(2.*PI)/sigmak/r
!tab & -V_soft(r,sigmak*SQ2)/r
!tab
!tab END DO

  RETURN
END SUBROUTINE initfunctions
!------------------------------------------------------------

!------------------------------------------------------------

!old c------------------------------------------------------------
!old FUNCTION V_soft_tab(r,ip)
!old c------------------------------------------------------------
!old USE params
!old
!old deltar=dx/1733
!old ir = nint(r/deltar) + 1
!old rn = ir*deltar
!old drn = r-rn
!old
!old
!old IF (ir.gt.kfermi) THEN
!old V_soft_tab=0.0
!old ELSE
!old V_soft_tab= V_soft_tab0(ir,ip)+drn*V_soft_tab1(ir,ip)
!old END IF
!old
!old
!old RETURN
!old END
!------------------------------------------------------------

!tabc------------------------------------------------------------
!tab FUNCTION V_soft_tabc(r)
!tabc------------------------------------------------------------
!tabUSE params
!tab
!tab deltar=dx/1733
!tab ir = nint(r/deltar) + 1
!tab rn = ir*deltar
!tab drn = r-rn
!tab
!tab
!tab IF (ir.gt.kfermi) THEN
!tab V_soft_tabc=0.0
!tab elseif(r.lt.8.0*sigmac*SQ2) THEN
!tab V_soft_tabc= V_soft_tabc0(ir)+drn*V_soft_tabc1(ir)
!tab ELSE
!tabc V_soft_tabc= V_soft_tabc0(ir)
!tab V_soft_tabc=1/r
!tab END IF
!tab
!tab
!tab RETURN
!tab END
!tabc------------------------------------------------------------
!tab
!tab
!tabc------------------------------------------------------------
!tab FUNCTION V_soft_tabe(r)
!tabc------------------------------------------------------------
!tabUSE params
!tab
!tab deltar=dx/1733
!tab ir = nint(r/deltar) + 1
!tab rn = ir*deltar
!tab drn = r-rn
!tab
!tab
!tab IF (ir.gt.kfermi) THEN
!tab V_soft_tabe=0.0
!tab elseif (r.lt.8.0*sigmav*SQ2) THEN
!tab V_soft_tabe= V_soft_tabe0(ir)+drn*V_soft_tabe1(ir)
!tab ELSE
!tabc V_soft_tabe= V_soft_tabe0(ir)
!tab V_soft_tabe=1/r
!tab END IF
!tab
!tab
!tab RETURN
!tab END
!tabc------------------------------------------------------------
!tab
!tabc------------------------------------------------------------
!tab FUNCTION V_soft_tabk(r)
!tabc------------------------------------------------------------
!tabUSE params
!tab
!tab deltar=dx/1733
!tab ir = nint(r/deltar) + 1
!tab rn = ir*deltar
!tab drn = r-rn
!tab
!tab
!tab IF (ir.gt.kfermi) THEN
!tab V_soft_tabk=0.0
!tab elseif (r.lt.8.0*sigmak*SQ2) THEN
!tab V_soft_tabk= V_soft_tabk0(ir)+drn*V_soft_tabk1(ir)
!tab ELSE
!tabc V_soft_tabk= V_soft_tabk0(ir)
!tab V_soft_tabk= 1/r
!tab END IF
!tab
!tab
!tab RETURN
!tab END
!tabc------------------------------------------------------------

!------------------------------------------------------------
!
REAL(DP) FUNCTION varelcore(r)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN)::r
  REAL(DP):: deltar
  REAL(DP):: rn, drn
  INTEGER:: ir
! linear interpolation of Fermi FUNCTION

! deltar = dx/1733*kxbox
  deltar = dx/1733D0

  ir = nint(r/deltar) + 1
  rn = ir*deltar
  drn = r - rn

  varelcore = varelcore0(ir) + drn*varelcore1(ir)

  RETURN
END FUNCTION varelcore
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION vfermi(r)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN)::r
  REAL(DP):: deltar
  REAL(DP):: rn, drn
  INTEGER:: ir
! linear interpolation of general Fermi FUNCTION

  deltar = dx/1733D0

  ir = nint(r/deltar) + 1
  rn = ir*deltar
  drn = r - rn

  vfermi = vfermi0(ir) + drn*vfermi1(ir)

  RETURN
END FUNCTION vfermi
!------------------------------------------------------------
#endif

