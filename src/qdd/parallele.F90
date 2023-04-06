
!-----init_parallele-------------------------------------------------

SUBROUTINE init_parallele()

! Initializes parallele computing and determines actual node
! NUMBER 'myn'. The actual 'myn' is communicated via module 'params'.

  USE params
  IMPLICIT NONE
#if(mpi)
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)
  INTEGER :: nb
#endif
  INTEGER :: nprocs

!------------------------------------------------------------------

  myn = 0

#if(mpi)
  WRITE (*, *) ' before mpi_init'
  CALL mpi_init(mpi_ierror)
  WRITE (*, *) ' before mpi_comm_size'
  CALL mpi_comm_size(mpi_comm_world, nprocs, mpi_ierror)
  WRITE (*, *) nprocs, knode
  knode = nprocs
  WRITE (*, *) ' before mpi_comm_rank'
  CALL mpi_comm_rank(mpi_comm_world, myn, mpi_ierror)
  nb = myn
#endif

! Set the number of OpenMP threads to use equal to the number of CPUs or,
! IF the OMP_NUM_THREADS env. var. is set, use that value.
! Also, the dynamic adjustment of the number of threads
! available for execution of parallel regions is enabled to prevent
! creating more threads than there are CPU's
#if(omp)
  IF (numthr > 0) setdyn = .FALSE.
!  CALL OMP_SET_DYNAMIC(setdyn)
  numthr = OMP_GET_MAX_THREADS()
!  CALL OMP_SET_NUM_THREADS(numthr) ! set the system wide number of OMP threads
  nthr = numthr - 1
#if(omp_debug)
  WRITE (*, *) 'setdyn = ', setdyn
  WRITE (*, *) 'numthr = ', numthr
  WRITE (*, *) 'OMP_GET_MAX_THREADS() = ', OMP_GET_MAX_THREADS()
  WRITE (*, *) 'OMP_GET_NUM_PROCS() = ', OMP_GET_NUM_PROCS()
  WRITE (*, *) 'OMP_GET_DYNAMIC() = ', OMP_GET_DYNAMIC()
  WRITE (*, *) 'OMP_GET_NESTED() = ', OMP_GET_NESTED()
!$OMP PARALLEL
  WRITE (*, *) 'OMP_GET_NUM_THREADS() = ', OMP_GET_NUM_THREADS()
  WRITE (*, *) 'OMP_GET_THREAD_NUM() = ', OMP_GET_THREAD_NUM()
!$OMP END PARALLEL
#endif
  WRITE (*, *) 'INIT_PARALLELE(): numthr=', numthr
  WRITE (*, *) 'INIT_PARALLELE(): nthr=', nthr
#else
  nthr = 0
#endif

  RETURN
END SUBROUTINE init_parallele

#if(mpi)
!-----comm_inputparams-------------------------------------------------

SUBROUTINE comm_inputparams()

! Communicates input parameters from node "0" to other nodes.

  USE params
  USE kinetic
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER :: is(mpi_status_size)
  INTEGER :: ii, ic, nod
#if(raregas)
  INTEGER :: i, j
#endif
!------------------------------------------------------------------

  CALL mpi_comm_rank(mpi_comm_world, myn, mpi_ierror)

  CALL i_scatter(kxbox)
  CALL i_scatter(kybox)
  CALL i_scatter(kzbox)
  CALL i_scatter(kstate)

  CALL pi_scatter(scaleionx)
  CALL pi_scatter(scaleiony)
  CALL pi_scatter(scaleionz)
  CALL pi_scatter(shiftionx)
  CALL pi_scatter(shiftiony)
  CALL pi_scatter(shiftionz)
  CALL pi_scatter(rotionx)
  CALL pi_scatter(rotiony)
  CALL pi_scatter(rotionz)
  CALL pi_scatter(shiftwfx)
  CALL pi_scatter(shiftwfy)
  CALL pi_scatter(shiftwfz)

  CALL i_scatter(ispinsep)
  CALL i_scatter(tshiftcmtoorigin)
  CALL i_scatter(tplotorbitals)
  CALL i_scatter(ievaluate)

  CALL i_scatter(tmoms_rel_cm)
  CALL i_scatter(ifspemoms)
  CALL i_scatter(iffastpropag)
  CALL i_scatter(ifexpevol)
#if(extended)
!  CALL i_scatter(ifcnevol)
#endif
  CALL i_scatter(iflocaliz)
  CALL i_scatter(jgeomion)
  CALL i_scatter(jangabso)
  CALL i_scatter(jgeomel)
  CALL i_scatter(jelf)
  CALL i_scatter(jstinf)
  CALL i_scatter(iangabso)
  CALL i_scatter(nangtheta)
  CALL i_scatter(nangphi)
  CALL i_scatter(tfreezekspot)
  CALL i_scatter(ispherabso)
  CALL i_scatter(tfixcmion)
  CALL i_scatter(nmptheta)
  CALL i_scatter(nmpphi)
  CALL i_scatter(jmp)
  CALL i_scatter(jnorms)
  CALL i_scatter(jdensitydiff)
  CALL i_scatter(jdensitydiff2d)
  CALL i_scatter(jdensity2d)
  CALL i_scatter(jcharges)

  CALL pi_scatter(drcharges)
  CALL pi_scatter(tempion)
  CALL pi_scatter(projcharge)
  CALL pi_scatter(delomega)
  CALL pi_scatter(angthetal)
  CALL pi_scatter(angthetah)
  CALL pi_scatter(angphil)
  CALL pi_scatter(angphih)
  CALL pi_scatter(powabso)

#if(extended)
  CALL pi_scatter(dinmargin)
  CALL i_scatter(izforcecorr)
  CALL i_scatter(ntref)
  CALL i_scatter(jescmask)
  CALL i_scatter(jescmaskorb)
  CALL i_scatter(iscatterelectron)
  CALL pi_scatter(trequest)
  CALL pi_scatter(timefrac)
  CALL pi_scatter(scatterelectronenergy)
  CALL pi_scatter(scatterelectronvxn)
  CALL pi_scatter(scatterelectronvyn)
  CALL pi_scatter(scatterelectronvzn)
  CALL pi_scatter(scatterelectronx)
  CALL pi_scatter(scatterelectrony)
  CALL pi_scatter(scatterelectronz)
  CALL pi_scatter(scatterelectronw)
#endif

#if(raregas)
  CALL pi_scatter(ekin0pp)
  CALL pi_scatter(vxn0)
  CALL pi_scatter(vyn0)
  CALL pi_scatter(vzn0)
#endif

  IF (myn == 0 .AND. knode /= 1) THEN
  DO nod = 1, knode - 1

    CALL mpi_send(outnam, 13, mpi_character, nod, 1, mpi_comm_world, ic)

    CALL mpi_send(nion2, 1, mpi_integer, nod, 1, mpi_comm_world, ic)

    CALL mpi_send(nabsorb, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(ifsicp, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(ionmdtyp, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(ifredmas, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(icooltyp, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(ipsptyp, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(numspin, 1, mpi_integer, nod, 1, mpi_comm_world, ic)

    CALL mpi_send(ismax, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(itmax, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(isaves, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(isaved, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(istinf, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(jstinf, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(dt1, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(irest, 1, mpi_integer, nod, 1, mpi_comm_world, ic)

    CALL mpi_send(centfx, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(centfy, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(centfz, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)

    CALL mpi_send(tspindip, 1, mpi_logical, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(irotat, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(phirot, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)

    CALL mpi_send(jpos, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(jposcm, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(jvel, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(jesc, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(jforce, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(tstat, 1, mpi_logical, nod, 1, mpi_comm_world, ic)

    CALL mpi_send(jlaser, 1, mpi_integer, nod, 1, mpi_comm_world, ic)

    CALL mpi_send(jdip, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(jdiporb, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(jquad, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(jang, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(jspdp, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(jinfo, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(jenergy, 1, mpi_integer, nod, 1, mpi_comm_world, ic)

    CALL mpi_send(itft, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(tnode, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(deltat, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(tpeak, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(omega, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(e0, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)

    CALL mpi_send(e1x, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(e1y, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(e1z, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(e2x, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(e2y, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(e2z, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)
    CALL mpi_send(phi, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)

    CALL mpi_send(nelect, 1, mpi_integer, nod, 2, mpi_comm_world, ic)

    CALL mpi_send(nion, 1, mpi_integer, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(scaleion, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)

    CALL mpi_send(charge, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(nspdw, 1, mpi_integer, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(temp, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(occmix, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(b2occ, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(gamocc, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(deocc, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(osfac, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(init_ao, 1, mpi_integer, nod, 2, mpi_comm_world, ic)

    CALL mpi_send(dx, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(dy, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(dz, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)

    CALL mpi_send(endcon, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(itback, 1, mpi_integer, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(epswf, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(e0dmp, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(epsoro, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(dpolx, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(dpoly, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
    CALL mpi_send(dpolz, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
#if(extended)
      CALL mpi_send(ipseudo, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
      CALL mpi_send(idenspl, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
      CALL mpi_send(i3dx, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
      CALL mpi_send(i3dz, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
      CALL mpi_send(i3dstate, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
      CALL mpi_send(eproj, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)
      CALL mpi_send(nproj, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
      CALL mpi_send(nproj_states, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
      WRITE (*, *) 'proj_states=', proj_states
      DO ii = 1, nproj_states
        CALL mpi_send(proj_states(ii), 1, mpi_integer, nod, 1, mpi_comm_world, ic)
      END DO
#endif
#if(raregas)
      CALL mpi_send(ivdw, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
      CALL mpi_send(ekmat, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)
      CALL mpi_send(nc, 1, mpi_integer, nod, 2, mpi_comm_world, ic)
      CALL mpi_send(NE, 1, mpi_integer, nod, 2, mpi_comm_world, ic)
      CALL mpi_send(nk, 1, mpi_integer, nod, 2, mpi_comm_world, ic)
      CALL mpi_send(surftemp, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
      CALL mpi_send(sigmac, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
      CALL mpi_send(sigmav, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
      CALL mpi_send(sigmak, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
      CALL mpi_send(chgc0, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
      CALL mpi_send(chge0, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
      CALL mpi_send(chgk0, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
      CALL mpi_send(cspr, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
      CALL mpi_send(mion, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
      CALL mpi_send(me, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
      CALL mpi_send(mkat, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)

      DO i = 1, 5
        DO j = 1, 5
          CALL mpi_send(isrtyp(i, j), 1, mpi_integer, nod, 2, mpi_comm_world, ic)
        END DO
      END DO

      CALL mpi_send(shiftx, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
      CALL mpi_send(shifty, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)
      CALL mpi_send(shiftz, 1, mpi_double_precision, nod, 2, mpi_comm_world, ic)

      CALL mpi_send(iaxis, 1, mpi_integer, nod, 2, mpi_comm_world, ic)
#endif
    END DO

  ELSE IF (myn /= 0 .AND. knode /= 1) THEN

    CALL mpi_recv(outnam, 13, mpi_character, 0, mpi_any_tag, mpi_comm_world, is, ic)

    CALL mpi_recv(nion2, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)

    CALL mpi_recv(nabsorb, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(ifsicp, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(ionmdtyp, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(ifredmas, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(icooltyp, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(ipsptyp, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(numspin, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)

    CALL mpi_recv(ismax, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(itmax, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(isaves, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(isaved, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(istinf, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(jstinf, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(dt1, 1, mpi_double_precision, 0, mpi_any_tag, &
                mpi_comm_world, is, ic)
    CALL mpi_recv(irest, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)

    CALL mpi_recv(centfx, 1, mpi_double_precision, 0, mpi_any_tag, &
                mpi_comm_world, is, ic)
    CALL mpi_recv(centfy, 1, mpi_double_precision, 0, mpi_any_tag, &
                mpi_comm_world, is, ic)
    CALL mpi_recv(centfz, 1, mpi_double_precision, 0, mpi_any_tag, &
                mpi_comm_world, is, ic)

    CALL mpi_recv(tspindip, 1, mpi_logical, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(irotat, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(phirot, 1, mpi_double_precision, 0, mpi_any_tag, &
                mpi_comm_world, is, ic)


    CALL mpi_recv(jpos, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(jposcm, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(jvel, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(jesc, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(jforce, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(tstat, 1, mpi_logical, 0, mpi_any_tag, mpi_comm_world, is, ic)

    CALL mpi_recv(jlaser, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)

    CALL mpi_recv(jdip, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(jdiporb, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(jquad, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(jang, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(jspdp, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(jinfo, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(jenergy, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)

    CALL mpi_recv(itft, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(tnode, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(deltat, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(tpeak, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(omega, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)

    CALL mpi_recv(e0, 1, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(e1x, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(e1y, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(e1z, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(e2x, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(e2y, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(e2z, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(phi, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)

    CALL mpi_recv(nelect, 1, mpi_integer, 0, 2, mpi_comm_world, is, ic)

    CALL mpi_recv(nion, 1, mpi_integer, 0, 2, mpi_comm_world, is, ic)
    CALL mpi_recv(scaleion, 1, mpi_double_precision, 0, 2, mpi_comm_world, is, ic)

    CALL mpi_recv(charge, 1, mpi_double_precision, 0, 2, mpi_comm_world, is, ic)
    CALL mpi_recv(nspdw, 1, mpi_integer, 0, 2, mpi_comm_world, is, ic)
    CALL mpi_recv(temp, 1, mpi_double_precision, 0, 2, mpi_comm_world, is, ic)
    CALL mpi_recv(occmix, 1, mpi_double_precision, 0, 2, mpi_comm_world, is, ic)
    CALL mpi_recv(b2occ, 1, mpi_double_precision, 0, 2, mpi_comm_world, is, ic)
    CALL mpi_recv(gamocc, 1, mpi_double_precision, 0, 2, mpi_comm_world, is, ic)
    CALL mpi_recv(deocc, 1, mpi_double_precision, 0, 2, mpi_comm_world, is, ic)
    CALL mpi_recv(osfac, 1, mpi_double_precision, 0, 2, mpi_comm_world, is, ic)
    CALL mpi_recv(init_ao, 1, mpi_integer, 0, 2, mpi_comm_world, is, ic)

    CALL mpi_recv(dx, 1, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(dy, 1, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(dz, 1, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, ic)

    CALL mpi_recv(endcon, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(itback, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(epswf, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(e0dmp, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(epsoro, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(dpolx, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(dpoly, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(dpolz, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
#if(extended)
    CALL mpi_recv(ipseudo, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(idenspl, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(i3dx, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(i3dz, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(i3dstate, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(eproj, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(nproj, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    WRITE (*, *) 'nproj=', nproj
    CALL mpi_recv(nproj_states, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    WRITE (*, *) 'nproj_states=', nproj_states
    IF (nproj_states > 0) THEN
      ALLOCATE (proj_states(nproj_states))
      DO ii = 1, nproj_states
        CALL mpi_recv(proj_states(ii), 1, mpi_integer, 0, mpi_any_tag, &
                      mpi_comm_world, is, ic)
        WRITE (*, *) 'proj_states(', ii, ')=', proj_states(ii)
      END DO
    END IF
#endif
#if(raregas)
    CALL mpi_recv(ivdw, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(nc, 1, mpi_integer, 0, 2, mpi_comm_world, is, ic)
    CALL mpi_recv(NE, 1, mpi_integer, 0, 2, mpi_comm_world, is, ic)
    CALL mpi_recv(nk, 1, mpi_integer, 0, 2, mpi_comm_world, is, ic)
    CALL mpi_recv(ekmat, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(surftemp, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(sigmac, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(sigmav, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(sigmak, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(chgc0, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(chge0, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(chgk0, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(cspr, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(mion, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(me, 1, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(mkat, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)

    DO i = 1, 5
      DO j = 1, 5
        CALL mpi_recv(isrtyp(i, j), 1, mpi_integer, 0, mpi_any_tag, &
                      mpi_comm_world, is, ic)
      END DO
    END DO

    CALL mpi_recv(shiftx, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(shifty, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)
    CALL mpi_recv(shiftz, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, ic)

    CALL mpi_recv(iaxis, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
#endif
  END IF

  RETURN
END SUBROUTINE comm_inputparams

!-----comm_ionconfig-------------------------------------------------

SUBROUTINE comm_ionconfig()

! Communicates ion configuration from node "0" to other nodes.

  USE params
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)
  INTEGER :: ic, ion, nod
!------------------------------------------------------------------

  IF (myn == 0 .AND. knode /= 1) THEN

    DO nod = 1, knode - 1
      CALL mpi_send(cx(1), nion, mpi_double_precision, nod, 1, mpi_comm_world, ic)
      CALL mpi_send(cy(1), nion, mpi_double_precision, nod, 1, mpi_comm_world, ic)
      CALL mpi_send(cz(1), nion, mpi_double_precision, nod, 1, mpi_comm_world, ic)
      CALL mpi_send(np(1), nion, mpi_integer, nod, 1, mpi_comm_world, ic)
      CALL mpi_send(cpx(1), nion, mpi_double_precision, nod, 1, mpi_comm_world, ic)
      CALL mpi_send(cpy(1), nion, mpi_double_precision, nod, 1, mpi_comm_world, ic)
      CALL mpi_send(cpz(1), nion, mpi_double_precision, nod, 1, mpi_comm_world, ic)
      IF (init_ao) THEN
        CALL mpi_send(radini(1), nion, mpi_double_precision, &
                      nod, 1, mpi_comm_world, ic)
        CALL mpi_send(ipol(1), nion, mpi_integer, &
                      nod, 1, mpi_comm_world, ic)
        CALL mpi_send(nmaxst(1), nion, mpi_integer, &
                      nod, 1, mpi_comm_world, ic)
        DO ion = 1, nion
          WRITE (6, *) 'myn', myn, 'ion', ion, 'ipol', ipol(ion)
          CALL mpi_send(initord(1, ion), 3, mpi_integer, nod, 1, mpi_comm_world, ic)
        END DO
      END IF
      CALL mpi_send(tblock, ng + 1, mpi_logical, nod, 1, mpi_comm_world, ic)
      IF (ipsptyp == 1) THEN
        CALL mpi_send(tnonloc, ng + 1, mpi_logical, nod, 1, mpi_comm_world, ic)
        CALL mpi_send(tnonlocany, 1, mpi_logical, nod, 1, mpi_comm_world, ic)
      END IF
    END DO

  ELSE IF (myn /= 0 .AND. knode /= 1) THEN

    CALL mpi_recv(cx(1), nion, mpi_double_precision, 0, &
                  mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(cy(1), nion, mpi_double_precision, 0, &
                  mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(cz(1), nion, mpi_double_precision, 0, &
                  mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(np(1), nion, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(cpx(1), nion, mpi_double_precision, 0, &
                  mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(cpy(1), nion, mpi_double_precision, 0, &
                  mpi_any_tag, mpi_comm_world, is, ic)
    CALL mpi_recv(cpz(1), nion, mpi_double_precision, 0, &
                  mpi_any_tag, mpi_comm_world, is, ic)
    IF (init_ao) THEN
      CALL mpi_recv(radini(1), nion, mpi_double_precision, &
                    0, mpi_any_tag, mpi_comm_world, is, ic)
      CALL mpi_recv(ipol(1), nion, mpi_integer, &
                    0, mpi_any_tag, mpi_comm_world, is, ic)
      DO ion = 1, nion
        WRITE (6, *) 'myn', myn, 'ion', ion, 'ipol', ipol(ion)
      END DO
      CALL mpi_recv(nmaxst(1), nion, mpi_integer, &
                    0, mpi_any_tag, mpi_comm_world, is, ic)
      DO ion = 1, nion
        CALL mpi_recv(initord(1, ion), 3, mpi_integer, 0, &
                      mpi_any_tag, mpi_comm_world, is, ic)
      END DO
    END IF
    CALL mpi_recv(tblock, ng + 1, mpi_logical, 0, &
                  mpi_any_tag, mpi_comm_world, is, ic)
    IF (ipsptyp == 1) THEN
      CALL mpi_recv(tnonloc, ng + 1, mpi_logical, 0, mpi_any_tag, mpi_comm_world, is, ic)
      CALL mpi_recv(tnonlocany, 1, mpi_logical, 0, mpi_any_tag, mpi_comm_world, is, ic)
    END IF

  END IF

  RETURN
END SUBROUTINE comm_ionconfig
!-----comm_periodic-------------------------------------------------

SUBROUTINE comm_periodic()

! Communicates table of PsP from node "0" to other nodes.

  USE params
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)
  INTEGER :: nod
!------------------------------------------------------------------

  IF (myn == 0 .AND. knode /= 1) THEN

    DO nod = 1, knode - 1
      CALL mpi_send(ch, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(amu, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(cc1, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(cc2, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(crloc, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(crs, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(chs, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(chg1, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(chg2, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(sgm1, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(sgm2, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(r0g, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(r1g, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(r2g, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(h0_11g, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(h0_22g, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(h0_33g, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(h1_11g, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(h1_22g, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(h2_11g, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(h0_12g, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(radiong, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(radjel, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(surjel, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(bbeta, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(gamma, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(beta4, 185, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
    END DO

  ELSE IF (myn /= 0 .AND. knode /= 1) THEN

    CALL mpi_recv(ch, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(amu, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(cc1, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(cc2, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(crloc, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(crs, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(chs, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(chg1, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(chg2, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(sgm1, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(sgm2, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(r0g, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(r1g, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(r2g, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(h0_11g, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(h0_22g, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(h0_33g, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(h1_11g, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(h1_22g, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(h2_11g, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(h0_12g, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(radiong, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(radjel, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(surjel, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(bbeta, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(gamma, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(beta4, 185, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)

  END IF

  RETURN
END SUBROUTINE comm_periodic

!-----comm_siman-------------------------------------------------

SUBROUTINE comm_simann()

! Communicates parameters for simulated annealing

  USE params
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)
  INTEGER :: n, nod
!------------------------------------------------------------------

  n = myn

  IF (n == 0 .AND. knode /= 1) THEN

    DO nod = 1, knode - 1
      CALL mpi_send(ionsin, 1, mpi_integer, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(iknow, 1, mpi_integer, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(facann, 1, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(nrun, 1, mpi_integer, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(nloop1, 1, mpi_integer, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(nloop2, 1, mpi_integer, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(cptemp, 1, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(delpos, 1, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(ERR, 1, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(errks0, 1, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(trfac2, 1, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(prfac2, 1, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(errsim, 1, 1, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(ncsim, 1, mpi_integer, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(errtot, 1, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(ncon, 1, mpi_integer, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(erfac1, 1, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(trfac1, 1, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(prfac1, 1, mpi_double_precision, nod, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(rand_seed, isize_seed, mpi_integer, nod, 1, mpi_comm_world, mpi_ierror)
    END DO

  ELSE IF (n /= 0 .AND. knode /= 1) THEN

    CALL mpi_recv(ionsin, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(iknow, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(facann, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(nrun, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(nloop1, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(nloop2, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(cptemp, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(delpos, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(ERR, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(errks0, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(trfac2, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(prfac2, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(errsim, 1, 1, mpi_double_precision, 0, &
                  mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(ncsim, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(errtot, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(ncon, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(erfac1, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(trfac1, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(prfac1, 1, mpi_double_precision, 0, mpi_any_tag, &
                  mpi_comm_world, is, mpi_ierror)
    CALL mpi_recv(rand_seed, isize_seed, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, mpi_ierror)
  END IF

  RETURN
END SUBROUTINE comm_simann

! *****************************

SUBROUTINE pi_allreduce(rin, rout, kdf, id, id1, id2, ic)

! *****************************

  USE params
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  REAL, INTENT(IN OUT) :: rin(kdf)
  REAL, INTENT(IN OUT) :: rout(kdf)
  INTEGER, INTENT(IN OUT) :: kdf
  INTEGER, INTENT(IN OUT) :: id
  INTEGER, INTENT(IN OUT) :: id1
  INTEGER, INTENT(IN OUT) :: id2
  INTEGER, INTENT(IN OUT) :: ic

  CALL mpi_allreduce(rin, rout, kdf, mpi_double_precision, &
                     mpi_sum, mpi_comm_world, ic)
  RETURN
END SUBROUTINE pi_allreduce

! *****************************

SUBROUTINE pi_allgather(rin, mint, nt, rout, kdf)

! *****************************

  USE params
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  REAL, INTENT(IN OUT) :: rin(kdf)
  INTEGER, INTENT(IN OUT) :: mint
  INTEGER, INTENT(IN OUT) :: nt
  REAL, INTENT(IN OUT) :: rout(kdf)
  INTEGER, INTENT(IN OUT) :: kdf

  CALL mpi_allgather(rin(mint), nt, mpi_double_precision, &
                     rout, nt, mpi_double_precision, mpi_comm_world, mpi_ierror)
  RETURN
END SUBROUTINE pi_allgather

! *****************************

SUBROUTINE i_allreduce(rin, rout, kdf, id, id1, id2, ic)

! *****************************

  USE params
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  REAL, INTENT(OUT) :: rin(kdf)
  REAL, INTENT(IN OUT) :: rout(kdf)
  INTEGER, INTENT(IN OUT) :: kdf
  INTEGER, INTENT(IN OUT) :: id
  INTEGER, INTENT(IN OUT) :: id1
  INTEGER, INTENT(IN OUT) :: id2
  INTEGER, INTENT(IN OUT) :: ic

  INTEGER :: is(mpi_status_size)
  INTEGER :: i, it, j, l, mx, n, nprocs, nt
  CALL mpi_comm_size(mpi_comm_world, nprocs, ic)
  IF (nprocs /= kstate) STOP
  CALL mpi_comm_rank(mpi_comm_world, n, ic)
  l = LOG(kstate*1D0)/LOG(2D0) - 1
  DO it = 0, l
    j = (n/(2D0**REAL(it, DP)))
    IF (MOD(j, 2) == 0) THEN
      nt = n + 2D0**REAL(it, DP)
    ELSE
      nt = n - 2D0**REAL(it, DP)
    END IF
    IF (MOD(j, 2) == 0) THEN
      CALL mpi_send(rin, mx, mpi_double_precision, nt, 0, mpi_comm_world, ic)
      CALL mpi_recv(rout, mx, mpi_double_precision, nt, 0, mpi_comm_world, is, ic)
    ELSE
      CALL mpi_recv(rout, mx, mpi_double_precision, nt, 0, mpi_comm_world, is, ic)
      CALL mpi_send(rin, mx, mpi_double_precision, nt, 0, mpi_comm_world, ic)
    END IF
    DO i = 1, mx
      rin(i) = rin(i) + rout(i)
    END DO
  END DO
  DO i = 1, mx
    rout(i) = rin(i)
  END DO
  RETURN
END SUBROUTINE i_allreduce

! *****************************

SUBROUTINE i_allgather(rin, mint, ntr, rout, kdf)

! *****************************

  USE params
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  REAL, INTENT(IN) :: rin(kdf)
  INTEGER, INTENT(IN) :: mint
  INTEGER, INTENT(IN) :: ntr
  REAL, INTENT(OUT) :: rout(kdf)
  INTEGER, INTENT(IN OUT) :: kdf
  INTEGER :: is(mpi_status_size)

  INTEGER :: i, it, j, l, maxr, minr, n, nax, nin, nt

  CALL mpi_comm_rank(mpi_comm_world, n, mpi_ierror)

  l = LOG(kstate*1D0)/LOG(2D0) - 1
  nax = ntr
  nin = mint
  DO i = nin, nin + nax
! dir$ unroll 8
    rout(i) = rin(i)
  END DO
  DO it = 0, l
    j = (n/(2D0**REAL(it, DP)))
    IF (MOD(j, 2) == 0) THEN
      nt = n + 2D0**REAL(it, DP)
    ELSE
      nt = n - 2D0**REAL(it, DP)
    END IF
    IF (MOD(j, 2) == 0) THEN
      CALL mpi_send(nax, 1, mpi_integer, nt, 1, mpi_comm_world, mpi_ierror)
      CALL mpi_send(nin, 1, mpi_integer, nt, 2, mpi_comm_world, mpi_ierror)
      CALL mpi_send(rout(nin), nax, mpi_double_precision, &
                    nt, 3, mpi_comm_world, mpi_ierror)
      CALL mpi_recv(maxr, 1, mpi_integer, nt, 4, mpi_comm_world, is, mpi_ierror)
      CALL mpi_recv(minr, 1, mpi_integer, nt, 5, mpi_comm_world, is, mpi_ierror)
      CALL mpi_recv(rout(minr), maxr, mpi_double_precision, &
                    nt, 6, mpi_comm_world, is, mpi_ierror)
    ELSE
      CALL mpi_recv(maxr, 1, mpi_integer, nt, 1, mpi_comm_world, is, mpi_ierror)
      CALL mpi_recv(minr, 1, mpi_integer, nt, 2, mpi_comm_world, is, mpi_ierror)
      CALL mpi_recv(rout(minr), maxr, mpi_double_precision, &
                    nt, 3, mpi_comm_world, is, mpi_ierror)
      CALL mpi_send(nax, 1, mpi_integer, nt, 4, mpi_comm_world, mpi_ierror)
      CALL mpi_send(nin, 1, mpi_integer, nt, 5, mpi_comm_world, mpi_ierror)
      CALL mpi_send(rout(nin), nax, mpi_double_precision, &
                    nt, 6, mpi_comm_world, mpi_ierror)
      nin = nin - maxr
    END IF
    nax = nax + maxr

  END DO

  RETURN
END SUBROUTINE i_allgather
! *****************************

SUBROUTINE pi_scatter(buf)

! *****************************

  USE params
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)

  REAL, INTENT(IN OUT) :: buf

  INTEGER :: ic, nod

  IF (myn == 0 .AND. knode /= 1) THEN
    DO nod = 1, knode - 1
      CALL mpi_send(buf, 1, mpi_double_precision, nod, 1, mpi_comm_world, ic)
    END DO
  ELSE IF (myn /= 0 .AND. knode /= 1) THEN
    CALL mpi_recv(buf, 1, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, is, ic)
  ELSE
    STOP ' PI_SCATTER1 ONLY for KNODE > 1'
  END IF

END SUBROUTINE pi_scatter

SUBROUTINE i_scatter(buf)

! *****************************

  USE params
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)

  INTEGER, INTENT(IN OUT) :: buf

  INTEGER :: ic, nod

  IF (myn == 0 .AND. knode /= 1) THEN
    DO nod = 1, knode - 1
      CALL mpi_send(buf, 1, mpi_integer, nod, 1, mpi_comm_world, ic)
    END DO
  ELSE IF (myn /= 0 .AND. knode /= 1) THEN
    CALL mpi_recv(buf, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, is, ic)
  ELSE
    STOP ' PI_SCATTER1 ONLY for KNODE > 1'
  END IF

END SUBROUTINE i_scatter

!-----prispe_parallele-------------------------------------------------

SUBROUTINE prispe_parallele(iunit, it)

! Collects single particle properties and writes it on
! UNIT 'iunit'

  USE params
  USE util, ONLY: safeopen, tmodbuf
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iunit

  INCLUDE 'mpif.h'

  INTEGER :: i, ic, it, k, n, nact, nb, nod, nod2
  INTEGER :: is(mpi_status_size)
  REAL(DP) :: prisav(kstate, 5)
  REAL(DP) :: prisav_all(ksttot, 5)
  INTEGER :: iprisav(kstate, 2)

  LOGICAL, PARAMETER :: ttest = .FALSE.

!------------------------------------------------------------------

! collect into aux fields
  DO nb = 1, nstate
    prisav(nb, 1) = occup(nb)
    prisav(nb, 2) = ekinsp(nb)
    prisav(nb, 3) = amoy(nb)
    prisav(nb, 4) = spvariance(nb)
    prisav(nb, 5) = spnorm(nb)
    iprisav(nb, 1) = nrel2abs(nb)
    iprisav(nb, 2) = 3 - 2*ispin(nrel2abs(nb))
  END DO

! communicate and print
  prisav_all = 0D0
  IF (myn /= 0) THEN
    nod = myn
    CALL mpi_send(prisav, 5*kstate, mpi_double_precision, 0, nod, mpi_comm_world, ic)
    CALL mpi_send(iprisav, 2*kstate, mpi_integer, 0, nod, mpi_comm_world, ic)
    IF (ttest) WRITE (*, *) ' INFO: sent at node:', myn
    DO n = 1, nstate
      nact = nrel2abs_other(n, 0)
      prisav_all(nact, :) = prisav(n, :)
    END DO
  ELSE
    DO nod2 = 0, knode - 1
      IF (nod2 > 0) THEN
        CALL mpi_recv(prisav, 5*kstate, mpi_double_precision, &
                      nod2, mpi_any_tag, mpi_comm_world, is, ic)
        CALL mpi_recv(iprisav, 2*kstate, mpi_integer, &
                      nod2, mpi_any_tag, mpi_comm_world, is, ic)
        IF (ttest) WRITE (*, *) ' INFO: recv from node=', nod2
      END IF
      DO nb = 1, nstate_node(nod2)
        nact = nrel2abs_other(nb, nod2)
        prisav_all(nact, :) = prisav(nb, :)
        IF (numspin == 2) THEN
          WRITE (iunit, '(a,i3,a,i3,3f9.5,1pg12.4)') 'level:', iprisav(nb, 1), &
            ' spin,occup,ekin,esp,variance =', iprisav(nb, 2), (prisav(nb, k), k=1, 4)
        ELSE
          WRITE (iunit, '(a,i3,a,3f9.5,1pg12.4)') &
            'level:', nrel2abs(nb), &
            ' occup,ekin,esp,variance=', (prisav(nb, k), k=1, 4)
        END IF
      END DO
    END DO

    IF (tmodbuf(it, jstinf)) THEN
      CALL safeopen(77, it, jstinf, 'pspenergies')
      WRITE (77, '(1f15.6,500f12.6)') tfs, (prisav_all(nb, 3), nb=1, nstate_all)
      CALL flush (77)
      CALL safeopen(76, it, jstinf, 'pspvariances')
      WRITE (76, '(1f15.6,500f12.6)') tfs, (prisav_all(nb, 4), nb=1, nstate_all)
      CALL flush (76)
    END IF

    IF (jesc > 0 .AND. tmodbuf(it, jnorms)) THEN
      CALL safeopen(806, it, jnorms, 'pescOrb')
      WRITE (806, '(1f15.6,500f12.8)') tfs, 1D0 - prisav_all(1:nstate_all, 5)
      CALL flush (806)
    END IF

  END IF

  RETURN
END SUBROUTINE prispe_parallele

!-----init_boxpara-------------------------------------------------

SUBROUTINE init_boxpara()

! Initializes the length of the subarrays used in pi_scatterv
! and pi_gatherv.

  USE params
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)

  INTEGER :: i, ntranche, rest, nprocs
!------------------------------------------------------------------

  CALL mpi_comm_rank(mpi_comm_world, myn, mpi_ierror)
  CALL mpi_comm_size(mpi_comm_world, nprocs, mpi_ierror)

  ALLOCATE (lengnod(nprocs))
  ALLOCATE (displ(0:nprocs - 1))

  IF (myn == 0) THEN
    ntranche = nxyz/nprocs
    rest = MOD(nxyz, nprocs)
    displ(0) = 0
    DO i = 1, nprocs - 1
      lengnod(i) = ntranche
      displ(i) = lengnod(i) + displ(i - 1)
    END DO
    lengnod(nprocs) = ntranche + rest
  END IF
  CALL mpi_bcast(lengnod, nprocs, mpi_integer, 0, mpi_comm_world, mpi_ierror)
  CALL mpi_bcast(displ, nprocs, mpi_integer, 0, mpi_comm_world, mpi_ierror)

  RETURN
END SUBROUTINE init_boxpara

!-----pi_scatterv-------------------------------------------------

SUBROUTINE pi_scatterv(sendbuf, N, recvbuf, sizenod, ic)

! Scatters from node 0 a partition of the array 'sendbuf' of
! size 'N'=(nxyz) into nprocs subarrays, locally called on the current
! process 'recvbuf', which size is stored in the array 'lengnod'.

  USE params
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)

  INTEGER, INTENT(IN OUT) :: ic
  INTEGER, INTENT(IN) :: sizenod
  INTEGER, INTENT(IN) :: N
  REAL(DP), INTENT(IN OUT) :: sendbuf(N)
  REAL(DP), INTENT(IN OUT) :: recvbuf(sizenod)

  IF (N /= sum(lengnod)) &
    STOP 'wrong SIZE for sender array IN pi_scatterv'

!ALLOCATE(recvbuf(sizenod))
!WRITE(*,*) 'allocation of recvbuf'

  CALL mpi_scatterv(sendbuf, lengnod, displ, mpi_double_precision, &
                    recvbuf, sizenod, mpi_double_precision, 0, mpi_comm_world, ic)
  RETURN
END SUBROUTINE pi_scatterv

!-----pi_allgatherv-------------------------------------------------

SUBROUTINE pi_allgatherv(sendbuf, sizenod, recvbuf, N, ic)

! Gathers all subarrays 'sendbuf' of size sizenod, into the
! array 'recvbf' of size 'N'(=nxyz), and communicates it to all
! processes.

  USE params
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)

  INTEGER, INTENT(IN OUT) :: ic
  INTEGER, INTENT(IN) :: sizenod
  INTEGER, INTENT(IN) :: N
  REAL(DP), INTENT(IN OUT) :: sendbuf(nxyz)
  REAL(DP), INTENT(IN OUT) :: recvbuf(N)

  CALL mpi_allgatherv(sendbuf, sizenod, mpi_double_precision, &
                      recvbuf, lengnod, displ, mpi_double_precision, mpi_comm_world, ic)
  RETURN
END SUBROUTINE pi_allgatherv
#endif
