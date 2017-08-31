SUBROUTINE rk

  USE mpi
  USE Euler
  USE Lagrange
  USE para

  DOUBLE PRECISION:: fac
  INTEGER:: krk

  ! krk=1:
  krk=1
  fac=0.5d0
  CALL Euler_reset(krk)
  IF(isingular.EQV..TRUE.) THEN
     IF(object_move.EQV..TRUE.) THEN
        CALL surface_property
        CALL Raycrossing
        CALL allocate_jcs(krk)
     ENDIF
     CALL jc_firstsecond(krk)
     CALL principle_jc_p
     CALL mac_distribute
     CALL jump_contribution(krk)
  ELSE
     CALL interpolateuv(0)
  ENDIF
  CALL uv_strain
  CALL old_save
  CALL rhs4p(fac)
  CALL solvePoisson(1)
  CALL rhs4uv(krk)
  CALL updateuv(krk,fac)
  CALL exchange4ghost
  IF(object_move.EQV..TRUE.) THEN
     CALL free_jcs
     CALL surface_move(krk,fac)
  ENDIF

  ! krk=2:
  krk=2
  fac=1.0d0
  !! if use rk4, fac=0.5d0 !!
  CALL Euler_reset(krk)
  IF(isingular.EQV..TRUE.) THEN
     IF(object_move.EQV..TRUE.) THEN
        CALL surface_property
        CALL Raycrossing
        CALL allocate_jcs(krk)
     ENDIF
     CALL jc_firstsecond(krk)
     CALL principle_jc_p
     CALL mac_distribute
     CALL jump_contribution(krk)
  ELSE
     CALL interpolateuv(0)
  ENDIF
  CALL uv_strain
  CALL rhs4p(fac)
  CALL solvePoisson(1)
  CALL rhs4uv(krk)
  CALL updateuv(krk,fac)
  CALL exchange4ghost
  IF(object_move.EQV..TRUE.) THEN
     CALL free_jcs
     CALL surface_move(krk,fac)
  ENDIF

  ! krk=3:
  krk=3
  fac=1.0d0
  CALL Euler_reset(krk)
  IF(isingular.EQV..TRUE.) THEN
     IF(object_move.EQV..TRUE.) THEN
        CALL surface_property
        CALL Raycrossing
        CALL allocate_jcs(krk)
     ENDIF
     CALL jc_firstsecond(krk)
     CALL principle_jc_p
     CALL mac_distribute
     CALL jump_contribution(krk)
  ELSE
     CALL interpolateuv(0)
  ENDIF
  CALL uv_strain
  CALL rhs4p(fac)
  CALL solvePoisson(1)
  CALL rhs4uv(krk)
  CALL updateuv(krk,fac)
  CALL exchange4ghost
  IF(object_move.EQV..TRUE.) THEN
     CALL free_jcs
     CALL surface_move(krk,fac)
  ENDIF

END SUBROUTINE rk
