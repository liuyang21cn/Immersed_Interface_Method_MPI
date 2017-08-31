SUBROUTINE free_jcs

  USE para
  USE Euler
  USE Lagrange

  IF(object_move.EQV..TRUE.) THEN

     ! Raycrossing
     DEALLOCATE(xc_int_vertex)
     DEALLOCATE(xf_int_vertex)
     DEALLOCATE(yc_int_vertex)
     DEALLOCATE(yf_int_vertex)
     DEALLOCATE(xc_int_info)
     DEALLOCATE(xf_int_info)
     DEALLOCATE(yc_int_info)
     DEALLOCATE(yf_int_info)

     ! jc_firstsecond
     DEALLOCATE(ujc)
     DEALLOCATE(vjc)
     DEALLOCATE(pjc)

     ! mac_distribute
     DEALLOCATE(ujcxf)
     DEALLOCATE(vjcxf)
     DEALLOCATE(vjcxc)
     DEALLOCATE(ujcxc)
     DEALLOCATE(pjcxc)
     DEALLOCATE(ujcyf)
     DEALLOCATE(vjcyf)
     DEALLOCATE(vjcyc)
     DEALLOCATE(ujcyc)
     DEALLOCATE(pjcyc)

  ENDIF


END SUBROUTINE free_jcs
