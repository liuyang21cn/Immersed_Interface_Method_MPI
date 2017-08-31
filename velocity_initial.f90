SUBROUTINE velocity_initial

  USE mpi
  USE para
  USE Euler
  USE Lagrange

  DOUBLE PRECISION::grad,ucpdx,vcpdy
  INTEGER:: ixc,jyc,myvertex,krk

  krk=-1

  ! initialize rhsp
  IF(isingular.EQV..TRUE.) THEN
     ! find intersection points and allocate related jcs
     CALL surface_property
     CALL Raycrossing
     CALL allocate_jcs(0)
     ! compute necessary jcs
     CALL jc_firstsecond(krk)
     CALL principle_jc_p
     CALL pressure_distribute
     ! compute jump contributions
     DO index = 1, n_yc_int
        CALL correction_pressure(index,3,ucpdx)
        jyc = yc_int_vertex(2, index)
        ixc = yc_int_vertex(3, index)
        IF(ixc<=iuend-1 .AND. ixc>= iubeg+1) THEN
           IF(jyc <=juend-1 .AND. jyc>=jubeg+1) u(ixc,jyc)=u(ixc,jyc)-ucpdx
        ENDIF
     ENDDO
     DO index = 1, n_xc_int
        CALL correction_pressure(index,1,vcpdy)
        ixc = xc_int_vertex(2, index)
        jyc = xc_int_vertex(3, index)
        IF(ixc<=ivend-1 .AND. ixc>= ivbeg+1) THEN
           IF(jyc <=jvend-1 .AND. jyc>=jvbeg+1) v(ixc,jyc)=v(ixc,jyc)-vcpdy
        ENDIF
     ENDDO
  ENDIF

  ! solve pressure
  CALL rhs4p_clean
  CALL solvePoisson(1)

  ! update u & v
  DO j=jubeg,juend
     DO i=iubeg,iuend
        ! DO j=jubeg,juend
        !  DO i=iubeg,iuend-1
        IF(i>0 .AND. i<nx) THEN
           IF(j>0 .AND. j<=ny) THEN
              grad=-xifx(i)*(p(i+1,j)-p(i,j))/dxi
              u(i,j)=u(i,j)+grad
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  DO j=jvbeg,jvend
     DO i=ivbeg,ivend
        ! DO j=jvbeg,jvend-1
        !  DO i=ivbeg,ivend
        IF(j>0 .AND. j<ny) THEN
           IF(i>0 .AND. i<=nx) THEN
              grad=-etafy(j)*(p(i,j+1)-p(i,j))/deta
              v(i,j)=v(i,j)+grad
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  IF(isingular.EQV..TRUE.) CALL velocity_reset
  CALL exchange4ghost
  CALL interpolateuv(0)
  IF(object_move.EQV..TRUE.) THEN
     CALL free_jcs
     CALL surface_move(5,0)
  ENDIF



END SUBROUTINE velocity_initial
