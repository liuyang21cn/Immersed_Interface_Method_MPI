SUBROUTINE updateuv(krk,fac)

  USE para
  USE euler
  DOUBLE PRECISION:: fac
  INTEGER:: krk

  IF(krk==1) THEN
     DO j=jubeg,juend
        DO i=iubeg,iuend
           u(i,j)=un(i,j)+fac*dt*rhsu1(i,j)
        ENDDO
     ENDDO
     DO j=jvbeg,jvend
        DO i=ivbeg,ivend
           v(i,j)=vn(i,j)+fac*dt*rhsv1(i,j)
        ENDDO
     ENDDO

  ELSEIF(krk==2) THEN
     DO j=jubeg,juend
        DO i=iubeg,iuend
           u(i,j)=un(i,j)+fac*dt*rhsu2(i,j)
        ENDDO
     ENDDO
     DO j=jvbeg,jvend
        DO i=ivbeg,ivend
           v(i,j)=vn(i,j)+fac*dt*rhsv2(i,j)
        ENDDO
     ENDDO

  ELSEIF(krk==3) THEN

     !!========= rk3 ============!!!
     DO j=jubeg,juend
        DO i=iubeg,iuend
           u(i,j)=un(i,j)+(dt/6.0d0)*(rhsu1(i,j)+4.0d0*rhsu2(i,j)+rhsu3(i,j))
        ENDDO
     ENDDO
     DO j=jvbeg,jvend
        DO i=ivbeg,ivend
           v(i,j)=vn(i,j)+(dt/6.0d0)*(rhsv1(i,j)+4.0d0*rhsv2(i,j)+rhsv3(i,j))
        ENDDO
     ENDDO


     !=========== rk4 ================!!
     !    DO j=jubeg,juend
     !       DO i=iubeg,iuend
     !          u(i,j)=un(i,j)+fac*dt*rhsu3(i,j)
     !       ENDDO
     !    ENDDO
     !    DO j=jvbeg,jvend
     !       DO i=ivbeg,ivend
     !          v(i,j)=vn(i,j)+fac*dt*rhsv3(i,j)
     !       ENDDO
     !    ENDDO
     !
     ! ELSEIF(krk==4) THEN
     !    DO j=jubeg,juend
     !       DO i=iubeg,iuend
     !          u(i,j)=un(i,j)+(dt/6.0d0)*(rhsu1(i,j)+2.0d0*(rhsu2(i,j)+rhsu3(i,j))+rhsu4(i,j))
     !       ENDDO
     !    ENDDO
     !    DO j=jvbeg,jvend
     !       DO i=ivbeg,ivend
     !          v(i,j)=vn(i,j)+(dt/6.0d0)*(rhsv1(i,j)+2.0d0*(rhsv2(i,j)+rhsv3(i,j))+rhsv4(i,j))
     !       ENDDO
     !    ENDDO
  ENDIF

  CALL velocity_reset
  CALL ubc
  CALL vbc

END SUBROUTINE updateuv
