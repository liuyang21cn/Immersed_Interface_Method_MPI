SUBROUTINE rhs4uv(krk)

  USE MPI
  USE para
  USE Euler

  INTEGER,INTENT(in):: krk
  DOUBLE PRECISION:: conv,visc,grad

  ! DO j=jubeg,juend
  !    DO i=iubeg,iuend-1
  DO j=jubeg,juend
     DO i=iubeg,iuend
        IF(i>0 .AND. i<nx) THEN
           IF(j>0 .AND. j<=ny) THEN
              conv=-xifx(i) *(ucc(i+1,j)*ucc(i+1,j)-ucc(i,j)*ucc(i,j))/dxi &
                   -etacy(j)*(uee(i,j)*vee(i,j)-uee(i,j-1)*vee(i,j-1))/deta
              grad=-xifx(i)*(p(i+1,j)-p(i,j))/dxi
              visc=Re1*xifx(i)*(xicx(i)*u(i-1,j)-(xicx(i)+xicx(i+1))*u(i,j)+xicx(i+1)*u(i+1,j))/dxi/dxi+&
                   Re1*etacy(j)*(etafy(j-1)*u(i,j-1)-(etafy(j-1)+etafy(j))*u(i,j)+etafy(j)*u(i,j+1))/deta/deta
              IF(krk.EQ.1) THEN
                 rhsu1(i,j) = rhsu1(i,j)+conv+grad+visc
              ELSEIF ( krk.EQ.2 ) THEN
                 rhsu2(i,j) = rhsu2(i,j)+conv+grad+visc
              ELSEIF ( krk.EQ.3 ) THEN
                 rhsu3(i,j) = rhsu3(i,j)+conv+grad+visc
              ELSE
                 rhsu4(i,j) = rhsu4(i,j)+conv+grad+visc
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDDO

  ! DO j=jvbeg,jvend-1
  !  DO i=ivbeg,ivend
  DO j=jvbeg,jvend
     DO i=ivbeg,ivend
        IF(j>0 .AND. j<ny) THEN
           IF(i>0 .AND. i<=nx) THEN
              conv=-etafy(j)*(vcc(i,j+1)*vcc(i,j+1)-vcc(i,j)*vcc(i,j))/deta &
                   -xicx(i) *(uee(i,j)*vee(i,j)-uee(i-1,j)*vee(i-1,j))/dxi
              grad=-etafy(j)*(p(i,j+1)-p(i,j))/deta
              visc=Re1*xicx(i)* (xifx(i-1)*v(i-1,j)-(xifx(i-1)+xifx(i))*v(i,j)  &
                   +xifx(i)*v(i+1,j))/dxi/dxi + &
                   Re1*etafy(j)*(etacy(j)*v(i,j-1)-(etacy(j)+etacy(j+1))*v(i,j) &
                   +etacy(j+1)*v(i,j+1))/deta/deta
              IF(krk.EQ.1) THEN
                 rhsv1(i,j)=rhsv1(i,j)+conv+grad+visc
              ELSEIF ( krk.EQ.2 ) THEN
                 rhsv2(i,j)=rhsv2(i,j)+conv+grad+visc
              ELSEIF ( krk.EQ.3 ) THEN
                 rhsv3(i,j)=rhsv3(i,j)+conv+grad+visc
              ELSE
                 rhsv4(i,j)=rhsv4(i,j)+conv+grad+visc
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDDO


END SUBROUTINE rhs4uv
