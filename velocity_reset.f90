SUBROUTINE velocity_reset

  USE para
  USE Euler
  USE Lagrange

  DOUBLE PRECISION xx,yy

  DO j=jubeg,juend
     DO i=iubeg,iuend
        DO myobj=1,nobj4proc
           xx=xf(i)
           yy=yc(j)
           IF(iou(i,j)==obj4proc(myobj)) THEN
              u(i,j)=xsct(myobj)-thetat(myobj)*yy
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  DO j=jvbeg,jvend
     DO i=ivbeg,ivend
        DO myobj=1,nobj4proc
           xx=xc(i)
           yy=yf(j)
           IF(iov(i,j)==obj4proc(myobj)) THEN
              v(i,j)=ysct(myobj)+thetat(myobj)*xx
           ENDIF
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE velocity_reset
