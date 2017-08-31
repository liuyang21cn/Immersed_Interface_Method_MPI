SUBROUTINE Divergence_reset

  USE para
  USE Euler
  USE Lagrange

  INTEGER:: ixc,ixc1,jyc,jyc1,index,myvertex
  DOUBLE PRECISION:: signx,signy
  DOUBLE PRECISION:: w1,w2

  DO j=jdbeg,jdend
     DO i=idbeg,idend
        DO myobj=1,nobj4proc
           IF(iop(i,j)==obj4proc(myobj)) THEN
              d(i,j)=0.0d0
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  DO index = 1, n_xc_int
     myvertex=xc_int_vertex(1,index)
     signy = xc_int_info(4, index)
     ixc = xc_int_vertex(2, index)
     jyc = xc_int_vertex(3, index)
     jyc1 = jyc+1
     IF(signy.GE.0.0d0) THEN
        w1=(yc(jyc1+1)-yc(jyc1)  )/(yc(jyc1+2)-yc(jyc1))
        w2=(yc(jyc1+2)-yc(jyc1+1))/(yc(jyc1+2)-yc(jyc1))
        d(ixc,jyc1)=(d(ixc,jyc1+1)-w2*d(ixc,jyc1+2))/w1
     ELSE
        w1=(yc(jyc-1)-yc(jyc-2))/(yc(jyc)-yc(jyc-2))
        w2=(yc(jyc)  -yc(jyc-1))/(yc(jyc)-yc(jyc-2))
        d(ixc,jyc)=(d(ixc,jyc-1)-w1*d(ixc,jyc-2))/w2
     ENDIF
  ENDDO

  DO index = 1, n_yc_int
     signx = yc_int_info(3, index)
     jyc = yc_int_vertex(2, index)
     ixc = yc_int_vertex(3, index)
     ixc1 = ixc+1
     IF(signx.GE.0.0d0) THEN
        w1=(xc(ixc1+1)-xc(ixc1)  )/(xc(ixc1+2)-xc(ixc1))
        w2=(xc(ixc1+2)-xc(ixc1+1))/(xc(ixc1+2)-xc(ixc1))
        d(ixc1,jyc)=(d(ixc1+1,jyc)-w2*d(ixc1+2,jyc))/w1
     ELSE
        w1=(xc(ixc-1)-xc(ixc-2))/(xc(ixc)-xc(ixc-2))
        w2=(xc(ixc)  -xc(ixc-1))/(xc(ixc)-xc(ixc-2))
        d(ixc,jyc)=(d(ixc-1,jyc)-w1*d(ixc-2,jyc))/w2
     ENDIF
  ENDDO

END SUBROUTINE Divergence_reset
