SUBROUTINE correction_pressure(index,flag,dpdx)

  USE para
  USE Euler
  USE Lagrange
  USE stretchxy

  INTEGER:: ixf,ixc,ixc1,jyf,jyc,jyc1,myvertex
  INTEGER, INTENT(in):: flag,index
  DOUBLE PRECISION:: xint,yint,pcdyy(2),vcpdy,pcdxx(2),ucpdx
  DOUBLE PRECISION:: pjc0,pjcdxi,pjcdxi2,pjcdeta,pjcdeta2
  DOUBLE PRECISION, INTENT(out):: dpdx

  dpdx = zero
  IF(flag == 1) THEN
     myvertex=xc_int_vertex(1,index)
     ixc = xc_int_vertex(2,index)
     jyc = xc_int_vertex(3,index)
     jyf = xc_int_vertex(4,index)
     yint = xc_int_info(2,index)
     jyc1=jyc+1

     pjc0 = pjcxc(6,index)
     pjcdeta = pjcxc(2,index)*dydeta(yint)
     pjcdeta2 = pjcxc(4,index)*dydeta(yint)*dydeta(yint) + pjcxc(2,index)*d2ydeta2(yint)

     IF((yc(jyc).LT.yint).AND.(yint.LT.yf(jyc))) THEN
        pcdyy(1)=-etacy(jyc)*etafy(jyc)/deta/deta*&
             (pjc0 + pjcdeta*(etac(jyc)-y2eta(yint)) + 0.5d0*pjcdeta2*(etac(jyc)-y2eta(yint))**2.0d0)&
             -etacy(jyc)/deta*(pjcxc(2,index)+pjcxc(4,index)*dydeta(yint)*(etaf(jyc)-y2eta(yint)))
        pcdyy(2)=etacy(jyc1)*etafy(jyc)/deta/deta*&
             (pjc0 + pjcdeta*(etac(jyc)-y2eta(yint)) + 0.5d0*pjcdeta2*(etac(jyc)-y2eta(yint))**2.0d0)
     ELSEIF((yf(jyc).LT.yint).AND.(yint.LT.yc(jyc1))) THEN
        pcdyy(1)=-etacy(jyc)*etafy(jyc)/deta/deta*&
             (pjc0 + pjcdeta*(etac(jyc1)-y2eta(yint)) + 0.5d0*pjcdeta2*(etac(jyc1)-y2eta(yint))**2.0d0)
        pcdyy(2)=etacy(jyc1)*etafy(jyc)/deta/deta*&
             (pjc0 + pjcdeta*(etac(jyc1)-y2eta(yint)) + 0.5d0*pjcdeta2*(etac(jyc1)-y2eta(yint))**2.0d0)&
             -etacy(jyc1)/deta*(pjcxc(2,index)+pjcxc(4,index)*dydeta(yint)*(etaf(jyc)-y2eta(yint)))
     ENDIF

     IF(jyc.EQ.jyf) THEN
        vcpdy=-1.0d0/deta/dydeta(yint)*(pjc0+ &
             pjcdeta *(etac(jyc1)-y2eta(yint))+ &
             0.5d0*pjcdeta2*(etac(jyc1)-y2eta(yint))*(etac(jyc1)-y2eta(yint)))
     ELSE
        vcpdy=-1.0d0/deta/dydeta(yint)*(pjc0+ &
             pjcdeta *(etac(jyc)- y2eta(yint))+ &
             0.5d0*pjcdeta2*(etac(jyc)- y2eta(yint))*(etac(jyc)- y2eta(yint)))
     ENDIF
     ! update rhsp
     dpdx = vcpdy
     IF(ixc<=ipend .AND. ixc>=ipbeg) THEN
        IF(jyc <=jpend .AND. jyc >=jpbeg) rhsp(ixc,jyc)  =rhsp(ixc,jyc)  -pcdyy(1)*dxi*dxi
        IF(jyc1<=jpend .AND. jyc1>=jpbeg) rhsp(ixc,jyc+1)=rhsp(ixc,jyc+1)-pcdyy(2)*dxi*dxi
     ENDIF
  ENDIF

  IF(flag == 3) THEN
     myvertex=yc_int_vertex(1,index)
     xint = yc_int_info(2, index)
     jyc = yc_int_vertex(2, index)
     ixc = yc_int_vertex(3, index)
     ixf = yc_int_vertex(4, index)
     ixc1=ixc+1

     pjc0 = pjcyc(5,index)
     pjcdxi = pjcyc(1,index)*dxdxi(xint)
     pjcdxi2 = pjcyc(3,index)*dxdxi(xint)*dxdxi(xint) + pjcyc(1,index)*d2xdxi2(xint)

     IF((xc(ixc).LT.xint).AND.(xint.LT.xf(ixc))) THEN
        pcdxx(1)=-xicx(ixc)*xifx(ixc)/dxi/dxi*&
             (pjc0 + pjcdxi*(xic(ixc)-x2xi(xint)) + 0.5d0*pjcdxi2*(xic(ixc)-x2xi(xint))**2.0d0)&
             -xicx(ixc)/dxi*(pjcyc(1,index)+pjcyc(3,index)*dxdxi(xint)*(xif(ixc)-x2xi(xint)))
        pcdxx(2)=xicx(ixc1)*xifx(ixc)/dxi/dxi*&
             (pjc0 + pjcdxi*(xic(ixc)-x2xi(xint)) + 0.5d0*pjcdxi2*(xic(ixc)-x2xi(xint))**2.0d0)
        IF(ixc+1 == 316 .AND. jyc==227) PRINT*,myid,ixc,jyc,pcdxx(1),pcdxx(2),'1'
        IF(ixc+1 == 316 .AND. jyc==255) PRINT*,myid,ixc,jyc,pcdxx(1),pcdxx(2),'1'
     ELSEIF((xf(ixc).LT.xint).AND.(xint.LT.xc(ixc1))) THEN
        pcdxx(1)=-xicx(ixc)*xifx(ixc)/dxi/dxi*&
             (pjc0 + pjcdxi*(xic(ixc1)-x2xi(xint)) + 0.5d0*pjcdxi2*(xic(ixc1)-x2xi(xint))**2.0d0)
        pcdxx(2)=xicx(ixc1)*xifx(ixc)/dxi/dxi*&
             (pjc0 + pjcdxi*(xic(ixc1)-x2xi(xint)) + 0.5d0*pjcdxi2*(xic(ixc1)-x2xi(xint))**2.0d0)&
             -xicx(ixc1)/dxi*(pjcyc(1,index)+pjcyc(3,index)*dxdxi(xint)*(xif(ixc)-x2xi(xint)))
     ENDIF

     IF(ixc.EQ.ixf) THEN
        ucpdx=-1.0d0/dxi/dxdxi(xint)*(pjc0+ &
             pjcdxi *(xic(ixc1)-x2xi(xint))+ &
             0.5d0*pjcdxi2*(xic(ixc1)-x2xi(xint))*(xic(ixc1)-x2xi(xint)))
     ELSE
        ucpdx=-1.0d0/dxi/dxdxi(xint)*(pjc0+ &
             pjcdxi *(xic(ixc)- x2xi(xint))+ &
             0.5d0*pjcdxi2*(xic(ixc)- x2xi(xint))*(xic(ixc)- x2xi(xint)))
     ENDIF

     ! update rhsp
     IF(jyc>=jpbeg .AND. jyc<=jpend) THEN
        IF(ixc <=ipend .AND. ixc >=ipbeg) rhsp(ixc,  jyc)=rhsp(ixc,  jyc)-pcdxx(1)*dxi*dxi
        IF(ixc1<=ipend .AND. ixc1>=ipbeg) rhsp(ixc+1,jyc)=rhsp(ixc+1,jyc)-pcdxx(2)*dxi*dxi
     ENDIF
     dpdx = ucpdx
  ENDIF

END SUBROUTINE correction_pressure
