SUBROUTINE correction_strain(index, flag)

  USE para
  USE Euler
  USE Lagrange
  USE stretchxy

  INTEGER:: flag,index,ixc,ixf,ixf1,jyc,jyf,jyf1
  DOUBLE PRECISION:: xint,yint,pcudy,pcvdy,pcudx,pcvdx
  DOUBLE PRECISION:: ujcdeta,ujcdeta2,ujcdxi,ujcdxi2
  DOUBLE PRECISION:: vjcdeta,vjcdeta2,vjcdxi,vjcdxi2

  ! xc
  IF(flag == 1) THEN
     yint = xc_int_info(2,index)
     ixc = xc_int_vertex(2,index)
     jyc = xc_int_vertex(3,index)
     jyf = xc_int_vertex(4,index)
     jyf1 = jyf+1
     myvertex=xc_int_vertex(1,index)

     ujcdeta=ujcxc(2,index)*dydeta(yint)
     vjcdeta=vjcxc(2,index)*dydeta(yint)
     ujcdeta2=ujcxc(4,index)*dydeta(yint)*dydeta(yint)+ujcxc(2,index)*d2ydeta2(yint)
     vjcdeta2=vjcxc(4,index)*dydeta(yint)*dydeta(yint)+vjcxc(2,index)*d2ydeta2(yint)

     IF(jyc==jyf) THEN
        pcudy=-1.0d0/deta*detady(yint)*(ujcdeta*(etaf(jyf)-y2eta(yint))+ &
             0.5d0*ujcdeta2*(etaf(jyf)-y2eta(yint))**2.0d0)
        pcvdy=-1.0d0/deta*detady(yint)*(vjcdeta*(etaf(jyf)-y2eta(yint))+ &
             0.5d0*vjcdeta2*(etaf(jyf)-y2eta(yint))**2.0d0)
     ELSE
        pcudy=-1.0d0/deta*detady(yint)*(ujcdeta*(etaf(jyf1)-y2eta(yint))+ &
             0.5d0*ujcdeta2*(etaf(jyf1)-y2eta(yint))**2.0d0)
        pcvdy=-1.0d0/deta*detady(yint)*(vjcdeta*(etaf(jyf1)-y2eta(yint))+ &
             0.5d0*vjcdeta2*(etaf(jyf1)-y2eta(yint))**2.0d0)
     ENDIF
     IF(ixc>=idbeg .AND. ixc<=idend) THEN
        IF(jyf1>=jdbeg .AND. jyf1<=jdend) THEN
           uy(ixc,jyf1) = uy(ixc,jyf1) + pcudy
           vy(ixc,jyf1) = vy(ixc,jyf1) + pcvdy
           d(ixc,jyf1)  = d(ixc,jyf1)  + pcvdy
        ENDIF
     ENDIF
     IF(ixc>=ipbeg .AND. ixc<=ipend) THEN
        IF(jyf1>=jpbeg .AND. jyf1<=jpend) THEN
          rhsph(ixc,jyf1)=rhsph(ixc,jyf1)-pcudy
        ENDIF
     ENDIF

  ENDIF

  ! yc
  IF(flag == 3) THEN
     xint = yc_int_info(2, index)
     jyc = yc_int_vertex(2, index)
     ixc = yc_int_vertex(3, index)
     ixf = yc_int_vertex(4, index)
     myvertex=yc_int_vertex(1,index)
     ixf1=ixf+1

     ujcdxi=ujcyc(1,index)*dxdxi(xint)
     vjcdxi=vjcyc(1,index)*dxdxi(xint)
     ujcdxi2=ujcyc(3,index)*dxdxi(xint)*dxdxi(xint)+ujcyc(1,index)*d2xdxi2(xint)
     vjcdxi2=vjcyc(3,index)*dxdxi(xint)*dxdxi(xint)+vjcyc(1,index)*d2xdxi2(xint)

     IF(ixc==ixf) THEN
        pcudx=-1.0d0/dxi*dxidx(xint)*(ujcdxi*(xif(ixf)-x2xi(xint))+ &
             0.5d0*ujcdxi2*(xif(ixf)-x2xi(xint))**2.0d0)
        pcvdx=-1.0d0/dxi*dxidx(xint)*(vjcdxi*(xif(ixf)-x2xi(xint))+ &
             0.5d0*vjcdxi2*(xif(ixf)-x2xi(xint))**2.0d0)
     ELSE
        pcudx=-1.0d0/dxi*dxidx(xint)*(ujcdxi*(xif(ixf1)-x2xi(xint))+ &
             0.5d0*ujcdxi2*(xif(ixf1)-x2xi(xint))**2.0d0)
        pcvdx=-1.0d0/dxi*dxidx(xint)*(vjcdxi*(xif(ixf1)-x2xi(xint))+ &
             0.5d0*vjcdxi2*(xif(ixf1)-x2xi(xint))**2.0d0)
     ENDIF
     IF(ixf1>=idbeg .AND. ixf1<=idend) THEN
        IF(jyc>=jdbeg .AND. jyc<=jdend) THEN
           ux(ixf1, jyc) = ux(ixf1,jyc) + pcudx
           vx(ixf1, jyc) = vx(ixf1,jyc) + pcvdx
           d(ixf1, jyc)  = d(ixf1, jyc) + pcudx
           rhsph(ixf1,jyc)=rhsph(ixf1,jyc)+pcvdx
        ENDIF
     ENDIF
     IF(ixf1>=ipbeg .AND. ixf1<=ipend) THEN
        IF(jyc>=jpbeg .AND. jyc<=jpend) THEN
           rhsph(ixf1,jyc)=rhsph(ixf1,jyc)+pcvdx
        ENDIF
     ENDIF
  ENDIF


END SUBROUTINE correction_strain
