SUBROUTINE correction_velocity(index,flag,krk,dpdx)

  USE stretchxy
  USE para
  USE Lagrange
  USE Euler

  INTEGER:: ixc,ixc1,ixf,ixf1,jyc,jyc1,jyf,jyf1,myvertex
  INTEGER,INTENT(in):: index,flag,krk
  DOUBLE PRECISION,INTENT(in):: dpdx
  DOUBLE PRECISION:: ucudx,vcdxx(2),vcdx,ucdxx(2)
  DOUBLE PRECISION:: vcvdy,vcdyy(2),ucdy,ucdyy(2)
  DOUBLE PRECISION:: xint,yint,signx,signy
  DOUBLE PRECISION:: ui,vi,djc,ddjc
  DOUBLE PRECISION:: ujcdxi,ujcdxi2,ujcdeta,ujcdeta2
  DOUBLE PRECISION:: vjcdxi,vjcdxi2,vjcdeta,vjcdeta2
  DOUBLE PRECISION:: ucpdx,vcpdy

  ! xc
  IF(flag == 1) THEN

     yint = xc_int_info(2, index)
     signx = xc_int_info(3, index)
     signy = xc_int_info(4, index)

     ixc = xc_int_vertex(2, index)
     jyc = xc_int_vertex(3, index)
     jyf = xc_int_vertex(4, index)

     ui=xc_int_info(5,index)
     vi=xc_int_info(6,index)

     jyf1=jyf+1
     jyc1=jyc+1

     vjcdeta=vjcxc(2,index)*dydeta(yint)
     vjcdeta2=vjcxc(4,index)*dydeta(yint)*dydeta(yint)+vjcxc(2,index)*d2ydeta2(yint)

     djc=2.0d0*vi*vjcdeta
     ddjc=2.0d0*vi*vjcdeta2 + &
          2.0d0*(xc_int_info(10,index)*xc_int_info(10,index)- &
          xc_int_info(11,index)*xc_int_info(11,index))*dydeta(yint)*dydeta(yint)

     IF(jyc == jyf) THEN
        vcvdy=-1.0d0/deta*detady(yint)*(djc*(etac(jyc1)-y2eta(yint))+ &
             0.5d0*ddjc*(etac(jyc1)-y2eta(yint))**2.0d0)
        vcdyy(1)=(-detady(yint)*detady(yint)/deta/deta-d2etady2(yint)/2.0d0/deta)*&
             (vjcdeta*(etaf(jyf1)-y2eta(yint))+0.5d0*vjcdeta2*(etaf(jyf1)-y2eta(yint))**2.0d0)
        vcdyy(2)=( detady(yint)*detady(yint)/deta/deta-d2etady2(yint)/2.0d0/deta)*&
             (vjcdeta*(etaf(jyf) -y2eta(yint))+0.5d0*vjcdeta2*(etaf(jyf) -y2eta(yint))**2.0d0)
     ELSE
        vcvdy=-1.0d0/deta*detady(yint)*(djc*(etac(jyc)-y2eta(yint))+ &
             0.5d0*ddjc*(etac(jyc)-y2eta(yint))**2.0d0)
        vcdyy(1)=(-detady(yint)*detady(yint)/deta/deta-d2etady2(yint)/2.0d0/deta)*&
             (vjcdeta*(etaf(jyf1)-y2eta(yint))+0.5d0*vjcdeta2*(etaf(jyf1)-y2eta(yint))**2.0d0)
        vcdyy(2)=( detady(yint)*detady(yint)/deta/deta-d2etady2(yint)/2.0d0/deta)*&
             (vjcdeta*(etaf(jyf) -y2eta(yint))+0.5d0*vjcdeta2*(etaf(jyf) -y2eta(yint))**2.0d0)
     ENDIF
     vcpdy = dpdx
      IF(ixc<=ivend-1 .AND. ixc>= ivbeg+1) THEN
         IF( krk.EQ.1 ) THEN
            IF(jyc <=jvend-1  .AND. jyc>=jvbeg+1)  rhsv1(ixc,jyc)  = rhsv1(ixc,jyc)  - vcvdy - vcpdy
            IF(jyf <=jvend-1  .AND. jyf>=jvbeg+1)  rhsv1(ixc,jyf)  = rhsv1(ixc,jyf)  + vcdyy(1)/Re
            IF(jyf1 <=jvend-1 .AND. jyf1>=jvbeg+1) rhsv1(ixc,jyf1) = rhsv1(ixc,jyf1) + vcdyy(2)/Re
         ELSEIF ( krk.EQ.2 ) THEN
            IF(jyc  <=jvend-1 .AND. jyc>=jvbeg+1)  rhsv2(ixc,jyc)  = rhsv2(ixc,jyc)  - vcvdy - vcpdy
            IF(jyf  <=jvend-1 .AND. jyf>=jvbeg+1)  rhsv2(ixc,jyf)  = rhsv2(ixc,jyf)  + vcdyy(1)/Re
            IF(jyf1 <=jvend-1 .AND. jyf1>=jvbeg+1) rhsv2(ixc,jyf1) = rhsv2(ixc,jyf1) + vcdyy(2)/Re
         ELSEIF ( krk.EQ.3 ) THEN
            IF(jyc <=jvend-1 .AND.  jyc>=jvbeg+1)  rhsv3(ixc,jyc)  = rhsv3(ixc,jyc)  - vcvdy - vcpdy
            IF(jyf <=jvend-1 .AND.  jyf>=jvbeg+1)  rhsv3(ixc,jyf)  = rhsv3(ixc,jyf)  + vcdyy(1)/Re
            IF(jyf1 <=jvend-1 .AND. jyf1>=jvbeg+1) rhsv3(ixc,jyf1) = rhsv3(ixc,jyf1) + vcdyy(2)/Re
         ELSE
            IF(jyc <=jvend-1 .AND.  jyc>=jvbeg+1)  rhsv4(ixc,jyc)  = rhsv4(ixc,jyc)  - vcvdy - vcpdy
            IF(jyf <=jvend-1 .AND.  jyf>=jvbeg+1)  rhsv4(ixc,jyf)  = rhsv4(ixc,jyf)  + vcdyy(1)/Re
            IF(jyf1 <=jvend-1 .AND. jyf1>=jvbeg+1) rhsv4(ixc,jyf1) = rhsv4(ixc,jyf1) + vcdyy(2)/Re
         ENDIF
      ENDIF
    !  IF( krk.EQ.1 ) THEN
    !     rhsv1(ixc,jyc)  = rhsv1(ixc,jyc)  - vcvdy - vcpdy
    !     rhsv1(ixc,jyf)  = rhsv1(ixc,jyf)  + vcdyy(1)/Re
    !     rhsv1(ixc,jyf1) = rhsv1(ixc,jyf1) + vcdyy(2)/Re
    !  ELSEIF ( krk.EQ.2 ) THEN
    !     rhsv2(ixc,jyc)  = rhsv2(ixc,jyc)  - vcvdy - vcpdy
    !     rhsv2(ixc,jyf)  = rhsv2(ixc,jyf)  + vcdyy(1)/Re
    !     rhsv2(ixc,jyf1) = rhsv2(ixc,jyf1) + vcdyy(2)/Re
    !  ELSEIF ( krk.EQ.3 ) THEN
    !     rhsv3(ixc,jyc)  = rhsv3(ixc,jyc)  - vcvdy - vcpdy
    !     rhsv3(ixc,jyf)  = rhsv3(ixc,jyf)  + vcdyy(1)/Re
    !     rhsv3(ixc,jyf1) = rhsv3(ixc,jyf1) + vcdyy(2)/Re
    !  ELSE
    !     rhsv4(ixc,jyc)  = rhsv4(ixc,jyc)  - vcvdy - vcpdy
    !     rhsv4(ixc,jyf)  = rhsv4(ixc,jyf)  + vcdyy(1)/Re
    !     rhsv4(ixc,jyf1) = rhsv4(ixc,jyf1) + vcdyy(2)/Re
    !  ENDIF
  ENDIF

  ! xf
  IF(flag == 2) THEN

     myvertex=xf_int_vertex(1,index)
     yint = xf_int_info(2, index)
     signx = xf_int_info(3, index)
     signy = xf_int_info(4, index)

     ixf = xf_int_vertex(2, index)
     jyc = xf_int_vertex(3, index)
     jyf = xf_int_vertex(4, index)

     ui=xf_int_info(5,index)
     vi=xf_int_info(6,index)

     jyf1=jyf+1
     jyc1=jyc+1

     ujcdeta=ujcxf(2,index)*dydeta(yint)
     vjcdeta=vjcxf(2,index)*dydeta(yint)
     ujcdeta2=ujcxf(4,index)*dydeta(yint)*dydeta(yint)+ujcxf(2,index)*d2ydeta2(yint)
     vjcdeta2=vjcxf(4,index)*dydeta(yint)*dydeta(yint)+vjcxf(2,index)*d2ydeta2(yint)

     djc=ui*vjcdeta+vi*ujcdeta
     ddjc=ui*vjcdeta2 + vi*ujcdeta2 + &
          2.0d0*(xf_int_info(7,index)*xf_int_info(9,index)- &
          xf_int_info(8,index)*xf_int_info(10,index))*dydeta(yint)*dydeta(yint)

     IF(jyc == jyf) THEN
        ucdy=-1.0d0/deta*detady(yint)*(djc*(etaf(jyf)-y2eta(yint))+ &
             0.5d0*ddjc*(etaf(jyf)-y2eta(yint))**2.0d0)
        ucdyy(1)=(-detady(yint)*detady(yint)/deta/deta-d2etady2(yint)/2.0d0/deta)*&
             (ujcdeta*(etac(jyc1)-y2eta(yint))+0.5d0*ujcdeta2*(etac(jyc1)-y2eta(yint))**2.0d0)
        ucdyy(2)=( detady(yint)*detady(yint)/deta/deta-d2etady2(yint)/2.0d0/deta)*&
             (ujcdeta*(etac(jyc) -y2eta(yint))+0.5d0*ujcdeta2*(etac(jyc) -y2eta(yint))**2.0d0)
     ELSE
        ucdy=-1.0d0/deta*detady(yint)*(djc*(etaf(jyf1)-y2eta(yint))+ &
             0.5d0*ddjc*(etaf(jyf1)-y2eta(yint))**2.0d0)
        ucdyy(1)=(-detady(yint)*detady(yint)/deta/deta-d2etady2(yint)/2.0d0/deta)*&
             (ujcdeta*(etac(jyc1)-y2eta(yint))+0.5d0*ujcdeta2*(etac(jyc1)-y2eta(yint))**2.0d0)
        ucdyy(2)=( detady(yint)*detady(yint)/deta/deta-d2etady2(yint)/2.0d0/deta)*&
             (ujcdeta*(etac(jyc) -y2eta(yint))+0.5d0*ujcdeta2*(etac(jyc) -y2eta(yint))**2.0d0)
     ENDIF
      IF(ixf<=iuend-1 .AND. ixf>= iubeg+1) THEN
         IF( krk.EQ.1 ) THEN
            IF(jyf1 <=juend-1 .AND. jyf1>=jubeg+1) rhsu1(ixf,jyf1) = rhsu1(ixf,jyf1) - ucdy
            IF(jyc <=juend-1 .AND. jyc>=jubeg+1)   rhsu1(ixf,jyc)  = rhsu1(ixf,jyc)  + ucdyy(1)*Re1
            IF(jyc1 <=juend-1 .AND. jyc1>=jubeg+1) rhsu1(ixf,jyc1) = rhsu1(ixf,jyc1) + ucdyy(2)*Re1
         ELSEIF ( krk.EQ.2 ) THEN
            IF(jyf1 <=juend-1 .AND. jyf1>=jubeg+1) rhsu2(ixf,jyf1) = rhsu2(ixf,jyf1) - ucdy
            IF(jyc <=juend-1 .AND. jyc>=jubeg+1)   rhsu2(ixf,jyc)  = rhsu2(ixf,jyc)  + ucdyy(1)*Re1
            IF(jyc1 <=juend-1 .AND. jyc1>=jubeg+1) rhsu2(ixf,jyc1) = rhsu2(ixf,jyc1) + ucdyy(2)*Re1
         ELSEIF ( krk.EQ.3 ) THEN
            IF(jyf1 <=juend-1 .AND. jyf1>=jubeg+1) rhsu3(ixf,jyf1) = rhsu3(ixf,jyf1) - ucdy
            IF(jyc <=juend-1 .AND. jyc>=jubeg+1)   rhsu3(ixf,jyc)  = rhsu3(ixf,jyc)  + ucdyy(1)*Re1
            IF(jyc1 <=juend-1 .AND. jyc1>=jubeg+1) rhsu3(ixf,jyc1) = rhsu3(ixf,jyc1) + ucdyy(2)*Re1
         ELSE
            IF(jyf1 <=juend-1 .AND. jyf1>=jubeg+1) rhsu4(ixf,jyf1) = rhsu4(ixf,jyf1) - ucdy
            IF(jyc <=juend-1 .AND. jyc>=jubeg+1)   rhsu4(ixf,jyc)  = rhsu4(ixf,jyc)  + ucdyy(1)*Re1
            IF(jyc1 <=juend-1 .AND. jyc1>=jubeg+1) rhsu4(ixf,jyc1) = rhsu4(ixf,jyc1) + ucdyy(2)*Re1
         ENDIF
      ENDIF
     !
    !  IF( krk.EQ.1 ) THEN
    !     rhsu1(ixf,jyf1) = rhsu1(ixf,jyf1) - ucdy
    !     rhsu1(ixf,jyc)  = rhsu1(ixf,jyc)  + ucdyy(1)*Re1
    !     rhsu1(ixf,jyc1) = rhsu1(ixf,jyc1) + ucdyy(2)*Re1
    !  ELSEIF ( krk.EQ.2 ) THEN
    !     rhsu2(ixf,jyf1) = rhsu2(ixf,jyf1) - ucdy
    !     rhsu2(ixf,jyc)  = rhsu2(ixf,jyc)  + ucdyy(1)*Re1
    !     rhsu2(ixf,jyc1) = rhsu2(ixf,jyc1) + ucdyy(2)*Re1
    !  ELSEIF ( krk.EQ.3 ) THEN
    !     rhsu3(ixf,jyf1) = rhsu3(ixf,jyf1) - ucdy
    !     rhsu3(ixf,jyc)  = rhsu3(ixf,jyc)  + ucdyy(1)*Re1
    !     rhsu3(ixf,jyc1) = rhsu3(ixf,jyc1) + ucdyy(2)*Re1
    !  ELSE
    !     rhsu4(ixf,jyf1) = rhsu4(ixf,jyf1) - ucdy
    !     rhsu4(ixf,jyc)  = rhsu4(ixf,jyc)  + ucdyy(1)*Re1
    !     rhsu4(ixf,jyc1) = rhsu4(ixf,jyc1) + ucdyy(2)*Re1
    !  ENDIF
  ENDIF

  ! yc
  IF(flag == 3) THEN

     xint = yc_int_info(2, index)
     signx = yc_int_info(3, index)
     signy = yc_int_info(4, index)

     jyc = yc_int_vertex(2, index)
     ixc = yc_int_vertex(3, index)
     ixf = yc_int_vertex(4, index)

     ui = yc_int_info(5,index)
     vi = yc_int_info(6,index)

     ixf1=ixf+1
     ixc1=ixc+1

     ujcdxi=ujcyc(1,index)*dxdxi(xint)
     ujcdxi2=ujcyc(3,index)*dxdxi(xint)*dxdxi(xint)+ujcyc(1,index)*d2xdxi2(xint)

     djc=2.0d0*ui*ujcdxi
     ddjc=2.0d0*ui*ujcdxi2 + &
          2.0d0*(yc_int_info(7,index)*yc_int_info(7,index)- &
          yc_int_info(8,index)*yc_int_info(8,index))*dxdxi(xint)*dxdxi(xint)

     IF(ixc == ixf) THEN
        ucudx=-1.0d0/dxi*dxidx(xint)*(djc*(xic(ixc1)-x2xi(xint))+ &
             0.5d0*ddjc*(xic(ixc1)-x2xi(xint))**2.0d0)
        ucdxx(1)=(-dxidx(xint)*dxidx(xint)/dxi/dxi-d2xidx2(xint)/2.0d0/dxi)*&
             (ujcdxi*(xif(ixf1)-x2xi(xint))+0.5d0*ujcdxi2*(xif(ixf1)-x2xi(xint))**2.0d0)
        ucdxx(2)=( dxidx(xint)*dxidx(xint)/dxi/dxi-d2xidx2(xint)/2.0d0/dxi)*&
             (ujcdxi*(xif(ixf) -x2xi(xint))+0.5d0*ujcdxi2*(xif(ixf) -x2xi(xint))**2.0d0)
     ELSE
        ucudx=-1.0d0/dxi*dxidx(xint)*(djc*(xic(ixc)-x2xi(xint))+ &
             0.5d0*ddjc*(xic(ixc)-x2xi(xint))**2.0d0)
        ucdxx(1)=(-dxidx(xint)*dxidx(xint)/dxi/dxi-d2xidx2(xint)/2.0d0/dxi)*&
             (ujcdxi*(xif(ixf1)-x2xi(xint))+0.5d0*ujcdxi2*(xif(ixf1)-x2xi(xint))**2.0d0)
        ucdxx(2)=( dxidx(xint)*dxidx(xint)/dxi/dxi-d2xidx2(xint)/2.0d0/dxi)*&
             (ujcdxi*(xif(ixf) -x2xi(xint))+0.5d0*ujcdxi2*(xif(ixf) -x2xi(xint))**2.0d0)
     ENDIF
     ucpdx = dpdx

     IF(jyc <=juend-1 .AND. jyc>=jubeg+1) THEN
        IF( krk.EQ.1 ) THEN
           IF(ixc<=iuend-1 .AND. ixc>= iubeg+1) rhsu1(ixc,jyc)  = rhsu1(ixc,jyc)  - ucudx - ucpdx
           IF(ixf<=iuend-1 .AND. ixf>= iubeg+1) rhsu1(ixf,jyc)  = rhsu1(ixf,jyc)  + ucdxx(1)*Re1
           IF(ixf1<=iuend-1 .AND. ixf1>= iubeg+1) rhsu1(ixf1,jyc) = rhsu1(ixf1,jyc) + ucdxx(2)*Re1
        ELSEIF ( krk.EQ.2 ) THEN
           IF(ixc<=iuend-1 .AND. ixc>= iubeg+1) rhsu2(ixc,jyc)  = rhsu2(ixc,jyc)  - ucudx - ucpdx
           IF(ixf<=iuend-1 .AND. ixf>= iubeg+1) rhsu2(ixf,jyc)  = rhsu2(ixf,jyc)  + ucdxx(1)*Re1
           IF(ixf1<=iuend-1 .AND. ixf1>= iubeg+1) rhsu2(ixf1,jyc) = rhsu2(ixf1,jyc) + ucdxx(2)*Re1
        ELSEIF ( krk.EQ.3 ) THEN
           IF(ixc<=iuend-1 .AND. ixc>= iubeg+1) rhsu3(ixc,jyc)  = rhsu3(ixc,jyc)  - ucudx - ucpdx
           IF(ixf<=iuend-1 .AND. ixf>= iubeg+1) rhsu3(ixf,jyc)  = rhsu3(ixf,jyc)  + ucdxx(1)*Re1
           IF(ixf1<=iuend-1 .AND. ixf1>= iubeg+1) rhsu3(ixf1,jyc) = rhsu3(ixf1,jyc) + ucdxx(2)*Re1
        ELSE
           IF(ixc<=iuend-1 .AND. ixc>= iubeg+1) rhsu4(ixc,jyc)  = rhsu4(ixc,jyc)  - ucudx - ucpdx
           IF(ixf<=iuend-1 .AND. ixf>= iubeg+1) rhsu4(ixf,jyc)  = rhsu4(ixf,jyc)  + ucdxx(1)*Re1
           IF(ixf1<=iuend-1 .AND. ixf1>= iubeg+1) rhsu4(ixf1,jyc) = rhsu4(ixf1,jyc) + ucdxx(2)*Re1
        ENDIF
     ENDIF
    !  IF( krk.EQ.1 ) THEN
    !     rhsu1(ixc,jyc)  = rhsu1(ixc,jyc)  - ucudx - ucpdx
    !     rhsu1(ixf,jyc)  = rhsu1(ixf,jyc)  + ucdxx(1)*Re1
    !     rhsu1(ixf1,jyc) = rhsu1(ixf1,jyc) + ucdxx(2)*Re1
    !  ELSEIF ( krk.EQ.2 ) THEN
    !     rhsu2(ixc,jyc)  = rhsu2(ixc,jyc)  - ucudx - ucpdx
    !     rhsu2(ixf,jyc)  = rhsu2(ixf,jyc)  + ucdxx(1)*Re1
    !     rhsu2(ixf1,jyc) = rhsu2(ixf1,jyc) + ucdxx(2)*Re1
    !  ELSEIF ( krk.EQ.3 ) THEN
    !     rhsu3(ixc,jyc)  = rhsu3(ixc,jyc)  - ucudx - ucpdx
    !     rhsu3(ixf,jyc)  = rhsu3(ixf,jyc)  + ucdxx(1)*Re1
    !     rhsu3(ixf1,jyc) = rhsu3(ixf1,jyc) + ucdxx(2)*Re1
    !  ELSE
    !     rhsu4(ixc,jyc)  = rhsu4(ixc,jyc)  - ucudx - ucpdx
    !     rhsu4(ixf,jyc)  = rhsu4(ixf,jyc)  + ucdxx(1)*Re1
    !     rhsu4(ixf1,jyc) = rhsu4(ixf1,jyc) + ucdxx(2)*Re1
    !  ENDIF
  ENDIF

  ! yf
  IF(flag == 4) THEN

     xint = yf_int_info(2, index)
     signx = yf_int_info(3, index)
     signy = yf_int_info(4, index)

     jyf = yf_int_vertex(2, index)
     ixc = yf_int_vertex(3, index)
     ixf = yf_int_vertex(4, index)

     ui = yf_int_info(5,index)
     vi = yf_int_info(6,index)

     ixf1=ixf+1
     ixc1=ixc+1

     ujcdxi=ujcyf(1,index)*dxdxi(xint)
     vjcdxi=vjcyf(1,index)*dxdxi(xint)
     ujcdxi2=ujcyf(3,index)*dxdxi(xint)*dxdxi(xint)+ujcyf(1,index)*d2xdxi2(xint)
     vjcdxi2=vjcyf(3,index)*dxdxi(xint)*dxdxi(xint)+vjcyf(1,index)*d2xdxi2(xint)

     djc=ui*vjcdxi+vi*ujcdxi
     ddjc=ui*vjcdxi2 + vi*ujcdxi2 + &
          2.0d0*(yf_int_info(7,index)*yf_int_info(9,index)- &
          yf_int_info(8,index)*yf_int_info(10,index))*dxdxi(xint)*dxdxi(xint)

     IF(ixc == ixf) THEN
        vcdx=-1.0d0/dxi*dxidx(xint)*(djc*(xif(ixf)-x2xi(xint)) + &
             0.5d0*ddjc*(xif(ixf)-x2xi(xint))**2.0d0)
        vcdxx(1)=(-dxidx(xint)*dxidx(xint)/dxi/dxi-d2xidx2(xint)/2.0d0/dxi)*&
             (vjcdxi*(xic(ixc1)-x2xi(xint))+0.5d0*vjcdxi2*(xic(ixc1)-x2xi(xint))**2.0d0)
        vcdxx(2)=( dxidx(xint)*dxidx(xint)/dxi/dxi-d2xidx2(xint)/2.0d0/dxi)*&
             (vjcdxi*(xic(ixc) -x2xi(xint))+0.5d0*vjcdxi2*(xic(ixc) -x2xi(xint))**2.0d0)
     ELSE
        vcdx=-1.0d0/dxi*dxidx(xint)*(djc*(xif(ixf1)-x2xi(xint)) + &
             0.5d0*ddjc*(xif(ixf1)-x2xi(xint))**2.0d0)
        vcdxx(1)=(-dxidx(xint)*dxidx(xint)/dxi/dxi-d2xidx2(xint)/2.0d0/dxi)*&
             (vjcdxi*(xic(ixc1)-x2xi(xint))+0.5d0*vjcdxi2*(xic(ixc1)-x2xi(xint))**2.0d0)
        vcdxx(2)=( dxidx(xint)*dxidx(xint)/dxi/dxi-d2xidx2(xint)/2.0d0/dxi)*&
             (vjcdxi*(xic(ixc) -x2xi(xint))+0.5d0*vjcdxi2*(xic(ixc) -x2xi(xint))**2.0d0)
     ENDIF
     IF(jyf <=jvend-1 .AND. jyf>=jvbeg+1) THEN
        IF( krk.EQ.1 ) THEN
           IF(ixf1<=ivend-1 .AND. ixf1>= ivbeg+1) rhsv1(ixf1,jyf) = rhsv1(ixf1,jyf) - vcdx
           IF(ixc<=ivend-1 .AND. ixc>= ivbeg+1)   rhsv1(ixc,jyf)  = rhsv1(ixc,jyf)  + vcdxx(1)*Re1
           IF(ixc1<=ivend-1 .AND. ixc1>= ivbeg+1) rhsv1(ixc1,jyf) = rhsv1(ixc1,jyf) + vcdxx(2)*Re1
        ELSEIF ( krk.EQ.2 ) THEN
           IF(ixf1<=ivend-1 .AND. ixf1>= ivbeg+1) rhsv2(ixf1,jyf) = rhsv2(ixf1,jyf) - vcdx
           IF(ixc<=ivend-1 .AND. ixc>= ivbeg+1)   rhsv2(ixc,jyf)  = rhsv2(ixc,jyf)  + vcdxx(1)*Re1
           IF(ixc1<=ivend-1 .AND. ixc1>= ivbeg+1) rhsv2(ixc1,jyf) = rhsv2(ixc1,jyf) + vcdxx(2)*Re1
        ELSEIF ( krk.EQ.3 ) THEN
           IF(ixf1<=ivend-1 .AND. ixf1>= ivbeg+1) rhsv3(ixf1,jyf) = rhsv3(ixf1,jyf) - vcdx
           IF(ixc<=ivend-1 .AND. ixc>= ivbeg+1)   rhsv3(ixc,jyf)  = rhsv3(ixc,jyf)  + vcdxx(1)*Re1
           IF(ixc1<=ivend-1 .AND. ixc1>= ivbeg+1) rhsv3(ixc1,jyf) = rhsv3(ixc1,jyf) + vcdxx(2)*Re1
        ELSE
           IF(ixf1<=ivend-1 .AND. ixf1>= ivbeg+1) rhsv4(ixf1,jyf) = rhsv4(ixf1,jyf) - vcdx
           IF(ixc<=ivend-1 .AND. ixc>= ivbeg+1)   rhsv4(ixc,jyf)  = rhsv4(ixc,jyf)  + vcdxx(1)*Re1
           IF(ixc1<=ivend-1 .AND. ixc1>= ivbeg+1) rhsv4(ixc1,jyf) = rhsv4(ixc1,jyf) + vcdxx(2)*Re1
        ENDIF
     ENDIF
    !  IF( krk.EQ.1 ) THEN
    !     rhsv1(ixf1,jyf) = rhsv1(ixf1,jyf) - vcdx
    !     rhsv1(ixc,jyf)  = rhsv1(ixc,jyf)  + vcdxx(1)*Re1
    !     rhsv1(ixc1,jyf) = rhsv1(ixc1,jyf) + vcdxx(2)*Re1
    !  ELSEIF ( krk.EQ.2 ) THEN
    !     rhsv2(ixf1,jyf) = rhsv2(ixf1,jyf) - vcdx
    !     rhsv2(ixc,jyf)  = rhsv2(ixc,jyf)  + vcdxx(1)*Re1
    !     rhsv2(ixc1,jyf) = rhsv2(ixc1,jyf) + vcdxx(2)*Re1
    !  ELSEIF ( krk.EQ.3 ) THEN
    !     rhsv3(ixf1,jyf) = rhsv3(ixf1,jyf) - vcdx
    !     rhsv3(ixc,jyf)  = rhsv3(ixc,jyf)  + vcdxx(1)*Re1
    !     rhsv3(ixc1,jyf) = rhsv3(ixc1,jyf) + vcdxx(2)*Re1
    !  ELSE
    !     rhsv4(ixf1,jyf) = rhsv4(ixf1,jyf) - vcdx
    !     rhsv4(ixc,jyf)  = rhsv4(ixc,jyf)  + vcdxx(1)*Re1
    !     rhsv4(ixc1,jyf) = rhsv4(ixc1,jyf) + vcdxx(2)*Re1
    !  ENDIF
  ENDIF


END SUBROUTINE correction_velocity
