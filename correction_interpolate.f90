SUBROUTINE correction_interpolate(index, flag)

  USE stretchxy
  USE para
  USE Euler
  USE Lagrange

  ! indices
  INTEGER, INTENT(in):: flag,index
  INTEGER:: is,js,myvertex,ixc,ixf,jyc,jyf
  INTEGER:: u_jc_index_y, u_jc_index_x, v_jc_index_y,v_jc_index_x
  ! jump contribution
  DOUBLE PRECISION:: vcixfyf,vciycxc,vciycxf
  DOUBLE PRECISION:: uciyfxc,uciyfxf,ucixcyc
  DOUBLE PRECISION:: xint,yint
  DOUBLE PRECISION:: ujcdeta,ujcdeta2,ujcdxi,ujcdxi2
  DOUBLE PRECISION:: vjcdeta,vjcdeta2,vjcdxi,vjcdxi2

  ! xc
  IF(flag == 1) THEN
     ixc = xc_int_vertex(2,index)
     jyc = xc_int_vertex(3,index)
     jyf = xc_int_vertex(4,index)
     yint = xc_int_info(2,index)
     myvertex=xc_int_vertex(1,index)

     ujcdeta =ujcxc(2,index)*dydeta(yint)
     vjcdeta =vjcxc(2,index)*dydeta(yint)
     ujcdeta2=ujcxc(4,index)*dydeta(yint)*dydeta(yint)+ujcxc(2,index)*d2ydeta2(yint)
     vjcdeta2=vjcxc(4,index)*dydeta(yint)*dydeta(yint)+vjcxc(2,index)*d2ydeta2(yint)

     IF( jyc == jyf ) THEN
        v_jc_index_y = jyc+1
        js = jyf
        vciycxc = 0.5d0*(vjcdeta*(etaf(js)-y2eta(yint))+ &
             0.5d0*vjcdeta2*(etaf(js)-y2eta(yint))**2.0d0)
        u_jc_index_y = jyf
        js = jyc+1
        uciyfxc =-0.5d0*(ujcdeta*(etac(js)-y2eta(yint))+ &
             0.5d0*ujcdeta2*(etac(js)-y2eta(yint))**2.0d0)
     ELSE
        v_jc_index_y = jyc
        js = jyf+1
        vciycxc =-0.5d0*(vjcdeta*(etaf(js)-y2eta(yint))+ &
             0.5d0*vjcdeta2*(etaf(js)-y2eta(yint))**2.0d0)
        u_jc_index_y = jyf+1
        js = jyc
        uciyfxc = 0.5d0*(ujcdeta*(etac(js)-y2eta(yint))+ &
             0.5d0*ujcdeta2*(etac(js)-y2eta(yint))**2.0d0)
     ENDIF
     IF(ixc<=ipend+1 .AND. ixc>= ipbeg-1) THEN
        IF(v_jc_index_y <=jvend .AND. v_jc_index_y>=jvbeg) THEN
           vcc(ixc, v_jc_index_y) = vcc(ixc, v_jc_index_y) + vciycxc
        ENDIF
     ENDIF
     IF(ixc<=ipend .AND. ixc>= ipbeg) THEN
        IF(u_jc_index_y <=jpend .AND. u_jc_index_y>=jpbeg-1) THEN
           uce(ixc, u_jc_index_y) = uce(ixc, u_jc_index_y) + uciyfxc
        ENDIF
     ENDIF
  ENDIF

  ! xf
  IF(flag == 2) THEN
     ixf = xf_int_vertex(2,index)
     jyc = xf_int_vertex(3,index)
     jyf = xf_int_vertex(4,index)
     yint = xf_int_info(2,index)
     myvertex=xf_int_vertex(1,index)

     ujcdeta=ujcxf(2,index)*dydeta(yint)
     vjcdeta=vjcxf(2,index)*dydeta(yint)
     ujcdeta2=ujcxf(4,index)*dydeta(yint)*dydeta(yint)+ujcxf(2,index)*d2ydeta2(yint)
     vjcdeta2=vjcxf(4,index)*dydeta(yint)*dydeta(yint)+vjcxf(2,index)*d2ydeta2(yint)

     IF( jyc == jyf) THEN
        u_jc_index_y = jyf
        js = jyc+1
        uciyfxf =-0.5d0*(ujcdeta*(etac(js)-y2eta(yint))+ &
             0.5d0*ujcdeta2*(etac(js)-y2eta(yint))**2.0d0)
        v_jc_index_y = jyc+1
        js = jyf
        vciycxf = 0.5d0*(vjcdeta*(etaf(js)-y2eta(yint))+ &
             0.5d0*vjcdeta2*(etaf(js)-y2eta(yint))**2.0d0)
     ELSE
        u_jc_index_y = jyf+1
        js = jyc
        uciyfxf = 0.5d0*(ujcdeta*(etac(js)-y2eta(yint))+ &
             0.5d0*ujcdeta2*(etac(js)-y2eta(yint))**2.0d0)
        v_jc_index_y =jyc
        js = jyf+1
        vciycxf =-0.5d0*(vjcdeta*(etaf(js)-y2eta(yint))+ &
             0.5d0*vjcdeta2*(etaf(js)-y2eta(yint))**2.0d0)
     ENDIF
     IF(ixf<=ipend .AND. ixf>= ipbeg-1) THEN
        IF(u_jc_index_y <=jpend .AND. u_jc_index_y>=jpbeg-1) THEN
           uee(ixf, u_jc_index_y) = uee(ixf, u_jc_index_y) + uciyfxf
        ENDIF
        IF(v_jc_index_y <=jpend .AND. v_jc_index_y>=jpbeg) THEN
           vec(ixf, v_jc_index_y) = vec(ixf, v_jc_index_y) + vciycxf
        ENDIF
     ENDIF
  ENDIF

  ! yc
  IF(flag == 3) THEN
     myvertex=yc_int_vertex(1,index)
     xint = yc_int_info(2, index)
     jyc = yc_int_vertex(2, index)
     ixc = yc_int_vertex(3, index)
     ixf = yc_int_vertex(4, index)

     ujcdxi=ujcyc(1,index)*dxdxi(xint)
     ujcdxi2=ujcyc(3,index)*dxdxi(xint)*dxdxi(xint)+ujcyc(1,index)*d2xdxi2(xint)

     IF(ixc == ixf) THEN
        u_jc_index_x = ixc+1
        is = ixf
        ucixcyc= 0.5d0*(ujcdxi*(xif(is)-x2xi(xint))+ &
             0.5d0*ujcdxi2*(xif(is)-x2xi(xint))**2.0d0)
     ELSE
        u_jc_index_x = ixc
        is = ixf+1
        ucixcyc=-0.5d0*(ujcdxi*(xif(is)-x2xi(xint))+ &
             0.5d0*ujcdxi2*(xif(is)-x2xi(xint))**2.0d0)
     ENDIF
     IF(u_jc_index_x<=iuend .AND. u_jc_index_x>= iubeg) THEN
        IF(jyc <=jpend+1 .AND. jyc>=jpbeg-1) THEN
           ucc(u_jc_index_x, jyc) = ucc(u_jc_index_x, jyc) + ucixcyc
        ENDIF
     ENDIF
  ENDIF

  ! yf
  IF(flag == 4) THEN
     xint = yf_int_info(2, index)
     jyf = yf_int_vertex(2, index)
     ixc = yf_int_vertex(3, index)
     ixf = yf_int_vertex(4, index)
     myvertex=yf_int_vertex(1,index)

     vjcdxi=vjcyf(1,index)*dxdxi(xint)
     vjcdxi2=vjcyf(3,index)*dxdxi(xint)*dxdxi(xint)+vjcyf(1,index)*d2xdxi2(xint)

     IF( ixc == ixf ) THEN
        v_jc_index_x = ixf
        is = ixc+1
        vcixfyf=-0.5d0*(vjcdxi*(xic(is)-x2xi(xint))+ &
             0.5d0*vjcdxi2*(xic(is)-x2xi(xint))**2.0d0)
     ELSE
        v_jc_index_x = ixf+1
        is = ixc
        vcixfyf= 0.5d0*(vjcdxi*(xic(is)-x2xi(xint))+ &
             0.5d0*vjcdxi2*(xic(is)-x2xi(xint))**2.0d0)
     ENDIF
     IF(v_jc_index_x<=ipend .AND. v_jc_index_x>= ipbeg-1) THEN
        IF(jyf <=jpend .AND. jyf>=jpbeg-1) THEN
           vee(v_jc_index_x,jyf) = vee(v_jc_index_x,jyf)+vcixfyf
        ENDIF
     ENDIF

  ENDIF

END SUBROUTINE correction_interpolate
