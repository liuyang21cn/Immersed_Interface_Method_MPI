SUBROUTINE jc_pressure(index,flag)

  USE para
  USE Euler
  USE Lagrange

  INTEGER:: index,flag,ixc,jyc,myvertex
  DOUBLE PRECISION:: signx,signy,r3
  DOUBLE PRECISION:: duyp,duyn,dvyp,dvyn,duxp,duxn,dvxp,dvxn

  ! xc
  IF(flag == 1) THEN
     myvertex=xc_int_vertex(1,index)
     signx = xc_int_info(3, index)
     signy = xc_int_info(4, index)

     duxp = xc_int_info(7, index)
     duxn = xc_int_info(8, index)
     dvxp = xc_int_info(9, index)
     dvxn = xc_int_info(10,index)
     duyp = xc_int_info(11,index)
     duyn = xc_int_info(12,index)
     dvyp = xc_int_info(13,index)
     dvyn = xc_int_info(14,index)

     r3=2.0d0*signy*((duxp*dvyp-duxn*dvyn)-(dvxp*duyp-dvxn*duyn))

     pjcxc(3,index)=pjcxc(7,index)+(pjcxc(3,index)-pjcxc(7,index))*r3
     pjcxc(4,index)=pjcxc(8,index)+(pjcxc(4,index)-pjcxc(8,index))*r3

     pjcxc(1,index)=signx*pjcxc(1,index)
     pjcxc(3,index)=signx*pjcxc(3,index)
     pjcxc(5,index)=signx*pjcxc(5,index)

     pjcxc(2,index)=signy*pjcxc(2,index)
     pjcxc(4,index)=signy*pjcxc(4,index)
     pjcxc(6,index)=signy*pjcxc(6,index)
  ENDIF

  ! yc
  IF(flag == 3) THEN
     myvertex=yc_int_vertex(1,index)

     signx = yc_int_info(3, index)
     signy = yc_int_info(4, index)

     duxp = yc_int_info(7, index)
     duxn = yc_int_info(8, index)
     dvxp = yc_int_info(9, index)
     dvxn = yc_int_info(10,index)
     duyp = yc_int_info(11,index)
     duyn = yc_int_info(12,index)
     dvyp = yc_int_info(13,index)
     dvyn = yc_int_info(14,index)

     r3=2.0d0*signx*((duxp*dvyp-duxn*dvyn)-(dvxp*duyp-dvxn*duyn))

     pjcyc(3,index)=pjcyc(7,index)+(pjcyc(3,index)-pjcyc(7,index))*r3
     pjcyc(4,index)=pjcyc(8,index)+(pjcyc(4,index)-pjcyc(8,index))*r3

     pjcyc(1,index)=signx*pjcyc(1,index)
     pjcyc(3,index)=signx*pjcyc(3,index)
     pjcyc(5,index)=signx*pjcyc(5,index)

     pjcyc(2,index)=signy*pjcyc(2,index)
     pjcyc(4,index)=signy*pjcyc(4,index)
     pjcyc(6,index)=signy*pjcyc(6,index)

  ENDIF


END SUBROUTINE jc_pressure
