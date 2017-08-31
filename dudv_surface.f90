SUBROUTINE dudv_surface(index, flag)

  USE stretchxy
  USE Euler
  USE para
  USE Lagrange

  INTEGER:: index,flag,ixc,ixf,jyc,jyf,myvertex,myobj
  DOUBLE PRECISION:: signx,signy
  DOUBLE PRECISION:: duxjc,duyjc,dvxjc,dvyjc
  DOUBLE PRECISION:: duyp,duyn,dvyp,dvyn,duxp,duxn,dvxp,dvxn
  DOUBLE PRECISION:: uu,vv,xx,yy

  ! xc
  IF(flag == 1) THEN

     signx = xc_int_info(3, index)
     signy = xc_int_info(4, index)
     myvertex = xc_int_vertex(1, index)
     myobj = xc_int_vertex(5, index)

     ! jump conditions
     duxjc=signy*signx*ujcxc(1,index)
     dvxjc=signy*signx*vjcxc(1,index)
     duyjc=ujcxc(2,index)
     dvyjc=vjcxc(2,index)

     duxn = zero
     dvxn =  thetat(myobj)
     duyn = -thetat(myobj)
     dvyn = zero

     duxp = duxjc + duxn
     dvxp = dvxjc + dvxn
     duyp = duyjc + duyn
     dvyp = dvyjc + dvyn

     xc_int_info(7, index)=duxp
     xc_int_info(8, index)=duxn
     xc_int_info(9, index)=dvxp
     xc_int_info(10,index)=dvxn
     xc_int_info(11,index)=duyp
     xc_int_info(12,index)=duyn
     xc_int_info(13,index)=dvyp
     xc_int_info(14,index)=dvyn
  ENDIF

  ! xf
  IF(flag == 2) THEN

     signx = xf_int_info(3, index)
     signy = xf_int_info(4, index)
     myvertex = xf_int_vertex(1, index)
     myobj = xf_int_vertex(5, index)

     ! jump conditions
     duyjc = ujcxf(2,index)
     dvyjc = vjcxf(2,index)

     duyn = -thetat(myobj)
     dvyn = zero

     duyn = duyp - duyjc
     dvyn = dvyp - dvyjc
     duyp = duyjc + duyn
     dvyp = dvyjc + dvyn

     xf_int_info(7, index)=duyp
     xf_int_info(8, index)=duyn
     xf_int_info(9, index)=dvyp
     xf_int_info(10,index)=dvyn

  ENDIF

  ! yc
  IF(flag == 3) THEN

     signx = yc_int_info(3, index)
     signy = yc_int_info(4, index)
     myvertex = yc_int_vertex(1, index)
     myobj = yc_int_vertex(5, index)

     ! jump conditions
     duxjc=ujcyc(1,index)
     dvxjc=vjcyc(1,index)

     duyjc=signx*signy*ujcyc(2,index)
     dvyjc=signx*signy*vjcyc(2,index)

     duxn = zero
     dvxn =  thetat(myobj)
     duyn = -thetat(myobj)
     dvyn = zero

     duxp = duxjc + duxn
     dvxp = dvxjc + dvxn
     duyp = duyjc + duyn
     dvyp = dvyjc + dvyn

     yc_int_info(7, index)=duxp
     yc_int_info(8, index)=duxn
     yc_int_info(9, index)=dvxp
     yc_int_info(10,index)=dvxn
     yc_int_info(11,index)=duyp
     yc_int_info(12,index)=duyn
     yc_int_info(13,index)=dvyp
     yc_int_info(14,index)=dvyn
  ENDIF

  ! yf
  IF(flag == 4) THEN

     signx = yf_int_info(3, index)
     signy = yf_int_info(4, index)
     myvertex = yf_int_vertex(1, index)
     myobj = yf_int_vertex(5, index)

     ! jump conditions
     duxjc = ujcyf(1,index)
     dvxjc = vjcyf(1,index)

     duxn = zero
     dvxn =  thetat(myobj)

     duxp = duxjc + duxn
     dvxp = dvxjc + dvxn

     yf_int_info(7, index)=duxp
     yf_int_info(8, index)=duxn
     yf_int_info(9, index)=dvxp
     yf_int_info(10,index)=dvxn
  ENDIF

END SUBROUTINE dudv_surface
