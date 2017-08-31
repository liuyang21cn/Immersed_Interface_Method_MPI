SUBROUTINE mac_distribute

  USE Euler
  USE para
  USE Lagrange

  DOUBLE PRECISION, DIMENSION(1:2)::A,B,corner0,corner1
  DOUBLE PRECISION:: signx,signy,f(22)
  DOUBLE PRECISION:: xsm0,xsm1,ysm0,ysm1,xx,yy
  INTEGER:: index, myvertex

  ! xc
  DO index = 1, n_xc_int
     myvertex = xc_int_vertex(1,index)
     xx = xc_int_info(1, index)
     xsm0 = xs(myvertex)
     xsm1 = xs(myvertex+1)
     ysm0 = ys(myvertex)
     ysm1 = ys(myvertex+1)

     DO n=1,6
        ujcxc(n,index)=ujc(n,myvertex)+(ujc(n,myvertex+1)-ujc(n,myvertex))*&
             (xx-xsm0)/(xsm1-xsm0)
        vjcxc(n,index)=vjc(n,myvertex)+(vjc(n,myvertex+1)-vjc(n,myvertex))*&
             (xx-xsm0)/(xsm1-xsm0)
     ENDDO
     DO n=1,8
        pjcxc(n,index)=pjc(n,myvertex)+(pjc(n,myvertex+1)-pjc(n,myvertex))*&
             (xx-xsm0)/(xsm1-xsm0)
     ENDDO
  ENDDO

  ! xf
  DO index = 1,n_xf_int
     myvertex = xf_int_vertex(1,index)
     xx = xf_int_info(1, index)
     signx = xf_int_info(3,i)
     signy = xf_int_info(4,i)
     xsm0 = xs(myvertex)
     xsm1 = xs(myvertex+1)
     ysm0 = ys(myvertex)
     ysm1 = ys(myvertex+1)
     DO n=1,6
        ujcxf(n,index) = ujc(n,myvertex)+(ujc(n,myvertex+1)-ujc(n,myvertex))*&
             (xx-xsm0)/(xsm1-xsm0)
        vjcxf(n,index) = vjc(n,myvertex)+(vjc(n,myvertex+1)-vjc(n,myvertex))*&
             (xx-xsm0)/(xsm1-xsm0)
     ENDDO
  ENDDO

  ! yc
  DO index = 1, n_yc_int
     myvertex = yc_int_vertex(1,index)
     yy = yc_int_info(1, index)
     xsm0 = xs(myvertex)
     xsm1 = xs(myvertex+1)
     ysm0 = ys(myvertex)
     ysm1 = ys(myvertex+1)
     DO n=1,6
        ujcyc(n,index)=ujc(n,myvertex)+(ujc(n,myvertex+1)-ujc(n,myvertex))*&
             (yy-ysm0)/(ysm1-ysm0)
        vjcyc(n,index)=vjc(n,myvertex)+(vjc(n,myvertex+1)-vjc(n,myvertex))*&
             (yy-ysm0)/(ysm1-ysm0)
     ENDDO
     DO n=1,8
        pjcyc(n,index)=pjc(n,myvertex)+(pjc(n,myvertex+1)-pjc(n,myvertex))*&
             (yy-ysm0)/(ysm1-ysm0)
     ENDDO
  ENDDO

  ! yf
  DO index = 1,n_yf_int
     myvertex = yf_int_vertex(1,index)
     yy = yf_int_info(1, index)
     xsm0 = xs(myvertex)
     xsm1 = xs(myvertex+1)
     ysm0 = ys(myvertex)
     ysm1 = ys(myvertex+1)
     DO n=1,6
        ujcyf(n,index)=ujc(n,myvertex)+(ujc(n,myvertex+1)-ujc(n,myvertex))*&
             (yy-ysm0)/(ysm1-ysm0)
        vjcyf(n,index)=vjc(n,myvertex)+(vjc(n,myvertex+1)-vjc(n,myvertex))*&
             (yy-ysm0)/(ysm1-ysm0)
     ENDDO
  ENDDO

  !xc
  DO i = 1, n_xc_int
     signx = xc_int_info(3,i)
     signy = xc_int_info(4,i)

     ujcxc(1,i)=signx*ujcxc(1,i)
     vjcxc(1,i)=signx*vjcxc(1,i)
     ujcxc(3,i)=signx*ujcxc(3,i)
     vjcxc(3,i)=signx*vjcxc(3,i)
     ujcxc(5,i)=signx*ujcxc(5,i)
     vjcxc(5,i)=signx*vjcxc(5,i)

     ujcxc(2,i)=signy*ujcxc(2,i)
     vjcxc(2,i)=signy*vjcxc(2,i)
     ujcxc(4,i)=signy*ujcxc(4,i)
     vjcxc(4,i)=signy*vjcxc(4,i)
     ujcxc(6,i)=signy*ujcxc(6,i)
     vjcxc(6,i)=signy*vjcxc(6,i)
  ENDDO

  ! xf
  DO i = 1, n_xf_int
     signx = xf_int_info(3,i)
     signy = xf_int_info(4,i)

     ujcxf(1,i)=signx*ujcxf(1,i)
     vjcxf(1,i)=signx*vjcxf(1,i)
     ujcxf(3,i)=signx*ujcxf(3,i)
     vjcxf(3,i)=signx*vjcxf(3,i)
     ujcxf(5,i)=signx*ujcxf(5,i)
     vjcxf(5,i)=signx*vjcxf(5,i)

     ujcxf(2,i)=signy*ujcxf(2,i)
     vjcxf(2,i)=signy*vjcxf(2,i)
     ujcxf(4,i)=signy*ujcxf(4,i)
     vjcxf(4,i)=signy*vjcxf(4,i)
     ujcxf(6,i)=signy*ujcxf(6,i)
     vjcxf(6,i)=signy*vjcxf(6,i)
  ENDDO

  ! yc
  DO j = 1, n_yc_int
     signx = yc_int_info(3,j)
     signy = yc_int_info(4,j)
     myvertex=yc_int_vertex(1,j)
     ujcyc(1,j)=signx*ujcyc(1,j)
     vjcyc(1,j)=signx*vjcyc(1,j)
     ujcyc(3,j)=signx*ujcyc(3,j)
     vjcyc(3,j)=signx*vjcyc(3,j)
     ujcyc(5,j)=signx*ujcyc(5,j)
     vjcyc(5,j)=signx*vjcyc(5,j)

     ujcyc(2,j)=signy*ujcyc(2,j)
     vjcyc(2,j)=signy*vjcyc(2,j)
     ujcyc(4,j)=signy*ujcyc(4,j)
     vjcyc(4,j)=signy*vjcyc(4,j)
     ujcyc(6,j)=signy*ujcyc(5,j)
     vjcyc(6,j)=signy*vjcyc(6,j)
  ENDDO

  ! yf
  DO j = 1, n_yf_int
     signx = yf_int_info(3,j)
     signy = yf_int_info(4,j)

     ujcyf(1,j)=signx*ujcyf(1,j)
     vjcyf(1,j)=signx*vjcyf(1,j)
     ujcyf(3,j)=signx*ujcyf(3,j)
     vjcyf(3,j)=signx*vjcyf(3,j)
     ujcyf(5,j)=signx*ujcyf(5,j)
     vjcyf(5,j)=signx*vjcyf(5,j)

     ujcyf(2,j)=signy*ujcyf(2,j)
     vjcyf(2,j)=signy*vjcyf(2,j)
     ujcyf(4,j)=signy*ujcyf(4,j)
     vjcyf(4,j)=signy*vjcyf(4,j)
     ujcyf(6,j)=signy*ujcyf(6,j)
     vjcyf(6,j)=signy*vjcyf(6,j)
  ENDDO

END SUBROUTINE mac_distribute
