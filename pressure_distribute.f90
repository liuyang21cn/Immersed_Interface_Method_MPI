SUBROUTINE pressure_distribute

  USE Euler
  USE para
  USE Lagrange

  DOUBLE PRECISION, DIMENSION(1:2)::A,B,corner0,corner1
  DOUBLE PRECISION:: signx,signy,f(22)
  DOUBLE PRECISION:: xsm0,xsm1,ysm0,ysm1,xx,yy
  INTEGER:: int_idx, myvertex

  ! xc
  DO int_idx = 1, n_xc_int
     myvertex = xc_int_vertex(1,int_idx)
     xx = xc_int_info(1, int_idx)
     signx = xc_int_info(3,int_idx)
     signy = xc_int_info(4,int_idx)
     xsm0 = xs(myvertex)
     xsm1 = xs(myvertex+1)
     ysm0 = ys(myvertex)
     ysm1 = ys(myvertex+1)

     DO n=1,6
        pjcxc(n,int_idx)=pjc(n,myvertex)+(pjc(n,myvertex+1)-pjc(n,myvertex))*&
             (xx-xsm0)/(xsm1-xsm0)
     ENDDO
     pjcxc(1,int_idx)=signx*pjcxc(1,int_idx)
     pjcxc(3,int_idx)=signx*pjcxc(3,int_idx)
     pjcxc(5,int_idx)=signx*pjcxc(5,int_idx)

     pjcxc(2,int_idx)=signy*pjcxc(2,int_idx)
     pjcxc(4,int_idx)=signy*pjcxc(4,int_idx)
     pjcxc(6,int_idx)=signy*pjcxc(6,int_idx)
  ENDDO

  ! yc
  DO int_idx = 1, n_yc_int
     myvertex = yc_int_vertex(1,int_idx)
     yy = yc_int_info(1, int_idx)
     signx = yc_int_info(3,int_idx)
     signy = yc_int_info(4,int_idx)
     xsm0 = xs(myvertex)
     xsm1 = xs(myvertex+1)
     ysm0 = ys(myvertex)
     ysm1 = ys(myvertex+1)

     DO n=1,6
        pjcyc(n,int_idx)=pjc(n,myvertex)+(pjc(n,myvertex+1)-pjc(n,myvertex))*&
             (yy-ysm0)/(ysm1-ysm0)
     ENDDO
     pjcyc(1,int_idx)=signx*pjcyc(1,int_idx)
     pjcyc(3,int_idx)=signx*pjcyc(3,int_idx)
     pjcyc(5,int_idx)=signx*pjcyc(5,int_idx)

     pjcyc(2,int_idx)=signy*pjcyc(2,int_idx)
     pjcyc(4,int_idx)=signy*pjcyc(4,int_idx)
     pjcyc(6,int_idx)=signy*pjcyc(6,int_idx)

  ENDDO

END SUBROUTINE pressure_distribute
