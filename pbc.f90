SUBROUTINE pbc(fac)

  USE Euler
  USE Lagrange
  USE para

  DOUBLE PRECISION:: conv,visc,dpdx,dpdy,fac

  !c  dirichlet
  IF(ipbcx0==1) THEN
     IF(ipbeg==1) THEN
        DO j=jpbeg,jpend
           rhsp(ipbeg,j)=0.0d0
        ENDDO
     ENDIF
  ENDIF

  IF(ipbcx1==1) THEN
     IF(ipend==nx) THEN
        DO j=jpbeg,jpend
           rhsp(ipend,j)=0.0d0
        ENDDO
     ENDIF
  ENDIF

  IF(jpbcy0==1) THEN
     IF(jpbeg==1) THEN
        DO i=ipbeg,ipend
           rhsp(i,jpbeg)=0.0d0
        ENDDO
     ENDIF
  ENDIF

  IF(jpbcy1==1) THEN
     IF(jpend==ny) THEN
        DO i=ipbeg,ipend
           rhsp(i,jpend)=0.0d0
        ENDDO
     ENDIF
  ENDIF

  !c  newmann

  IF(ipbcx0==2) THEN
     IF(ipbeg==1) THEN
        i=1
        DO j=jpbeg,jpend
           dpdx=  xicx(i)/dxi/dxi/Re*( xifx(i-1)*(ucc(i+1,j) &
                + 2.0d0*dxi/deta*etacy(j)/xicx(i)*(v(i,j)-v(i,j-1))) &
                + xifx(i)*ucc(i+1,j) - (xifx(i-1)+xifx(i))*ucc(i,j)) &
                + etacy(j)/deta/deta/Re*(etafy(j)*ucc(i,j+1)-(etafy(j)+etafy(j-1))*ucc(i,j)+etafy(j-1)*ucc(i,j-1)) &
                - xicx(i)/dxi*(u(1,j)*u(1,j)-u(0,j)*u(0,j)) &
                - (uce(1,j)*v(1,j)-uce(1,j-1)*v(1,j-1))*etacy(j)/deta
           rhsp(i,j)=rhsp(i,j)+2.0d0*xifx(i-1)*dxi*(dpdx)
        ENDDO
     ENDIF
  ENDIF

  IF(ipbcx1==2) THEN
     IF(ipend==nx) THEN
        i=nx
        DO j=jpbeg,jpend
           dpdx=  xicx(i)/dxi/dxi/Re*( xifx(i)*(ucc(i-1,j)  &
                - 2.0d0*dxi/deta*etacy(j)/xicx(i)*(v(i,j)-v(i,j-1))) &
                + xifx(i-1)*ucc(i-1,j) - (xifx(i-1)+xifx(i))*ucc(i,j)) &
                + etacy(j)/deta/deta/Re*(etafy(j)*ucc(i,j+1)-(etafy(j)+etafy(j-1))*ucc(i,j)+etafy(j-1)*ucc(i,j-1)) &
                - xicx(i)/dxi*(u(nx,j)*u(nx,j)-u(nx-1,j)*u(nx-1,j)) &
                - (uce(nx,j)*v(nx,j)-uce(nx,j-1)*v(nx,j-1))*etacy(j)/deta
           rhsp(i,j)=rhsp(i,j)+2.0d0*xifx(i)*dxi*(-dpdx)
        ENDDO
     ENDIF
  ENDIF

  IF(jpbcy0==2) THEN
     IF(jpbeg==1) THEN
        j=1
        DO i=ipbeg,ipend
           dpdy= etacy(j)/deta/deta/Re*(etafy(j-1)*(vcc(i,j+1) &
                + 2.0d0*deta/dxi*xicx(i)/etacy(j)*(u(i,j)-u(i-1,j))) &
                -(etafy(j)+etafy(j-1))*vcc(i,j)+etafy(j)*vcc(i,j+1)) &
                +xicx(i)/dxi/dxi/Re*(xifx(i)*vcc(i+1,j)-(xifx(i)+xifx(i-1))*vcc(i,j)+xifx(i-1)*vcc(i-1,j)) &
                -etacy(j)/deta*(v(i,j)*v(i,j)-v(i,j-1)*v(i,j-1)) &
                -xicx(i)/dxi *(vec(i,j)*u(i,j)-vec(i-1,j)*u(i-1,j))
           rhsp(i,j)=rhsp(i,j)+2.0d0*etafy(j-1)/deta*dxi*dxi*(dpdy)
        ENDDO
     ENDIF
  ENDIF

  IF(jpbcy1==2) THEN
     IF(jpend==ny) THEN
        j=ny
        DO i=ipbeg,ipend
           dpdy= etacy(j)/deta/deta/Re*(etafy(j)*(vcc(i,j-1)  &
                - 2.0d0*deta/dxi*xicx(i)/etacy(j)*(u(i,j)-u(i-1,j))) &
                -(etafy(j)+etafy(j-1))*vcc(i,j)+etafy(j-1)*vcc(i,j-1)) &
                +xicx(i)/dxi/dxi/Re*(xifx(i)*vcc(i+1,j)-(xifx(i)+xifx(i-1))*vcc(i,j)+xifx(i-1)*vcc(i-1,j)) &
                -etacy(j)/deta*(v(i,j)*v(i,j)-v(i,j-1)*v(i,j-1)) &
                -xicx(i)/dxi *(vec(i,j)*u(i,j)-vec(i-1,j)*u(i-1,j))
           rhsp(i,j)=rhsp(i,j)+2.0d0*etafy(j)/deta*dxi*dxi*(-dpdy)
        ENDDO
     ENDIF
  ENDIF

END SUBROUTINE pbc
