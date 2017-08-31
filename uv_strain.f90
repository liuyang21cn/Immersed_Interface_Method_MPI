SUBROUTINE uv_strain

  USE stretchxy
  USE para
  USE Euler
  USE Lagrange

  DO j=jdbeg,jdend
     DO i=idbeg,idend
        !dudx
        ux(i,j)=ux(i,j)+xicx(i)*(u(i,j)-u(i-1,j))/dxi
        ! dvdx
        vx(i,j)=vx(i,j)+xicx(i)*(vcc(i+1,j)-vcc(i-1,j))/twodxi
        ! dudy
        uy(i,j)=uy(i,j)+etacy(j)*(ucc(i,j+1)-ucc(i,j-1))/twodeta
        !dvdy
        vy(i,j)=vy(i,j)+etacy(j)*(v(i,j)-v(i,j-1))/deta
        ! divergence
        d(i,j)=d(i,j)+xicx(i)*(u(i,j)-u(i-1,j))/dxi+etacy(j)*(v(i,j)-v(i,j-1))/deta
     ENDDO
  ENDDO

END SUBROUTINE uv_strain
