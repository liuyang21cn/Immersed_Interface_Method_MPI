SUBROUTINE rhs4p(fac)

  USE mpi
  USE para
  USE Euler
  USE stretchxy

  DOUBLE PRECISION:: fac

  CALL Divergence_reset
  DO j=jpbeg,jpend
     DO i=ipbeg,ipend
        sp=TWO*(ux(i,j)*vy(i,j)-uy(i,j)*vx(i,j))
        rhsp(i,j)=rhsp(i,j)+(dn(i,j)/(fac*dt) + sp &
             -TWO*(d(i,j)*d(i,j)+ucc(i,j)*(d(i+1,j)-d(i-1,j))*xicx(i)/twodxi    &
             +vcc(i,j)*(d(i,j+1)-d(i,j-1))*etacy(j)/twodeta) &
             +Re1*(xicx(i) *(xifx(i-1) *d(i-1,j)-(xifx(i-1) +xifx(i)) *d(i,j)   &
             +xifx(i) *d(i+1,j))/dxi/dxi     &
             +etacy(j)*(etafy(j-1)*d(i,j-1)-(etafy(j-1)+etafy(j))*d(i,j)   &
             +etafy(j)*d(i,j+1))/deta/deta))*dxi*dxi
     ENDDO
  ENDDO
  ! if periodic condition, comment pbc out
  CALL pbc
  CALL rhs4p_clean

END SUBROUTINE rhs4p
