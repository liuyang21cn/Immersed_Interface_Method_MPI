SUBROUTINE streamFunction

  USE mpi
  USE para
  USE Euler
  USE Lagrange

  DOUBLE PRECISION::dpdx
  DOUBLE PRECISION::sum1,sum2,sum_zz,sum_bz

  ! rhsph is actually vorticity
  rhsph=zero
  CALL interpolateuv(0)
  DO j=jpbeg,jpend
     DO i=ipbeg,ipend
        rhsph(i,j)=rhsph(i,j) + xicx(i)/dxi*(vec(i,j)-vec(i-1,j)) &
             -(uce(i,j)-uce(i,j-1))*etacy(j)/deta
        ! rhsph(i,j) = rhsph(i,j) + xicx(i)/2.0d0/dxi*(v(i+1,j)-v(i-1,j)) &
        !  -(u(i,j+1)-u(i,j-1))*etacy(j)/2.0d0/deta
        ! rhsph(i,j)=rhsph(i,j) + xicx(i)*(vcc(i+1,j)-vcc(i-1,j))/twodxi &
        !  -etacy(j)*(ucc(i,j+1)-ucc(i,j-1))/twodeta
        rhsph(i,j)=rhsph(i,j)*dxi*dxi
     ENDDO
  ENDDO
  DO j=jpbeg,jpend
     DO i=ipbeg,ipend
        DO myobj=1,nobj4proc
           IF(iop(i,j)==obj4proc(myobj)) THEN
              rhsph(i,j)=2.0d0*thetat(myobj)*dxi*dxi
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  ! BCs for stream function
  IF(ipbcx0==2) THEN
     IF(ipbeg==1) THEN
        i=1
        DO j=jpbeg,jpend
           dpdx=vcc(1,j)
            ! rhsph(i,j)=rhsph(i,j)+2.0d0*xifx(i-1)*dxi*(dpdx)
        ENDDO
     ENDIF
  ENDIF

  IF(ipbcx1==2) THEN
     IF(ipend==nx) THEN
        i=nx
        DO j=jpbeg,jpend
           !  dpdx=vcc(nx,j)
           dpdx=xicx(i)*(xifx(i-1)*v(i-1,j)-(xifx(i-1)+xifx(i))*v(i,j)+xifx(i)*v(i+1,j))/dxi/dxi
          !  rhsph(i,j)=rhsph(i,j)+2.0d0*xifx(i)*dxi*(-dpdx)
        ENDDO
     ENDIF
  ENDIF

  IF(jpbcy0==2) THEN
     IF(jpbeg==1) THEN
        j=1
        DO i=ipbeg,ipend
           dpdy=-ucc(i,1)
            ! rhsph(i,j)=rhsph(i,j)+2.0d0*etafy(j-1)/deta*dxi*dxi*(dpdy)
        ENDDO
     ENDIF
  ENDIF

  IF(jpbcy1==2) THEN
     IF(jpend==ny) THEN
        j=ny
        DO i=ipbeg,ipend
           dpdy=-ucc(i,ny)
            ! rhsph(i,j)=rhsph(i,j)+2.0d0*etafy(j)/deta*dxi*dxi*(-dpdy)
        ENDDO
     ENDIF
  ENDIF

  ! clean rhs of stream function
  sum_zz=zero
  sum_bz=zero
  sum1=zero
  sum2=zero
  DO j=jpbeg,jpend
     DO i=ipbeg,ipend
        sum1 = sum1 + rhsph(i,j)*(ppp(i)*qqq(j))
        sum2 = sum2 + (ppp(i)*qqq(j))*(ppp(i)*qqq(j))
     ENDDO
  ENDDO
  CALL MPI_Allreduce(sum1,sum_bz,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm2d,ierr)
  CALL MPI_Allreduce(sum2,sum_zz,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm2d,ierr)
  DO j=jpbeg,jpend
     DO i=ipbeg,ipend
        rhsph(i,j)=rhsph(i,j)-sum_bz/sum_zz*(ppp(i)*qqq(j))
     ENDDO
  ENDDO
  CALL solvePoisson(0)

END SUBROUTINE streamfunction
