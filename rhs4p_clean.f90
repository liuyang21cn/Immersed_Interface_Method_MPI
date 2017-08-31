SUBROUTINE rhs4p_clean

  USE mpi
  USE para
  USE Euler
  USE Lagrange

  DOUBLE PRECISION::sum1,sum2,sum_zz,sum_bz

  ! compatibility condition enforcement
  sum_zz=zero
  sum_bz=zero
  sum1=zero
  sum2=zero
  DO j=jpbeg,jpend
     DO i=ipbeg,ipend
        sum1 = sum1 + rhsp(i,j)*(ppp(i)*qqq(j))
        sum2 = sum2 + (ppp(i)*qqq(j))*(ppp(i)*qqq(j))
     ENDDO
  ENDDO
  CALL MPI_Allreduce(sum1,sum_bz,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm2d,ierr)
  CALL MPI_Allreduce(sum2,sum_zz,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm2d,ierr)
  DO j=jpbeg,jpend
     DO i=ipbeg,ipend
        rhsp(i,j)=rhsp(i,j)-sum_bz/sum_zz*(ppp(i)*qqq(j))
     ENDDO
  ENDDO

END SUBROUTINE rhs4p_clean
