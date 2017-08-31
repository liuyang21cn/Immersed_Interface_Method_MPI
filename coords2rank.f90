SUBROUTINE coords2rank(mycoords,myrank)
  USE para
  USE MPI
  INTEGER, DIMENSION(1:2), INTENT(in):: mycoords
  INTEGER, INTENT(out):: myrank

  IF((ipbcx0/=0.AND.(mycoords(1)==-1.OR.mycoords(1)==px)).OR.&
       (jpbcy0/=0.AND.(mycoords(2)==-1.OR.mycoords(2)==py))) THEN
     myrank=MPI_PROC_NULL
     RETURN
  END IF

  CALL MPI_CART_RANK(comm2d,mycoords,myrank,ierr)

END SUBROUTINE coords2rank
